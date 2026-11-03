#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <limits>
#include <cctype>
#include <chrono>
#include <ctime>
#include <math.h>
#include <stdio.h>

#include "utility.h"
#include "functions.h"
#include "beamprune.h"
#include "ac_matcher.h"
#include "common.h"

using namespace std;
using namespace oligodesign;

double internal_time = 0.0;
double internal_time1 = 0.0;
bool show_processing = true;

bool F_speedup = false;

template <typename T>
bool update(std::unordered_map<int, State<T>> &stateMap, int index, double newscore, int index1, int index2, int index3, int q, Manner MANNER) {
    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        stateMap[index] = State<T>(newscore, index1, index2, index3, q, MANNER);
        return true;
    }
    return false;
}

template <typename T>
bool update(std::unordered_map<int, State<T>> &stateMap, int index, double newscore, int index1, int index2, int index3, int q, Manner MANNER, int first_pair_pos, int last_pair_pos) {
    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        stateMap[index] = State<T>(newscore, index1, index2, index3, q, MANNER, first_pair_pos, last_pair_pos);
        return true;
    }
    return false;
}

template <typename T>
bool update(std::array<State<T>, 257> &stateMap, std::array<bool, 257> &checkMap, int index, double newscore, int index1, int index2, int index3, int q, Manner MANNER) {
    if (!checkMap[index] || newscore > stateMap[index].score) {
        stateMap[index] = State<T>(newscore, index1, index2, index3, q, MANNER);
        checkMap[index] = true;
        return true;
    }
    return false;
}

/** *************************************************************************************************************** */

// 1. Calculate the N (unpaired)
template <typename T>
inline void compute_N(
    int i, int j, int nuci, int nucj, int index,
    double nuci_score, double nucj_score,
    const vector<int> *nucj_1_list,
    vector<int> &con_seq,
    std::unordered_map<int, State<T>> &bestN_j,
    const std::vector<std::unordered_map<int, State<T>>> &bestN)
{
    if (i == j) {
        double newscore = nuci_score;
        update(bestN_j, index, newscore, -1, -1, -1, -1, NONE);
    } else {
        // N(i,j) = N(i,j-1) + N(j,j)
        const std::unordered_map<int, State<T>> &bestN_j_1 = bestN[j - 1];
        for (auto &nucj_1 : *nucj_1_list) {
            if (!IsLegal(con_seq, nucj_1, nucj, j - 1)) {
                continue;
            }

            int index_j_1 = GetIndex(i, nuci, nucj_1); // the index of N(i, j-1)
            auto it = bestN_j_1.find(index_j_1);
            if (it != bestN_j_1.end()) {
                double newscore = it->second.score + nucj_score;
                update(bestN_j, index, newscore, index_j_1, -1, -1, -1, MANNER_N_N_F);
            }
        }
    }
}

// 2. Calculate the S (unpaired)
template <typename T>
inline void compute_S(
    int i, int j, int l, int nuci, int nucj, int index1,
    double nuci_score, double nucj_score,
    const vector<int> *nucj_1_list,
    vector<int> &con_seq,
    std::vector<std::array<State<T>, 257>> &bestSS_j,
    std::vector<std::array<bool, 257>> &bestSS_check_j,
    const std::vector<std::vector<std::array<State<T>, 257>>> &bestSS,
    const std::vector<std::vector<std::array<bool, 257>>> &bestSS_check)
{
    if (i == j) {
        double newscore = nuci_score;
        update(bestSS_j[1], bestSS_check_j[1], index1, newscore, -1, -1, -1, -1, NONE);
    } else {
        // S(i,j) = S(i,j-1) + S(j,j)
        if (l <= SINGLE_MAX_LEN - 1) { // -1是因為S的開頭會是 i或是q
            const std::array<State<T>, 257> &bestSS_j_1 = bestSS[j - 1][l];
            const std::array<bool, 257> &bestSS_check_j_1 = bestSS_check[j - 1][l];
            for (auto &nucj_1 : *nucj_1_list) {
                if (!IsLegal(con_seq, nucj_1, nucj, j - 1)) {
                    continue;
                }

                int index_j_1 = GetIndex(i, nuci, nucj_1); // the index of S(i, j-1)
                int index_ = GetIndex1(nuci, nucj_1);     // the index of S(i, j-1)
                if (bestSS_check_j_1[index_]) {
                    double newscore = bestSS_j_1[index_].score + nucj_score;
                    update(bestSS_j[l + 1], bestSS_check_j[l + 1], index1, newscore, index_j_1, -1, -1, -1, MANNER_S_S_F);
                }
            }
        }
    }
}

// 3. Hairpin
template <typename T>
inline void compute_hairpin(
    int i, int j, int l, int nuci, int nucj, int index,
    double nuci_score, double nucj_score,
    const vector<int> *nuci1_list,
    const vector<int> *nucj_1_list,
    vector<int> &con_seq,
    const std::vector<std::unordered_map<int, State<T>>> &bestN,
    std::unordered_map<int, State<T>> &bestSC_j_l1)
{
    if ((j - i > HAIRPIN_GAP) && _allowed_pairs[nuci][nucj]) {
        const std::unordered_map<int, State<T>> &bestN_j_1 = bestN[j - 1];
        for (auto &nuci1 : *nuci1_list) {
            if (!IsLegal(con_seq, nuci, nuci1, i)) {
                continue;
            }

            for (auto &nucj_1 : *nucj_1_list) {
                if (!IsLegal(con_seq, nucj_1, nucj, j - 1)) {
                    continue;
                }

                int index_i1_j_1 = GetIndex(i + 1, nuci1, nucj_1);
                auto it = bestN_j_1.find(index_i1_j_1);
                if (it != bestN_j_1.end()) {
                    double newscore = -MFE_scalar * v_score_hairpin(i, j, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj));
                    newscore += it->second.score + nuci_score + nucj_score;
                    update(bestSC_j_l1, index, newscore, index_i1_j_1, -1, -1, -1, MANNER_C_H);
                }
            }
        }
    }
}

// 4. Stacking
template <typename T>
inline void compute_stacking(
    int i, int j, int l, int nuci, int nucj, int index,
    double nuci_score, double nucj_score,
    const string &rna_seq,
    vector<int> &con_seq,
    const std::vector<std::vector<std::unordered_map<int, State<T>>>> &bestSC,
    std::unordered_map<int, State<T>> &bestSC_j_l1)
{
    if ((j - i > MIN_STACKPAIR_GAP) && _allowed_pairs[nuci][nucj]) {
        int p = i + 1;
        int q = j - 1;

        vector<int> &nucp_list = Base_table.at(rna_seq[p]);
        vector<int> &nucq_list = Base_table.at(rna_seq[q]);
        const std::unordered_map<int, State<T>> &bestSC_q = bestSC[q][q - p + 1];

        for (auto &nucp : nucp_list) {
            if (!IsLegal(con_seq, nuci, nucp, i)) {
                continue;
            }

            for (auto &nucq : nucq_list) {
                if (!_allowed_pairs[nucp][nucq]) {
                    continue;
                }

                if (!IsLegal(con_seq, nucq, nucj, q)) {
                    continue;
                }

                int index_p_q = GetIndex(p, nucp, nucq);
                auto it = bestSC_q.find(index_p_q);
                if (it == bestSC_q.end()) {
                    continue;
                }

                double newscore = -MFE_scalar * v_score_single_stacking(EnrBASE(nuci), EnrBASE(nucj), EnrBASE(nucp), EnrBASE(nucq));
                newscore += it->second.score + nuci_score + nucj_score;
                update(bestSC_j_l1, index, newscore, -1, index_p_q, -1, q, MANNER_C_I);
            }
        }
    }
}

// 5. Left Bulge
template <typename T>
inline void compute_left_bulge(
    int i, int j, int l, int nuci, int nucj, int index,
    double nuci_score, double nucj_score,
    const vector<int> *nucj_1_list,
    const string &rna_seq,
    vector<int> &con_seq,
    const std::vector<std::vector<std::array<State<T>, 257>>> &bestSS,
    const std::vector<std::vector<std::array<bool, 257>>> &bestSS_check,
    const std::vector<std::vector<std::unordered_map<int, State<T>>>> &bestSC,
    std::unordered_map<int, State<T>> &bestSC_j_l1)
{
    if ((j - i > MIN_STACKPAIR_GAP) && _allowed_pairs[nuci][nucj]) {
        int p = i + 1;
        vector<int> &nucp_list = Base_table.at(rna_seq[p]);

        int minq = MAX2(p + 1 + HAIRPIN_GAP, j - 1 - SINGLE_MAX_LEN);
        for (int q = minq; q < j; q++) {
            if (p == i + 1 && q == j - 1) { continue; }

            const std::array<State<T>, 257> &bestSS_j_1 = bestSS[j - 1][j - 1 - (q + 1) + 1];
            const std::array<bool, 257> &bestSS_check_j_1 = bestSS_check[j - 1][j - 1 - (q + 1) + 1];
            vector<int> &nucq_list = Base_table.at(rna_seq[q]);
            vector<int> &nucq1_list = Base_table.at(rna_seq[q + 1]);

            const std::unordered_map<int, State<T>> &bestSC_q = bestSC[q][q - p + 1];
            for (auto &nucp : nucp_list) {
                if (!IsLegal(con_seq, nuci, nucp, i)) {
                    continue;
                }

                for (auto &nucq : nucq_list) {
                    if (!_allowed_pairs[nucp][nucq]) {
                        continue;
                    }

                    if ((j - q == 1) && !IsLegal(con_seq, nucq, nucj, q)) {
                        continue;
                    }

                    int index_p_q = GetIndex(p, nucp, nucq);
                    auto it = bestSC_q.find(index_p_q);
                    if (it == bestSC_q.end()) {
                        continue;
                    }

                    double score_i1_p_1 = 0;
                    for (auto &nucj_1 : *nucj_1_list) {
                        if (!IsLegal(con_seq, nucj_1, nucj, j - 1)) {
                            continue;
                        }

                        for (auto &nucq1 : nucq1_list) {
                            if (!IsLegal(con_seq, nucq, nucq1, q)) {
                                continue;
                            }

                            int index_q1_j_ = GetIndex1(nucq1, nucj_1);
                            double score_q1_j_1 = 0;
                            if (j - q == 1) {
                                if ((nucj_1 != nucq) || (nucq1 != nucj)) {
                                    continue;
                                }
                                score_q1_j_1 = 0;
                            } else if (j - q == 2) {
                                if ((nucj_1 != nucq1) || !bestSS_check_j_1[index_q1_j_]) {
                                    continue;
                                }
                                score_q1_j_1 = bestSS_j_1[index_q1_j_].score;
                            } else {
                                if (!bestSS_check_j_1[index_q1_j_]) {
                                    continue;
                                }
                                score_q1_j_1 = bestSS_j_1[index_q1_j_].score;
                            }

                            int index_q1_j_1 = GetIndex(q + 1, nucq1, nucj_1);

                            double newscore = -MFE_scalar * v_score_single_bulge(i, j, p, q, EnrBASE(nuci), EnrBASE(nucp), EnrBASE(nucj_1), EnrBASE(nucj), EnrBASE(nuci), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucq1));

                            newscore += score_i1_p_1 + it->second.score + score_q1_j_1;
                            newscore += nuci_score + nucj_score;
                            update(bestSC_j_l1, index, newscore, -1, index_p_q, index_q1_j_1, q, MANNER_C_I);
                        }
                    }
                }
            }
        }
    }
}

// 6. Right Bulge
template <typename T>
inline void compute_right_bulge(
    int i, int j, int l, int nuci, int nucj, int index,
    double nuci_score, double nucj_score,
    const vector<int> *nuci1_list,
    const string &rna_seq,
    vector<int> &con_seq,
    const std::vector<std::vector<std::array<State<T>, 257>>> &bestSS,
    const std::vector<std::vector<std::array<bool, 257>>> &bestSS_check,
    const std::vector<std::vector<std::unordered_map<int, State<T>>>> &bestSC,
    std::unordered_map<int, State<T>> &bestSC_j_l1)
{
    if ((j - i > MIN_STACKPAIR_GAP) && _allowed_pairs[nuci][nucj]) {
        int q = j - 1;
        vector<int> &nucq_list = Base_table.at(rna_seq[q]);

        int maxp = MIN2(q - 1 - HAIRPIN_GAP, i + SINGLE_MAX_LEN + 1);
        for (int p = i + 1; p <= maxp; p++) {
            const std::array<State<T>, 257> &bestSS_p_1 = bestSS[p - 1][p - 1 - (i + 1) + 1];
            const std::array<bool, 257> &bestSS_check_p_1 = bestSS_check[p - 1][p - 1 - (i + 1) + 1];
            vector<int> &nucp_list = Base_table.at(rna_seq[p]);
            vector<int> &nucp_1_list = Base_table.at(rna_seq[p - 1]);

            const std::unordered_map<int, State<T>> &bestSC_q = bestSC[q][q - p + 1];
            for (auto &nucp : nucp_list) {
                if ((p - i == 1) && !IsLegal(con_seq, nuci, nucp, i)) {
                    continue;
                }

                for (auto &nucq : nucq_list) {
                    if (!_allowed_pairs[nucp][nucq]) {
                        continue;
                    }

                    if (!IsLegal(con_seq, nucq, nucj, q)) {
                        continue;
                    }

                    int index_p_q = GetIndex(p, nucp, nucq);
                    auto it = bestSC_q.find(index_p_q);
                    if (it == bestSC_q.end()) {
                        continue;
                    }

                    for (auto &nuci1 : *nuci1_list) {
                        if (!IsLegal(con_seq, nuci, nuci1, i)) {
                            continue;
                        }
                        for (auto &nucp_1 : nucp_1_list) {
                            if (!IsLegal(con_seq, nucp_1, nucp, p - 1)) {
                                continue;
                            }

                            int index_i1_p_ = GetIndex1(nuci1, nucp_1);
                            double score_i1_p_1 = 0;

                            if (p - i == 1) {
                                if ((nuci1 != nucp) || (nucp_1 != nuci)) {
                                    continue;
                                }
                                score_i1_p_1 = 0;
                            } else if (p - i == 2) {
                                if ((nuci1 != nucp_1) || !bestSS_check_p_1[index_i1_p_]) {
                                    continue;
                                }
                                score_i1_p_1 = bestSS_p_1[index_i1_p_].score;
                            } else {
                                if (!bestSS_check_p_1[index_i1_p_]) {
                                    continue;
                                }
                                score_i1_p_1 = bestSS_p_1[index_i1_p_].score;
                            }

                            double score_q1_j_1 = 0;
                            int index_i1_p_1 = GetIndex(i + 1, nuci1, nucp_1);

                            double newscore = -MFE_scalar * v_score_single_bulge(i, j, p, q, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucq), EnrBASE(nucj), EnrBASE(nucp_1), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucj));

                            newscore += score_i1_p_1 + it->second.score + score_q1_j_1;
                            newscore += nuci_score + nucj_score;
                            update(bestSC_j_l1, index, newscore, index_i1_p_1, index_p_q, -1, q, MANNER_C_I);
                        }
                    }
                }
            }
        }
    }
}

// 7. Internal Loop
template <typename T>
inline void compute_internal(
    int i, int j, int l, int nuci, int nucj, int index,
    double nuci_score, double nucj_score,
    const vector<int> *nuci1_list,
    const vector<int> *nucj_1_list,
    const string &rna_seq,
    vector<int> &con_seq,
    const std::vector<std::vector<std::array<State<T>, 257>>> &bestSS,
    const std::vector<std::vector<std::array<bool, 257>>> &bestSS_check,
    const std::vector<std::vector<std::unordered_map<int, State<T>>>> &bestSC,
    std::unordered_map<int, State<T>> &bestSC_j_l1)
{
    if ((j - i > MIN_STACKPAIR_GAP) && _allowed_pairs[nuci][nucj]) {
        int maxp = MIN2(j - 2 - HAIRPIN_GAP, i + SINGLE_MAX_LEN + 1);
        for (int p = i + 2; p <= maxp; p++) {
            const std::array<State<T>, 257> &bestSS_p_1 = bestSS[p - 1][p - 1 - (i + 1) + 1];
            const std::array<bool, 257> &bestSS_check_p_1 = bestSS_check[p - 1][p - 1 - (i + 1) + 1];

            vector<int> &nucp_list = Base_table.at(rna_seq[p]);
            vector<int> &nucp_1_list = Base_table.at(rna_seq[p - 1]);

            int minq = MAX2(p + 1 + HAIRPIN_GAP, j - i + p - SINGLE_MAX_LEN - 2);
            for (int q = minq; q < j - 1; q++) {
                const std::array<State<T>, 257> &bestSS_j_1 = bestSS[j - 1][j - 1 - (q + 1) + 1];
                const std::array<bool, 257> &bestSS_check_j_1 = bestSS_check[j - 1][j - 1 - (q + 1) + 1];

                vector<int> &nucq_list = Base_table.at(rna_seq[q]);
                vector<int> &nucq1_list = Base_table.at(rna_seq[q + 1]);

                const std::unordered_map<int, State<T>> &bestSC_q = bestSC[q][q - p + 1];

                for (auto &nucp : nucp_list) {
                    for (auto &nucq : nucq_list) {
                        if (!_allowed_pairs[nucp][nucq]) {
                            continue;
                        }

                        int index_p_q = GetIndex(p, nucp, nucq);
                        auto it = bestSC_q.find(index_p_q);
                        if (it == bestSC_q.end()) {
                            continue;
                        }

                        for (auto &nuci1 : *nuci1_list) {
                            if (!IsLegal(con_seq, nuci, nuci1, i))
                                continue;
                            for (auto &nucp_1 : nucp_1_list) {
                                if (!IsLegal(con_seq, nucp_1, nucp, p - 1))
                                    continue;

                                int index_i1_p_ = GetIndex1(nuci1, nucp_1);
                                double score_i1_p_1 = 0;

                                if (p - i == 2) {
                                    if ((nuci1 != nucp_1) || !bestSS_check_p_1[index_i1_p_]) {
                                        continue;
                                    }
                                    score_i1_p_1 = bestSS_p_1[index_i1_p_].score;
                                } else {
                                    if (!bestSS_check_p_1[index_i1_p_]) {
                                        continue;
                                    }
                                    score_i1_p_1 = bestSS_p_1[index_i1_p_].score;
                                }

                                for (auto &nucj_1 : *nucj_1_list) {
                                    if (!IsLegal(con_seq, nucj_1, nucj, j - 1)) {
                                        continue;
                                    }

                                    for (auto &nucq1 : nucq1_list) {
                                        if (!IsLegal(con_seq, nucq, nucq1, q)) {
                                            continue;
                                        }

                                        int index_q1_j_ = GetIndex1(nucq1, nucj_1);
                                        double score_q1_j_1 = 0;
                                        if (j - q == 2) {
                                            if ((nucj_1 != nucq1) || !bestSS_check_j_1[index_q1_j_]) {
                                                continue;
                                            }
                                            score_q1_j_1 = bestSS_j_1[index_q1_j_].score;
                                        } else {
                                            if (!bestSS_check_j_1[index_q1_j_]) {
                                                continue;
                                            }
                                            score_q1_j_1 = bestSS_j_1[index_q1_j_].score;
                                        }

                                        int index_i1_p_1 = GetIndex(i + 1, nuci1, nucp_1);
                                        int index_q1_j_1 = GetIndex(q + 1, nucq1, nucj_1);

                                        double newscore = -MFE_scalar * v_score_single_internal(i, j, p, q, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj), EnrBASE(nucp_1), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucq1));

                                        newscore += score_i1_p_1 + it->second.score + score_q1_j_1;
                                        newscore += nuci_score + nucj_score;
                                        update(bestSC_j_l1, index, newscore, index_i1_p_1, index_p_q, index_q1_j_1, q, MANNER_C_I);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// 8. Multi-loop Closing (C = M2 + ML_closing)
template <typename T>
inline void compute_multiloop_closed(
    int i, int j, int l, int nuci, int nucj, int index,
    double nuci_score, double nucj_score,
    const vector<int> *nuci1_list,
    const vector<int> *nucj_1_list,
    vector<int> &con_seq,
    const std::vector<std::unordered_map<int, State<T>>> &bestM2,
    std::unordered_map<int, State<T>> &bestSC_j_l1)
{
    if ((j - i > MIN_STACKPAIR_GAP * 2) && _allowed_pairs[nuci][nucj]) {
        const std::unordered_map<int, State<T>> &bestM2_j_1 = bestM2[j - 1];
        for (auto &nuci1 : *nuci1_list) {
            for (auto &nucj_1 : *nucj_1_list) {
                if (IsLegal(con_seq, nuci, nuci1, i) && IsLegal(con_seq, nucj_1, nucj, j - 1)) {
                    int index_i1_j_1 = GetIndex(i + 1, nuci1, nucj_1);
                    auto it = bestM2_j_1.find(index_i1_j_1);
                    if (it != bestM2_j_1.end()) {
                        double newscore = -MFE_scalar * v_score_multi(EnrBASE(nuci), EnrBASE(nucj));
                        newscore += it->second.score + nuci_score + nucj_score;
                        update(bestSC_j_l1, index, newscore, index_i1_j_1, -1, -1, -1, MANNER_C_M2);
                    }
                }
            }
        }
    }
}

// 9. Calculate M2 where i is unpaired and j is unpaired
template <typename T>
inline void compute_M2_unpaired(
    int i, int j, int nuci, int nucj, int index,
    double nuci_score, double nucj_score,
    const vector<int> *nuci1_list,
    const vector<int> *nucj_1_list,
    vector<int> &con_seq,
    const std::vector<std::unordered_map<int, State<T>>> &bestM2,
    std::unordered_map<int, State<T>> &bestM2_j)
{
    if (j - i > MIN_STACKPAIR_GAP * 2 - 1) {
        // M2(i,j) = M2(i,j-1) + F(j,j)
        const std::unordered_map<int, State<T>> &bestM2_j_1 = bestM2[j - 1];
        for (auto &nucj_1 : *nucj_1_list) {
            if (IsLegal(con_seq, nucj_1, nucj, j - 1)) {
                int index_j_1 = GetIndex(i, nuci, nucj_1);
                auto it = bestM2_j_1.find(index_j_1);
                if (it != bestM2_j_1.end()) {
                    int first_pair_pos = it->second.first_pair_pos;
                    int last_pair_pos = it->second.last_pair_pos;

                    if (j - last_pair_pos <= SINGLE_MAX_LEN) {
                        double newscore = it->second.score + nucj_score;
                        update(bestM2_j, index, newscore, index_j_1, -1, -1, -1, MANNER_M2_M2_F, first_pair_pos, last_pair_pos);
                    }
                }
            }
        }

        // M2(i,j) = F(i,i) + M2(i+1,j)
        for (auto &nuci1 : *nuci1_list) {
            if (IsLegal(con_seq, nuci, nuci1, i)) {
                int index_j = GetIndex(i + 1, nuci1, nucj);
                auto it = bestM2_j.find(index_j);
                if (it != bestM2_j.end()) {
                    int first_pair_pos = it->second.first_pair_pos;
                    int last_pair_pos = it->second.last_pair_pos;
                    if (first_pair_pos - i <= SINGLE_MAX_LEN) {
                        double newscore = it->second.score + nuci_score;
                        update(bestM2_j, index, newscore, index_j, -1, -1, -1, MANNER_M2_F_M2, first_pair_pos, last_pair_pos);
                    }
                }
            }
        }
    }
}

// 10. Calculate the M1
template <typename T>
inline void compute_M1(
    int i, int j, int l, int nuci, int nucj, int index,
    double nuci_score, double nucj_score,
    const vector<int> *nuci1_list,
    const vector<int> *nucj_1_list,
    vector<int> &con_seq,
    const std::vector<std::unordered_map<int, State<T>>> &bestM1,
    std::unordered_map<int, State<T>> &bestM1_j,
    const std::unordered_map<int, State<T>> &bestM2_j,
    const std::unordered_map<int, State<T>> &bestSC_j_l1,
    bool &is_M1_OCT,
    bool &is_M1_STEP,
    double &m1_coterminus_score)
{
    // M1(i,j) = M2(i,j)
    if (j - i > MIN_STACKPAIR_GAP * 2 - 2) {
        auto it = bestM2_j.find(index);
        if (it != bestM2_j.end()) {
            int first_pair_pos = it->second.first_pair_pos;
            int last_pair_pos = it->second.last_pair_pos;
            double newscore = it->second.score;
            update(bestM1_j, index, newscore, -1, -1, -1, -1, MANNER_M1_M2, first_pair_pos, last_pair_pos);
        }
    }

    // M1(i,j) = M1(i,j-1) + F(j,j)
    if (j - i > HAIRPIN_GAP + 1) {
        const std::unordered_map<int, State<T>> &bestM1_j_1 = bestM1[j - 1];
        for (auto &nucj_1 : *nucj_1_list) {
            if (IsLegal(con_seq, nucj_1, nucj, j - 1)) {
                int index_j_1 = GetIndex(i, nuci, nucj_1);
                auto it = bestM1_j_1.find(index_j_1);
                if (it != bestM1_j_1.end()) {
                    int first_pair_pos = it->second.first_pair_pos;
                    int last_pair_pos = it->second.last_pair_pos;
                    if (j - last_pair_pos <= SINGLE_MAX_LEN) {
                        double newscore = it->second.score + nucj_score;
                        update(bestM1_j, index, newscore, index_j_1, -1, -1, -1, MANNER_M1_M1_F, first_pair_pos, last_pair_pos);
                    }
                }
            }
        }
    }

    // M1(i,j) = F(i,i) + M1(i+1,j)
    if (j - i > HAIRPIN_GAP + 1) {
        for (auto &nuci1 : *nuci1_list) {
            if (IsLegal(con_seq, nuci, nuci1, i)) {
                int index_j = GetIndex(i + 1, nuci1, nucj);
                auto it = bestM1_j.find(index_j);
                if (it != bestM1_j.end()) {
                    int first_pair_pos = it->second.first_pair_pos;
                    int last_pair_pos = it->second.last_pair_pos;
                    if (first_pair_pos - i <= SINGLE_MAX_LEN) {
                        double newscore = it->second.score + nuci_score;
                        auto it1 = bestM1_j.find(index);
                        if ((it1 == bestM1_j.end()) || (it1->second.score <= newscore)) {
                            update(bestM1_j, index, newscore, index_j, -1, -1, -1, MANNER_M1_F_M1, first_pair_pos, last_pair_pos);
                            is_M1_STEP = false;
                        }
                    }
                }
            }
        }
    }

    // M1(i,j) = C(i,j) + ML_stem
    if ((j - i > HAIRPIN_GAP) && _allowed_pairs[nuci][nucj]) {
        auto it = bestSC_j_l1.find(index);
        if (it != bestSC_j_l1.end()) {
            m1_coterminus_score = it->second.score - MFE_scalar * v_score_M1(EnrBASE(nuci), EnrBASE(nucj));
            if (m1_coterminus_score > bestM1_j[index].score) {
                update(bestM1_j, index, m1_coterminus_score, -1, -1, -1, -1, MANNER_M1_C, i, j);
                is_M1_OCT = true;
                is_M1_STEP = true;
            }
        }
    }

    if (is_M1_STEP) {
        auto it = bestM1_j.find(index);
        if (it != bestM1_j.end()) {
            it->second.is_M1_STEP = true;
        }
    }
}

// 11. Forward Update M2(k,j), 0 <= k <= i-1
template <typename T>
inline void forward_update_M2(
    int i, int j, int nuci, int nucj, int index,
    double m1_coterminus_score,
    bool is_M1_OCT,
    vector<int> &con_seq,
    const std::vector<std::unordered_map<int, State<T>>> &bestM1,
    std::unordered_map<int, State<T>> &bestM2_j)
{
    if (is_M1_OCT && (i > MIN_STACKPAIR_GAP - 1)) {
        // M2(k,j) = M1(k,i-1) + M1(i,j) (STEP + OCT) Forward Update
        const std::unordered_map<int, State<T>> &bestM1_i_1 = bestM1[i - 1];
        for (auto &state_k : bestM1_i_1) {
            if (!state_k.second.is_M1_STEP) {
                continue;
            }
            int k, nuck, nuci_1;
            std::tie(k, nuck, nuci_1) = GetIndexTuple(state_k.first);
            if (!IsLegal(con_seq, nuci_1, nuci, i - 1)) {
                continue;
            }

            int first_pair_pos = state_k.second.first_pair_pos;
            int index_k_j = GetIndex(k, nuck, nucj);
            double newscore = state_k.second.score + m1_coterminus_score;
            update(bestM2_j, index_k_j, newscore, state_k.first, index, -1, i - 1, MANNER_M2_M1_M1, first_pair_pos, j);
        }
    }
}
// 12. Calculate the F (non-speedup mode)
template <typename T>
inline void compute_F(
    int i, int j, int l, int nuci, int nucj, int index,
    double nucj_score,
    const vector<int> *nuci1_list,
    vector<int> &con_seq,
    const std::vector<std::unordered_map<int, State<T>>> &bestF,
    std::unordered_map<int, State<T>> &bestF_j,
    const std::unordered_map<int, State<T>> &bestSC_j_l1)
{
        bool is_F_OCT = false;
        bool is_F_STEP = true;
        double f_coterminus_score = 0.0;
        double newscore;

        // F(i,j) = F(i,i) + F(i+1, j), the non-STEP case
        if (j > i + 1) {
            int index_i_i = GetIndex(i, nuci, nuci);
            for (auto &nuci1 : *nuci1_list) {
                if (IsLegal(con_seq, nuci, nuci1, i)) {
                    int index_i1 = GetIndex(i + 1, nuci1, nucj);

                    auto it = bestF_j.find(index_i1);
                    if (it != bestF_j.end()) {
                        newscore = bestF[i].at(index_i_i).score + it->second.score;

                        auto it1 = bestF_j.find(index);
                        if ((it1 == bestF_j.end()) || (it1->second.score <= newscore)) {
                            bool res = update(bestF_j, index, newscore, index_i_i, index_i1, -1, i, MANNER_F_F_F);
                            is_F_STEP = false;

                            if (res) {
                                bestF_j[index].outer_pairs = it->second.outer_pairs;
                            }
                        }
                    }
                }
            }
        } else if (j == i + 1) {
            if (IsLegal(con_seq, nuci, nucj, i)) {
                int index_i_i = GetIndex(i, nuci, nuci);
                int index_i1 = GetIndex(i + 1, nucj, nucj);
                newscore = bestF[i].at(index_i_i).score + bestF_j[index_i1].score;

                auto it1 = bestF_j.find(index);
                if ((it1 == bestF_j.end()) || (it1->second.score <= newscore)) {
                    update(bestF_j, index, newscore, index_i_i, index_i1, -1, i, MANNER_F_F_F);
                    is_F_STEP = false;
                }
            }
        }

        // F(i,j) = max( F(i,j), C(i,j) )
        if (i == j) {
            f_coterminus_score = nucj_score;
            update(bestF_j, index, f_coterminus_score, -1, -1, -1, -1, NONE);
            is_F_OCT = true;
            is_F_STEP = false;
        } else {
            auto it = bestSC_j_l1.find(index);
            if (it != bestSC_j_l1.end()) {
                f_coterminus_score = it->second.score + (-MFE_scalar * v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj)));

                if (f_coterminus_score > bestF_j[index].score) {
                    bool res = update(bestF_j, index, f_coterminus_score, -1, -1, -1, -1, MANNER_F_C);
                    is_F_OCT = true;
                    is_F_STEP = true;
                    if (res) {
                        bestF_j[index].outer_pairs.push_back(index);
                    }
                }
            }
        }

        if (is_F_STEP) {
            auto it = bestF_j.find(index);
            if (it != bestF_j.end()) {
                it->second.is_F_STEP = true;
            }
        }

        // Forward Update F(k,j), 0 <= k <= i-1
        if (is_F_OCT && (i != 0)) {
            const std::unordered_map<int, State<T>> &bestF_i_1 = bestF[i - 1];
            for (auto &state_k : bestF_i_1) {
                if (!state_k.second.is_F_STEP) {
                    continue;
                }
                int k, nuck, nuci_1;
                std::tie(k, nuck, nuci_1) = GetIndexTuple(state_k.first);

                if (!IsLegal(con_seq, nuci_1, nuci, i - 1)) {
                    continue;
                }

                int index_k_j = GetIndex(k, nuck, nucj);
                newscore = state_k.second.score + f_coterminus_score;

                bool res = update(bestF_j, index_k_j, newscore, state_k.first, index, -1, i - 1, MANNER_F_F_F);
                if (res) {
                    bestF_j[index_k_j].outer_pairs = state_k.second.outer_pairs;
                    if (i != j) {
                        bestF_j[index_k_j].outer_pairs.push_back(index);
                    }
                }
            }
        }
    }

    // 13. Calculate the F (speedup mode: i == 0)
template <typename T>
inline void compute_F_speedup(
    int i, int j,
    const vector<int> &nucj_list,
    vector<int> &con_seq,
    const string &ami_seq,
    bool is_DN,
    const std::vector<std::unordered_map<int, State<T>>> &bestF,
    std::unordered_map<int, State<T>> &bestF_j,
    const std::vector<std::unordered_map<int, State<T>>> &bestSC_j)
{
        if (i == j) { // F(0, 0)
            for (auto &nucj : nucj_list) {
                int index = GetIndex(i, nucj, nucj);
                update(bestF_j, index, 0.0, -1, -1, -1, -1, NONE);
            }
        } else {
            const std::unordered_map<int, State<T>> &bestF_j_1 = bestF[j - 1];

            // F(0, j) = F(0, j-1) + F(j, j)
            for (auto &nucj : nucj_list) {
                double nucj_score = isLastNuc(j) ? getCAI(ami_seq[j / 3], nucj, is_DN) : 0.0;

                for (auto &itemF_j_1 : bestF_j_1) {
                    const State<T> &stateF = itemF_j_1.second;
                    int i_, nuci_, nucj_1;
                    std::tie(i_, nuci_, nucj_1) = GetIndexTuple(itemF_j_1.first);

                    if (IsLegal(con_seq, nucj_1, nucj, j - 1)) {
                        double newscore = stateF.score + nucj_score;
                        int index = GetIndex(i_, nuci_, nucj);
                        int index_j_j = GetIndex(j, nucj, nucj);
                        update(bestF_j, index, newscore, itemF_j_1.first, index_j_j, -1, j - 1, MANNER_F_F_F);
                    }
                }
            }

            // F(0, j) = C(0, j)
            if (j - i > HAIRPIN_GAP) {
                int k = 0;
                for (auto &itemC_k : bestSC_j[j - k + 1]) {
                    const State<T> &stateC = itemC_k.second;
                    int k_, nuck, nucj_;
                    std::tie(k_, nuck, nucj_) = GetIndexTuple(itemC_k.first);

                    int index = itemC_k.first;
                    double newscore = stateC.score + (-MFE_scalar * v_score_external_paired(EnrBASE(nuck), EnrBASE(nucj_)));
                    update(bestF_j, index, newscore, -1, -1, -1, -1, MANNER_F_C);
                }
            }

            // F(0, j) = F(0, k-1) + C(k, j)
            for (int k = 1; j - k > HAIRPIN_GAP; k++) {
                const std::unordered_map<int, State<T>> &bestF_k_1 = bestF[k - 1];

                for (auto &itemC_k : bestSC_j[j - k + 1]) {
                    const State<T> &stateC = itemC_k.second;
                    int k_, nuck, nucj_;
                    std::tie(k_, nuck, nucj_) = GetIndexTuple(itemC_k.first);

                    for (auto &itemF_k_1 : bestF_k_1) {
                        const State<T> &stateF = itemF_k_1.second;
                        int i_, nuci_, nuck_1;
                        std::tie(i_, nuci_, nuck_1) = GetIndexTuple(itemF_k_1.first);

                        if (!IsLegal(con_seq, nuck_1, nuck, k - 1)) {
                            continue;
                        }

                        int index = GetIndex(i_, nuci_, nucj_);
                        double newscore = stateF.score + stateC.score + (-MFE_scalar * v_score_external_paired(EnrBASE(nuck), EnrBASE(nucj_)));

                        update(bestF_j, index, newscore, itemF_k_1.first, itemC_k.first, -1, k - 1, MANNER_F_F_C);
                    }
                }
            }
        }
    }
/** *************************************************************************************************************** */

template <typename T>
void SparseCDSfoldCAI(AllTables<T> &alltables, string &rna_seq, vector<int> &con_seq, string &ami_seq, int ires_len, bool is_DN) {
    int seq_length = rna_seq.size();
    double newscore;
    double nuci_score, nucj_score;
    double f_coterminus_score, m1_coterminus_score;
    int index, index_i_i, index_i1, index_j_1, index_j;
    int index1;
    int index_k_j, k, nuck, nuci_1;

    bool is_F_OCT, is_M1_OCT, is_F_STEP, is_M1_STEP;

    struct timespec start, end;
    struct timespec start1, end1;
    // double internal_time = 0.0;
    double multi_time = 0.0;
    double F_time = 0.0;

    State<T> state_i_1;

    std::vector<std::unordered_map<int, State<T>>> &bestN = alltables.bestN;
    std::vector<std::vector<std::array<State<T>, 257>>> &bestSS = alltables.bestSS;
    std::vector<std::vector<std::array<bool, 257>>> &bestSS_check = alltables.bestSS_check;
    std::vector<std::unordered_map<int, State<T>>> &bestF = alltables.bestF;
    std::vector<std::vector<std::unordered_map<int, State<T>>>> &bestSC = alltables.bestSC;
    std::vector<std::unordered_map<int, State<T>>> &bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>> &bestM2 = alltables.bestM2;

    int size_ = seq_length * (1 + seq_length) / 2;
    for (int j = 0; j < seq_length; j++) {
        std::unordered_map<int, State<T>> &bestN_j = bestN[j];
        std::vector<std::array<State<T>, 257>> &bestSS_j = bestSS[j];
        std::vector<std::array<bool, 257>> &bestSS_check_j = bestSS_check[j];
        std::unordered_map<int, State<T>> &bestF_j = bestF[j];
        std::vector<std::unordered_map<int, State<T>>> &bestSC_j = bestSC[j];
        std::unordered_map<int, State<T>> &bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>> &bestM1_j = bestM1[j];
        vector<int> &nucj_list = Base_table.at(rna_seq[j]);
        const std::vector<int> *nucj_1_list = nullptr;
        if (j > 0) {
            nucj_1_list = &Base_table.at(rna_seq[j - 1]);
        }

        for (int i = j; i >= 0; i--) {
            if (show_processing)
                Processing(size_);

            int l = j - i; // 兩端距離

            vector<int> &nuci_list = Base_table.at(rna_seq[i]);
            const std::vector<int> *nuci1_list = nullptr;
            if (i < seq_length - 1) {
                nuci1_list = &Base_table.at(rna_seq[i + 1]);
            }

            for (auto &nucj : nucj_list) {
                for (auto &nuci : nuci_list) {
                    if ((i == j) && (nuci != nucj)) { continue; }

                    index = GetIndex(i, nuci, nucj);
                    index1 = GetIndex1(nuci, nucj);

                    is_F_OCT = false;
                    is_M1_OCT = false;
                    is_F_STEP = true;
                    is_M1_STEP = true;

                    nuci_score = isLastNuc(i) ? getCAI(ami_seq[i / 3], nuci, is_DN) : 0.0;
                    nucj_score = isLastNuc(j) ? getCAI(ami_seq[j / 3], nucj, is_DN) : 0.0;

                    // 1. calculate the N, which is unpair
                    compute_N(i, j, nuci, nucj, index, nuci_score, nucj_score, nucj_1_list, con_seq, bestN_j, bestN);

                    // 2. calculate the S, which is unpair
                    compute_S(i, j, l, nuci, nucj, index1, nuci_score, nucj_score, nucj_1_list, con_seq, bestSS_j, bestSS_check_j, bestSS, bestSS_check);

                    // 3. hairpin
                    compute_hairpin(i, j, l, nuci, nucj, index, nuci_score, nucj_score, nuci1_list, nucj_1_list, con_seq, bestN, bestSC_j[l + 1]);

                    clock_gettime(CLOCK_MONOTONIC, &start);

                    // 4. stacking
                    compute_stacking(i, j, l, nuci, nucj, index, nuci_score, nucj_score, rna_seq, con_seq, bestSC, bestSC_j[l + 1]);

                    // 5. left bulge
                    compute_left_bulge(i, j, l, nuci, nucj, index, nuci_score, nucj_score, nucj_1_list, rna_seq, con_seq, bestSS, bestSS_check, bestSC, bestSC_j[l + 1]);

                    // 6. right bulge
                    compute_right_bulge(i, j, l, nuci, nucj, index, nuci_score, nucj_score, nuci1_list, rna_seq, con_seq, bestSS, bestSS_check, bestSC, bestSC_j[l + 1]);

                    // 7. internal
                    compute_internal(i, j, l, nuci, nucj, index, nuci_score, nucj_score, nuci1_list, nucj_1_list, rna_seq, con_seq, bestSS, bestSS_check, bestSC, bestSC_j[l + 1]);

                    clock_gettime(CLOCK_MONOTONIC, &end);
                    internal_time += calTime(start, end);

                    // 8. multi-loop closing
                    compute_multiloop_closed(i, j, l, nuci, nucj, index, nuci_score, nucj_score, nuci1_list, nucj_1_list, con_seq, bestM2, bestSC_j[l + 1]);

                    // 9. calculate the M2 which i is unpair and j is unpair
                    compute_M2_unpaired(i, j, nuci, nucj, index, nuci_score, nucj_score, nuci1_list, nucj_1_list, con_seq, bestM2, bestM2_j);

                    // 10. calculate the M1
                    compute_M1(i, j, l, nuci, nucj, index, nuci_score, nucj_score, nuci1_list, nucj_1_list, con_seq, bestM1, bestM1_j, bestM2_j, bestSC_j[l + 1], is_M1_OCT, is_M1_STEP, m1_coterminus_score);

                    // 11. Forward Update M2(k,j), 0 <= k <= i-1
                    forward_update_M2(i, j, nuci, nucj, index, m1_coterminus_score, is_M1_OCT, con_seq, bestM1, bestM2_j);

                    if (!F_speedup) {
                        compute_F(i, j, l, nuci, nucj, index, nucj_score, nuci1_list, con_seq, bestF, bestF_j, bestSC_j[l + 1]);
                    }
                } // nuci end
            } // nucj end

            if (i == 0 && F_speedup) {
                compute_F_speedup(i, j, nucj_list, con_seq, ami_seq, is_DN, bestF, bestF_j, bestSC_j);
            }

        } // i end
    } // j end

    return;
}
