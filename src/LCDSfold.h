#pragma once

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
#include <queue>

#include "functions.h"
#include "utility.h"
#include "beamprune.h"
#include "ac_matcher.h"
#include "common.h"

using namespace std;
using namespace oligodesign;
#define GET_ACGU_NUM_V(x) ((x=='A'? 1 : (x=='C'? 2 : (x=='G'? 3 : (x=='U'?4: 0)))))
bool debug = false;

// #define debug

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, Manner MANNER, int cross_pair_num = 0) {
    auto it = stateMap.find(index);

    if (it == stateMap.end()) {
        stateMap.emplace(index, State<T>(newscore, preindex_1, -1, -1, -1, MANNER, -1, -1, cross_pair_num));
    } else if (newscore > it->second.score) {
        it->second = State<T>(newscore, preindex_1, -1, -1, -1, MANNER, -1, -1, cross_pair_num);
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, int preindex_2, Manner MANNER, int cross_pair_num = 0) {
    auto it = stateMap.find(index);

    if (it == stateMap.end()) {
        stateMap.emplace(index, State<T>(newscore, preindex_1, preindex_2, -1, -1, MANNER, -1, -1, cross_pair_num));
    } else if (newscore > it->second.score) {
        it->second = State<T>(newscore, preindex_1, preindex_2, -1, -1, MANNER, -1, -1, cross_pair_num);
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, Manner MANNER, int last_pair_pos, int cross_pair_num) {
    auto it = stateMap.find(index);

    if (it == stateMap.end()) {
        stateMap.emplace(index, State<T>(newscore, preindex_1, -1, -1, -1, MANNER, -1, last_pair_pos, cross_pair_num));
    } else if (newscore > it->second.score) {
        it->second = State<T>(newscore, preindex_1, -1, -1, -1, MANNER, -1, last_pair_pos, cross_pair_num);
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, int preindex_2, Manner MANNER, int last_pair_pos, int cross_pair_num) {
    auto it = stateMap.find(index);

    if (it == stateMap.end()) {
        stateMap.emplace(index, State<T>(newscore, preindex_1, preindex_2, -1, -1, MANNER, -1, last_pair_pos, cross_pair_num));
    } else if (newscore > it->second.score) {
        it->second = State<T>(newscore, preindex_1, preindex_2, -1, -1, MANNER, -1, last_pair_pos, cross_pair_num);
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, int preindex_2, int preindex_3, Manner MANNER, int cross_pair_num = 0) {
    auto it = stateMap.find(index);

    if (it == stateMap.end()) {
        stateMap.emplace(index, State<T>(newscore, preindex_1, preindex_2, preindex_3, -1, MANNER, -1, -1, cross_pair_num));
    } else if (newscore > it->second.score) {
        it->second = State<T>(newscore, preindex_1, preindex_2, preindex_3, -1, MANNER, -1, -1, cross_pair_num);
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, double newscore, int preindex_1, int preindex_2, int preindex_3, Manner MANNER, int last_pair_pos, int cross_pair_num) {
    auto it = stateMap.find(index);

    if (it == stateMap.end()) {
        stateMap.emplace(index, State<T>(newscore, preindex_1, preindex_2, preindex_3, -1, MANNER, -1, last_pair_pos, cross_pair_num));
    } else if (newscore > it->second.score) {
        it->second = State<T>(newscore, preindex_1, preindex_2, preindex_3, -1, MANNER, -1, last_pair_pos, cross_pair_num);
    }
}

/** *************************************************************************************************************** */

/* 
 *  @brief 
 *      None -> S/N state
 *       N   ->  N  state
 *       N   ->  C  state
 */
template<typename T>
void init_and_N_state(std::vector<std::unordered_map<int, State<T>>>& bestN,
                      std::vector<std::unordered_map<int, State<T>>>& bestS_j,
                      std::unordered_map<int, State<T>>& bestC_j,
                      string& rna_seq,
                      vector<int>& con_seq,
                      string& ami_seq,
                      const vector<int>& nucj_list,
                      int j,
                      int ires_len,
                      bool is_DN,
                      double param) {
    int index;
    T newscore;
    int j_1 = j - 1;
    std::unordered_map<int, State<T>>& bestN_j = bestN[j];
    std::unordered_map<int, State<T>>& bestN_j_1 = bestN[j_1];

    for (auto &nucj : nucj_list) { // all possible nuc at position j
        // NONE->S
        // NONE->N
        index = GetIndex(j, nucj, nucj);
        newscore = 0;

        if (isLastNuc(j) && param)
            newscore += getCAI(ami_seq[j / 3], nucj, is_DN);

        update(bestS_j[1], index, newscore, -1, MANNER_NONEtoS, 0); // initialize bestS[j]
        update(bestN_j, index, newscore, -1, MANNER_NONEtoN, 0); // initialize bestN[j]

        /*all N -> another state*/
        for (auto &itemN_j_1 : bestN_j_1) { // all different states -> N at position j-1
            int index_j_1 = itemN_j_1.first; // key
            State<T> stateN_j_1 = itemN_j_1.second; // state
            int i, nucj_1, nuci;
            std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1); // from which position i and nuc of i and j-1
            if (IsLegal(con_seq, nucj_1, nucj, j_1)) {
                // because N -> N means Non-pairing to Non-pairing, only need to consider CAI(at position j)
                newscore = stateN_j_1.score;

                if (isLastNuc(j) && param)
                    newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                // N -> N
                index = GetIndex(i, nuci, nucj);
                update(bestN_j, index, newscore, index_j_1, MANNER_N_EtoN, stateN_j_1.cross_pair_num);

                int i_1 = i - 1; // the continuous Non-pairing is from i to j-1, so if we want to transfer N -> C ,the pairing must be i-1 and j (loop : i ~ j-1)
                if ((i_1 >= 0) && (j - i_1 > HAIRPIN_GAP)) { // consider hairpin(gap > 3)
                    vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); // i-1 list
                    for (auto &nuci_1 : nuci_1_list) {
                        if (IsLegal(con_seq, nuci_1, nuci, i_1)) {
                            if (_allowed_pairs[nuci_1][nucj]) {
                                // because N -> C means Non-pairing to Closing, need to consider hairpin structure
                                if (is_DN) {
                                    newscore = -lambda * v_score_hairpin(i_1, j, EnrBASE(nuci_1), EnrBASE(nuci), EnrBASE(nucj_1), EnrBASE(nucj)) + stateN_j_1.score;
                                } else {
                                    newscore = -v_score_hairpin(i_1, j, EnrBASE(nuci_1), EnrBASE(nuci), EnrBASE(nucj_1), EnrBASE(nucj)) + stateN_j_1.score;
                                }
                                if (isLastNuc(j) && param)
                                    newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                                if (isLastNuc(i_1) && param)
                                    newscore += getCAI(ami_seq[i_1 / 3], nuci_1, is_DN);

                                int cross_delta = isCrossTalk(i_1, j, ires_len);
                                int total_cross = stateN_j_1.cross_pair_num + cross_delta;
                                newscore -= cross_delta * crosspair_penalty;
                                
                                // N->C
                                index = GetIndex(i_1, nuci_1, nucj);
                                update(bestC_j, index, newscore, index_j_1, MANNER_NtoC, total_cross);
                            }
                        }
                    }
                }
            }
        }
    }
}

/*
 * @brief 
 *   S -> S state
 */
template<typename T>
void S_to_S_state(std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS,
                  std::vector<std::unordered_map<int, State<T>>>& bestS_j,
                  vector<int>& con_seq,
                  string& ami_seq,
                  const vector<int>& nucj_list,
                  int j,
                  bool is_DN,
                  double param) {
    int index;
    T newscore;
    int j_1 = j - 1;

    for (auto &nucj : nucj_list) {
        for (int l = 1; l <= min(SINGLE_MAX_LEN - 1, j); l++) {
            std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1][l];
            for (auto &itemS_j_1 : bestS_j_1) {
                int index_j_1 = itemS_j_1.first;
                State<T> stateS_j_1 = itemS_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if (IsLegal(con_seq, nucj_1, nucj, j_1) && (j - i) < SINGLE_MAX_LEN) {
                    // S means single length extend
                    newscore = stateS_j_1.score;
                    if (isLastNuc(j) && param)
                        newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                    // S -> S
                    index = GetIndex(i, nuci, nucj);
                    update(bestS_j[l + 1], index, newscore, index_j_1, MANNER_S_EtoS, stateS_j_1.cross_pair_num);
                }
            }
        }
    }
}

/*
 * @brief 
 *   F -> F state
 *   C -> F state
 *   F + C -> F state
 */
template<typename T>
void C_and_F_state(std::vector<std::unordered_map<int, State<T>>>& bestF,
                   std::unordered_map<int, State<T>>& bestF_j,
                   std::unordered_map<int, State<T>>& bestC_j,
                   vector<int>& con_seq,
                   string& ami_seq,
                   const vector<int>& nucj_list,
                   int j,
                   int seq_length,
                   bool is_DN,
                   double param) {
    int index;
    T newscore;
    int j_1 = j - 1;

    for (auto &nucj : nucj_list) {
        std::unordered_map<int, State<T>>& bestF_j_1 = bestF[j_1];
        for (auto &itemF_j_1 : bestF_j_1) {
            int index_j_1 = itemF_j_1.first;
            State<T> stateF_j_1 = itemF_j_1.second;
            int i, nucj_1, nuci;
            std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
            if (IsLegal(con_seq, nucj_1, nucj, j_1)) {
                newscore = stateF_j_1.score;
                if (isLastNuc(j) && param)
                    newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                // F->F
                index = GetIndex(i, nuci, nucj);
                update(bestF_j, index, newscore, index_j_1, MANNER_F_EtoF, stateF_j_1.cross_pair_num);
            }
        }
    }

    for (auto &itemC_j : bestC_j) {
        int index_j = itemC_j.first;
        State<T> stateC_j = itemC_j.second;
        int i, nucj, nuci;
        std::tie(i, nuci, nucj) = GetIndexTuple(index_j);

        if (i == 0) {
            if (is_DN) {
                newscore = stateC_j.score - lambda * v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
            } else {
                newscore = stateC_j.score - v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
            }
            // C->F
            update(bestF_j, index_j, newscore, index_j, MANNER_CtoF, stateC_j.cross_pair_num);
        }

        int i_1 = i - 1;
        if (i_1 >= 0) {
            std::unordered_map<int, State<T>>& bestF_i_1 = bestF[i_1];
            for (auto &itemF_i_1 : bestF_i_1) {
                int index_i_1 = itemF_i_1.first;
                State<T> stateF_i_1 = itemF_i_1.second;
                int k, nuck, nuci_1;
                std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                if ((k == 0 || j == seq_length - 1) && IsLegal(con_seq, nuci_1, nuci, i_1)) {
                    if (is_DN) {
                        newscore = stateF_i_1.score + stateC_j.score - lambda * v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                    } else {
                        newscore = stateF_i_1.score + stateC_j.score - v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                    }
                    // F+C->F
                    int total_cross = stateF_i_1.cross_pair_num + stateC_j.cross_pair_num;
                    index = GetIndex(k, nuck, nucj);
                    update(bestF_j, index, newscore, index_i_1, index_j, MANNER_F_CtoF, total_cross);
                }
            }
        }
    }
}

/*
 * @brief 
 *   C + S -> CS state
 */
template<typename T>
void C_to_CS_state(std::vector<std::unordered_map<int, State<T>>>& bestS_j,
                   std::vector<std::unordered_map<int, State<T>>>& bestC,
                   std::unordered_map<int, State<T>>& bestCS_j,
                   vector<int>& con_seq,
                   int j) {
    int index;
    T newscore;

    for (int l = min(SINGLE_MAX_LEN, j - 4); l >= 1; --l) {
        for (auto &itemS_j : bestS_j[l]) {
            int index_j = itemS_j.first;
            State<T> stateS_j = itemS_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            int i_1 = i - 1;

            if (i_1 >= 0) {
                std::unordered_map<int, State<T>>& bestC_i_1 = bestC[i_1];
                for (auto &itemC_i_1 : bestC_i_1) {
                    int index_i_1 = itemC_i_1.first;
                    State<T> stateC_i_1 = itemC_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);

                    if (k >= 0 && IsLegal(con_seq, nuci_1, nuci, i_1)) {
                        newscore = stateS_j.score + stateC_i_1.score;
                        int total_cross = stateS_j.cross_pair_num + stateC_i_1.cross_pair_num;

                        // C+S->CS
                        index = GetIndexCS(k, nuck, nuci_1, nucj, nuci, l);
                        update(bestCS_j, index, newscore, index_i_1, index_j, MANNER_C_StoCS, total_cross);
                    }
                }
            }
        }
    }
}

/*
 * @brief 
 *   CS  -> C state
 *  S+CS -> C state 
 */
template<typename T>
void CS_to_C_state(std::unordered_map<int, State<T>>& bestCS_j_1,
                   std::unordered_map<int, State<T>>& bestC_j,
                   std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS,
                   string& rna_seq,
                   vector<int>& con_seq,
                   string& ami_seq,
                   int nucj,
                   int j,
                   int ires_len,
                   bool is_DN,
                   double param) {
    int index;
    T newscore;
    int j_1 = j - 1;

    for (auto &itemCS_j_1 : bestCS_j_1) {
        int index_j_1 = itemCS_j_1.first;
        State<T> stateCS_j_1 = itemCS_j_1.second;
        int i1, nucj_1, nuci1, nuci1_pair, nucx, len;
        std::tie(i1, nuci1, nuci1_pair, nucj_1, nucx, len) = GetIndexTupleCS(index_j_1);
        int i = i1 - 1;
        if (i >= 0 && IsLegal(con_seq, nucj_1, nucj, j_1)) {
            vector<int> nuci_list = Base_table.at(rna_seq[i]);
            for (auto &nuci : nuci_list) {
                if (IsLegal(con_seq, nuci, nuci1, i) && _allowed_pairs[nuci][nucj]) {
                    int pre_index_2 = stateCS_j_1.index_2;
                    int pre_index_1 = stateCS_j_1.index_1;
                    int s, nucs, nucl;
                    std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                    int s_1 = s - 1;
                    int a, nuca, nucs_1;
                    std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);
                    if (is_DN) {
                        newscore = -lambda * v_score_single(i, j, i1, s_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score;
                    } else {
                        newscore = -v_score_single(i, j, i1, s_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score;
                    }

                    if (isLastNuc(j) && param)
                        newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                    if (isLastNuc(i) && param)
                        newscore += getCAI(ami_seq[i / 3], nuci, is_DN);

                    int cross_delta = isCrossTalk(i, j, ires_len);
                    int total_cross = stateCS_j_1.cross_pair_num + cross_delta;
                    newscore -= cross_delta * crosspair_penalty;

                    // CS->C
                    index = GetIndex(i, nuci, nucj);
                    update(bestC_j, index, newscore, index_j_1, MANNER_CStoC, total_cross);
                }
            }

            for (int l = 1; l <= min(SINGLE_MAX_LEN - len, i); l++) {
                std::unordered_map<int, State<T>>& bestS_i = bestS[i][l];
                for (auto &itemS_i : bestS_i) {
                    int index_i = itemS_i.first;
                    State<T> stateS_i = itemS_i.second;
                    int k1, nuck1, nuci;
                    std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                    int k = k1 - 1;

                    if (k >= 0 && IsLegal(con_seq, nuci, nuci1, i) && i - k + len <= SINGLE_MAX_LEN) {
                        vector<int> nuck_list = Base_table.at(rna_seq[k]);
                        for (auto &nuck : nuck_list) {
                            if (IsLegal(con_seq, nuck, nuck1, k) && _allowed_pairs[nuck][nucj]) {
                                int pre_index_2 = stateCS_j_1.index_2;
                                int pre_index_1 = stateCS_j_1.index_1;
                                int s, nucs, nucl;
                                std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                                int s_1 = s - 1;
                                int a, nuca, nucs_1;
                                std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);
                                if (is_DN) {
                                    newscore = -lambda * v_score_single(k, j, i1, s_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                        EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score + stateS_i.score;
                                } else {
                                    newscore = -v_score_single(k, j, i1, s_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                        EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score + stateS_i.score;
                                }

                                if (isLastNuc(j) && param)
                                    newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                                if (isLastNuc(k) && param)
                                    newscore += getCAI(ami_seq[k / 3], nuck, is_DN);

                                int cross_delta = isCrossTalk(k, j, ires_len);
                                int total_cross = stateCS_j_1.cross_pair_num + stateS_i.cross_pair_num + cross_delta;
                                newscore -= cross_delta * crosspair_penalty;

                                // S+CS->C
                                index = GetIndex(k, nuck, nucj);
                                update(bestC_j, index, newscore, index_i, index_j_1, MANNER_S_CStoC, total_cross);
                            }
                        }
                    }
                }
            }
        }
    }
}

/*
 * @brief 
 *   C + S -> C state
 *   S + C + S -> C state
 */
template<typename T>
void S_C_S_state(std::vector<std::unordered_map<int, State<T>>>& bestC,
                 std::unordered_map<int, State<T>>& bestC_j,
                 std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS,
                 std::vector<std::unordered_map<int, State<T>>>& bestF,
                 string& rna_seq,
                 vector<int>& con_seq,
                 string& ami_seq,
                 const vector<int>& nucj_list,
                 int j,
                 int ires_len,
                 bool is_DN,
                 double param,
                 bool mid_beam_prune = false) {
    int index;
    T newscore;
    int j_1 = j - 1;

    for (auto &nucj : nucj_list) {
        for (int l1 = min(SINGLE_MAX_LEN, j - 4); l1 >= 1; --l1) {
            std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1][l1];
            for (auto &itemS_j_1 : bestS_j_1) {
                int index_j_1 = itemS_j_1.first;
                State<T> stateS_j_1 = itemS_j_1.second;
                int q1, nucq1, nucj_1;
                std::tie(q1, nucq1, nucj_1) = GetIndexTuple(index_j_1);
                int q = q1 - 1;

                if (q >= 0 && IsLegal(con_seq, nucj_1, nucj, j_1)) {
                    std::unordered_map<int, State<T>>& bestC_q = bestC[q];
                    for (auto &itemC_q : bestC_q) {
                        int index_q = itemC_q.first;
                        State<T> stateC_q = itemC_q.second;
                        int p, nucp, nucq;
                        std::tie(p, nucp, nucq) = GetIndexTuple(index_q);
                        int p_1 = p - 1;

                        if (p_1 >= 0 && IsLegal(con_seq, nucq, nucq1, q)) {
                            // C_StoC (right bulge)
                            int i = p_1;
                            vector<int> nuci_list = Base_table.at(rna_seq[i]);
                            for (auto &nuci : nuci_list) {
                                if (IsLegal(con_seq, nuci, nucp, i) && _allowed_pairs[nuci][nucj]) {
                                    if (is_DN) {
                                        newscore = -lambda * v_score_single(i, j, p, q, EnrBASE(nuci), EnrBASE(nucp), EnrBASE(nucj_1), EnrBASE(nucj), EnrBASE(nuci), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucq1))
                                            + stateS_j_1.score + stateC_q.score;
                                    } else {
                                        newscore = -v_score_single(i, j, p, q, EnrBASE(nuci), EnrBASE(nucp), EnrBASE(nucj_1), EnrBASE(nucj), EnrBASE(nuci), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucq1))
                                            + stateS_j_1.score + stateC_q.score;
                                    }

                                    if (isLastNuc(j) && param)
                                        newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                                    if (isLastNuc(i) && param)
                                        newscore += getCAI(ami_seq[i / 3], nuci, is_DN);

                                    int cross_delta = isCrossTalk(i, j, ires_len);
                                    int total_cross = stateS_j_1.cross_pair_num + stateC_q.cross_pair_num + cross_delta;
                                    newscore -= cross_delta * crosspair_penalty;

                                    index = GetIndex(i, nuci, nucj);
                                    update(bestC_j, index, newscore, index_q, index_j_1, MANNER_C_StoC, total_cross);
                                }
                            }

                            if (mid_beam_prune && is_DN) {
                                BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
                            }

                            // S_C_StoC (internal loop)
                            for (int l2 = 1; l2 <= min(SINGLE_MAX_LEN - l1, p_1); ++l2) {
                                std::unordered_map<int, State<T>>& bestS_p_1 = bestS[p_1][l2];
                                for (auto &itemS_p_1 : bestS_p_1) {
                                    int index_p_1 = itemS_p_1.first;
                                    State<T> stateS_p_1 = itemS_p_1.second;
                                    int i1, nuci1, nucp_1;
                                    std::tie(i1, nuci1, nucp_1) = GetIndexTuple(index_p_1);
                                    int i = i1 - 1;
                                    if (i >= 0 && IsLegal(con_seq, nucp_1, nucp, p_1)) {
                                        vector<int> nuci_list = Base_table.at(rna_seq[i]);
                                        for (auto &nuci : nuci_list) {
                                            if (IsLegal(con_seq, nuci, nuci1, i) && _allowed_pairs[nuci][nucj]) {
                                                if (is_DN) {
                                                    newscore = -lambda * v_score_single(i, j, p, q, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj), EnrBASE(nucp_1), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucq1))
                                                        + stateS_j_1.score + stateC_q.score + stateS_p_1.score;
                                                } else {
                                                    newscore = -v_score_single(i, j, p, q, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj), EnrBASE(nucp_1), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucq1))
                                                        + stateS_j_1.score + stateC_q.score + stateS_p_1.score;
                                                }

                                                if (isLastNuc(j) && param)
                                                    newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                                                if (isLastNuc(i) && param)
                                                    newscore += getCAI(ami_seq[i / 3], nuci, is_DN);

                                                int cross_delta = isCrossTalk(i, j, ires_len);
                                                int total_cross = stateS_j_1.cross_pair_num + stateC_q.cross_pair_num + stateS_p_1.cross_pair_num + cross_delta;
                                                newscore -= cross_delta * crosspair_penalty;

                                                index = GetIndex(i, nuci, nucj);
                                                update(bestC_j, index, newscore, index_p_1, index_q, index_j_1, MANNER_S_C_StoC, total_cross);
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
    }
}

/*
 * @brief 
 *   C -> C state
 *   S + C -> C state
 */
template<typename T>
void C_state(std::vector<std::unordered_map<int, State<T>>>& bestC,
             std::unordered_map<int, State<T>>& bestC_j,
             std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS,
             string& rna_seq,
             vector<int>& con_seq,
             string& ami_seq,
             const vector<int>& nucj_list,
             int j,
             int ires_len,
             bool is_DN,
             double param) {
    int index;
    T newscore;
    int j_1 = j - 1;

    for (auto &nucj : nucj_list) {
        std::unordered_map<int, State<T>>& bestC_j_1 = bestC[j_1];
        for (auto &itemC_j_1 : bestC_j_1) {
            int index_j_1 = itemC_j_1.first;
            State<T> stateC_j_1 = itemC_j_1.second;
            int i1, nucj_1, nuci1;
            std::tie(i1, nuci1, nucj_1) = GetIndexTuple(index_j_1);
            int i = i1 - 1;

            if (IsLegal(con_seq, nucj_1, nucj, j_1) && (i >= 0) && (j - i > (HAIRPIN_GAP + 2))) {
                vector<int> nuci_list = Base_table.at(rna_seq[i]);
                for (auto &nuci : nuci_list) {
                    if (_allowed_pairs[nuci][nucj] && IsLegal(con_seq, nuci, nuci1, i) && IsLegal(con_seq, nucj_1, nucj, j_1)) {
                        if (is_DN) {
                            newscore = -lambda * v_score_single(i, j, i1, j_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                                EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score;
                        } else {
                            newscore = -v_score_single(i, j, i1, j_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                                EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score;
                        }

                        if (isLastNuc(j) && param)
                            newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                        if (isLastNuc(i) && param)
                            newscore += getCAI(ami_seq[i / 3], nuci, is_DN);

                        int cross_delta = isCrossTalk(i, j, ires_len);
                        int total_cross = stateC_j_1.cross_pair_num + cross_delta;
                        newscore -= cross_delta * crosspair_penalty;

                        // C->C
                        index = GetIndex(i, nuci, nucj);
                        update(bestC_j, index, newscore, index_j_1, MANNER_CtoC, total_cross);
                    }
                }

                int i_1 = i - 1;
                if (i_1 >= 0) {
                    for (int l = 1; l <= min(SINGLE_MAX_LEN, i); l++) {
                        std::unordered_map<int, State<T>>& bestS_i = bestS[i][l];
                        for (auto &itemS_i : bestS_i) {
                            int index_i = itemS_i.first;
                            State<T> stateS_i = itemS_i.second;
                            int k1, nuck1, nuci;
                            std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                            int k = k1 - 1;
                            if (k >= 0 && IsLegal(con_seq, nuci, nuci1, i)) {
                                vector<int> nuck_list = Base_table.at(rna_seq[k]);
                                for (auto &nuck : nuck_list) {
                                    if (_allowed_pairs[nuck][nucj] && IsLegal(con_seq, nuck, nuck1, k)) {
                                        if (is_DN) {
                                            newscore = -lambda * v_score_single(k, j, i1, j_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                                EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score + stateS_i.score;
                                        } else {
                                            newscore = -v_score_single(k, j, i1, j_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                                EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score + stateS_i.score;
                                        }

                                        if (isLastNuc(j) && param)
                                            newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                                        if (isLastNuc(k) && param)
                                            newscore += getCAI(ami_seq[k / 3], nuck, is_DN);

                                        int cross_delta = isCrossTalk(k, j, ires_len);
                                        int total_cross = stateC_j_1.cross_pair_num + stateS_i.cross_pair_num + cross_delta;
                                        newscore -= cross_delta * crosspair_penalty;

                                        // S+C->C
                                        index = GetIndex(k, nuck, nucj);
                                        update(bestC_j, index, newscore, index_i, index_j_1, MANNER_S_CtoC, total_cross);
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

/*
 * @brief 
 *   Multi -> Multi state
 *   Multi -> C state
 */
template<typename T>
void Multi_state(std::vector<std::unordered_map<int, State<T>>>& bestMulti,
                 std::unordered_map<int, State<T>>& bestMulti_j,
                 std::unordered_map<int, State<T>>& bestC_j,
                 string& rna_seq,
                 vector<int>& con_seq,
                 string& ami_seq,
                 const vector<int>& nucj_list,
                 int j,
                 int ires_len,
                 bool is_DN,
                 double param) {
    int index;
    T newscore;
    int j_1 = j - 1;

    for (auto &nucj : nucj_list) {
        std::unordered_map<int, State<T>>& bestMulti_j_1 = bestMulti[j_1];
        for (auto &itemMulti_j_1 : bestMulti_j_1) {
            int index_j_1 = itemMulti_j_1.first;
            State<T> stateMulti_j_1 = itemMulti_j_1.second;
            int i, nucj_1, nuci;
            std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
            if (IsLegal(con_seq, nucj_1, nucj, j_1)) {
                int last_pair_pos_Multi = stateMulti_j_1.last_pair_pos;
                if (j - last_pair_pos_Multi <= SINGLE_MAX_LEN) {
                    newscore = stateMulti_j_1.score;

                    if (isLastNuc(j) && param)
                        newscore += getCAI(ami_seq[j / 3], nucj, is_DN);

                    // Multi->Multi
                    index = GetIndex(i, nuci, nucj);
                    update(bestMulti_j, index, newscore, index_j_1, -1, stateMulti_j_1.index_3, MANNER_Multi_EtoMulti, last_pair_pos_Multi, stateMulti_j_1.cross_pair_num);
                }

                int i_1 = i - 1;
                if (i_1 >= 0) {
                    vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]);
                    for (auto &nuci_1 : nuci_1_list) {
                        if (IsLegal(con_seq, nuci_1, nuci, i_1) && _allowed_pairs[nuci_1][nucj]) {
                            if (is_DN) {
                                newscore = -lambda * v_score_multi(EnrBASE(nuci_1), EnrBASE(nucj)) + stateMulti_j_1.score;
                            } else {
                                newscore = -v_score_multi(EnrBASE(nuci_1), EnrBASE(nucj)) + stateMulti_j_1.score;
                            }

                            if (isLastNuc(j) && param)
                                newscore += getCAI(ami_seq[j / 3], nucj, is_DN);
                            if (isLastNuc(i_1) && param)
                                newscore += getCAI(ami_seq[i_1 / 3], nuci_1, is_DN);

                            int cross_delta = isCrossTalk(i_1, j, ires_len);
                            int total_cross = stateMulti_j_1.cross_pair_num + cross_delta;
                            newscore -= cross_delta * crosspair_penalty;

                            // Multi->C
                            index = GetIndex(i_1, nuci_1, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_MultitoC, total_cross);
                        }
                    }
                }
            }
        }
    }
}

/*
 * @brief 
 *   C -> M1 state
 *   M1 + C -> M2 state
 */
template<typename T>
void M1_C_to_M2_state(std::vector<std::unordered_map<int, State<T>>>& bestM1,
                      std::unordered_map<int, State<T>>& bestM1_j,
                      std::unordered_map<int, State<T>>& bestM2_j,
                      std::unordered_map<int, State<T>>& bestC_j,
                      vector<int>& con_seq,
                      int j,
                      bool is_DN) {
    int index;
    T newscore;

    for (auto &itemC_j : bestC_j) {
        int index_j = itemC_j.first;
        State<T> stateC_j = itemC_j.second;
        int i, nucj, nuci;
        std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
        if (is_DN) {
            newscore = -lambda * v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score;
        } else {
            newscore = -v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score;
        }

        // C->M1
        update(bestM1_j, index_j, newscore, index_j, -1, 1, MANNER_CtoM1, j, stateC_j.cross_pair_num);

        int i_1 = i - 1;
        if (i_1 >= HAIRPIN_GAP + 1) {
            std::unordered_map<int, State<T>>& bestM1_i_1 = bestM1[i_1];
            for (auto &itemM1_i_1 : bestM1_i_1) {
                int index_i_1 = itemM1_i_1.first;
                State<T> stateM1_i_1 = itemM1_i_1.second;
                int k, nuck, nuci_1;
                std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                if (IsLegal(con_seq, nuci_1, nuci, i_1)) {
                    if (is_DN) {
                        newscore = -lambda * v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score + stateM1_i_1.score;
                    } else {
                        newscore = -v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score + stateM1_i_1.score;
                    }
                    int total_cross = stateC_j.cross_pair_num + stateM1_i_1.cross_pair_num;

                    // M1+C->M2
                    index = GetIndex(k, nuck, nucj);
                    update(bestM2_j, index, newscore, index_i_1, index_j, stateM1_i_1.index_3 + 1, MANNER_M1_CtoM2, j, total_cross);
                }
            }
        }
    }
}

/*
 * @brief 
 *   M1 -> M1 state
 */
template<typename T>
void M1_to_M1_state(std::vector<std::unordered_map<int, State<T>>>& bestM1,
                    std::unordered_map<int, State<T>>& bestM1_j,
                    vector<int>& con_seq,
                    string& ami_seq,
                    const vector<int>& nucj_list,
                    int j,
                    bool is_DN,
                    double param) {
    int index;
    T newscore;
    int j_1 = j - 1;

    for (auto &nucj : nucj_list) {
        std::unordered_map<int, State<T>>& bestM1_j_1 = bestM1[j_1];
        for (auto &itemM1_j_1 : bestM1_j_1) {
            int index_j_1 = itemM1_j_1.first;
            State<T> stateM1_j_1 = itemM1_j_1.second;
            int i, nucj_1, nuci;
            std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);

            int last_pair_pos_M1 = stateM1_j_1.last_pair_pos;
            assert(last_pair_pos_M1 != -1);
            if (j - last_pair_pos_M1 <= SINGLE_MAX_LEN && IsLegal(con_seq, nucj_1, nucj, j_1)) {
                newscore = stateM1_j_1.score;

                if (isLastNuc(j) && param)
                    newscore += getCAI(ami_seq[j / 3], nucj, is_DN);

                // M1->M1
                index = GetIndex(i, nuci, nucj);
                update(bestM1_j, index, newscore, index_j_1, -1, stateM1_j_1.index_3, MANNER_M1_EtoM1, last_pair_pos_M1, stateM1_j_1.cross_pair_num);
            }
        }
    }
}

/*
 * @brief 
 *   M2 -> M1 state
 *   M2 -> Multi state
 *   S + M2 -> Multi state
 */
template<typename T>
void M2_state(std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS,
              std::unordered_map<int, State<T>>& bestM1_j,
              std::unordered_map<int, State<T>>& bestMulti_j,
              std::unordered_map<int, State<T>>& bestM2_j,
              vector<int>& con_seq) {
    int index;
    T newscore;

    for (auto &itemM2_j : bestM2_j) {
        int index_j = itemM2_j.first;
        State<T> stateM2_j = itemM2_j.second;
        int i, nucj, nuci;
        std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
        newscore = stateM2_j.score;
        int last_pair_pos_M2 = stateM2_j.last_pair_pos;

        // M2->M1
        update(bestM1_j, index_j, newscore, index_j, -1, stateM2_j.index_3, MANNER_M2toM1, last_pair_pos_M2, stateM2_j.cross_pair_num);
        // M2->Multi
        update(bestMulti_j, index_j, newscore, index_j, -1, stateM2_j.index_3, MANNER_M2toMulti, last_pair_pos_M2, stateM2_j.cross_pair_num);

        int i_1 = i - 1;
        if (i_1 >= 0) {
            for (int l = 1; l <= min(SINGLE_MAX_LEN, i_1); l++) {
                std::unordered_map<int, State<T>>& bestS_i_1 = bestS[i_1][l];
                for (auto &itemS_i_1 : bestS_i_1) {
                    int index_i_1 = itemS_i_1.first;
                    State<T> stateS_i_1 = itemS_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if (IsLegal(con_seq, nuci_1, nuci, i_1)) {
                        newscore = stateM2_j.score + stateS_i_1.score;
                        int total_cross = stateM2_j.cross_pair_num + stateS_i_1.cross_pair_num;
                        // S+M2->Multi
                        index = GetIndex(k, nuck, nucj);
                        update(bestMulti_j, index, newscore, index_i_1, index_j, stateM2_j.index_3, MANNER_S_M2toMulti, last_pair_pos_M2, total_cross);
                    }
                }
            }
        }
    }
}

/** *************************************************************************************************************** */

template<typename T>
void LCDSfoldCAI_beam(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq, int ires_len, bool is_DN, bool is_beam_speedup){
    //DERNA
    int seq_length = rna_seq.size();
    T newscore;
    int index;
    int preindex;
    double param = is_DN ? (1-lambda) : lambda;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;
    // is_beam_speedup
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;

    vector<int> nuc0_list = Base_table.at(rna_seq[0]); //j list
    for(size_t c0 = 0; c0 < nuc0_list.size(); ++c0){
        int nuc0 = nuc0_list[c0];
        index = GetIndex(0, nuc0, nuc0);
        update(bestN[0], index, 0, -1, MANNER_NONEtoN, 0);
        update(bestF[0], index, 0, -1, MANNER_NONEtoF, 0);
        update(bestS[0][1], index, 0, -1, MANNER_NONEtoS, 0);
    }
    
    for(int j = 1; j < seq_length; ++j){// 1 ~ seq_length
        Processing(seq_length);
        std::unordered_map<int, State<T>>& bestN_j = bestN[j];
        std::vector<std::unordered_map<int, State<T>>>& bestS_j = bestS[j];
        std::unordered_map<int, State<T>>& bestF_j = bestF[j];
        std::unordered_map<int, State<T>>& bestC_j = bestC[j];
        std::unordered_map<int, State<T>>& bestM1_j = bestM1[j];
        std::unordered_map<int, State<T>>& bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>>& bestMulti_j = bestMulti[j];
        // is_beam_speedup
        std::unordered_map<int, State<T>>& bestCS_j = bestCS[j];
        
        vector<int> nucj_list = Base_table.at(rna_seq[j]); //j list
        int j_1 = j - 1;

        // None -> S/N state
        //  N   ->  N  state
        //  N   ->  C  state
        init_and_N_state<T>(bestN, bestS_j, bestC_j, rna_seq, con_seq, ami_seq, nucj_list, j, ires_len, is_DN, param);
        
        BeamPrune(con_seq, rna_seq, bestN_j, bestF, !param);
        
        // S -> S state
        S_to_S_state<T>(bestS, bestS_j, con_seq, ami_seq, nucj_list, j, is_DN, param);

        if(is_beam_speedup){
            for (auto &nucj : nucj_list){
                std::unordered_map<int, State<T>>& bestCS_j_1 = bestCS[j_1];
                // CS -> C
                // S+CS -> C
                CS_to_C_state<T>(bestCS_j_1, bestC_j, bestS, rna_seq, con_seq, ami_seq, nucj, j, ires_len, is_DN, param);
                if(is_DN) {
                    BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
                }

                // C+S -> CS
                C_to_CS_state<T>(bestS_j, bestC, bestCS_j, con_seq, j);
                BeamPrune(con_seq, rna_seq, bestCS_j, bestF, false);
            }
        } else {
            //  C+S  -> C state
            // S+C+S -> C state
            S_C_S_state<T>(bestC, bestC_j, bestS, bestF, rna_seq, con_seq, ami_seq, nucj_list, j, ires_len, is_DN, param, true);
            BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
        }

        // C -> C state
        // S+C -> C state
        C_state<T>(bestC, bestC_j, bestS, rna_seq, con_seq, ami_seq, nucj_list, j, ires_len, is_DN, param);
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
        
        // Multi -> Multi state
        // Multi -> C state
        Multi_state<T>(bestMulti, bestMulti_j, bestC_j, rna_seq, con_seq, ami_seq, nucj_list, j, ires_len, is_DN, param);
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);

        // C -> M1 state
        // M1 + C -> M2 state
        M1_C_to_M2_state<T>(bestM1, bestM1_j, bestM2_j, bestC_j, con_seq, j, is_DN);
        BeamPrune(con_seq, rna_seq, bestM2_j, bestF, false);

        // M1 -> M1 state
        M1_to_M1_state<T>(bestM1, bestM1_j, con_seq, ami_seq, nucj_list, j, is_DN, param);

        // M2 -> M1 state
        // M2 -> Multi state
        // S + M2 -> Multi state
        M2_state<T>(bestS, bestM1_j, bestMulti_j, bestM2_j, con_seq);

        BeamPrune(con_seq, rna_seq, bestM1_j, bestF, false);
        BeamPrune(con_seq, rna_seq, bestMulti_j, bestF, false);

        /*conclusion*/
        //  F  -> F State
        //  C  -> F State
        // F+C -> F State
        C_and_F_state<T>(bestF, bestF_j, bestC_j, con_seq, ami_seq, nucj_list, j, seq_length, is_DN, param);

    } //j end
}

template<typename T>
void LCDSfoldCAI_exact(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq, int ires_len, bool is_DN){
    //DERNA
    int seq_length = rna_seq.size();
    T newscore;
    int index;
    int preindex;
    double param = is_DN ? (1-lambda) : lambda;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    vector<int> nuc0_list = Base_table.at(rna_seq[0]); //j list
    for(size_t c0 = 0; c0 < nuc0_list.size(); ++c0){
        int nuc0 = nuc0_list[c0];
        index = GetIndex(0, nuc0, nuc0);
        update(bestN[0], index, 0, -1, MANNER_NONEtoN, 0);
        update(bestF[0], index, 0, -1, MANNER_NONEtoF, 0);
        update(bestS[0][1], index, 0, -1, MANNER_NONEtoS, 0);
    }
    
    for(int j = 1; j < seq_length; ++j){// 1 ~ seq_length
        Processing(seq_length);
        std::unordered_map<int, State<T>>& bestN_j = bestN[j];
        std::vector<std::unordered_map<int, State<T>>>& bestS_j = bestS[j];
        std::unordered_map<int, State<T>>& bestF_j = bestF[j];
        std::unordered_map<int, State<T>>& bestC_j = bestC[j];
        std::unordered_map<int, State<T>>& bestM1_j = bestM1[j];
        std::unordered_map<int, State<T>>& bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>>& bestMulti_j = bestMulti[j];
        
        vector<int> nucj_list = Base_table.at(rna_seq[j]); //j list
        int j_1 = j - 1;

        // None -> S/N state
        //  N   ->  N  state
        //  N   ->  C  state
        init_and_N_state<T>(bestN, bestS_j, bestC_j, rna_seq, con_seq, ami_seq, nucj_list, j, ires_len, is_DN, param);

        // S -> S state
        S_to_S_state<T>(bestS, bestS_j, con_seq, ami_seq, nucj_list, j, is_DN, param);

        //  C+S  -> C state
        // S+C+S -> C state
        S_C_S_state<T>(bestC, bestC_j, bestS, bestF, rna_seq, con_seq, ami_seq, nucj_list, j, ires_len, is_DN, param);
        
        // C -> C state
        // S+C -> C state
        C_state<T>(bestC, bestC_j, bestS, rna_seq, con_seq, ami_seq, nucj_list, j, ires_len, is_DN, param);
        
        // Multi -> Multi state
        // Multi -> C state
        Multi_state<T>(bestMulti, bestMulti_j, bestC_j, rna_seq, con_seq, ami_seq, nucj_list, j, ires_len, is_DN, param);

        // C -> M1 state
        // M1 + C -> M2 state
        M1_C_to_M2_state<T>(bestM1, bestM1_j, bestM2_j, bestC_j, con_seq, j, is_DN);
        
        // M1 -> M1 state
        M1_to_M1_state<T>(bestM1, bestM1_j, con_seq, ami_seq, nucj_list, j, is_DN, param);

        // M2 -> M1 state
        // M2 -> Multi state
        // S + M2 -> Multi state
        M2_state<T>(bestS, bestM1_j, bestMulti_j, bestM2_j, con_seq);

        /* conclusion */
        //  F  -> F State
        //  C  -> F State
        // F+C -> F State
        C_and_F_state<T>(bestF, bestF_j, bestC_j, con_seq, ami_seq, nucj_list, j, seq_length, is_DN, param);

    } //j end
}

/** *************************************************************************************************************** */

/**
 * Examine the elements in bestTables by logging out the indexs and their corresponding nucleotides. 
 * This is for debugging and analysis purposes, to understand what kind of pairs or structures are being stored in the best tables.
 * 
 * Usage:
 *    examine_bestTable(alltables.bestC, "bestC");
 *    examine_bestTable(alltables.bestM1, "bestM1");
 *    examine_bestTable(alltables.bestM2, "bestM2");
 *    examine_bestTable(alltables.bestMulti, "bestMulti");
 *    examine_bestSTable(alltables.bestS, "bestS");
 */

template<typename T>
void examine_bestTable(std::vector<std::unordered_map<int, State<T>>>& bestTable, const string& tableName){
    std::string fileName = "result_" + tableName + ".csv";
    std::ofstream outputcsv_bestC(fileName,std::ios::app);
    outputcsv_bestC<<"score,index_i,nuci,index_j,nucj,manner"<<std::endl;
    for (size_t i = 0; i < bestTable.size(); ++i) {
        for (auto& [id, state] : bestTable[i]) {
            int pos, nuci, nucj;
            std::tie(pos, nuci, nucj) = GetIndexTuple(id);
            outputcsv_bestC<<state.score<<","<<pos<<","<<reBASE(nuci)<<","<<i<<","<<reBASE(nucj)<<","<<mannerToString(state.MANNER)<<std::endl;
        }
    }
    return;
}

template<typename T>
void examine_bestSTable(std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestTable, const string& tableName){
    std::string fileName = "result_" + tableName + ".csv";
    std::ofstream outputcsv_bestC(fileName,std::ios::app);
    outputcsv_bestC<<"score,index_i,nuci,index_j,nucj,manner"<<std::endl;
    for (size_t i = 0; i < bestTable.size(); ++i) {
        for (size_t l = 1; l <= SINGLE_MAX_LEN; ++l) {
            for (auto& [id, state] : bestTable[i][l]) {
                int pos, nuci, nucj;
                std::tie(pos, nuci, nucj) = GetIndexTuple(id);
                outputcsv_bestC<<state.score<<","<<pos<<","<<reBASE(nuci)<<","<<i<<","<<reBASE(nucj)<<","<<mannerToString(state.MANNER)<<std::endl;
            }
        }
    }
    return;
}

template<typename T>
void Pareto_solution(double threshold1, double threshold2, string& ires_seq, string& rna_seq, vector<int>& con_seq, string& ami_seq, vector<double>& cai_vector, string &output_txt, string &output_csv, bool show_score, string& PID, int protein_len, string& ires_id, int ires_len){
    double left_lambda,right_lambda;
    double left_cai,left_mfe,right_cai,right_mfe;
    timeval start,end;
    std::queue<pair<double,double>> lambda_vector;
    std::unordered_map<double,double> MFE_map,CAI_map;
    lambda_vector.push({0.00001,0.99999});
    while(true){
        std::cout<<"Remain: "<<lambda_vector.size()<<" different lambda"<<std::endl;
        if(!lambda_vector.empty()){ 
            left_lambda = lambda_vector.front().first;
            right_lambda = lambda_vector.front().second;
            lambda_vector.pop();
        }
        else break;// No lambda need to try
        
        if(MFE_map.find(left_lambda) != MFE_map.end()){//left lambda has tried before
            left_cai = CAI_map[left_lambda];
            left_mfe = MFE_map[left_lambda];
        }
        else{
            lambda = left_lambda;
            std::cout << "Protein ID: " << PID << std::endl;
            std::cout << "Protein length: " << protein_len << std::endl;
            std::cout <<"Lambda: " << std::fixed << std::setprecision(3) << lambda << std::endl;
            if(ires_id != ""){
                std::cout <<"IRES ID: " << ires_id << std::endl;
                std::cout <<"IRES length: " << ires_len << std::endl;
                std::cout <<"Crosspair Penalty: " << std::fixed << std::setprecision(3) << crosspair_penalty << std::endl;
            }
            AllTables<T> alltables(rna_seq, rna_seq.size());
            std::string rna_solution;
            std::string structure_solution;
            initialize_CAI_table(cai_vector,true);

            gettimeofday(&start,0);
            initialize_Special_HP<double>(alltables, rna_seq, con_seq, ami_seq, ires_seq.length(), true);
            LCDSfoldCAI_exact<double>(alltables, rna_seq, con_seq, ami_seq, ires_seq.length(), true);
            gettimeofday(&end,0);
            double TimeSpend = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);

            std::vector<pair<double,double>> result = result_output<double>(alltables, rna_seq, con_seq, ami_seq,output_txt,output_csv, PID, protein_len, ires_id, ires_len, TimeSpend, show_score,true);
            CAI_map[left_lambda] = result.front().first;
            MFE_map[left_lambda] = result.front().second;
            left_cai = CAI_map[left_lambda];
            left_mfe = MFE_map[left_lambda];
                
        }

        if(MFE_map.find(right_lambda) != MFE_map.end()){//right lambda has tried before
            right_cai = CAI_map[right_lambda];
            right_mfe = MFE_map[right_lambda];
        }
        else{
            lambda = right_lambda;
            std::cout << "Protein ID: " << PID << std::endl;
            std::cout << "Protein length: " << protein_len << std::endl;
            std::cout <<"Lambda: " << std::fixed << std::setprecision(3) << lambda << std::endl;
            if(ires_id != ""){
                std::cout <<"IRES ID: " << ires_id << std::endl;
                std::cout <<"IRES length: " << ires_len << std::endl;
                std::cout <<"Crosspair Penalty: " << std::fixed << std::setprecision(3) << crosspair_penalty << std::endl;
            }
            AllTables<T> alltables(rna_seq, rna_seq.size());
            std::string rna_solution;
            std::string structure_solution;
            initialize_CAI_table(cai_vector,true);

            gettimeofday(&start,0);
            initialize_Special_HP<double>(alltables, rna_seq, con_seq, ami_seq, ires_seq.length(), true);
            LCDSfoldCAI_exact<double>(alltables, rna_seq, con_seq, ami_seq, ires_seq.length(), true);
            gettimeofday(&end,0);
            double TimeSpend = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);

            std::vector<pair<double,double>> result = result_output<double>(alltables, rna_seq, con_seq, ami_seq,output_txt,output_csv, PID, protein_len, ires_id, ires_len, TimeSpend, show_score, true);
            CAI_map[right_lambda] = result.front().first;
            MFE_map[right_lambda] = result.front().second;
            right_cai = CAI_map[right_lambda];
            right_mfe = MFE_map[right_lambda];

        }
        if (!compare(left_cai, right_cai) && !compare(left_mfe, right_mfe)) {
            if (right_lambda < threshold1) {//0.0025
                if (!compare(left_lambda,right_lambda,threshold2)) {//0.00075
                    double m = (left_lambda + right_lambda) / 2;
                    lambda_vector.push({left_lambda, m});
                    lambda_vector.push({m, right_lambda});
                }
            } else {
                if (!compare(left_lambda,right_lambda,threshold1)) {//0.0025
                    double m = (left_lambda + right_lambda) / 2;
                    lambda_vector.push({left_lambda, m});
                    lambda_vector.push({m, right_lambda});
                }
            }

        }
        
    }
    
    
}

