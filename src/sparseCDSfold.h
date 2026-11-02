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
#include "LCDSfold.h"
using namespace std;
using namespace oligodesign;

double internal_time = 0.0;
double internal_time1 = 0.0;
bool show_processing = true;

bool F_speedup = false;

template <typename T>
void BacktrackStructure_sparse(string &rna_solution, string &structure_solution, int seq_length, int maxindex, State<T> state, AllTables<T> &alltables, int start_j = -1)
{
    std::vector<std::unordered_map<int, State<T>>> &bestN = alltables.bestN;
    std::vector<std::vector<std::array<State<T>, 257>>> &bestSS = alltables.bestSS;
    std::vector<std::vector<std::array<bool, 257>>> &bestSS_check = alltables.bestSS_check;

    std::vector<std::unordered_map<int, State<T>>> &bestF = alltables.bestF;
    std::vector<std::vector<std::unordered_map<int, State<T>>>> &bestSC = alltables.bestSC;
    std::vector<std::unordered_map<int, State<T>>> &bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>> &bestM2 = alltables.bestM2;

    stack<tuple<int, int, State<T>>> stk;
    int j_start = (start_j >= 0) ? start_j : (seq_length - 1);
    stk.push(make_tuple(j_start, maxindex, state));

    while (!stk.empty())
    {
        tuple<int, int, State<T>> top = stk.top();
        int j = get<0>(top);
        int cindex = get<1>(top);
        state = get<2>(top);
        int index1 = state.index_1;
        int index2 = state.index_2;
        int index3 = state.index_3;
        int q = state.index_q;

        int i, nuci, nucj;
        std::tie(i, nuci, nucj) = GetIndexTuple(cindex);
        Manner MANNER = state.MANNER;
        State<T> next_state1, next_state2, next_state3;
        stk.pop();

        int p, nucp, nucq;
        int i_, nuci_, nucj_;
        int index_;

        switch (MANNER)
        {
        case NONE:
            rna_solution[j] = reBASE(nucj);
            structure_solution[j] = '.';
            break;
        case MANNER_N_N_F:
            next_state1 = bestN[j - 1][index1];
            stk.push(make_tuple(j - 1, index1, next_state1));
            rna_solution[j] = reBASE(nucj);
            structure_solution[j] = '.';
            break;
        case MANNER_S_S_F:
            std::tie(i_, nuci_, nucj_) = GetIndexTuple(index1);
            index_ = GetIndex1(nuci_, nucj_);

            next_state1 = bestSS[j - 1][j - i][index_];
            stk.push(make_tuple(j - 1, index1, next_state1));
            rna_solution[j] = reBASE(nucj);
            structure_solution[j] = '.';
            break;
        case MANNER_C_H:
            next_state1 = bestN[j - 1][index1];
            stk.push(make_tuple(j - 1, index1, next_state1));
            rna_solution[i] = reBASE(nuci);
            rna_solution[j] = reBASE(nucj);
            structure_solution[i] = '(';
            structure_solution[j] = ')';
            break;
        case MANNER_C_I:
            std::tie(p, nucp, nucq) = GetIndexTuple(index2);
            if (p - i > 1)
            {
                std::tie(i_, nuci_, nucj_) = GetIndexTuple(index1);
                index_ = GetIndex1(nuci_, nucj_);
                next_state1 = bestSS[p - 1][p - 1 - (i + 1) + 1][index_];
                stk.push(make_tuple(p - 1, index1, next_state1));
            }

            next_state2 = bestSC[q][q - p + 1][index2];
            stk.push(make_tuple(q, index2, next_state2));

            if (j - q > 1)
            {
                std::tie(i_, nuci_, nucj_) = GetIndexTuple(index3);
                index_ = GetIndex1(nuci_, nucj_);
                next_state3 = bestSS[j - 1][j - 1 - (q + 1) + 1][index_];
                stk.push(make_tuple(j - 1, index3, next_state3));
            }

            rna_solution[i] = reBASE(nuci);
            rna_solution[j] = reBASE(nucj);
            structure_solution[i] = '(';
            structure_solution[j] = ')';
            break;
        case MANNER_C_M2:
            next_state1 = bestM2[j - 1][index1];
            stk.push(make_tuple(j - 1, index1, next_state1));
            rna_solution[i] = reBASE(nuci);
            rna_solution[j] = reBASE(nucj);
            structure_solution[i] = '(';
            structure_solution[j] = ')';
            break;
        case MANNER_M2_M2_F:
            next_state1 = bestM2[j - 1][index1];
            stk.push(make_tuple(j - 1, index1, next_state1));
            rna_solution[j] = reBASE(nucj);
            structure_solution[j] = '.';
            break;
        case MANNER_M2_M1_M1:
            next_state1 = bestM1[q][index1];
            stk.push(make_tuple(q, index1, next_state1));
            next_state2 = bestM1[j][index2];
            stk.push(make_tuple(j, index2, next_state2));
            break;
        case MANNER_M2_F_M2:
            next_state1 = bestM2[j][index1];
            stk.push(make_tuple(j, index1, next_state1));
            rna_solution[i] = reBASE(nuci);
            structure_solution[i] = '.';
            break;
        case MANNER_M1_M2:
            next_state1 = bestM2[j][cindex];
            stk.push(make_tuple(j, cindex, next_state1));
            break;
        case MANNER_M1_M1_F:
            next_state1 = bestM1[j - 1][index1];
            stk.push(make_tuple(j - 1, index1, next_state1));
            rna_solution[j] = reBASE(nucj);
            structure_solution[j] = '.';
            break;
        case MANNER_M1_F_M1:
            next_state1 = bestM1[j][index1];
            stk.push(make_tuple(j, index1, next_state1));
            rna_solution[i] = reBASE(nuci);
            structure_solution[i] = '.';
            break;
        case MANNER_M1_C:
            next_state1 = bestSC[j][j - i + 1][cindex];
            stk.push(make_tuple(j, cindex, next_state1));
            break;
        case MANNER_F_C:
            next_state1 = bestSC[j][j - i + 1][cindex];
            stk.push(make_tuple(j, cindex, next_state1));
            break;
        case MANNER_F_F_F:
            next_state1 = bestF[q][index1];
            stk.push(make_tuple(q, index1, next_state1));
            next_state2 = bestF[j][index2];
            stk.push(make_tuple(j, index2, next_state2));
            break;
        case MANNER_F_F_C:
            next_state1 = bestF[q][index1];
            stk.push(make_tuple(q, index1, next_state1));
            next_state2 = bestSC[j][j - (q + 1) + 1][index2];
            stk.push(make_tuple(j, index2, next_state2));
            break;
        case MANNER_S_HP:
            rna_solution[i] = reBASE(nuci);
            rna_solution[j] = reBASE(nucj);
            structure_solution[i] = '(';
            structure_solution[j] = ')';
            for (int x = i + 1; x < j; x++)
            {
                rna_solution[x] = sp_backtrack[cindex][x - i];
                structure_solution[x] = '.';
            }
            break;
        default:
            break;
        }
    }
}

template <typename T>
T BackTrack_sparse(string &rna_solution, string &structure_solution, string &rna_seq, AllTables<T> &alltables, vector<int> &con_seq)
{
    std::vector<std::unordered_map<int, State<T>>> &bestF = alltables.bestF;

    int seq_length = rna_seq.size();
    rna_solution.assign(seq_length, '.');
    structure_solution.assign(seq_length, '.');

    vector<int> first_nuclist = Base_table.at(rna_seq[0]);             // i list
    vector<int> last_nuclist = Base_table.at(rna_seq[seq_length - 1]); // j list

    State<T> state;
    T maxscore = VALUE_MIN<T>();
    int maxindex;

    for (int ci = 0; ci < first_nuclist.size(); ci++) {
        for (int cj = 0; cj < last_nuclist.size(); cj++) {
            int nuci = first_nuclist[ci];
            int nucj = last_nuclist[cj];
            int index = GetIndex(0, nuci, nucj);

            auto it = bestF[seq_length - 1].find(index);
            if (it != bestF[seq_length - 1].end()) {
                T newscore = it->second.score;
                if (maxscore < newscore) {
                    maxscore = newscore;
                    maxindex = index;
                    state = it->second;
                }
            }
        }
    }

    BacktrackStructure_sparse(rna_solution, structure_solution, seq_length, maxindex, state, alltables);

    return maxscore;
}

template <typename T>
bool update(std::unordered_map<int, State<T>> &stateMap, int index, double newscore, int index1, int index2, int index3, int q, Manner MANNER)
{
    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score)
    {
        stateMap[index] = State<T>(newscore, index1, index2, index3, q, MANNER);
        return true;
    }

    return false;
}

template <typename T>
bool update(std::unordered_map<int, State<T>> &stateMap, int index, double newscore, int index1, int index2, int index3, int q, Manner MANNER, int first_pair_pos, int last_pair_pos)
{
    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score)
    {
        stateMap[index] = State<T>(newscore, index1, index2, index3, q, MANNER, first_pair_pos, last_pair_pos);
        return true;
    }

    return false;
}

template <typename T>
bool update(std::array<State<T>, 257> &stateMap, std::array<bool, 257> &checkMap, int index, double newscore, int index1, int index2, int index3, int q, Manner MANNER)
{
    if (!checkMap[index] || newscore > stateMap[index].score)
    {
        stateMap[index] = State<T>(newscore, index1, index2, index3, q, MANNER);
        checkMap[index] = true;
        return true;
    }

    return false;
}

double calTime(const struct timespec start, const struct timespec end)
{
    struct timespec temp;
    double time_used = 0.0;

    if ((end.tv_nsec - start.tv_nsec) < 0)
    {
        temp.tv_sec = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    }
    else
    {
        temp.tv_sec = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }

    time_used = temp.tv_sec + (double)temp.tv_nsec / 1000000000.0;
    return time_used;
}


template <typename T>
void initialize_Special_HP_sparse(AllTables<T> &alltables, string &rna_seq, vector<int> &con_seq, string &ami_seq, int ires_len, bool is_DN)
{
    // std::vector<std::unordered_map<int, State>> &bestC = alltables.bestC;
    std::vector<std::vector<std::unordered_map<int, State<T>>>> &bestSC = alltables.bestSC;

    ACMatcher<char, 256> matcher;

    for (auto sp_loop : sp_loops_vec)
    {
        string s = sp_loop;
        matcher.AddString(s.c_str(), s.length());
    }

    if (!matcher.MakeTree())
    {
        std::cout << "Failed to build the Aho-Corasick tree" << std::endl;
        return;
    }

    unordered_map<int, unordered_set<string>> global_sp; //<index, all special hairpins could appear>
    for (int i = ires_len; i < ami_seq.length() - 3; i++)
    {
        vector<string> &Amino_1 = Amino_to_nucs[ami_seq[i]];
        vector<string> &Amino_2 = Amino_to_nucs[ami_seq[i + 1]];
        vector<string> &Amino_3 = Amino_to_nucs[ami_seq[i + 2]];
        vector<string> &Amino_4 = Amino_to_nucs[ami_seq[i + 3]];
        for (auto &A1 : Amino_1)
            for (auto &A2 : Amino_2)
                for (auto &A3 : Amino_3)
                {
                    for (auto &A4 : Amino_4)
                    {
                        std::string sub_seq = A1 + A2 + A3 + A4;
                        // cout<<sub_seq<<endl;
                        unordered_map<int, int> sp;
                        matcher.SearchPositions(sub_seq.c_str(), i * 3, 0, sub_seq.length() - 1, sp);
                        for (auto ss : sp)
                        {
                            string special_hp = sub_seq.substr(ss.first - i * 3 - ss.second + 1, ss.second);
                            if (global_sp[ss.first].count(special_hp) == 0)
                                global_sp[ss.first].insert(special_hp);
                            // cout<<"position "<<ss.first<<" have sp loop :"<<sub_seq.substr(ss.first - i*3 - ss.second +1,ss.second)<<endl;
                        }
                    }
                }
    }

    for (auto g : global_sp)
    {
        // cout<<"position "<<g.first<<" have sp loop(s) : ";
        for (auto h : g.second)
        {
            vector<int> last_nuc, first_nuc;
            int last_idx = g.first, first_idx = last_idx - h.length() + 1;
            unordered_map<int, double> first_CAI;
            // cout<<"{"<<h <<", "<<sp_loops[h]<< ", ";
            if (isLastNuc(last_idx))
            { // third place in a codon
                string codon_ = h.substr(h.length() - 3, 3);
                last_nuc.push_back(getLastExtendedNuc(re_Amino_label[codon_][0], codon_));
                // cout<<CheckReBASE(getLastExtendedNuc(re_Amino_label[codon_][0],codon_))<<endl;//Not finish
            }
            else
            {
                vector<int> &nuc_list = Base_table.at(rna_seq[last_idx]);
                vector<int> &nuc_list2 = Base_table.at(rna_seq[last_idx - 1]);
                for (auto &nuc : nuc_list)
                {
                    if (reBASE(nuc) == h[h.length() - 1])
                    {
                        for (auto &nuc2 : nuc_list2)
                        {
                            if (reBASE(nuc2) == h[h.length() - 2] && IsLegal(con_seq, nuc2, nuc, last_idx - 1))
                                last_nuc.push_back(nuc);
                        }
                    }
                }
            }

            if (first_idx % 3 == 1)
            { // second place in a codon
                vector<int> &nuc_list = Base_table.at(rna_seq[first_idx]);
                vector<int> &nuc_list2 = Base_table.at(rna_seq[first_idx + 1]);
                for (auto &nuc : nuc_list)
                {
                    if (reBASE(nuc) == h[0])
                    {
                        for (auto &nuc2 : nuc_list2)
                        {
                            if (reBASE(nuc2) == h[1] && IsLegal(con_seq, nuc, nuc2, first_idx))
                            {
                                first_nuc.push_back(nuc);
                                first_CAI[nuc] = getCAI(ami_seq[(first_idx + 1) / 3], nuc2, is_DN);
                            }
                        }
                    }
                }
            }
            else
            {
                vector<int> &nuc_list = Base_table.at(rna_seq[first_idx]);
                for (auto &nuc : nuc_list)
                {
                    if (reBASE(nuc) == h[0])
                    {
                        first_nuc.push_back(nuc);
                    }
                }
            }
            std::vector<std::unordered_map<int, State<T>>> &bestSC_j = bestSC[last_idx];
            for (auto &nuci : first_nuc)
            {
                for (auto &nucj : last_nuc)
                {
                    double newscore = -MFE_scalar * sp_loops[h];
                    for (int x = first_idx; x <= last_idx; x++)
                    {
                        if (isLastNuc(x) && (CAI_scalar))
                        {
                            if (x - first_idx == 0)
                            { // special case1
                                newscore += getCAI(ami_seq[x / 3], nuci, is_DN);
                            }
                            else if (x - first_idx == 1)
                            { // special case2
                                newscore += first_CAI[nuci];
                            }
                            else
                            {
                                string codon_ = h.substr(x - first_idx - 2, 3);
                                int nucx = getLastExtendedNuc(re_Amino_label[codon_][0], codon_);
                                newscore += getCAI(ami_seq[x / 3], nucx, is_DN);
                            }
                        }
                    }
                    int index = GetIndex(first_idx, nuci, nucj);
                    int l = last_idx - first_idx;
                    update(bestSC_j[l + 1], index, newscore, -1, -1, -1, -1, MANNER_S_HP);
                    if (bestSC_j[l + 1][index].score == newscore)
                        sp_backtrack[index] = h;
                }
            }
        }
    }
}


template <typename T>
void SparseCDSfoldCAI(AllTables<T> &alltables, string &rna_seq, vector<int> &con_seq, string &ami_seq, bool is_DN)
{
    int seq_length = rna_seq.size();
    double newscore, score_i1_p_1, score_q1_j_1;
    double nuci_score, nucj_score;
    double f_coterminus_score, m1_coterminus_score;
    int index, index_i_i, index_i1_j_1, index_i1, index_k_1, index_q_1, index_j_1, index_j;
    int index1;
    int index_p_q, index_i1_p_1, index_q1_j_1;
    int index_k_j, index_i_1, k, nuck, nuci_1;

    int index_p_q_, index_i1_p_, index_q1_j_;

    int maxp, minq, maxk;
    bool is_F_OCT, is_M1_OCT, is_F_STEP, is_M1_STEP;

    struct timespec start, end;
    struct timespec start1, end1;
    // double internal_time = 0.0;
    double multi_time = 0.0;
    double F_time = 0.0;

    State<T> state_i_1;

    std::vector<std::unordered_map<int, State<T>>> &bestN = alltables.bestN;
    // std::vector<std::vector<std::unordered_map<int, State>>> &bestS = alltables.bestS;
    std::vector<std::vector<std::array<State<T>, 257>>> &bestSS = alltables.bestSS;
    std::vector<std::vector<std::array<bool, 257>>> &bestSS_check = alltables.bestSS_check;
    std::vector<std::unordered_map<int, State<T>>> &bestF = alltables.bestF;
    std::vector<std::vector<std::unordered_map<int, State<T>>>> &bestSC = alltables.bestSC;
    // std::vector<std::vector<std::array<State, 257>>> &bestC = alltables.bestC;
    // std::vector<std::vector<std::array<bool, 257>>> &bestC_check = alltables.bestC_check;
    std::vector<std::unordered_map<int, State<T>>> &bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>> &bestM2 = alltables.bestM2;

    int size_ = seq_length * (1 + seq_length) / 2;
    for (int j = 0; j < seq_length; j++)
    {
        std::unordered_map<int, State<T>> &bestN_j = bestN[j];
        // std::vector<std::unordered_map<int, State>> &bestS_j = bestS[j];
        std::vector<std::array<State<T>, 257>> &bestSS_j = bestSS[j];
        std::vector<std::array<bool, 257>> &bestSS_check_j = bestSS_check[j];
        std::unordered_map<int, State<T>> &bestF_j = bestF[j];
        std::vector<std::unordered_map<int, State<T>>> &bestSC_j = bestSC[j];
        // std::vector<std::array<State, 257>> &bestC_j = bestC[j];
        // std::vector<std::array<bool, 257>> &bestC_check_j = bestC_check[j];
        std::unordered_map<int, State<T>> &bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>> &bestM1_j = bestM1[j];
        vector<int> &nucj_list = Base_table.at(rna_seq[j]); // j list
        const std::vector<int> *nucj_1_list = nullptr;
        if (j > 0)
            nucj_1_list = &Base_table.at(rna_seq[j - 1]);

        for (int i = j; i >= 0; i--)
        {
            if (show_processing)
                Processing(size_);

            int l = j - i; // 兩端距離

            vector<int> &nuci_list = Base_table.at(rna_seq[i]); // i list
            const std::vector<int> *nuci1_list = nullptr;
            if (i < seq_length - 1)
                nuci1_list = &Base_table.at(rna_seq[i + 1]);

            for (auto &nucj : nucj_list) // iterate all nuclitide in position j
            {
                for (auto &nuci : nuci_list) // iterate all nuclitide in position i
                {
                    /*the start of (i, j, nuci, nucj)*/

                    if ((i == j) && (nuci != nucj))
                    {
                        continue; // illegal single nucleotide
                    }

                    index = GetIndex(i, nuci, nucj);
                    index1 = GetIndex1(nuci, nucj); // 0-256

                    is_F_OCT = false;
                    is_M1_OCT = false;
                    is_F_STEP = true;
                    is_M1_STEP = true;

                    if (isLastNuc(i))
                    {
                        nuci_score = getCAI(ami_seq[i / 3], nuci, is_DN);
                    }
                    else
                    {
                        nuci_score = 0.0;
                    }

                    if (isLastNuc(j))
                    {
                        nucj_score = getCAI(ami_seq[j / 3], nucj, is_DN);
                    }
                    else
                    {
                        nucj_score = 0.0;
                    }

                    // calculate the N, which is unpair
                    if (i == j)
                    {
                        newscore = nuci_score;
                        update(bestN_j, index, newscore, -1, -1, -1, -1, NONE);
                    }
                    else
                    {
                        // N(i,j) = N(i,j-1) + N(j,j)
                        std::unordered_map<int, State<T>> &bestN_j_1 = bestN[j - 1];
                        for (auto &nucj_1 : *nucj_1_list)
                        {
                            if (!IsLegal(con_seq, nucj_1, nucj, j - 1))
                            {
                                continue;
                            }

                            index_j_1 = GetIndex(i, nuci, nucj_1); // the index of N(i, j-1)
                            auto it = bestN_j_1.find(index_j_1);
                            if (it != bestN_j_1.end())
                            {
                                newscore = it->second.score;
                                newscore += nucj_score;
                                update(bestN_j, index, newscore, index_j_1, -1, -1, -1, MANNER_N_N_F);
                            }
                        }
                    }

                    // calculate the S, which is unpair
                    if (i == j)
                    {
                        newscore = nuci_score;
                        // update(bestS_j[1], index, newscore, -1, -1, -1, -1, NONE);
                        update(bestSS_j[1], bestSS_check_j[1], index1, newscore, -1, -1, -1, -1, NONE);
                    }
                    else
                    {
                        // S(i,j) = S(i,j-1) + S(j,j)
                        if (l <= SINGLE_MAX_LEN - 1) //-1是因為S的開頭會是 i或是q
                        {
                            // std::unordered_map<int, State> &bestS_j_1 = bestS[j - 1][l];
                            std::array<State<T>, 257> &bestSS_j_1 = bestSS[j - 1][l];
                            std::array<bool, 257> &bestSS_check_j_1 = bestSS_check[j - 1][l];
                            for (auto &nucj_1 : *nucj_1_list)
                            {
                                if (!IsLegal(con_seq, nucj_1, nucj, j - 1))
                                {
                                    continue;
                                }

                                index_j_1 = GetIndex(i, nuci, nucj_1); // the index of S(i, j-1)
                                // auto it = bestSS_j_1.find(index_j_1);
                                // if (it != bestSS_j_1.end())
                                // {
                                //     newscore = it->second.score;
                                //     newscore += nucj_score;
                                //     update(bestSS_j[l + 1], index, newscore, index_j_1, -1, -1, -1, MANNER_S_S_F);
                                // }
                                int index_ = GetIndex1(nuci, nucj_1); // the index of S(i, j-1)
                                if (bestSS_check_j_1[index_])
                                {
                                    newscore = bestSS_j_1[index_].score;
                                    newscore += nucj_score;
                                    update(bestSS_j[l + 1], bestSS_check_j[l + 1], index1, newscore, index_j_1, -1, -1, -1, MANNER_S_S_F);
                                }
                            }
                        }
                    }

                    // hairpin
                    if ((j - i > HAIRPIN_GAP) && _allowed_pairs[nuci][nucj])
                    {
                        std::unordered_map<int, State<T>> &bestN_j_1 = bestN[j - 1];
                        for (auto &nuci1 : *nuci1_list)
                        {
                            if (!IsLegal(con_seq, nuci, nuci1, i))
                            {
                                continue;
                            }

                            for (auto &nucj_1 : *nucj_1_list)
                            {
                                if (!IsLegal(con_seq, nucj_1, nucj, j - 1))
                                {
                                    continue;
                                }

                                index_i1_j_1 = GetIndex(i + 1, nuci1, nucj_1);
                                auto it = bestN_j_1.find(index_i1_j_1);
                                if (it != bestN_j_1.end())
                                {
                                    newscore = -MFE_scalar * v_score_hairpin(i, j, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj));
                                    newscore += it->second.score;
                                    newscore += nuci_score;
                                    newscore += nucj_score;
                                    update(bestSC_j[l + 1], index, newscore, index_i1_j_1, -1, -1, -1, MANNER_C_H);
                                    // update(bestC_j[l + 1], bestC_check_j[l + 1], index1, newscore, index_i1_j_1, -1, -1, -1, MANNER_C_H);
                                }
                            }
                        }
                    }

                    clock_gettime(CLOCK_MONOTONIC, &start);

                    // stacking
                    if ((j - i > MIN_STACKPAIR_GAP) && _allowed_pairs[nuci][nucj])
                    {
                        int p = i + 1;
                        int q = j - 1;

                        std::array<State<T>, 257> &bestSS_p_1 = bestSS[p - 1][p - 1 - (i + 1) + 1];
                        std::array<bool, 257> &bestSS_check_p_1 = bestSS_check[p - 1][p - 1 - (i + 1) + 1];
                        vector<int> &nucp_list = Base_table.at(rna_seq[p]);

                        std::array<State<T>, 257> &bestSS_j_1 = bestSS[j - 1][j - 1 - (q + 1) + 1];
                        std::array<bool, 257> &bestSS_check_j_1 = bestSS_check[j - 1][j - 1 - (q + 1) + 1];
                        vector<int> &nucq_list = Base_table.at(rna_seq[q]);

                        std::unordered_map<int, State<T>> &bestSC_q = bestSC[q][q - p + 1];

                        for (auto &nucp : nucp_list)
                        {
                            if (!IsLegal(con_seq, nuci, nucp, i))
                            {
                                continue;
                            }

                            for (auto &nucq : nucq_list)
                            {
                                if (!_allowed_pairs[nucp][nucq])
                                {
                                    continue;
                                }

                                if (!IsLegal(con_seq, nucq, nucj, q))
                                {
                                    continue;
                                }

                                // for each intloops(nuci, nucp, nucq, nucj)
                                index_p_q = GetIndex(p, nucp, nucq);
                                index_p_q_ = GetIndex1(nucp, nucq);
                                auto it = bestSC_q.find(index_p_q);
                                if (it == bestSC_q.end())
                                {
                                    continue;
                                }

                                newscore = -MFE_scalar * v_score_single_stacking(EnrBASE(nuci), EnrBASE(nucj), EnrBASE(nucp), EnrBASE(nucq));

                                newscore += it->second.score;
                                newscore += nuci_score;
                                newscore += nucj_score;
                                update(bestSC_j[l + 1], index, newscore, -1, index_p_q, -1, q, MANNER_C_I);
                            }
                        }
                    }

                    // left bulge
                    if ((j - i > MIN_STACKPAIR_GAP) && _allowed_pairs[nuci][nucj])
                    {
                        int p = i + 1;

                        std::array<State<T>, 257> &bestSS_p_1 = bestSS[p - 1][p - 1 - (i + 1) + 1];
                        std::array<bool, 257> &bestSS_check_p_1 = bestSS_check[p - 1][p - 1 - (i + 1) + 1];
                        vector<int> &nucp_list = Base_table.at(rna_seq[p]);
                        vector<int> &nucp_1_list = Base_table.at(rna_seq[p - 1]);

                        minq = MAX2(p + 1 + HAIRPIN_GAP, j - 1 - SINGLE_MAX_LEN);
                        for (int q = minq; q < j; q++)
                        {
                            if (p == i + 1 && q == j - 1)
                                continue;

                            std::array<State<T>, 257> &bestSS_j_1 = bestSS[j - 1][j - 1 - (q + 1) + 1];
                            std::array<bool, 257> &bestSS_check_j_1 = bestSS_check[j - 1][j - 1 - (q + 1) + 1];
                            vector<int> &nucq_list = Base_table.at(rna_seq[q]);
                            vector<int> &nucq1_list = Base_table.at(rna_seq[q + 1]);

                            std::unordered_map<int, State<T>> &bestSC_q = bestSC[q][q - p + 1];
                            for (auto &nucp : nucp_list)
                            {
                                if (!IsLegal(con_seq, nuci, nucp, i))
                                {
                                    continue;
                                }

                                for (auto &nucq : nucq_list)
                                {
                                    if (!_allowed_pairs[nucp][nucq])
                                    {
                                        continue;
                                    }

                                    if ((j - q == 1) && !IsLegal(con_seq, nucq, nucj, q))
                                    {
                                        continue;
                                    }

                                    // for each intloops(nuci, nucp, nucq, nucj)
                                    index_p_q = GetIndex(p, nucp, nucq);
                                    index_p_q_ = GetIndex1(nucp, nucq);
                                    auto it = bestSC_q.find(index_p_q);
                                    if (it == bestSC_q.end())
                                    {
                                        continue;
                                    }

                                    score_i1_p_1 = 0;
                                    for (auto &nucj_1 : *nucj_1_list)
                                    {
                                        if (!IsLegal(con_seq, nucj_1, nucj, j - 1))
                                            continue;

                                        for (auto &nucq1 : nucq1_list)
                                        {
                                            if (!IsLegal(con_seq, nucq, nucq1, q))
                                                continue;

                                            index_q1_j_ = GetIndex1(nucq1, nucj_1);
                                            if (j - q == 1)
                                            {
                                                if ((nucj_1 != nucq) || (nucq1 != nucj))
                                                {
                                                    continue;
                                                }
                                                score_q1_j_1 = 0;
                                            }
                                            else if (j - q == 2)
                                            {
                                                if ((nucj_1 != nucq1) || !bestSS_check_j_1[index_q1_j_])
                                                {
                                                    continue;
                                                }
                                                score_q1_j_1 = bestSS_j_1[index_q1_j_].score;
                                            }
                                            else
                                            {
                                                if (!bestSS_check_j_1[index_q1_j_])
                                                {
                                                    continue;
                                                }
                                                score_q1_j_1 = bestSS_j_1[index_q1_j_].score;
                                            }

                                            index_q1_j_1 = GetIndex(q + 1, nucq1, nucj_1);

                                            newscore = -MFE_scalar * v_score_single_bulge(i, j, p, q, EnrBASE(nuci), EnrBASE(nucp), EnrBASE(nucj_1), EnrBASE(nucj), EnrBASE(nuci), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucq1));

                                            newscore += score_i1_p_1 + it->second.score + score_q1_j_1;
                                            newscore += nuci_score;
                                            newscore += nucj_score;
                                            update(bestSC_j[l + 1], index, newscore, -1, index_p_q, index_q1_j_1, q, MANNER_C_I);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // right bulge
                    if ((j - i > MIN_STACKPAIR_GAP) && _allowed_pairs[nuci][nucj])
                    {
                        int q = j - 1;
                        std::array<State<T>, 257> &bestSS_j_1 = bestSS[j - 1][j - 1 - (q + 1) + 1];
                        std::array<bool, 257> &bestSS_check_j_1 = bestSS_check[j - 1][j - 1 - (q + 1) + 1];
                        vector<int> &nucq_list = Base_table.at(rna_seq[q]);
                        vector<int> &nucq1_list = Base_table.at(rna_seq[q + 1]);

                        maxp = MIN2(q - 1 - HAIRPIN_GAP, i + SINGLE_MAX_LEN + 1);
                        for (int p = i + 1; p <= maxp; p++)
                        {
                            std::array<State<T>, 257> &bestSS_p_1 = bestSS[p - 1][p - 1 - (i + 1) + 1];
                            std::array<bool, 257> &bestSS_check_p_1 = bestSS_check[p - 1][p - 1 - (i + 1) + 1];
                            vector<int> &nucp_list = Base_table.at(rna_seq[p]);
                            vector<int> &nucp_1_list = Base_table.at(rna_seq[p - 1]);

                            std::unordered_map<int, State<T>> &bestSC_q = bestSC[q][q - p + 1];
                            for (auto &nucp : nucp_list)
                            {
                                if ((p - i == 1) && !IsLegal(con_seq, nuci, nucp, i))
                                {
                                    continue;
                                }

                                for (auto &nucq : nucq_list)
                                {
                                    if (!_allowed_pairs[nucp][nucq])
                                    {
                                        continue;
                                    }

                                    if (!IsLegal(con_seq, nucq, nucj, q))
                                    {
                                        continue;
                                    }

                                    // for each intloops(nuci, nucp, nucq, nucj)
                                    index_p_q = GetIndex(p, nucp, nucq);
                                    index_p_q_ = GetIndex1(nucp, nucq);
                                    auto it = bestSC_q.find(index_p_q);
                                    if (it == bestSC_q.end())
                                    {
                                        continue;
                                    }

                                    for (auto &nuci1 : *nuci1_list)
                                    {
                                        if (!IsLegal(con_seq, nuci, nuci1, i))
                                            continue;
                                        for (auto &nucp_1 : nucp_1_list)
                                        {
                                            if (!IsLegal(con_seq, nucp_1, nucp, p - 1))
                                                continue;

                                            index_i1_p_ = GetIndex1(nuci1, nucp_1);

                                            if (p - i == 1)
                                            {
                                                if ((nuci1 != nucp) || (nucp_1 != nuci))
                                                {
                                                    continue;
                                                }
                                                score_i1_p_1 = 0;
                                            }
                                            else if (p - i == 2)
                                            {
                                                if ((nuci1 != nucp_1) || !bestSS_check_p_1[index_i1_p_])
                                                {
                                                    continue;
                                                }
                                                score_i1_p_1 = bestSS_p_1[index_i1_p_].score;
                                            }
                                            else
                                            {
                                                if (!bestSS_check_p_1[index_i1_p_])
                                                {
                                                    continue;
                                                }
                                                score_i1_p_1 = bestSS_p_1[index_i1_p_].score;
                                            }

                                            score_q1_j_1 = 0;

                                            index_i1_p_1 = GetIndex(i + 1, nuci1, nucp_1);

                                            newscore = -MFE_scalar * v_score_single_bulge(i, j, p, q, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucq), EnrBASE(nucj), EnrBASE(nucp_1), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucj));

                                            newscore += score_i1_p_1 + it->second.score + score_q1_j_1;
                                            newscore += nuci_score;
                                            newscore += nucj_score;
                                            update(bestSC_j[l + 1], index, newscore, index_i1_p_1, index_p_q, -1, q, MANNER_C_I);
                                            // update(bestC_j[l + 1], bestC_check_j[l + 1], index1, newscore, index_i1_p_1, index_p_q, index_q1_j_1, q, MANNER_C_I);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // internal
                    // C(i,j) = S(i+1,p-1) + C(p,q) + S(q+1,j-1)
                    if ((j - i > MIN_STACKPAIR_GAP) && _allowed_pairs[nuci][nucj])
                    {
                        maxp = MIN2(j - 2 - HAIRPIN_GAP, i + SINGLE_MAX_LEN + 1);
                        for (int p = i + 2; p <= maxp; p++)
                        {
                            std::array<State<T>, 257> &bestSS_p_1 = bestSS[p - 1][p - 1 - (i + 1) + 1];
                            std::array<bool, 257> &bestSS_check_p_1 = bestSS_check[p - 1][p - 1 - (i + 1) + 1];

                            vector<int> &nucp_list = Base_table.at(rna_seq[p]);
                            vector<int> &nucp_1_list = Base_table.at(rna_seq[p - 1]);

                            minq = MAX2(p + 1 + HAIRPIN_GAP, j - i + p - SINGLE_MAX_LEN - 2);
                            for (int q = minq; q < j - 1; q++)
                            {
                                std::array<State<T>, 257> &bestSS_j_1 = bestSS[j - 1][j - 1 - (q + 1) + 1];
                                std::array<bool, 257> &bestSS_check_j_1 = bestSS_check[j - 1][j - 1 - (q + 1) + 1];

                                vector<int> &nucq_list = Base_table.at(rna_seq[q]);
                                vector<int> &nucq1_list = Base_table.at(rna_seq[q + 1]);

                                std::unordered_map<int, State<T>> &bestSC_q = bestSC[q][q - p + 1];

                                for (auto &nucp : nucp_list)
                                {
                                    for (auto &nucq : nucq_list)
                                    {
                                        if (!_allowed_pairs[nucp][nucq])
                                        {
                                            continue;
                                        }

                                        // for each intloops(nuci, nucp, nucq, nucj)
                                        index_p_q = GetIndex(p, nucp, nucq);
                                        index_p_q_ = GetIndex1(nucp, nucq);
                                        auto it = bestSC_q.find(index_p_q);
                                        if (it == bestSC_q.end())
                                        {
                                            continue;
                                        }

                                        for (auto &nuci1 : *nuci1_list)
                                        {
                                            if (!IsLegal(con_seq, nuci, nuci1, i))
                                                continue;
                                            for (auto &nucp_1 : nucp_1_list)
                                            {
                                                if (!IsLegal(con_seq, nucp_1, nucp, p - 1))
                                                    continue;

                                                index_i1_p_ = GetIndex1(nuci1, nucp_1);

                                                if (p - i == 2)
                                                {
                                                    if ((nuci1 != nucp_1) || !bestSS_check_p_1[index_i1_p_])
                                                    {
                                                        continue;
                                                    }
                                                    score_i1_p_1 = bestSS_p_1[index_i1_p_].score;
                                                }
                                                else
                                                {
                                                    if (!bestSS_check_p_1[index_i1_p_])
                                                    {
                                                        continue;
                                                    }
                                                    score_i1_p_1 = bestSS_p_1[index_i1_p_].score;
                                                }

                                                for (auto &nucj_1 : *nucj_1_list)
                                                {
                                                    if (!IsLegal(con_seq, nucj_1, nucj, j - 1))
                                                        continue;

                                                    for (auto &nucq1 : nucq1_list)
                                                    {
                                                        if (!IsLegal(con_seq, nucq, nucq1, q))
                                                            continue;

                                                        index_q1_j_ = GetIndex1(nucq1, nucj_1);
                                                        if (j - q == 2)
                                                        {
                                                            if ((nucj_1 != nucq1) || !bestSS_check_j_1[index_q1_j_])
                                                            {
                                                                continue;
                                                            }
                                                            score_q1_j_1 = bestSS_j_1[index_q1_j_].score;
                                                        }
                                                        else
                                                        {
                                                            if (!bestSS_check_j_1[index_q1_j_])
                                                            {
                                                                continue;
                                                            }
                                                            score_q1_j_1 = bestSS_j_1[index_q1_j_].score;
                                                        }

                                                        index_i1_p_1 = GetIndex(i + 1, nuci1, nucp_1);
                                                        index_q1_j_1 = GetIndex(q + 1, nucq1, nucj_1);

                                                        newscore = -MFE_scalar * v_score_single_internal(i, j, p, q, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj), EnrBASE(nucp_1), EnrBASE(nucp), EnrBASE(nucq), EnrBASE(nucq1));

                                                        newscore += score_i1_p_1 + it->second.score + score_q1_j_1;
                                                        // newscore += score_i1_p_1 + bestC_q[index_p_q_].score + score_q1_j_1;
                                                        newscore += nuci_score;
                                                        newscore += nucj_score;
                                                        update(bestSC_j[l + 1], index, newscore, index_i1_p_1, index_p_q, index_q1_j_1, q, MANNER_C_I);
                                                        // update(bestC_j[l + 1], bestC_check_j[l + 1], index1, newscore, index_i1_p_1, index_p_q, index_q1_j_1, q, MANNER_C_I);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    clock_gettime(CLOCK_MONOTONIC, &end);
                    internal_time += calTime(start, end);

                    // clock_gettime(CLOCK_MONOTONIC, &start);
                    // multi-loop
                    // C(i,j) = M2(i+1,j-1) + ML_closing
                    if ((j - i > MIN_STACKPAIR_GAP * 2) && _allowed_pairs[nuci][nucj])
                    {
                        std::unordered_map<int, State<T>> &bestM2_j_1 = bestM2[j - 1];
                        for (auto &nuci1 : *nuci1_list)
                        {
                            for (auto &nucj_1 : *nucj_1_list)
                            {
                                if (IsLegal(con_seq, nuci, nuci1, i) && IsLegal(con_seq, nucj_1, nucj, j - 1))
                                {
                                    index_i1_j_1 = GetIndex(i + 1, nuci1, nucj_1);
                                    auto it = bestM2_j_1.find(index_i1_j_1);
                                    if (it != bestM2_j_1.end())
                                    {
                                        newscore = -MFE_scalar * v_score_multi(EnrBASE(nuci), EnrBASE(nucj));
                                        newscore += it->second.score;
                                        newscore += nuci_score;
                                        newscore += nucj_score;
                                        update(bestSC_j[l + 1], index, newscore, index_i1_j_1, -1, -1, -1, MANNER_C_M2);
                                        // update(bestC_j[l + 1], bestC_check_j[l + 1], index1, newscore, index_i1_j_1, -1, -1, -1, MANNER_C_M2);
                                    }
                                }
                            }
                        }
                    }

                    // calculate the M2 which i is unpair and j is unpair
                    if (j - i > MIN_STACKPAIR_GAP * 2 - 1)
                    { // M2(i,j) = M2(i,j-1) + F(j,j)
                        std::unordered_map<int, State<T>> &bestM2_j_1 = bestM2[j - 1];
                        for (auto &nucj_1 : *nucj_1_list)
                        {
                            if (IsLegal(con_seq, nucj_1, nucj, j - 1))
                            {
                                index_j_1 = GetIndex(i, nuci, nucj_1);
                                auto it = bestM2_j_1.find(index_j_1);
                                if (it != bestM2_j_1.end())
                                {
                                    int first_pair_pos = it->second.first_pair_pos;
                                    int last_pair_pos = it->second.last_pair_pos;

                                    if (j - last_pair_pos <= SINGLE_MAX_LEN)
                                    {
                                        newscore = it->second.score;
                                        newscore += nucj_score;
                                        update(bestM2_j, index, newscore, index_j_1, -1, -1, -1, MANNER_M2_M2_F, first_pair_pos, last_pair_pos);
                                    }
                                }
                            }
                        }

                        // M2(i,j) = F(i,i) + M2(i+1,j)
                        for (auto &nuci1 : *nuci1_list)
                        {
                            if (IsLegal(con_seq, nuci, nuci1, i))
                            {
                                index_j = GetIndex(i + 1, nuci1, nucj);
                                auto it = bestM2_j.find(index_j);
                                if (it != bestM2_j.end())
                                {
                                    int first_pair_pos = it->second.first_pair_pos;
                                    int last_pair_pos = it->second.last_pair_pos;
                                    if (first_pair_pos - i <= SINGLE_MAX_LEN)
                                    {
                                        newscore = it->second.score;
                                        newscore += nuci_score;
                                        update(bestM2_j, index, newscore, index_j, -1, -1, -1, MANNER_M2_F_M2, first_pair_pos, last_pair_pos);
                                    }
                                }
                            }
                        }
                    }

                    // calculate the M1
                    // M1(i,j) = M2(i,j)
                    if (j - i > MIN_STACKPAIR_GAP * 2 - 2)
                    {
                        auto it = bestM2_j.find(index);
                        if (it != bestM2_j.end())
                        {
                            int first_pair_pos = it->second.first_pair_pos;
                            int last_pair_pos = it->second.last_pair_pos;
                            newscore = it->second.score;
                            update(bestM1_j, index, newscore, -1, -1, -1, -1, MANNER_M1_M2, first_pair_pos, last_pair_pos);
                        }
                    }

                    // M1(i,j) = M1(i,j-1)+F(j,j)
                    if ((j - i > HAIRPIN_GAP + 1))
                    {
                        std::unordered_map<int, State<T>> &bestM1_j_1 = bestM1[j - 1];
                        for (auto &nucj_1 : *nucj_1_list)
                        {
                            if (IsLegal(con_seq, nucj_1, nucj, j - 1))
                            {
                                index_j_1 = GetIndex(i, nuci, nucj_1);
                                auto it = bestM1_j_1.find(index_j_1);
                                if (it != bestM1_j_1.end())
                                {
                                    int first_pair_pos = it->second.first_pair_pos;
                                    int last_pair_pos = it->second.last_pair_pos;
                                    if (j - last_pair_pos <= SINGLE_MAX_LEN)
                                    {
                                        newscore = it->second.score;
                                        newscore += nucj_score;
                                        update(bestM1_j, index, newscore, index_j_1, -1, -1, -1, MANNER_M1_M1_F, first_pair_pos, last_pair_pos);
                                    }
                                }
                            }
                        }
                    }

                    // M1(i,j) = F(i,i)+M1(i+1,j)
                    if ((j - i > HAIRPIN_GAP + 1))
                    {
                        for (auto &nuci1 : *nuci1_list)
                        {
                            if (IsLegal(con_seq, nuci, nuci1, i))
                            {
                                index_j = GetIndex(i + 1, nuci1, nucj);
                                auto it = bestM1_j.find(index_j);
                                if (it != bestM1_j.end())
                                {
                                    int first_pair_pos = it->second.first_pair_pos;
                                    int last_pair_pos = it->second.last_pair_pos;
                                    if (first_pair_pos - i <= SINGLE_MAX_LEN)
                                    {
                                        newscore = it->second.score;
                                        newscore += nuci_score;
                                        auto it1 = bestM1_j.find(index);
                                        if ((it1 == bestM1_j.end()) || (it1->second.score <= newscore))
                                        {
                                            update(bestM1_j, index, newscore, index_j, -1, -1, -1, MANNER_M1_F_M1, first_pair_pos, last_pair_pos);
                                            is_M1_STEP = false;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // M1(i,j) = C(i,j)+ML_stem
                    if ((j - i > HAIRPIN_GAP) && _allowed_pairs[nuci][nucj])
                    {
                        auto it = bestSC_j[l + 1].find(index);
                        if (it != bestSC_j[l + 1].end())
                        // if (bestC_check_j[l + 1][index1])
                        {
                            m1_coterminus_score = it->second.score - MFE_scalar * v_score_M1(EnrBASE(nuci), EnrBASE(nucj));
                            // m1_coterminus_score = bestC_j[l + 1][index1].score - v_score_M1(EnrBASE(nuci), EnrBASE(nucj));
                            if (m1_coterminus_score > bestM1_j[index].score)
                            {
                                // update(bestM1_j, index, m1_coterminus_score, -1, -1, -1, -1, MANNER_M1_C);
                                update(bestM1_j, index, m1_coterminus_score, -1, -1, -1, -1, MANNER_M1_C, i, j);
                                is_M1_OCT = true;
                                is_M1_STEP = true;
                            }
                        }
                    }

                    if (is_M1_STEP)
                    {
                        auto it = bestM1_j.find(index);
                        if (it != bestM1_j.end())
                        {
                            it->second.is_M1_STEP = true;
                        }
                    }

                    // Forward Update M2(k,j),  0 <= k <= i-1
                    if (is_M1_OCT && (i > MIN_STACKPAIR_GAP - 1))
                    { // M2(k,j) = M1(k,i-1)+ M1(i,j)    ( STEP + OCT ) Forward Update
                        std::unordered_map<int, State<T>> &bestM1_i_1 = bestM1[i - 1];
                        for (auto &state_k : bestM1_i_1)
                        { // M1(1,i-1), M1(2,i-1),...,M1(k,i-1),...,M1(i-5,i-1)
                            if (!state_k.second.is_M1_STEP)
                                continue;
                            std::tie(k, nuck, nuci_1) = GetIndexTuple(state_k.first);
                            if (!IsLegal(con_seq, nuci_1, nuci, i - 1))
                                continue;

                            int first_pair_pos = state_k.second.first_pair_pos;

                            index_k_j = GetIndex(k, nuck, nucj);
                            newscore = state_k.second.score + m1_coterminus_score;
                            // update(bestM2_j, index_k_j, newscore, state_k.first, index, -1, i - 1, MANNER_M2_M1_M1);
                            update(bestM2_j, index_k_j, newscore, state_k.first, index, -1, i - 1, MANNER_M2_M1_M1, first_pair_pos, j);
                        }
                    }

                    // clock_gettime(CLOCK_MONOTONIC, &end);
                    // multi_time += calTime(start, end);

                    if (!F_speedup)
                    {
                        // clock_gettime(CLOCK_MONOTONIC, &start);
                        // calculate the F
                        // F(i,j) = F(i,i) + F(i+1, j), the non-STEP case
                        if (j > i + 1)
                        {
                            index_i_i = GetIndex(i, nuci, nuci);
                            for (auto &nuci1 : *nuci1_list)
                            {
                                if (IsLegal(con_seq, nuci, nuci1, i))
                                {
                                    index_i1 = GetIndex(i + 1, nuci1, nucj);

                                    auto it = bestF_j.find(index_i1);
                                    if (it != bestF_j.end())
                                    {
                                        newscore = bestF[i][index_i_i].score + it->second.score;

                                        auto it1 = bestF_j.find(index);
                                        if ((it1 == bestF_j.end()) || (it1->second.score <= newscore))
                                        {
                                            bool res = update(bestF_j, index, newscore, index_i_i, index_i1, -1, i, MANNER_F_F_F);
                                            is_F_STEP = false;

                                            if (res)
                                            {
                                                bestF_j[index].outer_pairs = it->second.outer_pairs;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else if (j == i + 1)
                        {
                            if (IsLegal(con_seq, nuci, nucj, i))
                            {
                                index_i_i = GetIndex(i, nuci, nuci);
                                index_i1 = GetIndex(i + 1, nucj, nucj);
                                newscore = bestF[i][index_i_i].score + bestF_j[index_i1].score;

                                auto it1 = bestF_j.find(index);
                                if ((it1 == bestF_j.end()) || (it1->second.score <= newscore))
                                {
                                    update(bestF_j, index, newscore, index_i_i, index_i1, -1, i, MANNER_F_F_F);
                                    is_F_STEP = false;
                                }
                            }
                        }

                        // F(i,j) = max( F(i,j), C(i,j) )
                        if (i == j)
                        { // F(i,i) = F(i,i)
                            f_coterminus_score = nucj_score;
                            update(bestF_j, index, f_coterminus_score, -1, -1, -1, -1, NONE);
                            is_F_OCT = true;
                            is_F_STEP = false;
                        }
                        else
                        {
                            auto it = bestSC_j[l + 1].find(index);
                            if (it != bestSC_j[l + 1].end())
                            {
                                f_coterminus_score = it->second.score;

                                // if (i == 0 || j == seq_length - 1) // 確定在外層 所以給分
                                // {
                                f_coterminus_score += -MFE_scalar * v_score_external_paired(EnrBASE(nuci), EnrBASE(nucj));
                                // }

                                if (f_coterminus_score > bestF_j[index].score)
                                {
                                    bool res = update(bestF_j, index, f_coterminus_score, -1, -1, -1, -1, MANNER_F_C);
                                    is_F_OCT = true;
                                    is_F_STEP = true;
                                    if (res)
                                    {
                                        bestF_j[index].outer_pairs.push_back(index);
                                    }
                                }
                            }
                        }

                        if (is_F_STEP)
                        {
                            auto it = bestF_j.find(index);
                            if (it != bestF_j.end())
                            {
                                it->second.is_F_STEP = true;
                            }
                        }

                        // Forward Update F(k,j),  0 <= k <= i-1
                        if (is_F_OCT && (i != 0))
                        { // F(k,j) = F(k,i-1)+ F(i,j)    ( STEP + OCT ) Forward Update
                            std::unordered_map<int, State<T>> &bestF_i_1 = bestF[i - 1];
                            for (auto &state_k : bestF_i_1)
                            { // F(1,i-1), F(2,i-1),...,F(k,i-1),...,F(i-1,i-1)
                                if (!state_k.second.is_F_STEP)
                                    continue;
                                std::tie(k, nuck, nuci_1) = GetIndexTuple(state_k.first);

                                if (!IsLegal(con_seq, nuci_1, nuci, i - 1))
                                    continue;

                                index_k_j = GetIndex(k, nuck, nucj);
                                newscore = state_k.second.score + f_coterminus_score;

                                bool res = update(bestF_j, index_k_j, newscore, state_k.first, index, -1, i - 1, MANNER_F_F_F);
                                if (res)
                                {
                                    bestF_j[index_k_j].outer_pairs = state_k.second.outer_pairs;
                                    if (i != j)
                                        bestF_j[index_k_j].outer_pairs.push_back(index);
                                }
                            }
                        }
                    }
                } // nuci end
            } // nucj end

            if (i == 0 && F_speedup)
            {
                if (i == j) // F(0, 0)
                {
                    for (auto &nucj : nucj_list)
                    {
                        index = GetIndex(i, nucj, nucj);
                        update(bestF_j, index, 0.0, -1, -1, -1, -1, NONE);
                    }
                }
                else
                {
                    std::unordered_map<int, State<T>> &bestF_j_1 = bestF[j - 1];

                    // F(0, j) = F(0, j-1) + F(j, j)
                    for (auto &nucj : nucj_list)
                    {
                        if (isLastNuc(j))
                        {
                            nucj_score = getCAI(ami_seq[j / 3], nucj, is_DN);
                        }
                        else
                        {
                            nucj_score = 0.0;
                        }

                        for (auto &itemF_j_1 : bestF_j_1)
                        {
                            State<T> stateF = itemF_j_1.second;
                            int i_, nuci_, nucj_1;
                            std::tie(i_, nuci_, nucj_1) = GetIndexTuple(itemF_j_1.first);

                            if (IsLegal(con_seq, nucj_1, nucj, j - 1))
                            {
                                newscore = stateF.score + nucj_score;
                                index = GetIndex(i_, nuci_, nucj);

                                int index_j_j = GetIndex(j, nucj, nucj);
                                update(bestF_j, index, newscore, itemF_j_1.first, index_j_j, -1, j - 1, MANNER_F_F_F);
                            }
                        }
                    }

                    // F(0, j) = C(0, j)
                    if (j - i > HAIRPIN_GAP)
                    {
                        int k = 0;
                        for (auto &itemC_k : bestSC_j[j - k + 1])
                        {
                            State<T> stateC = itemC_k.second;
                            int k, nuck, nucj_;
                            std::tie(k, nuck, nucj_) = GetIndexTuple(itemC_k.first);

                            index = itemC_k.first;

                            newscore = stateC.score;
                            newscore += -MFE_scalar * v_score_external_paired(EnrBASE(nuck), EnrBASE(nucj_));
                            update(bestF_j, index, newscore, -1, -1, -1, -1, MANNER_F_C);
                        }
                    }

                    // F(0, j) = F(0, k-1) + C(k, j)
                    for (int k = 1; j - k > HAIRPIN_GAP; k++)
                    {
                        std::unordered_map<int, State<T>> &bestF_k_1 = bestF[k - 1];

                        for (auto &itemC_k : bestSC_j[j - k + 1])
                        {
                            State<T> stateC = itemC_k.second;
                            int k_, nuck, nucj_;
                            std::tie(k_, nuck, nucj_) = GetIndexTuple(itemC_k.first);

                            for (auto &itemF_k_1 : bestF_k_1)
                            {
                                State<T> stateF = itemF_k_1.second;
                                int i_, nuci_, nuck_1;
                                std::tie(i_, nuci_, nuck_1) = GetIndexTuple(itemF_k_1.first);

                                if (!IsLegal(con_seq, nuck_1, nuck, k - 1))
                                    continue;

                                index = GetIndex(i_, nuci_, nucj_);
                                newscore = stateF.score + stateC.score;
                                newscore += -MFE_scalar * v_score_external_paired(EnrBASE(nuck), EnrBASE(nucj_));

                                update(bestF_j, index, newscore, itemF_k_1.first, itemC_k.first, -1, k - 1, MANNER_F_F_C);
                            }
                        }
                    }
                }
            }

        } // i end
    } // j end

    // std::cout << internal_time << std::endl;

    // std::cout << multi_time << std::endl;
    // std::cout << F_time << std::endl;

    return;
}


template<typename T>
std::vector<pair<T,T>> sparse_result_output(AllTables<T> &alltables,string& rna_seq, vector<int>& con_seq, string& ami_seq, string& output_txt, string& output_csv,  string& PID, int protein_length, string& ires_id, int ires_len, double TimeSpend, bool show_score,bool is_DN){//return {CAI,MFE}
    std::ofstream outputfile(output_txt,std::ios::app);
    std::ofstream outputcsv(output_csv,std::ios::app);

    std::string rna_solution;
    std::string structure_solution;

    outputfile << "===================================================================================================================" << std::endl;

    outputfile << "Protein ID: " << PID << std::endl;
    outputfile << "Protein length: " << protein_length << std::endl;
    outputfile << "Lambda: " << std::fixed << std::setprecision(3) << lambda << std::endl;
    if(ires_id != ""){
        outputfile << "IRES ID: " << ires_id << std::endl;
        outputfile << "IRES length: " << ires_len << std::endl;
        outputfile << "Crosspair Penalty: " << std::fixed << std::setprecision(3) << crosspair_penalty << std::endl;
    }

    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Linear RNA Reconstruction:" << std::endl;
    std::cout << "=== LinearRNA (Sparse) ===" << std::endl;

    T maxscore = BackTrack_sparse<T>(rna_solution, structure_solution, rna_seq, alltables, con_seq);
    T weighted_cai_score = GetCAIScore<T>(rna_solution,is_DN);
    double cai_value = GetUnweghtedCAIScore(rna_solution);
    cai_value = exp(cai_value/double(ami_seq.size()));
    //output score and results
    std::cout << "Coding sequence and its secondary structure:" << std::endl;
    outputfile << "Coding sequence and its secondary structure:" << std::endl;
    std::cout << rna_solution << std::endl;
    outputfile << rna_solution << std::endl;
    std::cout << structure_solution << std::endl;
    outputfile << structure_solution << std::endl;
    
    if (show_score){
        std::cerr << "Score: " << maxscore << std::endl;
        outputfile << "Score: " << maxscore << std::endl;
    }

    double mfe_value ;
    if(is_DN){
        mfe_value = -double(double(maxscore - weighted_cai_score)/(100.0*lambda));//*100
        std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
        outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    } else {
        if(lambda){
            mfe_value = -double(double(maxscore - weighted_cai_score)/(100.0));//*100
            std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
            outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;  
        } else {//only consider MFE
            mfe_value = -(maxscore/(100.0));//*100
            std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
            outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;  
        }
    }
    std::cout << "CAI: " << cai_value << std::endl;
    outputfile << "CAI: " <<  cai_value << std::endl;

    std::cout << "Total runtime: " << TimeSpend << " s" << std::endl;
    outputfile <<  "Total runtime: " << TimeSpend << " s" << std::endl;

    check_ami_solution(ami_seq, rna_solution, ires_len);

    std::vector<pair<T,T>> result_container;
    result_container.push_back({round_up(cai_value),round_up(mfe_value)});
    // output csv
    outputcsv<<std::endl<<PID<<","<<std::to_string(protein_length)<<","<<std::to_string(lambda)<<","<<std::to_string(mfe_value)<<","<<std::to_string(cai_value)<<","<<TimeSpend;

    outputcsv.close();
    outputfile.close();
    
    return result_container;

}
