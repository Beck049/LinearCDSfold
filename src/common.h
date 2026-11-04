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

using namespace std;
using namespace oligodesign;

inline bool isLastNuc(int pos){
    if(pos%3 == 2)
        return true;
    return false;
}

inline double getCAI(char Amino_label, int nuc, bool is_DN){
    if(is_DN) {
        return CodonSetCAIMap_DN[Amino_label][nuc];
    } else {
        return CodonSetCAIMap[Amino_label][nuc];
    }
}

inline bool compare(double x, double y, double epsilon = 1.0E-5) {
    if(fabs(x - y) < epsilon) return true;
    return false;
}

double calTime(const struct timespec start, const struct timespec end) {
    struct timespec temp;
    double time_used = 0.0;

    if ((end.tv_nsec - start.tv_nsec) < 0) {
        temp.tv_sec = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }

    time_used = temp.tv_sec + (double)temp.tv_nsec / 1000000000.0;
    return time_used;
}


inline int isCrossTalk(int i, int j, int ires_len){
    if(i < ires_len && j >= ires_len){
        return 1;
    } else if(i >= ires_len && j < ires_len) {
        return 1;
    } else {
        return 0;
    }
}

inline int getCrossPairNum(const std::string& structure, int ires_len) {
    if (ires_len <= 0) return 0;
    int cross_count = 0;
    std::vector<int> stack;
    for (int j = 0; j < (int)structure.length(); ++j) {
        if (structure[j] == '(') {
            stack.push_back(j);
        } else if (structure[j] == ')') {
            if (!stack.empty()) {
                int i = stack.back();
                stack.pop_back();
                if (isCrossTalk(i, j, ires_len)) {
                    cross_count++;
                }
            }
        }
    }
    return cross_count;
}

/** *************************************************************************************************************** */
// Init Special HP
template <typename T>
void initialize_Special_HP(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq, int ires_len, bool is_DN, bool is_sparse = false) {
    ACMatcher<char, 256> matcher;

    for (const auto& sp_loop : sp_loops_vec) {
        string s = sp_loop;
        matcher.AddString(s.c_str(), s.length());
    }

    if (!matcher.MakeTree()) {
        std::cout << "Failed to build the Aho-Corasick tree" << std::endl;
        return;
    }

    // 收集所有可能出現 special hairpin 的 <結尾位置, hairpin字串集合>
    unordered_map<int, unordered_set<string>> global_sp;
    for (int i = ires_len; i < (int)ami_seq.length() - 3; i++) {
        const auto& Amino_1 = Amino_to_nucs[ami_seq[i]];
        const auto& Amino_2 = Amino_to_nucs[ami_seq[i + 1]];
        const auto& Amino_3 = Amino_to_nucs[ami_seq[i + 2]];
        const auto& Amino_4 = Amino_to_nucs[ami_seq[i + 3]];

        for (const auto& A1 : Amino_1) {
            for (const auto& A2 : Amino_2) {
                for (const auto& A3 : Amino_3) {
                    for (const auto& A4 : Amino_4) {
                        std::string sub_seq = A1 + A2 + A3 + A4;
                        unordered_map<int, int> sp;
                        matcher.SearchPositions(sub_seq.c_str(), i * 3, 0, sub_seq.length() - 1, sp);
                        for (const auto& ss : sp) {
                            string special_hp = sub_seq.substr(ss.first - i * 3 - ss.second + 1, ss.second);
                            global_sp[ss.first].insert(special_hp);
                        }
                    }
                }
            }
        }
    }

    // 遍歷所有匹配到的 hairpin 並計算分數與更新狀態
    for (const auto& g : global_sp) {
        for (const auto& h : g.second) {
            vector<int> last_nuc, first_nuc;
            int last_idx = g.first;
            int first_idx = last_idx - h.length() + 1;
            unordered_map<int, T> first_CAI;

            // 1. 結尾核苷酸判定
            if (isLastNuc(last_idx)) {
                string codon_ = h.substr(h.length() - 3, 3);
                last_nuc.push_back(getLastExtendedNuc(re_Amino_label[codon_][0], codon_));
            } else {
                const auto& nuc_list  = Base_table.at(rna_seq[last_idx]);
                const auto& nuc_list2 = Base_table.at(rna_seq[last_idx - 1]);
                for (auto nuc : nuc_list) {
                    if (reBASE(nuc) == h[h.length() - 1]) {
                        for (auto nuc2 : nuc_list2) {
                            if (reBASE(nuc2) == h[h.length() - 2] && IsLegal(con_seq, nuc2, nuc, last_idx - 1)) {
                                last_nuc.push_back(nuc);
                            }
                        }
                    }
                }
            }

            // 2. 開頭核苷酸判定
            if (first_idx % 3 == 1) {
                const auto& nuc_list  = Base_table.at(rna_seq[first_idx]);
                const auto& nuc_list2 = Base_table.at(rna_seq[first_idx + 1]);
                for (auto nuc : nuc_list) {
                    if (reBASE(nuc) == h[0]) {
                        for (auto nuc2 : nuc_list2) {
                            if (reBASE(nuc2) == h[1] && IsLegal(con_seq, nuc, nuc2, first_idx)) {
                                first_nuc.push_back(nuc);
                                first_CAI[nuc] = getCAI(ami_seq[(first_idx + 1) / 3], nuc2, is_DN);
                            }
                        }
                    }
                }
            } else {
                const auto& nuc_list = Base_table.at(rna_seq[first_idx]);
                for (auto nuc : nuc_list) {
                    if (reBASE(nuc) == h[0]) {
                        first_nuc.push_back(nuc);
                    }
                }
            }

            // 3. 基礎能量計算
            T base_score = is_DN ? (-lambda * sp_loops[h]) : (-1 * sp_loops[h]);
            double param = is_DN ? (1 - lambda) : lambda;
            int cross_num = is_sparse ? 0 : isCrossTalk(first_idx, last_idx, ires_len);

            // 4. 更新 DP 表格
            for (auto nuci : first_nuc) {
                for (auto nucj : last_nuc) {
                    T newscore = base_score;
                    for (int x = first_idx; x <= last_idx; x++) {
                        if (isLastNuc(x) && param) {
                            if (x - first_idx == 0) {
                                newscore += getCAI(ami_seq[x / 3], nuci, is_DN);
                            } else if (x - first_idx == 1) {
                                newscore += first_CAI[nuci];
                            } else {
                                string codon_ = h.substr(x - first_idx - 2, 3);
                                int nucx = getLastExtendedNuc(re_Amino_label[codon_][0], codon_);
                                newscore += getCAI(ami_seq[x / 3], nucx, is_DN);
                            }
                        }
                    }

                    int index = GetIndex(first_idx, nuci, nucj);

                    if (is_sparse) {
                        int l = last_idx - first_idx;
                        auto& target_map = alltables.bestSC[last_idx][l + 1];
                        update(target_map, index, newscore, -1, -1, -1, -1, MANNER_S_HP);
                        if (target_map[index].score == newscore) {
                            sp_backtrack[index] = h;
                        }
                    } else {
                        newscore -= cross_num * crosspair_penalty;
                        auto& target_map = alltables.bestC[last_idx];
                        update(target_map, index, newscore, -1, MANNER_S_HP, cross_num);
                        if (target_map[index].score == newscore) {
                            sp_backtrack[index] = h;
                        }
                    }
                }
            }
        }
    }
}

/** *************************************************************************************************************** */
// Backtracking
template <typename T>
void BacktrackStructure(string& rna_solution, string& structure_solution, int seq_length, int maxindex, State<T> state, AllTables<T>& alltables, int start_j = -1){
    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    stack<tuple<int, int, State<T>>> stk;
    int j_start = (start_j >= 0) ? start_j : (seq_length - 1);
    stk.push(make_tuple(j_start, maxindex, state));
    while ( !stk.empty() ) {

        tuple<int, int, State<T>> top = stk.top();
        int j = get<0>(top);
        int cindex = get<1>(top);
        state = get<2>(top);
        int index1 = state.index_1;
        int index2 = state.index_2;
        int index3 = state.index_3;
        int i, nuci, nuci_pair, nucj, nucx, len;
        Manner MANNER = state.MANNER;
        State<T> next_state, next_state1, next_state2, next_state3;
        // cout<<mannerToString(MANNER)<<endl;
        if(MANNER == MANNER_C_StoCS)//make_tuple(i, nuci, nuci_pair,nucj, nucx, len);
            std::tie(i, nuci, nuci_pair, nucj, nucx, len) = GetIndexTupleCS(cindex);
        else{
            std::tie(i, nuci, nucj) = GetIndexTuple(cindex);
        }

        
        stk.pop();
        int k, nuck, nuck_pair, nucl, posi, s, nucs, nucs_end;
        int struct_score, p, q, q1, p_1, l1, l2;
        switch (MANNER){
            
            case NONE:
                break;
            case MANNER_NONEtoN:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                break;
            case MANNER_NONEtoF:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                break;
            case MANNER_NONEtoS:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                break;
            case MANNER_N_EtoN:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestN[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_NtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state = bestN[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_S_HP:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                for(int x = i + 1; x < j ; x ++ ){
                    rna_solution[x]=sp_backtrack[cindex][x - i];
                    structure_solution[x]='.';
                }
                break;
            case MANNER_S_EtoS:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                std::tie(k, nuck, nucl) = GetIndexTuple(index1);
                next_state = bestS[j-1][(j-1)-k+1][index1];//consider len of single
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_CtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state = bestC[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_S_CtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestS[k-1][(k-1)-(i+1)+1][index1];//consider len of single
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestC[j-1][index2];
                stk.push(make_tuple(j-1 ,index2, next_state2));
                break;
            case MANNER_CStoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state1 = bestCS[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state1));
                break;
            case MANNER_C_StoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                std::tie(q1, nucs, nuck) = GetIndexTuple(index2);
                q = q1-1;
                next_state1 = bestC[q][index1];
                stk.push(make_tuple(q ,index1, next_state1));
                next_state2 = bestS[j-1][(j-1)-(q1)+1][index2];//consider len of single
                stk.push(make_tuple(j-1 ,index2, next_state2));
                break;
            case MANNER_S_CStoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state2 = bestCS[j-1][index2];
                stk.push(make_tuple(j-1 ,index2, next_state2));
                std::tie(k, nuck, nuck_pair,nucl, nucx,posi) = GetIndexTupleCS(index2);//from CS to get S
                next_state1 = bestS[k-1][(k-1)-(i+1)+1][index1];//consider len of single
                stk.push(make_tuple(k-1 ,index1, next_state1));
                break;
            case MANNER_S_C_StoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                std::tie(q1,nuck,nucs) = GetIndexTuple(index3);//nuc is not important
                l1 = (j-1) - q1 + 1;
                next_state3 = bestS[j-1][l1][index3];//consider len of single
                stk.push(make_tuple(j-1 ,index3, next_state3));
                q = q1-1;

                next_state2 = bestC[q][index2];
                stk.push(make_tuple(q ,index2, next_state2));
                std::tie(p,nuck,nucs) = GetIndexTuple(index2);//nuc is not important
                p_1 = p-1;

                l2 = p_1 - (i+1)+ 1;
                next_state1 = bestS[p_1][l2][index1];//consider len of single
                stk.push(make_tuple(p_1 ,index1, next_state1));
                break;
            case MANNER_C_StoCS:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestC[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestS[j][len][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            case MANNER_CtoM1:
                next_state = bestC[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_M1_CtoM2:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestM1[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestC[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            case MANNER_M1_EtoM1:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestM1[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_M2toM1:
                next_state = bestM2[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_M2toMulti:
                next_state = bestM2[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_S_M2toMulti:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                std::tie(s, nucs, nucs_end) = GetIndexTuple(index1);
                next_state1 = bestS[k-1][(k-1)-s+1][index1];//consider len of single
                // next_state1 = bestN[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestM2[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            case MANNER_Multi_EtoMulti:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestMulti[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_MultitoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state = bestMulti[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));

                // if(debug){
                //     struct_score += - v_score_multi(GET_ACGU_NUM_V(rna_solution[i]),GET_ACGU_NUM_V(rna_solution[j]));
                //     printf("Multi loop ( %d, %d) %c%c : %.2f\n", i+1, j+1, rna_solution[i], rna_solution[j], struct_score / -100.0);
                // }
                break;
            case MANNER_F_EtoF:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestF[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_CtoF:
                next_state = bestC[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_F_CtoF:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestF[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestC[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            default:
                cout << "Something error" << endl;
            
        }
        //cout<<mannerToString(MANNER)<<", "<<state.score/lambda<<", "<<i+1<<", "<<j+1<<endl;
        
    }
}

template <typename T>
void BacktrackStructure_sparse(string &rna_solution, string &structure_solution, int seq_length, int maxindex, State<T> state, AllTables<T> &alltables, int start_j = -1) {
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

    while (!stk.empty()) {
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

        switch (MANNER) {
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
            if (p - i > 1) {
                std::tie(i_, nuci_, nucj_) = GetIndexTuple(index1);
                index_ = GetIndex1(nuci_, nucj_);
                next_state1 = bestSS[p - 1][p - 1 - (i + 1) + 1][index_];
                stk.push(make_tuple(p - 1, index1, next_state1));
            }

            next_state2 = bestSC[q][q - p + 1][index2];
            stk.push(make_tuple(q, index2, next_state2));

            if (j - q > 1) {
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
            for (int x = i + 1; x < j; x++) {
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
T BackTrack(string& rna_solution, string& structure_solution, string& rna_seq, AllTables<T>& alltables, vector<int>& con_seq, bool is_sparse){
    // cout<<"start back tracking"<<endl;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;

    int seq_length = rna_seq.size();
    rna_solution.assign(seq_length, '.');
    structure_solution.assign(seq_length, '.');

    vector<int> first_nuclist = Base_table.at(rna_seq[0]);//i list
    vector<int> last_nuclist = Base_table.at(rna_seq[seq_length-1]);//j list

    // 找到最佳結構分數的 index & MANNER
    // go through 所有 nuc[0] 和 nuc[n-1] 的組合，找到 bestF 中分數最高的組合
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

    if(is_sparse){
        BacktrackStructure_sparse(rna_solution, structure_solution, seq_length, maxindex, state, alltables);
    } else {
        BacktrackStructure(rna_solution, structure_solution, seq_length, maxindex, state, alltables);
    }

    return maxscore;
}

/** *************************************************************************************************************** */
// Result Output
template<typename T>
std::vector<pair<T,T>> result_output(AllTables<T> &alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq, string& output_txt, string& output_csv, string& PID, int protein_length, string& ires_id, int ires_len, double TimeSpend, bool show_score, bool is_DN, bool is_sparse) { // return {CAI, MFE}
    std::ofstream outputfile(output_txt, std::ios::app);
    std::ofstream outputcsv(output_csv, std::ios::app);

    std::string rna_solution;
    std::string structure_solution;

    outputfile << "===================================================================================================================" << std::endl;

    outputfile << "Protein ID: " << PID << std::endl;
    outputfile << "Protein length: " << protein_length << std::endl;
    outputfile << "Lambda: " << std::fixed << std::setprecision(3) << lambda << std::endl;

    if (ires_id != "") {
        outputfile << "IRES ID: " << ires_id << std::endl;
        outputfile << "IRES length: " << ires_len << std::endl;
        outputfile << "Crosspair Penalty: " << std::fixed << std::setprecision(3) << crosspair_penalty << std::endl;
    }
    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Linear RNA Reconstruction:" << std::endl;
    std::cout << "=== LinearRNA";
    std::cout << (is_sparse?"(sparse)":"");
    std::cout << " ===" << std::endl;

    T maxscore = BackTrack<T>(rna_solution, structure_solution, rna_seq, alltables, con_seq, is_sparse);
    T weighted_cai_score = GetCAIScore<T>(rna_solution, is_DN);
    double cai_value = GetUnweghtedCAIScore(rna_solution);
    cai_value = exp(cai_value / double(ami_seq.size()));

    // 1. 計算 Crosstalk 數量與懲罰分數
    int cross_pair_count = getCrossPairNum(structure_solution, ires_len);
    double cross_penalty_score = cross_pair_count * crosspair_penalty;

    // 輸出序列與結構
    std::cout << "Coding sequence and its secondary structure:" << std::endl;
    outputfile << "Coding sequence and its secondary structure:" << std::endl;
    std::cout << rna_solution << std::endl;
    outputfile << rna_solution << std::endl;
    std::cout << structure_solution << std::endl;
    outputfile << structure_solution << std::endl;

    // 2. 還原真實的 MFE（排除 CAI 分數與 Crosspair 懲罰）
    double mfe_value;
    double mfe_score_component; // 記錄對應在 maxscore 中的 MFE 成分

    if (is_DN) {
        // maxscore = mfe_score_component + weighted_cai_score - cross_penalty_score
        // 其中 mfe_score_component = -100.0 * lambda * mfe_value
        mfe_score_component = double(maxscore) - double(weighted_cai_score) + cross_penalty_score;
        mfe_value = -mfe_score_component / ((double)LDCONST * lambda);
    } else {
        if (lambda) {
            mfe_score_component = double(maxscore) - double(weighted_cai_score) + cross_penalty_score;
            mfe_value = -mfe_score_component / (double)LDCONST;
        } else { // only consider MFE
            mfe_score_component = double(maxscore) + cross_penalty_score;
            mfe_value = -mfe_score_component / (double)LDCONST;
        }
    }

    // 3. 輸出分數組成分解
    if (show_score) {
        std::cerr << "--- Score Breakdown ---" << std::endl;
        std::cerr << "Total Score: " << maxscore << std::endl;
        std::cerr << "  - MFE Score: " << mfe_score_component << std::endl;
        std::cerr << "  - CAI Score (weighted): " << weighted_cai_score << std::endl;
        std::cerr << "  - Crosspair Penalty Score: -" << cross_penalty_score 
                  << " (count: " << cross_pair_count << ", penalty: " << crosspair_penalty << ")" << std::endl;

        outputfile << "--- Score Breakdown ---" << std::endl;
        outputfile << "Total Score: " << maxscore << std::endl;
        outputfile << "  - MFE Score: " << mfe_score_component << std::endl;
        outputfile << "  - CAI Score (weighted): " << weighted_cai_score << std::endl;
        outputfile << "  - Crosspair Penalty Score: -" << cross_penalty_score 
                   << " (count: " << cross_pair_count << ", penalty: " << crosspair_penalty << ")" << std::endl;
    }

    // 4. 輸出個別生理/物理指標
    std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;

    std::cout << "CAI: " << cai_value << std::endl;
    outputfile << "CAI: " << cai_value << std::endl;

    if(ires_len != 0) {
        std::cout << "Crosspair count: " << cross_pair_count << std::endl;
        outputfile << "Crosspair count: " << cross_pair_count << std::endl;
    }

    std::cout << "Total runtime: " << TimeSpend << " s" << std::endl;
    outputfile << "Total runtime: " << TimeSpend << " s" << std::endl;

    check_ami_solution(ami_seq, rna_solution, ires_len);

    std::vector<pair<T,T>> result_container;
    result_container.push_back({round_up(cai_value), round_up(mfe_value)});

    // 5. CSV 加入 crosspair 欄位輸出
    // 格式：PID, protein_length, lambda, MFE, CAI, Crosspair_count, Crosspair_penalty, Runtime
    outputcsv << std::endl << PID << "," << std::to_string(protein_length) << "," << std::to_string(lambda) << "," ;
    if(ires_len != 0) {
        outputcsv << std::to_string(crosspair_penalty) << "," ;
    }
    outputcsv << std::to_string(mfe_value) << "," << std::to_string(cai_value) << "," ;
     if(ires_len != 0) {
        outputcsv << std::to_string(cross_pair_count) << "," ;
    }
    outputcsv << TimeSpend;

    outputcsv.close();
    outputfile.close();

    return result_container;
}
