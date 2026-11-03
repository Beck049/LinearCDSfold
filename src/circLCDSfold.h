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
#include "LCDSfold.h"

using namespace std;
using namespace oligodesign;

template<typename T>
struct CircularCandidate {
    State<T> state;
    T score = VALUE_MIN<T>();
    int index = -1;
    int i = -1;
    int j = -1;
    std::string codon_i = "";
    std::string codon_j = "";
    Manner manner;
    int cross_pair_num = 0;
};

/**
 * @brief 給定 structure[i, j]，計算 acid[i] 或 acid[j] 處的額外 CAI 分數與最佳 Codon
 * @return std::pair<std::string, double> 回傳 {codon_j, extra_cai}
 * 
 * 當 type == 'i' 時，僅回傳 codon 前 2 char
 * 僅 type == 'j' 時，extra_cai 才有效
 */
std::pair<std::string, double> compute_extra_cai(
    int nth, 
    int nucn, 
    const std::string& ami_seq, 
    bool is_DN,
    char type // ['i', 'j']
) {
    std::string codon_n = "";
    double extra_cai = INT_MIN;
    
    char target_ami = ami_seq[nth / 3];
    const std::vector<std::string>& all_nuc = Amino_to_nucs.at(target_ami);
    int pos_in_codon = nth % 3;

    std::map<char, std::vector<int>> pos_order_map = {
        {'j', {2, 0, 1}},
        {'i', {0, 2, 1}}
    };

    if (pos_in_codon == pos_order_map.at(type)[0]) {
        // do nothing
    } else if (pos_in_codon == pos_order_map.at(type)[1]) {
        for (const auto& codon : all_nuc) {
            if (codon[pos_in_codon] == reBASE(nucn)) {
                double score = getScoreByCodon(std::string(1, target_ami), codon, is_DN);
                if (score > extra_cai) {
                    extra_cai = score;
                    codon_n = codon;
                }
            }
        }
    } else if (pos_in_codon == pos_order_map.at(type)[2]) {
        bool is_LRS = (target_ami == 'L' || target_ami == 'R' || target_ami == 'S');

        if (!is_LRS) {
            for (const auto& codon : all_nuc) {
                if (codon[pos_in_codon] == reBASE(nucn)) {
                    double score = getScoreByCodon(std::string(1, target_ami), codon, is_DN);
                    if (score > extra_cai) {
                        extra_cai = score;
                        codon_n = codon;
                    }
                }
            }
        } else {
            char nucn_1 = '\0';
            if (target_ami == 'R') {
                if (nucn == 5)       nucn_1 = 'A';
                else if (nucn == 6)  nucn_1 = 'C';
                else if (nucn == 13) nucn_1 = 'C';
            } else if (target_ami == 'L') {
                if (nucn == 7)       nucn_1 = 'U';
                else if (nucn == 8)  nucn_1 = 'C';
                else if (nucn == 14)  nucn_1 = 'C';
            } else if (target_ami == 'S') {
                if (nucn == 3)       nucn_1 = 'A';
                else if (nucn == 2)  nucn_1 = 'U';
            }

            if(type == 'j') {
                for (const auto& codon : all_nuc) {
                    if (codon[pos_in_codon] == reBASE(nucn) && codon[pos_in_codon-1] == nucn_1) {
                        double score = getScoreByCodon(std::string(1, target_ami), codon, is_DN);
                        if (score > extra_cai) {
                            extra_cai = score;
                            codon_n = codon;
                        }
                    }
                }
            } else if (type == 'i') {
                if (nucn_1 != '\0') {
                    codon_n = std::string(1, nucn_1) + reBASE(nucn) + '.';
                } else {
                    // 若找不到前一個配對字元，補預設值或原始字元
                    codon_n = std::string(".") + reBASE(nucn) + '.';
                }
            }
        }
    }

    return {codon_n, extra_cai};
}

/**
 * @brief 計算 Circular RNA 在連接跨端區間形成的 Interior loop / External connection MFE
 * 
 * @tparam T 數值型別 (float, double)
 * @param alltables DP 狀態記錄表
 * @param rna_seq 原始輸入 RNA 序列
 * @param ami_seq 胺基酸序列
 * @param BestAminoMap 各胺基酸最佳 CAI 的 Codon 對應表
 * @param state 當前候選的狀態 (State<T>)
 * @param id 當前狀態的編號 index
 * @param ith 當前結構左端索引
 * @param jth 當前結構右端索引
 * @param codon_i 補齊 i 所屬 Codon 的字串
 * @param codon_j 補齊 j 所屬 Codon 的字串
 * @param seq_length 序列總長度
 * @param debug 是否印出偵錯資訊
 * @return double 計算得到的額外能量補償 (extra_mfe)
 */
template<typename T>
double calcInternalMFE(
    AllTables<T>& alltables,
    const std::string& rna_seq,
    const std::string& ami_seq,
    const std::unordered_map<char, std::pair<std::string, double>>& BestAminoMap,
    const State<T>& state,
    int id,
    int ith,
    int jth,
    const std::string& codon_i,
    const std::string& codon_j,
    int seq_length,
    bool debug = true
) {
    // 1. 回溯取得當前局部解的二級結構與序列
    std::string tmp_rna(seq_length, '.');
    std::string tmp_structure(seq_length, '.');
    BacktrackStructure(tmp_rna, tmp_structure, seq_length, id, state, alltables, jth);

    // 2. 找出在 [ith, jth] 範圍內 exterior loop 的所有 base pairs
    std::vector<std::pair<int, int>> exterior_pairs;
    std::stack<int> bp_stack;

    for (int pos = ith; pos <= jth; ++pos) {
        if (tmp_structure[pos] == '(') {
            bp_stack.push(pos);
        } else if (tmp_structure[pos] == ')') {
            if (!bp_stack.empty()) {
                int open_idx = bp_stack.top();
                bp_stack.pop();
                if (bp_stack.empty()) {
                    exterior_pairs.push_back({open_idx, pos});
                }
            }
        }
    }

    // 3. 補齊 ith 與 jth 所屬 Codon
    if(codon_i != ""){
        for(int k = 0; k < 3; ++k){
            if(codon_i[k] == '.') continue;
            tmp_rna[((ith/3)*3)+k] = codon_i[k];
        }
    }
    if(codon_j != ""){
        for(int k = 0; k < 3; ++k){
            tmp_rna[((jth/3)*3)+k] = codon_j[k];
        }
    }

    // 4. 補齊 ith-1 與 jth+1 的 codon
    const int ami_i = (ith / 3) - 1;
    if( ami_i >= 0 ){
        const std::string&  best_codon_i_1 = BestAminoMap.at(ami_seq[ami_i]).first;
        for (int k = 0; k < 3; ++k) {
            tmp_rna[ami_i * 3 + k] = best_codon_i_1[k];
        }
    }
    
    const int ami_j = (jth / 3) + 1;
    if( ami_j < ami_seq.size() ){
        const std::string& best_codon_j_1 = BestAminoMap.at(ami_seq[ami_j]).first;
        for (int k = 0; k < 3; ++k) {
            tmp_rna[ami_j * 3 + k] = best_codon_j_1[k];
        }
    } else if (ami_i != 0) {
        const std::string& best_codon_j_1 = BestAminoMap.at(ami_seq[0]).first;
        for (int k = 0; k < 3; ++k) {
            tmp_rna[k] = best_codon_j_1[k];
        }
    }

    // 5. 偵錯資訊輸出
    if (debug) {
        std::cout << "\n[Debug Interior Structure " << ith << " " << jth << "]: \n"
                  << tmp_rna << "\n" << tmp_structure << std::endl;
        std::cout << codon_i << ", " << codon_j << std::endl;
        std::cout << "Exterior pairs in range [" << ith << ", " << jth << "]: ";
        for (const auto& [u, v] : exterior_pairs) {
            std::cout << "(" << u + 1 << ", " << v + 1 << ") [" << tmp_rna[u] << "-" << tmp_rna[v] << "] ";
        }
        std::cout << std::endl;
    }

    // 6. 計算能量 (MFE)
    // 預設為閉合 Multiloop / Interior loop 扣除開環代價
    int i = exterior_pairs[0].first;
    int j = exterior_pairs[0].second;
    int p = exterior_pairs[1].first;
    int q = exterior_pairs[1].second;
    // std::cout << "(" << i << ", " << j << ")" << "NUC: " << tmp_rna[i] << ", " << tmp_rna[i+1] << ", " << tmp_rna[j-1] << ", " << tmp_rna[j] << std::endl;
    // std::cout << "(" << p << ", " << q << ")" << "NUC: " << tmp_rna[p-1] << ", " << tmp_rna[p] << ", " << tmp_rna[q] << ", " << tmp_rna[(q+1) % seq_length] << std::endl;

    auto to_nuc_id = [&](char c) -> int {
        char upper_c = std::toupper(static_cast<unsigned char>(c));
        switch (upper_c) {
            case 'A': return 1;
            case 'C': return 2;
            case 'G': return 3;
            case 'U': case 'T': return 4;
            default:
                std::cerr << "[Warning] Unrecognized nucleotide: '" << c << "'\n";
                return 0;
        }
    };

    double extra_mfe = -v_score_single(
        i, j, p, q,
        to_nuc_id(tmp_rna[i]), to_nuc_id(tmp_rna[i+1]), to_nuc_id(tmp_rna[j-1]), to_nuc_id(tmp_rna[j]),
        to_nuc_id(tmp_rna[p-1]), to_nuc_id(tmp_rna[p]), to_nuc_id(tmp_rna[q]), to_nuc_id(tmp_rna[(q+1) % seq_length])
    );
    return extra_mfe;
}


// 直接環化（用 vienna 環化方式）
template<typename T>
std::vector<pair<T,T>> circular_compose_v1(AllTables<T> &alltables,string& rna_seq, vector<int>& con_seq, string& ami_seq, string& output_txt, string& output_csv,  string& PID, int protein_length, double TimeSpend, bool show_score, bool is_beam, bool is_DN){//return {CAI,MFE}
    std::ofstream outputfile(output_txt,std::ios::app);
    std::ofstream outputcsv(output_csv,std::ios::app);

    double c_TimeSpend = 0;
    auto c_start_time = std::chrono::high_resolution_clock::now();

    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Circular RNA Reconstruction:" << std::endl;
    std::cout << "=== CircularRNA (v1) ===" << std::endl;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    int seq_length = rna_seq.size();

    // 用 CodonSetCAIMap or CodonSetCAIMap_DN
    std::unordered_map<char, std::pair<std::string, double>> BestAminoMap;

    for (const auto& [label, codon_list] : Amino_to_nucs) {
        // std::cout << "Label: " << label << std::endl;
        for (const auto& codon : codon_list) {
            // std::cout << "Codon: " << codon << " (" << unweighted_CAIMap.at(codon) << ")" << std::endl;
            if (unweighted_CAIMap.at(codon) == 0){
                BestAminoMap.insert({label, {codon, 0}});
            }
        }
    }

    // std::cout << "========= Best Amino Acid Codon Map =========\n";
    // for (const auto& [acid, info_pair] : BestAminoMap) {
    //     std::cout << "Amino: [" << acid << "] -> Best Codon: " << info_pair.first  << ", Score: " << info_pair.second << "\n";
    // }

    // 考慮 除了 dangling end 以外的 circular structure
    CircularCandidate<T> best_struct;

    /**
     * 只考慮頭尾皆 paired 的情況
     * 因為 C 屬於 M1，且 M2 也屬於 M1，所以僅考慮 M1
     * case 1 - CtoM1 (需考慮 Hairpin_Gap)
     * case 2 - M2toM1 (可能為 Multi)
     * 
     * 向前後延伸，並計算 CAI & MFE
     */
    // go through bestM1 中所有元素，找到分數最高的 circular structure
    for (size_t jth = 0; jth < bestM1.size(); ++jth) {
        for (auto& [id, state] : bestM1[jth]) {
            // 為 從 i 到 j 結構
            int ith, nuci, nucj, last_pair_pos;
            Manner manner = state.MANNER;
            std::tie(ith, nuci, nucj) = GetIndexTuple(id);
            // cout << "(" << nuci << "/" << nucj << ")";
            double cur_score = state.score;
            std::string codon_i = "";
            std::string codon_j = "";
            // 子結構為 C, M2 且 hairpin < MAXLOOP 才考慮
            //  nuc  為 ith ~ jth   // amino 為 ith/3 ~ jth/3
            int hairpin_len = seq_length - jth + ith;
            if ((manner == MANNER_CtoM1 ||  manner == MANNER_M2toM1) && hairpin_len < MAXLOOP) {
                // cout << "nuc (" << i << "~" << j << "); acid (" << i/3 << "~" << j/3 << ")" << std::endl;

                // CAI (只考慮 acid[j] 的 cai)
                auto [codon_j, extra_cai] = compute_extra_cai(jth, nucj, ami_seq, is_DN, 'j');

                // 只考慮 acid[i] 的 codon
                auto [codon_i, extra_cai_dump] = compute_extra_cai(ith, nuci, ami_seq, is_DN, 'i');

                // MFE
                double extra_mfe = 0;
                if (manner == MANNER_CtoM1) {   // hairpin
                    extra_mfe += hairpin37[min(hairpin_len, 30)];
                    extra_mfe += (hairpin_len > 30) ? (int)(lxc37*log((hairpin_len)/30.)) : 0;
                } else if (manner == MANNER_M2toM1) {
                    if(state.index_3 > 2) { // Multi
                        extra_mfe -= ML_intern37 + ML_closing37;
                    } else {    // Interior
                        extra_mfe += calcInternalMFE<T>(alltables, rna_seq, ami_seq, BestAminoMap, state,
                                                            id, ith, jth, codon_i, codon_j, seq_length, false );
                    }
                }
                if (is_DN) extra_mfe *= lambda;

                // compare
                if (best_struct.score < cur_score + extra_mfe - ((extra_cai == INT_MIN)?0:extra_cai) && (hairpin_len > HAIRPIN_GAP || hairpin_len == 0)) {
                    best_struct.score = cur_score + extra_mfe - ((extra_cai == INT_MIN)?0:extra_cai);
                    best_struct.index = id;
                    best_struct.i = ith;
                    best_struct.j = jth;
                    best_struct.codon_i = codon_i;
                    best_struct.codon_j = codon_j;
                    best_struct.state = state;
                    best_struct.manner = manner;
                }
            }
        }
    }

    if (best_struct.index == -1) {
        std::cout << "No circular structure found." << std::endl;
        outputfile << "No circular structure found." << std::endl;
        return {};
    }

    // 建構 Output
    double c_final_score = best_struct.score;
    std::string c_rna_solution(seq_length, 'X');
    std::string c_structure_solution(seq_length, '.');

    std::cout << "Coding sequence and its secondary structure (" << best_struct.i << "~" << best_struct.j
              << ", " << mannerToString(best_struct.manner) << "):" << std::endl;
    outputfile << "Coding sequence and its secondary structure :" << std::endl;

    // 取得 rna sequence & structure
    BacktrackStructure(c_rna_solution, c_structure_solution, seq_length, best_struct.index, best_struct.state, alltables, best_struct.j);
    // 補齊 i, j
    if(best_struct.codon_i != ""){
        for(int k = 0; k < 3; ++k){
            if(best_struct.codon_i[k] == '.') continue;
            c_rna_solution[((best_struct.i/3)*3)+k] = best_struct.codon_i[k];
        }
    }
    if(best_struct.codon_j != ""){
        for(int k = 0; k < 3; ++k){
            c_rna_solution[((best_struct.j/3)*3)+k] = best_struct.codon_j[k];
        }
    }
    // 補齊前後 sequence
    for(int ami = 0; ami < best_struct.i/3; ++ami){
        for(int k = 0; k < 3; ++k){
            c_rna_solution[((ami)*3)+k] = BestAminoMap.at(ami_seq[ami]).first[k];
        }
    }
    for(int ami = (best_struct.j/3)+1; ami < ami_seq.size(); ++ami){
        for(int k = 0; k < 3; ++k){
            c_rna_solution[((ami)*3)+k] = BestAminoMap.at(ami_seq[ami]).first[k];
        }
    }

    for (int pos = 0; pos < seq_length; ++pos) {
        if (c_rna_solution[pos] == '.') {
            c_rna_solution[pos] = rna_seq[pos];
        }
    }

    std::cout << c_rna_solution << std::endl;
    outputfile << c_rna_solution << std::endl;
    std::cout << c_structure_solution << std::endl;
    outputfile << c_structure_solution << std::endl;

    // MFE & CAI
    T weighted_cai_score = GetCAIScore<T>(c_rna_solution,is_DN);
    double cai_value = GetUnweghtedCAIScore(c_rna_solution);
    cai_value = exp(cai_value/double(ami_seq.size()));

    double mfe_value ;
    if(is_DN){
        mfe_value = -double(double(c_final_score - weighted_cai_score)/(100.0*lambda));//*100
        std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
        outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    } else {
        if(lambda){
            mfe_value = -double(double(c_final_score - weighted_cai_score)/(100.0));//*100
            std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
            outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;  
        } else {//only consider MFE
            mfe_value = -(c_final_score/(100.0));//*100
            std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
            outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;  
        }
    }
    std::cout << "CAI: " << std::round(cai_value * 1000.0) / 1000.0 << std::endl;
    outputfile << "CAI: " <<  std::round(cai_value * 1000.0) / 1000.0 << std::endl;

    auto c_end_time = std::chrono::high_resolution_clock::now();
    c_TimeSpend = std::chrono::duration_cast<std::chrono::duration<double>>(c_end_time - c_start_time).count();
    std::cout << "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;
    outputfile <<  "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;

    std::vector<pair<T,T>> result_container;
    result_container.push_back({round_up(cai_value),round_up(mfe_value)});

    // output csv
    outputcsv << "," << mfe_value << "," << cai_value << "," << c_TimeSpend;

    // PrintCalCAIMap(is_DN);
    // PrintUnweightedCAIMap();

    return result_container;
}

// 直接環化（用 mfold 環化方式）
template<typename T>
std::vector<pair<T,T>> circular_compose_m1(AllTables<T> &alltables,string& rna_seq, vector<int>& con_seq, string& ami_seq, string& output_txt, string& output_csv,  string& PID, int protein_length, double TimeSpend, bool show_score, bool is_beam, bool is_DN){//return {CAI,MFE}
    std::ofstream outputfile(output_txt,std::ios::app);
    std::ofstream outputcsv(output_csv,std::ios::app);

    double c_TimeSpend = 0;
    auto c_start_time = std::chrono::high_resolution_clock::now();

    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Circular RNA Reconstruction:" << std::endl;
    std::cout << "=== CircularRNA (m1) ===" << std::endl;

    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;

    int seq_length = rna_seq.size() / 2;

    // 考慮 除了 dangling end 以外的 circular structure
    CircularCandidate<T> best_struct;

    int best_i = -1;
    int best_j = -1;
    int best_id1 = -1;
    int best_id2 = -1;
    State<T> best_state1;
    State<T> best_state2;
    T best_score = -std::numeric_limits<T>::infinity();

    // 1. 遍歷 j from 0 to seq_length
    // 尋找 C(i, j) 與 C(j, i + seq_length) 的最優組合
    for (int j = 0; j < seq_length; ++j) {
        if (j >= bestC.size() || bestC[j].empty()) continue;

        for (const auto& [id1, state1] : bestC[j]) {
            int i, nuci1, nucj1;
            std::tie(i, nuci1, nucj1) = GetIndexTuple(id1);

            int second_end = i + seq_length;
            if (second_end >= bestC.size() || bestC[second_end].empty()) continue;

            // 尋找右側閉合區塊 C(j, i + seq_length)
            for (const auto& [id2, state2] : bestC[second_end]) {
                int start2, nucj2, nuci2;
                std::tie(start2, nucj2, nuci2) = GetIndexTuple(id2);

                // 確認起點與第一個區塊的終點 j 相接且核苷酸相容
                if (start2 == j && nucj1 == nucj2 && nuci1 == nuci2) {
                    T combined_score = state1.score + state2.score;

                    if (combined_score > best_score) { // 若分數越高越好用 >，反之用 <
                        best_score = combined_score;
                        best_i = i;
                        best_j = j;
                        best_id1 = id1;
                        best_id2 = id2;
                        best_state1 = state1;
                        best_state2 = state2;
                    }
                }
            }
        }
    }

    if (best_i == -1) {
        std::cout << "No circular structure found." << std::endl;
        outputfile << "No circular structure found." << std::endl;
        return {};
    }

    // 建構 Output
    // 2. 回溯雙倍長度的序列與結構 (長度 2 * seq_length)
    std::string full_rna(rna_seq.size(), '.');
    std::string full_struct(rna_seq.size(), '.');

    // 分別 Backtrack 兩個配對區塊 C(i, j) 與 C(j, i + seq_length)
    BacktrackStructure(full_rna, full_struct, rna_seq.size(), best_id1, best_state1, alltables, best_j);
    BacktrackStructure(full_rna, full_struct, rna_seq.size(), best_id2, best_state2, alltables, best_i + seq_length);

    // 未配對部分填入 rna_seq 原始核苷酸
    for (size_t pos = 0; pos < full_rna.size(); ++pos) {
        if (full_rna[pos] == '.' || full_rna[pos] == 'X') {
            full_rna[pos] = rna_seq[pos];
        }
    }

    // 3. 截取 [best_i, best_i + seq_length - 1] 作為單倍長度的 CDS 與二級結構輸出
    //    旋轉回到原本 0 ~ seq_length-1 (1 ~ seq_length) 的順序, 原區間為 [best_i, best_i + seq_length - 1]，對應的原座標為 k % seq_length
    std::string c_rna_solution(seq_length, ' ');
    std::string c_structure_solution(seq_length, '.');

    // 收集區間 [best_i, best_i + seq_length - 1] 內產生的配對清單 (u, v)
    std::vector<int> pair_table(2 * seq_length, -1);
    std::stack<int> st;
    for (int k = best_i; k < best_i + seq_length; ++k) {
        if (full_struct[k] == '(') {
            st.push(k);
        } else if (full_struct[k] == ')') {
            if (!st.empty()) {
                int left = st.top();
                st.pop();
                pair_table[left] = k;
                pair_table[k] = left;
            }
        }
    }

    // 將序列填回 0 ~ seq_length-1 (1-based 的 1 ~ seq_length)
    for (int offset = 0; offset < seq_length; ++offset) {
        int original_pos = (best_i + offset) % seq_length;
        c_rna_solution[original_pos] = full_rna[best_i + offset];
    }

    // 依據原始座標重構配對括號 (保持 '(' 在前、')' 在後)
    for (int k = best_i; k < best_i + seq_length; ++k) {
        if (pair_table[k] != -1 && k < pair_table[k]) {
            int pos1 = k % seq_length;
            int pos2 = pair_table[k] % seq_length;

            if (pos1 > pos2) std::swap(pos1, pos2);

            c_structure_solution[pos1] = '(';
            c_structure_solution[pos2] = ')';
        }
    }

    std::cout << "Optimal circular pair split at: i=" << best_i << ", j=" << best_j 
              << ", i+L=" << (best_i + seq_length) << std::endl;
    std::cout << "Coding sequence and its secondary structure (range [" << best_i << ", " << best_i + seq_length - 1 << "]):" << std::endl;
    outputfile << "Coding sequence and its secondary structure :" << std::endl;

    std::cout << c_rna_solution << std::endl;
    outputfile << c_rna_solution << std::endl;
    std::cout << c_structure_solution << std::endl;
    outputfile << c_structure_solution << std::endl;

    // 4. 計算 MFE 與 CAI
    T weighted_cai_score = GetCAIScore<T>(c_rna_solution, is_DN);
    double cai_value = GetUnweghtedCAIScore(c_rna_solution);
    cai_value = exp(cai_value / double(ami_seq.size() / 2)); // 此處 ami_seq 原為 2 倍長度，故除以原始長度

    double mfe_value;
    if (is_DN) {
        mfe_value = -double(double(best_score - weighted_cai_score) / (100.0 * lambda));
    } else {
        if (lambda) {
            mfe_value = -double(double(best_score - weighted_cai_score) / 100.0);
        } else {
            mfe_value = -(best_score / 100.0);
        }
    }

    std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    std::cout << "CAI: " << std::round(cai_value * 1000.0) / 1000.0 << std::endl;
    outputfile << "CAI: " << std::round(cai_value * 1000.0) / 1000.0 << std::endl;

    auto c_end_time = std::chrono::high_resolution_clock::now();
    c_TimeSpend = std::chrono::duration_cast<std::chrono::duration<double>>(c_end_time - c_start_time).count();
    std::cout << "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;
    outputfile << "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;

    std::vector<pair<T,T>> result_container;
    result_container.push_back({round_up(cai_value), round_up(mfe_value)});
    outputcsv << "," << mfe_value << "," << cai_value << "," << c_TimeSpend;

    return result_container;
}

// 包含 cross_pair_penalty 的環化重建函式
template<typename T>
std::vector<pair<T,T>> circular_compose_i1(
    AllTables<T> &alltables,
    string& rna_seq,
    vector<int>& con_seq,
    string& ami_seq,
    string& output_txt,
    string& output_csv,
    string& PID,
    int protein_length,
    int ires_len,
    double TimeSpend,
    bool show_score,
    bool is_beam,
    bool is_DN
) { // return {CAI, MFE}
    std::ofstream outputfile(output_txt, std::ios::app);
    std::ofstream outputcsv(output_csv, std::ios::app);

    double c_TimeSpend = 0;
    auto c_start_time = std::chrono::high_resolution_clock::now();

    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Circular RNA Reconstruction:" << std::endl;
    std::cout << "=== CircularRNA (i1) ===" << std::endl;

    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    int seq_length = rna_seq.size();

    // 建立 BestAminoMap (含安全 fallback，防止特定胺基酸未配對成功)
    std::unordered_map<char, std::pair<std::string, double>> BestAminoMap;
    for (const auto& [label, codon_list] : Amino_to_nucs) {
        std::string best_codon = "";
        double min_cai = std::numeric_limits<double>::max();
        for (const auto& codon : codon_list) {
            auto it = unweighted_CAIMap.find(codon);
            if (it != unweighted_CAIMap.end()) {
                if (it->second == 0) {
                    best_codon = codon;
                    min_cai = 0;
                    break;
                } else if (it->second < min_cai) {
                    min_cai = it->second;
                    best_codon = codon;
                }
            }
        }
        if (!best_codon.empty()) {
            BestAminoMap[label] = {best_codon, min_cai};
        } else if (!codon_list.empty()) {
            BestAminoMap[label] = {codon_list[0], 0.0};
        }
    }

    // 安全取得 Best Codon 的 Lambda 函式，避免 out_of_range
    auto safe_get_best_codon = [&](char aa) -> std::string {
        auto it = BestAminoMap.find(aa);
        if (it != BestAminoMap.end()) {
            return it->second.first;
        }
        auto it_amino = Amino_to_nucs.find(aa);
        if (it_amino != Amino_to_nucs.end() && !it_amino->second.empty()) {
            return it_amino->second[0];
        }
        if (aa == '*') return "UAA";
        return "NNN";
    };

    // 安全包裝 compute_extra_cai
    auto safe_compute_extra_cai = [&](int nth, int nucn, const std::string& seq, bool is_dn, char type) -> std::pair<std::string, double> {
        if (nth / 3 >= (int)seq.size()) {
            return {"", INT_MIN};
        }
        char aa = seq[nth / 3];
        if (Amino_to_nucs.find(aa) == Amino_to_nucs.end()) {
            return {"", INT_MIN};
        }
        return compute_extra_cai(nth, nucn, seq, is_dn, type);
    };

    CircularCandidate<T> best_struct;

    // 遍歷 bestM1 中所有狀態，搜尋分數最高的 circular 候選結構
    for (size_t jth = 0; jth < bestM1.size(); ++jth) {
        for (auto& [id, state] : bestM1[jth]) {
            int ith, nuci, nucj;
            Manner manner = state.MANNER;
            std::tie(ith, nuci, nucj) = GetIndexTuple(id);

            double cur_score = state.score;
            std::string codon_i = "";
            std::string codon_j = "";
            int hairpin_len = seq_length - jth + ith;

            if ((manner == MANNER_CtoM1 || manner == MANNER_M2toM1) && hairpin_len < MAXLOOP) {
                auto [c_j, extra_cai] = safe_compute_extra_cai(jth, nucj, ami_seq, is_DN, 'j');
                auto [c_i, extra_cai_dump] = safe_compute_extra_cai(ith, nuci, ami_seq, is_DN, 'i');
                codon_j = c_j;
                codon_i = c_i;

                double extra_mfe = 0;
                int cross_delta = 0;

                if (manner == MANNER_CtoM1) { // 跨端 Hairpin
                    extra_mfe += hairpin37[min(hairpin_len, 30)];
                    extra_mfe += (hairpin_len > 30) ? (int)(lxc37 * log((hairpin_len) / 30.)) : 0;
                    cross_delta = isCrossTalk(ith, jth, ires_len);
                } else if (manner == MANNER_M2toM1) {
                    if (state.index_3 > 2) { // Multi loop
                        extra_mfe -= ML_intern37 + ML_closing37;
                    } else { // Interior loop
                        extra_mfe += calcInternalMFE<T>(alltables, rna_seq, ami_seq, BestAminoMap, state,
                                                        id, ith, jth, codon_i, codon_j, seq_length, false);
                    }
                }

                if (is_DN) extra_mfe *= lambda;

                // 總分計算包含：原始 DP 分數 + 跨端閉合 MFE - CAI 補償 - 閉合產生的 Crosstalk 處罰
                double candidate_score = cur_score + extra_mfe 
                                       - ((extra_cai == INT_MIN) ? 0 : extra_cai)
                                       - (cross_delta * crosspair_penalty);

                if (best_struct.score < candidate_score && (hairpin_len > HAIRPIN_GAP || hairpin_len == 0)) {
                    best_struct.score = candidate_score;
                    best_struct.index = id;
                    best_struct.i = ith;
                    best_struct.j = jth;
                    best_struct.codon_i = codon_i;
                    best_struct.codon_j = codon_j;
                    best_struct.state = state;
                    best_struct.manner = manner;
                    best_struct.cross_pair_num = state.cross_pair_num + cross_delta;
                }
            }
        }
    }

    if (best_struct.index == -1) {
        std::cout << "No circular structure found." << std::endl;
        outputfile << "No circular structure found." << std::endl;
        return {};
    }

    double c_final_score = best_struct.score;
    std::string c_rna_solution(seq_length, 'X');
    std::string c_structure_solution(seq_length, '.');

    std::cout << "Coding sequence and its secondary structure (" << best_struct.i << "~" << best_struct.j
              << ", " << mannerToString(best_struct.manner) << "):" << std::endl;
    outputfile << "Coding sequence and its secondary structure :" << std::endl;

    BacktrackStructure(c_rna_solution, c_structure_solution, seq_length, best_struct.index, best_struct.state, alltables, best_struct.j);

    // 補齊 ith 與 jth 所屬 Codon
    if (best_struct.codon_i != "") {
        for (int k = 0; k < 3; ++k) {
            if (best_struct.codon_i[k] == '.') continue;
            int idx = ((best_struct.i / 3) * 3) + k;
            if (idx < seq_length) c_rna_solution[idx] = best_struct.codon_i[k];
        }
    }
    if (best_struct.codon_j != "") {
        for (int k = 0; k < 3; ++k) {
            int idx = ((best_struct.j / 3) * 3) + k;
            if (idx < seq_length) c_rna_solution[idx] = best_struct.codon_j[k];
        }
    }

    // 補齊前後序列
    for (int ami = 0; ami < best_struct.i / 3 && ami < (int)ami_seq.size(); ++ami) {
        std::string codon = safe_get_best_codon(ami_seq[ami]);
        for (int k = 0; k < 3; ++k) {
            int idx = ami * 3 + k;
            if (idx < seq_length) c_rna_solution[idx] = codon[k];
        }
    }
    for (int ami = (best_struct.j / 3) + 1; ami < (int)ami_seq.size(); ++ami) {
        std::string codon = safe_get_best_codon(ami_seq[ami]);
        for (int k = 0; k < 3; ++k) {
            int idx = ami * 3 + k;
            if (idx < seq_length) c_rna_solution[idx] = codon[k];
        }
    }

    // 未配對部分回填原始 rna_seq
    for (int pos = 0; pos < seq_length; ++pos) {
        if (c_rna_solution[pos] == '.' || c_rna_solution[pos] == 'X') {
            c_rna_solution[pos] = rna_seq[pos];
        }
    }

    std::cout << c_rna_solution << std::endl;
    outputfile << c_rna_solution << std::endl;
    std::cout << c_structure_solution << std::endl;
    outputfile << c_structure_solution << std::endl;

    // 1. 計算 Crosstalk 配對數量與懲罰分數
    int cross_pair_count = getCrossPairNum(c_structure_solution, ires_len);
    double cross_penalty_score = cross_pair_count * crosspair_penalty;

    // 2. MFE & CAI 分數分解
    T weighted_cai_score = GetCAIScore<T>(c_rna_solution, is_DN);
    double cai_value = GetUnweghtedCAIScore(c_rna_solution);
    cai_value = exp(cai_value / double(ami_seq.size()));

    double mfe_value;
    double mfe_score_component;

    if (is_DN) {
        mfe_score_component = double(c_final_score) - double(weighted_cai_score) + cross_penalty_score;
        mfe_value = -mfe_score_component / (100.0 * lambda);
    } else {
        if (lambda) {
            mfe_score_component = double(c_final_score) - double(weighted_cai_score) + cross_penalty_score;
            mfe_value = -mfe_score_component / 100.0;
        } else {
            mfe_score_component = double(c_final_score) + cross_penalty_score;
            mfe_value = -mfe_score_component / 100.0;
        }
    }

    // 3. 輸出分數細項 Breakdown
    if (show_score) {
        std::cerr << "--- Score Breakdown ---" << std::endl;
        std::cerr << "Total Score: " << c_final_score << std::endl;
        std::cerr << "  - MFE Score: " << mfe_score_component << std::endl;
        std::cerr << "  - CAI Score (weighted): " << weighted_cai_score << std::endl;
        std::cerr << "  - Crosspair Penalty Score: -" << cross_penalty_score 
                  << " (count: " << cross_pair_count << ", penalty: " << crosspair_penalty << ")" << std::endl;

        outputfile << "--- Score Breakdown ---" << std::endl;
        outputfile << "Total Score: " << c_final_score << std::endl;
        outputfile << "  - MFE Score: " << mfe_score_component << std::endl;
        outputfile << "  - CAI Score (weighted): " << weighted_cai_score << std::endl;
        outputfile << "  - Crosspair Penalty Score: -" << cross_penalty_score 
                   << " (count: " << cross_pair_count << ", penalty: " << crosspair_penalty << ")" << std::endl;
    }

    std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    std::cout << "CAI: " << std::round(cai_value * 1000.0) / 1000.0 << std::endl;
    outputfile << "CAI: " << std::round(cai_value * 1000.0) / 1000.0 << std::endl;
    std::cout << "Crosspair count: " << cross_pair_count << std::endl;
    outputfile << "Crosspair count: " << cross_pair_count << std::endl;

    auto c_end_time = std::chrono::high_resolution_clock::now();
    c_TimeSpend = std::chrono::duration_cast<std::chrono::duration<double>>(c_end_time - c_start_time).count();
    std::cout << "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;
    outputfile << "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;

    std::vector<pair<T,T>> result_container;
    result_container.push_back({round_up(cai_value), round_up(mfe_value)});

    // 格式保持與專案其他 CSV 輸出一致
    outputcsv << "," << mfe_value << "," << cai_value << "," << cross_pair_count << "," << c_TimeSpend;

    outputcsv.close();
    outputfile.close();

    return result_container;
}