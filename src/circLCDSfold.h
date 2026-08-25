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
                else if (nucn == 4)  nucn_1 = 'C';
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
                codon_n += nucn_1 + reBASE(nucn) + '.';
            }
        }
    }

    return {codon_n, extra_cai};
}

template<typename T>
std::vector<pair<T,T>> circular_compose_v1(AllTables<T> &alltables,string& rna_seq, vector<int>& con_seq, string& ami_seq, string& output_txt, string& output_csv,  string& PID, int protein_length, double TimeSpend, bool show_score, bool is_beam, bool is_DN){//return {CAI,MFE}
    std::ofstream outputfile(output_txt,std::ios::app);
    std::ofstream outputcsv(output_csv,std::ios::app);

    double c_TimeSpend = 0;
    auto c_start_time = std::chrono::high_resolution_clock::now();

    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Circular RNA Reconstruction:" << std::endl;
    std::cout << "=== CircularRNA ===" << std::endl;

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
                if (manner == MANNER_CtoM1) { // hairpin
                    extra_mfe += hairpin37[min(hairpin_len, 30)];
                    extra_mfe += (hairpin_len > 30) ? (int)(lxc37*log((hairpin_len)/30.)) : 0;
                } else if (manner == MANNER_M2toM1) { // multi
                    extra_mfe -= ML_intern37 + ML_closing37;
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