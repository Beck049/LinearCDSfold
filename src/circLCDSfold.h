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
std::vector<pair<T,T>> circular_compose_v1(AllTables<T> &alltables,string& rna_seq, vector<int>& con_seq, string& ami_seq, string& output_txt, string& output_csv,  string& PID, int protein_length, double TimeSpend, bool show_score, bool is_beam, bool is_DN){//return {CAI,MFE}
    std::ofstream outputfile(output_txt,std::ios::app);
    std::ofstream outputcsv(output_csv,std::ios::app);

    std::string rna_solution;
    std::string structure_solution;

    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Circular RNA Reconstruction:" << std::endl;
    std::cout << "=== CircularDesign ===" << std::endl;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::vector<std::unordered_map<int, State<T>>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    int seq_length = rna_seq.size();
    rna_solution.assign(seq_length, '.');
    structure_solution.assign(seq_length, '.');

    // 用 CodonSetCAIMap or CodonSetCAIMap_DN
    std::unordered_map<char, std::pair<std::string, double>> BestAminoMap;

    
    for (const auto& [acid_label, sub_map] : (is_DN ? CodonSetCAIMap_DN : CodonSetCAIMap)) {
        for (const auto& [last_nuc, new_score] : sub_map) {
            std::string acid_codon = reCodon(std::string(1, acid_label), last_nuc);
            if (BestAminoMap.count(acid_label)) {
                if (BestAminoMap.at(acid_label).second < new_score) {
                    BestAminoMap.at(acid_label) = {acid_codon, new_score};
                }
            } else {
                BestAminoMap.insert({acid_label, {acid_codon, new_score}});
            }
        }
    }
    // std::cout << "========= Best Amino Acid Codon Map =========\n";
    // for (const auto& [acid, info_pair] : BestAminoMap) {
    //     std::cout << "Amino: [" << acid << "] -> Best Codon: " << info_pair.first 
    //             << ", Score: " << info_pair.second << "\n";
    // }

    // 考慮 除了 dangling end 以外的 circular structure
    State<T> bestState;
    T max_score = VALUE_MIN<T>();
    int max_index = -1;
    int best_i, best_j;
    std::string best_codon_i = "";
    std::string best_codon_j = "";
    Manner max_manner;

    /**
     * 在 Exact 模式下，
     * 只考慮頭尾皆 paired 的情況
     * 因為 C 屬於 M1，且 M2 也屬於 M1，所以僅考慮 M1
     * case 1 - CtoM1 (需考慮 Hairpin_Gap)
     * case 2 - M2toM1
     * 
     * 向前後延伸，並計算 CAI & MFE
     */
    if(!is_beam) {
        // go through bestM1 中所有元素，找到分數最高的 circular structure
        for (size_t j = 0; j < bestM1.size(); ++j) {
            for (auto& [id, state] : bestM1[j]) {
                // 為 從 i 到 j 結構
                int i, nuci, nucj, last_pair_pos;
                Manner manner;
                std::tie(i, nuci, nucj) = GetIndexTuple(id);
                manner = state.MANNER;
                double cur_score = state.score;
                std::string codon_j = "";
                // 子結構為 C, M2
                //  nuc  為 i ~ j
                // amino 為 i/3 ~ j/3
                if (manner == MANNER_CtoM1 ||  manner == MANNER_M2toM1) {
                    int hairpin_len = seq_length - j + i;
                    // cout << "nuc (" << i << "~" << j << "); acid (" << i/3 << "~" << j/3 << ")" << std::endl;

                    // MFE
                    double extra_mfe = 0;
                    if (manner == MANNER_CtoM1) { // hairpin
                        extra_mfe += hairpin37[min(hairpin_len, 30)];
                        extra_mfe += (hairpin_len > 30) ? (int)(lxc37*log((hairpin_len)/30.)) : 0;
                    } else if (manner == MANNER_M2toM1) { // multi
                        extra_mfe -= ML_intern37 + ML_closing37;
                    }
                    if (is_DN) extra_mfe *= lambda;

                    // CAI (只考慮 acid[j] 的 cai)
                    double extra_cai = -9999;
                    if(j%3 == 2) {
                        extra_cai = 0;
                    } else if(j%3 == 0) {
                        // std::cout << ami_seq[j/3] << "(" << reBASE(nucj) << "," << j%3 << ")" << " (";
                        std::vector<std::string> all_nuc = Amino_to_nucs.at(ami_seq[j/3]);
                        for (size_t i = 0; i < all_nuc.size(); ++i) {
                            if(all_nuc[i][j%3] == reBASE(nucj)){
                                if (extra_cai < getScoreByCodon(std::string(1, ami_seq[j/3]), all_nuc[i], is_DN)){
                                    extra_cai = getScoreByCodon(std::string(1, ami_seq[j/3]), all_nuc[i], is_DN);
                                    codon_j = all_nuc[i];
                                }
                            }
                        }
                        // std::cout << codon_j << " " << getScoreByCodon(std::string(1, ami_seq[j/3]), codon_j, is_DN) << ")" << std::endl;
                    } else if(j%3 == 1) {
                        // 加 backtrack
                        State<T> c_state;
                        int nucj_1;
                        if(manner == MANNER_CtoM1){
                            c_state = bestC[j][state.index_1];
                        } else if (manner == MANNER_M2toM1) {
                            State<T> m2_state = bestM2[j][state.index_1];
                            c_state = bestC[j][m2_state.index_2];
                        }
                        Manner preManner = c_state.MANNER;
                        int i, nuci;
                        switch (preManner){
                            case MANNER_NtoC:
                                std::tie(i, nuci, nucj_1) = GetIndexTuple(c_state.index_1);
                                break;
                            case MANNER_CtoC:
                                std::tie(i, nuci, nucj_1) = GetIndexTuple(c_state.index_1);
                                break;
                            case MANNER_S_CtoC:
                                std::tie(i, nuci, nucj_1) = GetIndexTuple(c_state.index_2);
                                break;
                            case MANNER_CStoC:
                                std::tie(i, nuci, nucj_1) = GetIndexTuple(c_state.index_1);
                                break;
                            case MANNER_C_StoC:
                                std::tie(i, nuci, nucj_1) = GetIndexTuple(c_state.index_2);
                                break;
                            case MANNER_S_CStoC:
                                std::tie(i, nuci, nucj_1) = GetIndexTuple(c_state.index_2);
                                break;
                            case MANNER_S_C_StoC:
                                std::tie(i, nuci, nucj_1) = GetIndexTuple(c_state.index_3);
                                break;
                            case MANNER_MultitoC:
                                std::tie(i, nuci, nucj_1) = GetIndexTuple(c_state.index_1);
                                break;                            
                            default:
                                break;
                        }

                        // std::cout << ami_seq[j/3] << "(" << reBASE(nucj) << "," << j%3 << ")" << " (";
                        std::vector<std::string> all_nuc = Amino_to_nucs.at(ami_seq[j/3]);
                        for (size_t i = 0; i < all_nuc.size(); ++i) {
                            if(all_nuc[i][j%3] == reBASE(nucj) && all_nuc[i][(j%3)-1] == reBASE(nucj_1)){
                                if (extra_cai < getScoreByCodon(std::string(1, ami_seq[j/3]), all_nuc[i], is_DN)){
                                    extra_cai = getScoreByCodon(std::string(1, ami_seq[j/3]), all_nuc[i], is_DN);
                                    codon_j = all_nuc[i];
                                }
                            }
                        }
                        // std::cout << codon_j << " " << getScoreByCodon(std::string(1, ami_seq[j/3]), codon_j, is_DN) << ")" << std::endl;
                    }

                    // compare
                    if (max_score < cur_score + extra_mfe && (hairpin_len > HAIRPIN_GAP || hairpin_len == 0)) {
                        max_score = cur_score + extra_mfe;
                        max_index = id;
                        best_i = i;
                        best_j = j;
                        best_codon_j = codon_j;
                        bestState = state;
                        max_manner = manner;
                    }
                }
            }
        }
    } else {
    /** 
     * ToDo: Beam 模式
     */

    }

    if (max_index == -1) {
        std::cout << "No circular structure found." << std::endl;
        outputfile << "No circular structure found." << std::endl;
        return {};
    }

    std::cout << "Coding sequence and its secondary structure (" << best_i << "~" << best_j
              << ", " << mannerToString(max_manner) << "):" << std::endl;
    outputfile << "Coding sequence and its secondary structure :" << std::endl;

    // 取得 rna sequence & structure
    std::string c_rna_solution(seq_length, 'X');
    std::string c_structure_solution(seq_length, '.');
    BacktrackStructure(c_rna_solution, c_structure_solution, seq_length, max_index, bestState, alltables, best_j);
    // 補齊 i, j
    std::vector<std::string> all_nuci = Amino_to_nucs.at(ami_seq[best_i/3]);
    double max_nuci_score = -9999;
    for (size_t i = 0; i < all_nuci.size(); ++i) {
        if((best_i%3 == 1 && all_nuci[i][1] == c_rna_solution[best_i] && all_nuci[i][2] == c_rna_solution[best_i+1]) ||
            (best_i%3 == 2 && all_nuci[i][2] == c_rna_solution[best_i])){
            if (max_nuci_score < getScoreByCodon(std::string(1, ami_seq[best_i/3]), all_nuci[i], is_DN)){
                max_nuci_score = getScoreByCodon(std::string(1, ami_seq[best_i/3]), all_nuci[i], is_DN);
                best_codon_i = all_nuci[i];
            }
        }
    }
    if(best_codon_i != ""){
        for(int k = 0; k < 3; ++k){
            c_rna_solution[((best_i/3)*3)+k] = best_codon_i[k];
        }
    }
    if(best_codon_j != ""){
        for(int k = 0; k < 3; ++k){
            c_rna_solution[((best_j/3)*3)+k] = best_codon_j[k];
        }
    }
    // 補齊前後 sequence
    for(int ami = 0; ami < best_i/3; ++ami){
        for(int k = 0; k < 3; ++k){
            c_rna_solution[((ami)*3)+k] = BestAminoMap.at(ami_seq[ami]).first[k];
        }
    }
    for(int ami = (best_j/3)+1; ami < ami_seq.size(); ++ami){
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
    T weighted_cai_score = GetCAIScore<T>(rna_solution,is_DN);
    double cai_value = GetUnweghtedCAIScore(rna_solution);
    cai_value = exp(cai_value/double(ami_seq.size()));

    double mfe_value ;
    if(is_DN){
        mfe_value = -double(double(max_score - weighted_cai_score)/(100.0*lambda));//*100
        std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
        outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    } else {
        if(lambda){
            mfe_value = -double(double(max_score - weighted_cai_score)/(100.0));//*100
            std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
            outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;  
        } else {//only consider MFE
            mfe_value = -(max_score/(100.0));//*100
            std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
            outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;  
        }
    }
    std::cout << "CAI: " << cai_value << std::endl;
    outputfile << "CAI: " <<  cai_value << std::endl;

    

    return {};
}