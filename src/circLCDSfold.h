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
std::pair<std::string, double> compute_extra_cai(int nth, int nucn, const std::string& ami_seq, bool is_DN, char type) {  // ['i', 'j']
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
    int id, int ith, int jth, 
    const std::string& codon_i, 
    const std::string& codon_j, 
    int seq_length, 
    bool is_sparse = false, 
    bool debug = true,
    int ires_len = 0
) {
    // 1. 回溯取得當前局部解的二級結構與序列
    std::string tmp_rna(seq_length, '.');
    std::string tmp_structure(seq_length, '.');
    if (is_sparse) {
        BacktrackStructure_sparse(tmp_rna, tmp_structure, seq_length, id, state, alltables, jth);
    } else {
        BacktrackStructure(tmp_rna, tmp_structure, seq_length, id, state, alltables, jth);
    }

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

    // 驗證 exterior_pairs 是否至少有 2 對 base pair，若不足無法構成封閉 internal loop
    if (exterior_pairs.size() < 2) {
        if (debug) {
            std::cerr << "[Warning] exterior_pairs size < 2 in calcInternalMFE (found: " 
                      << exterior_pairs.size() << ") for [" << ith << ", " << jth << "]\n";
        }
        return 0.0;
    }

    // 3. 補齊 ith 與 jth 所屬 Codon
    if (!codon_i.empty()) {
        int base_i = (ires_len > 0 && ith >= ires_len) 
                     ? (ires_len + ((ith - ires_len) / 3) * 3) 
                     : ((ith / 3) * 3);
        for (int k = 0; k < 3; ++k) {
            if (codon_i[k] == '.') continue;
            if (base_i + k < seq_length) tmp_rna[base_i + k] = codon_i[k];
        }
    }
    if (!codon_j.empty()) {
        int base_j = (ires_len > 0 && jth >= ires_len) 
                     ? (ires_len + ((jth - ires_len) / 3) * 3) 
                     : ((jth / 3) * 3);
        for (int k = 0; k < 3; ++k) {
            if (codon_j[k] == '.') continue;
            if (base_j + k < seq_length) tmp_rna[base_j + k] = codon_j[k];
        }
    }

    // 4. 補齊 ith-1 與 jth+1 的 codon（僅在 CDS 區域查表）
    int cds_offset_i = ith - ires_len;
    if (cds_offset_i >= 3) {
        const int ami_i = (cds_offset_i / 3) - 1;
        if (ami_i >= 0 && ami_i < (int)ami_seq.size() && BestAminoMap.count(ami_seq[ami_i])) {
            const std::string& best_codon_i_1 = BestAminoMap.at(ami_seq[ami_i]).first;
            int base_idx = ires_len + ami_i * 3;
            for (int k = 0; k < 3; ++k) {
                if (base_idx + k < seq_length) tmp_rna[base_idx + k] = best_codon_i_1[k];
            }
        }
    }

    int cds_offset_j = jth - ires_len;
    if (cds_offset_j >= 0) {
        const int ami_j = (cds_offset_j / 3) + 1;
        if (ami_j < (int)ami_seq.size() && BestAminoMap.count(ami_seq[ami_j])) {
            const std::string& best_codon_j_1 = BestAminoMap.at(ami_seq[ami_j]).first;
            int base_idx = ires_len + ami_j * 3;
            for (int k = 0; k < 3; ++k) {
                if (base_idx + k < seq_length) tmp_rna[base_idx + k] = best_codon_j_1[k];
            }
        }
    }

    // ★ 關鍵修復：將仍為 '.' 或 'X' 的未配對位置以原始 rna_seq 補齊
    for (int pos = 0; pos < seq_length; ++pos) {
        if (tmp_rna[pos] == '.' || tmp_rna[pos] == 'X') {
            if (pos < (int)rna_seq.size()) {
                tmp_rna[pos] = rna_seq[pos];
            } else {
                tmp_rna[pos] = 'A'; // 極端防護，避免殘留 '.'
            }
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
    int i = exterior_pairs[0].first;
    int j = exterior_pairs[0].second;
    int p = exterior_pairs[1].first;
    int q = exterior_pairs[1].second;

    auto to_nuc_id = [&](char c) -> int {
        char upper_c = std::toupper(static_cast<unsigned char>(c));
        switch (upper_c) {
            case 'A': return 1;
            case 'C': return 2;
            case 'G': return 3;
            case 'U': case 'T': return 4;
            default:
                if (debug) {
                    std::cerr << "[Warning] Unrecognized nucleotide: '" << c << "'\n";
                }
                return 1; // 遭遇未知字元時回傳預設鹼基(如 A=1)，避免能量表越界崩潰
        }
    };

    // 環狀邊界安全索引
    int idx_ip1 = (i + 1) % seq_length;
    int idx_jm1 = (j - 1 + seq_length) % seq_length;
    int idx_pm1 = (p - 1 + seq_length) % seq_length;
    int idx_qp1 = (q + 1) % seq_length;

    double extra_mfe = -v_score_single(
        i, j, p, q,
        to_nuc_id(tmp_rna[i]), 
        to_nuc_id(tmp_rna[idx_ip1]), 
        to_nuc_id(tmp_rna[idx_jm1]), 
        to_nuc_id(tmp_rna[j]),
        to_nuc_id(tmp_rna[idx_pm1]), 
        to_nuc_id(tmp_rna[p]), 
        to_nuc_id(tmp_rna[q]), 
        to_nuc_id(tmp_rna[idx_qp1])
    );

    return extra_mfe;
}

// 直接環化（用 vienna 環化方式）(可包含 IRES)
template<typename T>
std::vector<pair<T,T>> circular_compose_v2(AllTables<T> &alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq, string& ires_seq, string& output_txt, string& output_csv, string& PID, int protein_length, double TimeSpend, bool show_score, bool is_beam, bool is_DN) { // return {CAI, MFE}
    std::ofstream outputfile(output_txt, std::ios::app);
    std::ofstream outputcsv(output_csv, std::ios::app);

    double c_TimeSpend = 0;
    auto c_start_time = std::chrono::high_resolution_clock::now();

    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    int seq_length = rna_seq.size();
    int ires_len = ires_seq.length();

    // 將 ires_seq 中的 T/t 轉換為 U，並確保大小寫統一
    for (char &base : ires_seq) {
        if (base == 'T' || base == 't') {
            base = 'U';
        } else if (base == 'u') {
            base = 'U';
        }
    }

    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Circular RNA Reconstruction:" << std::endl;
    if (ires_len == 0) {
        std::cout << "=== CircularRNA (v1) ===" << std::endl;
    } else {
        std::cout << "=== CircularRNA (v1 - IRES) ===" << std::endl;
    }

    // 建立 BestAminoMap (含安全 fallback，防止特定胺基酸未找到對應密碼子)
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

    // 安全取得 Best Codon 的 Lambda 函式
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

    // 安全包裝 compute_extra_cai (考慮 ires_len 造成的座標平移)
    auto safe_compute_extra_cai = [&](int nth, int nucn, const std::string& seq, bool is_dn, char type) -> std::pair<std::string, double> {
        if (nth < ires_len) {
            return {"", INT_MIN}; // 位於 IRES 區段，不參與 CDS CAI 計算
        }
        int cds_nth = nth - ires_len;
        if (cds_nth / 3 >= (int)seq.size()) {
            return {"", INT_MIN};
        }
        char aa = seq[cds_nth / 3];
        if (Amino_to_nucs.find(aa) == Amino_to_nucs.end()) {
            return {"", INT_MIN};
        }
        return compute_extra_cai(cds_nth, nucn, seq, is_dn, type);
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
                    if (ires_len > 0) {
                        cross_delta = isCrossTalk(ith, jth, ires_len);
                    }
                } else if (manner == MANNER_M2toM1) {
                    if (state.index_3 > 2) { // Multi loop
                        extra_mfe -= ML_intern37 + ML_closing37;
                    } else { // Interior loop
                        extra_mfe += calcInternalMFE<T>(alltables, rna_seq, ami_seq, BestAminoMap, state,
                                                        id, ith, jth, codon_i, codon_j, seq_length, false, false, ires_len);
                    }
                }

                if (is_DN) extra_mfe *= lambda;

                // 總分計算：若 ires_len > 0 則納入 crosstalk 處罰
                double candidate_score = cur_score + extra_mfe - ((extra_cai == INT_MIN) ? 0 : extra_cai);
                if (ires_len > 0) {
                    candidate_score -= (cross_delta * crosspair_penalty);
                }

                if (best_struct.score < candidate_score && (hairpin_len > HAIRPIN_GAP || hairpin_len == 0)) {
                    best_struct.score = candidate_score;
                    best_struct.index = id;
                    best_struct.i = ith;
                    best_struct.j = jth;
                    best_struct.codon_i = codon_i;
                    best_struct.codon_j = codon_j;
                    best_struct.state = state;
                    best_struct.manner = manner;
                    if (ires_len > 0) {
                        best_struct.cross_pair_num = state.cross_pair_num + cross_delta;
                    }
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

    // 補齊 ith 與 jth 所屬 Codon (注意加上 ires_len 偏移量)
    if (best_struct.codon_i != "" && best_struct.i >= ires_len) {
        int cds_offset = best_struct.i - ires_len;
        int codon_start = ires_len + (cds_offset / 3) * 3;
        for (int k = 0; k < 3; ++k) {
            if (best_struct.codon_i[k] == '.') continue;
            int idx = codon_start + k;
            if (idx < seq_length) c_rna_solution[idx] = best_struct.codon_i[k];
        }
    }
    if (best_struct.codon_j != "" && best_struct.j >= ires_len) {
        int cds_offset = best_struct.j - ires_len;
        int codon_start = ires_len + (cds_offset / 3) * 3;
        for (int k = 0; k < 3; ++k) {
            if (best_struct.codon_j[k] == '.') continue;
            int idx = codon_start + k;
            if (idx < seq_length) c_rna_solution[idx] = best_struct.codon_j[k];
        }
    }

    // 1. 若有 IRES，前 ires_len 個鹼基直接填入 ires_seq
    if (ires_len > 0) {
        for (int pos = 0; pos < ires_len && pos < seq_length; ++pos) {
            c_rna_solution[pos] = ires_seq[pos];
        }
    }

    // 2. 補齊結構前段 CDS: 從 ires_len 到 best_struct.i 之前
    int cds_start_ami = 0;
    int cds_end_ami_front = (best_struct.i > ires_len) ? (best_struct.i - ires_len) / 3 : 0;

    for (int ami = cds_start_ami; ami < cds_end_ami_front && ami < (int)ami_seq.size(); ++ami) {
        std::string codon = safe_get_best_codon(ami_seq[ami]);
        int base_idx = ires_len + ami * 3;
        for (int k = 0; k < 3; ++k) {
            int idx = base_idx + k;
            if (idx < seq_length) c_rna_solution[idx] = codon[k];
        }
    }

    // 3. 補齊結構後段 CDS: 從 best_struct.j 之後到 CDS 結尾
    int cds_start_ami_back = 0;
    if (best_struct.j >= ires_len) {
        cds_start_ami_back = ((best_struct.j - ires_len) / 3) + 1;
    }

    for (int ami = cds_start_ami_back; ami < (int)ami_seq.size(); ++ami) {
        std::string codon = safe_get_best_codon(ami_seq[ami]);
        int base_idx = ires_len + ami * 3;
        for (int k = 0; k < 3; ++k) {
            int idx = base_idx + k;
            if (idx < seq_length) c_rna_solution[idx] = codon[k];
        }
    }

    // 4. 未配對部分回填原始 rna_seq (再次確保前 ires_len 鎖定為 ires_seq)
    for (int pos = 0; pos < seq_length; ++pos) {
        if (c_rna_solution[pos] == '.' || c_rna_solution[pos] == 'X') {
            if (ires_len > 0 && pos < ires_len) {
                c_rna_solution[pos] = ires_seq[pos];
            } else {
                c_rna_solution[pos] = rna_seq[pos];
            }
        }
    }

    std::cout << c_rna_solution << std::endl;
    outputfile << c_rna_solution << std::endl;
    std::cout << c_structure_solution << std::endl;
    outputfile << c_structure_solution << std::endl;

    // 計算 MFE 與 CAI
    T weighted_cai_score = GetCAIScore<T>(c_rna_solution, is_DN);
    double cai_value = GetUnweghtedCAIScore(c_rna_solution);
    cai_value = exp(cai_value / double(ami_seq.size()));

    double mfe_value;
    int cross_pair_count = 0;
    double cross_penalty_score = 0;
    double mfe_score_component = 0;

    if (ires_len == 0) { // 完全遵循 v1 的 MFE 計算方式
        if (is_DN) {
            mfe_value = -double(double(c_final_score - weighted_cai_score) / (100.0 * lambda));
        } else {
            if (lambda) {
                mfe_value = -double(double(c_final_score - weighted_cai_score) / 100.0);
            } else {
                mfe_value = -(c_final_score / 100.0);
            }
        }
    } else { // 遵循 i1 包含 crosspair 扣除的 MFE 計算方式
        cross_pair_count = getCrossPairNum(c_structure_solution, ires_len);
        cross_penalty_score = cross_pair_count * crosspair_penalty;

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
    }

    std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    std::cout << "CAI: " << std::round(cai_value * 1000.0) / 1000.0 << std::endl;
    outputfile << "CAI: " << std::round(cai_value * 1000.0) / 1000.0 << std::endl;

    if (ires_len > 0) {
        std::cout << "Crosspair count: " << cross_pair_count << std::endl;
        outputfile << "Crosspair count: " << cross_pair_count << std::endl;
    }

    auto c_end_time = std::chrono::high_resolution_clock::now();
    c_TimeSpend = std::chrono::duration_cast<std::chrono::duration<double>>(c_end_time - c_start_time).count();
    std::cout << "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;
    outputfile << "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;

    std::vector<pair<T,T>> result_container;
    result_container.push_back({round_up(cai_value), round_up(mfe_value)});

    // 格式保持與對應分支一致
    if (ires_len == 0) {
        outputcsv << "," << mfe_value << "," << cai_value << "," << c_TimeSpend;
    } else {
        outputcsv << "," << mfe_value << "," << cai_value << "," << cross_pair_count << "," << c_TimeSpend;
    }

    outputcsv.close();
    outputfile.close();

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
    std::cout << "=== CircularRNA (mfold) ===" << std::endl;

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

// 直接環化（用 Vienna 環化方式 - Sparse 結構版本，支援 IRES 計算）
template<typename T>
std::vector<pair<T,T>> circular_compose_v2_sparse(AllTables<T> &alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq, string& ires_seq, string& output_txt, string& output_csv, string& PID, int protein_length, double TimeSpend, bool show_score, bool is_beam, bool is_DN) { // return {CAI, MFE}
    std::ofstream outputfile(output_txt, std::ios::app);
    std::ofstream outputcsv(output_csv, std::ios::app);

    double c_TimeSpend = 0;
    auto c_start_time = std::chrono::high_resolution_clock::now();

    // 將 ires_seq 中的 T/t 轉換為 U，並確保大小寫統一
    for (char &base : ires_seq) {
        if (base == 'T' || base == 't') {
            base = 'U';
        } else if (base == 'u') {
            base = 'U';
        }
    }

    int seq_length = rna_seq.size();
    int ires_len = ires_seq.length();

    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Circular RNA Reconstruction:" << std::endl;
    if (ires_len == 0) {
        std::cout << "=== CircularRNA (v1 - Sparse) ===" << std::endl;
    } else {
        std::cout << "=== CircularRNA (v1 - Sparse IRES) ===" << std::endl;
    }

    std::vector<std::unordered_map<int, State<T>>> &bestN = alltables.bestN;
    std::vector<std::vector<std::array<State<T>, 257>>> &bestSS = alltables.bestSS;
    std::vector<std::vector<std::array<bool, 257>>> &bestSS_check = alltables.bestSS_check;
    std::vector<std::unordered_map<int, State<T>>> &bestF = alltables.bestF;
    std::vector<std::vector<std::unordered_map<int, State<T>>>> &bestSC = alltables.bestSC;
    std::vector<std::unordered_map<int, State<T>>> &bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>> &bestM2 = alltables.bestM2;

    // 建立 BestAminoMap (含安全 fallback，防止特定胺基酸未找到對應密碼子)
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

    // 安全取得 Best Codon 的 Lambda 函式
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

    // 安全包裝 compute_extra_cai (考慮 ires_len 造成的座標平移)
    auto safe_compute_extra_cai = [&](int nth, int nucn, const std::string& seq, bool is_dn, char type) -> std::pair<std::string, double> {
        if (nth < ires_len) {
            return {"", INT_MIN}; // 位於 IRES 區段，不參與 CDS CAI 計算
        }
        int cds_nth = nth - ires_len;
        if (cds_nth / 3 >= (int)seq.size()) {
            return {"", INT_MIN};
        }
        char aa = seq[cds_nth / 3];
        if (Amino_to_nucs.find(aa) == Amino_to_nucs.end()) {
            return {"", INT_MIN};
        }
        return compute_extra_cai(cds_nth, nucn, seq, is_dn, type);
    };

    CircularCandidate<T> best_struct;

    // ==========================================
    // 1. 遍歷 bestSC(i, j) 尋找外側閉合為 Hairpin 的環化結構
    //    bestSC[jth][span_len]，其中起點 ith = jth - span_len
    // ==========================================
    for (size_t jth = 0; jth < bestSC.size(); ++jth) {
        for (size_t span_len = 0; span_len < bestSC[jth].size(); ++span_len) {
            if (bestSC[jth][span_len].empty()) continue;

            int ith = (int)jth - (int)span_len;
            if (ith < 0) continue;

            int hairpin_len = seq_length - jth + ith;
            if (hairpin_len <= MAXLOOP && (hairpin_len > HAIRPIN_GAP || hairpin_len == 0)) {
                for (const auto& [id, state] : bestSC[jth][span_len]) {
                    int ith_id, nuci, nucj;
                    std::tie(ith_id, nuci, nucj) = GetIndexTuple(id);

                    double cur_score = state.score;

                    // CAI 修正
                    auto [codon_j, extra_cai_j] = safe_compute_extra_cai(jth, nucj, ami_seq, is_DN, 'j');
                    auto [codon_i, extra_cai_i] = safe_compute_extra_cai(ith, nuci, ami_seq, is_DN, 'i');
                    double total_extra_cai = ((extra_cai_j == INT_MIN) ? 0 : extra_cai_j);

                    // MFE: 計算外側 Hairpin Loop 能量
                    double extra_mfe = hairpin37[std::min(hairpin_len, 30)];
                    extra_mfe += (hairpin_len > 30) ? (int)(lxc37 * log((hairpin_len) / 30.0)) : 0;
                    
                    int cross_delta = 0;
                    if (ires_len > 0) {
                        cross_delta = isCrossTalk(ith, jth, ires_len);
                    }

                    if (is_DN) extra_mfe *= lambda;

                    double total_score = cur_score + extra_mfe - total_extra_cai;
                    if (ires_len > 0) {
                        total_score -= (cross_delta * crosspair_penalty);
                    }

                    if (best_struct.score < total_score) {
                        best_struct.score = total_score;
                        best_struct.index = id;
                        best_struct.i = ith;
                        best_struct.j = jth;
                        best_struct.codon_i = codon_i;
                        best_struct.codon_j = codon_j;
                        best_struct.state = state;
                        best_struct.manner = state.MANNER;
                        if (ires_len > 0) {
                            best_struct.cross_pair_num = state.cross_pair_num + cross_delta;
                        }
                    }
                }
            }
        }
    }

    // ==========================================
    // 2. 遍歷 bestM2 尋找閉合為 Internal Loop 或 Multiloop 的環化結構
    //    if M2.pair_num == 2 => Internal loop
    //    if M2.pair_num > 2  => Multiloop
    // ==========================================
    for (size_t jth = 0; jth < bestM2.size(); ++jth) {
        for (const auto& [id, state] : bestM2[jth]) {
            int ith, nuci, nucj;
            std::tie(ith, nuci, nucj) = GetIndexTuple(id);

            int loop_len = seq_length - jth + ith;
            if (loop_len > MAXLOOP) continue;

            double cur_score = state.score;

            // CAI 修正
            auto [codon_j, extra_cai_j] = safe_compute_extra_cai(jth, nucj, ami_seq, is_DN, 'j');
            auto [codon_i, extra_cai_i] = safe_compute_extra_cai(ith, nuci, ami_seq, is_DN, 'i');
            double total_extra_cai = ((extra_cai_j == INT_MIN) ? 0 : extra_cai_j);

            double extra_mfe = 0;

            if (state.pair_num == 2) {
                // Internal Loop
                extra_mfe += calcInternalMFE<T>(alltables, rna_seq, ami_seq, BestAminoMap, state,
                                                id, ith, jth, codon_i, codon_j, seq_length, true, false, ires_len);
            } else if (state.pair_num > 2) {
                // Multiloop 閉合
                extra_mfe -= (ML_intern37 + ML_closing37);
            } else {
                continue;
            }

            if (is_DN) extra_mfe *= lambda;

            double total_score = cur_score + extra_mfe - total_extra_cai;
            // 注意：M2 狀態下的 cross_pair_num 已在 DP 過程中累加，若閉合處無額外 crosstalk 判定則直接繼承
            if (best_struct.score < total_score) {
                best_struct.score = total_score;
                best_struct.index = id;
                best_struct.i = ith;
                best_struct.j = jth;
                best_struct.codon_i = codon_i;
                best_struct.codon_j = codon_j;
                best_struct.state = state;
                best_struct.manner = state.MANNER;
                if (ires_len > 0) {
                    best_struct.cross_pair_num = state.cross_pair_num;
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

    // 回溯結構與序列
    BacktrackStructure_sparse(c_rna_solution, c_structure_solution, seq_length, best_struct.index, best_struct.state, alltables, best_struct.j);

    // 補齊 ith 與 jth 所屬 Codon (注意加上 ires_len 偏移量)
    if (best_struct.codon_i != "" && best_struct.i >= ires_len) {
        int cds_offset = best_struct.i - ires_len;
        int codon_start = ires_len + (cds_offset / 3) * 3;
        for (int k = 0; k < 3; ++k) {
            if (best_struct.codon_i[k] == '.') continue;
            int idx = codon_start + k;
            if (idx < seq_length) c_rna_solution[idx] = best_struct.codon_i[k];
        }
    }
    if (best_struct.codon_j != "" && best_struct.j >= ires_len) {
        int cds_offset = best_struct.j - ires_len;
        int codon_start = ires_len + (cds_offset / 3) * 3;
        for (int k = 0; k < 3; ++k) {
            if (best_struct.codon_j[k] == '.') continue;
            int idx = codon_start + k;
            if (idx < seq_length) c_rna_solution[idx] = best_struct.codon_j[k];
        }
    }

    // 1. 若有 IRES，前 ires_len 個鹼基直接填入 ires_seq
    if (ires_len > 0) {
        for (int pos = 0; pos < ires_len && pos < seq_length; ++pos) {
            c_rna_solution[pos] = ires_seq[pos];
        }
    }

    // 2. 補齊結構前段 CDS: 從 ires_len 到 best_struct.i 之前
    int cds_start_ami = 0;
    int cds_end_ami_front = (best_struct.i > ires_len) ? (best_struct.i - ires_len) / 3 : 0;

    for (int ami = cds_start_ami; ami < cds_end_ami_front && ami < (int)ami_seq.size(); ++ami) {
        std::string codon = safe_get_best_codon(ami_seq[ami]);
        int base_idx = ires_len + ami * 3;
        for (int k = 0; k < 3; ++k) {
            int idx = base_idx + k;
            if (idx < seq_length) c_rna_solution[idx] = codon[k];
        }
    }

    // 3. 補齊結構後段 CDS: 從 best_struct.j 之後到 CDS 結尾
    int cds_start_ami_back = 0;
    if (best_struct.j >= ires_len) {
        cds_start_ami_back = ((best_struct.j - ires_len) / 3) + 1;
    }

    for (int ami = cds_start_ami_back; ami < (int)ami_seq.size(); ++ami) {
        std::string codon = safe_get_best_codon(ami_seq[ami]);
        int base_idx = ires_len + ami * 3;
        for (int k = 0; k < 3; ++k) {
            int idx = base_idx + k;
            if (idx < seq_length) c_rna_solution[idx] = codon[k];
        }
    }

    // 4. 未配對部分回填原始 rna_seq (再次確保前 ires_len 鎖定為 ires_seq)
    for (int pos = 0; pos < seq_length; ++pos) {
        if (c_rna_solution[pos] == '.' || c_rna_solution[pos] == 'X') {
            if (ires_len > 0 && pos < ires_len) {
                c_rna_solution[pos] = ires_seq[pos];
            } else {
                c_rna_solution[pos] = rna_seq[pos];
            }
        }
    }

    std::cout << c_rna_solution << std::endl;
    outputfile << c_rna_solution << std::endl;
    std::cout << c_structure_solution << std::endl;
    outputfile << c_structure_solution << std::endl;

    // 計算 CDS 的 CAI（截取 pos >= ires_len 的序列）
    std::string cds_rna = (ires_len > 0 && seq_length >= ires_len) ? c_rna_solution.substr(ires_len) : c_rna_solution;
    T weighted_cai_score = GetCAIScore<T>(cds_rna, is_DN);
    double cai_value = GetUnweghtedCAIScore(cds_rna);
    cai_value = exp(cai_value / double(ami_seq.size()));

    double mfe_value;
    int cross_pair_count = 0;
    double cross_penalty_score = 0;
    double mfe_score_component = 0;

    if (ires_len == 0) { // 完全遵循無 IRES 的 MFE 計算方式
        if (is_DN) {
            mfe_value = -double(double(c_final_score - weighted_cai_score) / (100.0 * lambda));
        } else {
            if (lambda) {
                mfe_value = -double(double(c_final_score - weighted_cai_score) / 100.0);
            } else {
                mfe_value = -(c_final_score / 100.0);
            }
        }
    } else { // 遵循含 IRES 的 Crosspair 扣除計算方式
        cross_pair_count = getCrossPairNum(c_structure_solution, ires_len);
        cross_penalty_score = cross_pair_count * crosspair_penalty;

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
    }

    std::cout << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    outputfile << "Folding free energy: " << mfe_value << " kcal/mol" << std::endl;
    std::cout << "CAI: " << std::round(cai_value * 1000.0) / 1000.0 << std::endl;
    outputfile << "CAI: " << std::round(cai_value * 1000.0) / 1000.0 << std::endl;

    if (ires_len > 0) {
        std::cout << "Crosspair count: " << cross_pair_count << std::endl;
        outputfile << "Crosspair count: " << cross_pair_count << std::endl;
    }

    auto c_end_time = std::chrono::high_resolution_clock::now();
    c_TimeSpend = std::chrono::duration_cast<std::chrono::duration<double>>(c_end_time - c_start_time).count();
    std::cout << "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;
    outputfile << "Total runtime: " << std::round(c_TimeSpend * 1000.0) / 1000.0 << " s" << std::endl;

    std::vector<pair<T,T>> result_container;
    result_container.push_back({round_up(cai_value), round_up(mfe_value)});

    // 格式保持與對應分支一致
    if (ires_len == 0) {
        outputcsv << "," << mfe_value << "," << cai_value << "," << c_TimeSpend;
    } else {
        outputcsv << "," << mfe_value << "," << cai_value << "," << cross_pair_count << "," << c_TimeSpend;
    }

    outputcsv.close();
    outputfile.close();

    return result_container;
}

// 直接環化（用 mfold 環化方式）
template<typename T>
std::vector<pair<T,T>> circular_compose_m1_sparse(AllTables<T> &alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq, string& output_txt, string& output_csv, string& PID, int protein_length, double TimeSpend, bool show_score, bool is_beam, bool is_DN) { // return {CAI, MFE}
    std::ofstream outputfile(output_txt, std::ios::app);
    std::ofstream outputcsv(output_csv, std::ios::app);

    double c_TimeSpend = 0;
    auto c_start_time = std::chrono::high_resolution_clock::now();

    outputfile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    outputfile << "Circular RNA Reconstruction:" << std::endl;
    std::cout << "=== CircularRNA (mfold - Sparse) ===" << std::endl;

    // 型態為 std::vector<std::vector<std::unordered_map<int, State<T>>>>
    auto& bestSC = alltables.bestSC;

    int seq_length = rna_seq.size() / 2;

    int best_i = -1;
    int best_j = -1;
    int best_id1 = -1;
    int best_id2 = -1;
    State<T> best_state1;
    State<T> best_state2;
    T best_score = -std::numeric_limits<T>::infinity();

    // 1. 遍歷起點 i 與分割點 j
    // 區間 1: C(i, j)                 => 終點為 j,               跨度長度 len1 = j - i
    // 區間 2: C(j, i + seq_length)   => 終點為 i + seq_length,  跨度長度 len2 = (i + seq_length) - j
    for (int i = 0; i < seq_length; ++i) {
        int second_end = i + seq_length;
        if (second_end >= (int)bestSC.size()) continue;

        for (int j = i + 1; j < second_end; ++j) {
            int len1 = j - i;
            int len2 = second_end - j;

            // 邊界檢查
            if (j >= (int)bestSC.size() || len1 >= (int)bestSC[j].size()) continue;
            if (len2 >= (int)bestSC[second_end].size()) continue;

            // 取得對應的 unordered_map<int, State<T>>
            const auto& map1 = bestSC[j][len1];
            const auto& map2 = bestSC[second_end][len2];

            if (map1.empty() || map2.empty()) continue;

            // 遍歷 map 尋找最優配對
            for (const auto& [id1, state1] : map1) {
                int ith1, nuci1, nucj1;
                std::tie(ith1, nuci1, nucj1) = GetIndexTuple(id1);

                for (const auto& [id2, state2] : map2) {
                    int ith2, nuci2, nucj2;
                    std::tie(ith2, nuci2, nucj2) = GetIndexTuple(id2);

                    // 核苷酸相容性驗證：
                    // state1 是 C(i, j)，其核苷酸分別為 nuci1, nucj1
                    // state2 是 C(j, i+L)，其核苷酸分別為 nucj2, nuci2
                    if (nuci1 == nuci2 && nucj1 == nucj2) {
                        T combined_score = state1.score + state2.score;

                        if (combined_score > best_score) {
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
    }

    if (best_i == -1) {
        std::cout << "No circular structure found." << std::endl;
        outputfile << "No circular structure found." << std::endl;
        return {};
    }

    // 2. 回溯雙倍長度的序列與結構 (長度 2 * seq_length)
    std::string full_rna(rna_seq.size(), '.');
    std::string full_struct(rna_seq.size(), '.');

    // 分別 Backtrack 兩個配對區塊 C(i, j) 與 C(j, i + seq_length)
    BacktrackStructure_sparse(full_rna, full_struct, rna_seq.size(), best_id1, best_state1, alltables, best_j);
    BacktrackStructure_sparse(full_rna, full_struct, rna_seq.size(), best_id2, best_state2, alltables, best_i + seq_length);

    // 未配對部分填入 rna_seq 原始核苷酸
    for (size_t pos = 0; pos < full_rna.size(); ++pos) {
        if (full_rna[pos] == '.' || full_rna[pos] == 'X') {
            full_rna[pos] = rna_seq[pos];
        }
    }

    // 3. 截取 [best_i, best_i + seq_length - 1] 作為單倍長度的 CDS 與二級結構輸出
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

    // 將序列填回 0 ~ seq_length-1
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
    cai_value = exp(cai_value / double(ami_seq.size() / 2));

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
