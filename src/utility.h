#pragma once

#include "energy_parameter.h" // energy_parameter stuff

//  nuc-pairs : 0:NP 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA 7:NN
//  Vienna nuc: 0:N  1:A  2:C  3:G  4:U  5:u  6:g

#include "intl11.h"
#include "intl21.h"
#include "intl22.h"

/**
 * @brief Utility function on MFE 計算
 */

#define MAXLOOP 30

inline int MIN2(int a, int b) {if (a <= b)return a;else return b;}
inline int MAX2(int a, int b) {if (a >= b)return a;else return b;}

// 給定兩個 nucleotide x, y ，回傳它們的 pairing type
#define NUC_TO_PAIR(x,y) (x==1? (y==4?5:0) : (x==2? (y==3?1:0) : (x==3 ? (y==2?2:(y==4?3:0)) : (x==4 ? (y==3?4:(y==1?6:0)) : 0))))

/**
 * @brief 計算 Hairpin Loop 的自由能
 * @param i 髮夾環左側閉合鹼基的索引
 * @param j 髮夾環右側閉合鹼基的索引
 * @param base_i 閉合鹼基 i 的類型
 * @param base_ni 髮夾環內第一個鹼基 (i+1) 的類型 (Mismatch 用)
 * @param base_pj 髮夾環內最後一個鹼基 (j-1) 的類型 (Mismatch 用)
 * @param base_j 閉合鹼基 j 的類型
 * @return int 返回髮夾環的自由能值
 */
inline int v_score_hairpin(int i, int j, int base_i, int base_ni, int base_pj, int base_j) {
    int size = j-i-1;
    int type = NUC_TO_PAIR(base_i, base_j);
    int energy;

    if(size <= 30) {
        energy = hairpin37[size];
    } else {
        energy = hairpin37[30] + (int)(lxc37*log((size)/30.));
    }
    
    if(size == 3) return energy + (type>2?TerminalAU37:0); /* should only be the case when folding alignments */

    energy += mismatchH37[type][base_ni][base_pj];

    return energy;
}

/**
 * @brief 輔助功能：計算多分支環 (Multi-loop) 中單個螺旋臂 (Stem) 的能量貢獻
 * @param type 該螺旋臂閉合配對的類型
 */
inline int E_MLstem(int type) {
    int energy = 0;
    if(type > 2) {  energy += TerminalAU37; }
    energy += ML_intern37;
    return energy;
}

/**
 * @brief 計算多分支環 (Multi-loop) 閉合配對 (Closing Pair) 本身的能量
 * @param nuci 閉合配對左側的鹼基
 * @param nucj 閉合配對右側的鹼基
 * @note 因為是整個多分支環的「外圍閉合對」，視角反轉，故傳入時轉為 (nucj, nuci)
 */
inline int v_score_multi(int nuci, int nucj) {
	int tt = NUC_TO_PAIR(nucj, nuci); // : closing pair in multi: reversed
  return E_MLstem(tt) + ML_closing37;
}

/**
 * @brief 計算多分支環內部單個子螺旋臂 (Sub-stem) 的起始/局部能量
 * @param nuci 子螺旋左側鹼基
 * @param nuck 子螺旋右側鹼基
 */
inline int v_score_M1(int nuci, int nuck) {
  int tt = NUC_TO_PAIR(nuci, nuck);
  return E_MLstem(tt);
}

/**
 * @brief 計算外部環 (Exterior Loop) 中已配對螺旋末端的能量
 * @param nuci 外部配對左側鹼基
 * @param nucj 外部配對右側鹼基
 */
inline int v_score_external_paired(int nuci, int nucj) {
    int type = NUC_TO_PAIR(nuci, nucj);
    int energy = 0;

    if(type > 2) {  energy += TerminalAU37; }
    return energy;
}

/**
 * @brief 計算外部環中未配對鹼基的能量（Turner 模型中通常為 0）
 */
inline int v_score_external_unpaired(int i, int j) {
    return 0;
}

/**
 * @brief 計算內環 (Single/Internal Loop)、凸起環 (Bulge) 或堆疊 (Stack) 的綜合能量
 * @details 此函式涵蓋了所有由兩個配對對 (i,j) 和 (p,q) 圍成的局部結構（其中 i < p < q < j）
 * @param i 外層左鹼基位置, j 外層右鹼基位置
 * @param p 內層左鹼基位置, q 內層右鹼基位置
 * @param nuci~nucq1 各個關鍵位置的鹼基類型（用於特殊小環查表與 Mismatch 計算）
 */
inline int v_score_single(int i, int j, int p, int q,
                        int nuci, int nuci1, int nucj_1, int nucj,
                        int nucp_1, int nucp, int nucq, int nucq1){
  int type = NUC_TO_PAIR(nuci, nucj);
  int type_2 = NUC_TO_PAIR(nucq, nucp);
  int n1 = p-i-1;
  int n2 = j-q-1;
  int nl, ns, u, energy;
  energy = 0;

  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}

  if (nl == 0)
    return stack37[type][type_2];  /* stack */

  if (ns==0) {                      /* bulge */
      energy = (nl<=MAXLOOP)?bulge37[nl]:
      (bulge37[30]+(int)(lxc37*log(nl/30.)));
    if (nl==1) { 
      energy += stack37[type][type_2];
    } else {
      if (type>2) energy += TerminalAU37;
      if (type_2>2) energy += TerminalAU37;
    }
    return energy;
  } else {                            /* interior loop */
    if (ns==1) {
      if (nl==1)                    /* 1x1 loop */
        return int11_37[type][type_2][nuci1][nucj_1];
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1) {
          energy = int21_37[type][type_2][nuci1][nucq1][nucj_1];
        } else {
          energy = int21_37[type_2][type][nucq1][nuci1][nucp_1];
        }
        return energy;
      } else {  /* 1xn loop */
        energy = (nl+1<=MAXLOOP)?(internal_loop37[nl+1]) : (internal_loop37[30]+(int)(lxc37*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);
        energy += mismatch1nI37[type][nuci1][nucj_1] + mismatch1nI37[type_2][nucq1][nucp_1];
        return energy;
      }
    } else if (ns==2) {
      if(nl==2) {              /* 2x2 loop */
        return int22_37[type][type_2][nuci1][nucp_1][nucq1][nucj_1];
      } else if (nl==3){              /* 2x3 loop */
        energy = internal_loop37[5]+ninio37;
        energy += mismatch23I37[type][nuci1][nucj_1] + mismatch23I37[type_2][nucq1][nucp_1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      u = nl + ns;
      energy = (u <= MAXLOOP) ? (internal_loop37[u]) : (internal_loop37[30]+(int)(lxc37*log((u)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);

      energy += mismatchI37[type][nuci1][nucj_1] + mismatchI37[type_2][nucq1][nucp_1];
    }
  }
  return energy;
}

/**
 * @brief 計算內環外側閉合處 (Exterior closing pair) 的錯配能量
 * @details 用於模組化計算通用內環時，單獨提取外側的錯配交互作用能
 */
inline int v_score_exterior_mismatch(int i, int j, int nuci, int nuci1, int nucj_1, int nucj){
  int type = NUC_TO_PAIR(nuci, nucj);
  return mismatchI37[type][nuci1][nucj_1];
}

/**
 * @brief 計算內環內側閉合處 (Interior closing pair) 的錯配能量
 * @details 用於模組化計算通用內環時，單獨提取內側的錯配交互作用能
 */
inline int v_score_interior_mismatch(int p, int q, int nucp_1, int nucp, int nucq, int nucq1){
  int type_2 = NUC_TO_PAIR(nucq, nucp);
  return mismatchI37[type_2][nucq1][nucp_1];
}

/**
 * @brief 計算內環兩側長度不對稱造成的處罰能量 (Ninio Penalty)
 * @param n1 左側未配對長度, n2 右側未配對長度
 * @return int 返回不對稱處罰能量值，最高不超過 MAX_NINIO
 */
inline int v_score_asymmetry(int i, int j, int p, int q){
  int n1 = p-i-1;
  int n2 = j-q-1;
  int nl, ns;

  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}
  return MIN2(MAX_NINIO, (nl-ns)*ninio37);
}

/**
 * @brief 單純依據常規內環總長度 `l` 計算其基礎尺寸能量 (Size Energy)
 * @param l 內環中所有未配對鹼基的總數量 (l = nl + ns)
 */
inline int v_score_size(int l){
  int energy = (l <= MAXLOOP) ? (internal_loop37[l]) : (internal_loop37[30]+(int)(lxc37*log((l)/30.)));
  return energy;
}