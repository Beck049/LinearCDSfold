# 測試流程

## 參數簡介
| | Param | Description |
|---|---|---|
|`./LinearCDSfold`|`-l <lambda>`|object function 的權重|
||`-O <obj_func>`|`DN`: DERNA </br>`LD`: Linear Design|
||`-m <method>`|`exact`: 最佳解</br>`beam `: 近似解 (慢)</br>`beam2`: 近似解 (快)</br>`sparse`: Sparse|
||`-b <beam_size>`||
||`-C <circ_method>`|`vienna`: 用 vienna 直接環化</br>`mfold`: 用 mfold 直接環化</br>`IRES`: 用 vienna 加上 IRES 環化|
||`-i <IRES_file>`|IRES input file|
||`-k <crosspair_penalty>`|object function 中 crosspair 的權重|
||`-c <codon_usage_csv>`||
||`-o <output_txt>`||
||`-f <output_csv>`||
||`<input_fasta>`|輸入 蛋白質序列|

## 測試流程
分為 Long (大規模), Short (小規模) 測試

### 測試項目

```shell
./LinearCDSfold -l 2 -O LD -m exact -o LD_exact_target.txt -f LD_exact_target.csv dataset/test_target_${1}.fasta
```

```shell
./LinearCDSfold -l 2 -O DN -m beam -b 100 -o DN_beam_target.txt -f DN_beam_target.csv dataset/test_target_${1}.fasta
```

```shell
./LinearCDSfold -l 2 -O DN -m beam2 -b 100 -o DN_beam2_target.txt -f DN_beam2_target.csv dataset/test_target_${1}.fasta
```

### Output File Example
#### CSV
|UniPort ID|Protein length|Lambda|MFE|CAI|Run time|cMFE|cCAI|cRun time|
|------|---|-------|-----------|--------|-------|--------|--------|----------|
|A9JR22|95|2.000000|-165.100000|0.895415|4.35439|-172.197|0.886991|0.00242767|

#### TXT
```
Amino acid file: dataset/A9JR22.fasta
Codon usage table: codon_usage_freq_table_human.csv
Objective function: LinearDesign
Search mode: Exact search
===================================================================================================================
Protein ID: A9JR22
Protein length: 95
Lambda: 2.000
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Linear RNA Reconstruction:
Coding sequence and its secondary structure:
AUGGGGCUGCGGUACAGCAAGGAGGUGAGGGAUCGUCAUGGGGAUAAGGAUCCUGAGGGCAGGAUCCCUAUCACCCAGACGAUGCCUCAGACACUGUACGGCAGAUAUAAUUGUAAGAGCUGUUGGUUCGCCAACAAGGGGCUGAUCAAGUGCUCUAACCAUUAUAUCUGCCUGCGGUGUCUGACCUCCAUGCUGAACCGCACUGACUACUGUGAGAUCUGUGGGGAGGUGCUGCCCAAGCGGCUCACCUUUGAGACCACUCCCACAGCACCUCCCUACACCCCA
.(((((.((((((.(((((.((((((((((.((((((.((((((((.((((((((....)))))))).)))).)))))))))).))))((((((((((.((((((((((((.((.(((((..((((((.(((.......))).))))))..))))).)).)))))))))))))))))))))).)))))).))))).))))))((.((....)).))...(((((((((((((((....((.((((((....)))).)).))....))))))))))))))))))))
Folding free energy: -165.100 kcal/mol
CAI: 0.895
Total runtime: 4.354 s
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Circular RNA Reconstruction:
Coding sequence and its secondary structure :
AUGGGGCUGCGGUACAGCAAGGAGGUGAGGGAUCGUCAUGGGGAUAAGGAUCCUGAGGGCAGGAUCCCUAUCACCCAGACGAUGCCUCAGACACUGUACGGCAGAUAUAAUUGUAAGAGCUGUUGGUUCGCCAACAAGGGGCUGAUCAAGUGCUCUAACCAUUAUAUCUGCCUGCGGUGUCUGACCUCCAUGCUGAACCGCACUGACUACUGUGAGAUCUGUGGGGAGGUGCUGCCCAAGCGGCUCACCUUUGAGACCACUCCCACAGCACCUCCCUACACCCCG
..((((.((((((.(((((.((((((((((.((((((.((((((((.((((((((....)))))))).)))).)))))))))).))))((((((((((.((((((((((((.((.(((((..((((((.(((.......))).))))))..))))).)).)))))))))))))))))))))).)))))).))))).))))))((.((....)).))...(((((((((((((((....((.((((((....)))).)).))....))))))))))))))))))).
Folding free energy: -172.197 kcal/mol
CAI: 0.887
Total runtime: 0.002 s
```

## ToDo
1. 若 `-C IRES`，則
   1. 不輸出 Linear 結果！！
2. 實作 IRES Method
   1. IRES no CAI
3. Sparse Circular (IRES)