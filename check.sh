#!/bin/bash

# Usage: ./check.sh {long/short}

# read dataset: 
#   dataset/test_target_${long}.fasta
# generate result: 
#   - {LD, DN}_{exact, beam}_{vienna, mfold}_target.txt
#   - {LD, DN}_{exact, beam}_{vienna, mfold}_target.csv
# check result with answer: 
#   - {LD, DN}_{exact, beam}_{vienna, mfold}_{long, short}_answer.csv

if [ -f LinearCDSfold ]; then
    echo "clean up previous build..."
    make clean
fi

echo "building LinearCDSfold..."
make

if [ -z "$1" ]; then
    echo "Usage: $0 {long/short}"
    exit 1
fi

# Default Methods
case1="LD_exact_vienna"
case2="DN_beam_vienna"
case3="DN_beam2_vienna"

# Linear need to exec for 2n CDS
case4="LD_exact_mfold"
case5="DN_beam_mfold"
case6="DN_beam2_mfold"

# need IRES file + kappa
case7="LD_exact_IRES"
case8="DN_beam_IRES"
case9="DN_beam2_IRES"

# Default + Sparse
case10="LD_exact_vienna_sparse"
case11="DN_beam_vienna_sparse"
case12="DN_beam2_vienna_sparse"

./LinearCDSfold -l 2 -O LD -m exact -C vienna -o ${case1}_target.txt -f ${case1}_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam -b 100 -C vienna -o ${case2}_target.txt -f ${case2}_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam2 -b 100 -C vienna -o ${case3}_target.txt -f ${case3}_target.csv dataset/test_target_${1}.fasta

./LinearCDSfold -l 2 -O LD -m exact -C mfold -o ${case4}_target.txt -f ${case4}_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam -b 100 -C mfold -o ${case5}_target.txt -f ${case5}_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam2 -b 100 -C mfold -o ${case6}_target.txt -f ${case6}_target.csv dataset/test_target_${1}.fasta

python3 testing/Check.py ${case1}_target.csv testing/${case1}_${1}_answer.csv
python3 testing/Check.py ${case2}_target.csv testing/${case2}_${1}_answer.csv
python3 testing/Check.py ${case3}_target.csv testing/${case3}_${1}_answer.csv

python3 testing/Check.py ${case4}_target.csv testing/${case4}_${1}_answer.csv
python3 testing/Check.py ${case5}_target.csv testing/${case5}_${1}_answer.csv
python3 testing/Check.py ${case6}_target.csv testing/${case6}_${1}_answer.csv

# ./LinearCDSfold -l 2 -O LD -m exact dataset/A9JR22.fasta