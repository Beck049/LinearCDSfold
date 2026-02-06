#!/bin/bash

# Usage: ./check.sh {long/short}

# read dataset: 
#   dataset/test_target_${long}.fasta
# generate result: 
#   - LD_exact_target.txt
#   - DN_beam_target.txt
# check result with answer: 
#   - LD_exact_${long}_answer.csv
#   - DN_beam_${long}_answer.csv

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

case1="LD_exact"
case2="DN_beam"

./LinearCDSfold -l 2 -O LD -m exact -o ${case1}_target.txt -f ${case1}_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam -b 100 -o ${case2}_target.txt -f ${case2}_target.csv dataset/test_target_${1}.fasta

python3 testing/Check.py ${case1}_target.txt testing/${case1}_${1}_answer.csv
python3 testing/Check.py ${case2}_target.txt testing/${case2}_${1}_answer.csv