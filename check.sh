#!/bin/bash

# Usage: ./check.sh {long/short/ID}

# read dataset: 
#   dataset/test_target_${long}.fasta
# generate result: 
#   - {LD, DN}_{exact, beam, sparse}_{vienna, mfold}_target.txt
#   - {LD, DN}_{exact, beam, sparse}_{vienna, mfold}_target.csv
# check result with answer: 
#   - {LD, DN}_{exact, beam, sparse}_{vienna, mfold}_{long, short}_answer.csv

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
./LinearCDSfold -l 2 -O LD -m exact -C vienna -o LD_exact_vienna_target.txt -f LD_exact_vienna_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam -b 100 -C vienna -o DN_beam_vienna_target.txt -f DN_beam_vienna_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam2 -b 100 -C vienna -o DN_beam2_vienna_target.txt -f DN_beam2_vienna_target.csv dataset/test_target_${1}.fasta

python3 testing/Check.py LD_exact_vienna_target.csv testing/LD_exact_vienna_${1}_answer.csv
python3 testing/Check.py DN_beam_vienna_target.csv testing/DN_beam_vienna_${1}_answer.csv
python3 testing/Check.py DN_beam2_vienna_target.csv testing/DN_beam2_vienna_${1}_answer.csv

# Linear need to exec for 2n CDS
./LinearCDSfold -l 2 -O LD -m exact -C mfold -o LD_exact_mfold_target.txt -f LD_exact_mfold_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam -b 100 -C mfold -o DN_beam_mfold_target.txt -f DN_beam_mfold_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam2 -b 100 -C mfold -o DN_beam2_mfold_target.txt -f DN_beam2_mfold_target.csv dataset/test_target_${1}.fasta

python3 testing/Check.py LD_exact_mfold_target.csv testing/LD_exact_mfold_${1}_answer.csv
python3 testing/Check.py DN_beam_mfold_target.csv testing/DN_beam_mfold_${1}_answer.csv
python3 testing/Check.py DN_beam2_mfold_target.csv testing/DN_beam2_mfold_${1}_answer.csv

# Default + Sparse
./LinearCDSfold -l 2 -O LD -m sparse -C vienna -o LD_sparse_vienna_target.txt -f LD_sparse_vienna_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O LD -m sparse -C mfold -o LD_sparse_mfold_target.txt -f LD_sparse_mfold_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m sparse -C vienna -o DN_sparse_vienna_target.txt -f DN_sparse_vienna_target.csv dataset/test_target_${1}.fasta

python3 testing/Check.py LD_sparse_vienna_target.csv testing/LD_sparse_vienna_${1}_answer.csv
python3 testing/Check.py LD_sparse_mfold_target.csv testing/LD_sparse_mfold_${1}_answer.csv
python3 testing/Check.py DN_sparse_vienna_target.csv testing/DN_sparse_vienna_${1}_answer.csv

# need IRES file + kappa
./LinearCDSfold -l 2 -O LD -m exact -C IRES -k 100 -i vir_ires_00036.fasta -o LD_exact_IRES_target.txt -f LD_exact_IRES_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam -b 100 -C IRES -k 100 -i vir_ires_00036.fasta -o DN_beam_IRES_target.txt -f DN_beam_IRES_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O DN -m beam2 -b 100 -C IRES -k 100 -i vir_ires_00036.fasta -o DN_beam2_IRES_target.txt -f DN_beam2_IRES_target.csv dataset/test_target_${1}.fasta
./LinearCDSfold -l 2 -O LD -m sparse -C IRES -k 100 -i vir_ires_00036.fasta -o LD_sparse_IRES_target.txt -f LD_sparse_IRES_target.csv dataset/test_target_${1}.fasta

python3 testing/Check.py LD_exact_IRES_target.csv testing/LD_exact_IRES_${1}_answer.csv
python3 testing/Check.py DN_beam_IRES_target.csv testing/DN_beam_IRES_${1}_answer.csv
python3 testing/Check.py DN_beam2_IRES_target.csv testing/DN_beam2_IRES_${1}_answer.csv
python3 testing/Check.py LD_sparse_IRES_target.csv testing/LD_sparse_IRES_${1}_answer.csv