#!/bin/bash

# Usage: ./check.sh {long/short/ID}

# read dataset: 
#   dataset/test_target_${long}.fasta
# generate result: 
#   - {LD, DN}_{exact, beam, beam2, sparse}_{vienna, mfold, IRES}_target.txt
#   - {LD, DN}_{exact, beam, beam2, sparse}_{vienna, mfold, IRES}_target.csv
# check result with answer: 
#   - {LD, DN}_{exact, beam, sparse}_{vienna, mfold, IRES}_{long, short}_answer.csv

if [ -f LinearCDSfold ]; then
    echo "clean up previous build..."
    make clean
fi

echo "building LinearCDSfold..."
make

if [ -z "$1" ]; then
    echo "Usage: {long/short/ProtinId}"
    exit 1
fi

if [ "$1" = "long" ] || [ "$1" = "short" ]; then
    TARGET_NAME="test_target_${1}"
else
    TARGET_NAME="${1}"
fi

MAX_JOBS=2
pids=()
run_job() {
    # 若目前背景執行的任務數 >= MAX_JOBS，就等其中一個跑完再繼續
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 0.3
    done
    "$@" > /dev/null &
    pids+=($!)
}
echo "Running all LinearCDSfold jobs in parallel (Max: $MAX_JOBS)..."

# Default Methods
run_job ./LinearCDSfold -l 2 -O LD -m exact -C vienna -o LD_exact_vienna_target.txt -f LD_exact_vienna_target.csv dataset/${TARGET_NAME}.fasta
run_job ./LinearCDSfold -l 2 -O DN -m beam -b 100 -C vienna -o DN_beam_vienna_target.txt -f DN_beam_vienna_target.csv dataset/${TARGET_NAME}.fasta
run_job ./LinearCDSfold -l 2 -O DN -m beam2 -b 100 -C vienna -o DN_beam2_vienna_target.txt -f DN_beam2_vienna_target.csv dataset/${TARGET_NAME}.fasta
# Linear need to exec for 2n CDS
run_job ./LinearCDSfold -l 2 -O LD -m exact -C mfold -o LD_exact_mfold_target.txt -f LD_exact_mfold_target.csv dataset/${TARGET_NAME}.fasta
run_job ./LinearCDSfold -l 2 -O DN -m beam -b 100 -C mfold -o DN_beam_mfold_target.txt -f DN_beam_mfold_target.csv dataset/${TARGET_NAME}.fasta
run_job ./LinearCDSfold -l 2 -O DN -m beam2 -b 100 -C mfold -o DN_beam2_mfold_target.txt -f DN_beam2_mfold_target.csv dataset/${TARGET_NAME}.fasta
# Default + Sparse
run_job ./LinearCDSfold -l 2 -O LD -m sparse -C vienna -o LD_sparse_vienna_target.txt -f LD_sparse_vienna_target.csv dataset/${TARGET_NAME}.fasta
run_job ./LinearCDSfold -l 2 -O LD -m sparse -C mfold -o LD_sparse_mfold_target.txt -f LD_sparse_mfold_target.csv dataset/${TARGET_NAME}.fasta
run_job ./LinearCDSfold -l 2 -O DN -m sparse -C vienna -o DN_sparse_vienna_target.txt -f DN_sparse_vienna_target.csv dataset/${TARGET_NAME}.fasta
# need IRES file + kappa
run_job ./LinearCDSfold -l 2 -O LD -m exact -C IRES -k 0 -i vir_ires_00036.fasta -o LD_exact_IRES_target.txt -f LD_exact_IRES_target.csv dataset/${TARGET_NAME}.fasta
run_job ./LinearCDSfold -l 2 -O DN -m beam -b 100 -C IRES -k 0 -i vir_ires_00036.fasta -o DN_beam_IRES_target.txt -f DN_beam_IRES_target.csv dataset/${TARGET_NAME}.fasta
run_job ./LinearCDSfold -l 2 -O DN -m beam2 -b 100 -C IRES -k 0 -i vir_ires_00036.fasta -o DN_beam2_IRES_target.txt -f DN_beam2_IRES_target.csv dataset/${TARGET_NAME}.fasta
run_job ./LinearCDSfold -l 2 -O LD -m sparse -C IRES -k 0 -i vir_ires_00036.fasta -o LD_sparse_IRES_target.txt -f LD_sparse_IRES_target.csv dataset/${TARGET_NAME}.fasta

# 等待上面所有背景程序執行完畢
for pid in "${pids[@]}"; do
    wait "$pid" || echo "Warning: Job PID $pid failed!"
done
echo "All LinearCDSfold jobs finished. Starting Check.py..."

# Default Methods
python3 testing/Check.py    LD_exact_vienna_target.csv    testing/LD_exact_vienna_answer.csv
python3 testing/Check.py    DN_beam_vienna_target.csv     testing/DN_beam_vienna_answer.csv
python3 testing/Check.py    DN_beam2_vienna_target.csv    testing/DN_beam2_vienna_answer.csv
# Linear need to exec for 2n CDS
python3 testing/Check.py    LD_exact_mfold_target.csv     testing/LD_exact_mfold_answer.csv
python3 testing/Check.py    DN_beam_mfold_target.csv      testing/DN_beam_mfold_answer.csv
python3 testing/Check.py    DN_beam2_mfold_target.csv     testing/DN_beam2_mfold_answer.csv
# Default + Sparse
python3 testing/Check.py    LD_sparse_vienna_target.csv   testing/LD_sparse_vienna_answer.csv
python3 testing/Check.py    LD_sparse_mfold_target.csv    testing/LD_sparse_mfold_answer.csv
python3 testing/Check.py    DN_sparse_vienna_target.csv   testing/DN_sparse_vienna_answer.csv
# need IRES file + kappa
python3 testing/Check.py    LD_exact_IRES_target.csv      testing/LD_exact_IRES_answer.csv
python3 testing/Check.py    DN_beam_IRES_target.csv       testing/DN_beam_IRES_answer.csv
python3 testing/Check.py    DN_beam2_IRES_target.csv      testing/DN_beam2_IRES_answer.csv
python3 testing/Check.py    LD_sparse_IRES_target.csv     testing/LD_sparse_IRES_answer.csv