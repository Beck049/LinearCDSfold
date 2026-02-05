if(LinearCDSfold exist )then make clean
make

 ./LinearCDSfold -l 2 -O LD -m exact -o LD_exact_target.txt -f LD_exact_target.csv dataset/test_target.fasta
 ./LinearCDSfold -l 2 -O DN -m beam -b 100 -o DN_exact_target.txt -f DN_beam_target.csv dataset/test_target.fasta

python3 Check.py LD_exact_target.txt LD_exact_answer.txt
python3 Check.py DN_beam_target.txt DN_beam_answer.txt

