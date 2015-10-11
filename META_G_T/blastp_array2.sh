#!/bin/bash
#$ -cwd 
#$ -j y
#$ -pe threaded 4
#$ -V

${BLASTP} -query out.${SGE_TASK_ID} -task blastp -db ${BLASTP_DB} -out blastp_${SGE_TASK_ID}.txt -evalue 1e-5 -outfmt '6 std qcovs salltitles' -max_target_seqs 10 -num_threads 4
