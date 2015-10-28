#!/bin/bash
#$ -cwd 
#$ -j y
#$ -pe threaded 4
#$ -V

${PEAR} -j 4 -v 10 -n 350 -m 500 -f ./Trimmed/${SGE_TASK_ID}"_ptrim_R1.fastq" -r ./Trimmed/${SGE_TASK_ID}"_ptrim_R2.fastq" -o ./Merged/${SGE_TASK_ID}  > ./Merged/${SGE_TASK_ID}"_merged.log"

