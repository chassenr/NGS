#!/bin/bash
#$ -cwd 
#$ -j y
#$ -pe threaded 1
#$ -V

${REFORMAT} in=./Merged/${SGE_TASK_ID}".assembled.fastq" out=./Swarm/${SGE_TASK_ID}"_good.fasta" fastawrap=1000 
