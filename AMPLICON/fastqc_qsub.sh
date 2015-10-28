#!/bin/bash
#$ -cwd 
#$ -j y
#$ -pe threaded 1
#$ -V

cd ./FastQC/
${FASTQC} -o . ../Merged/$SGE_TASK_ID".assembled.fastq"
unzip $SGE_TASK_ID".assembled_fastqc.zip"
cd ..
