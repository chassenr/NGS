#!/bin/bash
#$ -cwd 
#$ -j y
#$ -pe threaded 4
#$ -V

java -jar ${TRIMMOMATIC} PE -threads 4 -trimlog ./Trimmed/$SGE_TASK_ID"_trim.log" ./Clipped/${SGE_TASK_ID}"_clip_R1.fastq" ./Clipped/${SGE_TASK_ID}"_clip_R2.fastq" ./Trimmed/${SGE_TASK_ID}"_ptrim_R1.fastq" ./Trimmed/${SGE_TASK_ID}"_strim_R1.fastq" ./Trimmed/${SGE_TASK_ID}"_ptrim_R2.fastq" ./Trimmed/${SGE_TASK_ID}"_strim_R2.fastq" SLIDINGWINDOW:4:15 MINLEN:100

