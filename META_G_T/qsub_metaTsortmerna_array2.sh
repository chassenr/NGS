#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe threaded 6
#$ -V

./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db --reads ${FILELOCATION}/sort${i}/out.${SGE_TASK_ID} --aligned ${FILELOCATION}/sort${i}/rRNA.${SGE_TASK_ID} --other ${FILELOCATION}/sort${i}/non_rRNA.${SGE_TASK_ID} --sam --fastx --blast 3 --best 1 --log -v -a 6
