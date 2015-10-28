#!/bin/bash
#$ -cwd 
#$ -j y
#$ -pe threaded 6
#$ -V

${SINA} -i out.$SGE_TASK_ID --intype fasta -o sina_out.$SGE_TASK_ID --outtype fasta --search --meta-fmt csv --overhang remove --insertion forbid --filter none --fs-kmer-no-fast --fs-kmer-len 10 --fs-req 2 --fs-req-full 1 --fs-min 40 --fs-max 40 --fs-weight 1 --fs-full-len 1400 --fs-msc 0.7 --match-score 1 --mismatch-score -1 --pen-gap 5 --pen-gapext 2 --search-cover query --search-iupac optimistic --search-min-sim 0.9 --turn all --lca-quorum 0.7 --search-db ${SINA_PT} --ptdb ${SINA_PT} --lca-fields tax_slv
