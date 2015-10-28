#!/bin/bash
#$ -cwd 
#$ -j y
#$ -pe threaded 1
#$ -V

cd ./Swarm/
grep -v "^>" ${SGE_TASK_ID}"_good.fasta" | \
grep -v [^ACGTacgt] | sort -d | uniq -c | \
while read abundance sequence ; do
    hash=$(printf "${sequence}" | sha1sum)
    hash=${hash:0:40}
    printf ">%s_%d_%s\n" "${hash}" "${abundance}" "${sequence}"
done | sort -t "_" -k2,2nr -k1.2,1d | \
sed -e 's/\_/\n/2' > ${SGE_TASK_ID}"_dereplicated.fasta"
cd ..
