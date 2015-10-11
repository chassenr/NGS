#!/bin/bash
#$ -cwd 
#$ -j y
#$ -pe threaded 1
#$ -V

if [ -f ${output}.efetch ]
then
 rm ${output}.efetch
fi

while read line
do
 echo ${line} > ${line}.uid.tmp
 efetch -db protein -id ${line} -mode xml | grep "<OrgName_lineage>" | head -1 > ${line}.efetch.tmp
 paste ${line}.efetch.tmp ${line}.uid.tmp >> ${output}.efetch
 rm ${line}.efetch.tmp
 rm ${line}.uid.tmp
done < ${input}.uid
sed -e 's/^ *//' -e 's/<OrgName_lineage>//' -e 's/<\/OrgName_lineage>//' ${output}.efetch > ${output}.txt


if [ $(cut -f1 ${output}.txt | sed '/^$/d' | wc -l) -ne $(cut -f2 ${output}.txt | sed '/^$/d' | wc -l) ]
then
 mv ${output}.txt ${output}.missing.txt
 awk 'BEGIN {FS="\t"}  $1!="" {print}' ${output}.missing.txt > ${output}.txt
 
 until [ $(ls -l ${output}.missing.txt | cut -d' ' -f5) -eq $(ls -l ${output}.txt | cut -d' ' -f5) ]
 do
  if [ -f tmp.efetch ]
  then
   rm tmp.efetch
  fi

  awk 'BEGIN {FS="\t"}  $1=="" {print}' ${output}.missing.txt | cut -f2 > repeat.uid
  
  while read line
  do
   echo ${line} > ${line}.uid.tmp
   efetch -db protein -id ${line} -mode xml | grep "<OrgName_lineage>" | head -1 > ${line}.efetch.tmp
   paste ${line}.efetch.tmp ${line}.uid.tmp >> tmp.efetch
   rm ${line}.efetch.tmp
   rm ${line}.uid.tmp
  done < repeat.uid
  
  sed -e 's/^ *//' -e 's/<OrgName_lineage>//' -e 's/<\/OrgName_lineage>//' tmp.efetch >> ${output}.txt
 done
fi



