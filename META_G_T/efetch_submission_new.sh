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
  
  until [ $(cut -f1 ${output}.txt | sed '/^$/d' | wc -l) -eq $(cut -f2 ${output}.txt | sed '/^$/d' | wc -l) ]
  do

    mv ${output}.txt ${output}.tmp.txt
    awk 'BEGIN {FS="\t"}  $1=="" {print}' ${output}.tmp.txt | cut -f2 > repeat.uid
    awk 'BEGIN {FS="\t"}  $1!="" {print}' ${output}.tmp.txt > ${output}.txt

    if [ -f repeat.efetch ]
    then
      rm repeat.efetch
    fi

    while read line
    do
      echo ${line} > ${line}.uid.tmp
      efetch -db protein -id ${line} -mode xml | grep "<OrgName_lineage>" | head -1 > ${line}.efetch.tmp
      paste ${line}.efetch.tmp ${line}.uid.tmp >> repeat.efetch
      rm ${line}.efetch.tmp
      rm ${line}.uid.tmp
    done < repeat.uid
  
    sed -e 's/^ *//' -e 's/<OrgName_lineage>//' -e 's/<\/OrgName_lineage>//' repeat.efetch > repeat.txt
    cat repeat.txt >> ${output}.txt

    awk 'BEGIN {FS="\t"}  $1=="" {print}' repeat.txt | cut -f2 > left.uid
    
    if [ $(cat left.uid | wc -l ) -eq $(cat repeat.uid | wc -l) ]
    then
      break
    fi

  done
fi
