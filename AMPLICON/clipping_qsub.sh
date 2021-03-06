#!/bin/bash
#$ -cwd 
#$ -j y
#$ -pe threaded 1
#$ -V

# Let's process the reads in the FR orientation...
## Look for forward primer in R1. R2 is processed at the same time to make sure these files go 
## through the same processes and any sequences discarded match between outputs.
## If primer is not found and no trimming is done, both R1 and R2 are discarded
## temporary, log, and information files are placed in the sample directory... 
${CUTADAPT} --no-indels -O ${OFWD} -g ${FBC} -e ${ERROR} --info-file ./Clipped/${SGE_TASK_ID}"_clip_fr.R1.info" -o ./Clipped/${SGE_TASK_ID}"_clip_R1.TMP.fastq" -p ./Clipped/${SGE_TASK_ID}"_clip_R2.TMP.fastq" ./Renamed/${SGE_TASK_ID}"_R1.fastq" ./Renamed/${SGE_TASK_ID}"_R2.fastq" --untrimmed-o /dev/null --untrimmed-p /dev/null > ./Clipped/${SGE_TASK_ID}"_clip_fr.R1.cutadapt.log" 2>&1
  
## Using the output of the previous command (the TMP files), we'll now look for the reverse primer in R2, bringing R1 along for the ride.
## note inverse orientation of input files. As before, if primer is not found, remove the pair.
## this step generates the final results for the FWD-REV orienatation, accepting those seqs that passed the first stage.
${CUTADAPT} --no-indels -O ${OREV} -g ${RBC} -e ${ERROR} --info-file ./Clipped/${SGE_TASK_ID}"_clip_fr.R2.info" -o ./Clipped/${SGE_TASK_ID}"_clip_R2.fastq" -p ./Clipped/${SGE_TASK_ID}"_clip_R1.fastq" ./Clipped/${SGE_TASK_ID}"_clip_R2.TMP.fastq" ./Clipped/${SGE_TASK_ID}"_clip_R1.TMP.fastq" --untrimmed-o /dev/null --untrimmed-p /dev/null > ./Clipped/${SGE_TASK_ID}"_clip_fr.R2.cutadapt.log" 2>&1

## clean up TMP files...
rm ./Clipped/${SGE_TASK_ID}"_clip_R1.TMP.fastq" 
rm ./Clipped/${SGE_TASK_ID}"_clip_R2.TMP.fastq"

# Now to process the reads with the RF orientation...
## First, we search for the reverse primer in R1...
${CUTADAPT} --no-indels -O ${OREV} -g ${RBC} -e ${ERROR} --info-file ./Clipped/${SGE_TASK_ID}"_clip_rf.R1.info" -o ./Clipped/${SGE_TASK_ID}"_clip_R1.TMP.fastq" -p ./Clipped/${SGE_TASK_ID}"_clip_R2.TMP.fastq" ./Renamed/${SGE_TASK_ID}"_R1.fastq" ./Renamed/${SGE_TASK_ID}"_R2.fastq" --untrimmed-o /dev/null --untrimmed-p /dev/null > ./Clipped/${SGE_TASK_ID}"_clip_rf.R1.cutadapt.log" 2>&1
  
## As before, we search for the forward primer in R2, only processing the output of the previous command
${CUTADAPT} --no-indels -O ${OFWD} -g ${FBC} -e ${ERROR} --info-file ./Clipped/${SGE_TASK_ID}"_clip_rf.R2.info" -o ./Clipped/${SGE_TASK_ID}"_clip_rf.R2.fastq" -p ./Clipped/${SGE_TASK_ID}"_clip_rf.R1.fastq" ./Clipped/${SGE_TASK_ID}"_clip_R2.TMP.fastq" ./Clipped/${SGE_TASK_ID}"_clip_R1.TMP.fastq" --untrimmed-o /dev/null --untrimmed-p /dev/null > ./Clipped/${SGE_TASK_ID}"_clip_rf.R2.cutadapt.log" 2>&1

## remove the temp files, as -o ${RFR2} -p ${RFR1} have this step's results
rm ./Clipped/${SGE_TASK_ID}"_clip_R1.TMP.fastq" 
rm ./Clipped/${SGE_TASK_ID}"_clip_R2.TMP.fastq"

# Reorient REV-FWD output to FWD-REV
# Change the read id in the headers "@MISEQ:41:000000000-A9A9U:1:1101:17488:1966 1:N:0:AGGCAGAAAGAGTAGA"
awk '{if (NR%4==1){gsub("^1:","2:",$2); print $0}else{print $0}}' ./Clipped/${SGE_TASK_ID}"_clip_rf.R1.fastq" > ./Clipped/${SGE_TASK_ID}"_rf2frR2.fastq"
awk '{if (NR%4==1){gsub("^2:","1:",$2); print $0}else{print $0}}' ./Clipped/${SGE_TASK_ID}"_clip_rf.R2.fastq" > ./Clipped/${SGE_TASK_ID}"_rf2frR1.fastq"
    
# rename the reorient. rev-fwd files
mv ./Clipped/${SGE_TASK_ID}"_rf2frR1.fastq" ./Clipped/${SGE_TASK_ID}"_clip_rf.R1.fastq"
mv ./Clipped/${SGE_TASK_ID}"_rf2frR2.fastq" ./Clipped/${SGE_TASK_ID}"_clip_rf.R2.fastq"

# add all corrected seqs to the fwd-rev fastqs
cat ./Clipped/${SGE_TASK_ID}"_clip_rf.R1.fastq" >> ./Clipped/${SGE_TASK_ID}"_clip_R1.fastq"
cat ./Clipped/${SGE_TASK_ID}"_clip_rf.R2.fastq" >> ./Clipped/${SGE_TASK_ID}"_clip_R2.fastq"


