#!/bin/bash

while getopts s:t:d:n:x option
do
 case "${option}"
in
s) SPACERFILE=${OPTARG};;
t) TARGETFILE=${OPTARG};;
d) FILEDATE=${OPTARG};;
n) THREADS=${OPTARG};;
x) NONE=${OPTARG};;
esac
done

echo "Converting spacer search file : ${PROJDIR}/output/${SPACERFILE}"
# Convert 2 column fasta-like list ex R to fasta single column
	awk '{for(i=1;i<=NF;i++) printf "%s\n",$i}' ${PROJDIR}/output/${SPACERFILE} > ${PROJDIR}/tmp/spacers_to_search.fasta

# Run GASSST
echo "Running GASSST for ${SPACERFILE} against ${TARGETFILE} target database"
Gassst -d ${PROJDIR}/phage_sequences/${TARGETFILE} -i ${PROJDIR}/tmp/spacers_to_search.fasta \
-o ${PROJDIR}/output/search_results_${SPACERFILE}_${TARGETFILE}.txt -p 80 -w 7 -n ${THREADS} -t 200 -g 5 -s 4

#rm ${PROJDIR}/tmp/spacers_to_search.fasta

# remove the header lines
grep -v "@" ${PROJDIR}/output/search_results_${SPACERFILE}_${TARGETFILE}.txt > ${PROJDIR}/tmp/tmp.txt
cp ${PROJDIR}/tmp/tmp.txt ${PROJDIR}/output/search_results_${SPACERFILE}_${TARGETFILE}.txt

rm ${PROJDIR}/tmp/tmp.txt

echo "Search done"
# Shuffle target database using MEME on the server...
echo "Run against the shuffled database"

# Run shuffled dataset for comparison....
echo "Running GASSST for ${SPACERFILE} against shuffled ${TARGETFILE} target database"
Gassst -d ${PROJDIR}/phage_sequences/shuffled_${TARGETFILE} -i ${PROJDIR}/tmp/spacers_to_search.fasta \
-o ${PROJDIR}/output/search_results_${SPACERFILE}_${TARGETFILE}_shuffled.txt -p 80 -w 7 -n ${THREADS} -t 200 -g 5 -s 4

rm ${PROJDIR}/tmp/spacers_to_search.fasta

# remove the header lines
grep -v "@" ${PROJDIR}/output/search_results_${SPACERFILE}_${TARGETFILE}_shuffled.txt > ${PROJDIR}/tmp/tmp.txt
cp ${PROJDIR}/tmp/tmp.txt ${PROJDIR}/output/search_results_${SPACERFILE}_${TARGETFILE}_shuffled.txt

rm ${PROJDIR}/tmp/tmp.txt

echo "Shuffled search done - all complete"








