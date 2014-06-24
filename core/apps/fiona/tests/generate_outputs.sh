#!/bin/sh
#
FIONA=../../../../build/gcc47/bin/fiona

# ============================================================
# Run on DNA (Adenoviruses).
# ============================================================

#echo ${FIONA} -g 10000 ../example/simulatedReads_Length75.fasta simulatedReads_corrected.fasta
#${FIONA} -nt 4 -v --super-packages 1 -i 2 -g 10000 ../example/simulatedReads_Length75.fasta simulatedReads_corrected.fasta

echo ${FIONA} -g 10000 ../reads.fa corrected.fa
${FIONA} -nt 4 -v -i 1 -g 10000 --correction-infos ../example/reads.fa corrected.fa

