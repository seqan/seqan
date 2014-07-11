#!/bin/sh
#
SAMCAT=../../../../build/clang/bin/samcat

# ============================================================
# Run on DNA (Adenoviruses).
# ============================================================

echo ${SAMCAT} ex1_a1.sam ex1_a2.sam ex1_a3.sam -o ex1_merged.sam
${SAMCAT} ex1_a1.sam ex1_a2.sam ex1_a3.sam -o ex1_merged.sam

echo ${SAMCAT} ex1_a1.sam ex1_a2.sam ex1_a3.sam -o ex1_merged.bam
${SAMCAT} ex1_a1.sam ex1_a2.sam ex1_a3.sam -o ex1_merged.bam

