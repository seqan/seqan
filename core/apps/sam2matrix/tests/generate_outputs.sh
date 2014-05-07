#!/bin/sh
#
# Output generation script for sam2matrix

SAM2MATRIX=../../../../build/Debug/core/apps/sam2matrix/sam2matrix

# ============================================================
# First Section
# ============================================================

${SAM2MATRIX} -m ecoli.sam -m ehec.sam -r ecoli_0.50_ehec_0.50.fastq -rf
ecoli.fasta -rf ehec.fasta -o test_sam2matrix.csv > out.stdout
