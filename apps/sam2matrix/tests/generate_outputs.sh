#!/bin/sh
#
# Output generation script for sam2matrix

SAM2MATRIX=../../../../build/Debug/apps/sam2matrix/sam2matrix

# ============================================================
# First Section
# ============================================================

${SAM2MATRIX} -m ecoli.sam -m ehec.sam -r ecoli_0.50_ehec_0.50.fq -rf
ecoli.fa -rf ehec.fa -o test_sam2matrix.tsv > out.stdout
