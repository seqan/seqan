#!/bin/sh
#
# Output generation script for sam2matrix

SAM2MATRIX=../../../../build/Debug/extras/apps/sam2matrix/sam2matrix

# ============================================================
# First Section
# ============================================================

${SAM2MATRIX} -sf ecoli.sam -sf ehec.sam -rf ecoli_0.50_ehec_0.50.fastq -gf ecoli.fasta -gf ehec.fasta -o test_sam2matrix.csv > out.stdout
