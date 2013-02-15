#!/bin/sh
#
# Output generation script for insegt

INSEGT=../../../../../build/Debug/bin/insegt

# ============================================================
# First Section
# ============================================================

${INSEGT} -p default_ alignments.sam annotations.gff

${INSEGT} -c 2 -p threshold-count2_ alignments.sam annotations.gff

${INSEGT} -n 3 -p ntuple3_ alignments.sam annotations.gff

${INSEGT} -m -p max-tuple_ alignments.sam annotations.gff

${INSEGT} -e  -p exact-ntuple_ alignments.sam annotations.gff


