#!/bin/sh
#
# Output generation script for insegt

INSEGT=../../../../../build/Debug/bin/insegt

# ============================================================
# First Section
# ============================================================

${INSEGT} -ro default_readOutput.gff -ao default_annoOutput.gff -to default_tupleOutput.gff alignments.sam annotations.gff

${INSEGT} -c 2 -ro threshold-count2_readOutput.gff -ao threshold-count2_annoOutput.gff -to threshold-count2_tupleOutput.gff alignments.sam annotations.gff

${INSEGT} -n 3 -ro ntuple3_readOutput.gff -ao ntuple3_annoOutput.gff -to ntuple3_tupleOutput.gff alignments.sam annotations.gff

${INSEGT} -m -ro max-tuple_readOutput.gff -ao max-tuple_annoOutput.gff -to max-tuple_tupleOutput.gff alignments.sam annotations.gff

${INSEGT} -e  -ro exact-ntuple_readOutput.gff -ao exact-ntuple_annoOutput.gff -to exact-ntuple_tupleOutput.gff alignments.sam annotations.gff


