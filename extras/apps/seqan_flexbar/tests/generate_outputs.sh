#!/bin/sh
#
# Output generation script for seqFlexQC

BASE=../../../../build/Debug/extras/apps/seqan_flexbar/

# ============================================================
# First Section
# ============================================================

${BASE}sflexQC
../../../seqan-git/extras/apps/seqan_flexbar/tests/testsample.fq -q 20
-o qc_test.fa -t -ni > qc_test.stdout

${BASE}sflexFilter
../../../seqan-git/extras/apps/seqan_flexbar/tests/testsample.fq -tl 3 -tr 4 -ml 70 -u 1 -s A -fl 70 -ni -o filter_test.fq > filter_test.stdout

${BASE}sflexAR
./../../seqan-git/extras/apps/seqan_flexbar/tests/testsample.fq -a ../../../seqan-git/extras/apps/seqan_flexbar/tests/adapter.fasta -o ar_test.fq -ni > ar_test.stdout
