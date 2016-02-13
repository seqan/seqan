#!/bin/sh
#
# Output generation script for seqFlexQC

BASE=../../../../build/Debug/apps/seqan_flexbar/

# ============================================================
# First Section
# ============================================================

${BASE}sflexQC
../../../seqan-git/apps/seqan_flexbar/tests/testsample.fq -q 20
-o qc_test.fa -t -ni > qc_test.stdout

${BASE}sflexFilter
../../../seqan-git/apps/seqan_flexbar/tests/testsample.fq -tl 3 -tr 4 -ml 70 -u 1 -s A -fl 70 -ni -o filter_test.fq > filter_test.stdout

${BASE}sflexAR
./../../seqan-git/apps/seqan_flexbar/tests/testsample.fq -a ../../../seqan-git/apps/seqan_flexbar/tests/adapter.fa -o ar_test.fq -ni > ar_test.stdout

${BASE}sflexDMulti
../../../seqan-git/apps/seqan_flexbar/tests/testsample_multiplex.fq
-b ../../../seqan-git/apps/seqan_flexbar/tests/barcodes.fa -o
test_de_multi.fq -ni > test_de_multi.stdout

