#!/bin/sh
#
# Output generation script for seqFlexQC

SEQFLEXQC=../../../../build/Debug/extras/apps/seqan_flexbar/seqFlexQC
SEQFLEXFILTER=../../../../build/Debug/extras/apps/seqan_flexbar/seqFlexFilter
SEQFLEXAR=../../../../build/Debug/extras/apps/seqan_flexbar/seqFlexAR

# ============================================================
# First Section
# ============================================================

${SEQFLEXQC}
../../../seqan-git/extras/apps/seqan_flexbar/tests/testsample.fq -q 20
-out qc_test.fa -t -ni > qc_test.stdout

${SEQFLEXFILTER}
../../../seqan-git/extras/apps/seqan_flexbar/tests/testsample.fq -tl 3 -tr 4 -ml 70 -u 1 -s A -fl 70 -ni -out filter_test.fq > filter_test.stdout

${SEQFLEXAR}
./../../seqan-git/extras/apps/seqan_flexbar/tests/testsample.fq -a ../../../seqan-git/extras/apps/seqan_flexbar/tests/adapter.fasta -out ar_test.fq -ni > ar_test.stdout
