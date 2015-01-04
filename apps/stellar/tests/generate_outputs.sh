#!/bin/sh
#
# Output generation script for stellar

STELLAR=../../../../build/release/bin/stellar

# ============================================================
# Varying error rates
# ============================================================

eps="e-1"
errRate=0.1
${STELLAR} -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -o $eps.gff 512_simSeq1_$eps.fa 512_simSeq2_$eps.fa > $eps.stdout

eps="75e-3"
errRate=0.075
${STELLAR} -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -o $eps.gff 512_simSeq1_$eps.fa 512_simSeq2_$eps.fa > $eps.stdout

eps="5e-2"
errRate=0.05
${STELLAR} -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -o $eps.gff 512_simSeq1_$eps.fa 512_simSeq2_$eps.fa > $eps.stdout

eps="25e-3"
errRate=0.025
${STELLAR} -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -o $eps.gff 512_simSeq1_$eps.fa 512_simSeq2_$eps.fa > $eps.stdout

eps="e-4"
errRate=0.0001
${STELLAR} -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -o $eps.gff 512_simSeq1_$eps.fa 512_simSeq2_$eps.fa > $eps.stdout


# ============================================================
# Varying minimal lengths
# ============================================================

minLen="20"
${STELLAR} -e 0.05 -l $minLen -x 10 -k 7 -n 5000 -s 10000 -f -v -t -o minLen$minLen.gff 512_simSeq1_5e-2.fa 512_simSeq2_5e-2.fa > minLen$minLen.stdout

minLen="150"
${STELLAR} -e 0.05 -l $minLen -x 10 -k 7 -n 5000 -s 10000 -f -v -t -o minLen$minLen.gff 512_simSeq1_5e-2.fa 512_simSeq2_5e-2.fa > minLen$minLen.stdout


# ============================================================
# Output format
# ============================================================

eps="5e-2"
errRate=0.05
${STELLAR} -e $errRate -l 50 -x 10 -k 7 -n 5000 -s 10000 -f -v -t -o $eps.txt 512_simSeq1_$eps.fa 512_simSeq2_$eps.fa > $eps"txt.stdout"  
