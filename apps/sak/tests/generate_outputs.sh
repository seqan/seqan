#!/bin/sh

# We slice and dice the file "adeno.fa".

# sak build from 2010-05-25
#SAK=../../seqan-trunk-build/Debug/bin/sak -ll 1000
SAK=../../../../../seqan-git-build/debug/bin/sak

# ============================================================
# Run on DNA (Adenoviruses).
# ============================================================

echo ${SAK} adeno.fa -o adeno.all.fa
${SAK} adeno.fa -o adeno.all.fa

echo ${SAK} adeno.fa -s 1 -o adeno.seq1.fa
${SAK} adeno.fa -s 1 -o adeno.seq1.fa

echo ${SAK} adeno.fa -ss 1-2 -o adeno.seq1-2.fa
${SAK} adeno.fa -ss 1-2 -o adeno.seq1-2.fa

echo ${SAK} adeno.fa -s 3 -o adeno.seq3.fa
${SAK} adeno.fa -s 3 -o adeno.seq3.fa

echo ${SAK} adeno.fa -sn 'gi|9626621' -o adeno.sn.fa
${SAK} adeno.fa -sn 'gi|9626621' -o adeno.sn.fa

echo ${SAK} adeno.fa -s 1 -i 5-25 -o adeno.s1i5-25.fa
${SAK} adeno.fa -s 1 -i 5-25 -o adeno.s1i5-25.fa

echo ${SAK} adeno.fa -ss 1-2 -i 5-25 -o adeno.s1-2i5-25.fa
${SAK} adeno.fa -ss 1-2 -i 5-25 -o adeno.s1-2i5-25.fa

echo ${SAK} adeno.fa -s 1 -rc -o adeno.s1rc.fa
${SAK} adeno.fa -s 1 -rc -o adeno.s1rc.fa
