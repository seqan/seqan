#!/bin/sh
#
# We slice and dice the file "adeno.fasta".

# sak build from 2010-05-25
#SAK=../../seqan-trunk-build/Debug/core/apps/sak/sak -ll 1000
SAK=../../../../../seqan-trunk-build/Debug/core/apps/sak/sak

# ============================================================
# Run on DNA (Adenoviruses).
# ============================================================

echo ${SAK} adeno.fasta -o adeno.all.out
${SAK} adeno.fasta -o adeno.all.out

echo ${SAK} adeno.fasta -s 1 -o adeno.seq1.out
${SAK} adeno.fasta -s 1 -o adeno.seq1.out

echo ${SAK} adeno.fasta -ss 1-2 -o adeno.seq1-2.out
${SAK} adeno.fasta -ss 1-2 -o adeno.seq1-2.out

echo ${SAK} adeno.fasta -s 3 -o adeno.seq3.out
${SAK} adeno.fasta -s 3 -o adeno.seq3.out

echo ${SAK} adeno.fasta -sn 'gi|9626621' -o adeno.sn.out
${SAK} adeno.fasta -sn 'gi|9626621' -o adeno.sn.out

echo ${SAK} adeno.fasta -s 1 -i 5-25 -o adeno.s1i5-25.out
${SAK} adeno.fasta -s 1 -i 5-25 -o adeno.s1i5-25.out

echo ${SAK} adeno.fasta -ss 1-2 -i 5-25 -o adeno.s1-2i5-25.out
${SAK} adeno.fasta -ss 1-2 -i 5-25 -o adeno.s1-2i5-25.out

echo ${SAK} adeno.fasta -s 1 -rc -o adeno.s1rc.out
${SAK} adeno.fasta -s 1 -rc -o adeno.s1rc.out
