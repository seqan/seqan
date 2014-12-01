#!/bin/bash

MASON_GENOME=../../../../../seqan-clean-build/release/bin/mason_genome
MASON_SIMULATOR=../../../../../seqan-clean-build/release/bin/mason_simulator
FIONA=../../../../../seqan-clean-build/release/bin/fiona
FIONA_ILLUMINA=../../../../../seqan-clean-build/release/bin/fiona_illumina

# ============================================================
# Simulate Genome and reads.
# ============================================================

echo "${MASON_GENOME} -l 10000 -o genome.10k.fa &>/dev/null"
${MASON_GENOME} -l 10000 -o genome.10k.fa &>/dev/null

echo "${MASON_SIMULATOR} -ir genome.10k.fa -n 300 -o reads.illumina.fq -oa reads.illumina.sam &>/dev/null"
${MASON_SIMULATOR} -ir genome.10k.fa -n 300 -o reads.illumina.fq -oa reads.illumina.sam &>/dev/null

echo "${MASON_SIMULATOR} --seq-technology 454 --fragment-mean-size 800 -ir genome.10k.fa -n 300 -o reads.454.fq -oa reads.454.sam &>/dev/null"
${MASON_SIMULATOR} --seq-technology 454 --fragment-mean-size 800 -ir genome.10k.fa -n 300 -o reads.454.fq -oa reads.454.sam &>/dev/null

# ============================================================
# Run Fiona.
# ============================================================

# Illumina Mode

echo "${FIONA_ILLUMINA} -nt 1 -i 1 -g 10000 reads.illumina.fq reads.illumina.corrected.i1.fa >reads.illumina.fq.i1.stdout 2>reads.illumina.fq.i1.stderr"
${FIONA_ILLUMINA} -nt 1 -i 1 -g 10000 reads.illumina.fq reads.illumina.corrected.i1.fa >reads.illumina.fq.i1.stdout 2>reads.illumina.fq.i1.stderr
echo "  => $?"

echo "${FIONA_ILLUMINA} -nt 1 -i 2 -g 10000 reads.illumina.fq reads.illumina.corrected.i2.fa >reads.illumina.fq.i2.stdout 2>reads.illumina.fq.i2.stderr"
${FIONA_ILLUMINA} -nt 1 -i 2 -g 10000 reads.illumina.fq reads.illumina.corrected.i2.fa >reads.illumina.fq.i2.stdout 2>reads.illumina.fq.i2.stderr
echo "  => $?"

# 454 Mode

echo "${FIONA} -nt 1 -i 1 -g 10000 reads.454.fq reads.454.corrected.i1.fa >reads.454.fq.i1.stdout 2>reads.454.fq.i1.stderr"
${FIONA} -nt 1 -i 1 -g 10000 reads.454.fq reads.454.corrected.i1.fa >reads.454.fq.i1.stdout 2>reads.454.fq.i1.stderr
echo "  => $?"

echo "${FIONA} -nt 1 -i 2 -g 10000 reads.454.fq reads.454.corrected.i2.fa >reads.454.fq.i2.stdout 2>reads.454.fq.i2.stderr"
${FIONA} -nt 1 -i 2 -g 10000 reads.454.fq reads.454.corrected.i2.fa >reads.454.fq.i2.stdout 2>reads.454.fq.i2.stderr
echo "  => $?"
