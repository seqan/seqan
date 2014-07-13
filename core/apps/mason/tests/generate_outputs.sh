#!/bin/sh
#
# Expected output generation for mason.

# We use version 13383 for building the test data.
MASON=${MASON:-../../../../../seqan-align-build/release/bin/mason}

# ============================================================
# Simulate 454 Reads, Single-End
# ============================================================

${MASON} 454 -N 100 -nm 400 -ne 40 -s 0 -rnp read -o 454-se-random-N100-nm400-ne40-s0.fasta random.fasta > 454-se-random-N100-nm400-ne40-s0.stdout
${MASON} 454 -N 100 -nm 400 -ne 40 -s 0 -sq -rnp read -o 454-se-random-N100-nm400-ne40-s0-sq.fastq random.fasta > 454-se-random-N100-nm400-ne40-s0-sq.stdout
${MASON} 454 -N 100 -nm 200 -ne 20 -s 0 -rnp read -o 454-se-random-N100-nm200-ne20-s0.fasta random.fasta > 454-se-random-N100-nm200-ne20-s0.stdout
${MASON} 454 -N 100 -nm 200 -ne 20 -s 0 -sq -rnp read -o 454-se-random-N100-nm200-ne20-s0-sq.fastq random.fasta > 454-se-random-N100-nm200-ne20-s0-sq.stdout

# ============================================================
# Simulate Illumina Reads, Single-End, Existing Reference
# ============================================================

${MASON} illumina -N 100 -n 36 -s 0 -rnp read -o illumina-se-adeno-N100-n36-s0.fasta adeno-genome.fa > illumina-se-adeno-N100-n36-s0.stdout
${MASON} illumina -N 100 -n 36 -s 0 -sq -rnp read -o illumina-se-adeno-N100-n36-s0-sq.fastq adeno-genome.fa > illumina-se-adeno-N100-n36-s0-sq.stdout
${MASON} illumina -N 100 -n 100 -s 0 -rnp read -o illumina-se-adeno-N100-n100-s0.fasta adeno-genome.fa > illumina-se-adeno-N100-n100-s0.stdout
${MASON} illumina -N 100 -n 100 -s 0 -sq -rnp read -o illumina-se-adeno-N100-n100-s0-sq.fastq adeno-genome.fa > illumina-se-adeno-N100-n100-s0-sq.stdout

# ============================================================
# Simulate Illumina Reads, Single-End
# ============================================================

${MASON} illumina -N 100 -n 36 -s 0 -rnp read -o illumina-se-random-N100-n36-s0.fasta random.fasta > illumina-se-random-N100-n36-s0.stdout
${MASON} illumina -N 100 -n 36 -s 0 -sq -rnp read -o illumina-se-random-N100-n36-s0-sq.fastq random.fasta > illumina-se-random-N100-n36-s0-sq.stdout
${MASON} illumina -N 100 -n 100 -s 0 -rnp read -o illumina-se-random-N100-n100-s0.fasta random.fasta > illumina-se-random-N100-n100-s0.stdout
${MASON} illumina -N 100 -n 100 -s 0 -sq -rnp read -o illumina-se-random-N100-n100-s0-sq.fastq random.fasta > illumina-se-random-N100-n100-s0-sq.stdout

# ============================================================
# Simulate Illumina Reads, Paired-End
# ============================================================

${MASON} illumina -mp -N 100 -n 36 -s 0 -rnp read -o illumina-pe-random-N100-n36-s0.fasta random.fasta > illumina-pe-random-N100-n36-s0.stdout
${MASON} illumina -mp -N 100 -n 36 -s 0 -sq -rnp read -o illumina-pe-random-N100-n36-s0-sq.fastq random.fasta > illumina-pe-random-N100-n36-s0-sq.stdout
${MASON} illumina -mp -N 100 -n 100 -s 0 -rnp read -o illumina-pe-random-N100-n100-s0.fasta random.fasta > illumina-pe-random-N100-n100-s0.stdout
${MASON} illumina -mp -N 100 -n 100 -s 0 -sq -rnp read -o illumina-pe-random-N100-n100-s0-sq.fastq random.fasta > illumina-pe-random-N100-n100-s0-sq.stdout

# ============================================================
# Simulate Sanger Reads, Single-End
# ============================================================

${MASON} sanger -N 100 -nm 400 -ne 40 -s 0 -rnp read -o sanger-se-random-N100-nm400-ne40-s0.fasta random.fasta > sanger-se-random-N100-nm400-ne40-s0.stdout
${MASON} sanger -N 100 -nm 400 -ne 40 -s 0 -sq -rnp read -o sanger-se-random-N100-nm400-ne40-s0-sq.fastq random.fasta > sanger-se-random-N100-nm400-ne40-s0-sq.stdout
${MASON} sanger -N 100 -nm 200 -ne 20 -s 0 -rnp read -o sanger-se-random-N100-nm200-ne20-s0.fasta random.fasta > sanger-se-random-N100-nm200-ne20-s0.stdout
${MASON} sanger -N 100 -nm 200 -ne 20 -s 0 -sq -rnp read -o sanger-se-random-N100-nm200-ne20-s0-sq.fastq random.fasta > sanger-se-random-N100-nm200-ne20-s0-sq.stdout
