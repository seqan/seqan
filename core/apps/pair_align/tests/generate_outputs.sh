#!/bin/sh
#
# The files adeno{1,2,3}.fa are fasta files with each two sequences, each a
# different pair.  They contain the first 280 base pairs of Adenoviruses.
#
# The files adeno{1,2,3}-rna.fa are the same as adeno{1,2,3}.fa but all
# occurences of T have been replaced by a U.
#
# The files 1ad2.fa, 1aab.fa and 2trx.fa contain pairs from the Balibase
# benchmark.
#
# The script will run the pair_align program as configured
# below with different parameters on all files and generate output
# files for them.
#
# They can be used to compare against the output of new pair_align programs.

# Revision 13682 was used for generating the output.
PAIR_ALIGN=../../../../../build/release/bin/pair_align.exe

# ============================================================
# Run on Proteins (Balibase).
# ============================================================

# Run with defaults for all non-mandatory options.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -s ${file}.fa -o ${file}_out.fasta
    ${PAIR_ALIGN} -s ${file}.fa -o ${file}_out.fasta > ${file}.stdout
done

# Run with explicit alphabet.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -a protein -s ${file}.fa -o ${file}.protein_out.fasta
    ${PAIR_ALIGN} -a protein -s ${file}.fa -o ${file}.protein_out.fasta > ${file}.protein.stdout
done

# Run with different alignment methods.
for file in 1aab 1ad2 2trx; do
    for m in nw gotoh sw lcs; do
        echo ${PAIR_ALIGN} -m ${m} -s ${file}.fa -o ${file}.m${m}_out.fasta
        ${PAIR_ALIGN} -m ${m} -s ${file}.fa -o ${file}.m${m}_out.fasta > ${file}.m${m}.stdout
    done
done

# Run with different scoring options.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -g -20 -s ${file}.fa -o ${file}.g-20_out.fasta
    ${PAIR_ALIGN} -g -20 -s ${file}.fa -o ${file}.g-20_out.fasta > ${file}.g-20.stdout
    echo ${PAIR_ALIGN} -e -5 -s ${file}.fa -o ${file}.e-5_out.fasta
    ${PAIR_ALIGN} -e -5 -s ${file}.fa -o ${file}.e-5_out.fasta > ${file}.e-5.stdout
    echo ${PAIR_ALIGN} -ms 10 -s ${file}.fa -o ${file}.ms10_out.fasta
    ${PAIR_ALIGN} -ms 10 -s ${file}.fa -o ${file}.ms10_out.fasta > ${file}.ms10.stdout
    echo ${PAIR_ALIGN} -mm -8 -s ${file}.fa -o ${file}.mm-8_out.fasta
    ${PAIR_ALIGN} -mm -8 -s ${file}.fa -o ${file}.mm-8_out.fasta > ${file}.mm-8.stdout
done

# Run with matrix file.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200_out.fasta
    ${PAIR_ALIGN} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200_out.fasta > ${file}.maVTML200.stdout
done

# Run with different banded alignment options.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -lo 5 -s ${file}.fa -o ${file}.lo5_out.fasta
    ${PAIR_ALIGN} -lo 5 -s ${file}.fa -o ${file}.lo5_out.fasta > ${file}.lo5.stdout
    echo ${PAIR_ALIGN} -hi 5 -s ${file}.fa -o ${file}.hi5_out.fasta
    ${PAIR_ALIGN} -hi 5 -s ${file}.fa -o ${file}.hi5_out.fasta > ${file}.hi5.stdout
done

# Run with different matrix configuration options
for file in 1aab 1ad2 2trx; do
    for c in ffff tttt ffft fftf ftff tfff fftt fttf ttff tfft; do
        echo ${PAIR_ALIGN} -c ${c} -s ${file}.fa -o ${file}.c${c}_out.fasta
        ${PAIR_ALIGN} -c ${c} -s ${file}.fa -o ${file}.c${c}_out.fasta > ${file}.c${c}.stdout
    done
done

# ============================================================
# Run on DNA (Adenoviruses).
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 1 2 3; do
    echo ${PAIR_ALIGN} -a dna -s adeno${i}.fa -o adeno${i}_out.fasta
    ${PAIR_ALIGN} -a dna -s adeno${i}.fa -o adeno${i}_out.fasta > adeno${i}.stdout
done

# ============================================================
# Run on RNA.
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 1 2 3; do
    echo ${PAIR_ALIGN} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna_out.fasta
    ${PAIR_ALIGN} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna_out.fasta > adeno${i}-rna.stdout
done
