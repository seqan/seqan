#!/bin/sh
#
# The files adeno{1,2,3}.fa are fa files with each two sequences, each a
# different pair.  They contain the first 280 base pairs of Adenoviruses.
#
# The files adeno{1,2,3}-rna.fa are the same as adeno{1,2,3}.fa but all
# occurrences of T have been replaced by a U.
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
    echo ${PAIR_ALIGN} -s ${file}.fa -o ${file}_out.fa
    ${PAIR_ALIGN} -s ${file}.fa -o ${file}_out.fa > ${file}.stdout
done

# Run with explicit alphabet.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -a protein -s ${file}.fa -o ${file}.protein_out.fa
    ${PAIR_ALIGN} -a protein -s ${file}.fa -o ${file}.protein_out.fa > ${file}.protein.stdout
done

# Run with different alignment methods.
for file in 1aab 1ad2 2trx; do
    for m in nw gotoh sw lcs; do
        echo ${PAIR_ALIGN} -m ${m} -s ${file}.fa -o ${file}.m${m}_out.fa
        ${PAIR_ALIGN} -m ${m} -s ${file}.fa -o ${file}.m${m}_out.fa > ${file}.m${m}.stdout
    done
done

# Run with different scoring options.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -g -20 -s ${file}.fa -o ${file}.g-20_out.fa
    ${PAIR_ALIGN} -g -20 -s ${file}.fa -o ${file}.g-20_out.fa > ${file}.g-20.stdout
    echo ${PAIR_ALIGN} -e -5 -s ${file}.fa -o ${file}.e-5_out.fa
    ${PAIR_ALIGN} -e -5 -s ${file}.fa -o ${file}.e-5_out.fa > ${file}.e-5.stdout
    echo ${PAIR_ALIGN} -ms 10 -s ${file}.fa -o ${file}.ms10_out.fa
    ${PAIR_ALIGN} -ms 10 -s ${file}.fa -o ${file}.ms10_out.fa > ${file}.ms10.stdout
    echo ${PAIR_ALIGN} -mm -8 -s ${file}.fa -o ${file}.mm-8_out.fa
    ${PAIR_ALIGN} -mm -8 -s ${file}.fa -o ${file}.mm-8_out.fa > ${file}.mm-8.stdout
done

# Run with matrix file.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200_out.fa
    ${PAIR_ALIGN} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200_out.fa > ${file}.maVTML200.stdout
done

# Run with different banded alignment options.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -lo 5 -s ${file}.fa -o ${file}.lo5_out.fa
    ${PAIR_ALIGN} -lo 5 -s ${file}.fa -o ${file}.lo5_out.fa > ${file}.lo5.stdout
    echo ${PAIR_ALIGN} -hi 5 -s ${file}.fa -o ${file}.hi5_out.fa
    ${PAIR_ALIGN} -hi 5 -s ${file}.fa -o ${file}.hi5_out.fa > ${file}.hi5.stdout
done

# Run with different matrix configuration options
for file in 1aab 1ad2 2trx; do
    for c in ffff tttt ffft fftf ftff tfff fftt fttf ttff tfft; do
        echo ${PAIR_ALIGN} -c ${c} -s ${file}.fa -o ${file}.c${c}_out.fa
        ${PAIR_ALIGN} -c ${c} -s ${file}.fa -o ${file}.c${c}_out.fa > ${file}.c${c}.stdout
    done
done

# ============================================================
# Run on DNA (Adenoviruses).
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 1 2 3; do
    echo ${PAIR_ALIGN} -a dna -s adeno${i}.fa -o adeno${i}_out.fa
    ${PAIR_ALIGN} -a dna -s adeno${i}.fa -o adeno${i}_out.fa > adeno${i}.stdout
done

# ============================================================
# Run on RNA.
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 1 2 3; do
    echo ${PAIR_ALIGN} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna_out.fa
    ${PAIR_ALIGN} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna_out.fa > adeno${i}-rna.stdout
done
