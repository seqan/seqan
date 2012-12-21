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

# Revision 13383 was used for generating the output.
PAIR_ALIGN=../../../../../seqan-align-build/release/bin/pair_align

# ============================================================
# Run on Proteins (Balibase).
# ============================================================

# Run with defaults for all non-mandatory options.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -s ${file}.fa -o ${file}.out
    ${PAIR_ALIGN} -s ${file}.fa -o ${file}.out > ${file}.stdout
done

# Run with explicit alphabet.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -a protein -s ${file}.fa -o ${file}.protein.out
    ${PAIR_ALIGN} -a protein -s ${file}.fa -o ${file}.protein.out > ${file}.protein.stdout
done

# Run with different alignment methods.
for file in 1aab 1ad2 2trx; do
    for m in nw gotoh sw lcs; do
        echo ${PAIR_ALIGN} -m ${m} -s ${file}.fa -o ${file}.m${m}.out
        ${PAIR_ALIGN} -m ${m} -s ${file}.fa -o ${file}.m${m}.out > ${file}.m${m}.stdout
    done
done

# Run with different scoring options.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -g -20 -s ${file}.fa -o ${file}.g-20.out
    ${PAIR_ALIGN} -g -20 -s ${file}.fa -o ${file}.g-20.out > ${file}.g-20.stdout
    echo ${PAIR_ALIGN} -e -5 -s ${file}.fa -o ${file}.e-5.out
    ${PAIR_ALIGN} -e -5 -s ${file}.fa -o ${file}.e-5.out > ${file}.e-5.stdout
    echo ${PAIR_ALIGN} -ms 10 -s ${file}.fa -o ${file}.ms10.out
    ${PAIR_ALIGN} -ms 10 -s ${file}.fa -o ${file}.ms10.out > ${file}.ms10.stdout
    echo ${PAIR_ALIGN} -mm -8 -s ${file}.fa -o ${file}.mm-8.out
    ${PAIR_ALIGN} -mm -8 -s ${file}.fa -o ${file}.mm-8.out > ${file}.mm-8.stdout
done

# Run with matrix file.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200.out
    ${PAIR_ALIGN} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200.out > ${file}.maVTML200.stdout
done

# Run with different banded alignment options.
for file in 1aab 1ad2 2trx; do
    echo ${PAIR_ALIGN} -lo 5 -s ${file}.fa -o ${file}.lo5.out
    ${PAIR_ALIGN} -lo 5 -s ${file}.fa -o ${file}.lo5.out > ${file}.lo5.stdout
    echo ${PAIR_ALIGN} -hi 5 -s ${file}.fa -o ${file}.hi5.out
    ${PAIR_ALIGN} -hi 5 -s ${file}.fa -o ${file}.hi5.out > ${file}.hi5.stdout
done

# Run with different matrix configuration options
for file in 1aab 1ad2 2trx; do
    for c in ffff tttt ffft fftf ftff tfff fftt fttf ttff tfft; do
        echo ${PAIR_ALIGN} -c ${c} -s ${file}.fa -o ${file}.c${c}.out
        ${PAIR_ALIGN} -c ${c} -s ${file}.fa -o ${file}.c${c}.out > ${file}.c${c}.stdout
    done
done

# ============================================================
# Run on DNA (Adenoviruses).
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 1 2 3; do
    echo ${PAIR_ALIGN} -a dna -s adeno${i}.fa -o adeno${i}.out
    ${PAIR_ALIGN} -a dna -s adeno${i}.fa -o adeno${i}.out > adeno${i}.stdout
done

# ============================================================
# Run on RNA.
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 1 2 3; do
    echo ${PAIR_ALIGN} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna.out
    ${PAIR_ALIGN} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna.out > adeno${i}-rna.stdout
done
