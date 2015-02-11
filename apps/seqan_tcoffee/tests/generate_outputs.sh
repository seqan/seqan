#!/bin/sh
#
# The files adeno{2,3,4}.fa contain fasta files with up to four sequences,
# each sequence has the first 280 base pairs of an Adenovirus.
#
# The files adeno{2,3,4}-rna.fa contain fasta files with up to four sequences,
# each sequence has the first 280 base pairs of an Adenovirus, however all
# T characters are replaced by an U.
#
# The files 1ad2.fa, 1aab.fa and 2trx.fa are sequences from the balibase
# benchmark set in FASTA format.
#
# The script will run the seqan_tcoffee program as configured
# below with different parameters on all files and generate fastaput
# files for them.
#
# They can be used to compare against the fastaput of new seqan_tcoffee programs.
#
# TODO(holtgrew): Run with different match file variants.

# fastaput was generated on 2011-10-20 (r10627).
# Some test files were renamed manually on 2012-07-27 but not rebuild (r12420)
# gernerate_fastaputs.sh reflects these name changes (r12421) 

#TCOFFEE=../../../../../build/Debug/apps/seqan_tcoffee/seqan_tcoffee
TCOFFEE=../../../../../seqan-trunk-build/debug/bin/seqan_tcoffee
# ============================================================
# Run on Proteins (Balibase).
# ============================================================

# Run with defaults for all non-mandatory options.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -s ${file}.fa -o ${file}.fasta
    ${TCOFFEE} -s ${file}.fa -o ${file}.fasta
done

# Run with explicit alphabet.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -a protein -s ${file}.fa -o ${file}.protein.fasta
    ${TCOFFEE} -a protein -s ${file}.fa -o ${file}.protein.fasta
done

# Run with different segment match generation options.  We run with with single
# values and combinations of neighbours
for file in 1aab 1ad2 2trx; do
    for m in global local overlap lcs; do
        echo ${TCOFFEE} -m ${m} -s ${file}.fa -o ${file}.m${m}.fasta
        ${TCOFFEE} -m ${m} -s ${file}.fa -o ${file}.m${m}.fasta
    done
    m1=global
    m2=local
    echo ${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.fasta
	${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.fasta
	
    m1=local
    m2=overlap	
    echo ${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.fasta
	${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.fasta
	
    m1=overlap
    m2=lcs
    echo ${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.fasta
	${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.fasta
	
    m1=global
    m2=lcs
    echo ${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.fasta
	${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.fasta
done

# Run with different match files variations.

# Run with different scoring options.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -g -20 -s ${file}.fa -o ${file}.g-20.fasta
    ${TCOFFEE} -g -20 -s ${file}.fa -o ${file}.g-20.fasta
    echo ${TCOFFEE} -e -5 -s ${file}.fa -o ${file}.e-5.fasta
    ${TCOFFEE} -e -5 -s ${file}.fa -o ${file}.e-5.fasta
    echo ${TCOFFEE} -ms 10 -s ${file}.fa -o ${file}.ms10.fasta
    ${TCOFFEE} -ms 10 -s ${file}.fa -o ${file}.ms10.fasta
    echo ${TCOFFEE} -mm -8 -s ${file}.fa -o ${file}.mm-8.fasta
    ${TCOFFEE} -mm -8 -s ${file}.fa -o ${file}.mm-8.fasta
done

# Run with matrix file.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200.fasta
    ${TCOFFEE} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200.fasta
done

# Run with manual guide tree.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -u ${file}.newick -s ${file}.fa -o ${file}.u.fasta
    ${TCOFFEE} -u ${file}.newick -s ${file}.fa -o ${file}.u.fasta
done

# Run with different guide tree building options.
for file in 1aab 1ad2 2trx; do
    for b in nj min max avg wavg; do
        echo ${TCOFFEE} -b ${b} -s ${file}.fa -o ${file}.b${b}.fasta
        ${TCOFFEE} -b ${b} -s ${file}.fa -o ${file}.b${b}.fasta
    done
done

# Run alignment evaluation.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -i ${file}.fasta "(redirecting to ${file}.i.fasta)"
    ${TCOFFEE} -i ${file}.fasta > ${file}.i.fasta
done

# ============================================================
# Run on DNA (Adenoviruses).
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 2 3 4; do
    echo ${TCOFFEE} -a dna -s adeno${i}.fa -o adeno${i}.fasta
    ${TCOFFEE} -a dna -s adeno${i}.fa -o adeno${i}.fasta
done

# ============================================================
# Run on RNA.
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 2 3 4; do
    echo ${TCOFFEE} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna.fasta
    ${TCOFFEE} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna.fasta
done
