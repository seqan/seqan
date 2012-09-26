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
# below with different parameters on all files and generate output
# files for them.
#
# They can be used to compare against the output of new seqan_tcoffee programs.
#
# TODO(holtgrew): Run with different match file variants.

# Output was generated on 2011-10-20 (r10627).
# Some test files were renamed manually on 2012-07-27 but not rebuild (r12420)
# gernerate_outputs.sh reflects these name changes (r12421) 

TCOFFEE=../../../../build/Debug/core/apps/seqan_tcoffee/seqan_tcoffee

# ============================================================
# Run on Proteins (Balibase).
# ============================================================

# Run with defaults for all non-mandatory options.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -s ${file}.fa -o ${file}.out
    ${TCOFFEE} -s ${file}.fa -o ${file}.out
done

# Run with explicit alphabet.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -a protein -s ${file}.fa -o ${file}.protein.out
    ${TCOFFEE} -a protein -s ${file}.fa -o ${file}.protein.out
done

# Run with different segment match generation options.  We run with with single
# values and combinations of neighbours
for file in 1aab 1ad2 2trx; do
    for m in global local overlap lcs; do
        echo ${TCOFFEE} -m ${m} -s ${file}.fa -o ${file}.m${m}.out
        ${TCOFFEE} -m ${m} -s ${file}.fa -o ${file}.m${m}.out
    done
    m1=global
    m2=local
    echo ${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.out
	${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.out
	
    m1=local
    m2=overlap	
    echo ${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.out
	${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.out
	
    m1=overlap
    m2=lcs
    echo ${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.out
	${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.out
	
    m1=global
    m2=lcs
    echo ${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.out
	${TCOFFEE} -m ${m1} -m ${m2} -s ${file}.fa -o ${file}.m${m1}.m${m2}.out
done

# Run with different match files variations.

# Run with different scoring options.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -g -20 -s ${file}.fa -o ${file}.g-20.out
    ${TCOFFEE} -g -20 -s ${file}.fa -o ${file}.g-20.out
    echo ${TCOFFEE} -e -5 -s ${file}.fa -o ${file}.e-5.out
    ${TCOFFEE} -e -5 -s ${file}.fa -o ${file}.e-5.out
    echo ${TCOFFEE} -ms 10 -s ${file}.fa -o ${file}.ms10.out
    ${TCOFFEE} -ms 10 -s ${file}.fa -o ${file}.ms10.out
    echo ${TCOFFEE} -mm -8 -s ${file}.fa -o ${file}.mm-8.out
    ${TCOFFEE} -mm -8 -s ${file}.fa -o ${file}.mm-8.out
done

# Run with matrix file.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200.out
    ${TCOFFEE} -ma VTML200I -s ${file}.fa -o ${file}.maVTML200.out
done

# Run with manual guide tree.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -u ${file}.newick -s ${file}.fa -o ${file}.u.out
    ${TCOFFEE} -u ${file}.newick -s ${file}.fa -o ${file}.u.out
done

# Run with different guide tree building options.
for file in 1aab 1ad2 2trx; do
    for b in nj min max avg wavg; do
        echo ${TCOFFEE} -b ${b} -s ${file}.fa -o ${file}.b${b}.out
        ${TCOFFEE} -b ${b} -s ${file}.fa -o ${file}.b${b}.out
    done
done

# Run alignment evaluation.
for file in 1aab 1ad2 2trx; do
    echo ${TCOFFEE} -i ${file}.out "(redirecting to ${file}.i.out)"
    ${TCOFFEE} -i ${file}.out > ${file}.i.out
done

# ============================================================
# Run on DNA (Adenoviruses).
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 2 3 4; do
    echo ${TCOFFEE} -a dna -s adeno${i}.fa -o adeno${i}.out
    ${TCOFFEE} -a dna -s adeno${i}.fa -o adeno${i}.out
done

# ============================================================
# Run on RNA.
# ============================================================

# Run with defaults for all non-mandatory options.
for i in 2 3 4; do
    echo ${TCOFFEE} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna.out
    ${TCOFFEE} -a rna -s adeno${i}-rna.fa -o adeno${i}-rna.out
done
