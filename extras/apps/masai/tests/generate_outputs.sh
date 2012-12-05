#!/bin/sh
#
# Masai reference output generation. Version v0.6.1 was used.

PATH=~/Documents/Code/SeqAn-Trunk/build/Release/bin/
#PATH=../../../../../seqan-trunk-build/Release/extras/apps/masai
INDEXER=$PATH/masai_indexer
MAPPER=$PATH/masai_mapper
SINGLE=$PATH/masai_output_se
PAIRED=$PATH/masai_output_pe

# ============================================================
# Run Indexer
# ============================================================

# Run with different indices.
for i in sa esa fm qgram; do
    ${INDEXER} --index ${i} --index-prefix adeno-index-${i}.out adeno-genome.fa &> adeno-index-${i}.stdout
done

# ============================================================
# Run Mapper Single-End
# ============================================================

for rl in 100; do #36 100; do
    # Run with default parameters.
    ${MAPPER} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1.out &> se-adeno-reads${rl}_1.stdout

    # Run with different seed length.
    for sl in 16 50; do
        ${MAPPER} --seed-length ${sl} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-sl${sl}.out &> se-adeno-reads${rl}_1-sl${sl}.stdout
    done

    # Run with different indices.
#    for i in in esa fm qgram; do
#        ${MAPPER} --index ${i} --index-prefix adeno-index-${i}.out adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-i${i}.out &> se-adeno-reads${rl}_1-i${i}.stdout
#    done

    # Run with different mapping modes.
    for mm in all; do
        ${MAPPER} --mapping-mode ${mm} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-mm${mm}.out &> se-adeno-reads${rl}_1-mm${mm}.stdout
    done

    # Run with different absolute number of errors.
    for e in 0 1 2 3 4; do
        ${MAPPER} --errors ${e} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-e${e}.out &> se-adeno-reads${rl}_1-e${e}.stdout
    done

    # Run with different output formats.
    for of in sam sam-no-cigar; do
        ${MAPPER} --output-format ${of} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-of${of}.out &> se-adeno-reads${rl}_1-of${of}.stdout
    done

    # Run without gaps.
    ${MAPPER} --no-gaps adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-nogaps.out &> se-adeno-reads${rl}_1-nogaps.stdout
done

# ============================================================
# Run Mapper Paired-End
# ============================================================

# TODO(esiragusa): Write Paired-End tests.
