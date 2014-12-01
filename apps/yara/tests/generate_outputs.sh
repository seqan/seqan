#!/bin/bash
#
# Yara output generation.

PATH=~/Code/seqan-builds/Release-Gcc/bin
INDEXER=$PATH/yara_indexer
MAPPER=$PATH/yara_mapper

# ============================================================
# Run Indexer
# ============================================================

# Run with different organisms.
for organism in adeno; do
    ${INDEXER} input/$organism-genome.fa -o gold/$organism-genome
done

# ============================================================
# Run Single-End Mapper
# ============================================================

# Run with different threads.
MAPPER_ARGS=("--threads 1") # "--threads 8")
MAPPER_SUFFIX=("t1") # "t8")

# Run with different organismnisms.
for organism in adeno; do
    # Run with different arguments.
    for ((a=0; a<${#MAPPER_ARGS[@]}; a++)); do
        ${MAPPER} gold/$organism-genome input/$organism-reads_1.fa -o gold/$organism-reads_1.${MAPPER_SUFFIX[$a]}.sam ${MAPPER_ARGS[$a]}
    done
done
