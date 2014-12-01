#!/bin/sh
#
# Output generation for RazerS.

# We use the current trunk version of 2011-10-18 (r10612) for building the
# reference.
RAZERS="../../../../../seqan-knime-build/release/bin/razers --low-memory"

# ============================================================
# Run Adeno Single-End Example
# ============================================================

for rl in 36 100; do
    # Run with defaults for everything.
    ${RAZERS} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1.razers > se-adeno-reads${rl}_1.stdout

    # Allow indels.
    ${RAZERS} -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id.razers > se-adeno-reads${rl}_1-id.stdout

    # Compute forward/reverse maches only.
    ${RAZERS} -f -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id-f.razers > se-adeno-reads${rl}_1-id-f.stdout
    ${RAZERS} -r -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id-r.razers > se-adeno-reads${rl}_1-id-r.stdout

    # Compute with different identity rates.
    for i in 90 91 92 93 94 95 96 97 98 99 100; do
        ${RAZERS} -i ${i} -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id-i${i}.razers > se-adeno-reads${rl}_1-id-i${i}.stdout
    done

    # Run with different output formats.
    for pair in 0.razers 1.fa 2.eland 3.gff; do
        of=$(echo $pair | sed 's/\..*$//g')
        suffix=$(echo $pair | sed 's/^.\.//g')
        ${RAZERS} -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id-of${of}.${suffix} > se-adeno-reads${rl}_1-id-of${of}.stdout
    done

    # Run with different match ordering.
    for so in 0 1; do
        ${RAZERS} -so ${so} -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id-so${so}.razers > se-adeno-reads${rl}_1-id-so${so}.stdout
    done
done

# ============================================================
# Run Adeno Paired-End Example
# ============================================================

for rl in 36 100; do
    # Run with defaults for everything.
    ${RAZERS} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2.razers > pe-adeno-reads${rl}_2.stdout

    # Allow indels.
    ${RAZERS} -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id.razers > pe-adeno-reads${rl}_2-id.stdout

    # Compute forward/backward maches only.
    ${RAZERS} -f -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id-f.razers > pe-adeno-reads${rl}_2-id-f.stdout
    ${RAZERS} -r -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id-r.razers > pe-adeno-reads${rl}_2-id-r.stdout

    # Compute with different identity rates.
    for i in 90 91 92 93 94 95 96 97 98 99 100; do
        ${RAZERS} -i ${i} -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id-i${i}.razers > pe-adeno-reads${rl}_2-id-i${i}.stdout
    done

    # Run with different output formats.
    for pair in 0.razers 1.fa 2.eland 3.gff; do
        of=$(echo $pair | sed 's/\..*$//g')
        suffix=$(echo $pair | sed 's/^.\.//g')
        ${RAZERS} -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id-of${of}.${suffix} > pe-adeno-reads${rl}_2-id-of${of}.stdout
    done

    # Run with different match ordering.
    for so in 0 1; do
        ${RAZERS} -so ${so} -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id-so${so}.razers > pe-adeno-reads${rl}_2-id-so${so}.stdout
    done
done
