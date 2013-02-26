#!/bin/sh
#
# Currently, we run RazerS only on

# We use revision 13383 for building the gold standard.
RAZERS=../../../../build/make/bin/razers3

# ============================================================
# Run Adeno Single-End Example
# ============================================================

for rl in 36 100; do
    # Run with defaults for everything.
    ${RAZERS} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1.razers > se-adeno-reads${rl}_1.stdout

    # No gaps.
    ${RAZERS} --no-gaps adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-ng.razers > se-adeno-reads${rl}_1-ng.stdout

    # Compute forward/reverse maches only.
    ${RAZERS} -f adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-f.razers > se-adeno-reads${rl}_1-f.stdout
    ${RAZERS} -r adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-r.razers > se-adeno-reads${rl}_1-r.stdout

    # Compute with different identity rates.
    for i in 90 91 92 93 94 95 96 97 98 99 100; do
        ${RAZERS} -i ${i} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-i${i}.razers > se-adeno-reads${rl}_1-i${i}.stdout
    done

    # Run with different output formats.
    for pair in 0.razers 1.fa 2.eland 3.gff 4.sam 5.afg; do
        of=$(echo $pair | sed 's/\..*$//g')
        suffix=$(echo $pair | sed 's/^.\.//g')
        ${RAZERS} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-of${of}.${suffix} > se-adeno-reads${rl}_1-of${of}.stdout
    done

    # Run with different match ordering.
    for so in 0 1; do
        ${RAZERS} -so ${so} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-so${so}.razers > se-adeno-reads${rl}_1-so${so}.stdout
    done
done

# ============================================================
# Run Adeno Paired-End Example
# ============================================================

for rl in 36 100; do
    # Run with defaults for everything.
    ${RAZERS} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2.razers > pe-adeno-reads${rl}_2.stdout

    # Allow indels.
    ${RAZERS} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2.razers > pe-adeno-reads${rl}_2.stdout

    # Compute forward/backward maches only.
    ${RAZERS} -f adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-f.razers > pe-adeno-reads${rl}_2-f.stdout
    ${RAZERS} -r adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-r.razers > pe-adeno-reads${rl}_2-r.stdout

    # Compute with different identity rates.
    for i in 90 91 92 93 94 95 96 97 98 99 100; do
        ${RAZERS} -i ${i} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-i${i}.razers > pe-adeno-reads${rl}_2-i${i}.stdout
    done

    # Run with different output formats.
    for pair in 0.razers 1.fa 2.eland 3.gff 4.sam 5.afg; do
        of=$(echo $pair | sed 's/\..*$//g')
        suffix=$(echo $pair | sed 's/^.\.//g')
        ${RAZERS} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-of${of}.${suffix} > pe-adeno-reads${rl}_2-of${of}.stdout
    done

    # Run with different match ordering.
    for so in 0 1; do
        ${RAZERS} -so ${so} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-so${so}.razers > pe-adeno-reads${rl}_2-so${so}.stdout
    done
done
