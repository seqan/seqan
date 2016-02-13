#!/bin/sh
#
# Currently, we run RazerS only on

# We use revision 13383 for building the gold standard.
RAZERS=${RAZERS:-../../../../../seqan-knime-build/release/bin/razers2}

# ============================================================
# Run Adeno Single-End Example
# ============================================================

for rl in 36 100; do
    # Run with defaults for everything.
    ${RAZERS} --low-memory adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1.razers > se-adeno-reads${rl}_1.stdout

    # Allow indels.
    ${RAZERS} --low-memory -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id.razers > se-adeno-reads${rl}_1-id.stdout

    # Compute forward/reverse maches only.
    ${RAZERS} --low-memory -f -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id-f.razers > se-adeno-reads${rl}_1-id-f.stdout
    ${RAZERS} --low-memory -r -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id-r.razers > se-adeno-reads${rl}_1-id-r.stdout

    # Compute with different identity rates.
    for i in 90 91 92 93 94 95 96 97 98 99 100; do
        ${RAZERS} --low-memory -i ${i} -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id-i${i}.razers > se-adeno-reads${rl}_1-id-i${i}.stdout
    done

    # Run with different output formats.
    for suf in razers fa eland gff sam afg; do
        ${RAZERS} --low-memory -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id.${suf} > se-adeno-reads${rl}_1-id.${suf}.stdout
    done

    # Run with different match ordering.
    for so in 0 1; do
        ${RAZERS} --low-memory -so ${so} -id adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-id-so${so}.razers > se-adeno-reads${rl}_1-id-so${so}.stdout
    done
done

# ============================================================
# Run Adeno Paired-End Example
# ============================================================

for rl in 36 100; do
    # Run with defaults for everything.
    ${RAZERS} --low-memory adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2.razers > pe-adeno-reads${rl}_2.stdout

    # Allow indels.
    ${RAZERS} --low-memory -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id.razers > pe-adeno-reads${rl}_2-id.stdout

    # Compute forward/backward maches only.
    ${RAZERS} --low-memory -f -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id-f.razers > pe-adeno-reads${rl}_2-id-f.stdout
    ${RAZERS} --low-memory -r -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id-r.razers > pe-adeno-reads${rl}_2-id-r.stdout

    # Compute with different identity rates.
    for i in 90 91 92 93 94 95 96 97 98 99 100; do
        ${RAZERS} --low-memory -i ${i} -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id-i${i}.razers > pe-adeno-reads${rl}_2-id-i${i}.stdout
    done

    # Run with different output formats.
    for suf in razers fa eland gff sam afg; do
        ${RAZERS} --low-memory -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id.${suf} > pe-adeno-reads${rl}_2-id.${suf}.stdout     
    done

    # Run with different match ordering.
    for so in 0 1; do
        ${RAZERS} --low-memory -so ${so} -id adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-id-so${so}.razers > pe-adeno-reads${rl}_2-id-so${so}.stdout
    done
done
