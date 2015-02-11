#!/bin/sh
#
# Currently, we run RazerS only on

# We use revision 13816 for building the gold standard.
RAZERS=${RAZERS:-../../../../../seqan-trunk-build/release/bin/razers3}

for tc in 0 1; do  # number of threads, 0 and >0 have different code paths

# ============================================================
# Run Adeno Single-End Example
# ============================================================

for rl in 36 100; do
    # Run with defaults for everything.
    ${RAZERS} -tc ${tc} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-tc${tc}.razers > se-adeno-reads${rl}_1-tc${tc}.stdout

    # No gaps.
    ${RAZERS} -tc ${tc} --no-gaps adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-ng-tc${tc}.razers > se-adeno-reads${rl}_1-ng-tc${tc}.stdout

    # Compute forward/reverse maches only.
    ${RAZERS} -tc ${tc} -f adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-f-tc${tc}.razers > se-adeno-reads${rl}_1-f-tc${tc}.stdout
    ${RAZERS} -tc ${tc} -r adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-r-tc${tc}.razers > se-adeno-reads${rl}_1-r-tc${tc}.stdout

    # Compute with different identity rates.
    for i in 90 91 92 93 94 95 96 97 98 99 100; do
        ${RAZERS} -tc ${tc} -i ${i} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-i${i}-tc${tc}.razers > se-adeno-reads${rl}_1-i${i}-tc${tc}.stdout
    done

    # Run with different output formats.
    for pair in 0.razers 1.fa 2.eland 3.gff 4.sam 5.afg; do
        of=$(echo $pair | sed 's/\..*$//g')
        suffix=$(echo $pair | sed 's/^.\.//g')
        ${RAZERS} -tc ${tc} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-of${of}-tc${tc}.${suffix} > se-adeno-reads${rl}_1-of${of}-tc${tc}.stdout
    done

    # Run with different match ordering.
    for so in 0 1; do
        ${RAZERS} -tc ${tc} -so ${so} adeno-genome.fa adeno-reads${rl}_1.fa -o se-adeno-reads${rl}_1-so${so}-tc${tc}.razers > se-adeno-reads${rl}_1-so${so}-tc${tc}.stdout
    done
done

# ============================================================
# Run Adeno Paired-End Example
# ============================================================

for rl in 36 100; do
    # Run with defaults for everything.
    ${RAZERS} -tc ${tc} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-tc${tc}.razers > pe-adeno-reads${rl}_2-tc${tc}.stdout

    # Allow indels.
    ${RAZERS} -tc ${tc} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-tc${tc}.razers > pe-adeno-reads${rl}_2-tc${tc}.stdout

    # Compute forward/backward maches only.
    ${RAZERS} -tc ${tc} -f adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-f-tc${tc}.razers > pe-adeno-reads${rl}_2-f-tc${tc}.stdout
    ${RAZERS} -tc ${tc} -r adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-r-tc${tc}.razers > pe-adeno-reads${rl}_2-r-tc${tc}.stdout

    # Compute with different identity rates.
    for i in 90 91 92 93 94 95 96 97 98 99 100; do
        ${RAZERS} -tc ${tc} -i ${i} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-i${i}-tc${tc}.razers > pe-adeno-reads${rl}_2-i${i}-tc${tc}.stdout
    done

    # Run with different output formats.
    for pair in 0.razers 1.fa 2.eland 3.gff 4.sam 5.afg; do
        of=$(echo $pair | sed 's/\..*$//g')
        suffix=$(echo $pair | sed 's/^.\.//g')
        ${RAZERS} -tc ${tc} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-of${of}-tc${tc}.${suffix} > pe-adeno-reads${rl}_2-of${of}-tc${tc}.stdout
    done

    # Run with different match ordering.
    for so in 0 1; do
        ${RAZERS} -tc ${tc} -so ${so} adeno-genome.fa adeno-reads${rl}_1.fa adeno-reads${rl}_2.fa -o pe-adeno-reads${rl}_2-so${so}-tc${tc}.razers > pe-adeno-reads${rl}_2-so${so}-tc${tc}.stdout
    done
done

done  # for nt
