#!/bin/sh
#
# Output generation for variant_comp.

# We use the current trunk version (r13805) for building the reference.
GUSTAF="../../../../../build/debug/bin/gustaf"
STELLAR="../../../../../build/debug/bin/stellar"

# ============================================================
# Creating Stellar output file
# ============================================================

    ${STELLAR} adeno.fa adeno_modified_reads.fa -l 30 -o stellar.gff > stellar.stdout 2> stellar.stderr

# ============================================================
# Sanity check with default values and empty output file
# ============================================================

out="st2_l100"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -bpo ${out}.gff > ${out}.stdout 2> ${out}.stderr

# ============================================================
# -st 1 -l 30
# ============================================================

out="st1_l30"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -st 1 -l 30 -bpo ${out}.gff > ${out}.stdout 2> ${out}.stderr

# ============================================================
# -st 1 -m stellar.gff
# ============================================================

out="st1_l30_m"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -st 1 -m stellar.gff -bpo ${out}.gff > ${out}.stdout 2> ${out}.stderr

# ============================================================
# -st 1 -l 30 -ith 5
# ============================================================

out="st1_l30_ith5"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -st 1 -l 30 -ith 5 -bpo ${out}.gff > ${out}.stdout 2> ${out}.stderr

# ============================================================
# -st 1 -l 30 -gth 3
# ============================================================

out="st1_l30_gth3"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -st 1 -l 30 -gth 3 -bpo ${out}.gff > ${out}.stdout 2> ${out}.stderr
