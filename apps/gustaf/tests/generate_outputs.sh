#!/bin/sh
#
# Output generation for variant_comp.

# We use the current trunk version (r13805) for building the reference.
GUSTAF="../../../../build/debug/bin/gustaf"
JOINMATES="../../../../build/debug/bin/gustaf_mate_joining"
STELLAR="../../../../build/debug/bin/stellar"

# ============================================================
# Creating Stellar input files for paired-end mode using gustaf_mate_joining
# ============================================================

    ${JOINMATES} adeno_modified_reads_mates1.fa adeno_modified_reads_mates2.fa \
    -o adeno_modified_reads_joinedMates.fa -rc \
    > gustaf_mate_joining.stdout 2> gustaf_mate_joining.stderr

# ============================================================
# Gustaf_mate_joining output files
# ============================================================

    ${JOINMATES} reads_simulated_mates1_gold.fa reads_simulated_mates2_gold.fa \
    -o reads_simulated_joined_rc.fa -rc 1 \
    > gustaf_mate_joining.stdout 2> gustaf_mate_joining.stderr

    ${JOINMATES} reads_simulated_mates1_gold.fa reads_simulated_mates2_gold.fa \
    -o reads_simulated_joined.fa \
    > gustaf_mate_joining.stdout 2> gustaf_mate_joining.stderr

    ${JOINMATES} reads_simulated_joined_gold.fa \
    -o reads_simulated_mates1_rc.fa -o reads_simulated_mates2_rc.fa -rc 1 \
    > gustaf_mate_joining.stdout 2> gustaf_mate_joining.stderr

    ${JOINMATES} reads_simulated_joined_gold.fa \
    -o reads_simulated_mates1.fa -o reads_simulated_mates2.fa \
    > gustaf_mate_joining.stdout 2> gustaf_mate_joining.stderr

# ============================================================
# Creating Stellar output files
# ============================================================

#single-end
    ${STELLAR} adeno.fa adeno_modified_reads.fa -l 30 -o stellar.gff > stellar.stdout 2> stellar.stderr
#paired-end
    ${STELLAR} adeno.fa adeno_modified_reads_joinedMates.fa -l 30 -o stellar_joinedMates_l30.gff > \
    stellar_joinedMates.stdout 2> stellar_joinedMates.stderr

# ============================================================
# Sanity check with default values and empty output file
# ============================================================

out="st2_l100"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -gff ${out}.gff -vcf ${out}.vcf > ${out}.stdout 2> ${out}.stderr

# ============================================================
# -st 1 -l 30
# ============================================================

out="st1_l30"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -st 1 -l 30 -gff ${out}.gff -vcf ${out}.vcf > ${out}.stdout 2> ${out}.stderr

# ============================================================
# -st 1 -m stellar.gff
# ============================================================

out="st1_l30_m"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -st 1 -m stellar.gff -gff ${out}.gff -vcf ${out}.vcf > ${out}.stdout 2> ${out}.stderr

# ============================================================
# -st 1 -l 30 -ith 5
# ============================================================

out="st1_l30_ith5"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -st 1 -l 30 -ith 5 -bth 5 -gff ${out}.gff -vcf ${out}.vcf > ${out}.stdout 2> ${out}.stderr

# ============================================================
# -st 1 -l 30 -gth 3
# ============================================================

out="st1_l30_gth3"
    ${GUSTAF} adeno.fa adeno_modified_reads.fa -st 1 -l 30 -gth 3 -gff ${out}.gff -vcf ${out}.vcf > ${out}.stdout 2> ${out}.stderr

# ============================================================
# paired-end
# -st 1 -m stellar_joinedMates_l30.gff
# ============================================================

out="pairedEnd_st1_l30"
    ${GUSTAF} adeno.fa adeno_modified_reads_mates1.fa adeno_modified_reads_mates2.fa -m stellar_joinedMates_l30.gff -st 1 \
    -ll 1000 -le 30 -rc -gff ${out}.gff -vcf ${out}.vcf > ${out}.stdout 2> ${out}.stderr

# ============================================================
# Sanity check vcf out multiple references
# -st 1 -l 30
# ============================================================

out="reference2_st1_l30"
    ${GUSTAF} adeno.fa read_reference2.fa -st 1 \
    -l 30 -gff ${out}.gff -vcf ${out}.vcf > ${out}.stdout 2> ${out}.stderr

# ============================================================
# paired-end no artificial breakpoint
# -st 1 -m stellar_joinedMates_l30.gff -ll 800 -le 30 -rc
# ============================================================
#
#out="pairedEnd_st1_l30_ll800_gold"
#    ${GUSTAF} adeno.fa adeno_modified_reads_mates1.fa adeno_modified_reads_mates2.fa -m stellar_joinedMates_l30.gff -st 1 \
#    -ll 800 -le 30 -rc -gff ${out}.gff -vcf ${out}.vcf -do -j ${out} > ${out}.stdout 2> ${out}.stderr
#
