#!/bin/sh
#
# Expected output generation for rabema.

# We use the current trunk version of 2011-10-13 (r10463) for building the
# reference.
RAZERS=../../../../build/Release/apps/razers2/razers2
RABEMA_PREPARE=../../../../build/Release/apps/rabema/rabema_prepare_sam
RABEMA_BUILD=../../../../build/Release/apps/rabema/rabema_build_gold_standard
RABEMA_EVALUATE=../../../../build/Release/apps/rabema/rabema_evaluate

# ============================================================
# Map reads for gold standard and with too low error rate.
# ============================================================

${RAZERS} -m 10000 -vv -of 4 -ds -i 92     -o gold-adeno-hamming-08.sam adeno-genome.fa reads.fasta
${RAZERS} -m 10000 -vv -of 4 -ds -i 92 -id -o gold-adeno-edit-08.sam    adeno-genome.fa reads.fasta

${RAZERS} -vv -of 4 -ds -i 92     -o razers2-adeno-hamming-08.sam adeno-genome.fa reads.fasta
${RAZERS} -vv -of 4 -ds -i 92 -id -o razers2-adeno-edit-08.sam    adeno-genome.fa reads.fasta

${RAZERS} -vv -of 4 -ds -i 96     -o razers2-adeno-hamming-04.sam adeno-genome.fa reads.fasta
${RAZERS} -vv -of 4 -ds -i 96 -id -o razers2-adeno-edit-04.sam    adeno-genome.fa reads.fasta

# ============================================================
# Prepare SAM.
# ============================================================

${RABEMA_PREPARE} -i gold-adeno-hamming-08.sam -o gold-adeno-hamming-08.by_qname.sam
samtools view -Sb gold-adeno-hamming-08.by_qname.sam > gold-adeno-hamming-08.by_qname.bam
samtools sort gold-adeno-hamming-08.by_qname.bam gold-adeno-hamming-08.by_coordinate
samtools view gold-adeno-hamming-08.by_coordinate.bam > gold-adeno-hamming-08.by_coordinate.sam
${RABEMA_PREPARE} -i gold-adeno-edit-08.sam > gold-adeno-edit-08.by_qname.sam
samtools view -Sb gold-adeno-edit-08.by_qname.sam -o gold-adeno-edit-08.by_qname.bam
samtools sort gold-adeno-edit-08.by_qname.bam gold-adeno-edit-08.by_coordinate
samtools view gold-adeno-edit-08.by_coordinate.bam > gold-adeno-edit-08.by_coordinate.sam

# ============================================================
# Build Gold Standard
# ============================================================

${RABEMA_BUILD} --distance-metric hamming -e 8 -o gold-adeno-hamming-08.gsi --reference adeno-genome.fa --in-bam gold-adeno-hamming-08.by_coordinate.sam > gold-adeno-hamming-08.stdout
${RABEMA_BUILD} --distance-metric edit    -e 8 -o gold-adeno-edit-08.gsi    --reference adeno-genome.fa --in-bam gold-adeno-edit-08.by_coordinate.sam > gold-adeno-edit-08.stdout

# ============================================================
# Compare Against Gold Standard
# ============================================================

${RABEMA_EVALUATE} --distance-metric hamming -e 8 --reference adeno-genome.fa --in-bam razers2-adeno-hamming-08.sam --in-gsi gold-adeno-hamming-08.gsi > razers2-adeno-hamming-08.stdout
${RABEMA_EVALUATE} --distance-metric hamming -e 8 --reference adeno-genome.fa --in-bam razers2-adeno-hamming-04.sam --in-gsi gold-adeno-hamming-08.gsi > razers2-adeno-hamming-04.stdout
${RABEMA_EVALUATE} --distance-metric edit    -e 8 --reference adeno-genome.fa --in-bam razers2-adeno-edit-08.sam    --in-gsi gold-adeno-edit-08.gsi > razers2-adeno-edit-08.stdout
${RABEMA_EVALUATE} --distance-metric edit    -e 8 --reference adeno-genome.fa --in-bam razers2-adeno-edit-04.sam    --in-gsi gold-adeno-edit-08.gsi > razers2-adeno-edit-04.stdout
