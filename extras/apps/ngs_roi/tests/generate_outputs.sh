#!/bin/bash
#
# Output generation for roi_intersect

# Configure Bash.
set -x  # print executed commands
set -e  # exit immediately on errors

# Make execution errors visible.
function check_result {
  [ "$?" == "0" ] || echo "HALTED: program execution failed" 1>&2
}
trap check_result EXIT

# Shortcut for binary bas path.
BINDIR=../../../../../seqan-stream-build/debug/bin

PROJECT=${BINDIR}/roi_feature_projection
BAM2ROI=${BINDIR}/bam2roi

# ----------------------------------------------------------------------------
# Projection to BED or GFF/GTF in BED style.
# ----------------------------------------------------------------------------

for mode in intersection projection union difference; do
    for format in bed gff gtf; do
        ${PROJECT} --mode ${mode} -ir small.roi -if small.${format} -or out_small_${format}_m${mode}.roi >out_small_${format}_m${mode}.stdout 2>out_small_${format}_m${mode}.stderr
        ${PROJECT} --strand-specific --mode ${mode} -ir small.roi -if small.${format} -or out_small_${format}_m${mode}_ss.roi >out_small_${format}_m${mode}_ss.stdout 2>out_small_${format}_m${mode}_ss.stderr
    done
done

# TODO(holtgrew): App test with projection to transcripts from GFF/GTF

# ----------------------------------------------------------------------------
# Projection to BED or GFF/GTF in BED style.
# ----------------------------------------------------------------------------

${BAM2ROI} -if micro_rna_sorted_2l.bam -of out_mrna_2l.roi >out_mrna_2l.roi.stdout 2>out_mrna_2l.roi.stderr
${BAM2ROI} --strand-specific -if micro_rna_sorted_2l.bam -of out_mrna_2l_ss.roi >out_mrna_2l_ss.roi.stdout 2>out_mrna_2l_ss.roi.stderr
