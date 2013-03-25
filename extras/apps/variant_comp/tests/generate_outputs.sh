#!/bin/sh
#
# Output generation for variant_comp.

# We use the current trunk version (r13805) for building the reference.
VARIANT_COMP="../../../../../seqan-trunk-build/release/bin/variant_comp"

for f in predicted_full predicted_wrong gold_standard; do
    ${VARIANT_COMP} adeno.fa -ir gold_standard.gff -ip ${f}.gff -pt 10 -st 10 -v -o ${f}_pt10_st10_out.txt -on ${f}_pt10_st10_outFN.txt > ${f}_pt10_st10.stdout 2> ${f}_pt10_st10.stderr

    ${VARIANT_COMP} adeno.fa -ir gold_standard.gff -ip ${f}.gff -pt 1 -st 10 -v -o ${f}_pt1_st10_out.txt -on ${f}_pt1_st10_outFN.txt > ${f}_pt1_st10.stdout 2> ${f}_pt1_st10.stderr
    
    ${VARIANT_COMP} adeno.fa -ir gold_standard.gff -ip ${f}.gff -pt 10 -st 1 -v -o ${f}_pt10_st1_out.txt -on ${f}_pt10_st1_outFN.txt > ${f}_pt10_st1.stdout 2> ${f}_pt10_st1.stderr
    
    ${VARIANT_COMP} adeno.fa -ir gold_standard.gff -ip ${f}.gff -pt 1 -st 1 -v -o ${f}_pt1_st1_out.txt -on ${f}_pt1_st1_outFN.txt > ${f}_pt1_st1.stdout 2> ${f}_pt1_st1.stderr
done
