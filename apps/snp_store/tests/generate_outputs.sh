#!/bin/sh
#
# Output generation script for snp_store

#SNP_STORE=../../.. ../../../../seqan-build/release/bin/snp_store

# ============================================================
# First Section
# ============================================================

genome=human-chr22-inf2.fa
readsGff=human-reads2.gff
readsSam=human-reads2.sam

#echo "snp_store $genome $readsGff -o snps_default.vcf -id indels_default.vcf > snp_store_default.stdout"
snp_store $genome $readsGff -o snps_default.vcf -id indels_default.vcf > snp_store_default.stdout

# command options often in use: sam input and do realignment, show qualities
#echo "snp_store $genome $readsSam -re -sq -o snps_realign.vcf -id indels_realign.vcf > snp_store_realign.stdout"
snp_store $genome $readsSam -re -sq -o snps_realign.vcf -id indels_realign.vcf > snp_store_realign.stdout

# orientation aware and pile up correction, threshold method, indel thershold 1
#echo "snp_store $genome $readsSam -it 1 -re -oa -mp 1 -m maq -o snps_realign_m1mp1oa.vcf -id indels_realign_m1mp1oa.vcf > snp_store_realign_m1mp1oa.stdout"
snp_store $genome $readsSam -it 1 -re -oa -mp 1 -m maq -o snps_realign_m1mp1oa.vcf -id indels_realign_m1mp1oa.vcf > snp_store_realign_m1mp1oa.stdout

# orientation aware and pile up correction, maq method, indel thershold 2
#echo "snp_store $genome $readsSam  -it 2 -re -oa -mp 1 -o snps_realign_m0mp1oa.vcf -id indels_realign_m0mp1oa.vcf > snp_store_realign_m0mp1oa.stdout"
snp_store $genome $readsGff -it 2 -re -oa -mp 1 -o snps_realign_m0mp1oa.vcf -id indels_realign_m0mp1oa.vcf > snp_store_realign_m0mp1oa.stdout

# orientation aware, maq method, indel thershold 1, indel percentage threshold 0.1
#echo "snp_store $genome $readsSam -it 1 -ipt 0.1 -osc -re -oa -o snps_realign_m0mp1oa_it1ipt01.vcf -id indels_realign_m0mp1oa_it1ipt01.vcf > snp_store_realign_m0mp1oa_it1ipt01.stdout"
snp_store $genome $readsSam -it 1 -ipt 0.1 -osc -re -oa -o snps_realign_m0mp1oa_it1ipt01.vcf -id indels_realign_m0mp1oa_it1ipt01.vcf > snp_store_realign_m0mp1oa_it1ipt01.stdout

# the following line replaces the path to the reference genome with the path, the snp-store would use when calle via run_test.py to avoid failing of the tests du to path differences.
sed -i '4s@.*@##reference=../../../apps/snp_store/tests/human-chr22-inf2.fa@' *.vcf

