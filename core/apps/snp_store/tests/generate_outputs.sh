#!/bin/sh
#
# Output generation script for snp_store

SNP_STORE=../../../../../seqan-dev/build/clang/bin/snp_store

# ============================================================
# First Section
# ============================================================

genome=human-chr22-inf2.fa
readsGff=human-reads2.gff
readsSam=human-reads2.sam

#echo "${SNP_STORE} $genome $readsGff -o snps_default.vcf -id indels_default.gff > snp_store_default.stdout"
${SNP_STORE} $genome $readsGff -o snps_default.vcf -id indels_default.gff > snp_store_default.stdout

# command options often in use: sam input and do realignment
#echo "${SNP_STORE} $genome $readsSam -re -o snps_realign.vcf -id indels_realign.gff > snp_store_realign.stdout"
${SNP_STORE} $genome $readsSam -re -o snps_realign.vcf -id indels_realign.gff > snp_store_realign.stdout

# orientation aware and pile up correction, threshold method, hide qualites, indel thershold 1
#echo "${SNP_STORE} $genome $readsSam -it 1 -re -oa -mp 1 -m maq -hq -o snps_realign_m1mp1oa.vcf -id indels_realign_m1mp1oa.gff > snp_store_realign_m1mp1oa.stdout"
${SNP_STORE} $genome $readsSam -it 1 -re -oa -mp 1 -m maq -hq -o snps_realign_m1mp1oa.vcf -id indels_realign_m1mp1oa.gff > snp_store_realign_m1mp1oa.stdout

# orientation aware and pile up correction, maq method, hide qualites, indel thershold 2
#echo "${SNP_STORE} $genome $readsSam  -it 2 -re -oa -mp 1 -hq -o snps_realign_m0mp1oa.vcf -id indels_realign_m0mp1oa.gff > snp_store_realign_m0mp1oa.stdout"
${SNP_STORE} $genome $readsGff -it 2 -re -oa -mp 1 -hq -o snps_realign_m0mp1oa.vcf -id indels_realign_m0mp1oa.gff > snp_store_realign_m0mp1oa.stdout

# orientation aware, maq method, hide qualites, indel thershold 1, indel percentage threshold 0.1
#echo "${SNP_STORE} $genome $readsSam -it 1 -ipt 0.1 -osc -re -oa -hq -o snps_realign_m0mp1oa_it1ipt01.vcf -id indels_realign_m0mp1oa_it1ipt01.gff > snp_store_realign_m0mp1oa_it1ipt01.stdout"
${SNP_STORE} $genome $readsSam -it 1 -ipt 0.1 -osc -re -oa -hq -o snps_realign_m0mp1oa_it1ipt01.vcf -id indels_realign_m0mp1oa_it1ipt01.gff > snp_store_realign_m0mp1oa_it1ipt01.stdout


