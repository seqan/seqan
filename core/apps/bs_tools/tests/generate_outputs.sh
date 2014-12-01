#!/bin/sh
#
# Output generation script for bs_tools
# prepare data... createTestData.sh 


##############################################################
# Bisar
##############################################################
BISAR=../../../../../build/Release/bin/bisar

# ============================================================
# se 
# ============================================================

# default
${BISAR} -e3 4 -e4 5 -o reads_se_N6000_0.CT_GA.verified.sam reads_se_N6000.CT_GA.sam hg18_chr21_3000.fa reads_se_N6000.fastq > bisar_out_se_0.txt

# same as above currently
${BISAR} -gas -4.5 -ges -2.0 -der 0.001 -bsc 0.99 -gmr 0.5 -i 0.8 -rn 0.001 -pms 0.9 -e3 4 -e4 5 -o reads_se_N6000_1.CT_GA.verified.sam reads_se_N6000.CT_GA.sam hg18_chr21_3000.fa reads_se_N6000.fastq > bisar_out_se_1.txt

# non-uniform model 
${BISAR} -nse -nsi -nsd -gas -4.5 -ges -2.0 -der 0.001 -bsc 0.99 -gmr 0.5 -i 0.8 -rn 0.001 -pms 0.9 -e3 4 -e4 5 -o reads_se_N6000_2.CT_GA.verified.sam reads_se_N6000.CT_GA.sam hg18_chr21_3000.fa reads_se_N6000.fastq  > bisar_out_se_2.txt

# global methylation rate ...
${BISAR} -nse -nsi -nsd -gas -4.5 -ges -2.0 -der 0.001 -bsc 0.99 -gmr 0.2 -i 0.8 -rn 0.001 -pms 0.9 -e3 4 -e4 5 -o reads_se_N6000_3.CT_GA.verified.sam reads_se_N6000.CT_GA.sam hg18_chr21_3000.fa reads_se_N6000.fastq  > bisar_out_se_3.txt
${BISAR} -nse -nsi -nsd -gas -4.5 -ges -2.0 -der 0.001 -bsc 0.99 -gmr 0.8 -i 0.8 -rn 0.001 -pms 0.9 -e3 4 -e4 5 -o reads_se_N6000_4.CT_GA.verified.sam reads_se_N6000.CT_GA.sam hg18_chr21_3000.fa reads_se_N6000.fastq  > bisar_out_se_4.txt

# ============================================================
# pe
# ============================================================

# default
${BISAR} -e3 4 -e4 5 -o reads_pe_N6000_0.CT_GA.verified.sam reads_pe_N6000.CT_GA.sam hg18_chr21_3000.fa reads_pe_N6000.L.fastq reads_pe_N6000.R.fastq > bisar_out_pe_0.txt



# sort by coordinate before calling
# java -Xmx2g -jar /home/takifugu/sabrina7/picard-tools-1.84/SortSam.jar INPUT=reads_se_N6000_2.CT_GA.verified.sam OUTPUT=reads_se_N6000_2.CT_GA.verified.pos_so.sam SORT_ORDER=coordinate 
# java -Xmx2g -jar /home/takifugu/sabrina7/picard-tools-1.84/SortSam.jar INPUT=reads_pe_N6000_0.CT_GA.verified.sam OUTPUT=reads_pe_N6000_0.CT_GA.verified.pos_so.sam SORT_ORDER=coordinate 
##############################################################
# Casbar
##############################################################
CASBAR=../../../../../build/Release/bin/casbar

# ============================================================
# se
# ============================================================

${CASBAR} -nec -mc 6 -msc 5 -mpc 0.5 -hes 0.005 -o snps_se_0.vcf -b meths_se_0.bed  hg18_chr21_3000.fa reads_se_N6000_2.CT_GA.verified.pos_so.sam > casbar_out_se_0.stdout

${CASBAR} -nec -mc 2 -msc 3 -mpc 0.5 -hes 0.005 -o snps_se_1.vcf -b meths_se_1.bed  hg18_chr21_3000.fa reads_se_N6000_2.CT_GA.verified.pos_so.sam > casbar_out_se_1.stdout

# ============================================================
# pe
# ============================================================

${CASBAR} -nec -mc 6 -msc 5 -mpc 0.5 -hes 0.005 -o snps_pe_0.vcf -b meths_pe_0.bed  hg18_chr21_3000.fa reads_pe_N6000_0.CT_GA.verified.pos_so.sam > casbar_out_pe_0.stdout




