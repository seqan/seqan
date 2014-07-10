#!/bin/sh
#
# We create variants from a randomly generated genome.

SEQCONS=../../../../../seqan-git-build/debug-clang/bin/seqcons2

# ============================================================
# seqcons2
# ============================================================

# MSA Consensus and NOP for sequences.

echo "${SEQCONS} -m msa_consensus -i seqs1.fa -oa seqs1.msa_consensus.sam -oc seqs1.msa_consensus.fa >seqs1.msa_consensus.stdout 2>seqs1.msa_consensus.stderr"
${SEQCONS} -m msa_consensus -i seqs1.fa -oa seqs1.msa_consensus.sam -oc seqs1.msa_consensus.fa >seqs1.msa_consensus.stdout 2>seqs1.msa_consensus.stderr
echo "=> $?"
echo "${SEQCONS} -m nop -i seqs1.fa -oa seqs1.nop.sam -oc seqs1.nop.fa >seqs1.nop.stdout 2>seqs1.nop.stderr"
${SEQCONS} -m nop -i seqs1.fa -oa seqs1.nop.sam -oc seqs1.nop.fa >seqs1.nop.stdout 2>seqs1.nop.stderr
echo "=> $?"

# All consensus variants for alignments.

echo "${SEQCONS} -m msa_consensus -i alns1.sam -oa alns1.msa_consensus.sam -oc alns1.msa_consensus.fa >alns1.msa_consensus.sam.stdout 2>alns1.msa_consensus.sam.stderr"
${SEQCONS} -m msa_consensus -i alns1.sam -oa alns1.msa_consensus.sam -oc alns1.msa_consensus.fa >alns1.msa_consensus.sam.stdout 2>alns1.msa_consensus.sam.stderr
echo "=> $?"
echo "${SEQCONS} -m contig_consensus -i alns1.sam -oa alns1.contig_consensus.sam -oc alns1.contig_consensus.fa >alns1.contig_consensus.sam.stdout 2>alns1.contig_consensus.sam.stderr"
${SEQCONS} -m contig_consensus -i alns1.sam -oa alns1.contig_consensus.sam -oc alns1.contig_consensus.fa >alns1.contig_consensus.sam.stdout 2>alns1.contig_consensus.sam.stderr
echo "=> $?"
echo "${SEQCONS} -m pos_consensus -i alns1.sam -oa alns1.pos_consensus.sam -oc alns1.pos_consensus.fa >alns1.pos_consensus.sam.stdout 2>alns1.pos_consensus.sam.stderr"
${SEQCONS} -m pos_consensus -i alns1.sam -oa alns1.pos_consensus.sam -oc alns1.pos_consensus.fa >alns1.pos_consensus.sam.stdout 2>alns1.pos_consensus.sam.stderr
echo "=> $?"
echo "${SEQCONS} -m realign -i alns1.sam -oa alns1.realign.sam -oc alns1.realign.fa >alns1.realign.sam.stdout 2>alns1.realign.sam.stderr"
${SEQCONS} -m realign -i alns1.sam -oa alns1.realign.sam -oc alns1.realign.fa >alns1.realign.sam.stdout 2>alns1.realign.sam.stderr
echo "=> $?"
echo "${SEQCONS} -m nop -i alns1.sam -oa alns1.nop.sam -oc alns1.nop.fa >alns1.nop.sam.stdout 2>alns1.nop.sam.stderr"
${SEQCONS} -m nop -i alns1.sam -oa alns1.nop.sam -oc alns1.nop.fa >alns1.nop.sam.stdout 2>alns1.nop.sam.stderr
echo "=> $?"

echo "${SEQCONS} -m msa_consensus -i alns1.sam -oa alns1.msa_consensus.txt >alns1.msa_consensus.txt.stdout 2>alns1.msa_consensus.txt.stderr"
${SEQCONS} -m msa_consensus -i alns1.sam -oa alns1.msa_consensus.txt >alns1.msa_consensus.txt.stdout 2>alns1.msa_consensus.txt.stderr
echo "=> $?"
echo "${SEQCONS} -m contig_consensus -i alns1.sam -oa alns1.contig_consensus.txt >alns1.contig_consensus.txt.stdout 2>alns1.contig_consensus.txt.stderr"
${SEQCONS} -m contig_consensus -i alns1.sam -oa alns1.contig_consensus.txt >alns1.contig_consensus.txt.stdout 2>alns1.contig_consensus.txt.stderr
echo "=> $?"
echo "${SEQCONS} -m pos_consensus -i alns1.sam -oa alns1.pos_consensus.txt >alns1.pos_consensus.txt.stdout 2>alns1.pos_consensus.txt.stderr"
${SEQCONS} -m pos_consensus -i alns1.sam -oa alns1.pos_consensus.txt >alns1.pos_consensus.txt.stdout 2>alns1.pos_consensus.txt.stderr
echo "=> $?"
echo "${SEQCONS} -m realign -i alns1.sam -oa alns1.realign.txt >alns1.realign.txt.stdout 2>alns1.realign.txt.stderr"
${SEQCONS} -m realign -i alns1.sam -oa alns1.realign.txt >alns1.realign.txt.stdout 2>alns1.realign.txt.stderr
echo "=> $?"
echo "${SEQCONS} -m nop -i alns1.sam -oa alns1.nop.txt >alns1.nop.txt.stdout 2>alns1.nop.txt.stderr"
${SEQCONS} -m nop -i alns1.sam -oa alns1.nop.txt >alns1.nop.txt.stdout 2>alns1.nop.txt.stderr
echo "=> $?"

