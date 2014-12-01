                               Mason Simulator
                         Methylation Level Simulator

                      http://www.seqan.de/projects/mason
                              Version 2.0.0-beta1
                                February, 2014

                               Manuel Holtgrewe

------------------------------------------------------------------------------
Table of Contents
------------------------------------------------------------------------------

  1. Overview
  2. Examples
  3. Reference and Contact
  4. File Formats

------------------------------------------------------------------------------
1. Overview
------------------------------------------------------------------------------

Mason Simulator is the read simulator from the Mason package.

The simulator has the following features:

  * Single-end and paired-end data,
  * Illumina, 454, and Sanger error models,
  * simulate haplotypes directly from VCF file,
  * BS-Seq simulation, and
  * low main memory usage

------------------------------------------------------------------------------
1.1 Simulation Process
------------------------------------------------------------------------------

The simulation process is depicted in the following picture:

   [Reference]
       +--------<----- (Optional VCF FILE)
       |
       V
   [Haplotype]
       |
       V
   [Fragment]
       |
       V
   [Single or Paired-End Read] + [Alignment]

The reference is loaded and optionally, variants from a VCF file can be
applied.  You can simulate such variants using mason_variator.  This results
in a haplotype that thus can differ from the reference sequence.

Next, fragments are simulated from the sequence (the molecules simulated by
the "Random Fragmentation" step in NGS protocols).  The fragments can then be
sequenced from either side or both side (single-end or paired-end reads).

When using BS-seq simulation, methylation levels are loaded together with the
reference and the variants are also applied to the VCF file (structural
variants might move parts of the level strings and at SNPs and around
breakpoints, the levels are recomputed since the context may have changed.

------------------------------------------------------------------------------
2. Examples
------------------------------------------------------------------------------

You can find the binary "mason_simulator" in the directory "bin" and the
example file in the directory "examples".

------------------------------------------------------------------------------
2.1 Help
------------------------------------------------------------------------------

The command:

  $ mason_simulator --help

prints the help for Mason Simulator.

------------------------------------------------------------------------------
2.2 Single-End Illumina Simulation
------------------------------------------------------------------------------

Simulation of 1000 single-end Illumina reads of length 150bp (default is
100bp) from a genome.  Write out a FASTA file for the reads an a SAM file for
the alignments.

  mason_simulator -ir adeno_vrius.fa -n 1000 --illumina-read-length 150 \
    -o reads.fa -oa alignments.sam

------------------------------------------------------------------------------
2.3 Paired-End Illumina Simulation, With Variants
------------------------------------------------------------------------------

Simulation of paired-end Illumina sequencing (1000 read pairs) from a genome,
including variation from a VCF file.  The result is written to FASTQ files and
alignments are written to a SAM file.

  mason_simulator -ir genome.fa -iv variants.vcf -n 1000 -o reads_1.fq \
    -or reads_2.fq -oa alignments.sam

------------------------------------------------------------------------------
2.4 Paired-End Illumina Simulation, Without Variants
------------------------------------------------------------------------------

Simulation of mate-pair Illumina sequencing (1000 read pairs) from a genome
without a VCF file.  Write out results as FASTQ file and alignments as SAM
file.

  mason_simulator -ir genome.fa -n 1000 -o reads_1.fq -or reads_2.fq \
    -oa alignments.sam

------------------------------------------------------------------------------
2.5 Single-End 454 Read Simulation
------------------------------------------------------------------------------

Simulation of single-end 454 sequencing (1000 reads) from a genome without a
VCF file, write out alignments as SAM file and reads as FASTQ file.  The read
length follows a normal distribution with mean 800 and standard deviation 200.

Note that we explicitely increas the fragment length otherwise the reads will
be longer than the simulated fragments (which by default follow a normal
distribution with mean 300 and deviation 30).

  mason_simulator -ir genome.fa -n 1000 --seq-technology 454 \
    --454-read-length-mean 800 --454-read-length-stddev 200  \
    --fragment-mean-size 2000 --fragment-size-std-dev 200 \
    -o reads.fq -oa alignments.sam

------------------------------------------------------------------------------
3. Reference and Contact
------------------------------------------------------------------------------

In case of questions and problems please contact the mailing list

  https://lists.fu-berlin.de/listinfo/seqan-dev

or file a bug at

  http://trac.seqan.de/newticket

