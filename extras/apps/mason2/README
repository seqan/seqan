                                    Mason
                   Tools for Biological Sequence Simulation

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

------------------------------------------------------------------------------
1. Overview
------------------------------------------------------------------------------

Mason is a collection of tools for the simulation of biological sequences.
The collection consists of the following tools.  For each tool, there is file
README.${tool_name} that contains detailed documentation of the tool itself.
Also you can get detailed documentation of a tool's parameters by calling it
with the "--help" argument.

 * mason_frag_sequencing
       Simulation of fragment sampling from a genome.
 * mason_genome
       Simulation of random genomic sequences.
 * mason_materializer
       Apply the variation from a VCF file to a genome in a FASTA file.
 * mason_methylation
       Simulate methylation levels for a genome dependent on the context for
       each possible site.
 * mason_simulator
       Simulate NGS reads given a genome and optionally also a VCF file with
       variants for a given donor to use as the source.
 * mason_splicing
       Compute the transcriptome from a genome FASTA file and a GFF file with
       the genes.
 * mason_variator
       Simulate SNPs, small indels, and structural variants for genomic data.
       The result is written out as a VCF file.  Optionally, the resulting
       modified sequence can also be written out.

------------------------------------------------------------------------------
2. Examples
------------------------------------------------------------------------------

The binaries mason_* can be found in the directory "bin".

------------------------------------------------------------------------------
2.1 Help
------------------------------------------------------------------------------

Each program has a verbose built-in help that you can view using the "--help"
option.  For example,

  $ mason_genome --help

prints the help for Mason Genome.

------------------------------------------------------------------------------
2.1 Simulation of Variants and Reads
------------------------------------------------------------------------------

In the following, we will give a quick example of how to simulate variants
into a genome, write the variants into a VCF file and then use this VCF file
as the input of the read simulator to simulate reads of the genome with
variants.  Note that there are separated README files for the mason_variator
and the mason_simulator programs.

Let us first simulate variants into the FASTA sequence adeno_virus.fa (from
the "example" folder).  The variants are written to adeno_out.vcf.

  $ mason_variator -ir adeno_virus.fa -ov adeno_out.vcf

Next, we simulate 100 read pairreads from the input genome and the variants
VCF file:

  $ mason_simulator -ir adeno_virus.fa -n 100 -o left.fq -or right.fq

------------------------------------------------------------------------------
3. Reference and Contact
------------------------------------------------------------------------------

In case of questions and problems please contact the mailing list

  https://lists.fu-berlin.de/listinfo/seqan-dev

or file a bug at

  http://trac.seqan.de/newticket


