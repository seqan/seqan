                                Mason Variator
                         Genomic Variation Simulator

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

Mason Variator is a stand-alone program for the simulation of SNPs, small
indels, and structural variants into genomes.  SNPs and small indels are
simulated on a per-position probability.  Structural variants can be simulated
based on such rates as well, but also given as pairs of types and lengths from
a TSV file.

The simulation result can be written out in a VCF file.  Optionally, the
program can also write out a FASTA file with sequence that includes the
variants and the breakpoints of the variants in the simulated sequence.

------------------------------------------------------------------------------
2. Example
------------------------------------------------------------------------------

You can find the binary "mason_variator" in the directory "bin" and the
example file in the directory "examples".

------------------------------------------------------------------------------
2.1 Help
------------------------------------------------------------------------------

The command:

  $ mason_variator --help

prints the help for Mason Variator

------------------------------------------------------------------------------
2.2 Simple SV Simulation
------------------------------------------------------------------------------

As an introduction, let us simulate some SVs given by their type and length
from a TSV file (sv_sizes.tsv, see Section 4.2) into a small genome
(adeno_virus.fa).  We will write the variants into the VCF file adeno_out.vcf
(see Section 4.1), the simulated sequence into the FASTA file adeno_out.fa,
and the positions of the breakpoints in adeno_out.fa into the TSV file
adeno_out.tsv (see Section 4.3).

Note that for the human genome, you will want to use the option -n /
--num-haplotypes to select the number of haplotypes/allelles to simulate.

  $ mason_variator -ir adeno_virus.fa -ov adeno_out.vcf -of adeno_out.fa \
      -it sv_sizes.tsv --out-breakpoints adeno_out.tsv
  ...

The output to the console should show you that 3 SNPs and 6 structural
variants have been simulated.  Now, have a look at the generated files.

We can now use the script breakpoint_contigs.awk (in folder "bin") to cut out
the regions around the breakpoints in adeno_out.fa and write them into the
file adeno_breakpoints.fa.  Note that this requires you to have samtools
installed.

We cut out the sequence around the breakpoints within a radius of 200bp
towards each side:

  $ awk -f breakpoint_contigs.awk -v context=200 -v ref=adeno_out.fa \
      adeno_out.tsv > breakpoint_contigs.fa
  ...

These contigs around breakpoints could now for example be aligned to the
original genome using STELLAR and the result then analyzed using GUSTAF:

  $ stellar adeno_virus.fa breakpoint_contigs.fa -o adeno_out.gff
  ...
  $ gustaf -st 1 -m adeno_out.gff adeno_virus.fa breakpoint_contigs.fa
  ...

------------------------------------------------------------------------------
3. Reference and Contact
------------------------------------------------------------------------------

In case of questions and problems please contact the mailing list

  https://lists.fu-berlin.de/listinfo/seqan-dev

or file a bug at

  http://trac.seqan.de/newticket

------------------------------------------------------------------------------
4. File Formats
------------------------------------------------------------------------------

The program mason_variator uses the FASTA file for input and output sequences.

It writes out the variation in the VCF format with almost no special tags or
values.  See 4.1 for a description of the INFO fields that are not alread
described in the VCF 4.1 standard.

Sections 4.2 and 4.3 describe the TSV file format for specifying variations to
simulate and the TSV format that gives the breakpoints in haplotype
coordinates.

------------------------------------------------------------------------------
4.1 VCF Format
------------------------------------------------------------------------------

mason_variator creates VCF 4.1 files for the simulated variants.  The format
is described on the website of the 1000 genomes project [1].  mason_variator
introduces an INFO key called TARGETPOS that allows compact representations of
duplications.

SNPs are simulated such that all simulated haplotypes can contain the
alternative base or the original base but at least one haplotype must have the
SNP of course.

Small indels and structural variants are simulated such that only one
haplotype can have the indel.  No SNP or structural variant can overlap or be
within one base with an SV breakpoint.  SVs cannot overlap each other or be
nested.  Recombination can only happen within one contig/chromosome.  These
restrictions make simulation easier and yet are sufficient for the
experimental evaluation of SV tools.

SNPs, small indels, and SV indels are given directly in the REF/ALT columns.
SV indels are also marked as SVTYPE={INS,DEL}.  Duplications are given as
<DUP> in the ALT column, inversions as <INV>.  Translocations are given as 6
breakend BND entries, 2 for each of the breakpoints.

[1] http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41

------------------------------------------------------------------------------
4.2 Variation Selection TSV File
------------------------------------------------------------------------------

The user can give the path to a TSV file with the --in-variant-tsv option.
The first two columns are mandatory.  The first column gives the variation
type to simulate and the second column gives the length of the variant.
Comment lines begin with a hash (#).  For insertions, you can also give a
third column with the sequence to insert.  Note that if this sequence is given
then the length of the sequence takes precedence over the the length in the
second column.

The following variation types are interpreted:

 * INS: An insertion.
 * DEL: A deletion.
 * INV: An inversion.
 * CTR: An intra-chromosomal translocation.
 * DUP: A duplication

An example for the content of such a file is given below:

  #type length
  INS   100
  DEL   1000
  CTR   400
  INS   20     CGATTTAGATATACGTATAC

This file makes the program simulate a 100bp insertion with random characters,
a 1000bp deletion, a 400bp intra-chromosomal translocation, and a 20bp
insertion with specified sequence.  Note that when the --in-variant-tsv
parameter is given, structural variants are not simulated on a per-position
rate.

------------------------------------------------------------------------------
4.3 Breakpoint TSV File
------------------------------------------------------------------------------

The VCF file written by mason_variator gives the coordinates in the reference.
However, it might also be of interest to easily retrieve the coordinates of
breakpoints on the haplotypes with variants.

For this, you can give the parameter --out-breakpoints with a path to a TSV
file to write.  This file contains three columns, one for the contig/reference
name, one for the variant id, and one for the (1-based) position on the
contig.

The following shows an example with 4 breakpoints, two on each haplotype,
spread over three chromosomes.  Note that there are two breakpoints for the
insertion on chrI, also two for the inversion on chrII and one breakpoint for
the deletion on chrIII.  Missing variant names are given with dots (".").

  #ref     id             pos
  chrI/1   sim_sv_indel_1 2051
  chrI/1   sim_sv_indel_1 2351
  chrII/0  sim_inv_1      2451
  chrII/0  sim_inv_2      3573
  chrIII/1 sim_sv_indel_2 9297
