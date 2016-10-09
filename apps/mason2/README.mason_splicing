                                Mason Splicing
                        Transcript Splicing Simulator

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

Mason Splicing allows one to construct transcripts from a GFF/GTF file, a
reference sequence and an optional variant file.

------------------------------------------------------------------------------
2. Examples
------------------------------------------------------------------------------

You can find the binary "mason_splicing" in the directory "bin" and the
example file in the directory "examples".

------------------------------------------------------------------------------
2.1 Help
------------------------------------------------------------------------------

The command:

  $ mason_splicing --help

prints the help for Mason Splicing.

------------------------------------------------------------------------------
2.2 Transcript Simulation
------------------------------------------------------------------------------

We can use mason_splicing to simulate the transcript from a GFF file and its
reference.  Note that this would also work with a GTF file:

  $ mason_splicing -ir adeno_virus_b.fa -ig adeno_virus_b.gff -o out.fa \
      --gff-type cds --gff-group-by translation_id
  ...
  $ head out.fa
  >pro.1
  TTATGGCCTGGGGCGTTTACAGCTCAAGTCCAAAGGTTGCCCAGACTCGTTAAGCAAGTCCTCGATACAT
  TCCACAGCCTGGCGACGCCCACCAACTCTCACGGCAACTGGTTTAATGGGGCACAGCGGGACCACCGGGT
  GTATCTCAGGAGGTGTGTTAGAAGGACCGGAGTCACAGCTATCCGTACTACTATTGCATTCTCTAGACAC
  AGGTGATGTCGGGCGTCTCAGGATAGCAGGCACCAATTTAGGACGCCGGGTAGGTCTTGCAGGCTCCGGT
  TCTGGCTCGGGCTCAGGCTCAGGTTCAGACACAGGACCTTTTAAAAAATCACAATACAAAATTCTTTAAA
  CCACAAAACTGTAAAAATTAAAAAAAAATTACCACACCAAACCCACCACTCTATCACCGACTGCCCATAA
  TTTTCACTTACTGTAGACAAACATGCCACAGGTCCTCATATAGCAAAGCGAACACATAATATCTGGGTCC
  CCCGTATTCCTCCGGTGATAATGACAAGACCTGCAACCGTGCCCGGGGTGCTCCACATAATCTAACACAA
  ACTCCTCACCCTCTTCATCCTCGTCGTCACTGGGTGGAAAGCCAGCCTCGTGGCAGGTAAGATCGATCAC

Besides the path to the input FASTA and GFF file and the output FASTA file, we
also give the parameters --gff-type and --gff-group-by.

The GFF type is the value of the third column (in both GFF and GTF) to filter
for.  Here, we use "cds", another value could be "exon".

The GFF group-by key is the name of a key/tag from the last column of the
GTF/GTF file.  Here, this is transcript_id, for your data this could also be
"Parent", for example.

mason_splicing will read in the reference FASTA file and the GFF file.  For
each group of GFF/GTF records with the same "group-by" key, it will generate
one sequence in the output FASTA file.  For this, it will take all records
with the type given by --gff-type and concatenate the sequence from these
features in the order that they occur in in the genome.

------------------------------------------------------------------------------
3. Reference and Contact
------------------------------------------------------------------------------

In case of questions and problems please contact the mailing list

  https://lists.fu-berlin.de/listinfo/seqan-dev

or file a bug at

  http://trac.seqan.de/newticket

