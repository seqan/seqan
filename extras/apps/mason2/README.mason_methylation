                              Mason Methylation
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

Mason Methylation is a program for the simulation of methylation levels of
genomes.  The methylation levels are simulated dependin on the context of the
cytamine.  The levels have a resolution of 2.5% and are encoded as ASCII in
extra methylation FASTA files.

------------------------------------------------------------------------------
2. Examples
------------------------------------------------------------------------------

You can find the binary "mason_methylation" in the directory "bin" and the
example file in the directory "examples".

------------------------------------------------------------------------------
2.1 Help
------------------------------------------------------------------------------

The command:

  $ mason_methylation --help

prints the help for Mason Methylation.

------------------------------------------------------------------------------
2.2 Generating Methylation Levels
------------------------------------------------------------------------------

Let us generate some methylation levels for the file adeno_virus.fasta in the
example directory:

  $ mason_methylation -i adeno_virus.fa -o adeno_virus_m.fa
  ...

Let us look at the result:

  $ head adeno_virus_m.fa
  >gi|56160436|ref|AC_000005.1|/TOP
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!r!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!r!!!!!!!!!!!!!!!!!r!!!!!!!!!!!!!!!!r!r!!!!!!!!
  !!!r!!!r!r!!r!!!!!!!!r!!!!!!!!r!!!!!!!!!!!r!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!r!r!!!!!!!!!!!!!r!!!!!!!!!!!!!!!!!!q!!!!!!!!!!!!!
  !!!!r!r!!!!!!!!!!!r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!f!!!!!r!
  !!!!!!!!!!!!!!r!!!!!r!!!!!!!!!!r!!!!!!!!!!!!!!!!!!!!!!!!r!!!!!!!!!!!!!
  !!!rr!!p!!!!!r!!!!!!!!!!!!!!!!!!!!!!!r!!r!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!&!!X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!;!!!Z!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!-!!!!!!!r!!!!!!!!!!!!!!!!!!!!!!!r!!!!!!!!!r!!
  $ grep '^>gi' adeno_virus_m.fa
  >gi|56160436|ref|AC_000005.1|/TOP
  >gi|56160436|ref|AC_000005.1|/BOT

As you can see, we generated a FASTA file with with two entries for the
sequence in the input FASTA file.  The sequence identifiers are the one from
the input file with the suffixes "/TOP" and "/BOT".

See 4.1 for details on the file format.

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

------------------------------------------------------------------------------
4.1 Methylation Levels FASTA
------------------------------------------------------------------------------

Methylation levles are stored in FASTA files.  For each sequence entry "SEQ"
in the reference, there are two entries "SEQ/TOP" and "SEQ/BOT" in the
methylation levels FASTA file.  The "/TOP" entry provides the levels for the
forward/top strand, the "/BOT" entry provides the levels for the
bottom/reverse strand.

The levels are encoded in ASCII values.  Given a level L between 0 and 100
the ASCII character for the encoding is given by:

  if (L / 80 + '!' < '>')
      return L / 80 + '!';
  else
      return L / 80 + '!' + 1;

Which means the values 0..100 are projected to the ASCII characters '!'...'p'
but the character '>' is skipped so there is no ambiguity when reading the
FASTA file.
