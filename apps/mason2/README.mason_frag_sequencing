                          Mason Fragment Sequencing

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

Mason Fragment Sequencing is a simple sequence simulator.  Instead of the
end-to-end read simulator "mason_simulator", the program
"mason_frag_sequencing" expects as the input DNA fragment.  For these
fragments, sequencing is then simulated either from one or both ends.

The program thus exposes the NGS sequencing engine of the Mason package
without fragment selection and variant simulation.

------------------------------------------------------------------------------
2. Examples
------------------------------------------------------------------------------

You can find the binary "mason_genome" in the directory "bin".

------------------------------------------------------------------------------
2.1 Help
------------------------------------------------------------------------------

The command:

  $ mason_genome --help

prints the help for Mason Genome.

------------------------------------------------------------------------------
2.1 Simulating Paired-End Reads from Fragments
------------------------------------------------------------------------------

We call mason_frag_sequencing on the file adeno_fragments.fa to simulate
paired-end Illumina reads (default) of length 100.

  $ head adeno_fragments.fa
  >fragment0
  CCTATCTAATAATATACCTTATACTGGACTAGTGCCAATATTAAAATGAAGTGGGCGTAGTGTGTAATTT
  GATTGGGTGGAGGTGTGGCTTTGGCGTGCTTGTAAGTTTGGGCGGATGAGGAAGTGGGGCGCGGCGTGGG
  AGCCGGGCGCGCCGGATGTGACGTTTTAGACGCCATTTTACACGGAAATGATGTTTTTTGGGCGTTGTTT
  GTGCAAATTTTGTGTTTTAGGCGCGAAAACTGAAATGCGGAAGTGAAAATTGATGACGGCAATTTTATTA
  TAGGC
  >fragment1
  ACGGGGAAACTCCACGTTGGCGCTCAAAGGGCGCGTTTATTGTTCTGTCAGCTGATCGTTTGGGTATTTA
  ATGCCGCCGTGTTCGTCAAGAGGCCACTCTTGAGTGCCAGCGAGAAGAGTTTTCTCTGCCAGCTCATTTT
  CACGGCGCCATTATGAGAACTGAAATGACTCCCTTGGTCCTGTCGTATCAGGAAGCTGACGACATATTGG
  $ mason_frag_sequencing -i adeno_fragments.fa -o left.fq -or right.fq
  ...
  $ head left.fq 
  @simulated.1
  GCCTATAATAAAATTGCCGTCATCAATTTTCACTTCCGCATTTCAGTTTTCGCGCCTAAAACACAAAATTTGCACAAACAACGCCCAAAAAACATCATTT
  +
  IIIHHIHIHIIIIIHHHHIIHGIGGIIIHIGIIFGEIIICGIHHIFFIEIFIHCHFFEIIEHFHIIHIAIIDDIGIHIIIHIIHFIIDBFHDHIIIAIDI
  @simulated.2
  GGCAACTCTGGATGGTCTAACTGAAACTCCTCACGTTCCCTATCAGCGGCAGCAGCAGCTGCGGATGCAGAAACATGCGCCATTCCGTTCTCGTCTTGCT
  +
  IHIIHIHHIHIIIHHHHIIHFHHGFIIHGIIIGHIIEGIHIHHHEGGIIIGGIIIIDGBIIIIICIDDEICGIFII=IIEHFDDHIIIIIII;DAGEBI?
  @simulated.3
  AGTCCTGTGAGCACCACCGGAATAGTACTGGAAATACTGACTTAATGTGCTCTTTGTGCTATCTGCGAGCCTACAACATGTTCATTTACAGTAAGTGTGC
  $ head right.fq 
  @simulated.1
  CCTATCTAATAATATACCTTATACTGGACTAGTGCCAATATTAAAATGAAGTGGGCGTAGTGTGTAATTTGATTGGGTGGAGGTGTGGCTTTGGCGTGCT
  +
  HHHHHIIHIHIHHIHHHIHHIIIFIHHHHIIIIGIGHIIIIEHIIIHHIIIBHHIIIHFIG?IIFIHGIGGIGDEBGIE?IHDEIIII?IIIGIC@IAIB
  @simulated.2
  ACGGGGAAACTCCACGTTGGCGCTCAAAGGGCGCGTTTATTGTTCTGTCAGCTGATCGTTTGGGTATTTAATGCCGCCGTGTTCGTCAAGAGGCCACTCT
  +
  HHHIIHHHHIHHHIHGFHIGHIEGGHGIFIIFHHHIIHIGCFEEIGIIIHIIECGFIDIDIEDABFGIIEFHIHAI<IIIHAIEHIH>ICIHIGGIIIIG
  @simulated.3
  TGCTATTAGGTTCAGGCTCATTATCGGAAACAGGACCTAAAAACAACAAAATATTATTTTTCACTGCTTAAGAAAAAAAAATCACCTCCCACCTCCCATA

------------------------------------------------------------------------------
3. Reference and Contact
------------------------------------------------------------------------------

In case of questions and problems please contact the mailing list

  https://lists.fu-berlin.de/listinfo/seqan-dev

or file a bug at

  http://trac.seqan.de/newticket

