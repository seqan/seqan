.. We create this roles for putting the "Introduction: etc. headings
   on this page without them displaying in the ToC.  This would break
   rendering the ToC correctly on readthedocs style.  The rubric
   headings are formatted using CSS.
.. role:: rubric-heading1
   :class: rubric-heading1
.. role:: rubric-heading2
   :class: rubric-heading2

.. _tutorial:

Tutorial
--------

.. The inclusion order in the following determines the order in the table
   of contents to the left.
.. toctree::
   :hidden:
   :maxdepth: 2

   Tutorial/GettingStarted
   Tutorial/FirstStepsInSeqAn
   Tutorial/BackgroundAndMotivation

   Tutorial/Sequences
   Tutorial/Alphabets
   Tutorial/StringSets
   Tutorial/SequencesInDepth

   Tutorial/Iterators

   Tutorial/AlignmentRepresentation
   Tutorial/PairwiseSequenceAlignment
   Tutorial/MultipleSequenceAlignment

   Tutorial/Indices
   Tutorial/IndexIterators
   Tutorial/IndexQGram

   Tutorial/PatternMatching

   Tutorial/Graphs

   Tutorial/BasicSequenceIO
   Tutorial/IndexedFastaIO
   Tutorial/BasicSamBamIO
   Tutorial/VcfIO
   Tutorial/BedIO
   Tutorial/GffGtfIO

   Tutorial/Modifiers

   Tutorial/Randomness

   Tutorial/SeedAndExtend

   Tutorial/ParsingCommandLineArguments

   Tutorial/AnnotationStore

   Tutorial/InputOutputOverview
   Tutorial/SequenceFileIO
   Tutorial/SamBamIO
   Tutorial/FileIO
   Tutorial/LexicalCasting
   Tutorial/Parsing

   Tutorial/FragmentStore
   Tutorial/ConsensusAlignment
   Tutorial/Realignment

   Tutorial/SimpleRnaSeq
   Tutorial/SimpleReadMapping
   Tutorial/MiniBowtie
   Tutorial/JournalSet
   Tutorial/KnimeNode

   Tutorial/BasicTechniques
   Tutorial/Metafunctions
   Tutorial/TemplateSubclassing
   Tutorial/GlobalFunctionInterface

   Tutorial/Basics
   Tutorial/WritingTests

The SeqAn tutorials are the best way to get started with learning how to develop using SeqAn.
In contrast, the `API Documentation <http://docs.seqan.de/>`_ gives more comprehensive but less verbose documentation about the library while the How-Tos are strictly task driven and narrower in scope.

The main audience of the tutorials are graduate students and professionals who want to learn how to use SeqAn.
Previous programming knowledge is required, knowledge of C++ is recommended.

.. rubric:: :rubric-heading1:`Introduction`

These tutorials show you how to get started with SeqAn, including the :ref:`installation <tutorial-getting-started>`.
Then, you can learn about the background and motivation of SeqAn.
You should then definitely start your engines and read the :ref:`tutorial-first-steps-in-seqan` tutorial to see an example highlighting many important concepts in the SeqAn library.

:ref:`tutorial-getting-started`
  This tutorial will walk you through the installation of SeqAn and its dependencies.
  Then, you will create your first minimal SeqAn application!

:ref:`tutorial-first-steps-in-seqan`
  This tutorial gives practical examples and applications of the most important basic techniques.
  You should read this tutorial if you are starting out with SeqAn.

:ref:`tutorial-background-and-motivation`
  This tutorial gives an overview over the design aims and principles of SeqAn and a motivation for the employed mechanisms.

We highly recommend you to follow the Getting Started instructions if you are starting out with SeqAn.
Note that it is also possible to use SeqAn strictly as a library with your own build system.
The article :ref:`build-manual-integration-with-your-own-build-system` contains detailed information about this.

.. rubric:: :rubric-heading1:`A Stroll Through SeqAn`

.. rubric:: :rubric-heading2:`Sequences`

:ref:`tutorial-sequences`
  This tutorial introduces you to the basics of fundamental concept of sequences, namely Strings and Segments.

:ref:`tutorial-alphabets`
  This tutorial introduces you to SeqAn's alphabets, or in other words, the contained types of sequences.

:ref:`tutorial-string-sets`
  StringSets This tutorial introduces you to SeqAn's ``StringSet``, an efficient data structure to store a set of sequences.

:ref:`tutorial-sequences-in-depth`
  In this tutorial you will learn how to optimize the work with sequences, using different specializations of Strings and different overflow strategies for capacity changes.

.. rubric:: :rubric-heading2:`Iterators`

:ref:`tutorial-iterators`
  This tutorial explains how to use iterators in SeqAn, illustrated on containers.

.. rubric:: :rubric-heading2:`Alignments`

:ref:`tutorial-alignment-representation`
  This section of the tutorial introduces you to the data structures that are used to represent alignments in SeqAn.

:ref:`tutorial-pairwise-sequence-alignment`
  In this part of the tutorial we demonstrate how to compute pairwise sequence alignments in SeqAn.
  It shows the use of different scoring schemes, and which parameters can be used to customize the alignment algorithms.

:ref:`tutorial-multiple-sequence-alignment`
  In the last section of this tutorial we show how to compute multiple sequence alignments in SeqAn using a scoring matrix.

.. rubric:: :rubric-heading2:`Indices`

:ref:`tutorial-indices`
  This tutorial introduces you to the various indices in SeqAn like extended suffix arrays or k-mer indices.

:ref:`tutorial-index-iterators`
  This tutorial introduces you to the various index iterators with which you can use indices as if traversing search trees or tries.

:ref:`tutorial-q-gram-index`
  This tutorial introduces you to SeqAn's q-gram index.

.. rubric:: :rubric-heading2:`Pattern Matching`

:ref:`tutorial-pattern-matching`
  This section of the tutorial introduces you to the algorithms in SeqAn for exact and approximate pattern matching.

.. rubric:: :rubric-heading2:`Graphs`

:ref:`tutorial-graphs`
  This section of the tutorial introduces you to the graph type in SeqAn.
  We will discuss the various graph specializations and show you how to create directed and undirected graphs as well as HMMs, how to store additional information for edges and vertices and last but not least how to apply standard algorithms to the graphs.

.. rubric:: :rubric-heading2:`I/O Basics`

:ref:`tutorial-basic-sequence-io`
  Basic Sequence I/O This tutorial explains how to use the high-level API for reading and writing sequence files.

:ref:`tutorial-indexed-fasta-io`
  Indexed FASTA I/O This tutorial explains how to use FASTA index files for quick random access within FASTA files: read contigs or just sections without having to read through whole FASTA file.

:ref:`tutorial-basic-sam-bam-io`
  Basic SAM and BAM I/O This tutorial explains how to use the high-level API for reading and writing SAM and BAM files.

:ref:`tutorial-vcf-io`
  VCF I/O This tutorial explains how to use the high-level API for reading and writing VCF files.

:ref:`tutorial-bed-io`
  BED I/O This tutorial explains how to use the high-level API for reading and writing BED files.

:ref:`tutorial-gff-and-gtf-io`
  GFF and GTF I/O This tutorial explains how to use the high-level API for reading and writing GFF and GTF files.

.. rubric:: :rubric-heading2:`Modifiers`

:ref:`tutorial-modifiers`
  Modifiers Modifiers can be used to change the elements of a container without touching them.
  Here you will see, what modifiers are available in SeqAn.

.. rubric:: :rubric-heading2:`Randomness`

:ref:`tutorial-randomness`
  This chapter shows module random that provides pseudo random number generation functionality.

.. rubric:: :rubric-heading2:`Seed-And-Extend`

:ref:`tutorial-seed-and-extend`
  In this part of the tutorial we will introduce SeqAn's seed class, demonstrate seed extension and banded alignment with seeds, and finally show the usage of seed chaining algorithms.

.. rubric:: :rubric-heading2:`Parsing Command Line Arguments`

:ref:`tutorial-parsing-command-line-arguments`
  Parsing Command Line Arguments In this tutorial, you will learn how to use the  ArgumentParser class for parsing command line arguments.

.. rubric:: :rubric-heading2:`Genome Annotations`

:ref:`tutorial-genome-annotations`
  You will learn how to work with annotations in SeqAn and analyzing them, using the :dox:`FragmentStore::annotationStore` which is part of SeqAn's :dox:`FragmentStore`.

.. rubric:: :rubric-heading2:`More I/O`

These tutorials explain how to use the I/O functionality in SeqAn beyond the basic sequence, SAM/BAM and indexed FASTA I/O from above.
The tutorials are targeted at developers that either want to use the lower level I/O routines in SeqAn or write their own parsers.
We recommended to start out reading the I/O Overview and then jump to the chapter that interests you most.
In this Section we introduce the three main techniques of programming in SeqAn, namely the ''global function interface'', the use of
''Metafunctions'', and the concept of  ''Template subclassing''.

:ref:`tutorial-input-output-overview`
  This article gives an overview of the I/O functionality in SeqAn.

After reading, you will have a better understanding of the different bits in this section of the library.
The following tutorials introduce the lower level I/O routines for specific file formats.

:ref:`tutorial-sequence-file-io`
  This tutorial explains the RecordReader- and Stream-based interface for reading sequence files.

:ref:`tutorial-sam-bam-io`
  This tutorial explains the lower level API for reading and writing SAM and BAM files.

Read the following tutorials to learn how to write your own I/O routines.

:ref:`tutorial-file-io`
  This chapter shows how to use the file I/O facilities of SeqAn, including streams, compressed streams and memory mapped files.

:ref:`tutorial-lexical-casting`
  This tutorial explains the :dox:`lexicalCast` and :dox:`lexicalCast2` functions that allow to convert strings representing numbers into their numeric values.

:ref:`tutorial-parsing`
  In this part of the tutorial, you will be introduced to the parsing and tokenizing functionality using the RecordReader class.
  You will get the necessary information to write your own file parsers.

.. rubric:: :rubric-heading1:`Advanced Tutorials`

:ref:`tutorial-fragment-store`
  This tutorial shows how to use the fragment store which is a database for read mapping, sequence assembly or gene annotation.
  It supports to read/write multiple read alignments in SAM or AMOS format and access and modify them.
  It supports to read/write gene annotations in GFF/GTF and UCSC format, to create custom annotation types, and to traverse and modify the annotation tree.

:ref:`tutorial-consensus-alignment`
  This tutorial describes how to compute consensus alignments from NGS reads or other nucleic sequence, such as transcripts.
  The DNA sequences are stored in a fragment store, such that rough alignment information is available.

:ref:`tutorial-realignment`
  This tutorial describes how to use SeqAn's realignment module for refining multi-read alignment (or other sequences) stored in a fragment store.

:ref:`tutorial-simple-rna-seq`
  In this tutorial you will learn how to implement a simple RNA-Seq based gene quantification tool, that computes RPKM expression levels based on a given genome annotation and RNA-Seq read alignments.

:ref:`tutorial-simple-read-mapping`
  This tutorial shows how to implement a simple read mapping program based on the SWIFT filter and online Hamming finder for verification.

:ref:`tutorial-mini-bowtie`
  Mini-Bowtie is a very basic read aligner that is inspired by the well known Bowtie program :cite:`Langmead2009`.
  It serves as an example to show that you can write sophisticated programs with SeqAn using few lines of code.

:ref:`tutorial-data-journaling`
  In this tutorial we demonstrate how you can handle multiple large sequence in main memory while the data structures themself support a certain parallel sequence analysis.

:ref:`tutorial-knime-nodes`
  Here you can learn how to use SeqAn apps in KNIME.

.. rubric:: :rubric-heading1:`Developer's Corner`

First, congratulations on becoming an offical SeqAn developer!
After you went through the tutorials and before you actually start to develop your own application with SeqAn you might want to learn :ref:`how-to-write-tests` and read about the :ref:`API documentation <style-guide-dox-api-docs>`.
In addition, we follow a SeqAn specific :ref:`style-guide`.
Information like this can be found on the section site.
There are plenty of information completing your knowledge about SeqAn so have a look!

.. rubric:: :rubric-heading2:`Frequently used Software Techniques`

We assume that the user is acquainted with the basic data types of SeqAn, the introductory example and the demo programs.
Also you should be acquainted with the STL and template programming.
In this Section we introduce the three main techniques of programming in SeqAn, namely the *global function interface*, the use of
*Metafunctions*, and the concept of *Template subclassing*.

:ref:`tutorial-basic-techniques`
  Here we remind you of the basics of template programming and the use of the STL.

:ref:`tutorial-metafunctions`
  In this section you find an introductory explanation how Metafunctions are used in SeqAn to obtain information about data types used which will only be instantiated at compile time.

:ref:`tutorial-template-subclassing`
  In this section you find a short example that illustrates the power of template subclassing.

:ref:`tutorial-global-function-interface`
  In this section you find a useful piece of code that shows you the flexibility of the global function interface.
