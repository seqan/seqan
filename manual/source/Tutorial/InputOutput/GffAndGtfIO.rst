.. sidebar:: ToC

    .. contents::

.. _tutorial-io-gff-and-gtf-io:

GFF and GTF I/O
===============

Learning Objective
  In this tutorial, you will learn how to read and write GFF and GTF files.

Difficulty
  Average

Duration
 45 min

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-io-input-output-overview`, `GFF Format Specification <http://genome.ucsc.edu/FAQ/FAQformat.html#format3>`_

This tutorial shows how to read and write GFF and GTF files using the :dox:`GffFileIn` and :dox:`GffFileOut` classes.
It starts out with a quick reminder on the structure of GFF and GTF files and will then continue with how to read and write GFF and GTF files.

The GFF and GTF formats are used for annotating genomic intervals (an interval with begin/end position on a contig/chromosome).
GFF exists in versions 2 and 3 and GTF is sometimes called "GFF 2.5".
There are specifications for `GFF 2 <http://www.sanger.ac.uk/resources/software/gff/spec.html>`_, `GFF 3 <http://www.sequenceontology.org/gff3.shtml>`_, and `GTF <http://mblab.wustl.edu/GTF22.html>`_ available elsewhere.
GFF and GTF are TSV-based formats and in general have the same structure.
The main difference is the underlying system/ontology for the annotation but also smaller differences in the format.

In this tutorial, we will focus on the format GFF 3 since it is the most current one with most complete tool support.
The information of this tutorial can easily be translated to the other two formats.

The SeqAn module ``gff_io`` allows the reading and writing of the GFF and GTF formats.

.. tip::

    Format Version Support in SeqAn

    :dox:`GffFileIn` allows to read GFF files in version 2 and 3 and GTF files.
    For writing, :dox:`GffFileOut` supports only GFF 3 and GTF.

GFF Format
----------

The following is an example of a GFF 3 file:

.. includefrags:: demos/tutorial/gff_and_gtf_io/example.gff

The meaning of the columns are as follows:

seq id (1)
  Name of the reference sequence.

source (2)
  Free text field describing the source of the annotation, such as a software (e.g. "Genescan") or a a database (e.g. "Genebank"), "``.``" for none.

type (3)
  The type of the annotation.

start (4)
  The 1-based begin position of the annotation.

end (5)
  The 1-based end position of the annotation.

score (6)
  The score of the annotation, "``.``" for none.

strand (7)
  The strand of the annotation, "``+``" and "``-``" for forward and reverse strand, "``.``" for features that are not stranded.

phase (8)
  Shift of the feature regarding to the reading frame, one of "``0``", "``1``", "``2``", and "``.``" for missing/dont-care.

attributes (9)
  A list of key/value attributes.
  For GFF 3, this is a list of ``key=value`` pairs, separated by semicolons (e.g. ``ID=cds00003;Parent=mRNA00003``).
  For GTF and GFF 2, this is a list of tuples, separated by semicolon.
  The first entry gives the key, the following entries are values.
  Strings are generally enclosed in quotes (e.g. ``Target "HBA_HUMAN" 11 55 ; E_value 0.0003``)

.. tip::

   1-based and 0-based positions.

   There are two common ways of specifying intervals.

   #. Start counting positions at 1 and give intervals by the first and last position that are part of the interval (closed intervals).
      For example, the interval ``[1000; 2000]`` starts at character 1,000 and ends at character 2000 and includes it.
      This way is natural to non-programmers and used when giving coordinates in GFF files or genome browsers such as UCSC Genome Browser and IGV.
   #. Start counting positions at 0 and give intervals by the first position that is part of the interval and giving the position behind the last position that is part of the interval.
      The interval from above would be ``[999; 2000)`` in this case.

   In text representations, such as GFF and GTF, 1-based closed intervals are used whereas in the internal binary data structures, SeqAn uses 0-based half-open intervals.

A First Working Example
-----------------------

The following example shows an example of a program that reads the file ``example.gff`` and prints its contents back to the user on standard output.

.. includefrags:: demos/tutorial/gff_and_gtf_io/example1.cpp

The program first opens a :dox:`GffFileIn` for reading and a :dox:`GffFileOut` for writing.
The GFF records are read into :dox:`GffRecord` objects which we will focus on below.

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Reproduction

   Objective
     Create a file with the sample GFF content from above and adjust the path ``"example.gff"`` to the path to your GFF file (e.g. ``"/path/to/my_example.gff"``).

   Solution
      .. container:: foldable

         .. includefrags:: demos/tutorial/gff_and_gtf_io/solution1.cpp

         .. includefrags:: demos/tutorial/gff_and_gtf_io/solution1.cpp.stdout

Accessing the Records
---------------------

The class :dox:`GffRecord` stores one record in a Gff file.

.. includefrags:: demos/tutorial/gff_and_gtf_io/base.cpp
      :fragment: GffRecord

The static members ``INVALID_POS`` and ``INVALID_REFID`` store sentinel values for marking positions and reference sequence ids as invalid.
The static funtion ``INVALID_SCORE()`` returns the IEEE float "NaN" value.

Assignment 2
""""""""""""

.. container:: assignment

   Counting Records

   Type
     Review

   Objective
     Change the result of `Assignment 1`_ by counting the number of variants for each chromosome/contig instead of writing out the records.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/gff_and_gtf_io/solution2.cpp

        The output is

        .. includefrags:: demos/tutorial/gff_and_gtf_io/solution2.cpp.stdout


Creating a New File
-------------------

Assignment 3
""""""""""""

.. container:: assignment

   Generating GFF From Scratch

   Type
     Application

   Objective
     Write a program that prints the following GFF file.
     Create ``GffRecord`` objects and write them to a ``GffFileOut`` using ``writeRecord()``.

     .. includefrags::demos/tutorial/gff_and_gtf_io/solution3.cpp.stdcout

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/gff_and_gtf_io/solution3.cpp

Next Steps
----------

* Continue with the :ref:`tutorial`.
