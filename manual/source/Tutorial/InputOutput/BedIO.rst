.. sidebar:: ToC

    .. contents::

.. _tutorial-io-bed-io:

BED I/O
=======

Learning Objective
  In this tutorial, you will learn how to read and write BED files.

Difficulty
  Average

Duration
  45min

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-io-input-output-overview`, `BED Format Specification <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_

This tutorial shows how to read and write BED files using the :dox:`BedFileIn` and :dox:`BedFileOut` classes.
It starts out with a quick reminder on the structure of BED files and will then continue with how to read and write BED files.

Originally, the BED format was designed for storing annotation tracks in the UCSC genome browser.
Such an annotation track consists of multiple annotation records.
Each annotation adds some meta information to a genomic interval (an interval with begin/end position on a contig/chromosome) The original specification of the format can be found in the `UCSC Genome Browser FAQ <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_.

The BED format is a TSV format and contains 12 columns.
The first three columns specify a genomic region (contig/chromsome name, begin, and end position) and the remaining columns contain additional information.
The full format will be described below.

Since genomic intervals are very useful and because there were many tools for manipulating BED files (sorting, intersecting intervals etc.), many other authors and projects created variants of the BED format.
Usually, three or more columns have the same meaning as in BED and are followed by other, arbitrary tab-separated columns with additional annotation information.
The "full" BED format is then called BED12. BED3, BED4, BED5, and BED6 use the first 3-6 columns and keep the remaining information as data.

BED files can be manipuluated using standard Unix tools such as ``sed``, ``awk``, and ``sort``.
There also is the `bedtools <https://code.google.com/p/bedtools/>`_ suite with additional functionality.

The SeqAn module ``bed_io`` allows the reading and writing of BED files.

BED Format
----------

The following is an example of a BED file:

    .. literalinclude:: ../../../../demos/tutorial/bed_io/example.bed

The meaning of the columns are as follows:

ref (1)
  Name of the reference sequence.

beginPos (2)
  Begin position of the interval.

endPos (3)
  End position of the interval.

name (4)
  Name of the interval.

score (5)
  A score, could also be in scientific notation or several values in a comma/colon-separated list.

strand (6)
  The strand of the feature, ``+`` for forward, ``-`` for reverse, ``.`` for unknown/dont-care.

thickBegin (7)
  Begin position where the feature is drawn thick in the UCSC browser.

thickEnd (8)
  End position where the feature is drawn thick in the UCSC browser.

itemRgb (9)
  Comma-separated triple with RGB values (0..255 each)

blockCount (10)
  The number of blocks (exons) in the BED line (for the UCSC browser).

blockStarts (11)
  Comma-separated list with begin positions of exons (for the UCSC browser, should be consistent with ``blockCount``).

blockSizes (12)
  Comma-separated list with exon lists (for the UCSC browser, should be consistent with ``blockCount``).

.. tip::

   1-based and 0-based positions.

   There are two common ways of specifying intervals.

   #. Start counting positions at 1 and give intervals by the first and last position that are part of the interval (closed intervals).
      For example, the interval ``[1000; 2000]`` starts at character 1000 and ends at character 2000 and includes it.
      This way is natural to non-programmers and used when giving coordinates in GFF files or genome browsers such as UCSC Genome Browser and IGV.
   #. Start counting positions at 0 and give intervals by the first position that is part of the interval and giving the position behind the last position that is part of the interval.
      The interval from above would be ``[999; 2000)`` in this case.

   In text representations, such as GFF and GTF, 1-based closed intervals are used whereas in the internal binary data structures, SeqAn uses 0-based half-open intervals.
   BED is a text format using 0-based positions.

A First Working Example
-----------------------

The following example shows an example of a program that reads the file ``example.bed`` and prints its contents back to the user on standard output.

.. includefrags:: demos/tutorial/bed_io/example1.cpp

The program first opens a :dox:`BedFileIn` for reading and a :dox:`BedFileOut` for writing.
The BED records are read into :dox:`BedRecord` objects which we will focus on below.
In this case, we use the :dox:`Bed3Record` specialization of the :dox:`BedRecord` class.

The output of the example program looks as follows:

.. includefrags:: demos/tutorial/bed_io/example1.cpp.stdout

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Reproduction

   Objective
     Create a file with the sample BED content from above and adjust the path ``"example.bed"`` to the path to your BED file (e.g. ``"/path/to/my_example.bed"``).

   Solution
      .. container:: foldable

         .. includefrags:: demos/tutorial/bed_io/solution1.cpp


Accessing the Records
---------------------

The class :dox:`BedRecord` stores one record in a BED file.
Note that there are various specializations, each storing a different number of fields.
We show the quasi-definition of :dox:`BedRecord` below.
The other specializations have less fields.

.. includefrags:: demos/tutorial/bed_io/base.cpp
      :fragment: bedRecord

The static members ``INVALID_POS`` and ``INVALID_REFID`` store sentinel values for marking positions and reference sequence ids as invalid.

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

        .. includefrags:: demos/tutorial/bed_io/solution2.cpp

        The output is

        .. includefrags:: demos/tutorial/bed_io/solution2.cpp.stdout

Creating a New File
-------------------

Assignment 3
""""""""""""

.. container:: assignment

   Generating BED From Scratch

   Type
     Application

   Objective
     Write a program that prints the following BED file.
     Create ``BedRecord<Bed6>`` objects and write them to a ``BedFileOut`` using ``writeRecord()``.

     .. includefrags:: demos/tutorial/bed_io/solution3.cpp.stdout

   Solution
    .. container:: foldable

       .. includefrags:: demos/tutorial/bed_io/solution3.cpp

Next Steps
----------

* Continue with the :ref:`tutorial`.
