.. sidebar:: ToC

   .. contents::


.. _tutorial-bed-io:

BED I/O
=======

Learning Objective
  In this tutorial, you will learn how to use the high-level interface :dox:`BedStream` class to read and write BED files.

Difficulty
  Average

Duration
  45min

Prerequisites
  Exposure to the BED format is useful.

This tutorial deals with how to easily read and write BED files using the :dox:`BedStream` class.
It starts out with a quick reminder on the structure of BED files and will then continue with how to read and write BED files.

Originally, the BED format was designed for storing annotation tracks in the UCSC genome browser.
Such an annotation track consists of multiple annotation records.
Each annotation adds some meta information to a genomic interval (an interval with begin/end position on a contig/chromosome) The original specification of the format can be found in the `UCSC Genome Browser FAQ <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_.

The BED format is a TSV format and contains 12 columns.
The first three column specify a genomic region (contig/chromsome name, begin, and end position) and the remaining columns contain additional information.
The full format will be described below.

Since genomic intervals are very useful and because there were many tools for manipulating BED files (sorting, intersecting intervals etc.), many other authors and projects created variants of the BED format.
Usually, three or more columns have the same meaning as in BED and are followed by other, arbitrary tab-separated columns with additional annotation information.
The "full" BED format is then called BED12, and BED3, BED4, BED5, and BED6 use the first 3-6 columns and keep the remaining information as data.

BED files can be manipuluated using standard Unix tools such as ``sed``, ``awk``, and ``sort``.
There also is the `bedtools <https://code.google.com/p/bedtools/>`_ suite with additional functionality.

The SeqAn module ``bed_io`` allows the reading and writing of BED files.

BED File Structure
------------------

The following is an example of a BED file:

::

    chr1    66999824    67210768    NM_032291   0   +   6700004167208778    0   25  227,64,25,72,57,55,176,12,12,25,52,86,93,75,501,128,127,60,112,156,133,203,65,165,2013, 0,91705,98928,101802,105635,108668,109402,126371,133388,136853,137802,139139,142862,145536,147727,155006,156048,161292,185152,195122,199606,205193,206516,207130,208931,
    chr1    48998526    50489626    NM_032785   0   -   4899984450489468    0   14  1439,27,97,163,153,112,115,90,40,217,95,125,123,192,    0,2035,6787,54149,57978,101638,120482,130297,334336,512729,712915,1164458,1318541,1490908,
    chr1    16767166    16786584    NM_018090   0   +   1676725616785385    0   8   182,101,105,82,109,178,76,1248, 0,2960,7198,7388,8421,11166,15146,18170,
    chr1    33546713    33585995    NM_052998   0   +   3354785033585783    0   12  182,121,212,177,174,173,135,166,163,113,215,351,0,275,488,1065,2841,10937,12169,13435,15594,16954,36789,38931,
    chr1    16767166    16786584    NM_001145278    0   +   1676725616785385    0   8   104,101,105,82,109,178,76,1248, 0,2960,7198,7388,8421,11166,15146,18170,

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
      For example, the interval ``[1,000; 2,000]`` starts at character 1,000 and ends at character 2,000 and includes it.
      This way is natural to non-programmers and used when giving coordinates in GFF files or genome browsers such as UCSC Genome Browser and IGV.
   #. Start counting positions at 0 and give intervals by the first position that is part of the interval and giving the position behind the last position that is part of the interval.
      The interval from above would be ``[999; 2,000)`` in this case.

   In text representations, such as GFF and GTF, 1-based closed intervals are used whereas in the internal binary data structures, SeqAn uses 0-based half-open intervals.
   BED is a text format using 0-based positions.

A First Working Example
-----------------------

The following example shows an example of a program that reads the file with the path ``example.bed`` and prints its contents back to the user on stdout.
If you want to try out this program then create a file with the sample BED content from above and adjust the path ``"example.bed"`` in the program below to the path to your BED file (e.g. ``"path/to/my_example.bed"``).

.. includefrags:: extras/demos/tutorial/bed_io/example1.cpp

The program first opens a :dox:`BedStream` for reading, then one for writing.
You can read from stdin and write to stdout using ``"-"`` as the file name.

The member ``sequenceNames`` of your :dox:`BedStream` object ``bedIn`` contains the names of the reference sequences that have been seen in records so far.
This :dox:`StringSet` of :dox:`CharString` thus gets new elements as you read the BED file.
For the translation between reference names and numeric ids, a cache is used.
The function :dox:`BedStream#addSequenceName` can be used to register the sequence name with the ``bedOut`` stream.
This will also update the cache.

The BED records are read into :dox:`BedRecord` objects which we will focus on below.
In this case, we use the :dox:`Bed3Record` specialization of the :dox:`BedRecord` class.

.. tip::

   BED records and additional data.

   This means that the first three columns are read and interpreted and available in the class members.
   The remaining data is stored in the ``data`` member variable of the record.
   This means that the data stored after the first three columns could be empty or of an arbitrary format.

Note that the example above is missing error handling.
This means that if the input format is ill-formed, error return codes are not handled appropriately and the program might do something unexpected in the case of an error.
We will fix this in `Assignment 1`_.

You can see the output of the program below when called with the input file from above.

.. code-block:: console

   chr1    66999824    67210768    NM_032291   0   +   6700004167208778    0   25  227,64,25,72,57,55,176,12,12,25,52,86,93,75,501,128,127,60,112,156,133,203,65,165,2013, 0,91705,98928,101802,105635,108668,109402,126371,133388,136853,137802,139139,142862,145536,147727,155006,156048,161292,185152,195122,199606,205193,206516,207130,208931,
   chr1    48998526    50489626    NM_032785   0   -   4899984450489468    0   14  1439,27,97,163,153,112,115,90,40,217,95,125,123,192,    0,2035,6787,54149,57978,101638,120482,130297,334336,512729,712915,1164458,1318541,1490908,
   chr1    16767166    16786584    NM_018090   0   +   1676725616785385    0   8   182,101,105,82,109,178,76,1248, 0,2960,7198,7388,8421,11166,15146,18170,
   chr1    33546713    33585995    NM_052998   0   +   3354785033585783    0   12  182,121,212,177,174,173,135,166,163,113,215,351,0,275,488,1065,2841,10937,12169,13435,15594,16954,36789,38931,
   chr1    16767166    16786584    NM_001145278    0   +   1676725616785385    0   8   104,101,105,82,109,178,76,1248, 0,2960,7198,7388,8421,11166,15146,18170,

To add error handling, we have to check return values.
The :dox:`BedStream#readRecord` call returns a status code different from ``0``, indicating an error.

In `Assignment 1`_, we will add error handling to the program.

Assignment 1
""""""""""""

.. container:: assignment

   Adding Error Handling

   Type
     Review

   Objective
     Add error handling using the hints below.

   Hints
     The functions :dox:`BedStream#readRecord` and :dox:`BedStream#writeRecord` return a status code ``int``, ``0`` on success, ``1`` on errors.
     The function :dox:`BedStream#isGood` checks whether the state of a :dox:`BedStream` is errorneous.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/bed_io/solution1.cpp

The Class :dox:`BedRecord`
--------------------------

The class :dox:`BedRecord` stores one record in a BED file.
Note that there are various specializations, each storing a different number of fields.
We show the quasi-definition of :dox:`Bed12Record` below.
The other specializations have less fields.

.. code-block:: cpp

   namespace seqan {

   class BedRecord
   {
   public:
       CharString ref;      // reference name
       __int32 rID;         // index in sequenceNames of BedStream
       __int32 beginPos;    // begin position of the interval
       __int32 endPos;      // end position of the interval
       CharString name;     // name of the interval
       CharString score;    // score of the interval
       char strand;         // strand of the interval

       __int32 thickBegin;  // begin position for drawing thickly
       __int32 thickEnd;    // end position for drawing thickly
       BedRgb itemRgb;      // color for the item
       __int32 blockCount;  // number of blocks/exons
       String<__int32> blockSizes;   // block sizes
       String<__int32> blockBegins;  // block begin positions

       CharString data;    // any data not fitting into other members

       // Constants for marking reference id and position as invalid.
       static const __int32 INVALID_REFID = -1;
       static const __int32 INVALID_POS = -1;
   };

    }  // namespace seqan

The static members ``INVALID_POS``, ``INVALID_REFID`` store sentinel values for marking positions and reference sequence ids as invalid.

The member ``ref`` stores the contig/reference name of the genomic interval.
This information is somewhat redundant with the ``rID`` member that is filled automatically when reading from a :dox:`BedStream` such that the BedStream's ``sequenceNames[record.rID] == record.ref``.
Translating reference names to integers is useful in many applications.

When writing and ``record.rID == INVALID_REFID`` then ``record.ref`` is written out as the reference name and ``sequenceNames[record.rID]`` is written out otherwise.
The user has to take care that ``record.rID`` is a valid reference id in this case.

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

        .. includefrags:: extras/demos/tutorial/bed_io/solution2.cpp

        The output is

        .. code-block:: console

           RECORDS ON CONTIGS
           chr1    5

Assignment 3
""""""""""""

.. container:: assignment

   Generating BED From Scratch

   Type
     Application

   Objective
     Write a program that prints the following BED file.
     Create ``BedRecord<Bed6>`` objects and write them to a ``BedStream`` using ``writeRecord()``.

     .. code-block:: console

        chr7    127471196   127472363   Pos1    0   +
        chr7    127472363   127473530   Pos2    0   +


   Solution
    .. container:: foldable

       .. includefrags:: extras/demos/tutorial/bed_io/solution3.cpp

Next Steps
----------

* Continue with the :ref:`tutorial`.
