.. sidebar:: ToC

    .. contents::

.. _tutorial-io-sam-bam-io:

SAM and BAM I/O
===============

Learning Objective
  In this tutorial, you will learn how to read and write SAM and BAM files.

Difficulty
  Average

Duration
  1 h (45 min if you know the SAM format)

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-io-input-output-overview`, `SAM Format Specification <http://samtools.sourceforge.net/SAM1.pdf>`_

Overview
--------

.. warning::

    Before you can read/write BAM files (bgzf compressed SAM files) you need to make sure that your program is linked against the zlib library.
    When you build your application within the SeqAn build infrastructure, the zlib library is automatically located by looking at the standard places for the library.
    Also have a look at :ref:`tutorial-io-input-output-overview-formatted-files` to read more about support of compressed file I/O.
    If the macro ``SEQAN_HAS_ZLIB`` is set to ``0`` then reading/writing BAM file format is disabled.
    It is set to ``1`` if the zlib could be found and reading/writing of compressed files is enabled automatically.
    You can read :ref:`infra-use-cmake`, :ref:`infra-use-custom` and :ref:`infra-use-install-dependencies` for further notes about using the zlib and libbz2 in your build infrastructure.

This tutorial shows how to read and write SAM and BAM files using the :dox:`BamFileIn` and :dox:`BamFileOut` classes.
It starts out with a quick reminder on the structure of SAM (and also BAM) files and continues with how to read and write SAM/BAM files and access the tags of a record.

.. important::

    Note that this tutorial is targeted at readers that already know about the SAM format.
    If you do not know about the SAM format yet, then this tutorial will be harder for your to understand.

Both SAM and BAM files store multi-read alignments.
Storing alignments of longer sequences such as contigs from assemblies is also possible, but less common.
Here, we will focus on multi-read alignments.

SAM files are text files, having one record per line.
BAM files are just binary, compressed versions of SAM files that have a stricter organization and aim to be more efficiently usable by programs and computers.
The nuts and bolts of the formats are described in the `SAM Format Specification <http://samtools.sourceforge.net/SAM1.pdf>`_.

The SAM and BAM related I/O functionality in SeqAn focuses on allowing access to these formats in SeqAn with thin abstractions.
The :ref:`tutorial-datastructures-store-fragment-store` Tutorial shows how to get a more high-level abstraction for multi-read alignments.

.. important::

    SAM/BAM I/O vs. Fragment Store

    The :ref:`tutorial-datastructures-store-fragment-store` provides a high-level view of multi-read alignments.
    This is very useful if you want to do SNP or small indel detection because you need to access the alignment of the reads around your candidate regions.
    However, storing the whole alignment of a 120GB BAM file obviously is not a good idea.

    The SAM/BAM I/O functionality in SeqAn is meant for sequentially reading through SAM and BAM files.
    Jumping within BAM files using BAI indices is described in the `Using BAM Indices`_ section of this tutorial.


SAM / BAM Format
----------------

The following shows an example of a SAM file.

.. includefrags:: demos/tutorial/sam_and_bam_io/example.sam

SAM files are TSV (tab-separated-values) files and begin with an optional header.
The header consists of multiple lines, starting with an ``'@'`` character, each line is a record.
Each record starts with its identifier and is followed by tab-separated tags.
Each tag in the header consists of a two-character identifier, followed by ``':'``, followed by the value.

If present, the ``@HD`` record must be the first record which specifies the SAM version (tag ``VN``) used in this file and the sort order (``SO``).
The optional ``@SQ`` header records give the reference sequence names (tag ``SN``) and lengths (tag ``LN``).
There also are other header record types.

The optional header section is followed by the alignment records.
The alignment records are again tab-separated.
There are 11 mandatory columns.

+-----------+-------------+--------------+-----------------+-------------------------------------------+
| Col       | Field       | Type         | N/A Value       | Description                               |
+===========+=============+==============+=================+===========================================+
| 1         | QNAME       | string       | mandatory       | The query/read name.                      |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 2         | FLAG        | int          | mandatory       | The record's flag.                        |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 3         | RNAME       | string       | ``*``           | The reference name.                       |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 4         | POS         | 32-bit int   | ``0``           | 1-based position on the reference.        |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 5         | MAPQ        | 8-bit int    | ``255``         | The mapping quality.                      |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 6         | CIGAR       | string       | ``*``           | The CIGAR string of the alignment.        |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 7         | RNEXT       | string       | ``*``           | The reference of the next mate/segment.   |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 8         | PNEXT       | string       | ``0``           | The position of the next mate/seqgment.   |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 9         | TLEN        | string       | ``0``           | The observed length of the template.      |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 10        | SEQ         | string       | ``*``           | The query/read sequence.                  |
+-----------+-------------+--------------+-----------------+-------------------------------------------+
| 11        | QUAL        | string       | ``*``           | The ASCII PHRED-encoded base qualities.   |
+-----------+-------------+--------------+-----------------+-------------------------------------------+

Notes:

* The SAM standard talks about "queries".
  In the context of read mapping, where the format originates, queries are reads.
* The SAM standard talks about "templates" and "segments".
  In the case of paired-end and mate-pair mapping the template consists of two segments, each is one read.
  The template length is the insert size.
* Paired-end reads are stored as two alignments records with the same QNAME.
  The first and second mate are discriminated by the FLAG values.
* When the FLAG indicates that SEQ is reverse-complemented, then QUAL is reversed.
* Positions in the SAM file are 1-based.
  When read into a :dox:`BamAlignmentRecord` (see below), the positions become 0-based.
* The qualities must be stored as ASCII PHRED-encoded qualities.
* The query and reference names must not contain whitespace.
  It is common to trim query and reference ids at the first space.

There are many ambiguities, recommendations, and some special cases in the formats that we do not describe here.
We recommend that you follow this tutorial, start working with the SAM and BAM formats and later read the SAM specification "on demand" when you need it.

The 11 mandatory columns are followed by an arbitrary number of optional tags.
Tags have a two-character identifier followed by ``":${TYPE}:"``, followed by the tag's value.

BAM files store their header as plain-text SAM headers.
However, they additionally store the name and length information about the reference sequences.
This information is mandatory since in BAM, the alignment records only contain the numeric ids of the reference sequences.
Thus, the name is stored outside the record in the header.

A First Working Example
-----------------------

The following program reads a file named ``example.sam`` and prints its contents back to the user on standard output.

.. includefrags:: demos/tutorial/sam_and_bam_io/solution1.cpp

.. includefrags:: demos/tutorial/sam_and_bam_io/solution1.cpp.stdout

We instantiate a :dox:`BamFileIn` object for reading and a :dox:`BamFileOut` object for writing.
First, we read the BAM header with :dox:`FormattedFileIn#readRecord` and we write it with :dox:`FormattedFileOut#writeRecord`.
Then, we read each record from the input file and print it back on standard output.
The alignment records are read into :dox:`BamAlignmentRecord` objects, which we will focus on below.

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Reproduction

   Objective
     Create a file with the sample SAM content from above and adjust the path ``"example.sam"`` to the path to your SAM file (e.g. ``"/path/to/my_example.sam"``).

   Solution
      .. container:: foldable

         .. includefrags:: demos/tutorial/sam_and_bam_io/solution1.cpp


Accessing the Header
--------------------

Sequence information (i.e. @SQ records) from the BAM header is stored in the :dox:`BamIOContext`.
All remaining BAM header information is stored in the class :dox:`BamHeader`.

.. important::
   The header is not mandatory in SAM files and might be missing.

The following program accesses the :dox:`BamIOContext` of its :dox:`BamFileIn` and prints the reference sequence names and lengths present in the BAM header.

.. includefrags:: demos/tutorial/sam_and_bam_io/example2.cpp

The output looks like this:

.. includefrags:: demos/tutorial/sam_and_bam_io/example2.cpp.stdout

Accessing the Records
---------------------

The class :dox:`BamAlignmentRecord` stores one alignment record of a SAM or BAM file.
The class gives an in-memory representation that

#. is independent of whether it comes from/goes to a SAM or BAM file,
#. at the same time follows both formats closely,
#. allows for efficient storage and usage in C++ and
#. integrates well with the rest of the SeqAn library.

The following definition gives an overview of the available fields, their types, and how they map to the SAM and BAM fields.
Note that we use the :dox:`CigarElement` class to store entries in the CIGAR string.

.. includefrags:: demos/tutorial/sam_and_bam_io/base.cpp
      :fragment: bamRecord

The static members ``INVALID_POS``, ``INVALID_REFID``, and ``INVALID_LEN`` store sentinel values for marking positions, reference sequence ids, and lengths as invalid or N/A.

.. tip::
   A :dox:`BamAlignmentRecord` is linked to a reference sequence by the field ``rID``.
   The reference sequence information is stored in the BAM header and kept in the :dox:`BamIOContext`.
   To easily access reference sequence name and and length relative to a given :dox:`BamAlignmentRecord` within a :dox:`BamFileIn`, use functions :dox:`BamAlignmentRecord#getContigName` and :dox:`BamAlignmentRecord#getContigLength`.

An important related type is the enum :dox:`BamFlags` that provides constants for bit operations on the ``flag`` field.
The functions :dox:`BamAlignmentRecord#hasFlagAllProper`, :dox:`BamAlignmentRecord#hasFlagDuplicate`, :dox:`BamAlignmentRecord#hasFlagFirst`, :dox:`BamAlignmentRecord#hasFlagLast`, :dox:`BamAlignmentRecord#hasFlagMultiple`, :dox:`BamAlignmentRecord#hasFlagNextRC`, :dox:`BamAlignmentRecord#hasFlagNextUnmapped`, :dox:`BamAlignmentRecord#hasFlagQCNoPass`, :dox:`BamAlignmentRecord#hasFlagRC`, :dox:`BamAlignmentRecord#hasFlagSecondary`, :dox:`BamAlignmentRecord#hasFlagUnmapped`, and :dox:`BamAlignmentRecord#hasFlagSupplementary` allow for easy reading of flags.



Assignment 2
""""""""""""

.. container:: assignment

   Counting Records

   Type
     Review

   Objective
     Count the number of unmapped reads.

   Hints
     Use the function :dox:`BamAlignmentRecord#hasFlagUnmapped`.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/sam_and_bam_io/solution2.cpp

        .. includefrags:: demos/tutorial/sam_and_bam_io/solution2.cpp.stdout


Accessing the Records' Tags
---------------------------

You can use the :dox:`BamTagsDict` class to access the the tag list of a record in a dictionary-like fashion.
This class also performs the necessary casting when reading and writing tag list entries.

:dox:`BamTagsDict` acts as a wrapper around the raw ``tags`` member of a :dox:`BamAlignmentRecord`, which is of type :dox:`CharString`:

.. includefrags:: demos/tutorial/sam_and_bam_io/base.cpp
      :fragment: BamTagsDict

We can add a tag using the function :dox:`BamTagsDict#setTagValue`.
When setting an already existing tag's value, its value will be overwritten.
Note that in the following, we give the tags value in SAM format because it is easier to read, although they are stored in BAM format internally.

.. includefrags:: demos/tutorial/sam_and_bam_io/base.cpp
      :fragment: addTag

The first parameter to :dox:`BamTagsDict#setTagValue` is the :dox:`BamTagsDict`, the second one is a two-character string with the key, and the third one is the value.
Note that the type of tag entry will be taken automatically from the type of the third parameter.

Reading values is slightly more complex because we have to handle the case that the value is not present.
First, we get the index of the tag in the tag list.

.. includefrags:: demos/tutorial/sam_and_bam_io/base.cpp
      :fragment: getIndex

Then, we can read the value from the :dox:`BamTagsDict` using the function :dox:`BamTagsDict#extractTagValue`.

.. includefrags:: demos/tutorial/sam_and_bam_io/base.cpp
      :fragment: extractValue

The function returns a ``bool`` that is ``true`` on success and ``false`` otherwise.
The extraction can fail if the index is out of bounds or the value in the dictionary cannot be cast to the type of the first parameter.

The value in the tags dictionary will be casted to the type of the first parameter of :dox:`BamTagsDict#extractTagValue`:

.. includefrags:: demos/tutorial/sam_and_bam_io/base.cpp
      :fragment: cast

Assignment 3
""""""""""""

.. container:: assignment

   Reading Tags

   Type
     Review

   Objective
     Modify the solution of Assignment 2 to count the number of records having the ``"XX"`` tag.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/sam_and_bam_io/solution3.cpp

        .. includefrags:: demos/tutorial/sam_and_bam_io/solution3.cpp.stdout

Using BAM Indices
-----------------

SeqAn also contains features for reading BAM indices with the format ``.bai``. These indices can be built using the ``samtools index`` command. In the near future we plan to support building the bam index with SeqAn as well.

You can read indices into a :dox:`BaiBamIndex` object with the function :dox:`BamIndex#open`. Then, you can use the function :dox:`BamFileIn#jumpToRegion` to jump to a specific position within BAM files. After jumping, the next record to be read is before the given region. Therefore, you have to skip records until you access the one you are looking for.

.. includefrags:: demos/tutorial/sam_and_bam_io/example7.cpp

.. includefrags:: demos/tutorial/sam_and_bam_io/example7.cpp.stdout

Next Steps
----------

* Read the `SAM Format Specification <http://samtools.sourceforge.net/SAM1.pdf>`_.
* Continue with the :ref:`tutorial`.
