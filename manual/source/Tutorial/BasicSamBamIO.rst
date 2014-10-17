.. sidebar:: ToC

   .. contents::


.. _tutorial-basic-sam-bam-io:

Basic SAM and BAM I/O
=====================

Learning Objective
  In this tutorial, you will learn how to use the high-level interface :dox:`BamStream` class to read and write SAM and BAM files.

Difficulty
  Average

Duration
  1 h (45 min if you know the SAM format)

Prerequisites
  :ref:`tutorial-sequences`, :ref:`tutorial-basic-sequence-io`, Exposure to the SAM format

This tutorial deals with how to easily read and write SAM and BAM files using the :dox:`BamStream` class.
It starts out with a quick reminder on the structure of SAM (and also BAM) files and will then continue with how to read and write SAM/BAM files and access the tags of a record.

.. important::

    Note that this tutorial is targeted at readers that already know about the SAM format.
    If you do not know about the SAM format yet then this tutorial will be harder for your to understand.
    Teaching the ins and outs of SAM is out of the scope of such a tutorial.

Both SAM and BAM file store multi-read alignments.
Storing alignments of longer sequences such as contigs from assemblies is also possible, but less common.
Here, we will focus on multi-read alignments.

SAM files are text files, having one record per line.
BAM files are just binary, compressed versions of SAM files that have a stricter organization and aim to be more efficiently useable by programs and computers.
The nuts and bolts of the formats are described in the `SAM Format Specification <http://samtools.sourceforge.net/SAM1.pdf>`_.

The SAM and BAM related I/O functionality in SeqAn focuses on allowing access to these formats in SeqAn with thin abstractions.
The :ref:`tutorial-fragment-store` Tutorial shows how to get a more high-level abstraction for multi-read alignments.

.. important::

    SAM/BAM I/O vs. Fragment Store

    The :ref:`tutorial-fragment-store` provides a high-level view of multi-read alignments.
    This is very useful if you want to do SNP or small indel detection because you need to access the alignment of the reads around your candidate regions.
    However, storing the whole alignment of a 120GB BAM file obviously is not a good idea.

    The SAM/BAM I/O functionaliy in SeqAn is meant for sequentially reading through SAM and BAM files.
    Jumping within BAM files using BAI indices is described in the :ref:`tutorial-sam-bam-io` tutorial.

SAM / BAM File Structure
------------------------

We will give an quick overview of the SAM and BAM formats here.
Note that this overview serves more to remind you what the formats are about and are not meant to teach how to use the SAM and BAM format.

The following shows an example of a SAM file.

::

    @HD VN:1.3  SO:coordinate
    @SQ SN:ref  LN:45
    @SQ SN:ref2 LN:40
    r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
    r002    0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0   AAAAGATAAGGGATAAA   *
    r003    0   ref 9   30  5H6M    *   0   0   AGCTAA  *
    r004    0   ref 16  30  6M14N1I5M   *   0   0   ATAGCTCTCAGC    *
    r003    16  ref 29  30  6H5M    *   0   0   TAGGC   *
    r001    83  ref 37  30  9M  =   7   -39 CAGCGCCAT   *

SAM files are TSV (tab-separated-values) files and begin with an optional header.
The header consists of multiple lines, starting with an ``'@'`` character, each line is a record.
Each record starts with its identifier and is followed by tab-separated tags.
Each tag in the header consists of a two-character identifier, followed by ``':'``, followed by the value.

If present, the ``@HD`` record must be the first record and specifies the SAM version (tag ``VN``) used in this file and the sort order (``SO``).
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
* The qualities must be stored as ASCII PRED-encoded qualities.
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

The following example shows an example of a program that reads the file with the path ``example.sam`` and prints its contents back to the user on stdout.
If you want to try out this program then create a file with the sample SAM content from above and adjust the path ``"example.sam"`` in the program below to the path to your SAM file (e.g. ``"path/to/my_example.sam"``).

.. includefrags:: extras/demos/tutorial/basic_sam_bam_io/example1.cpp

The program first opens a :dox:`BamStream` for reading, then one for writing.
Note that :dox:`BamStream` automatically guesses the file type from the file contents when reading and from the file name when writing.
You can also force a format using :dox:`BamStream::BamStream BamStream's constructor`.
You can read from stdin and write to stdout using ``"-"`` as the file name.

The header is automatically read when a :dox:`BamStream` is opened.
After the header has been read, it is copied over into the output stream.
Then, the input stream is read record by record and written out to the output stream.
Note that the header is written out automatically before the first alignment record is written.

The alignment records are read into :dox:`BamAlignmentRecord` objects which we will focus on below.

Note that the example above is missing error handling.
This means that if the input format is ill-formed, error return codes are not handled appropriately and the program might do something unexpected in the case of an error.

For example, if the file contains trailing empty lines, the program will loop indefinitely as can be seen in the shell output below:

.. code-block:: console

   # tutorial_basic_sam_bam_io_example1
   @HD     VN:1.3  SO:coordinate
   @SQ     SN:ref  LN:45
   @SQ     SN:ref2 LN:40
   r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
   r002    0       ref     9       30      1S2I6M1P1I1P1I4M2I      *       0       0       AAAAGATAAGGGATAAA       *
   r003    0       ref     9       30      5H6M    *       0       0       AGCTAA  *
   r004    0       ref     16      30      6M14N1I5M       *       0       0       ATAGCTCTCAGC    *
   r003    16      ref     29      30      6H5M    *       0       0       TAGGC   *
   r001    83      ref     37      30      9M      =       7       -39     CAGCGCCAT       *
	   83      *       *       *       *       *       0       *       *       *
	   83      *       *       *       *       *       0       *       *       *
   ...

We can fix this problem by introducing error handling.
The :dox:`BamStream#readRecord` call returns a status code different from ``0``, indicating an error because an empty line does not form a valid SAM record line.
However, it stops processing as soon as an errernous record is detected which makes the call to :dox:`BamStream#atEnd` return ``false`` and run in an infinite loop

In Assignment 1, we will add error handling to the program.

Assignment 1
""""""""""""

.. container:: assignment

   Adding Error Handling

   Type
     Review

   Objective
     Add error handling using the hints below.

   Hints
      The functions :dox:`BamStream#readRecord` and :dox:`BamStream#writeRecord` return a status code ``int``, ``0`` on success, ``1`` on errors.
      The function :dox:`BamStream#isGood` checks whether the state of a :dox:`BamStream` is errorneous.

   Solution
      .. container:: foldable

         .. includefrags:: extras/demos/tutorial/basic_sam_bam_io/solution1.cpp

The Class :dox:`BamAlignmentRecord`
-----------------------------------

The class :dox:`BamAlignmentRecord` stores one alignment record in a SAM or BAM file.
The class gives a in-memory representation that (1) is independent of whether it comes from/goes to a SAM or BAM file, (2) at the same time follows both formats closely, (3) allows for efficient storage and usage in C++, and (4) integrates well with the rest of the SeqAn library.

The following definition gives an overview that annotate which fields are available, the field types, and how they map to the SAM and BAM fields.
Note that we use the :dox:`CigarElement` class to store entries in the CIGAR string.

.. code-block:: cpp

   namespace seqan {

   class BamAlignmentRecord
   {
   public:
       CharString qName;               // QNAME
       __uint16 flag;                  // FLAG
       __int32 rID;                    // REF
       __int32 beginPos;               // POS
       __uint8 mapQ;                   // MAPQ mapping quality, 255 for */invalid
       __uint16 bin;                   // bin for indexing
       String<CigarElement<> > cigar;  // CIGAR string
       __int32 rNextId;                // RNEXT (0-based)
       __int32 pNext;                  // PNEXT (0-based)
       __int32 tLen;                   // TLEN
       CharString seq;                 // SEQ, as in SAM/BAM file.
       CharString qual;                // Quality string as in SAM (Phred).
       CharString tags;                // Tags, raw as in BAM.

       // Constants for marking pos, reference id and length members invalid (== 0/*).
       static __int32 const INVALID_POS = -1;
       static __int32 const INVALID_REFID = -1;
       static __int32 const INVALID_LEN = 0;
   };

   }  // namespace seqan

The static members ``INVALID_POS``, ``INVALID_REFID``, and ``INVALID_LEN`` store sentinel values for marking positions, reference sequence ids, and lengths as invalid or N/A.

An important related type is the enum :dox:`BamFlags` that provides constants for bit operations on the ``flag`` field.
The functions :dox:`BamAlignmentRecord#hasFlagAllProper`, :dox:`BamAlignmentRecord#hasFlagDuplicate`, :dox:`BamAlignmentRecord#hasFlagFirst`, :dox:`BamAlignmentRecord#hasFlagLast`, :dox:`BamAlignmentRecord#hasFlagMultiple`, :dox:`BamAlignmentRecord#hasFlagNextRC`, :dox:`BamAlignmentRecord#hasFlagNextUnmapped`, :dox:`BamAlignmentRecord#hasFlagQCNoPass`, :dox:`BamAlignmentRecord#hasFlagRC`, :dox:`BamAlignmentRecord#hasFlagSecondary`, :dox:`BamAlignmentRecord#hasFlagUnmapped`, and :dox:`BamAlignmentRecord#hasFlagSupplementary` allow for easy reading of flags.

For example, the following loop sums up the length of the sequences that did not align:

.. code-block:: cpp

   seqan::BamAlignmentRecord record;
   unsigned lenSum = 0;
   while (atEnd(bamStreamIn))
       if (hasFlagUnmapped(record))
	   lenSum += length(record.seq);

.. container:: assignment

   Counting Records

   Type
     Review

   Objective
     Extend the result of Assignment 1 by counting the number of unmapped reads.

   Hints
     Use the function :dox:`BamAlignmentRecord#hasFlagUnmapped`.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/basic_sam_bam_io/solution2.cpp

The Classes :dox:`BamHeader` and :dox:`BamHeaderRecord`
-------------------------------------------------------

The header information is stored in the class :dox:`BamHeader`.
This class gives a unified in-memory representation for SAM and BAM files.

The class has two members: ``records`` and ``sequenceInfos``.
We will focus on ``sequenceInfos`` here.
``sequenceInfos`` is a :dox:`String` of :dox:`Pair` objects.
The first entry of the pair is a :dox:`CharString` with the sequence name and the second entry is a ``_int32`` with the sequence length.
Note that the ``@SQ`` header lines in the header and the ``sequenceInfos`` fields are not kept in sync automatically.

The following example program prints the sequences and lengths from a BAM file.

.. includefrags:: extras/demos/tutorial/basic_sam_bam_io/example2.cpp

Note that this is only guaranteed to work for BAM files because this information is not mandatory in SAM files and might be missing.
When writing files, you have to fill the ``sequenceInfos`` string appropriately before writing any record.

.. tip::

    Building Ref-ID Mappings Using ``sequenceInfos``.

    The following example gives a typical example for using the ``sequenceInfos`` member:
    You want to post-process a BAM file together with the reference FASTA file.
    The sequences in the FASTA file are the same but their order may have changed.
    For example, because the FASTA file from the mapping step has been generated from the chromosomes by concatenation in a different order than the currently present one.

    .. includefrags:: extras/demos/tutorial/basic_sam_bam_io/example3.cpp

Assignment 3
""""""""""""

.. container:: assignment

   Generating SAM From Scratch

   Type
     Application

   Objective
     Write a program that prints a SAM file, including headers ``@HD`` and ``@SQ``.
     The content should be all 12-mers of the reference sequence ``"CCCGATGAGCACACGATCACACGATGACA"``, called ``"REF"``.
     The name should be ``"REF_${begin pos}_${end pos}"``.
     You only have to fill the members ``qId``, ``rID``, ``beginPos``, ``cigar``, and ``flag`` (set ``flag`` to ``0``).

   Hints
     You can convert integers into strings using the ``<sstream>`` STL header.

     .. code-block: cpp

	#include <sstream>
	// ...
	std::stringstream ss;
	ss << 10;
	seqan::CharString str = ss.str();  // => == "10"
	// To reset ss, we need two calls:
	ss.str("");  // Remove contents.
	ss.clear();  // Reset any error bits.

     The first lines of the result should read as follows:

     ::

	 @HD VN:1.4
	 @SQ SN:REF  LN:29
	 REF_0_12    0   REF 1   *   12= *   0   *   CCCGATGAGCAC    *
	 REF_1_13    0   REF 2   *   12= *   0   *   CCGATGAGCACA    *
	 REF_2_14    0   REF 3   *   12= *   0   *   CGATGAGCACAC    *
	 REF_3_15    0   REF 4   *   12= *   0   *   GATGAGCACACG    *


   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/basic_sam_bam_io/solution3.cpp


Accessing the Tags
------------------

As seen above, accessing the header tags is simple since it is a string of tag/value pairs.
The whole header is completely read, parsed, and converted into this structure when the file is opened.
The header is expected to be small, especially when compared to the rest of the file, and thus the time and memory spent is neglectable.

The alignment record tags are a different story.
The tags only contain auxiliary information that are not of interest for all use cases.
Always parsing the tags would not be in agreement with C++'s and SeqAn's device "you only pay for what you use", especially for BAM files that are expected to contain millions of records.
Also, the tags of the alignment records are typed, e.g. ``NM:i:10`` is an integer tag named ``"NM"`` with the value ``10``.

Thus, the following strategy is used.
Alignment record tags from BAM files are copied byte-wise into the ``tag`` member of :dox:`BamAlignmentRecord` in a verbatim fashion.
When reading from SAM, the tags are converted into format used by BAM tags.

Then, you can use the :dox:`BamTagsDict` class to access the the tag list of a record in a dictionary-like fashion.
This class also performs the necessary casting when reading and writing tag list entries.

:dox:`BamTagsDict` acts as a wrapper around the ``tags`` member (which is of type :dox:`CharString`) of a :dox:`BamAlignmentRecord`:

.. code-block:: cpp

   seqan::BamAlignmentRecord record;
   seqan::BamTagsDict tagsDict(record.tags);

We can add a tag using the function :dox:`BamTagsDict#setTagValue`.
When setting an already existing tag's value, its value will be overwritten.
Note that in the following, we give the tags value in SAM format because it is easier to read, although they are stored in BAM format internally.

.. code-block:: cpp

   setTagValue(tagsDict, "NM", 2);
   // => tags: "NM:i:2"
   setTagValue(tagsDict, "NH", 1);
   // => tags: "NM:i:2 NH:i:1"
   setTagValue(tagsDict, "NM", 3);
   // => tags: "NM:i:3 NH:i:1"

The first parameter to :dox:`BamTagsDict#setTagValue` is the :dox:`BamTagsDict`, the second one is a two-character string with the key, and the third one is the value.
Note that the type of tag entry will be taken automatically from the type of the third parameter.

Reading values is slightly more complex because we have to handle the case that the value is not present.
First, we get the index of the tag in the tag list.

.. code-block:: cpp

   unsigned myIdx = 0;
   bool keyFound = findTagKey(myIdx, tagsDict, "NH");
   if (keyFound)
       std::cerr << "ERROR: Unknown key!\n";

Then, we can read the value from the :dox:`BamTagsDict` using the function :dox:`BamTagsDict#extractTagValue`.

.. code-block:: cpp

   int valInt = 0;
   bool ok = extractTagValue(valInt, tagsDict, myIdx);
   if (ok)
       std::cerr << "ERROR: There was an error extracting NH from tags!\n";

The function returns a ``bool`` that is ``true`` on success and ``false`` otherwise.
The extraction can fail if the index is out of bounds or the value in the dictionary cannot be cast to the type of the first parameter.

The value in the tags dictionary will be casted to the type of the first parameter (result parameter) of :dox:`BamTagsDict#extractTagValue`:

.. code-block:: cpp

   short valShort = 0;
   extractTagValue(valShort, tagsDict, myIdx);

Assignment 4
""""""""""""

.. container:: assignment

   Writing Tags

   Type
     Review

   Objective
     Modify the solution of Assignment 3 to also write the ``"NH"`` tag.
     This tag stores an ``int`` value that is the number of records for this query.
     In our case, the value is always ``1``.

     The first lines of the result should read as follows:

     ::

         @HD VN:1.4
         @SQ SN:REF  LN:29
         REF_0_12    0   REF 1   *   12= *   0   *   CCCGATGAGCAC    *   NH:i:1
         REF_1_13    0   REF 2   *   12= *   0   *   CCGATGAGCACA    *   NH:i:1
         REF_2_14    0   REF 3   *   12= *   0   *   CGATGAGCACAC    *   NH:i:1
         REF_3_15    0   REF 4   *   12= *   0   *   GATGAGCACACG    *   NH:i:1


   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/basic_sam_bam_io/solution4.cpp

Congratulations, you have now learned to read and write SAM and BAM files.

Next Steps
~~~~~~~~~~

* Read the `SAM Specification (pdf) <http://samtools.sourceforge.net/SAM1.pdf>`_.
* Continue with the :ref:`tutorial`.
