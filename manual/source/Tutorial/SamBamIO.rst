.. sidebar:: ToC

   .. contents::


.. _tutorial-sam-bam-io:

SAM and BAM I/O
===============

Learning Objective
  This tutorial will explain how to use the lower-level API for reading from and writing to SAM and BAM files.

Difficulty
  Advanced

Duration
  20 min

Prerequisites
  :ref:`tutorial-basic-sam-bam-io`, :ref:`tutorial-input-output-overview`

After reading the :ref;`tutorial-basic-sam-bam-io` tutorial, you should have a good understanding of the representation of SAM/BAM records in SeqAn library.
This tutorial will explain how the lower-level implementation of SAM and BAM I/O works and which design decisions lead to this.
You will also learn about the classes :dox:`BamIOContext` and :dox:`NameStoreCache`.
This article also explains how to write template-based code that can read both SAM files from :dox:`RecordReader RecordReaders` and BAM files from compressed streams.
Furthermore, you will learn how to use BAI indices for jumping in BAM files.

This tutorial will mostly work by showing full source code examples and explaining them.
Do not feel intimidated by its length, it is inflated by the amount of embedded source code.

The usage of SAM and BAM tags is explained in the :ref:`tutorial-basic-sam-bam-io` tutorial and will not be discussed here further.

.. important::

   Note that this tutorial is targeted at readers that already know about the SAM format.
   If you do not know about the SAM format yet then this tutorial will be harder for your to understand.
   Teaching the ins and outs of SAM is out of the scope of such a tutorial.

``BamAlignmentRecord`` Design Overview
--------------------------------------

Besides the fact that BAM is a compressed, binary format, the main difference to the plain-text format SAM is that BAM has additional header information.
The binary header of BAM files contains a plain-text SAM header with the same information as a SAM file.
However, BAM files contain an additional binary header section that stores all reference sequence names and their length.

This is equivalent to the ``@SQ`` entries in a SAM header but the ``@SQ`` entries are optional in both SAM and the SAM header in a BAM file.
Also, the ``@SQ`` entries in the SAM header of a BAM file can be out of sync with the binary reference information entries.
The authorative source for this information (most prominently the names and the order of the reference sequences).
BAM records then only contain an integer ``rID`` for each reference that referes to the ``rID``-th reference sequence in the binary header.

Because SAM and BAM files store the same information, there is only one record type to store it in SeqAn.
The record type is the class :dox:`BamAlignmentRecord` that was already introduced in the :ref:`tutorial-basic-sam-bam-io` tutorial tutorial.
The structure of this class is designed after the BAM file so conversion from the BAM file on disk to the in-memory representation is fast.
For example, the tags are stored as a binary string, the same as in a BAM file.
When reading SAM, the tags are converted into the BAM representation whereas BAM tags can be copied verbatim.
There is one important deviation, though: The qualities are stored using a phred-style ASCII encoding for consistencies with the rest of the SeqAn library.

As a reminder, here is the synopsis of the :dox:`BamAlignmentRecord` class again.

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

       // Constants for marking pos, reference id and length members invalid (== */0).
       static __int32 const INVALID_POS = -1;
       static __int32 const INVALID_REFID = -1;
       static __int32 const INVALID_LEN = 0;
   };

   }  // namespace seqan

Name Stores and Name Store Caches
---------------------------------

In order to translate from numeric reference id (``rID``) to text reference sequence name, the names have to be stored in a :dox:`StringSet` which we will call a **name store**.
For being able to translate back from a textual name (stored as a :dox:`CharString`, for example), we need a :dox:`NameStoreCache` that allows the fast lookup of numeric ids from textual names.
Both the name store and the cache are then wrapped by a :dox:`BamIOContext`.
This context object is used to prescient from the differences of SAM and BAM files when reading and writing.

For example, when writing out a :dox:`BamAlignmentRecord` to a SAM file, we need to look up the name of the reference from its numeric id to write it out as a string.
When reading a record from a SAM file, we have to translate its name string into a numeric id.
Even more, if the sequence is not know yet (remember, the ``@SQ`` headers are optional), we have to append it to the name store and register it with the cache.

Here is a minimal example of setting up a name store, name store cache, and a :dox:`BamIOContext`.
We will build upon this example below when showing how to read and write SAM and BAM files.

.. includefrags:: extras/demos/tutorial/bam_io/example1.cpp

BGZF Files / Stream
-------------------

By default, the BAM format is compressed using the BGZF compression scheme (originating from `Tabix <http://samtools.sourceforge.net/swlist.shtml>`_, but also described in the `SAM standard <http://samtools.sourceforge.net/SAM1.pdf>`_).
You can read BGZF files with tools for processing ``.gz`` files, e.g. ``gzip`` and ``zcat``.

However, there is a big difference between files written in BGZF and ``.gz`` files.
BGZF is a sequence of compressed blocks.
If the offset of a block is known, it can be decompressed independent of the rest of the file.
This information can then be used together with indices.

SeqAn provides the :dox:`BgzfStream BGZF Stream` class in the module ``<seqan/stream.h>`` to access such streams.
Here is an example for using a :dox:`BgzfStream Stream` for reading:

.. includefrags:: extras/demos/tutorial/bam_io/example2.cpp

Using a :dox:`BgzfStream BGZF Stream` for writing:

.. includefrags:: extras/demos/tutorial/bam_io/example3.cpp

Assignment 1
""""""""""""

.. container:: assignment

   Uncompressing a BGZF file.

   Type
     Review

   Objective
     Write a program that reads in a BGZF compressed file using :dox:`BgzfStream BGZF Stream` and writes the uncompressed data out again.

   Hint
     Use the function :dox:`StreamConcept#streamReadBlock` and :dox:`StreamConcept#streamWriteBlock` for reading and writing data into and from buffers.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/bam_io/solution1.cpp

Reading and Writing Headers
---------------------------

The data structure :dox:`BamHeader` has already been described in the :ref:`tutorial-basic-sam-bam-io` so we will not repeat that here.
Instead, we will focus on how to read headers from SAM and BAM files.

Here is a minimal example of reading and writing a header from and to a SAM file.
The example contains the creation of a :dox:`BamIOContext`, the necessary :dox:`RecordReader` and full error handling.

.. includefrags:: extras/demos/tutorial/bam_io/example4.cpp

Reading and writing headers from and to BAM files is simple.
We simply replace ``seqan::Sam()`` by ``seqan::Bam()`` and use :dox:`BgzfStream BGZF Stream` objects instead of uncompressed streams.
Also, we do not need a :dox:`RecordReader` any more.

.. includefrags:: extras/demos/tutorial/bam_io/example5.cpp

Note that except for the types, the signatures of the functions ``readRecord()`` and ``write()`` are the same.
Thus, we can make copying of the header a template function ``copyHeader()``.
This function can now be used for both BAM and SAM.

.. includefrags:: extras/demos/tutorial/bam_io/example6.cpp

Assignment 2
""""""""""""

.. container:: assignment

   Converting BAM header to SAM.

   Type
     Application

   Objective
     Write a program that reads the header from a BAM file and writes it out as a SAM header to ``std::cout``.

   Solution
     .. container:: foldable

         .. includefrags:: extras/demos/tutorial/bam_io/solution2.cpp

Reading and Writing Alignment Records
-------------------------------------

:dox:`BamAlignmentRecord BamAlignmentRecords` can be read and written the same way as :dox:`BamHeader` objects.
Here is an example for reading and writing of alignment records from SAM and to files.

.. code-block:: cpp

   // Copy over records.
   seqan::BamAlignmentRecord record;
   while (atEnd(reader))
   {
       if (readRecord(record, context, reader, seqan::Sam()) != 0)
       {
           std::cerr << "ERROR: Could not read record from SAM file " << argv[1] << "\n";
           return 1;
       }

       if (write2(outStream, record, context, seqan::Sam()) != 0)
       {
           std::cerr << "ERROR: Could not write record to SAM file " << argv[2] << "\n";
           return 1;
       }
   }

And here is the modified version for the BAM format.
The only changes are that

* we do not read from a :dox:`RecordReader` but a :dox:`BgzfStream BGZF Stream` instead,
* we need to write to a :dox:`BgzfStream BGZF Stream`, and
* we need to use the tag ``seqan::Bam()`` instead of ``seqan::Sam()``.

.. code-block:: cpp

   // Copy over records.
   seqan::BamAlignmentRecord record;
   while (atEnd(reader))
   {
       if (readRecord(record, context, inStream, seqan::Bam()) != 0)
       {
           std::cerr << "ERROR: Could not read record from BAM file " << argv[1] << "\n";
           return 1;
       }

       if (write2(outStream, record, context, seqan::Bam()) != 0)
       {
           std::cerr << "ERROR: Could not write record to BAM file " << argv[2] << "\n";
           return 1;
       }
    }

Assignment 3
""""""""""""

.. container:: assignment

   Converting whole BAM files to SAM.

   Type
      Application

   Objective
      Modify the solution of Assignment 2 to not only convert the header to BAM but also the alignment records.

   Solution
      .. container:: foldable

         .. includefrags:: extras/demos/tutorial/bam_io/solution3.cpp

Using Indices
-------------

SeqAn also contains support for reading BAM indices with the format ``.bai``.
These indices can be built using the ``samtools index`` command.

You can read such indices into a :dox:`BaiBamIndex BAI BamIndex` object with the function :dox:`BamIndex#read`.
Then, you can use the function seqan:"Function.BamIndex#jumpToRegion" to jump within BAM files.

After jumping, the next record that is read is before at the given position.
This means, you have to manually read as many records up until the one you are looking for is found.
The reason for this is that the function :dox:`BamIndex#jumpToRegion` would have to read until it finds the first record that is right from or at the given position.
This would lead to this record being lost.

.. includefrags:: extras/demos/tutorial/bam_io/example7.cpp

Next Steps
----------

* Read the `SAM Specification (pdf) <http://samtools.sourceforge.net/SAM1.pdf>`_.
* Continue with the :ref:`tutorial`.
