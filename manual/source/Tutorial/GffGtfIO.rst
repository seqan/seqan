.. sidebar:: ToC

   .. contents::


.. _tutorial-gff-and-gtf-io:

GFF and GTF I/O
===============

Learning Objective
  In this tutorial, you will learn how to use the high-level interface :dox:`GffStream` class to read and write GFF and GTF files.

Difficulty
  Average

Duration
 45 min

Prerequisites
  Exposure to the GFF and GTF formats is useful.

This tutorial deals with how to easily read and write GFF and GTF files using the :dox:`GffStream` class.
It starts out with a quick reminder on the structure of GFF and GTF files and will then continue with how to read and write GFF and GTF files.

The GFF and GTF formats are used for annotating genomic intervals (an interval with begin/end position on a contig/chromosome).
GFF exist in versions 2 and 3 and GTF is sometimes called "GFF 2.5".
There are specifications for `GFF 2 <http://www.sanger.ac.uk/resources/software/gff/spec.html>`_, `GFF 3 <http://www.sequenceontology.org/gff3.shtml>`_, and `GTF <http://mblab.wustl.edu/GTF22.html>`_ available elsewhere.
GFF and GTF are TSV-based formats and in general have the same structure.
The main difference is the underlying system/ontology for the annotation but also smaller differences in the format.

In this tutorial, we will focus on the format GFF 3 since it is the most current one with most complete tool support.
The information of this tutorial can easily be translated to the other two formats.

The SeqAn module ``gff_io`` allows the reading and writing of the GFF and GTF formats.

.. tip::

    Format Version Support in SeqAn

    :dox:`GffStream` allows to read GFF files in version 2 and 3 and GTF files.
    For writing, only GFF 3 and GTF are supported at the moment.

GFF File Structure
------------------

The following is an example of a GFF 3 file:

::

    ctg123  .   gene    1000    9000    .   +   .   ID=gene00001;Name=EDEN
    ctg123  .   TF_binding_site 1000    1012    .   +   .   Parent=gene00001
    ctg123  .   mRNA    1050    9000    .   +   .   ID=mRNA00001;Parent=gene00001
    ctg123  .   mRNA    1050    9000    .   +   .   ID=mRNA00002;Parent=gene00001
    ctg123  .   mRNA    1300    9000    .   +   .   ID=mRNA00003;Parent=gene00001
    ctg123  .   exon    1300    1500    .   +   .   Parent=mRNA00003
    ctg123  .   exon    1050    1500    .   +   .   Parent=mRNA00001,mRNA00002
    ctg123  .   exon    3000    3902    .   +   .   Parent=mRNA00001,mRNA00003
    ctg123  .   exon    5000    5500    .   +   .   Parent=mRNA00001,mRNA00002,mRNA00003
    ctg123  .   exon    7000    9000    .   +   .   Parent=mRNA00001,mRNA00002,mRNA00003
    ctg123  .   CDS 1201    1500    .   +   0   ID=cds00001;Parent=mRNA00001
    ctg123  .   CDS 3000    3902    .   +   0   ID=cds00001;Parent=mRNA00001
    ctg123  .   CDS 5000    5500    .   +   0   ID=cds00001;Parent=mRNA00001
    ctg123  .   CDS 7000    7600    .   +   0   ID=cds00001;Parent=mRNA00001
    ctg123  .   CDS 1201    1500    .   +   0   ID=cds00002;Parent=mRNA00002
    ctg123  .   CDS 5000    5500    .   +   0   ID=cds00002;Parent=mRNA00002
    ctg123  .   CDS         7000    7600    .   +   0   ID=cds00002;Parent=mRNA00002
    ctg123  .   CDS 3301    3902    .   +   0   ID=cds00003;Parent=mRNA00003
    ctg123  .   CDS         5000    5500    .   +   1   ID=cds00003;Parent=mRNA00003
    ctg123  .   CDS         7000    7600    .   +   1   ID=cds00003;Parent=mRNA00003
    ctg123  .   CDS 3391    3902    .   +   0   ID=cds00004;Parent=mRNA00003
    ctg123  .   CDS         5000    5500    .   +   1   ID=cds00004;Parent=mRNA00003
    ctg123  .   CDS         7000    7600    .   +   1   ID=cds00004;Parent=mRNA00003

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
      For example, the interval ``[1,000; 2,000]`` starts at character 1,000 and ends at character 2,000 and includes it.
      This way is natural to non-programmers and used when giving coordinates in GFF files or genome browsers such as UCSC Genome Browser and IGV.
   #. Start counting positions at 0 and give intervals by the first position that is part of the interval and giving the position behind the last position that is part of the interval.
      The interval from above would be ``[999; 2,000)`` in this case.

   In text representations, such as GFF and GTF, 1-based closed intervals are used whereas in the internal binary data structures, SeqAn uses 0-based half-open intervals.

A First Working Example
-----------------------

The following example shows an example of a program that reads the file with the path ``example.gff`` and prints its contents back to the user on stdout.
If you want to try out this program then create a file with the sample GFF content from above and adjust the path ``"example.gff"`` in the program below to the path to your GFF file (e.g. ``"path/to/my_example.gff"``).

.. includefrags:: extras/demos/tutorial/gff_io/example1.cpp

The program first opens a :dox:`GffStream` for reading, then one for writing.
You can read from stdin and write to stdout using ``"-"`` as the file name.

The member ``sequenceNames`` of your :dox:`GffStream` object ``gffIn`` contains the names of the reference sequences that have been seen in records so far.
This :dox:`StringSet` of :dox:`CharString` thus gets new elements as you read the Gff file.
For the translation between reference names and numeric ids, a cache is used.
The function [dox:GffStream#addSequenceName addSequenceName can be used to register the sequence name with the ``gffOut`` stream.
This will also update the cache.

Note that the example above is missing error handling.
This means that if the input format is ill-formed, error return codes are not handled appropriately and the program might do something unexpected in the case of an error.
We will fix this in `Assignment 1`_.

You can see the output of the program below when called with the input file from above.

.. code-block:: console

   ctg123  .   gene    1000    9000    .   +   .   ID=gene00001;Name=EDEN
   ctg123  .   TF_binding_site 1000    1012    .   +   .   Parent=gene00001
   ctg123  .   mRNA    1050    9000    .   +   .   ID=mRNA00001;Parent=gene00001
   ctg123  .   mRNA    1050    9000    .   +   .   ID=mRNA00002;Parent=gene00001
   ctg123  .   mRNA    1300    9000    .   +   .   ID=mRNA00003;Parent=gene00001
   ctg123  .   exon    1300    1500    .   +   .   Parent=mRNA00003
   ctg123  .   exon    1050    1500    .   +   .   Parent=mRNA00001,mRNA00002
   ctg123  .   exon    3000    3902    .   +   .   Parent=mRNA00001,mRNA00003
   ctg123  .   exon    5000    5500    .   +   .   Parent=mRNA00001,mRNA00002,mRNA00003
   ctg123  .   exon    7000    9000    .   +   .   Parent=mRNA00001,mRNA00002,mRNA00003
   ctg123  .   CDS 1201    1500    .   +   0   ID=cds00001;Parent=mRNA00001
   ctg123  .   CDS 3000    3902    .   +   0   ID=cds00001;Parent=mRNA00001
   ctg123  .   CDS 5000    5500    .   +   0   ID=cds00001;Parent=mRNA00001
   ctg123  .   CDS 7000    7600    .   +   0   ID=cds00001;Parent=mRNA00001
   ctg123  .   CDS 1201    1500    .   +   0   ID=cds00002;Parent=mRNA00002
   ctg123  .   CDS 5000    5500    .   +   0   ID=cds00002;Parent=mRNA00002
   ctg123  .   CDS         7000    7600    .   +   0   ID=cds00002;Parent=mRNA00002
   ctg123  .   CDS 3301    3902    .   +   0   ID=cds00003;Parent=mRNA00003
   ctg123  .   CDS         5000    5500    .   +   1   ID=cds00003;Parent=mRNA00003
   ctg123  .   CDS         7000    7600    .   +   1   ID=cds00003;Parent=mRNA00003
   ctg123  .   CDS 3391    3902    .   +   0   ID=cds00004;Parent=mRNA00003
   ctg123  .   CDS         5000    5500    .   +   1   ID=cds00004;Parent=mRNA00003
   ctg123  .   CDS         7000    7600    .   +   1   ID=cds00004;Parent=mRNA00003

To add error handling, we have to check return values.
The :dox:`GffStream#readRecord` call returns a status code different from ``0``, indicating an error.

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
     The functions :dox:`GffStream#readRecord` and :dox:`GffStream#writeRecord` return a status code ``int``, ``0`` on success, ``1`` on errors.
     The function :dox:`GffStream#isGood` checks whether the state of a :dox:`GffStream` is errorneous.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/gff_io/solution1.cpp

The Class :dox:`GffRecord`
--------------------------

The class :dox:`GffRecord` stores one record in a Gff file.

.. code-block:: cpp

   namespace seqan {

   class GffRecord
   {
   public:
       CharString ref;      // reference name
       __int32 rID;         // index in sequenceNames of GffStream
       CharString source;   // source free text descriptor
       CharString type;     // type of the feature
       __int32 beginPos;    // begin position of the interval
       __int32 endPos;      // end position of the interval
       float score;         // score of the annotation
       char strand;         // the strand
       char phase;          // one of '0', '1', '2', and '.'

       // The key/value list, split into a list of keys and values.
       StringSet<CharString> tagNames;
       StringSet<CharString> tagValues;

       // Returns float value for an invalid score.
       static float INVALID_SCORE();

       // Constants for marking reference id and position as invalid.
       static const __int32 INVALID_IDX = -1;
       static const __int32 INVALID_POS = -1;
   };

   }  // namespace seqan

The static members ``INVALID_POS``, ``INVALID_REFID`` store sentinel values for marking positions and reference sequence ids as invalid.
The static funtion ``INVALID_SCORE()`` returns the IEEE float "NaN" value.
In C++11, there will be a ``std::nan()`` function but for now, we need this here.

The member ``ref`` stores the contig/reference name of the genomic interval.
This information is somewhat redundant with the ``rID`` member that is filled automatically when reading from a :dox:`GffStream` such that the GffStream's ``sequenceNames[record.rID] == record.ref``.
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

        .. includefrags:: extras/demos/tutorial/gff_io/solution2.cpp

        The output is

        .. code-block:: console

           RECORDS ON CONTIGS
           ctg123  23

Assignment 3
""""""""""""

.. container:: assignment

   Generating GFF From Scratch

   Type
     Application

   Objective
     Write a program that prints the following GFF file.
     Create ``GffRecord`` objects and write them to a ``GffStream`` using ``writeRecord()``.

     .. code-block:: console

        ctg123  .   gene    1000    9000    .   +   .   ID=gene00001;Name=EDEN
        ctg123  .   TF_binding_site 1000    1012    .   +   .   Parent=gene00001

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/gff_io/solution3.cpp

GFF and GTF
-----------

The class :dox:`GffStream` transparently reads files in both GFF and GTF format.
When writing, you can select the output format with the third parameter to the constructor :dox:`GffStream` or the function :dox:`GffStream#open`.
When using ``GffStream::GFF``, GFF 3 is used, when using ``GffStream::GTF``, the GTF format.
The default is to use ``GffStream::GFF``.

The following program converts a GFF file into a GTF file.

.. includefrags:: extras/demos/tutorial/gff_io/example2.cpp

Given the GFF file at the top, the output of the example above will look as follows:

.. code-block:: console

   ctg123  .   gene    1000    9000    .   +   .   ID "gene00001"; Name "EDEN";
   ctg123  .   TF_binding_site 1000    1012    .   +   .   Parent "gene00001";
   ctg123  .   mRNA    1050    9000    .   +   .   ID "mRNA00001"; Parent "gene00001";
   ctg123  .   mRNA    1050    9000    .   +   .   ID "mRNA00002"; Parent "gene00001";
   ctg123  .   mRNA    1300    9000    .   +   .   ID "mRNA00003"; Parent "gene00001";
   ctg123  .   exon    1300    1500    .   +   .   Parent "mRNA00003";
   ctg123  .   exon    1050    1500    .   +   .   Parent "mRNA00001,mRNA00002";
   ctg123  .   exon    3000    3902    .   +   .   Parent "mRNA00001,mRNA00003";
   ctg123  .   exon    5000    5500    .   +   .   Parent "mRNA00001,mRNA00002,mRNA00003";
   ctg123  .   exon    7000    9000    .   +   .   Parent "mRNA00001,mRNA00002,mRNA00003";
   ctg123  .   CDS 1201    1500    .   +   0   ID "cds00001"; Parent "mRNA00001";
   ctg123  .   CDS 3000    3902    .   +   0   ID "cds00001"; Parent "mRNA00001";
   ctg123  .   CDS 5000    5500    .   +   0   ID "cds00001"; Parent "mRNA00001";
   ctg123  .   CDS 7000    7600    .   +   0   ID "cds00001"; Parent "mRNA00001";
   ctg123  .   CDS 1201    1500    .   +   0   ID "cds00002"; Parent "mRNA00002";
   ctg123  .   CDS 5000    5500    .   +   0   ID "cds00002"; Parent "mRNA00002";
   ctg123  .   CDS 7000    7600    .   +   0   ID "cds00002"; Parent "mRNA00002";
   ctg123  .   CDS 3301    3902    .   +   0   ID "cds00003"; Parent "mRNA00003";
   ctg123  .   CDS 5000    5500    .   +   1   ID "cds00003"; Parent "mRNA00003";
   ctg123  .   CDS 7000    7600    .   +   1   ID "cds00003"; Parent "mRNA00003";
   ctg123  .   CDS 3391    3902    .   +   0   ID "cds00004"; Parent "mRNA00003";
   ctg123  .   CDS 5000    5500    .   +   1   ID "cds00004"; Parent "mRNA00003";
   ctg123  .   CDS 7000    7600    .   +   1   ID "cds00004"; Parent "mRNA00003";

Next Steps
----------

* Continue with the :ref:`tutorial`.
