.. sidebar:: ToC

    .. contents::

.. _tutorial-io-vcf-io:

VCF I/O
=======

Learning Objective
  In this tutorial, you will learn how to read and write VCF files.

Difficulty
  Average

Duration
  1 h (45 min if you know the VCF format)

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-io-input-output-overview`, `VCF Format Specification (v4.2) <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_

Overview
--------

This tutorial shows how to read and write VCF files using the :dox:`VcfFileIn` and :dox:`VcfFileOut` classes.
It starts out with a quick reminder on the structure of VCF files and will then continue with how to read and write VCF files and access the tags of a record.

.. important::

   Note that this tutorial is targeted at readers that already know about the VCF format.
   If you do not know about the VCF format yet, then this tutorial will be harder for your to understand.

The VCF format allows storing genomic variants of individuals with respect to a reference.
The general file structure starts with 

#. several meta-information lines starting with ``##``,
#. one header line giving the names of the individuals, and
#. an arbitrary number of records.

The information of (1) and (2) will be read and written together as the "header" of the file.
For simple variants such as SNPs and small indels, each record corresponds to a variant.
More complex variants can be stored in multiple records (see the VCF standard on "breakends" for more information).

The ``vcf_io`` module of SeqAn allows to read and write VCF files record-wise.
Since the structure of the fields in the VCF format often is very complex and the format undergoes changes in this respect, SeqAn only offers basic parsing functionality: The position is stored as a 0-based integer, reference names are stored in a reference name store (similar as in the :ref:`tutorial-io-sam-bam-io` Tutorial), and the quality is stored as a ``float`` value.

The remaining fields have to be parsed from and composed as strings in the user's application.

VCF Format
----------

This section gives a very brief overview of the VCF file structure.
For more details, see the `VCF Format Specification (v4.2) <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_.

The following is an example of a VCF file:

.. includefrags:: demos/tutorial/vcf_io/example.vcf

The file starts with meta information lines (starting with ``##``) with a key/value structure.
The most important lines have the keys **contig**, **INFO**, **FILTER**, and **FORMAT**.

contig
  Lines with this key list the contigs of the reference genome.``

INFO
  These lines give valid keys (and the format of the values) for the INFO column.

FILTER
  Valid values of the FILTER column.

FORMAT
  Valid entries for the INFO column.

The meta information lines are followed by the header line which gives the names of the first 9 columns which are always the same (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) and a non-empty list of sample names.
The columns are separated by spaces.

The header line is followed by the records which contains a value for each column in the header.

CHROM
  Name of the chromosome/reference sequence that the variant lies on.

POS
  The 1-based position of the variant.

ID
  A name of the variant.
  ``.`` is used if no name is available.

REF
  The value of the reference allele.

ALT
  The alternate allele values (multiple values are comma-separated).

QUAL
  Quality value of the call (float).

FILTER
  A value for the filter result (given in a ``FILTER`` meta information line).

INFO
  Information about a variant.

FORMAT
  Colon-separated list of entries that are found for each variant.

The 9 mandatory columns are followed by as many columns as there are individuals.
For each individual, there is a colon-separated list of values in the order given in the FORMAT cell.

.. tip::

    1-based and 0-based positions.

    There are two common ways of specifying intervals.

     #. Start counting positions at 1 and give intervals by the first and last position that are part of the interval (closed intervals).
        For example, the interval ``[1000; 2000]`` starts at character 1000 and ends at character 2000 and includes it.
        This way is natural to non-programmers and used when giving coordinates in GFF files or genome browsers such as UCSC Genome Browser and IGV.
     #. Start counting positions at 0 and give intervals by the first position that is part of the interval and giving the position behind the last position that is part of the interval.
        The interval from above would be ``[999; 2000)`` in this case.

    In text representations, such as VCF, 1-based closed intervals are used whereas in the internal binary data structures, SeqAn uses 0-based half-open intervals.
    When fields are read as text, numbers are not translated, of course.

A First Working Example
-----------------------

The following example shows a program that reads the file ``example.vcf`` and prints its contents back to the user on standard output.

.. includefrags:: demos/tutorial/vcf_io/example1.cpp

The program first opens a :dox:`VcfFileIn` for reading the file, then a :dox:`VcfFileOut` for writing it.
First, the header is copied by means of a :dox:`VcfHeader` object that we will see below.
Then, the input file is read record by record and written out to standard output.
The alignment records are read into :dox:`VcfRecord` objects which we will focus on below.

The output of the example program looks as follows:

.. includefrags:: demos/tutorial/vcf_io/example1.cpp.stdout

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Reproduction

   Objective
     Create a file with the sample VCF content from above and adjust the path ``"example.vcf"`` to the path to your SAM file (e.g. ``"/path/to/my_example.sam"``).

   Solution
      .. container:: foldable

         .. includefrags:: demos/tutorial/vcf_io/solution1.cpp

Accessing the Header
--------------------

Sequence information from the VCF header is stored in the :dox:`VcfIOContext`.
The :dox:`VcfIOContext` of a :dox:`VcfFileIn` can be accessed through the function :dox:`FormattedFile#context`.
The VCF sequence informations can be in turn accessed through functions :dox:`VcfIOContext#contigNames` and :dox:`VcfIOContext#sampleNames`.
All remaining VCF header information is stored in the class :dox:`VcfHeader`.

Accessing the Records
---------------------

The class :dox:`VcfRecord` stores one record in a VCF file.
It is best explained by its definition.
Note how most fields are represented by :dox:`CharString Strings`:

.. includefrags:: demos/tutorial/vcf_io/base.cpp
      :fragment: vcfRecord

The static members ``INVALID_POS`` and ``INVALID_REFID`` store sentinel values for marking positions and reference sequence ids as invalid.
The static funtion ``MISSING_QUAL()`` returns the IEEE float "NaN" value.

.. tip::
   A :dox:`VcfRecord` is linked to a reference sequence by the field ``rID`` and to samples by ``genotypeInfos``.
   The sequence information is stored in the VCF header and kept in the :dox:`VcfIOContext`.


Assignment 2
""""""""""""

.. container:: assignment

   Counting Records

   Type
     Review

   Objective
     Change the result of `Assignment 1`_ by counting the number of variants for each chromosome/contig instead of writing out the records.

   Hints
     The reference sequence information from the VCF header is stored inside the :dox:`VcfIOContext` of its :dox:`VcfFileIn`.
     You can obtain the number of contigs from the :dox:`ContainerConcept#length` of the :dox:`VcfIOContext#contigNames`.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/vcf_io/solution2.cpp

        The output is

        .. includefrags:: demos/tutorial/vcf_io/solution2.cpp.stdout

Creating a New File
-------------------

Assignment 3
""""""""""""


.. container:: assignment

   Generating VCF From Scratch

   Type
     Application

   Objective
     Write a program that manually creates the VCF file from above and than prints it back on standard output.

   Hint
     * use :dox:`VcfHeaderRecord` to create a header record that can be appended to the :dox:`VcfHeader`
     * use ``appendValue()`` to append information to the :dox:`VcfIOContext` or the :dox:`VcfHeader`
     * use the direct member access operator ``.`` if you want to access the information of a :dox:`VcfRecord`

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/vcf_io/solution3.cpp

Next Steps
----------

* Continue with the :ref:`tutorial`
