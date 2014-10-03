.. sidebar:: ToC

   .. contents::


.. _tutorial-vcf-io:

VCF I/O
=======

Learning Objective
  In this tutorial, you will learn how to use the high-level interface :dox:`VcfStream` class to read and write VCF files.

Difficulty
  Average

Duration
  1 h (45 min if you know the VCF format)

Prerequisites
  :ref:`tutorial-sequences`, exposure to the VCf format

This tutorial deals with how to easily read and write VCF files using the :dox:`VcfStream` class.
It starts out with a quick reminder on the structure of VCF files and will then continue with how to read and write VCF files and access the tags of a record.

.. important::

   Note that this tutorial is targeted at readers that already know about the VCF format.
   If you do not know about the VCF format yet then this tutorial will be harder for your to understand.
   The 1000 Genomes project hosts the `VCF format specification (v.4.1) <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_.

The VCF format allows storing genomic variants of individuals with respect to a reference.
The general file structure starts with (1) meta-information lines starting with ``##``, one (2) header line giving the names of the individuals, and (3) an arbitrary number of records.

The information of (1) and (2) will be read and written together as the "header" of the file.
For simple variants such as SNPs and small indels, each record corresponds to a variant.
More complex variants can be stored in multiple records (see the VCF standard on "breakends" for more information).

The ``vcf_io`` module of SeqAn allows the record-wise reading and writing to VCF files.
Since the structure of the fields in the VCF format often is very complex and the format undergoes changes in this respect, SeqAn only offers basic parsing functionality: The position is stored as a 0-based integer, reference names are stored in a reference name store (similar as in the :ref:`tutorial-basic-sam-bam-io` Tutorial), and the quality is stored as a ``float`` value.

The remaining fields have to be parsed from and composed as strings in the user's application.

VCF File Structure
------------------

This section gives a very brief overview of the VCF file structure.
For more details, see the `VCF format specification (v4.1) <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_.

The following is an example of a VCF file:

::

    ##fileformat=VCFv4.1
    ##fileDate=20090805
    ##source=myImputationProgramV3.1
    ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
    ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
    ##phasing=partial
    ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
    ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
    ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
    ##FILTER=<ID=q10,Description="Quality below 10">
    ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
    #CHROM POS     ID        REF    ALT     QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
    20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
    20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
    20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
    20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
    20     1234567 microsat1 GTC    G,GTCT  50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3

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

The 9 mandatory columns are followed by as many columns as there are individual.
For each individual, there is a colon-separated list of values in the order given in the FORMAT cell.

.. tip::

    1-based and 0-based positions.

    There are two common ways of specifying intervals.

     #. Start counting positions at 1 and give intervals by the first and last position that are part of the interval (closed intervals).
        For example, the interval ``[1,000; 2,000]`` starts at character 1,000 and ends at character 2,000 and includes it.
        This way is natural to non-programmers and used when giving coordinates in GFF files or genome browsers such as UCSC Genome Browser and IGV.
     #. Start counting positions at 0 and give intervals by the first position that is part of the interval and giving the position behind the last position that is part of the interval.
        The interval from above would be ``[999; 2,000)`` in this case.

    In text representations, such as VCF, 1-based closed intervals are used whereas in the internal binary data structures, SeqAn uses 0-based half-open intervals.
    When fields are reads as text, numbers are not translated, of course.

A First Working Example
-----------------------

The following example shows an example of a program that reads the file with the path ``example.vcf`` and prints its contents back to the user on stdout.
If you want to try out this program then create a file with the sample VCF content from above and adjust the path ``"example.vcf"`` in the program below to the path to your VCF file (e.g.  ``"path/to/my_example.vcf"``).

.. includefrags:: extras/demos/tutorial/vcf_io/example1.cpp

The program first opens a :dox:`VcfStream` for reading, then one for writing.
You can read from stdin and write to stdout using ``"-"`` as the file name.

The header is automatically read when a :dox:`VcfStream` is opened.
After the header has been read, it is copied over into the output stream.
Then, the input stream is read record by record and written out to the output stream.
Note that the header is written out automatically before the first variant record is written.

The alignment records are read into :dox:`VcfRecord` objects which we will focus on below.

Note that the example above is missing error handling.
This means that if the input format is ill-formed, error return codes are not handled appropriately and the program might do something unexpected in the case of an error.
We will fix this in `Assignment 1`_.

You can see the output of the program below when called with the input file from above.

::

    ##fileformat=VCFv4.1
    ##fileDate=20090805
    ##source=myImputationProgramV3.1
    ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
    ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
    ##phasing=partial
    ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
    ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
    ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
    ##FILTER=<ID=q10,Description="Quality below 10">
    ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
    #CHROM POS     ID        REF    ALT     QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
    20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
    20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
    20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
    20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
    20     1234567 microsat1 GTC    G,GTCT  50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3

To add error handling, we have to check return values.
The :dox:`VcfStream#readRecord readRecor` call returns a status code different from ``0``, indicating an error.

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
      The functions :dox:`VcfStream#readRecord` and :dox:`VcfStream#writeRecord` return a status code ``int``, ``0`` on success, ``1`` on errors.
      The function :dox:`VcfStream#isGood` checks whether the state of a :dox:`VcfStream` is errorneous.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/vcf_io/solution1.cpp

The Class :dox:`VcfRecord`
--------------------------

The class :dox:`VcfRecord` stores one record in a VCF file.
It is best explained by its definition.
Note how most fields are represented by strings:

.. code-block:: cpp

   namespace seqan {

   class VcfRecord
   {
   public:
       __int32 rID;                          // CHROM
       __int32 beginPos;                     // POS
       CharString id;                        // ID
       CharString ref;                       // REF
       CharString alt;                       // ALT
       float qual;                           // QUAL
       CharString filter;                    // FILTER
       CharString info;                      // INFO
       CharString format;                    // FORMAT
       StringSet<CharString> genotypeInfos;  // <individual1> <individual2> ..

       // Constants for marking reference id and position as invalid.
       static const __int32 INVALID_REFID = -1;
       static const __int32 INVALID_POS = -1;
       // This function returns the float value for "invalid quality".
       static float MISSING_QUAL();
   };

   }  // namespace seqan

The static members ``INVALID_POS``, ``INVALID_REFID`` store sentinel values for marking positions and reference sequence ids as invalid.
The static funtion ``MISSING_QUAL()`` returns the IEEE float "NaN" value.
In C++11, there will be a ``std::nan()`` function but for now, we need this here.

Assignment 2
""""""""""""

.. container:: assignment

   Counting Records

   Type
     Review

   Objective
     Change the result of `Assignment 1`_ by counting the number of variants for each chromosome/contig instead of writing out the records.

   Hints
     The header contains the sequence names in ``vcfIn.header.sequenceNames``.
     You can use the length of this :dox:`StringSet` of :dox:`CharString` to get the number of contigs.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/vcf_io/solution2.cpp

        The output is

        .. code-block:: console

           VARIANTS ON CONTIGS
           20  5

The Classes :dox:`VcfHeader` and :dox:`VcfHeaderRecord`
-------------------------------------------------------

The header information is stored in the class :dox:`VcfHeader`.
Objects of this class store the information present in the VCF meta information and header lines.

The class has three members: ``sequenceNames``, ``sampleNames``, and ``headerRecords``.
``sequenceNames`` and ``sampleNames`` are :dox:`StringSet StringSets` of :dox:`CharString CharStrings`.
The member ``rID`` of ``VcfRecord`` points into ``sequenceNames`` and gives the reference sequence.
The ``genotypeInfos`` member of ``VcfRecord`` has the same number of entires as ``sampleNames`` and ``record.genotypeInfos[i]`` contains the variant information for ``sampleNames[i]``.

When writing VCF files, you have to fill these three members of :dox:`VcfHeader` before writing any record.

Assignment 3
""""""""""""


.. container:: assignment

   Generating VCF From Scratch

   Type
     Application

   Objective
     Write a program that prints the VCF file from above.

   Hints
     You can convert integers into strings using the ``<sstream>`` STL header.

     .. code-block:: cpp

        #include <sstream>
        // ...
        std::stringstream ss;
        ss << 10;
        seqan::CharString str = ss.str();  // => == "10"
        // To reset ss, we need two calls:
        ss.str("");  // Remove contents.
        ss.clear();  // Reset any error bits.


   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/vcf_io/solution3.cpp

Next Steps
----------

* Continue with the :ref:`tutorial`
