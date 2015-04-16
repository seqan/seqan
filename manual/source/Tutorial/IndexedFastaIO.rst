.. sidebar:: ToC

   .. contents::


.. _tutorial-indexed-fasta-io:

Indexed FASTA I/O
=================

Learning Objective ::
  In this tutorial, you will learn how to use a FASTA Index file (``.fai``) for indexed random-access to FASTA files.
  This is useful for retrieving regions (e.g. ``chr1:123-10004``) or single sequences (e.g. ``chr1``) from FASTA files quickly.

Difficulty
  Average

Duration
  30 min

Prerequisites
  :ref:`tutorial-sequences`

The idea of FASTA index files (*FAI*) comes from the `samtools <http://samtools.sourceforge.net/samtools.shtml>`_ program by Heng Li.
The program provides a command ``samtools faidx`` for rapidly accessing parts of a large FASTA file (e.g. querying for the first chromosome by the identifier "chr1" or querying for 900 characters starting from character 100 (1-based) by ``chr1:100-1,000``).
To do this, the program creates an index file that contains one entry for each sequence.
If the FASTA file is named ``path/sequence.fasta``, the index file is usually named ``path/sequence.fasta.fai``.

Using such index files, it is possible to rapidly read parts of the given sequence file.
The module ``<seqan/seq_io.h>`` allows to create and read such ``.fai`` index files and exposes an API to read parts randomly of FASTA file.

.. note::

    FASTA/FASTQ Meta Data and Sequence Ids

    FASTA and FASTQ files have one metadata record for each sequence.
    This usually contains the sequence name but sometimes a lot of additional information is stored.
    There is no consensus for the metadata.

    However, it is common to store the sequence identifier (*id*) at the beginning of the metadata field before the first space.
    The id is unique to the whole file and often identifies the associated sequence uniquely in a database (see section Sequence Identifiers on the `Wikipedia FASTA format <http://en.wikipedia.org/wiki/FASTA_format>`_ page).

    While not documented anywhere explicitely, **only the characters up to the first space are used as identifiers** by widely used tools such as `BWA <http://bio-bwa.sourceforge.net/>`_.
    Only the identifier is carried over into files generated from the input files (BWA uses the sequence id from the genome FASTA to identify the contig/chromosome and the read id as the read name in the SAM output).

How Does It Work?
-----------------

There are two requirements that a FASTA file has to fulfill to work with the FAI scheme.
For each sequence in the FASTA file, the number of characters stored per line has to be the same and the number of bytes per line has to be the same.
The first restriction speaks for itself, the second restriction means that the same line ending character has to be used and no line should contain any additional spaces.

The index file then stores records of sequence identifier, length, the offset of the first sequence character in the file, the number of characters per line, and the number of bytes per line.
With this information, we can easily compute the byte offset of the i-th character of a sequence in a file by looking at its index record.
We skip to this byte offset in the file and from there, we can read the necessary sequence characters.

Building the Index
------------------

The class :dox:`FaiIndex` allows for building and loading FAI indices.
fo build such an index, we use the function :dox:`FaiIndex#build` of the class :dox:`FaiIndex`.
The first parameter is the :dox:`FaiIndex` object, the second is the path to the FASTA file.
The function returns a ``bool`` indicating whether the mapping was successful (``true`` on success, ``false`` on failure).

.. code-block:: cpp

   #include <seqan/seq_io.h>

   seqan::FaiIndex faiIndex;
   if (!build(faiIndex, "path/to/file.fasta"))
       std::cerr << "ERROR: Could not build the index!\n";

There is an alternative variant of this function where you can pass the path to the FAI file that is to be built as third parameter.
The FAI file name will be stored in the :dox:`FaiIndex`.

.. code-block:: cpp

   #include <seqan/seq_io.h>

   seqan::FaiIndex faiIndex;
   if (!build(faiIndex, "path/to/file.fasta", "another/path/file.fasta.fai"))
       std::cerr << "ERROR: Could not build the index!\n";

We can write out the index after building it using the function :dox:`FaiIndex#save`:

.. code-block:: cpp

   #include <seqan/seq_io.h>

   seqan::FaiIndex faiIndex;
   // ... index building here ...

   if (!save(faiIndex, "path/to/file.fasta.fai"))
       std::cerr << "ERROR: Could not save the index to file!\n";

Assignment 1
""""""""""""

.. container:: assignment

   Building a FAI index

   Type
     Application

   Objective
      Write a small program ``build_fai`` that takes one parameter from the command line, the path to a FASTA file.
      The program should then build a FAI index and write it out.

   Hints
     .. container:: foldable

	Using the two-parameter variant of :dox:`FaiIndex#build` is good enough.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/indexed_fasta_io/solution1.cpp

Using the Index
---------------

To load a FAI file, we use the function :dox:`FaiIndex#open`: We pass the :dox:`FaiIndex` object as the first and the path to the FASTA file as the second parameter.
The function returns a ``bool`` indicating whether the mapping was successful (``true`` on success, ``false`` on failure).

.. code-block:: cpp

   #include <seqan/seq_io.h>

   seqan::FaiIndex faiIndex;
   if (!open(faiIndex, "path/to/file.fasta"))
       std::cerr << "ERROR: Could not open FAI index path/to/file.fasta.fai\n";

In the example above, the FAI file ``"path/to/file.fasta.fai"`` would be
loaded. Optionally, we can specify an extra path to the FAI file:

.. code-block:: cpp

   #include <seqan/seq_io.h>

   seqan::FaiIndex faiIndex;
   if (!open(faiIndex, "path/to/file.fasta", "path/to/index.fai"))
       std::cerr << "ERROR: Could not load FAI index path/to/index.fai\n";

After loading the index, we can then use the index to map a sequence id to its (zero-based) position (a position *i* meaning that it is the *i*-th sequence) in the FASTA file using :dox:`FaiIndex#getIdByName`.
The function gets the :dox:`FaiIndex` to use, the id of the sequence, and an ``unsigned`` position as parameters.
It returns a ``bool`` indicating whether the mapping was successful (``true`` on success, ``false`` on failure).

.. code-block:: cpp

   unsigned idx = 0;
   if (!getIdByName(faiIndex, "chr1", idx))
       std::cerr << "ERROR: FAI index has no entry for chr1.\n";

Once we have the index for the sequence in the FASTA file, we can then query the :dox:`FaiIndex` for the length of the sequence using :dox:`FaiIndex#sequenceLength`, get the whole sequence using :dox:`FaiIndex#readSequence`, or get just a part of the sequence using :dox:`FaiIndex#readRegion`.

.. code-block:: cpp

   unsigned seqLength = sequenceLength(faiIndex, idx);

   // Load first 1000 characters of chr1.
   seqan::CharString seqChr1Prefix;
   readRegion(seqChr1Prefix, faiIndex, idx, 0, 1000);

   // Load all of chr1.
   seqan::CharString seqChr1;
   readSequence(seqChr1, faiIndex, idx);

The sequence length can be determined by only looking at the index.
When loading the sequence or a sequence infix, only the relevant part of the file will be touched.
Thus, only the minimal amount of memory, time, and disk I/O is used.

Assignment 2
""""""""""""

.. container:: assignment

   Using the FAI index

   Type
     Application

   Objective
     Write a small program ``query_fai`` that takes four parameters from the command line:
     A path to a FASTA file, the id of the sequence, a begin and an end position.
     The program should then read the given infix of the given sequence from the file and print it to stdout.

  Hint
    .. container:: foldable

       Use the function :dox:`lexicalCast` to convert strings of numbers into integers.

  Solution
    .. container:: foldable

       The program appears to be very long, but most is error handling, as usual with robust I/O code.

       .. includefrags:: demos/tutorial/indexed_fasta_io/solution2.cpp


Next Steps
----------

* Read the Wikipedia articles about the `FASTA file format <http://en.wikipedia.org/wiki/FASTA_format>`_ and the `FASTQ file format and quality values <http://en.wikipedia.org/wiki/FASTQ_format>`_ to refresh your knowledge.
* Read the API documentation of the :dox:`GenomicRegion` class for storing regions (sequence identifier, start and end position).
  There also is functionality for parsing strings like ``chr1:2,032-3,212`` into :dox:`GenomicRegion` objects.
* Continue with the :ref:`tutorial`.
