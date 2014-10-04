.. sidebar:: ToC

   .. contents::


.. _tutorial-basic-sequence-io:

Basic Sequence I/O
==================

Learning Objective
  You will learn how to read and write sequence files (FASTA and FASTQ) using a simple, high-level API in the SeqAn library.
  This includes reading and writing compressed files.

Difficulty
  Basic

Duration
  30 min

Prerequisites
  :ref:`tutorial-sequences`

This tutorial explains how to read and write sequence files using the :dox:`SequenceStream` class.

The SeqAn library's I/O functionality (i.e. reading from and writing to files) is organized in layers and :dox:`SequenceStream` is in the highest layer: It provides an easy-to-use API for reading and writing sequence files, automatically compressing and decompressing data, and support for different sequence formats.
This flexibility comes at a slight cost of performance, compared to using the more low-level APIs.

The lower layers are responsible for providing raw file I/O functionality and adding parsing functionality.
Their usage is not part of this tutorial and is explainend in the Tutorials :ref:`tutorial-file-io`, :ref:`tutorial-sequence-file-io`, and :ref:`tutorial-parsing`.

After completing the tutorial, you will be able to read and write sequence files in the formats supported by the :dox:`SequenceStream` class.

A First Working Example
-----------------------

Let us start out with a minimal working example.
The following small program will read the file ``example.fa`` (which we will create later) from the current directory and print out the identifier and the sequence of the first record.

.. includefrags:: extras/demos/tutorial/basic_sequence_io/example1.cpp

We use the :dox:`SequenceStream::SequenceStream SequenceStream constructor` with one parameter, the path to the file we want to read.
This will open the file and guess the file format from the file contents.
The class will read the first some thousand characters and try to guess what format the file is.
Note that you can also pass ``"-"`` as the file name.
This will open the standard input when reading and the standard output when writing.

After construction, the ``seqStream`` object is ready for reading.
We use the function :dox:`SequenceStream#readRecord` to read the first record from the file.
:dox:`SequenceStream#readRecord` always reads the next record in the file.

.. tip::

   FASTA/FASTQ and Record-Based Files

   Most files in bioinformatics have a record-based structure.
   Often, a file format requires or allows for a header that contains information about the file format.
   Then, the file contains a list of records, one after another.

   The FASTA and FASTQ formats do not have a header but only contain lists of records.
   For example, a FASTQ record contains the sequence id, the sequence characters, and a quality value for each character.

Note that we do not have to close the file manually.
The :dox:`SequenceStream` object will automatically close any open files when it goes out of scope and it is destructred.
If you want to force a file to be closed, you can use the function :dox:`SequenceStream#close`.

Adding Error Handling
---------------------

Now, create a new FASTA file named ``example.fa`` in a directory of your choice with the following content:

::

    >seq1
    CCCCCCCCCCCCCCC
    >seq2
    CGATCGATC
    >seq3
    TTTTTTT

Then, copy the program above into new application ``basic_seq_io_example``, adjust the path ``"example.fa"`` to the just created FASTA file, compile the program, and run it.
For example, if you stored the file ``example.fa`` in ``/home/username/example.fa``, you replace the line ``seqan::SequenceStream seqStream("example.fa");`` from above with ``seqan::SequenceStream seqStream("/home/username/example.fa");``.
You should see the following output:

.. code-block:: console

   # basic_seq_io
   seq1    CCCCCCCCCCCCCCC

Assignment 1
""""""""""""

.. container:: assignment

   Type ::
     Review
   Objective ::
     Adjust the program above to use the first command line parameter ``argv[1]``, i.e. the first argument.
     Check that there actually is such an argument (``argc >= 2``) and let ``main()`` return ``1`` otherwise.
   Solution ::
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/basic_sequence_io/solution1.cpp

Our program is very simple but there is one large problem.
Anything can go wrong during file I/O and have not used any means to handle such errors.
Possible errors include: the file permissions forbid a certain operations, the file does not exist, there is a disk reading error, a file read from a remote location gets deleted while we are reading from it, or there is a physical error in the hard disk.

Let us add some error handling.
At the very least, we should detect errors.
If possible, we should try to recover from the error (sometimes it is possible to return default values instead of loading values from a file) or otherwise stop the current task in an organized fashion and notify the user about the problem.

We can use the Function :dox:`SequenceStream#isGood` to check whether the :dox:`SequenceStream` object is ready for any more reading.
After the creation of the object, this function indicates whether the file could be opened successfully by returning ``true``.
The function :dox:`SequenceStream#readRecord` returns an ``int`` that indicates whether the reading was successful.
If everything went fine, it returns ``0``, and a different value otherwise.

Note that :dox:`SequenceStream#isGood` queries the state of the stream and returns a ``bool`` indicating whether the stream is ready for reading/writing (``true`` for "is good" and ``false`` for "is not good").
:dox:`SequenceStream#readRecord`, on the other hand, returns an ``int`` indicating whether there was any error (``0`` for "is good" and a non-\ ``0`` value for "is not good", as it is customary in Unix programming).

The program will now read as follows:

.. includefrags:: extras/demos/tutorial/basic_sequence_io/example2.cpp

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     Change your program from above to perform these checks, too.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/basic_sequence_io/solution2.cpp

Assignment 3
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Change your program from above to loop over all sequences and print them in the same fashion.

   Hint
     You can use the function :dox:`SequenceStream#atEnd` to check whether a :dox:`SequenceStream` object is at the end of the file.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/basic_sequence_io/solution3.cpp

After completing Assignment 3, you should be able to run your program on the example file we created above and see the following output:

.. code-block:: console

    # basic_seq_io_example example.fa
    seq1    CCCCCCCCCCCCCCC
    seq2    CGATCGATC
    seq3    TTTTTTT

The Interface for Reading
-------------------------

There are three major usage patterns for sequence I/O:

#. We want to read **all records** from the file into memory, for example for building an index.
#. We want to read the file into memory **record by record**, so the memory usage is minimal.
   We could then perform some computation on each record, e.g. search it in an index.
#. We want to read a **batch of records** into memory, e.g. 100k records at a time.
   Then, we perform some computation on the records, for example in parallel with 4 threads on 25k records each.

These use cases are supported by the functions :dox:`SequenceStream#readAll`, :dox:`SequenceStream#readRecord`, and :dox:`SequenceStream#readBatch`.

Each of these functions is available in two variants.
The first accepting only the sequence identifier and sequence characters besides the :dox:`SequenceStream` object and the second also accepting the a :dox:`CharString` for the PHRED base qualities.
If a file does not contain any qualities and the function variant with quality values is used then the quality strings are returned as empty.
When writing a file with qualities and the function variant without quality values is used then the qualities are written out as ``'I'``, i.e. PHRED score 40.

When :dox:`DnaQ` or :dox:`Dna5Q` are used, then you should use the function variant without a parameter for qualities.
The qualities are simply stored directly in the sequence characters.

As to be expected, when there are characters in the file that are not valid characters in the :dox:`String` then the alphabet-dependent conversion is performed.
For example, for :dox:`Dna` and :dox:`Rna` this means a conversion of the invalid character to ``'A'``, and for :dox:`Dna5 Dna5 and [dox:Rna5 Rna5` this means a conversion to ``'N'``.

Here is an example for using :dox:`SequenceStream#readRecord`:

.. code-block:: cpp

   seqan::CharString id;
   seqan::Dna5String seq;
   seqan::CharString qual;
   int res = 0;

   seqan::SequenceStream seqStream("in.fq");

   res = readRecord(id, seq, seqStream);
   res = readRecord(id, seq, qual, seqStream);

The functions :dox:`SequenceStream#readAll` and :dox:`SequenceStream#readBatch` use :dox:`StringSet` instead of :dox:`String`.
The function :dox:`SequenceStream#readBatch` reads up to the given number of records.
It is not an error if there are less records.

.. code-block:: cpp

   seqan::StringSet<seqan::CharString> ids;
   seqan::StringSet<seqan::Dna5String> seqs;
   seqan::StringSet<seqan::CharString> quals;
   int res = 0;

   seqan::SequenceStream seqStream("in.fq");

   res = readAll(ids, seqs, seqStream);
   res = readAll(ids, seqs, quals, seqStream);

   res = readBatch(ids, seqs, seqStream, 10);
   res = readBatch(ids, seqs, quals, seqStream, 10);

Assignment 4
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Change your result of Assignment 3 to use the variant of :dox:`SequenceStream#readRecord` that also reads in the qualities and writes them next to the sequences.
     Create the following FASTQ file ``example.fq``.

     ::

         @seq1
         CCCCCCCCCCCCCCC
         +
         IIIIIHIIIIIIIII
         @seq2
         CGATCGATC
         +
         IIIIIIIII
         @seq3
         TTTTTTT
         +
         IIIIHHG

     When your program is called on this file, the result should look as follows.

     .. code-block:: console

        # basic_seq_io_example example.fq
        seq1    CCCCCCCCCCCCCCC    IIIIIHIIIIIIIII
        seq2    CGATCGATC    IIIIIIIII
        seq3    TTTTTTT      IIIIHHG

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/basic_sequence_io/solution4.cpp

The Interface for Writing
-------------------------

Now that you know how to read sequence files, writing them will come easy to you.
We can open files for writing by giving ``seqan::SequenceStream::WRITE`` as the second parameter to the :dox:`SequenceStream::SequenceStream SequenceStream constructor`.
Create a new SeqAn app ``basic_seq_io_example2`` in your sandbox and change the C++ file ``basic_seq_io_example2.cpp`` in this application to have the content below.
This program already has all the bells and whistles for error checking.

.. includefrags:: extras/demos/tutorial/basic_sequence_io/example3.cpp

The first lines are similar to those in the solution to Assignment 4.
However, instead of opening the file using ``seqan::SequenceStream seqStream(argv[1]);``, we use ``seqan::SequenceStream seqStream(argv[1], seqan::SequenceStream::WRITE);``.
this opens the file with the name in ``argv[1]`` for writing instead of for reading.
Also, instead of reading records, we write one record.

The program writes out one sequence with id "seq1" and the contents "CGAT" to the file given on the command line.
Note that :dox:`SequenceStream` will guess the format from the file name.
A file ending in ``.fa`` and ``.fasta`` mean FASTA, ``.fq`` and ``.fastq`` means FASTQ.
Optionally, you can force to use any file format with the third parameter to the :dox:`SequenceStream::SequenceStream SequenceStream constructor`.

Let us try out the program from above:

.. code-block:: console

   # basic_seq_io_example2 out.fa
   # cat out.fa
   >seq1
   CGAT
   # basic_seq_io_example2 out.fq
   # cat out.fq
   @seq
   CGAT
   +
   IIII

Assignment 5
""""""""""""

.. container:: assignment

   Type
     Reproduction

   Objective
     Change the program from above to write out a second sequence.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/basic_sequence_io/solution5.cpp

There are two functions for writing to sequence files using :dox:`SequenceStream`.
One, :dox:`SequenceStream#writeRecord`, for writing one sequence record from :dox:`String Strings`, and another one, :dox:`SequenceStream#writeAll`, for writing all sequences from :dox:`StringSet StringSets`.

Again, they come in one variant with and another variant without base qualities.
When writing to a FASTQ file using the function without qualities, the PHRED score 40 is written for each character (``'I'``) and when writing to a FASTA file with the variant with qualities, the qualities are ignored.
When using :dox:`DnaQ` or :dox:`Dna5Q`, the variant without qualities parameter writes out the qualities stored in the sequence characters themselves.

Here is an example for using :dox:`SequenceStream#writeRecord`:

.. code-block:: cpp

   seqan::CharString id;
   seqan::Dna5String seq;
   seqan::CharString qual;

   seqan::SequenceStream seqStream("out.fq", seqan::SequenceStream::WRITE);

   res = writeRecord(seqStream, id, seq);
   res = writeRecord(seqStream, id, seq, qual);

And here is an example for using :dox:`SequenceStream#writeAll`:

.. code-block:: cpp

   seqan::StringSet<seqan::CharString> ids;
   seqan::StringSet<seqan::Dna5String> seqs;
   seqan::StringSet<seqan::CharString> quals;

   seqan::SequenceStream seqStream("out.fq", seqan::SequenceStream::WRITE);

   res = writeAll(seqStream, ids, seqs);
   res = writeAll(seqStream, ids, seqs, quals);

Assignment 6
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Change the result of Assignment 5 to store the data for the two records in :dox:`StringSet StringSets` and write them out using :dox:`SequenceStream#writeAll`.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/basic_sequence_io/solution6.cpp

Compressed Files
----------------

Using compressed files is simple.
When opening a file for reading, :dox:`SequenceStream` will automatically detect whether the file is compressed or not, the same it detects the sequence file format for you.
If you run into problems here, make sure that you have zlib and/or libbz2 installed (see `Dependencies on Compression Libraries`_ below).

When opening a file for writing, :dox:`SequenceStream` will infer the compression type (gzip, bzip2, or plain text only) and the file format (FASTA or FASTQ) from the file ending.
First, the file type is guessed: A file ending in ``.gz`` means "gzip-compressed", one ending in ``.bz2`` means "bzip2-compressed".
Then, the ``.gz`` or ``.bz2`` suffix is ignored when guessing the file format.
A path ending in ``.fa`` and ``.fasta`` mean FASTA, ``.fq`` and ``.fastq`` mean FASTQ.
Since the suffixes ``.gz`` and ``.bz2`` are ignored, ``.fa.gz``, ``.fa.bz2``, ... mean FASTA too and ``.fq.gz``, .\ ``fq.bz2``, ... mean FASTQ.

File type detection from standard input is currently limited to either gzip-compressed or plain-text data.

Note that you can also use additional parameters in the :dox:`SequenceStream::SequenceStream SequenceStream constructor` to force a certain file type and file format when writing.
You can also force a certain file type and format when reading but this is only helpful in the few instances where the automatic detection fails.

This means that all the examples and your solutions to the assignments from above **already have compression support built-in**, if the compression libraries are available.

Dependencies on Compression Libraries
-------------------------------------

For accessing compressed files, you need to have zlib installed for reading ``.gz`` files and libbz2 for reading ``.bz2`` files.

If you are using Linux or Mac Os X and you followed the :ref:`tutorial-getting-started` tutorial closely then you should have already installed the necessary libraries.
On Windows, you will need to follow :ref:`how-to-install-contribs-on-windows` to get the necessary libraries.

You can check whether you have installed the libraries to use zlib and libbz2 by running CMake again.
Simply call ``cmake .`` in your build directory.
At the end of the output, there will be a section "SeqAn Features".
If you can read ``ZLIB - FOUND`` and ``BZIP2 - FOUND`` then you can use zlib and libbz2 in your programs.

Congratulations, you have now learned to write simple and robust sequence I/O code using SeqAn!

Next Steps
----------

* Read the Wikipedia articles about the `FASTA file format <http://en.wikipedia.org/wiki/FASTA_format>`_ and the `FASTQ file format and quality values <http://en.wikipedia.org/wiki/FASTQ_format>`_ to refresh your knowledge.
* Read the :ref:`tutorial-indexed-fasta-io` tutorial to learn how to read FASTA files efficiently in a random-access fashion.
* Continue with the :ref:`tutorial`.
