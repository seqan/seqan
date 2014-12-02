.. sidebar:: ToC

   .. contents::


.. _tutorial-sequence-io:

Sequence I/O
============

Learning Objective
  You will learn how to read and write sequence files in FASTA, FASTQ, EMBL or GenBank format.

Difficulty
  Basic

Duration
  30 min

Prerequisites
  :ref:`tutorial-sequences`, :ref:`tutorial-input-output-overview`

This tutorial explains how to read and write sequence files using the :dox:`SeqFileIn` and :dox:`SeqFileOut` classes.
These classes provide an API for accessing sequence files in different file formats, either compressed or uncompressed.


FASTA/FASTQ Format
------------------

FASTA/FASTQ are record-based files.
For example, a FASTQ record contains the sequence id, the sequence characters, and a quality value for each character.

Create a new FASTA file named ``example.fa`` in a directory of your choice with the following content:

::

    >seq1
    CCCCCCCCCCCCCCC
    >seq2
    CGATCGATC
    >seq3
    TTTTTTT


A First Working Example
-----------------------

Let us start out with a minimal working example.
The following small program will read the file ``example.fa`` just created and print out the identifier and the sequence of the first record.

.. includefrags:: demos/tutorial/seq_io/example1.cpp

We call the :dox:`SeqFileIn::SeqFileIn SeqFileIn constructor` with the path to the file to read.
Successively, we call the function :dox:`SeqFileIn#readRecord` to read the first record from the file.
Note that :dox:`SeqFileIn#readRecord` always moves ahead to the next record in the file.

We do not have to close the file manually.
The :dox:`SeqFileIn` object will automatically close any open files when it goes out of scope and it is destructred.
If you want to force a file to be closed, you can use the function :dox:`SeqFileIn#close`.

Copy the program above into new application ``basic_seq_io_example``, adjust the path ``"example.fa"`` to the just created FASTA file, compile the program, and run it.
For example, if you stored the file ``example.fa`` in ``/home/username/example.fa``, you replace the line ``seqan::SeqFileIn seqFileIn("example.fa");`` from above with ``seqan::SeqFileIn seqFileIn("/home/username/example.fa");``.

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

        .. includefrags:: demos/tutorial/seq_io/solution1.cpp

Handling Errors
---------------

The program will now read as follows:

.. includefrags:: demos/tutorial/seq_io/example2.cpp

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     Change the above program to catch IOError exceptions.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/seq_io/solution2.cpp


Reading Files
-------------

These use cases are supported by the functions :dox:`SeqFileIn#readRecord` and :dox:`SeqFileIn#readRecords`.

These functions are available in two variants:
#. the first variant accepts only the sequence identifier and sequence characters, besides the :dox:`SeqFileIn` object;
#. the second variant accepts an additional :dox:`CharString` for a PHRED base quality string.

If the second variant is used on a file not containing any qualities, the quality strings are returned empty.
Note that invalid characters in the file will be signaled by :dox:`SeqFileIn#readRecord` via parsing exceptions.

.. tip::

    When :dox:`DnaQ` or :dox:`Dna5Q` :dox:`String Strings` are used, then you should use the variant without qualities.
    The qualities are simply stored directly in the sequence characters.

Here is an example for using :dox:`SeqFileIn#readRecord`:

.. code-block:: cpp

   seqan::CharString id;
   seqan::Dna5String seq;
   seqan::CharString qual;

   seqan::SeqFileIn seqFileIn("in.fq");

   readRecord(id, seq, seqFileIn);
   readRecord(id, seq, qual, seqFileIn);

The function :dox:`SeqFileIn#readRecords` use :dox:`StringSet` instead of :dox:`String`.
By default, it reads all remaining records.
Optionally, one can specify a batch of records to be read, e.g. 10 records.

.. code-block:: cpp

   seqan::StringSet<seqan::CharString> ids;
   seqan::StringSet<seqan::Dna5String> seqs;
   seqan::StringSet<seqan::CharString> quals;

   seqan::SeqFileIn seqFileIn("in.fq");

   readRecords(ids, seqs, seqFileIn, 10);
   readRecords(ids, seqs, quals, seqFileIn, 10);

   readRecords(ids, seqs, seqFileIn);
   readRecords(ids, seqs, quals, seqFileIn);


Assignment 3
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Change your program from above to loop over all sequences and print them in the same fashion.

   Hint
     You can use the function :dox:`SeqFileIn#atEnd` to check whether a :dox:`SeqFileIn` object is at the end of the file.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/seq_io/solution3.cpp

After completing Assignment 3, you should be able to run your program on the example file we created above and see the following output:

.. code-block:: console

    # basic_seq_io_example example.fa
    seq1    CCCCCCCCCCCCCCC
    seq2    CGATCGATC
    seq3    TTTTTTT


Assignment 4
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Change your result of Assignment 3 to use the variant of :dox:`SeqFileIn#readRecord` that also reads in the qualities and writes them next to the sequences.
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

        .. includefrags:: demos/tutorial/seq_io/solution4.cpp

Writing Files
-------------

We can write sequence files with the :dox:`SeqFileOut` class.

Create a new SeqAn app ``basic_seq_io_example2`` in your sandbox and change the C++ file ``basic_seq_io_example2.cpp`` in this application to have the content below.
This program already has all the bells and whistles for error checking.

.. includefrags:: demos/tutorial/seq_io/example3.cpp

The first lines are similar to those in the solution to Assignment 4.
However, instead of reading records, we write one record.

The program writes out one sequence with id "seq1" and the contents "CGAT" to the file given on the command line.
Note that :dox:`SeqFileOut` will guess the format from the file name.
A file ending in ``.fa`` and ``.fasta`` mean FASTA, ``.fq`` and ``.fastq`` means FASTQ.

.. COMMENT Optionally, you can force to use any file format with the third parameter to the :dox:`SequenceStream::SequenceStream SequenceStream constructor`.
.. COMMENT When writing a file with qualities and the function variant without quality values is used then the qualities are written out as ``'I'``, i.e. PHRED score 40.

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

        .. includefrags:: demos/tutorial/seq_io/solution5.cpp

As for reading, there are two functions for writing sequence files: :dox:`SeqFileOut#writeRecord` and :dox:`SeqFileOut#writeRecords`.

Again, both functions come in two variants: with or without base qualities.
When writing to a FASTQ file using the function without qualities, the PHRED score 40 is written for each character (``'I'``) and when writing to a FASTA file with the variant with qualities, the qualities are ignored.

When using :dox:`DnaQ` or :dox:`Dna5Q`, the variant without qualities parameter writes out the qualities stored in the sequence characters themselves.

Here is an example for using :dox:`SeqFileOut#writeRecord`:

.. code-block:: cpp

   seqan::CharString id;
   seqan::Dna5String seq;
   seqan::CharString qual;

   seqan::SeqFileOut seqFileOut("out.fq");

   writeRecord(seqFileOut, id, seq);
   writeRecord(seqFileOut, id, seq, qual);

And here is an example for using :dox:`SeqFileOut#writeRecords`:

.. code-block:: cpp

   seqan::StringSet<seqan::CharString> ids;
   seqan::StringSet<seqan::Dna5String> seqs;
   seqan::StringSet<seqan::CharString> quals;

   seqan::SeqFileOut seqFileOut("out.fq");

   writeRecords(seqFileOut, ids, seqs);
   writeRecords(seqFileOut, ids, seqs, quals);

Assignment 6
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Change the result of Assignment 5 to store the data for the two records in :dox:`StringSet StringSets` and write them out using :dox:`SequenceStream#writeAll`.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/seq_io/solution6.cpp

Next Steps
----------

* Read the Wikipedia articles about the `FASTA file format <http://en.wikipedia.org/wiki/FASTA_format>`_ and the `FASTQ file format and quality values <http://en.wikipedia.org/wiki/FASTQ_format>`_ to refresh your knowledge.
* Read the :ref:`tutorial-indexed-fasta-io` tutorial to learn how to read FASTA files efficiently in a random-access fashion.
* Continue with the :ref:`tutorial`.
