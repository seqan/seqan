.. sidebar:: ToC

    .. contents::

.. _tutorial-io-sequence-io:

Sequence I/O
============

Learning Objective
  You will learn how to read and write sequence files in FASTA, FASTQ, EMBL or GenBank format.

Difficulty
  Basic

Duration
  20 min

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-io-input-output-overview`

This tutorial explains how to read and write sequence files using the :dox:`SeqFileIn` and :dox:`SeqFileOut` classes.
These classes provide an API for accessing sequence files in different file formats, either compressed or uncompressed.


FASTA/FASTQ Format
------------------

FASTA/FASTQ are record-based files.
A FASTA record contains the sequence id and the sequence characters.
Here is an example of FASTA file:

.. includefrags:: demos/tutorial/sequence_io/example.fa

In addition to that, a FASTQ record contains also a quality value for each sequence character.
Here is an example of FASTQ file:


.. includefrags:: demos/tutorial/sequence_io/example.fq


SeqFile Formats
---------------

We can read sequence files with the :dox:`SeqFileIn` class and write them with the :dox:`SeqFileOut` class.
These classes support files in FASTA, FASTQ, EMBL or GenBank format.

Note that :dox:`SeqFileOut` will guess the format from the file name.

+--------------+----------------------+
| File Format  | File Extension       |
+==============+======================+
| FASTA        | ``.fa``, ``.fasta``  |
+--------------+----------------------+
| FASTQ        | ``.fq``, ``.fastq``  |
+--------------+----------------------+
| EMBL         | ``.embl``            |
+--------------+----------------------+
| GenBank      | ``.gbk``             |
+--------------+----------------------+

A First Working Example
-----------------------

Let us start out with a minimal working example.
The following program reads a FASTA file called ``example.fa`` and prints out the identifier and the sequence of the first record.

.. includefrags:: demos/tutorial/sequence_io/example1.cpp

We call the :dox:`FormattedFile#FormattedFile SeqFileIn constructor` with the path to the file to read.
Successively, we call the function :dox:`SeqFileIn#readRecord` to read the first record from the file.
Note that, differently from all others :dox:`FormattedFileIn` classes, :dox:`SeqFileIn#readRecord` accepts **separate** identifier and sequence :dox:`String Strings` rather than one single record object.

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Reproduction

   Objective
     Copy the above example of a FASTA file in a new file ``example.fa`` in a directory of your choice.

     Copy the program above into a new application ``basic_seq_io_example``, adjust the path ``"example.fa"`` to the just created FASTA file, compile the program, and run it.

     You should see the following output:

     .. includefrags:: demos/tutorial/sequence_io/example1.cpp.stdout

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/sequence_io/solution1.cpp


Handling Errors
---------------

As explained in the :ref:`tutorial-io-input-output-overview` tutorial, :dox:`SeqFileIn` and :dox:`SeqFileOut` throw exceptions to signal eventual errors.
Invalid characters inside an input file will be signaled by :dox:`SeqFileIn#readRecord` via parsing exceptions.

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Improve the above program to handle errors.

   Hint
     You can use the generic class :dox:`Exception` to catch both low-level and high-level I/O errors.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/sequence_io/solution2.cpp


Accessing Records in Batches
----------------------------

There are three use cases for reading or writing record-based files:

#. read or write the file **record by record**;
#. read or write a **batch of records**, e.g. 100k records at a time;
#. read or write **all records** from or to the file.

The class :dox:`SeqFileIn` provides the functions :dox:`SeqFileIn#readRecord` and :dox:`SeqFileIn#readRecords`, while the class :dox:`SeqFileOut` provides the functions :dox:`SeqFileOut#writeRecord` and :dox:`SeqFileOut#writeRecords`.

.. tip::

    Reading records in batches is more efficient than reading single records.


Note that the function :dox:`SeqFileIn#readRecords` uses :dox:`StringSet` instead of :dox:`String`.
By default, :dox:`SeqFileIn#readRecords` reads **all** remaining records.
Optionally, one can specify a batch of records to be read.

.. includefrags:: demos/tutorial/sequence_io/base.cpp
      :fragment: batch


Assignment 3
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Change your program from above to load all sequences and print them in the same fashion.

     You should be able to run your program on the example file we created above and see the following output:

     .. includefrags:: demos/tutorial/sequence_io/solution3.cpp.stdout

   Hint
     You can use the function :dox:`SeqFileIn#readRecords` to load all records at once.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/sequence_io/solution3.cpp


Accessing Qualities
-------------------

Functions :dox:`SeqFileIn#readRecord`, :dox:`SeqFileIn#readRecords`, :dox:`SeqFileOut#writeRecord` and :dox:`SeqFileOut#writeRecords` are available in two variants:

#. the first variant accepts only the sequence identifier and sequence characters, besides the :dox:`SeqFileIn` object;
#. the second variant accepts an additional :dox:`CharString` for a PHRED base quality string.

If the first variant is used on an output file containing qualities, e.g. a FASTQ file, then :dox:`SeqFileOut#writeRecord` writes qualities as ``'I'``, i.e. PHRED score 40.
If the second variant is used on an input file containing no qualities, e.g. a FASTA file, then :dox:`SeqFileIn#readRecord` returns **empty** quality strings.

Here is an example for the second variant of :dox:`SeqFileIn#readRecord`:

.. includefrags:: demos/tutorial/sequence_io/base.cpp
      :fragment: qual

.. tip::

    When :dox:`DnaQ` or :dox:`Dna5Q` :dox:`String Strings` are used, then you should use the second variant.
    The qualities are simply stored directly in the sequence characters.


Assignment 4
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Copy the above example of FASTQ file in a new file ``example.fq`` in a directory of your choice.

     Change your result of Assignment 3 to use the variant of :dox:`SeqFileIn#readRecord` that also reads in the qualities and writes them next to the sequences.

     When your program is called on this file, the result should look as follows.

     .. includefrags:: demos/tutorial/sequence_io/solution4.cpp.stdout

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/sequence_io/solution4.cpp


Next Steps
----------

* Read the Wikipedia articles about the `FASTA file format <http://en.wikipedia.org/wiki/FASTA_format>`_ and the `FASTQ file format and quality values <http://en.wikipedia.org/wiki/FASTQ_format>`_ to refresh your knowledge.
* Read the :ref:`tutorial-io-indexed-fasta-io` tutorial to learn how to read FASTA files efficiently in a random-access fashion.
* Continue with the :ref:`tutorial`.
