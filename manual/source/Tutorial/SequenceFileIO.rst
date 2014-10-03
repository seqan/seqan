.. sidebar:: ToC

   .. contents::


.. _tutorial-sequence-file-io:

Sequence File I/O
=================

Learning Objective
  This tutorial will explain how to read and write sequence files in SeqAn using the :dox:`RecordReader` class and the :dox:`MultiSeqFile` API.

Difficulty
  Advanced

Duration
  30 min

Prerequisites
  :ref:`tutorial-input-output-overview`

This tutorial explains how to read using the :dox:`RecordReader` class and the :dox:`MultiSeqFile MultipleSeqFile` API and how to write sequence files using :dox:`StreamConcept Streams`.
Currently, SeqAn supports reading sequences (and qualities where it applies) in FASTA, FASTQ, EMBL and GenBank format.

This tutorial does not explain how to use the :dox:`SequenceStream` class.
This class is covered in the :ref:`tutorial-basic-sequence-io` tutorial.

Reading Sequence Files
----------------------

Because sequences play a central role in the SeqAn, the library provides a comprehensive implementation of input routines.
The sequence I/O facilities are available through the ``seq_io`` module.
You get access to this module by ``#include <seqan/seq_io.h>``.
We will first describe the API for record-wise reading of sequence files and then go over to reading all records in a file at once.

We start out by creating a :dox:`RecordReader` object.
For most cases, the single-pass specialization is the most appropriate one.
Then, we can read the file record-by-record using the function :dox:`RecordReader#readRecord`.

To this function, we pass the buffers we want to read the record identifier and the sequence into (``id`` and ``seq``).
This is followed by the ``reader`` object we create with the ``std::fstream`` previously opened.
The last argument is the tag ``seqan::Fasta()`` which uses the variant of :dox:`RecordReader#readRecord` for reading FASTA files.

.. includefrags:: core/demos/tutorial/sequence_io/example1.cpp

When reading FASTQ files, we use the tag ``seqan::Fastq()``.
For storing the qualities, we can pass an optional third parameter of type :dox:`CharString` which stores the quality values.

.. includefrags:: core/demos/tutorial/sequence_io/example2.cpp

Optionally, we can also read the sequence into a string of [dox:Dna5Q
Dna5Q] characters which will store the qualities directly in the
string's characters.

.. includefrags:: core/demos/tutorial/sequence_io/example3.cpp

.. important::

    Sequence Parsing Behaviour

    * When using :dox:`Dna5` or :dox:`Dna5Q` as the sequence's alphabet type, the parsing routine will allow the characters ``'C'``, ``'G'``, ``'A'``, ``'T'``, and ``'N'`` in the sequences of the file.
      This can make problems if the sequenc contains different characters, for example when it contains IUPAC characters.
      In this case, you can simply use :dox:`CharString` as the ``seq`` parameter and then assign them to a :dox:`Dna5String`.
    * Accordingly, when using :dox:`Dna` or :dox:`DnaQ`, only the characters ``'C'``, ``'G'``, ``'A'``, and ``'T'`` are allowed.
    * When omitting the ``qual`` parameter when reading FASTQ, the quality values from the file will be ignored.

Assignment 1
""""""""""""

.. container:: assignment

   Record-Wise Reading Sequences into :dox:`CharString`

   Type
     Review

   Objective
     Modify the example above to read the sequence into a :dox:`CharString` instead of a :dox:`Dna5String`.

   Solution
     .. container:: foldable

        .. includefrags:: core/demos/tutorial/sequence_io/solution1.cpp

When we want to read a whole sequence (e.g. FASTA or FASTQ) file into memory then we only have to slightly adjust the example from above.
For example, here is how we can read a whole FASTQ file into memory using the function :dox:`RecordReader#read` into :dox:`StringSet StringSets` of :dox:`CharString CharStrings` and :dox:`Dna5String Dna5Strings`.

.. warning::

   For a short time, ``read()`` will still be called ``read2()`` because of name clashes with the old I/O system.

.. includefrags:: core/demos/tutorial/sequence_io/example4.cpp

Assignment 2
""""""""""""

.. container:: assignment

   Document-Wise Reading Sequences into :dox:`CharString`

   Type
     Review

   Objective
     Modify the example above to read the sequence into a :dox:`StringSet` of :dox:`CharString CharStrings` instead of a :dox:`Dna5String Dna5Strings`.

   Solution
     .. container:: foldable

        .. includefrags:: core/demos/tutorial/sequence_io/solution2.cpp

Choice of Record Reader
-----------------------

In most cases, you will want to use a :dox:`SinglePassRecordReader Single-Pass RecordReader` for reading files.
Mostly, it is the fastest and best way to read files and also all file formats have a single-pass implementation.

Using a double-pass record reader almost only makes sense if read a whole file into main memory using the document reading API.
The file is read twice.
In the first pass, the total length of ids and sequence characters is determined.
When reading sequences into :dox:`StringSet StringSets`, the exact number of elements can be reserved.
Even more, when using :dox:`ConcatDirectStringSet Concat-Direct StringSet`, no superflous memory has to be allocated at all.
The string sets are then filled in the second pass.

Using double-pass I/O also only makes sense for document reading when used in conjunction with :dox:`MMapString MMap Strings`.
When using streams, the :dox:`RecordReader` has to buffer the read data in memory because not all stream implementation allow for jumping.
In the case of :dox:`MMapString MMap Strings`, no buffer is used because the record reader directly operates on the memory mapped file (and thus directly on the disk buffers of the kernel).

Assignment 3
""""""""""""

.. container:: assignment

   Using a :dox:`DoublePassRecordReader Double-Pass RecordReader` with a :dox:`MMapString MMap String`.

   Type
     Application

   Objective
     Change solution of Assignment 2 such that a :dox:`DoublePassRecordReader Double-Pass RecordReader` is used with a :dox:`MMapString MMap String`.

   Hint
     You can open files into MMap Strings as follows (include the ``<seqan/file.h>`` header):

      .. code-block:: cpp

         typedef seqan::String<char, seqan::MMap<> > TMMapString;
         TMMapString mmapString;
         bool success = open(mmapString, "filename.fa", seqan::OPEN_RDONLY);


     You can then define a :dox:`DoublePassRecordReader` wrapping the just opened ``mmapString`` as follows:

     .. code-block:: cpp

        typedef seqan::RecordReader<
                TMMapString,
                seqan::DoublePass<seqan::StringReader> > TReader;
        TReader reader(mmapString);


   Solution
     .. container:: foldable

        .. includefrags:: core/demos/tutorial/sequence_io/solution3.cpp

Auto Format Detection
---------------------

Passing the format as the tag is appropriate when the format is known beforehand.
Otherwise, you can use a variable of type :dox:`AutoSeqStreamFormat` instead of the tag.

:dox:`AutoSeqStreamFormat`\ t objects can be first passed to the function :dox:`guessStreamFormat`.
This function tries to parse the file as different formats on the first some thousand bytes.
When this succeeds, the successfully recognized file type is stored in the object.

You can then subsequently use the :dox:`AutoSeqStreamFormat` instead of a tag to the functions :dox:`RecordReader#readRecord` or :dox:`RecordReader#read`.

.. includefrags:: core/demos/tutorial/sequence_io/example9.cpp

Assignment 4
""""""""""""

.. container:: assignment

   Using :dox:`AutoSeqStreamFormat`

   Type
     Application

   Objective
     Adjust the solution of Assignment 3 to use a :dox:`AutoSeqStreamFormat` for format detection.

   Solution
     .. container:: foldable

        .. includefrags:: core/demos/tutorial/sequence_io/solution6.cpp

.. note::

    Qualities and FASTA files.

    When passing a ``qual`` parameter to :dox:`RecordReader#readRecord` or :dox:`RecordReader#read` then this cannot be filled with qualities from the file since FASTA files do not contain any.
    Instead, the ``qual`` string will be empty after the call to :dox:`RecordReader#readRecord` and after the call to :dox:`RecordReader#read`, it will be a string set with empty entries.
    The string set will have a size that is equal to the number of records in the file.

Writing Sequence Files
----------------------

Similar to reading, sequence files can be written record-by-record or as a whole.

For record-wise writing, we use the function :dox:`StreamConcept#writeRecord`.
This function expects as parameters, the :dox:`StreamConcept` to write to, the data to write, followed by the format tag.
The following example writes an identifier and a sequence :dox:`StringSet` record-by-record to stdout.

.. includefrags:: core/demos/tutorial/sequence_io/example6.cpp

The result on the console looks like this:

.. code-block:: console

    >id1
    CGATCGATCGAT
    >id2
    AAAAAAAAAAAA

Assignment 5
""""""""""""

.. container:: assignment

   Writing out FASTQ.

   Type
     Application

   Objective
     Change the example above such that the two sequences are written as FASTQ with qualities.
     Use the quality strings ``"IIIIIIIIIHII"`` and ``"IIIIIIIIIIII"``.

   Hint
     Simply use a new :dox:`StringSet` ``quals`` of :dox:`CharString`, append the quality strings, and modify the line with the ``writeRecord()`` call.

   Solution
     .. container:: foldable

        .. includefrags:: core/demos/tutorial/sequence_io/solution5.cpp

        The output looks as follows:

        .. code-block:: console

            @id1
            CGATCGATCGAT
            +
            IIIIIIIIIHII
            @id2
            AAAAAAAAAAAA
            +
            IIIIIIIIIIII

For writing out whole string sets at once, we use the function :dox:`StreamConcept#write2 write`.
The transition from record-wise writing to writing whole string sets is of similar simplicity as for reading:

.. warning::

   For a short time, ``write()`` will still be called ``write2()`` because of name clashes with the old I/O system.

.. includefrags:: core/demos/tutorial/sequence_io/example8.cpp

Using MultiSeqFile
------------------

.. warning::

   Deprecate ``MultiSeqFile`` in favour of ``FaiIndex``?

The class :dox:`MultiSeqFile` (which actually is a shortcut to a memory mapped string set) allows to read sequence files in a two-pass approach.
First, the file is read and the start positions of each sequence record in the file is stored in memory.
The file is kept open as a memory mapped file.

Then, we can access the identifier, sequence, and quality string of a record using functions such as :dox:`assignSeqId`.

Indexed reading can be done through :dox:`MultiSeqFile` which is a shortcut to a memory mapped string set.
We open the file using :dox:`File#open` on its ``concat`` member (which is a :dox:`MMapString MMap String`).
The function :dox:`split` then parses the file contents and sets the separating indexes of the :dox:`StringSet`.
For this, we need the file format. We could give a specify format in the tag (e.g. ``seqan::Fastq()``) or use :dox:`AutoSeqFormat` together with :dox:`guessFormat`.

The following example demonstrates how to use :dox:`MultiSeqFile` to read sequence files.
First, we include the necessary headers and start our ``main()`` function.

.. includefrags:: core/demos/tutorial/sequence_io/example5.cpp
   :fragment: includes-main

Then, we declare the :dox:`MultiSeqFile` object and open it with the value of ``argv[1]``.
If no parameters are given then we exit the program with status code ``1``.

.. includefrags:: core/demos/tutorial/sequence_io/example5.cpp
   :fragment: open

This is followed by using :dox:`AutoSeqFormat` for guessing the sequence file type.

.. includefrags:: core/demos/tutorial/sequence_io/example5.cpp
   :fragment: guess

After guessing the file type, we can now use this knowledge to compute the start positions of each record using the function :dox:`split`.

.. includefrags:: core/demos/tutorial/sequence_io/example5.cpp
   :fragment: split

After the call to :dox:`split`, we can get the number of sequences in the file using the function :dox:`ContainerConcept#length`.
We declare the :dox:`StringSet StringSets` for storing the sequences and sequence ids and reserve the exact space for the number of elements we need.

.. includefrags:: core/demos/tutorial/sequence_io/example5.cpp
   :fragment: load

Then, we declare some buffers for storing the sequence id, characters, and the quality values.

.. includefrags:: core/demos/tutorial/sequence_io/example5.cpp
   :fragment: buffers

Now, we can access the sequence, qualities and ids using the functions :dox:`assignSeq`, :dox:`assignQual`, and :dox:`assignSeqId`.
Note that these functions still have to do some parsing of the input file.
The number of sequences is the same as the number of entries in the ``MultiSeqFile`` ``StringSet`` as returned by :dox:`ContainerConcept#length`.

In the following loop, we first extract the sequences, qualities, and the sequence id.
Then, the qualities are stored in the :dox:`Dna5Q` letters of the string.
The sequence with qualities and the sequence ids are then stored in the variables ``seqs`` and ``seqIDs`` we allocated above.

.. includefrags:: core/demos/tutorial/sequence_io/example5.cpp
   :fragment: output

Finally, we return the status code ``0`` at the end of our ``main()`` function.

.. includefrags:: core/demos/tutorial/sequence_io/example5.cpp
   :fragment: return

Indexed reading has multiple advantages.

*  Its performance is only slightly worse than when reading sequentially
   with a double-pass String RecordReader.
*  The input file is mapped into main memory and otherwise complicated page-wise memory management is done by the operating system and does not have to be implemented by the user.
   The user can access the file almost at random and only the used parts will be loaded into main memory.
   This is quite efficient when only few sequences are needed.

If you need to have fast random access to all sequences in a file then loading it into a :dox:`ConcatDirectStringSet Concat-Direct StringSet` with the batch-reading API is faster than using :dox:`MultiSeqFile`.

.. container:: assignment

   MultiSeqFile Review

   Type
     Review

   Objective
      Change the example above, so the sequence file that is read is written to the user in a TSV format.
      For each record in the input file with id ``${ID}``, sequence ``${SEQ}``, and quality string ``${QUAL}``, write out a line ``${ID}\t${SEQ}\t${QUAL}``.

   Solution
     .. container:: foldable

        .. code-block:: cpp

           #include <seqan/file.h>
           #include <iostream>

           int main (int argc, char const ** argv)
           {
               seqan::MultiSeqFile multiSeqFile;
               if (argc < 2 || !open(multiSeqFile.concat, argv[1], seqan::OPEN_RDONLY))
                   return 1;

               seqan::AutoSeqFormat format;
               guessFormat(multiSeqFile.concat, format);
               split(multiSeqFile, format);

               seqan::String<seqan::Dna5> seq;
               seqan::CharString qual;
               seqan::CharString id;

               for (unsigned i = 0; i < seqCount; ++i)
               {
                   assignSeq(seq, multiSeqFile[i], format);    // read sequence
                   assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
                   assignSeqId(id, multiSeqFile[i], format);   // read sequence id

                   std::cout << id << '\t' << seq << '\t' << qual << '\n';
               }

               return 0;
           }

Next Steps
----------

* Read the Wikipedia articles about the `FASTA file format <http://en.wikipedia.org/wiki/FASTA_format>`_ and the `FASTQ file format and quality values <http://en.wikipedia.org/wiki/FASTQ_format>`_ to refresh your knowledge.
* Read the :ref:`tutorial-basic-sequence-io` tutorial to learn how to use the :dox:`SequenceStream` class.
* Read the :ref:`tutorial-indexed-fasta-io` tutorial tutorial to learn how to read FASTA files efficiently in a random-access fashion.
* Continue with :ref:`tutorial`.
