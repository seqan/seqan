.. sidebar:: ToC

   .. contents::


.. _tutorial-parsing:

Parsing
=======

Learning Objective
  This tutorial will explain how to use the class :dox:`RecordReader` for the parsing of text file fo
  You will see in-detail examples for parsing TSV-based formats such as GFF and BLAST tabular output and also for parsing the recursive Newick tree format.

Difficulty
  Advanced

Duration
  40 min

Prerequisites
  :ref:`tutorial-input-output-overview`, :ref:`tutorial-lexical-casting`

In this tutorial, you will learn how to use the :dox:`RecordReader` functions to easily create parsers for structured text file formats.
We will first give a quick example for parsing a simple TSV format.
Then, single-pass parsing will be explained (which is the most important variant) and the interface of the :dox:`RecordReader` class and the ``skip*()`` and ``read*()`` functions will be described.
This is followed by extensive examples on parsing the GFF and BLAST tabular format and an example on how to parse the non-linear Newick format for phylogenetic trees.
The tutorial closes with an explanation of how to write double-pass I/O code and in which situations it is useful to do so.

A First Example
---------------

Let us start off with a quick first example.
The following program reads a two-column TSV file from the standard input.
The first column contains keys, the second one contains values.
The program prints the data as ``${key} -> ${value}`` to stdout.

.. includefrags:: extras/demos/tutorial/parsing/example1.cpp

As you can see, using the :dox:`RecordReader` is straightforward.
First, we construct the :dox:`RecordReader` to wrap ``std::cin`` as also described in the :ref:`tutorial-input-output-overview` tutorial.

Each iteration of the loop loads one record/line from standard input and writes out the record.
We use ``atEnd()`` to check whether we are at the end of the file and loop.
The function ``readUntilChar()`` reads the characters from the underlying file into a buffer ``key`` until a given character occurs, here the character is ``'\t'``.
The reader will not copy the tabulator into ``key`` and stop on the character.
The function ``goNext()`` can be used to go to the next character in the current file.
The call to the function ``readLine()`` copies the data into ``value`` until the end of line, skipping the end-of-line marker (``'\n'`` or ``'\r\n'``) and does not copy the end-of-line marker to the ``value``.
Finally, we print the key/value pair to stdout.

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Reading CSV instead of TSV.

   Type
     Review

   Objective
     Modify the example above to use a comma (``','``) instead of a tab character for separating columns.

   Hint
     Yes, it is very easy.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/parsing/solution1.cpp

The Single-Pass ``RecordReader`` Class Interface
------------------------------------------------

Single-pass record readers can simply be seen and used as an abstraction of streams.
Read the file character-wise, from beginning to the end.

The low-level API for the single-pass reader is as follows:

+-------------------------------------------+-------------------------------------------------------------------------------------+
| **Function**                              | **Description**                                                                     |
+===========================================+=====================================================================================+
| :dox:`RecordReader#atEnd`                 | Return ``true`` if the reader is at the end of the file, ``false`` otherwise.       |
+-------------------------------------------+-------------------------------------------------------------------------------------+
| :dox:`RecordReader#goNext`                | Advance reader in file, return ``true`` if at end of file, ``false`` otherwise.     |
+-------------------------------------------+-------------------------------------------------------------------------------------+
| :dox:`RecordReader#value`                 | Return the character the reader points to at the moment.                            |
+-------------------------------------------+-------------------------------------------------------------------------------------+
| :dox:`RecordReader#resultCode resutlCode` | Return ``int`` with I/O status. 0 for no error, non-0 value for error when reading. |
+-------------------------------------------+-------------------------------------------------------------------------------------+

The following program shows another example of single-pass I/O.
We read a text file line-by-line and append the results to a :dox:`String` of :dox:`CharString CharStrings`.

.. includefrags:: extras/demos/tutorial/parsing/reader_single_demo.cpp

Character Classes and the ``read*`` and ``skip*`` Functions
-----------------------------------------------------------

Character Classes And ``is*``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In SeqAn, the same character classes are used as in the POSIX standard.
See `this list of character classes <http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/ctype.h.html>`_ for a comprehensive list and description.

For example:

.. code-block:: cpp

   printf("isdigit('a') == %d\n", isdigit('a'));  // => "isdigit('a') == 0"
   printf("isdigit('0') == %d\n", isdigit('0'));  // => "isdigit('0') == 1"
   printf("isblank(' ') == %d\n", isdigit(' '));  // => "isdigit(' ') == 0"

The ``read*`` And ``skip*`` Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parsing functionality in SeqAn built on top of the :dox:`StreamConcept` concept and :dox:`RecordReader` class is optimized for reading bioinformatics text file formats.

These formats mostly consist of fairly flat data files, i.e. a sequence of records, each having very few levels of subrecords.
A typical example are FASTQ files where one record consists of adjacent lines, containing the identifier, sequence, and qualities.
Another example are TSV (tab-separated-values) files where each record spans a line and there possibly is a header.
SAM is an example for a TSV file with a header at the top of the file.

The main challenge in reading bioinformatics files is their size.
When parsing a word processor document file, a HTML document, or a computer program, the input file is typically not larger than some MB.
In bioinformatics, files having multiple GB are not uncommon, e.g. NGS data or the sequence of the human genome.

Thus, in SeqAn, the files are parsed "on the fly" as they are read.
Using compiler nomenclauture, bioinformatics parsers often only have to be `tokenizers <http://en.wikipedia.org/wiki/Tokenizing>`_.
Making writing such simple parsers easy is the main aim of the ``read*`` and ``skip*`` functions in SeqAn.
NB: By using :dox:`CharArrayStream Char Array Streams`, you can also use the parsing infrastructure on in-memory data.

For each considered class of characters, there often is a read and a skip function.
There are two big types of classes: White-listing/inclusion (``read*X*``) of certain characters and black-listing/exclusion (``readUntil*X*``) of certain characters.
The inclusion functions stop after the last read/included character, the exclusion functions stop on the first excluded/not read character.

Most functions have the following interface.
Note that all functions only **append** to the ``buffer`` argument, so you have to call :dox:`StringSet#clear` yourself.
This facilitates optimized reading into :dox:`ConcatDirectStringSet Concat Direct StringSets`.

.. code-block:: cpp

   int readUntilXXX (TBuffer & buffer, RecordReader<TStream, TPass> & reader);
   int readXXX      (TBuffer & buffer, RecordReader<TStream, TPass> & reader);
   int skipUntilXXX (RecordReader<TStream, TPass> & reader);
   int skipXXX      (RecordReader<TStream, TPass> & reader);

.. tip::

    I/O Return Values and EOF_BEFORE_SUCCESS

    The ``read*()`` and ``skip*()`` functions return an ``int`` value.
    Consistent with C return codes, the return value is ``== 0`` in case that the reading/skipping was successful and ``!= 0`` if reading/skipping was not successful.

    The cases of unsuccessful reading/skipping include real errors (e.g. hardware problems) but also that the reader is at the end of the file.
    In this case ``seqan::EOF_BEFORE_SUCCESS`` is returned.
    This behaviour is required for file format guessing where a return value of ``seqan::EOF_BEFORE_SUCCESS`` is interpreted as success.

    There are three cases in how code can handle the value ``seqan::EOF_BEFORE_SUCCESS``: (1) interpret it as an error, (2) return ``seqan::EOF_BEFORE_SUCCESS`` itself, or (3) interpret it as "success".

    Here are some examples:

    '''(1) Interpret as Error'''

    Naively, one would assume that this is the correct treatment.
    However, (2) is the right choice for most cases.

    .. code-block:: cpp

       // TRecordReader reader created above.
       seqan::CharString buffer;
       while (atEnd(reader))
       {
           if (readLine(buffer, read) != 0)
               return 1;  // handle as error
       }

**(2) Interpret as ``seqan::EOF_BEFORE_SUCCESS``**

Returning this code gives the caller the opportunity to handle end-of-file different from any other error.
For example, a file format guesser can try to parse the first thousand bytes of a file and see whether they parse as valid.
When ``EOF_BEFORE_SUCCESS`` is returned, it would count this as an access.
Any other non-0 return code would be an error.

.. code-block:: cpp

   // TRecordReader reader created above.
   seqan::CharString buffer;
   int res = 0;
   while (atEnd(reader))
   {
       if ((res = readLine(buffer, read)) != 0)
           return res;  // handle as error or EOF_BEFORE_SUCCESS
   }

**(3) Interpret as Success**

In some cases, EOF is a valid event.
For example, if you have a line-based file format such as TSV, the last line could end with an EOF instead of a line break.

.. code-block:: cpp

   // TRecordReader reader created above.
   seqan::CharString buffer;
   int res = 0;
   while (atEnd(reader))
   {
       if ((res = readLine(buffer, read)) != 0 &&
           res != seqan::EOF_BEFORE_SUCCESS)
           return res;  // line not reached in case of EOF
   }

The following functions are available:

+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| **Function**                                                      | **Description**                                                               |
+===================================================================+===============================================================================+
| :dox:`FileFormatTokenization#readDigits`                          | Read digit characters.                                                        |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#readDna5IgnoringWhitespaces`         | Read DNA 5 characters, ignore whitespace.                                     |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#readLetters`                         | Read letter characters.                                                       |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#readLine`                            | Read whole line, line break is not written into buffer.                       |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#readLineStripTrailingBlanks`         | Read whole line, trailing blanks are not written into buffer.                 |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#readNChars`                          | Read a fixed number of characters.                                            |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#readNCharsIgnoringWhitespace`        | Read a fixed number of characters, whitespace is not written into the buffer. |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#readUntilBlank`                      | Read until a blank character occurs.                                          |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#readUntilChar`                       | Read until the given character occurs.                                        |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#readUntilWhitespace`                 | Read until a whitespace character occurs.                                     |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipBlanks`                          | Skip blank characters.                                                        |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipChar`                            | Skip one given character.                                                     |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipLine`                            | Skip from the current position to the end of the line.                        |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipNChars`                          | Skip a fixed number of characters.                                            |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipNCharsIgnoringWhitespace`        | Skip a fixed number of characters, ignore whitespace.                         |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipUntilBlank`                      | Skip until a blank character occurs.                                          |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipUntilChar`                       | Skip until a certain character occurs                                         |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipUntilGraph`                      | Skip until a graph character occurs.                                          |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipUntilLineBeginsWithChar`         | Skip until a line begins with a certain character.                            |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipUntilLineBeginsWithOneCharOfStr` | Skip until a line begins with one character of a given string/list.           |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipUntilLineBeginsWithStr`          | Skip until a line begins with a certain string.                               |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipUntilString`                     | Skip until a certain string is found.                                         |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipUntilWhitespace`                 | Skip until a whitespace character is found.                                   |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+
| :dox:`FileFormatTokenization#skipWhitespaces`                     | Skip until a non-whitespace character is found.                               |
+-------------------------------------------------------------------+-------------------------------------------------------------------------------+

In the following example, we read the first two fields of a TSV file from stdin and dump them to stdout.

.. code-block:: cpp

   seqan::RecordReader<std::istream, seqan::SinglePass<> > reader(std::cin);
   seqan::CharString buffer;

   while (atEnd(reader))
   {
       clear(buffer);
       int res = readUntilChar(buffer, reader, '\t');
       if (res != 0)
           return res;
       std::cout << buffer;

       if (goNext(reader))
           return seqan::EOF_BEFORE_SUCCESS;

       clear(buffer);
       res = readUntilChar(buffer, reader, '\t');
       if (res != 0)
           return res;
       std::cout << buffer << std::endl;

       res = skipLine(reader);
       if (res != 0 && res != seqan::EOF_BEFORE_SUCCESS)
           return 1;
   }

Writing Your Own ``read*`` and ``skip*`` Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Writing your own reading/skipping function is easy.
As an example, we write functions for reading and skipping the characters from the set *{x, y, z}*.
The functions follow the same pattern and use the functions ``_readHelper()`` and ``_skipHelper()``.

These functions read/skip characters as long as a specific overload of the predicate function ``_charCompare()`` (in the ``seqan`` namespace) returns ``true``.
The ``_charCompare()`` function gets two parameters: The character to test and a tag for selecting the specific ``_charCompare()`` overload.
The caracter to test is of type ``int``.
The tag is defined by you as a developer and the tag given to ``_charCompare()`` is the same as given to ``_readHelper()`` and ``_skipHelper()``.

For good examples, you can look at the file ``core/include/seqan/stream/tokenize.h`` to see how the rest of the ``read*`` and ``skip*`` functions from above are implemented.

.. code-block:: cpp

   struct Xyz_;
   typedef seqan::Tag<Xyz_> Xyz;

   inline int
   _charCompare(int const c, Xyz const & /* tag*/)
   {
       return c == 'x' || c == 'y' || c == 'z';
   }

   template <typename TStream, typename TPass, typename TBuffer>
   inline int
   readXyz(TBuffer & buffer, seqan::RecordReader<TStream, TPass> & reader)
   {
       return seqan::_readHelper(buffer, reader, Xyz(), false);
   }

   template <typename TBuffer, typename TStream, typename TPass>
   inline int
   readUntilXyz(TBuffer & buffer, seqan::RecordReader<TStream, TPass> & reader)
   {
       return seqan::_readHelper(buffer, reader, Xyz(), true);
   }

   template <typename TStream, typename TPass>
   inline int
   skipXyz(seqan::RecordReader<TStream, TPass> & reader)
   {
       return seqan::_skipHelper(reader, Xyz(), false);
   }

   template <typename TStream, typename TPass>
   inline int
   skipUntilXyz(seqan::RecordReader<TStream, TPass> & reader)
   {
       return seqan::_skipHelper(reader, Xyz(), true);
   }

Assignment 2
""""""""""""

.. container:: assignment

   Writing ``readHexNumber()``.

   Type
     Review

   Objective
     Write your own read and skip routines for hexadecimal numbers.
     Such numbers can only contain digits ``0-9`` and the characters ``a-f`` and ``A-F``.

   Solution
     .. container:: foldable

        The following program reads from stdin as long as the input forms a valid hexadecimal number.
        Note that you can send an end-of-file character to your application by pressing ``Ctrl + d``.

        .. includefrags:: extras/demos/tutorial/parsing/solution2.cpp

        An example session.
        The ``Ctrl + d`` is shown as ``^D``.

        .. code-block:: console

           # tutorial_parsing_solution2
           foo
           10
           20
           2a^D
           RECOGNIZED f
           RECOGNIZED 10
           RECOGNIZED 20
           RECOGNIZED 2a

Assignment 3
""""""""""""

.. container:: assignment

   Writing ``readPunctuation()``.

   Type
     Review

   Objective
     Modify the example above to read a sequence of punctuation characters in a function called ``readPunctuation()``.

   Hint
     You can use the function ``ispunct()``.

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/parsing/solution3.cpp

        An example session.
        The ``Ctrl + d`` is shown as ``^D``.

        .. code-block:: console

           ...
           asdf
           !!@#%%^
           RECOGNIZED ...
           RECOGNIZED !!
           RECOGNIZED !!@#%%^

File Parsing Practice
---------------------

This section will walk you through a parser for GFF, tabular BLAST output, and the Newick tree format.

Common Patterns
^^^^^^^^^^^^^^^

In order to support a new file format, you usually (1) introduce a ``struct`` type for storing records, (2) create tags for the file type and the records, and (3) provide overloads of the functions ``nextIs()`` and ``readRecord()``.
For example, for the GFF format, we

* create a ``struct GffRecord`` (1)
* create the tag ``Gff`` (2)
* create overloads of ``nextIs`` and ``readRecord`` for ``Gff`` (3).

A Simple GFF2 Parser
^^^^^^^^^^^^^^^^^^^^

We will implement a simple parser for the `GFF file format version 2 <http://www.sanger.ac.uk/resources/software/gff/spec.html>`_.
For the sake of simplicity, will not implement parsing of ``##`` and will read the whole *attributes* field as one and not subdivide it further.
Here, GFF2 files are TSV files with the following fields.

::

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

The following example shows a GFF2 parser.
First, include the necessary headers.

.. includefrags:: extras/demos/tutorial/parsing/parse_gff2.cpp
   :fragment: includes

Then, define ``Gff2`` tag and record struct.

.. includefrags:: extras/demos/tutorial/parsing/parse_gff2.cpp
   :fragment: tags-structs

We then implement a parser function for GFF records.
Note that most of the code is error handling.

.. includefrags:: extras/demos/tutorial/parsing/parse_gff2.cpp
   :fragment: read-record

On top of the record-reading routine, we implement reading of whole documents.
This is quite simple.

.. includefrags:: extras/demos/tutorial/parsing/parse_gff2.cpp
   :fragment: read-batch

Finally, some driver code to open a file and call the parser routine.
In the end, we dump some of the information we just read.

.. includefrags:: extras/demos/tutorial/parsing/parse_gff2.cpp
   :fragment: main

Let's look at an example run of the program.

.. code-block:: console

    # cat extras/demos/tutorial/parsing /gff2_example.txt
    IV     curated  mRNA   5506800 5508917 . + .   Transcript B0273.1; Note "Zn-Finger"
    IV     curated  5'UTR  5506800 5508999 . + .   Transcript B0273.1
    IV     curated  exon   5506900 5506996 . + .   Transcript B0273.1
    IV     curated  exon   5506026 5506382 . + .   Transcript B0273.1
    IV     curated  exon   5506558 5506660 . + .   Transcript B0273.1
    IV     curated  exon   5506738 5506852 . + .   Transcript B0273.1
    IV     curated  3'UTR  5506852 5508917 . + .   Transcript B0273.1
    # ./extras/demos/tutorial/parsing/tutorial_parse_gff2 extras/demos/tutorial/parsing/gff2_example.txt
    IV  +   0   5508917
    IV  +   0   5508999
    IV  +   0   5506996
    IV  +   0   5506382
    IV  +   0   5506660
    IV  +   0   5506852
    IV  +   0   5508917

Newick Tree Parsing (Recursion Example)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The newick tree format is used for encoding phylogenetic trees (see `Newick Tree Format Standard <http://evolution.genetics.washington.edu/phylip/newick_doc.html>`_ for a formal specification).
We will write a parser that reads Newick forest files (without allowing for comments).

Here is an example for the Newick format:

::

    (((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;

A file with this content encodes the following tree:

::

               +-+ One
            +--+
            |  +--+ Two
         +--+
         |  | +----+ Three
         |  +-+
         |    +--+ Four
         +
         +------+ Five

And here is the grammar of the Newick format in EBNF.

::

    forest        = tree+;
    tree          = node, ";";
    node          = children, label?, distance?
                  | children?, label, distance?;
    children      = "(", node, (",",node)*, ")";
    label         = quoted-list
                  | unquoted-list;
    distance      = ":", number;
    quoted-list   = "'", (qchar escaped-quote)*, "'";
    escaped-quote = "''";
    unquoted-list = uqchar;

The following demo shows the parsers, code to dump the tree from the internal data structures and a small driver program for the routines.

First, the necessary includes.

.. includefrags:: extras/demos/tutorial/parsing/parse_newick.cpp
   :fragment: includes

Then, we define a ``Newick`` tag and a struct for branch labels.

.. includefrags:: extras/demos/tutorial/parsing/parse_newick.cpp
   :fragment: tags-structs

In a next step, we write a ``readFloatLiteral()`` helper function that is reusable.

.. includefrags:: extras/demos/tutorial/parsing/parse_newick.cpp
   :fragment: read-float

The code for reading a Newick forest is recursive and a bit lengthy but not too complex.
We load such forests into strings of :dox:`Tree` objects.
Additionally, we have a vertex map for the branch distances and the vertex labels for each tree.

.. includefrags:: extras/demos/tutorial/parsing/parse_newick.cpp
   :fragment: reading

The code for dumping a Newick forest is also quite simple, if lengthy because of error checks.

.. includefrags:: extras/demos/tutorial/parsing/parse_newick.cpp
   :fragment: writing

Finally, the ``main()`` routine.

.. includefrags:: extras/demos/tutorial/parsing/parse_newick.cpp
   :fragment: main

Let's look at an example run.
Note that the children in SeqAn trees do not have a specific order and the Newick format does not introduce any normalized order.
In the written result, the order of the children has changed.

.. code-block:: console

    # cat extras/demos/tutorial/parsing/newick_example.txt
    (a,('Darwin''s Bulldog (Huxley)',c):-1.92e19)'The ''Root''':5;
    ((a_node,
      'another node',
      bird:0.3134)higher_node:4.5,
     c):1.03e10;
    ((<sub>),(,(</sub>,),));
    # tutorial_parse_newick extras/demos/tutorial/parsing/newick_example.txt
    ((c,'Darwin''s Bulldog (Huxley)'):-1.92e+19,a)'The ''Root''':5;
    (c,(bird:0.3134,'another node',a_node)higher_node:4.5):1.03e+10;
    ((,(<sub>,),),(</sub>));

Parsing Tabular BLAST
^^^^^^^^^^^^^^^^^^^^^

The program *BLASTN* can be given an ``-outfmt`` parameter that makes it generate tabular output.
This output is quite easy to parse (much easier than the human-readable BLAST reports) and looks as follows:

.. code-block:: console

    # blastn -subject NC_001405.fasta -query NC_001460.fasta -outfmt 7 > blast_example.txt
    # cat blast_example.txt
    # BLASTN 2.2.25+
    # Query: gi|9626621|ref|NC_001460.1| Human adenovirus A, complete genome
    # Subject: gi|9626158|ref|NC_001405.1| Human adenovirus C, complete genome
    # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    # 3 hits found
    gi|9626621|ref|NC_001460.1| gi|9626158|ref|NC_001405.1| 81.13   408 66  11  17730   18131   18827   19229   5e-87    316
    gi|9626621|ref|NC_001460.1| gi|9626158|ref|NC_001405.1| 81.63   98  12  6   383 476 433 528 9e-15   76.8
    gi|9626621|ref|NC_001460.1| gi|9626158|ref|NC_001405.1| 76.27   118 22  6   25147   25261   26644   26758   3e-09   58.4
    # BLAST processed 1 queries

The following example program takes the name of such a blastn output, reads it into record data structures and then prints it out in a different format again.
To do this, we will first implement a record-reading API that allows streaming through the file.
Then, we build a batch-reading API that reads such a file into a sequence of records that are all kept in main memory.

The program starts with including the required headers.

.. includefrags:: extras/demos/tutorial/parsing/parse_blastn.cpp
   :fragment: includes

Then, we define a record for the file format ``BlastnTab`` and tabs for the comment and alignment record types.

.. includefrags:: extras/demos/tutorial/parsing/parse_blastn.cpp
   :fragment: tags

Next, we define a record type.
Note that this record type is very specialized to the ``blastn -outfmt 7`` format.
When writing I/O code for multiple format for similar data, you might want to consider writing one record type for all of them.
See the (upcoming, TODO) SAM record I/O for the implementation of one record type for the SAM and then BAM format.

We also create a simple function to dump the record to a stream.

.. includefrags:: extras/demos/tutorial/parsing/parse_blastn.cpp
   :fragment: record

Then, we define :dox:`RecordReader#nextIs` functions for the ``BlastnTabComment`` and ``BlastnTabAlignment`` tags, and their represented record types.

.. includefrags:: extras/demos/tutorial/parsing/parse_blastn.cpp
   :fragment: next-is

Then, we implement a record-reading API on top of the ``skip*`` and ``read*`` functions.
Note that the error handling bloats up the number of required lines but is necessary.

.. includefrags:: extras/demos/tutorial/parsing/parse_blastn.cpp
   :fragment: read-record

On top of the record-reading API, we implement a batch-reading function.
This function turns out to be fairly simple.

.. includefrags:: extras/demos/tutorial/parsing/parse_blastn.cpp
   :fragment: batch-read

In the ``main()`` routine, we can then simply open a ``std::fstream``, create a :dox:`RecordReader`.
Then, use the batch-reading API to read the whole file into main memory and write it to *stdout* again.

.. includefrags:: extras/demos/tutorial/parsing/parse_blastn.cpp
   :fragment: main

The program's output looks as follow:

.. code-block:: console

    # ./extras/demos/tutorial/parsing/tutorial_parse_blastn ../../extras/demos/tutorial/parsing/blast_example.txt
    query name: gi|9626621|ref|NC_001460.1|
    subject name: gi|9626158|ref|NC_001405.1|
    identity: 81.13
    alignment length: 408
    mismatches: 66
    gap opens: 11
    query begin: 17730
    query end: 18131
    subject begin: 18827
    subject end: 19229
    evalue: 5e-87
    bit score: 316

    query name: gi|9626621|ref|NC_001460.1|
    subject name: gi|9626158|ref|NC_001405.1|
    identity: 81.63
    alignment length: 98
    mismatches: 12
    gap opens: 6
    query begin: 383
    query end: 476
    subject begin: 433
    subject end: 528
    evalue: 9e-15
    bit score: 76.8

    query name: gi|9626621|ref|NC_001460.1|
    subject name: gi|9626158|ref|NC_001405.1|
    identity: 76.27
    alignment length: 118
    mismatches: 22
    gap opens: 6
    query begin: 25147
    query end: 25261
    subject begin: 26644
    subject end: 26758
    evalue: 3e-09
    bit score: 58.4

Double-Pass I/O Using the ``RecordReader``
------------------------------------------

The :dox:`DoublePassRecordReader Double-Pass RecordReader` reader's API extends the function described above for the :dox:`SinglePassRecordReader Single-Pass RecordReader`.
It provides the following additional global interface functions.

+-----------------------------------------------+------------------------------+
| **Function**                                  | **Description**              |
+===============================================+==============================+
| :dox:`DoublePassRecordReader#startFirstPass`  | Start first pass of reading. |
+-----------------------------------------------+------------------------------+
| :dox:`DoublePassRecordReader#startSecondPass` | Second pass of reading.      |
+-----------------------------------------------+------------------------------+

It is used as follows: For each section of the file that is to be read in the next step (one or multiple records), you first call :dox:`DoublePassRecordReader#startFirstPass`.
This memoizes the current position in the file.
Then, you use the same API as for the single-pass reader to read the file.
When you are done with this section, you call :dox:`DoublePassRecordReader#startSecondPass`.
This will reset the position of the reader to the one where :dox:`DoublePassRecordReader#startFirstPass` was called.

Here is an example for using double-pass I/O:

.. includefrags:: extras/demos/tutorial/parsing/reader_double_demo.cpp

Note that all file contents read in the first pass are buffered when operating on streams.
Thus, double-pass I/O can have a high memory usage on streams when having large passes.
In this case, using memory mapped strings to read from can be more efficient.
However, in order to allow double-pass I/O when reading from compressed streams or stdin, this buffering is designed to lead to better performance or is even required.

Double-pass I/O has the advantage that the exact amount of memory can be allocated for the target data structures.
This can lead to reduced memory usage since no memory is pre-allocated and then left unused.
Thus, this is useful if the life span of your target data structures is long and a lot of memory is saved.

The disadvantage is the higher memory usage when reading the file itself.
All data read in the first pass has to be buffered if using streams.

So, **when should you use double-pass I/O?** A good **rule of thumb** is: *If you need to read a whole large file into main memory (e.g. NGS read set or a genome) and it is uncompressed then use a double-pass record reader with a memory mapped string. Otherwise, use single-pass I/O.*
