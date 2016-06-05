.. sidebar:: ToC

    .. contents::

.. _tutorial-io-blast-io:

Blast I/O
=========

Learning Objective
  In this tutorial, you will learn about different the Blast file formats and how to interact with them in SeqAn.

Difficulty
  Average

Duration
  1h30min - 2h30min

Prerequisite Tutorials
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-io-input-output-overview`, :ref:`tutorial-datastructures-alignment`, :ref:`tutorial-algorithms-alignment-pairwise-sequence-alignment`

Other recommended reading
  `Basics of Blast Statistics <http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html>`_

Technical requirements
  Full C++11 support required in compiler
  (GCC :math:`\ge` 4.9, Clang :math:`\ge` 3.4 or MSVC :math:`\ge` 2015)


The Basic local alignment search tool (BLAST) by the NCBI is one of the most widely used tools in Bioinformatics.
It supports a variety of formats which, although widely used, are not standardized and partly even declared as
"subject to unannounced and undocumented change". This makes implementing them very hard and is one of the reasons
dealing with Blast IO is more difficult than with other file formats.

Furthermore it is important to distinguish between the Blast version written in C (known by its `blastall` executable)
and the C++ version (individual `blastn`, `blastp`... executables), also known as BLAST+.
The formats and their identifiers changed between versions.
The following table gives an overview:

+---------------------------------------------------+------------------+-------------+
| Description                                       | `blastall`       |  `blast*`   |
+===================================================+==================+=============+
| pairwise                                          |  -m 0            |  -outfmt 0  |
+---------------------------------------------------+------------------+-------------+
| query-anchored showing identities                 |  -m 1            |  -outfmt 1  |
+---------------------------------------------------+------------------+-------------+
| query-anchored no identities                      |  -m 2            |  -outfmt 2  |
+---------------------------------------------------+------------------+-------------+
| flat query-anchored, show identities              |  -m 3            |  -outfmt 3  |
+---------------------------------------------------+------------------+-------------+
| flat query-anchored, no identities                |  -m 4            |  -outfmt 4  |
+---------------------------------------------------+------------------+-------------+
| query-anchored no identities and blunt ends       |  -m 5            |             |
+---------------------------------------------------+------------------+-------------+
| flat query-anchored, no identities and blunt ends |  -m 6            |             |
+---------------------------------------------------+------------------+-------------+
| XML Blast output                                  |  -m 7            |  -outfmt 5  |
+---------------------------------------------------+------------------+-------------+
| tabular                                           |  -m 8            |  -outfmt 6  |
+---------------------------------------------------+------------------+-------------+
| tabular with comment lines                        |  -m 9            |  -outfmt 7  |
+---------------------------------------------------+------------------+-------------+
| Text ASN.1                                        |  -m 10           |  -outfmt 8  |
+---------------------------------------------------+------------------+-------------+
| Binary ASN.1                                      |  -m 11           |  -outfmt 9  |
+---------------------------------------------------+------------------+-------------+
| Comma-separated values                            |                  |  -outfmt 10 |
+---------------------------------------------------+------------------+-------------+
| BLAST archive format (ASN.1)                      |                  |  -outfmt 11 |
+---------------------------------------------------+------------------+-------------+

The files written by `blastall` are considered the
"legacy"-format in SeqAn. SeqAn has support for the following formats:

+-----------------------------+----------------+----------------+--------------------+
| Format                      | read support   |  write support | based on version   |
+=============================+================+================+====================+
| pairwise                    |                |       ✓        |  Blast-2.2.26+     |
+-----------------------------+----------------+----------------+--------------------+
| tabular                     |      ✓         |       ✓        |  Blast-2.2.26+     |
+-----------------------------+----------------+----------------+--------------------+
| tabular (legacy)            |      ✓         |       ✓        |  Blast-2.2.26      |
+-----------------------------+----------------+----------------+--------------------+
| tabular w comments          |      ✓         |       ✓        |  Blast-2.2.26+     |
+-----------------------------+----------------+----------------+--------------------+
| tabular w comments (legacy) |      ✓         |       ✓        |  Blast-2.2.26      |
+-----------------------------+----------------+----------------+--------------------+

.. caution::

    Please note that *Blast-2.2.26+* is **not the same** as *Blast-2.2.26*! One is version 2.2.26 of the C++
    application suite (BLAST+) and the other is version 2.2.26 of the legacy application suite. There still are software
    releases for both generations.

There are different *program modes* in Blast which also influence the file format. In the legacy application suite
these where specified with the ``-p`` parameter, in the BLAST+ suite they each have their own executable.

+-----------------------------+------------------+--------------------+
| Program mode                | query alphabet   |  subject alphabet  |
+=============================+==================+====================+
| BlastN                      | nucleotide       |  nucleotide        |
+-----------------------------+------------------+--------------------+
| BlastP                      | protein          |  protein           |
+-----------------------------+------------------+--------------------+
| BlastX                      | translated nucl. |  protein           |
+-----------------------------+------------------+--------------------+
| TBlastN                     | protein          |  translated nucl.  |
+-----------------------------+------------------+--------------------+
| TBlastX                     | translated nucl. |  translated nucl.  |
+-----------------------------+------------------+--------------------+

Tabular formats
---------------

The tabular formats are tab-seperated-value formats (TSV), with twelve columns by default.
Each line represents one match (or *high scoring pair* in Blast terminology).
The twelve default columns are:

    1. Query sequence ID (truncated at first whitespace)
    2. Subject sequence ID (truncated at first whitespace)
    3. Percentage of identical positions
    4. Alignment length
    5. Number of mismatches
    6. Number of gap openings
    7. Start position of alignment on query sequence
    8. End position of alignment on query sequence
    9. Start position of alignment on subject sequence
    10. End position of alignment on subject sequence
    11. Expect value (length normalized bit score)
    12. Bit score (statistical significance indicator)

.. note:: Alignment positions in Blast

   #. **Interval notation:** Blast uses 1-based closed intervals for positions, i.e. a match from the 100th position to
      the 200th position of a sequence will be shown as ``100  200`` in the file. SeqAn internally uses 0-based
      half open intervals, i.e. it starts counting at position 0 and stores the first position behind the sequence
      as "end", e.g. position ``99`` and ``200`` for our example.
   #. **Reverse strands:** For matches found on the reverse complement strand the positions are counted backwards from
      the end of the sequence, e.g. a match from the 100th position to the 200th position on a reverse complement strand
      of a sequence of length 500 will be shown as ``400 300`` in the file.
   #. **Translation frames:** Positions given in the file are always on the original untranslated sequence!

   The ``writeRecord()`` function automatically does all of these conversions!


A **tabular** file could look like this (matches per query are sorted by e-value):

.. literalinclude:: ../../../../tests/blast/defaultfields.m8

The **tabular with comment lines** format additionally prefixes every block belonging to one query sequence with
comment lines that include the program version, the database name and column labels. The above example would look
like this:

.. literalinclude:: ../../../../tests/blast/defaultfields.m9
    :lines: 5-37

As you can see, comment lines are also printed for query sequences which don't have any matches.
The major
difference of these formats in BLAST+ vs the legacy application are that the *mismatches* column used to
include the number of gap characters, but it does not in BLAST+. The comments also look slightly different
in the **tabular with comment lines (legacy)** format:

.. literalinclude:: ../../../../tests/blast/defaultfields_legacy.m9
    :lines: 5-35



Pairwise format
---------------

The pairwise format is the default format in Blast. It is more verbose than the tabular formats and very human readable.

This is what the last record from above would look like (the other queries are omitted):

.. literalinclude:: ../../../../demos/tutorial/blast_io/plus_sub.m0


Blast formats in SeqAn
----------------------

There are three blast format related tags in SeqAn:

  #. :dox:`BlastReport` for the pairwise format with the :dox:`FormattedFile` output specialization
     :dox:`BlastReportFileOut`.
  #. :dox:`BlastTabular` for the tabular formats with the :dox:`FormattedFile` output and input specializations
     :dox:`BlastTabularFileOut` and :dox:`BlastTabularFileIn`.
  #. :dox:`BlastTabularLL` which provides light-weight, but very basic tabular IO.

The third tag can be used for simple file manipulation, e.g. filtering or column rearrangement, it is not covered in
this tutorial (see the dox for a simple example).

To work with the first two formats you need to understand at least the following data structures:
  * :dox:`BlastRecord`: the record covers all :dox:`BlastMatch` es belonging to one query sequence.
  * :dox:`FormattedFile`: one of :dox:`BlastReportFileOut`, :dox:`BlastTabularFileOut` and :dox:`BlastTabularFileIn`.
  * :dox:`BlastIOContext`: the context of the FormattedFile.

The context contains file-global data like the name of the database and can also be used to read/write certain file
format properties, e.g. "with comment lines" or "legacyFormat".

.. caution::
    Due to the structure of blast tabular files lots of information is repeated in every block of comment lines, e.g.
    the database name. Because it is expected that these stay the same they are saved in the context and not the record.
    You may still, however, check every time you ``readRecord()`` if you want to make sure.

File reading example
--------------------

Only tabular formats are covered in this example, because no input support is available for the pairwise format.

Copy the contents of the **tabular with comment lines** example above into a file and give it to the following
program as the only parameter. Please use ``.m9`` as file type extension.

.. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
    :language: c++
    :lines: 1-18, 74-84

Assignment 1
""""""""""""

.. container:: assignment

  Objective
    Complete the above example by reading the file according to :dox:`BlastTabularFileIn`.
    For every record print the query ID, the number of contained matches and the bit-score of the best match.

  Solution
      Top
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 1-18

      New code
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 19-34, 57-60

      Bottom
          .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 74-84

Assignment 2
""""""""""""

.. container:: assignment

  Objective
    Study the documentation of :dox:`BlastIOContext`. How can you adapt the previous program to check if there were any
    problems reading a record? If you have come up with a solution, try to read the file at
    ``tests/blast/defaultfields.m9``. What does the program print and why?

  Solution
      Top
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 1-34

      New code
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 41-60

            The program will print conformancyErrors for the last record, because there is a typo in the file
            ( ``Datacase`` instead of ``Database`` ).

      Bottom
          .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 74-84

Assignment 3
""""""""""""

.. container:: assignment

  Objective
    Now that you have a basic understanding of :dox:`BlastIOContext`, also print the following information after
    reading the records:

      * file format (with comment lines or without, BLAST+ or legacy?)
      * blast program and version
      * name of database

    Verify that the results are as expected on the files ``tests/blast/defaultfields.m8``,
    ``tests/blast/defaultfields.m9`` and ``tests/blast/defaultfields_legacy.m9``.

  Solution
      Top
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 1-34, 41-60

      New code
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 61-73

      Bottom
          .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 74-84

Assignment 4
""""""""""""

As previously mentioned, twelve columns are printed by default.
This can be changed in BLAST+, also by means of the ``--outfmt`` parameter.
A standards compliant **file with comment lines** and custom column composition can be read
without further configuration in SeqAn.

.. tip::
    Don't believe it? Look at ``tests/blast/customfields.m9``, as as you can see the bit score is in the 13th column
    (instead of the twelfth). If you run your program on this file, it should still print the correct bit-scores!

.. container:: assignment

  Objective
    Read :dox:`BlastIOContext` again focusing on :dox:`BlastIOContext::fields` and also read :dox:`BlastMatchField`.
    Now adapt the previous program to print for every record the ``optionLabel`` of each field used.

    Verify that the results are as expected on the files ``tests/blast/defaultfields.m9`` and
    ``tests/blast/customfields.m9``.

  Solution
      Top
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 1-34

      New code
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 35-40

      Bottom
          .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/read_assignment.cpp
                :language: c++
                :lines: 41-84

      If this was too easy, you can also try the same for tabular files without comment lines!

File writing example
--------------------

The following program stub creates three query sequences and two subject sequences in amino acid alphabet.
We will later generate records with matches and print these to disk.

.. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
    :language: c++
    :lines: 1-49, 105-114

Assignment 5
""""""""""""

.. container:: assignment

  Objective
    Before we can begin to align, certain properties of the :dox:`BlastIOContext` have to be set.
    Which ones? Add a block with the necessary initialization to the above code, use blast defaults where possible.
    Although not strictly required at this point: include a call to ``writeHeader()``.

    .. caution::
        Alignment score computation works slightly different in Blast and in SeqAn, please have a look at
        :dox:`BlastScoringScheme`. For the task at hand it should suffice to simply use the corresponding
        ``set*()`` functions on :dox:`BlastIOContext::scoringScheme`.

  Solution
      Top
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 1-49

      New code
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 50-64

      Bottom
          .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 106-116

Assignment 6
""""""""""""

.. container:: assignment

  Objective
        Next create a record for every query sequence, and in each record a match for every query-subject pair.
        Compute the local alignment for each of those matches. Use the align member of the match object.

  Solution
      Top
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 1-64

      New Code
       .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 65-82, 97-99, 103

      Bottom
       .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 106-116

Assignment 7
""""""""""""

.. container:: assignment

  Objective
        Now that you have the align member computed for every match, also save the begin and end positions, as well
        as the lengths. Blast Output needs to now about the number of gaps, mismatches... of every match, how can
        they be computed?

  Solution
      Top
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 1-82

      New Code
       .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 83-92

      Bottom
       .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 97-99, 103, 106-116

Assignment 8
""""""""""""

.. container:: assignment

  Objective
        Finally add e-value statistics and print the results to a file:
         * compute the bit score and e-value for every match
         * discard matches with an e-value greater than 1
         * for every record, sort the matches by bit-score
         * write each record
         * write the footer before exiting the program

  Solution
      Top
       .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 1-92

      New Code
       .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 93-96

      Bottom
       .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.cpp
                :language: c++
                :lines: 97-116

      Your output file should look like this:
        .. container:: foldable

            .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.m9

Assignment 9
""""""""""""

.. container:: assignment

  Objective
    Up until now you have only printed the **tabular with comment lines** format. What do you have to do to print
    without comment lines? In legacy format? What about the pairwise format?

  Solution
    Tabular without comment lines
      .. container:: foldable

        Add ``context(outfile).tabularSpec = BlastTabularSpec::NO_COMMENTS`` at l.53.
        Remember to use ``.m8`` as file extension!

        The result should look like this:

        .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.m8

    Tabular without comment lines (legacy)
      .. container:: foldable

        To print in legacy tabular format (with or without comment lines), add ``context(outfile).legacyFormat = true`` at l.53.

        The result should look like this(legacy and NO_COMMENTS):

        .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment_legacy.m8

    Pairwise format
      .. container:: foldable

        To print in the pairwise format replace l.52 with ``BlastReportFileOut<TContext> outfile(argv[1]);``.

        Remember to use ``.m0`` as file extension!

        The result should look like this:

        .. literalinclude:: ../../../../demos/tutorial/blast_io/write_assignment.m0
