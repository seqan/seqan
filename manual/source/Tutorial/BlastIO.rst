.. sidebar:: ToC

   .. contents::


.. _tutorial-blast-io:

Blast I/O
=========

Learning Objective
  In this tutorial, you will learn about different the Blast file formats and how to interact with them in SeqAn.

Difficulty
  Average

Duration
  999min

Prerequisite Tutorials
  :ref:`tutorial-sequences`, :ref:`tutorial-input-output-overview`, :ref:`tutorial-alignment-representation`, :ref:`tutorial-pairwise-sequence-alignment`

Other recommended reading
  `Basics of Blast Statistics <http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html>`_

Technical requirements
  Full C++11 support required in compiler (GCC>=4.9, Clang>=3.4 or MSVC>=2015)


The Basic local alignment search tool (BLAST) by the NCBI is one of the most widely used tools in Bioinformatics.
It supports a variety of formats which, although widely used, are not standardized and partly even declared as
"subject to unannounced and undocumented change". This makes implementing them very hard and is one of the reasons
dealing with Blast IO is more difficult than with other file formats.

Furthermore it is important to distinguish between the Blast version written in C (known by its `blastall` executable)
and the C++ version (individual `blastn`, `blastp`... executables), also known as BLAST+.
The formats and their identifiers changed
between versions, the following table gives an overview:

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

.. tip::

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

.. tip::

   Alignment positions in Blast.

   #. **Interval notation:** Blast uses 1-based closed intervals for positions, i.e. a match from the 100th position to
      the 200th position of a sequence will be shown as ``100  200`` in the file. SeqAn internally uses 0-based
      half open intervals, i.e. it starts counting at positition 0 and stores the first position behind the sequence
      as "end", e.g. position ``99`` and ``200`` for our example. More on how this is converted, later.
   #. **Reverse strands:** For matches found on the reverse complement strand the positions are counted backwards from
      the end of the sequence, e.g. a match from the 100th position to the 200th position on a reverse complement strand
      of a sequence of length 500 will be shown as ``400 300`` in the file.
   #. **Translation frames:** Positions given in the file are always on the original untranslated sequence!


A **tabular** file could look like this (matches per query are sorted by e-value):

.. literalinclude:: ../../../tests/blast/nocomments_defaults.m8

The **tabular with comment lines** format additionally prefixes every block belonging to one query sequence with
comment lines that include the program version, the database name and column labels. The above example would look
like this:

.. literalinclude:: ../../../tests/blast/plus_comments_defaults.m9
    :lines: 5-37

As you can see, comment lines are also printed for query sequences which don't have any matches.
The major
difference of these formats in BLAST+ vs the legacy application are that the *mismatches* column used to
include the number of gap characters, but it does not in BLAST+. The comments also look slightly different
in the **tabular with comment lines (legacy)** format:

.. literalinclude:: ../../../tests/blast/legacy_comments_defaults.m9
    :lines: 5-35



Pairwise format
---------------

The pairwise format is the default format in Blast, it is more verbose than the tabular formats and very human readable.

This is what the last record from above would look like (the other queries are omitted):

.. literalinclude:: ../../../demos/tutorial/blast/plus_sub.m0


Blast Formats in SeqAn
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
  * :dox:`BlastIOContext`: the context of the FormatteFile.

The context contains file-global data like the name of the database and can also be used to read/write certain file
format properties, e.g. "with comment lines" or "legacyFormat".

File reading
------------

The only

