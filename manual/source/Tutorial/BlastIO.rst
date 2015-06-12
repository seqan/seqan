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

Other recommended reading:
  `Basics of Blast Statistics <http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html>`_


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

+-----------------------------+------------------+----------------+--------------------+
| Format                      | read support     |  write support | based on version   |
+=============================+==================+================+====================+
| pairwise                    |                  |  ✓             |  Blast-2.2.26+     |
+-----------------------------+------------------+----------------+--------------------+
| tabular                     |  ✓               |  ✓             |  Blast-2.2.26+     |
+-----------------------------+------------------+----------------+--------------------+
| tabular (legacy)            |  ✓               |  ✓             |  Blast-2.2.26      |
+-----------------------------+------------------+----------------+--------------------+
| tabular w comments          |  ✓               |  ✓             |  Blast-2.2.26+     |
+-----------------------------+------------------+----------------+--------------------+
| tabular w comments (legacy) |  ✓               |  ✓             |  Blast-2.2.26      |
+-----------------------------+------------------+----------------+--------------------+

*Please note the missing + in the last column. There are still software releases for both Blast generations, in our
case we chose the same version number of each, but this does not imply that it is the same program.*

Tabular formats
---------------

TODO

