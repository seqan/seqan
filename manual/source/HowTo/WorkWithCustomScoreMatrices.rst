.. sidebar:: ToC

   .. contents::


.. _how-to-work-with-custom-score-matrices:

Working With Custom Score Matrices
==================================

This How To describes how to create new scoring matrices for Amino Acids and DNA alphabets and how to load score matrices from files.

Creating A New Built-In Score Matrix
------------------------------------

The following program demonstrate how to implement a new built-in score matrix.

.. includefrags:: demos/howto/scores/init_score.cpp
   :fragment: includes

Then, we perform the necessary definitions for the matrix.
This consists of three steps:

* defining a tag struct
* specializing the class ``ScoringMatrixData_`` with your tag

Note how we use enum values to compute the matrix size which itself is retrieved from the :dox:`FiniteOrderedAlphabetConcept#ValueSize` metafunction.

.. includefrags:: demos/howto/scores/init_score.cpp
   :fragment: user-defined-matrix

We define a function ``showScoringMatrix`` for displaying a matrix.

.. includefrags:: demos/howto/scores/init_score.cpp
   :fragment: show-scoring-matrix

Finally, the function ``main`` function demostrates some of the things you can do with scores:

* Construct empty score matrix object (2.)
* Programmatically fill the score matrix in a loop (3.1)
* Programatically fill the matrix with the user-defined matrix values (3.2)
* Directly create a score matrix with the user-defined matrix values (4)

.. includefrags:: demos/howto/scores/init_score.cpp
   :fragment: main

Here is the output of the program:

.. code-block:: console

    $ make tutorial_init_score
    $ ./demos/tutorial_init_score
    Coordinate Products
    	A	C	G	T	N
    A	0	0	0	0	0
    C	0	1	2	3	4
    G	0	2	4	6	8
    T	0	3	6	9	12
    N	0	4	8	12	16
    User defined matrix (also Dna5 scoring matrix)...
    	A	C	G	T	N
    A	1	0	0	0	0
    C	0	1	0	0	0
    G	0	0	1	0	0
    T	0	0	0	1	0
    N	0	0	0	0	0
    User DNA scoring scheme...
    	A	C	G	T	N
    A	1	0	0	0	0
    C	0	1	0	0	0
    G	0	0	1	0	0
    T	0	0	0	1	0
    N	0	0	0	0	0

Loading Score Matrices From File
------------------------------------

This small demo program shows how to load a score matrix from a file.
Examples for score file are ``demos/howto/scores/dna_example.txt`` for DNA alphabets and ``tests/sPAM250`` for amino acids.

Include the necessary headers.

.. includefrags:: demos/howto/scores/load_score.cpp
   :fragment: includes

We define a function that can show a scoring matrix.

.. includefrags:: demos/howto/scores/load_score.cpp
   :fragment: show-scoring-matrix

Finally, the main program loads the scoring matrix from the file given on the command line and then shows it.

.. includefrags:: demos/howto/scores/load_score.cpp
   :fragment: main

Here's the program output.

.. code-block:: console

   $ make tutorial_load_score
   $ ./demos/tutorial_load_score ../../demos/howto/scores/dna_example.txt
       A   C   G   T
   A   1   -1  -1  -1
   C   -1  1   -1  -1
   G   -1  -1  1   -1
   T   -1  -1  -1  1

