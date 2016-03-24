.. sidebar:: ToC

    .. contents::

.. _how-to-recipes-work-with-custom-score-matrices:

Working With Custom Score Matrices
==================================

This How To describes how to create new scoring matrices for Amino Acids and DNA alphabets and how to load score matrices from files.

Creating A New Built-In Score Matrix
------------------------------------

The following program demonstrates how to implement a new built-in score matrix.

.. includefrags:: demos/howto/scores/init_score.cpp
   :fragment: includes

We need to perform the necessary definitions for the matrix.
This consists of two steps:

#. Defining a tag struct.
#. Specializing the class ``ScoringMatrixData_`` with your tag.

Note how we use enum values to compute the matrix size which itself is retrieved from the :dox:`FiniteOrderedAlphabetConcept#ValueSize` metafunction.

.. includefrags:: demos/howto/scores/init_score.cpp
   :fragment: user-defined-matrix

Now we define a function ``showScoringMatrix`` for displaying a matrix.

.. includefrags:: demos/howto/scores/init_score.cpp
   :fragment: show-scoring-matrix

Finally, the ``main`` function demostrates some of the things you can do with scores:

* Construct an empty score matrix object (2.)
* Fill the score matrix in a loop (3.1)
* Fill the matrix with the user-defined matrix values (3.2)
* Directly create a score matrix with the user-defined matrix values (4)

.. includefrags:: demos/howto/scores/init_score.cpp
   :fragment: main

Here is the output of the program:

.. includefrags:: demos/howto/scores/init_score.cpp.stdout

Loading Score Matrices From File
------------------------------------

This small demo program shows how to load a score matrix from a file.
Examples for a score file are ``demos/howto/scores/dna_example.txt`` for DNA alphabets and ``tests/sPAM250`` for amino acids.

Include the necessary headers.

.. includefrags:: demos/howto/scores/load_score.cpp
   :fragment: includes

We define a function that can show a scoring matrix.

.. includefrags:: demos/howto/scores/load_score.cpp
   :fragment: show-scoring-matrix

Finally, the main program loads the scoring matrix and then shows it.

.. includefrags:: demos/howto/scores/load_score.cpp
   :fragment: main

Here's the program output.

.. includefrags:: demos/howto/scores/load_score.cpp.stdout
