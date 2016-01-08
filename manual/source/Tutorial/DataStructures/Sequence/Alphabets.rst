.. sidebar:: ToC

    .. contents::

.. _tutorial-datastructures-sequences-alphabets:

Alphabets
=========

Learning Objective
  You will learn the details about the alphabets in SeqAn.

Difficulty
  Basic

Duration
  15 min

Prerequisites
  :ref:`tutorial-getting-started-first-steps-in-seqan`

This tutorial will describe the different alphabets used in SeqAn, or in other words, you will learn about the contained types of a SeqAn :dox:`String`.
To continue with the other tutorials, it would be enough to know, that in SeqAn several standard alphabets are already predefined, e.g. :dox:`Dna`, :dox:`Dna5`, :dox:`Rna`, :dox:`Rna5`, :dox:`Iupac`, :dox:`AminoAcid`.

Types
-----

Any type that provides a default constructor, a copy constructor and an assignment operator can be used as the alphabet / contained type of a :dox:`String` (see also the tutorial :ref:`tutorial-datastructures-sequences`).
This includes the C++ `POD types <http://www.parashift.com/c++-faq-lite/intrinsic-types.html#faq-26.7>`_, e.g. ``char``, ``int``, ``double`` etc.
In addition you can use more complex types like :dox:`String` as the contained type of strings, e.g. ``String<String<char> >``.

SeqAn also provides the following types that are useful in bioinformatics.
Each of them is a specialization of the class :dox:`SimpleType`.

+------------------+-------------------------------------------------------------+
| Specialization   | Description                                                 |
+==================+=============================================================+
| :dox:`AminoAcid` | Amino Acid Alphabet                                         |
+------------------+-------------------------------------------------------------+
| :dox:`Dna`       | DNA alphabet                                                |
+------------------+-------------------------------------------------------------+
| :dox:`Dna5`      | ``N`` alphabet including ``N`` character                    |
+------------------+-------------------------------------------------------------+
| :dox:`DnaQ`      | ``N`` alphabet plus phred quality                           |
+------------------+-------------------------------------------------------------+
| :dox:`Dna5Q`     | ``N`` alphabet plus phred quality including ``N`` character |
+------------------+-------------------------------------------------------------+
| :dox:`Finite`    | Finite alphabet of fixed size.                              |
+------------------+-------------------------------------------------------------+
| :dox:`Iupac`     | ``N`` Iupac code.                                           |
+------------------+-------------------------------------------------------------+
| :dox:`Rna`       | ``N`` alphabet                                              |
+------------------+-------------------------------------------------------------+
| :dox:`Rna5`      | ``N`` alphabet including ``N`` character                    |
+------------------+-------------------------------------------------------------+

Functionality
-------------

In SeqAn, alphabets are value types that can take a limited number of values and which hence can be mapped to a range of natural numbers.
We can retrieve the number of different values of an alphabet, the alphabet size, by the metafunction :dox:`FiniteOrderedAlphabetConcept#ValueSize`.

.. includefrags:: demos/tutorial/alphabets/example_size.cpp
    :fragment: main

.. includefrags:: demos/tutorial/alphabets/example_size.cpp.stdout

Another useful metafunction called :dox:`AlphabetConcept#BitsPerValue` can be used to determine the number of bits needed to store a value of a given alphabet.

.. includefrags:: demos/tutorial/alphabets/example_bitsPerValue.cpp
    :fragment: main

.. includefrags:: demos/tutorial/alphabets/example_bitsPerValue.cpp.stdout

The order of a character in the alphabet (i.e. its corresponding natural number) can be retrieved by calling the function :dox:`FiniteOrderedAlphabetConcept#ordValue`.
See each specialization's documentation for the ordering of the alphabet's values.

.. includefrags:: demos/tutorial/alphabets/example_ordValue.cpp
    :fragment: main

.. includefrags:: demos/tutorial/alphabets/example_ordValue.cpp.stdout

.. tip::

    The return value of the :dox:`FiniteOrderedAlphabetConcept#ordValue` function is determined by the metafunction :dox:`FiniteOrderedAlphabetConcept#ValueSize`.
    :dox:`FiniteOrderedAlphabetConcept#ValueSize` returns the type which uses the least amount of memory while being able to represent all possible values.
    E.g. :dox:`FiniteOrderedAlphabetConcept#ValueSize` of :dox:`Dna` returns an ``_uint8`` which is able to represent 256 different characters.
    However, note that ``std::cout`` has no visible symbol for printing all values on the screen, hence a cast to ``unsigned`` might be necessary.

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     In this task you will learn how to access all the letters of an alphabet.
     Use the piece of code from below and adjust the function ``showAllLettersOfMyAlphabet()`` to go through all the characters of the current alphabet and print them.

     .. includefrags:: demos/tutorial/alphabets/assignment_1.cpp

   Hints
     You will need the Metafunction :dox:`FiniteOrderedAlphabetConcept#ValueSize`.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/alphabets/assignment_1_solution.cpp
