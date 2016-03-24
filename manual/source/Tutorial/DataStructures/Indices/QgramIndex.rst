.. sidebar:: ToC

    .. contents::

.. _tutorial-datastrucures-indices-q-gram-index:

Q-gram Index
============

Learning Objective
  You will know the features of the q-gram Index, how it can be used for searching and how to access the different fibres.

Difficulty
  Average

Duration
  1 h

Prerequisites
  :ref:`tutorial-datastructures-sequences`

The Q-gram Index
----------------

A q-gram index can be used to efficiently retrieve all occurrences of a certain q-gram in the text.
It consists of various tables, called fibres (see :ref:`how-to-recipes-access-index-fibres`), to retrieve q-gram positions, q-gram counts, etc.
However, it has no support for suffix tree iterators.
A q-gram index must be specialized with a :dox:`Shape` type.
A :dox:`Shape` defines q, the number of characters in a q-gram and possibly gaps between these characters.
There are different specializations of :dox:`Shape` available:

+-----------------------+--------------------+----------------------+
| Specialization        | Modifiable         | Number of Gaps       |
+=======================+====================+======================+
| :dox:`UngappedShape`  | \-                 | 0                    |
+-----------------------+--------------------+----------------------+
| :dox:`SimpleShape`    | \+                 | 0                    |
+-----------------------+--------------------+----------------------+
| :dox:`OneGappedShape` | \+                 | 0/1                  |
+-----------------------+--------------------+----------------------+
| :dox:`GappedShape`    | \-                 | any                  |
+-----------------------+--------------------+----------------------+
| :dox:`GenericShape`   | \+                 | any                  |
+-----------------------+--------------------+----------------------+

* \- *fixed at compile time*, \+ *can be changed at runtime*

Each shape evaluates a gapped or ungapped sequence of q characters to a hash value by the Functions :dox:`Shape#hash`, :dox:`Shape#hashNext`, etc.
For example, the shape ``1101`` represents a 3-gram with one gap of length 1.
This shape overlayed with the :dox:`Dna` text ``"GATTACA"`` at the third position corresponds to ``"TT-C"``.
The function :dox:`Shape#hash` converts this 3-gram into :math:`61 = (\mathbf{3} \cdot 4 + \mathbf{3}) \cdot 4 + 1`.
4 is the alphabet size in this example (see :dox:`FiniteOrderedAlphabetConcept#ValueSize`).

With :dox:`Shape#hash` and :dox:`Shape#hash hashNext`, we can compute the hash values of arbitrary / adjacent q-grams and a loop that outputs the hash values of all overlapping ungapped 3-grams could look as follows:

.. includefrags:: demos/tutorial/q_gram_index/index_qgram_hash.cpp
   :fragment: hash_loop1

Note that the shape not only stores the length and gaps of a q-gram shape but also stores the hash value returned by the last hash/hashNext call.
This hash value can be retrieved by calling :dox:`Shape#value` on the shape.
However, one drawback of the example loop above is that the first hash value must be computed with :dox:`Shape#hash` while the hash values of the following overlapping q-grams can more efficiently be computed by :dox:`Shape#hashNext`.
This complicates the structure of algorithms that need to iterate all hash values, as they have to handle this first hash differently.
As a remedy, the :dox:`Shape#hashInit` function can be used first and then :dox:`Shape#hashNext` on the first and all following text positions in the same way:

.. includefrags:: demos/tutorial/q_gram_index/index_qgram_hash.cpp
   :fragment: hash_loop2

The q-gram index offers different functions to search or count occurrences of q-grams in an indexed text, see :dox:`IndexQGram#getOccurrences`, :dox:`IndexQGram#countOccurrences`.
A q-gram index over a :dox:`StringSet` stores occurrence positions in the same way and in the same fibre (FibreSA) as the ESA index.
If only the number of q-grams per sequence are needed the QGramCounts and QGramCountsDir fibres can be used.
They store pairs ``(seqNo, count)``, ``count``>0, for each q-gram that occurs ``counts`` times in sequence number ``seqNo``.

To efficiently retrieve all occurrence positions or all pairs ``(seqNo, count)`` for a given q-gram, these positions or pairs are stored in contiguous blocks (in QGramSA, QGramCounts fibres), called buckets.
The begin position of bucket i is stored in directory fibres (QGramDir, QGramCountsDir) at position i, the end position is the begin positions of the bucket i+1.
The default implementation of the :dox:`IndexQGram` index maps q-gram hash values 1-to-1 to bucket numbers.
For large q or large alphabets the :dox:`OpenAddressingQGramIndex Open Addressing QGram Index` can be more appropriate as its directories are additionally bound by the text length.
This is realized by a non-trivial mapping from q-gram hashes to bucket numbers that requires an additional fibre (QGramBucketMap).

For more details on q-gram index fibres see :ref:`how-to-recipes-access-index-fibres` or :dox:`QGramIndexFibres QGram Index Fibres`.

Example
-------

We want to construct the q-gram index of the string ``"CATGATTACATA"`` and output the occurrences of the ungapped 3-gram ``"CAT"``.
As 3 is fixed at compile-time and the shape has no gaps we can use an :dox:`UngappedShape` which is the first template argument of :dox:`IndexQGram`, the second template argument of :dox:`Index`.
Next we create the string ``"CATGATTACATA"`` and specialize the first index template argument with the type of this string.
The string can be given to the index constructor.

.. includefrags:: demos/tutorial/q_gram_index/index_qgram.cpp
   :fragment: initialization

To get all occurrences of a q-gram, we first have to hash it with a shape of the same type as the index shape (we can even use the index shape returned by :dox:`IndexQGram#indexShape`).
The hash value returned by :dox:`Shape#hash` or :dox:`Shape#hashNext` is also stored in the shape and is used by the function :dox:`IndexQGram#getOccurrences` to retrieve all occurrences of our 3-gram.

.. includefrags:: demos/tutorial/q_gram_index/index_qgram.cpp
   :fragment: output

Program output:

.. includefrags:: demos/tutorial/q_gram_index/index_qgram.cpp.stdout

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Write a program that outputs all occurrences of the gapped q-gram "AT-A" in "CATGATTACATA".

   Solution
     .. container:: foldable

	Before we can create a :dox:`DnaString` index of "CATGATTACATA", we have to choose an appropriate :dox:`Shape`.
	Because our shape ``1101`` is known at compile-time and contains only one gap we could choose :dox:`OneGappedShape`, :dox:`GappedShape`, or :dox:`GenericShape` (see the commented-out code).
	Although the :dox:`GenericShape` could be used for every possible shape, it is a good idea to choose a :dox:`Shape` with restrictions as its :dox:`Shape#hash` functions are more efficient in general.

	.. includefrags:: demos/tutorial/q_gram_index/index_assignment5.cpp
	   :fragment: initialization

	Please note that the :dox:`Shape` object that corresponds to the :dox:`IndexQGram` index is empty initially and has to be set by :dox:`Shape#stringToShape` or :dox:`Shape#resize`.
	This initialization is not necessary for :dox:`Shape` that are defined at compile-time, i.e. :dox:`UngappedShape` and :dox:`GappedShape`.
	To search for "AT-A" we first have to hash it with the index shape or any other :dox:`Shape` with the same bitmap.
	The we can use :dox:`IndexQGram#getOccurrences` to output all matches.

	.. includefrags:: demos/tutorial/q_gram_index/index_assignment5.cpp
          :fragment: output

	.. tip::

	   Instead of ``length(getOccurrences(...))`` we could have used :dox:`IndexQGram#countOccurrences`.
	   But beware that :dox:`IndexQGram#countOccurrences` requires only the ``QGram_Dir`` fibre, whereas :dox:`IndexQGram#getOccurrences` requires both ``QGram_Dir`` and  ``QGram_SA``, see :ref:`how-to-recipes-access-index-fibres`.
	   Because ``QGram_SA`` can be much more efficiently constructed during the construction of ``QGram_Dir``, ``QGram_Dir`` would be constructed twice.

	Program output:

	.. includefrags:: demos/tutorial/q_gram_index/index_assignment5.cpp.stdout

Assignment 2
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Create and output a matrix M where M(i,j) is the number of common ungapped 5-grams between sequence i and sequence j for 3 random :dox:`Dna` sequences, each not longer than 200 characters.
     Optional: Run the matrix calculation twice, once for an :dox:`IndexQGram` and once for an :dox:`OpenAddressingQGramIndex Open Addressing QGram Index` and output the directory sizes (QGram_Dir, QGram_CountsDir fibre).

   Hint
     A common q-gram that occurs :math:`a` times in one and :math:`b` times in the other sequence counts for :math:`\min(a,b)`.

   Solution
     .. container:: foldable

        For generating random numbers we use the `std::mt19937 <http://www.cplusplus.com/reference/random/mt19937/>`_.
        The random numbers returned by the random number engine are arbitrary ``unsigned int`` values which we downscale to values between 0 and 3 and convert into :dox:`Dna` characters.
        The 3 generated strings are of random length and appended to a :dox:`StringSet`.
        The main algorithm is encapsulated in a template function ``qgramCounting`` to easily switch between the two :dox:`IndexQGram` specializations.

        .. includefrags:: demos/tutorial/q_gram_index/index_assignment6.cpp
           :fragment: initialization

        The main function expects the :dox:`StringSet` and the :dox:`Index` specialization as a tag.
        First, we define lots of types we need to iterate and access the fibres directly.
        We then notify the index about the fibres we require.
        For storing the common q-grams we use a 2-dimensional :dox:`Matrix` object whose lengths have to be set with ``setLength`` for each dimension.
        The matrix is initialized with zeros by :dox:`Matrix#resize`.

        .. includefrags:: demos/tutorial/q_gram_index/index_assignment6.cpp
           :fragment: matrix_init

        The main part of the function iterates over the CountsDir fibre.
        Each entry in this directory represents a q-gram bucket, a contiguous interval in the Counts fibre storing for every sequence the q-gram occurs in the number of occurrences in pairs (seqNo,count).
        The interval begin of each bucket is stored in the directory and the interval end is the begin of the next bucket.
        So the inner loops iterate over all non-empty buckets and two pairs (seqNo1,count1) and (seqNo2,count2) indicate that seqNo1 and seqNo2 have a common q-gram.
        At the end the matrix can simply be output by shifting it to the ``cout`` stream.

        .. includefrags:: demos/tutorial/q_gram_index/index_assignment6.cpp
           :fragment: matrix_calculation

        Please note that the :dox:`OpenAddressingQGramIndex open addressing` q-gram index directories are smaller than the :dox:`IndexQGram` index directories.

        Program output:

        .. includefrags:: demos/tutorial/q_gram_index/index_assignment6.cpp.stdout
