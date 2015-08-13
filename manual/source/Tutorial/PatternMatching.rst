.. sidebar:: ToC

   .. contents::


.. _tutorial-pattern-matching:

Pattern Matching
================

Learning Objective
  In this tutorial you will learn how to use the SeqAn classes finder and pattern to search a known pattern in a string or :dox:`StringSet`.

Difficulty
  Average

Duration
  40 min

Prerequisites
  :ref:`tutorial-sequences`, :ref:`tutorial-indices`

Pattern matching is about searching a known string or :dox:`StringSet` (``needle``) in another string or :dox:`StringSet` (``haystack``).
This tutorial will introduce you into the SeqAn classes finder and pattern.
It will demonstrate how to use the spezializations of the class finder to perform either an online search or an index based seach.
And you will learn how to specify the search algorithm, which can be exact or approximate.

Overview
--------

In the case of approximate searching errors are allowed, which are either only mismatches or also indels.
Additionally there are filtration algorithms which return potential matches, i.e. ``haystack`` segments possibly containing a pattern match.
All searching is done by calling the function :dox:`Finder#find`, which takes at least two arguments:

#. A :dox:`Finder` that stores all necessary information about the ``haystack`` and the last found position of the ``needle`` within the haystack.
#. A :dox:`Pattern` that stores all information about the ``needle``.
   Some variants of :dox:`Finder#find` support further arguments.
   The :dox:`Finder` and :dox:`Pattern` classes expect the underlying ``haystack`` and ``needle`` types as first template arguments.
   In addition, a second template argument specifies the search algorithm.

Each call of :dox:`Finder#find` finds only one match (or potential match) of the ``needle`` within the haystack.
The :dox:`Finder` can be asked for the begin and end position of the last found match.
The :dox:`Pattern` can be asked for the number of the found sequence if the ``needle`` is a :dox:`StringSet`.
Subsequent calls of find can be used to find more occurrences of the ``needle``, until no more occurrences can be found and find returns ``false``.

In general, search algorithms can be divided into algorithms that preprocess the ``needle`` (online search) or preprocess the ``haystack`` (index search).

Online Search
-------------

For all online search algorithms, the :dox:`Finder` is an iterator that scans over the ``haystack``.
The :dox:`Pattern` is a search algorithm dependent data structure preprocessed from the ``needle``.
The second template argument of the :dox:`Pattern` selects the search algorithm.

Exact Search
^^^^^^^^^^^^

The following code snippet illustrates the usage of online search algorithms in SeqAn using the example of the Hoorspool algorithm :cite:`Horspool1980`.
We begin by creating two strings of type ``char`` containing the ``haystack`` and the ``needle``.

.. includefrags:: demos/tutorial/find/find_exact.cpp
   :fragment: initialization

We then create :dox:`Finder` and :dox:`Pattern` objects of these strings and choose :dox:`HorspoolPattern Horspool` as the specialization in the second template argument of :dox:`Pattern`.

.. includefrags:: demos/tutorial/find/find_exact.cpp
   :fragment: output

Program output:

.. code-block:: console

   [2,4)   mo
   [12,14) mo
   [17,19) mo

Currently the following exact online algorithms for searching a single sequence are implemented in SeqAn:

:dox:`SimplePattern Simple`
  Brute force algorithm

:dox:`HorspoolPattern Horspool`
  :cite:`Horspool1980`

:dox:`BfamPattern Bfam`
  Backward Factor Automaton Matching

:dox:`BndmAlgoPattern BndmAlgo`
  Backward Nondeterministic DAWG Matching

:dox:`ShiftAndPattern ShiftAnd`
  Exact string matching using bit parallelism

:dox:`ShiftOrPattern ShiftOr`
  Exact string matching using bit parallelism

... and for multiple sequences:

:dox:`WuManberPattern WuManber`
  Extension of :dox:`HorspoolPattern Horspool`.

:dox:`MultiBfamPattern MultiBfam`
  Multiple version of :dox:`BfamPattern Bfam`, uses an automaton of reversed needles.

:dox:`SetHorspoolPattern SetHorspool`
  Another extension of :dox:`HorspoolPattern Horspool` using a trie of reversed needles.

:dox:`AhoCorasickPattern AhoCorasick`
  :cite:`Aho1975`

:dox:`MultipleShiftAndPattern MultipleShiftAnd`
  Extension of :dox:`ShiftAndPattern ShiftAnd`, should only be used if the sum of needle lengths doesn't exceed the machine word size.

Assignment 1
""""""""""""

.. container:: assignment

   Type
    Review

   Objective
    Use the given code example from below.
    Extend the code to search the given ``haystack`` simultaneously for "mo", "send" and "more".
    For every match output the begin and end position in the ``haystack`` and which ``needle`` has been found.

   Hint
     Online search algorithms for multiple sequences simply expect needles of type ``String<String<...> >``.

     .. includefrags:: demos/tutorial/find/find_assignment1.cpp

     You can use the specialization :dox:`WuManberPattern WuManber`.

   Solution
      Click **more...** to see the solution.

      .. container:: foldable

         .. includefrags:: demos/tutorial/find/find_assignment1_solution.cpp

	 We use a :dox:`Pattern` specialized with the :dox:`WuManberPattern WuManber` algorithm for the search and initialize it with our ``needles`` string.
	 For every match found by :dox:`Finder#find` we output the begin and end position and the match region in the ``haystack`` as well as the index of the found ``needle`` which is returned by ``position(pattern)``.

         .. code-block:: console

	    [2,4)   0   mo
	    [7,11)  1   send
	    [12,14) 0   mo
	    [12,16) 2   more
	    [17,19) 0   mo

Approximate Search
^^^^^^^^^^^^^^^^^^

The approximate search can be used to find segments in the ``haystack`` that are similar to a ``needle`` allowing errors, such as mismatches or indels.
Note that if only mismatches are allowed, the difference of the end and begin position of a match is the length of the found ``needle``.
However, in the case of indels this difference may vary and is only a rough estimate for the length.
Therefore, to find a begin position for a certain end position the :dox:`Finder#findBegin` interface should be used.
The usage is similar to :dox:`Finder#find` and is shown in the next example.
We want to find all semi-global alignments of a ``needle`` "more" with a :dox:`SimpleScore` of at least -2 using the scoring scheme (0,-2,-1) (match,mismatch,gap).

Again, we create ``haystack`` and ``needle`` strings first:

.. includefrags:: demos/tutorial/find/find_approx.cpp
   :fragment: initialization

We then create :dox:`Finder` and :dox:`Pattern` objects of these strings and choose :dox:`DPSearchPattern DPSearch` as the specialization in the second template argument of :dox:`Pattern`.
:dox:`DPSearchPattern DPSearch` expects the scoring function as the first template argument which is :dox:`SimpleScore` in our example.
The pattern is constructed using the ``needle`` as a template and our scoring object is initialized with the appropriate scores for match, mismatch and gap.
As in the previous example, the main iteration uses :dox:`Finder#find` to iterate over all end positions with a minimum best score of -2.
If such a semi-global alignment end position is found the begin position is searched via :dox:`Finder#findBegin`.
Please note that we have to set the minimum score to the score of the match found (:dox:`LocalAlignmentEnumerator#getScore`) in order to find the begin of a best match.
We then output all begin and end positions and the corresponding ``haystack`` segment for each match found.

.. includefrags:: demos/tutorial/find/find_approx.cpp
   :fragment: output

Program output:

.. code-block:: console

   [2,4)   mo
   [12,14) mo
   [12,15) mor
   [12,16) more
   [12,17) more
   [12,18) more m
   [17,19) mo
   [17,21) mone

The following specializations are available:

Specialization :dox:`DPSearchPattern DPSearch`
  Dynamic programming algorithm for many kinds of scoring scheme

Specialization :dox:`MyersPattern Myers`
  :cite:`Myers1999`, :cite:`Ukkonen1985`

Specialization :dox:`PexPattern Pex`
  :cite:`BaezaYates1999`

Specialization :dox:`AbndmAlgoPattern AbndmAlgo`
  Approximate Backward Nondeterministic DAWG Matching, adaption of :dox:`AbndmAlgoPattern AbndmAlgo`

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the example from above.
     Modify the code to search with the :dox:`MyersPattern Myers` algorithm for matches of ``"more"`` with an edit distance of at most 2.

   Solution
      Click **more...** to see the solution.

      .. container:: foldable

         .. includefrags:: demos/tutorial/find/find_assignment2_solution.cpp

	 We again set the ``needle`` to ``"more"``.
	 We then change the specialization tag of the :dox:`Pattern` to :dox:`MyersPattern Myers` with default arguments.
	 As :dox:`MyersPattern Myers` algorithm is only applicable to edit distance searches it cannot be specialized or initialized with a scoring scheme.
	 In SeqAn, edit distance corresponds to the scoring scheme (0,-1,-1) (match, mismatch, gap) and an edit distance of 2 corresponds to a minimum score of -2 given to the :dox:`Finder#find` function.

	 The program's output is as follows.

         .. code-block:: console

	    [2,4)   mo
	    [2,5)   mon
	    [2,6)   mon,
	    [12,14) mo
	    [12,15) mor
	    [12,16) more
	    [12,17) more
	    [12,18) more m
	    [17,19) mo
	    [17,20) mon
	    [17,21) mone
	    [17,22) money

Index Search
------------

Exact Search
^^^^^^^^^^^^

For the index based search the :dox:`Finder` needs to be specialized with an :dox:`Index` of the ``haystack`` in the first template argument.
The index itself requires two template arguments, the ``haystack`` type and a index specialization.
In contrast, since the ``needle`` is not preprocessed the second template argument of the :dox:`Pattern` has to be omitted.
The following source illustrates the usage of an index based search in SeqAn using the example of the :dox:`IndexEsa` index (an enhanced suffix array index).
This is the default index specialization if no second template argument for the index is given.
We begin to create an index object of our ``haystack`` ``"tobeornottobe"`` and a ``needle`` ``"be"``.

.. includefrags:: demos/tutorial/find/find_index.cpp
   :fragment: initialization

We proceed to create a :dox:`Pattern` of the needle and conduct the search in the usual way.

.. includefrags:: demos/tutorial/find/find_index.cpp
   :fragment: output

Instead of creating and using a pattern solely storing the ``needle`` we can pass the needle directly to :dox:`Finder#find`.
Please note that an :dox:`Index` based :dox:`Finder` has to be reset with :dox:`Finder#clear` before conducting another search.

.. includefrags:: demos/tutorial/find/find_index.cpp
   :fragment: output_short

Program output:

.. code-block:: console

    [11,13) be
    [2,4)   be
    [11,13) be
    [2,4)   be

All indices also support :dox:`StringSet` texts and can therefore be used to search multiple ``haystacks`` as the following example shows.
We simply exchange the :dox:`CharString` of the haystack with a :dox:`StringSet` of :dox:`CharString` and append some strings to it.

.. includefrags:: demos/tutorial/find/find_index_multiple.cpp
   :fragment: initialization

The rest of the program remains unchanged.

.. includefrags:: demos/tutorial/find/find_index_multiple.cpp
   :fragment: output

.. code-block:: console

   [< 0 , 11 >,< 0 , 13 >) be
   [< 1 , 3 >,< 1 , 5 >)   be
   [< 2 , 0 >,< 2 , 2 >)   be
   [< 0 , 2 >,< 0 , 4 >)   be

The following index specializations support the :dox:`Finder` interface as described above.

Specialization :dox:`IndexEsa`
  Enhanced suffix array based index.
  Supports arbitrary needles.

Specialization :dox:`IndexQGram`
  Q-gram index.
  Needle mustn't exceed the size of the q-gram.

Specialization :dox:`OpenAddressingQGramIndex Open Adressing QGram Index`
  Q-gram index with open addressing.
  Supports larger q-grams.
  Needle and q-gram must have the same size.

Besides the :dox:`Finder#find` interface there is another interface for indices using suffix tree iterators to search exact ``needle`` occurrences described in the tutorial :ref:`tutorial-indices`.

Assignment 3
""""""""""""

.. container:: assignment

     Type
       Application

     Objective
       Modify the example above to search with a :dox:`OpenAddressingQGramIndex Open Adressing QGram Index` q-gram index for matches of "tobe" in "tobeornottobe".

     Solution
      Click **more...** to see the solution.

      .. container:: foldable

         .. includefrags:: demos/tutorial/find/find_assignment3_solution.cpp

	 We simply add a second template argument to the definition of the :dox:`Index` as described in the documentation of the :dox:`OpenAddressingQGramIndex Open Adressing QGram Index`.
	 As shape we can use an :dox:`UngappedShape` of length 4.

	 Program output:

         .. code-block:: console

	    [0,4)   tobe
	    [9,13)  tobe

Approximate Filtration
^^^^^^^^^^^^^^^^^^^^^^

Currently there are no indices directly supporting an approximate search.
But nevertheless, there are approximate search filters available that can be used to filter out regions of the ``haystack`` that do not contain an approximate match, see :dox:`SwiftFinder` and :dox:`SwiftPattern`.
The regions found by these filters potentially contain a match and must be verified afterwards.
:dox:`Finder#beginPosition`, :dox:`Finder#endPosition` and :dox:`Finder#infix` can be used to return the boundaries or sequence of such a potential match.
For more details on using filters, see the article :ref:`how-to-filter-similar-sequences`.
