.. sidebar:: ToC

   .. contents::

.. _tutorial-algorithms-pattern-matching-online:

Online Pattern Matching
=======================

Learning Objective
  In this tutorial you will learn how to use the SeqAn classes finder and pattern to search a known pattern in a string or :dox:`StringSet`.

Difficulty
  Average

Duration
  40 min

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-datastructures-indices`

For all online search algorithms, the :dox:`Finder` is an iterator that scans over the ``haystack``.
The :dox:`Pattern` is a search algorithm dependent data structure preprocessed from the ``needle``.
The second template argument of the :dox:`Pattern` selects the search algorithm.

Exact Search
------------

The following code snippet illustrates the usage of online search algorithms in SeqAn using the example of the Hoorspool algorithm :cite:`Horspool1980`.
We begin by creating two strings of type ``char`` containing the ``haystack`` and the ``needle``.

.. includefrags:: demos/tutorial/pattern_matching/find_exact.cpp
   :fragment: initialization

We then create :dox:`Finder` and :dox:`Pattern` objects of these strings and choose :dox:`HorspoolPattern Horspool` as the specialization in the second template argument of :dox:`Pattern`.

.. includefrags:: demos/tutorial/pattern_matching/find_exact.cpp
   :fragment: output

Program output:

.. includefrags:: demos/tutorial/pattern_matching/find_exact.cpp.stdout

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
^^^^^^^^^^^^

.. container:: assignment

   Type
    Review

   Objective
    Use the given code example from below.
    Extend the code to search the given ``haystack`` simultaneously for "mo", "send" and "more".
    For every match output the begin and end position in the ``haystack`` and which ``needle`` has been found.

   Hint
     Online search algorithms for multiple sequences simply expect needles of type ``String<String<...> >``.

     .. includefrags:: demos/tutorial/pattern_matching/assignment1.cpp

     You can use the specialization :dox:`WuManberPattern WuManber`.

   Solution
      Click **more...** to see the solution.

      .. container:: foldable

         .. includefrags:: demos/tutorial/pattern_matching/assignment1_solution.cpp

	 We use a :dox:`Pattern` specialized with the :dox:`WuManberPattern WuManber` algorithm for the search and initialize it with our ``needles`` string.
	 For every match found by :dox:`Finder#find` we output the begin and end position and the match region in the ``haystack`` as well as the index of the found ``needle`` which is returned by ``position(pattern)``.

         .. includefrags:: demos/tutorial/pattern_matching/assignment1_solution.cpp.stdout

Approximate Search
------------------

The approximate search can be used to find segments in the ``haystack`` that are similar to a ``needle`` allowing errors, such as mismatches or indels.
Note that if only mismatches are allowed, the difference of the end and begin position of a match is the length of the found ``needle``.
However, in the case of indels this difference may vary and is only a rough estimate for the length.
Therefore, to find a begin position for a certain end position the :dox:`Finder#findBegin` interface should be used.
The usage is similar to :dox:`Finder#find` and is shown in the next example.
We want to find all semi-global alignments of a ``needle`` "more" with a :dox:`SimpleScore` of at least -2 using the scoring scheme (0,-2,-1) (match,mismatch,gap).

Again, we create ``haystack`` and ``needle`` strings first:

.. includefrags:: demos/tutorial/pattern_matching/find_approx.cpp
   :fragment: initialization

We then create :dox:`Finder` and :dox:`Pattern` objects of these strings and choose :dox:`DPSearchPattern DPSearch` as the specialization in the second template argument of :dox:`Pattern`.
:dox:`DPSearchPattern DPSearch` expects the scoring function as the first template argument which is :dox:`SimpleScore` in our example.
The pattern is constructed using the ``needle`` as a template and our scoring object is initialized with the appropriate scores for match, mismatch and gap.
As in the previous example, the main iteration uses :dox:`Finder#find` to iterate over all end positions with a minimum best score of -2.
If such a semi-global alignment end position is found the begin position is searched via :dox:`Finder#findBegin`.
Please note that we have to set the minimum score to the score of the match found (:dox:`LocalAlignmentEnumerator#getScore`) in order to find the begin of a best match.
We then output all begin and end positions and the corresponding ``haystack`` segment for each match found.

.. includefrags:: demos/tutorial/pattern_matching/find_approx.cpp
   :fragment: output

Program output:

.. includefrags:: demos/tutorial/pattern_matching/find_approx.cpp.stdout


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
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     Use the example from above.
     Modify the code to search with the :dox:`MyersPattern Myers` algorithm for matches of ``"more"`` with an edit distance of at most 2.

   Solution
      Click **more...** to see the solution.

      .. container:: foldable

         .. includefrags:: demos/tutorial/pattern_matching/assignment2_solution.cpp

	 We again set the ``needle`` to ``"more"``.
	 We then change the specialization tag of the :dox:`Pattern` to :dox:`MyersPattern Myers` with default arguments.
	 As :dox:`MyersPattern Myers` algorithm is only applicable to edit distance searches it cannot be specialized or initialized with a scoring scheme.
	 In SeqAn, edit distance corresponds to the scoring scheme (0,-1,-1) (match, mismatch, gap) and an edit distance of 2 corresponds to a minimum score of -2 given to the :dox:`Finder#find` function.

	 The program's output is as follows.

         .. includefrags:: demos/tutorial/pattern_matching/assignment2_solution.cpp.stdout
