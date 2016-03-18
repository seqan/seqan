.. sidebar:: ToC

    .. contents::

.. _tutorial-algorithms-pattern-matching-indexed:

Indexed Pattern Matching
========================

Learning Objective
  In this tutorial you will learn how to use the SeqAn classes finder and pattern to search a known pattern in a string or :dox:`StringSet`.

Difficulty
  Average

Duration
  30 min

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-datastructures-indices`, :ref:`tutorial-algorithms-pattern-matching-online`

Overview
--------

The :dox:`Finder` is an object that stores all necessary information for searching for a pattern. We have learned how to use it in :ref:`tutorial-algorithms-pattern-matching-online` tutorial.

The following line of code shows how the :dox:`Finder` is initialized with an index. In this example, we search for the pattern ``ACGT``.

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: finder

Calling the function :dox:`Finder#find` invokes the localization of all occurrences of a given pattern.
It works by modifying pointers of the ``Finder`` to tables of the index.
For example, the :dox:`Finder` of ``esaIndex`` stores two pointers, pointing to the first and last suffix array entry that stores an occurrence of the pattern. The return value of the :dox:`Finder#find` function tells us whether or not a given pattern occurs in the text. Furthermore, if there are several instances of a pattern, consecutive calls of :dox:`Finder#find` will modify the :dox:`Finder` such that it points to the next occurrence after each call:

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: finder_multiple

The above code is not very useful, since we do not know the locations of the first, second or third pattern occurrence.
The function :dox:`Finder#position` will help here.
:dox:`Finder#position` called on a finder returns the location of the ``x``\ th pattern, where ``x`` can be the first, second, or any other occurrence of the pattern.

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: finder_position

.. tip::

   Indices in SeqAn are built on demand.
   That means that the index tables are not build when the constructor is called, but when we search for a pattern for the first time.

Exact Search
------------

For the index based search the :dox:`Finder` needs to be specialized with an :dox:`Index` of the ``haystack`` in the first template argument.
The index itself requires two template arguments, the ``haystack`` type and a index specialization.
In contrast, since the ``needle`` is not preprocessed the second template argument of the :dox:`Pattern` has to be omitted.
The following source illustrates the usage of an index based search in SeqAn using the example of the :dox:`IndexEsa` index (an enhanced suffix array index).
This is the default index specialization if no second template argument for the index is given.
We begin to create an index object of our ``haystack`` ``"tobeornottobe"`` and a ``needle`` ``"be"``.

.. includefrags:: demos/tutorial/pattern_matching/find_index.cpp
   :fragment: initialization

We proceed to create a :dox:`Pattern` of the needle and conduct the search in the usual way.

.. includefrags:: demos/tutorial/pattern_matching/find_index.cpp
   :fragment: output

Instead of creating and using a pattern solely storing the ``needle`` we can pass the needle directly to :dox:`Finder#find`.
Please note that an :dox:`Index` based :dox:`Finder` has to be reset with :dox:`Finder#clear` before conducting another search.

.. includefrags:: demos/tutorial/pattern_matching/find_index.cpp
   :fragment: output_short

Program output:

.. includefrags:: demos/tutorial/pattern_matching/find_index.cpp.stdout


All indices also support :dox:`StringSet` texts and can therefore be used to search multiple ``haystacks`` as the following example shows.
We simply exchange the :dox:`CharString` of the haystack with a :dox:`StringSet` of :dox:`CharString` and append some strings to it.

.. includefrags:: demos/tutorial/pattern_matching/find_index_multiple.cpp
   :fragment: initialization

The rest of the program remains unchanged.

.. includefrags:: demos/tutorial/pattern_matching/find_index_multiple.cpp
   :fragment: output

.. includefrags:: demos/tutorial/pattern_matching/find_index_multiple.cpp.stdout


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

Besides the :dox:`Finder#find` interface there is another interface for indices using suffix tree iterators to search exact ``needle`` occurrences described in the tutorial :ref:`tutorial-datastructures-indices`.

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

     Type
       Application

     Objective
       Modify the example above to search with a :dox:`OpenAddressingQGramIndex Open Adressing QGram Index` q-gram index for matches of "tobe" in "tobeornottobe".

     Solution
      Click **more...** to see the solution.

      .. container:: foldable

         .. includefrags:: demos/tutorial/pattern_matching/assignment3_solution.cpp

	 We simply add a second template argument to the definition of the :dox:`Index` as described in the documentation of the :dox:`OpenAddressingQGramIndex Open Adressing QGram Index`.
	 As shape we can use an :dox:`UngappedShape` of length 4.

	 Program output:

         .. includefrags:: demos/tutorial/pattern_matching/assignment3_solution.cpp.stdout

Approximate Filtration
----------------------

Currently there are no indices directly supporting an approximate search.
But nevertheless, there are approximate search filters available that can be used to filter out regions of the ``haystack`` that do not contain an approximate match, see :dox:`SwiftFinder` and :dox:`SwiftPattern`.
The regions found by these filters potentially contain a match and must be verified afterwards.
:dox:`Finder#beginPosition`, :dox:`Finder#endPosition` and :dox:`Finder#infix` can be used to return the boundaries or sequence of such a potential match.
For more details on using filters, see the article :ref:`how-to-recipes-filter-similar-sequences`.
