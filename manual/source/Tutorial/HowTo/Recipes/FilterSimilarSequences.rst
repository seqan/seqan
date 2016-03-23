.. sidebar:: ToC

    .. contents::

.. _how-to-recipes-filter-similar-sequences:

Filtering Similar Sequences
===========================

Using Swift
-----------

In the next example we are going to use the **Swift** filter to efficiently find pairs of similar reads.
The Swift algorithm searches for so-called epsilon matches (local alignments) of two sequences with an error rate below a certain epsilon threshold.

The Swift implementation in SeqAn provides a :dox:`Finder#find` interface and requires the :dox:`Finder` and :dox:`Pattern` to be specialized with ``Swift<..>``.
Millions of sequences can be searched simultaneously with one :dox:`SwiftPattern Swift Pattern` in a :dox:`SwiftFinder Swift Finder` of a single haystack sequence.
The error rate of a local alignment is the number of errors divided by the length of the needle sequence part of the match.
There are currently two version of the Swift algorithm implemented in SeqAn: ``SwiftSemiGlobal`` and ``SwiftLocal``.
Both can be used to search epsilon-matches of a certain minimum length.

.. hint::

   ``SwiftSemiGlobal`` should only be used for short needles (sequenced reads) as it always returns potential epsilon matches spanning a whole needle sequence.
   ``SwiftLocal`` should be preferred for large needles as it returns needle sequences potentially having an intersection with an epsilon match.

The following program searches for semi-global alignments between pairs of reads with a maximal error rate of 10%.

.. includefrags:: demos/howto/filter_similar_sequences.cpp
   :fragment: includes

First we loads reads from a file into a :dox:`FragmentStore` with :dox:`FragmentStore#loadReads`.

.. includefrags:: demos/howto/filter_similar_sequences.cpp
   :fragment: load_reads

Swift uses a q-gram index of the needle sequences.
Thus, we have to specialize the :dox:`SwiftSemiGlobalPattern Swift Semi Global Pattern` with a :dox:`IndexQGram` index of the needle :dox:`StringSet` in the first template argument, create the index over the :dox:`FragmentStore::readSeqStore` and pass the index to the :dox:`SwiftSemiGlobalPattern Pattern` constructor.
:dox:`SwiftSemiGlobalFinder Swift Semi Global Finder` and :dox:`SwiftSemiGlobalPattern Swift Semi Global Pattern` classes have to be specialized with ``SwiftSemiGlobal`` in the second template argument.

.. note::

   Note, to use the local swift filter you simply switch the specialization tag to ``SwiftLocal``: :dox:`SwiftLocalFinder Swift Local Finder` and :dox:`SwiftLocalPattern Swift Local Pattern`.

The main loop iterates over all potential matches which can be further processed, e.g. by a semi-global or overlap aligner.

.. includefrags:: demos/howto/filter_similar_sequences.cpp
   :fragment: filter

