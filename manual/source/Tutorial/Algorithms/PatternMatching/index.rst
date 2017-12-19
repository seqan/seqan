.. _tutorial-algorithms-pattern-matching:

Pattern Matching
================

..  toctree::
    :hidden:
    :titlesonly:

    OnlinePatternMatching
    IndexedPatternMatching
    OptimalSearchSchemes

Pattern matching is about searching a known string or :dox:`StringSet` (``needle``) in another string or :dox:`StringSet` (``haystack``).
This tutorial will introduce you into the SeqAn classes :dox:`Finder` and :dox:`Pattern`.
It will demonstrate how to use the spezializations of the class finder to perform either an online search or an index based seach.
And you will learn how to specify the search algorithm, which can be either exact or approximate.

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
The additional section on Optimal Search Schemes is also an indexed search algorithm but due to a different interface separated into its own section.
