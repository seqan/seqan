.. sidebar:: ToC

    .. contents::

.. _tutorial-algorithms-optimal-search-schemes:

Optimal Search Schemes
======================

Learning Objective
  In this tutorial you will learn how to search a known pattern in a string or :dox:`StringSet` using a bidirectional index.
  The implemented algorithm in SeqAn is based on Optimal Search Schemes from the `FAMOUS paper <https://arxiv.org/abs/1711.02035>`_ which allows for searching for up to 4 errors based on Hamming distance or Edit distance.

Difficulty
  Average

Duration
  15 min

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-datastructures-indices`

Overview
--------

Searching in an index with more than two errors is computationally very expensive and thus trivial approaches such as simple backtracking are not recommended.
Optimal Search Schemes make use of a bidirectional index, split the pattern to be searched into multiple pieces and enumerate it in a more efficient fashion.
To use this algorithm, you can simply call the ``find()`` method.

.. includefrags:: demos/tutorial/indices/find2_index_approx.cpp
      :fragment: SinglePattern

The first argument is the delegate function that will be called if a match of the pattern is found.
It gets an iterator of the bidirectional index and you can choose how to proceed with the hits.
Additionally it passes a reference to the original pattern searched as well as the number of errors for the current hit in the index.

.. includefrags:: demos/tutorial/indices/find2_index_approx.cpp
      :fragment: Delegate

The second argument has to be a bidirectional index, currently SeqAn offers only bidirectional FM indices.
(Please check the corresponding tutorial on :ref:`tutorial-datastructures-indices-fm-index` for more information how FM indices can be configured and constructed.)

The third argument is a :dox:`String` that you want to search.

The number of allowed errors (lower and upper bounds) are passed as template arguments.
Please note that these parameters have to be ``constexpr``.
The distance metric is passed as fourth argument. You can choose between ``HammingDistance()`` and ```EditDistance()`` (also known as Levenshtein distance).

You also have to possibility to pass a :dox:`StringSet` of patterns to search instead of a single :dox:`String`.
The search can then be parallized by passing a fifth parameter tag ``Parallel()`` instead of ``Serial()``.

.. includefrags:: demos/tutorial/indices/find2_index_approx.cpp
      :fragment: MultiplePatterns

.. warning::

      Please be aware that the same hits might be reported multiple times and you might have to filter these duplicated depending on your application, especially with larger errors numbers.

.. note::

      When searching a StringSet in parallel mode the delegate is likely being executed by different threads at the same time.
      You have to take care by yourself that the threads do not interfere, e.g. write not into a variable simultaneously.
      One way of achieving this is by using lock guards.
      The following example buffers all hits in a set (to filter duplicates) and prints them afterwards.

      .. includefrags:: demos/tutorial/indices/find2_index_approx.cpp
            :fragment: ParallelMode

Here is a complete example for searching a string or a stringset (in serial mode):

.. includefrags:: demos/tutorial/indices/find2_index_approx.cpp
      :fragment: Complete
