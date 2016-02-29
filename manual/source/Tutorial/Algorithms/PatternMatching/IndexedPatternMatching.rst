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
40 min

Prerequisites
:ref:`tutorial-datastructures-sequences`, :ref:`tutorial-datastructures-indices`

Pattern matching is about searching a known string or :dox:`StringSet` (``needle``) in another string or :dox:`StringSet` (``haystack``).
This tutorial will introduce you into the SeqAn classes finder and pattern.
It will demonstrate how to use the spezializations of the class finder to perform either an online search or an index based seach.
And you will learn how to specify the search algorithm, which can be exact or approximate.

Overview
--------

.. image:: ../../../under_construction.jpg
