.. sidebar:: ToC

    .. contents::

.. _tutorial-datastructures-indices-fm-index:

FMIndex
=======

Learning Objective
  You will know the features of the FMIndex, how you can optimize it and how to search bidirectionally.

Difficulty
  Average

Duration
  0.5 h

Prerequisites
  :ref:`tutorial-datastructures-sequences`

The FMIndex
-----------

An FMIndex is a space-efficient and fast index structure that can be used to search arbitrary patterns.
In contrast to other indices like the :dox:`IndexSa`, the FMIndex can only be traversed as a prefix trie.
Patterns cannot be searched from left to right but from right to left.
This means that you either have to reverse the text before constructing the index or reverse the patterns before searching.
Also hits of a pattern do not refer to the starting position of the pattern but to the ending position.
This has to be taken into account when retrieving hits.

Different implementations
-------------------------

The FMIndex is based on the Burrow-Wheeler-Transform (BWT) which uses a rank dictionary.
We currently have three different implementations to choose from with different running time and space trade-offs (the space consumption only refers to the rank dictionary, not the entire index).
:math:`n` is the length of the text, :math:`\sigma` the alphabet size.
The running time refers to a single character search.

+----------------------+-----------------------------+----------------------------------+---------------------------------------------------------------------------------+
+ Rank Dictionary Type | Specialization              | Running time                     | Space consumption                                                               |
+======================+=============================+==================================+=================================================================================+
| :dox:`Levels`        | :dox:`LevelsRDConfig`       | :math:`\mathcal{O}(\log \sigma)` | :math:`\mathcal{O}(\log \sigma \cdot n) + o (\log \sigma \cdot \sigma \cdot n)` |
+----------------------+-----------------------------+----------------------------------+---------------------------------------------------------------------------------+
| :dox:`Levels`        | :dox:`LevelsPrefixRDConfig` | :math:`\mathcal{O}(1)`           | :math:`\mathcal{O}(\log \sigma \cdot n) + o (\log \sigma \cdot \sigma \cdot n)` |
+----------------------+-----------------------------+----------------------------------+---------------------------------------------------------------------------------+
| :dox:`WaveletTree`   | :dox:`WTRDConfig`           | :math:`\mathcal{O}(\log \sigma)` | :math:`\mathcal{O}(\log \sigma \cdot n)`                                        |
+----------------------+-----------------------------+----------------------------------+---------------------------------------------------------------------------------+

Generally speaking, :dox:`WaveletTree WaveletTrees` are recommended for larger alphabets, while :dox:`Levels` offer a faster runtime at the expense of a higher space consumption.
:dox:`LevelsRDConfig LevelsRDConfigs` are in practice noticeably faster than WaveletTrees since they only perform :math:`\mathcal{O}(\log \sigma)` logical operations instead of going down a tree of height :math:`\mathcal{O}(\log \sigma)` resulting in more cache misses.
:dox:`LevelsPrefixRDConfig` are recommended for smaller alphabets and bidirectional support on FM indices.
The running time is constant for both backward and forward searches.

The WaveletTree and the PrefixLevels implementations are already wrapped in configuration objects, namely :dox:`FMIndexConfig` and :dox:`FastFMIndexConfig`.

.. includefrags:: demos/tutorial/indices/fm_index.cpp
   :fragment: FMIndexConfigs

Optimizing Space or Running Time
--------------------------------

All three rank dictionaries are based on bit vectors with constant-time rank support.
We skip at this point details on how rank queries work and focus on the parameters and their effects.
To reduce the space consumption there are three parameters that you can adjust:

First of all *TSize* should be the smallest data type that can store the length of the indexed text (e.g. *uint8_t*, *uint16_t*, *uint32_t*, *uint64_t*).
For a :dox:`StringSet` the length is defined as the the sum of lengths of strings plus the number of strings.

The bit vector is clustered into blocks that store the precomputed rank up to every block.
To reduce the space of these blocks, one can add an additional level of blocks on top which group blocks together to a superblock.
We support up to 3 levels.
The effect on the running time is minimal since for every additional level only one more array lookup is conducted.
In practice two levels are in most cases faster than one level due to smaller tables and thus less cache misses.

To reduce the space consumption even further or to improve the running time, one can change the size of a block.
Each block contains :math:`64 \cdot WPB` (WORDS_PER_BLOCK) bits and its rank can thus be computed on a 64 bit machine by :math:`WPB` popcount operations.
By reducing (increasing) the block size, the running time (space consumption) can be improved noticeably.
By default, WPB is set to 1.

For texts with less than :math:`2^{32}` characters, two levels and a fast query time (only one popcount operation per rank query), the FMIndex can be configured using

.. includefrags:: demos/tutorial/indices/fm_index.cpp
   :fragment: FMIndexConfigs2

Bidirectional FMIndex
---------------------

The bidirectional FMIndex can be used to search a pattern into both directions, i.e. a string can be extend by a character to the left or right in an arbitrary manner.
For information on how to search in a bidirectional FMIndex, please check the Tutorial on :ref:`tutorial-datastructures-indices-index-iterators`.

.. includefrags:: demos/tutorial/indices/fm_index.cpp
   :fragment: BidirectionalIndex
