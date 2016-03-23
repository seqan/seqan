.. We create this roles for putting the "Introduction: etc. headings
    on this page without them displaying in the ToC.  This would break
    rendering the ToC correctly on readthedocs style.  The rubric
    headings are formatted using CSS.

.. role:: rubric-heading1
    :class: rubric-heading1
.. role:: rubric-heading2
    :class: rubric-heading2

.. _tutorial-datastructures:

Data Structures
===============

.. toctree::
    :hidden:
    :titlesonly:

    Sequence/index
    Indices/index
    Alignment/index
    Store/index
    Graphs
    Seeds
    Modifiers
    JournalStringTree

SeqAn has numerous data structures that are helpful for analyzing biological sequences. Those range from simple containers for strings that can be saved in different ways, to collection of strings or compressed strings.

The Journaled string tree, for example allows the user to traverse all sequence contexts, given a window of a certain size, that are present in a set of sequences.
Similar sequences are hence only traversed once, and the coordinate bookkeeping is all within the data  structure. This allows for example speedup of up to a 100x given sequences from the 1000 Genome project and compared to traversing the sequences one after another.

Another strong side of SeqAn are its generic string indices. You can think “suffix tree” but the implementations range from an enhanced suffix array to (bidirectional) FM-indices.

In this section you find tutorials addressing the most common of SeqAn's data structures.

The tutorials under :ref:`tutorial-datastructures-sequences` will introduce you to alphabets, sequence containers, iterators and various kinds of string sets, among them also comrpessed, reference based representations.
The tutorials under :ref:`tutorial-datastructures-indices` will introduce you to the interfaces and implementations of SeqAn's string and q-gram indices and the corresponding interators.

The tutorials under :ref:`tutorial-datastructures-alignment-alignment-gaps` will introduce you to how SeqAn implements alignment objects (e.g. gaps in sequences).
The tutorials under :ref:`tutorial-datastructures-store` will introduce you to SeqAn's store data structures for NGS read mapping and annotation

The tutorials under :ref:`tutorial-datastructures-graphs` will introduce you to SeqAn's graph type. Its simple and we provide in the algorithms section various standard graph algorithms for the datatype.
The tutorials under :ref:`tutorial-datastructures-seeds` will introduce you to seeds in SeqAn.
Its what you need if you think 'seed-and-extend'.

The tutorials under :ref:`tutorial-datastructures-modifiers` will introduce you SeqAn's modifier concept. If you want the reverse complement of a string, there is no need to explicitely allocate memory for it. Modifiers leave the sequence as it is but provide a different access to it.

The tutorials under :ref:`tutorial-datastructures-journaledstringtree` will introduce you to a data structure which allows you to have data parallel access to a collection of similar strings. If you have an algorithm that scans collections of sequences one after another and left to right, this will be of interest for you, since it avoids a lot of redundant work if the sequences in the collection are similar.


