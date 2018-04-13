.. sidebar:: ToC

    .. contents::

.. _tutorial-datastructures-indices-string-indices:

String Indices
==============

Learning Objective
  You will get an overview of the different kinds of indices in SeqAn and how they are used.

Difficulty
  Average

Duration
  1 h
Prerequisites
  :ref:`tutorial-datastructures-sequences`

Indices in SeqAn
----------------

Indices in SeqAn are substring indices, meaning that they allow efficient pattern queries in strings or sets of strings.
In contrast to, e.g., online-search algorithms that search through the text in :math:`\mathcal{O}(n)`, substring indices find a pattern in sublinear time :math:`o(n)`.

You can find the following indices in SeqAn:

:dox:`IndexSa`
  Suffix Array :cite:`Manber1993`
:dox:`IndexEsa`
  Extended Suffix Array :cite:`Abouelhoda2004`
:dox:`IndexWotd`
  Lazy suffix tree :cite:`Giegerich2003`
:dox:`IndexDfi`
  Deferred Frequency Index :cite:`Weese2008`
:dox:`IndexQGram`
  Q-gram index (see `here <QgramIndex.html>`_ for more details)
:dox:`FMIndex`
  Full-text minute index (see the `FMIndex <FMIndex.html>`_ for more details) :cite:`Ferragina2001`

Index Construction
------------------

We will now show how we can create the different indices in SeqAn before we show how they are used for pattern search.

All the mentioned indices belong to the generic :dox:`Index` class.
A SeqAn index needs two pieces of information: the type of the :dox:`String` or :dox:`StringSet` to be indexed and the index specialization, such as :dox:`IndexEsa` or :dox:`FMIndex`.

.. important::

    Indices based on suffix arrays (also including the FM index) are built using secondary memory.
    When building large indices, it is therefore possible to run out of disk space (in which case an exception will be
    thrown).
    To circumvent this, the directory used for temporary storage can be changed by specifying the TMPDIR environment variable (on UNIX)
    respectively TEMP environment variable (on Windows):

    .. code-block:: console

       # export TMPDIR=/somewhere/else/with/more/space

    .. code-block:: console

       # SET TEMP=C:\somewhere\else\with\more\space

The following code snippet creates an enhanced suffix array index of a string of type :dox:`Dna5`.

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: esa

In contrast, the next code snipped creates a FM index over a set of amino acid sequences:

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: fm

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Copy the code below and

     #. change it to build an :dox:`IndexEsa` over a string of type :dox:`Dna`,
     #. add an :dox:`IndexEsa` over a :dox:`StringSet` of :dox:`String Strings` of type :dox:`Dna`.

     .. includefrags:: demos/tutorial/indices/base.cpp
            :fragment: assignment1

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/indices/assignment_1.cpp

Assignment 2
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     Write a small program that prints the locations of all occurrences of ``"TATAA"`` in ``"TTATTAAGCGTATAGCCCTATAAATATAA"``.

   Hints
    Use the :dox:`Finder#find` function as the conditional instruction of a <tt>while</tt> loop.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/indices/assignment_2.cpp

You might have noticed that we only applied the :dox:`FMIndex` and :dox:`IndexEsa` in the examples.
The reason for this is that even though everything stated so far is true for the other indices as well, :dox:`IndexWotd` and :dox:`IndexDfi` are more useful when used with iterators as explained in the tutorial :ref:`tutorial-datastructures-indices-index-iterators` and the :dox:`IndexQGram` uses :dox:`Shape Shapes` which is also explained in another tutorial.

One last remark is necessary.

.. important::

    If you search for two different patterns with the same :dox:`Finder` object, you have to call the :dox:`Finder#clear` function of the finder between the search for the two patterns.
    Otherwise the behavior is undefined.

Handling Multiple Sequences (StringSets)
----------------------------------------

The previous sections already described how an index of a set of strings can be instantiated.
A character position of a :dox:`StringSet` can be one of the following:

#. A local position (default), i.e. a :dox:`Pair` (seqNo, seqOfs) where seqNo identifies the string within the :dox:`StringSet` and the seqOfs identifies the position within this string.
#. A global position, i.e. a single integer value between 0 and the sum of string lengths minus 1.
   This integer is the position in the gapless concatenation of all strings in the :dox:`StringSet` to a single string.

For indices, the meta-function :dox:`SAValue` determines, which position type (local or global) will be used for internal index tables (suffix array, q-gram array) and what type of position is returned by functions like :dox:`Finder#position` of a :dox:`Finder`.
:dox:`SAValue` returns a :dox:`Pair` (local position) by default, but could be specialized to return an integer type (global position) for some applications.
If you want to write algorithms for both variants you should use the functions :dox:`TextConcept#posLocalize`, :dox:`TextConcept#posGlobalize`, :dox:`TextConcept#getSeqNo`, and :dox:`TextConcept#getSeqOffset`.

Storing and Loading
-------------------

Storing and loading an index can be done with:

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: save

or

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: open

If you have built your q-gram index with variable shapes (i.e. :dox:`SimpleShape` :dox:`GenericShape`), you have to keep in mind that q or the shape is not stored or loaded.
This must be done manually and directly before or after loading with :dox:`Shape#resize` or :dox:`Shape#stringToShape`.

A newly instantiated index is initially empty.
If you assign a text to be indexed, solely the text fibre is set.
All other fibres are empty and created on demand.
Normally, a full created index should be saved to disk.
Therefore, you have to create the required fibres explicitly by hand.

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: require

For the :dox:`IndexEsa` index you could do:

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: require2

Indexes based on external strings, e.g.  ``Index<String<Dna,External<> >,IndexEsa<> >`` or ``Index<String<Dna,MMap<> >,IndexEsa<> >`` cannot be saved, as they are persistent implicitly.
The first thing after instantiating such an index should be associating it to a file with:

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: external

The file association implies that any change on the index, e.g. fibre construction, is synchronized to disk.
When instantiating and associating the index the next time, the index contains its previous state and all yet to be constructed fibres.

Reducing the memory consumption
-------------------------------

All :dox:`Index Indices` in SeqAn are capable of indexing :dox:`String Strings` or :dox:`StringSet StringSets` of arbitrary sizes, i.e. up to 2^64 characters.
This always comes at a cost in terms of memory consumption, as any :dox:`Index` has to represent 64 bit positions in the underlying text.
However, in many practical instances, the text to be indexed is shorter, e.g. it does not exceed 4.29 billion (2^32) characters.
In this case, one can reduce the memory consumption of an :dox:`Index` by changing its internal data types, with no drawback concerning running time.

SA Fibre
^^^^^^^^

All :dox:`Index Indices` in SeqAn internally use the :dox:`Fibre FibreSA`, i.e. some sort of suffix array.
For :dox:`String Strings`, each suffix array entry consumes 64 bit of memory per default, where 32 bit would be sufficient if the text size is appropriate.
In order to change the size type of the suffix array entry we simply have to overload the metafunction :dox:`SAValue`.

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: SAValue

If your text is a :dox:`StringSet`, then :dox:`SAValue` will return a :dox:`Pair` that can be overloaded in the same way.

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: SAValue2

The first type of the pair is used as the type for the index of a string in the string set.
So if you only have a few strings you could save even more memory like this.

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: SAValue3

FMIndex Fibres
^^^^^^^^^^^^^^

The size of a generalized :dox:`FMIndex` depends also on the total number of characters in a :dox:`StringSet` (see :dox:`StringSet#lengthSum`).
This trait can be configured via the :dox:`FMIndexConfig` object. For more information, see the :ref:`tutorial-datastructures-indices-fm-index` section.

.. includefrags:: demos/tutorial/indices/base.cpp
      :fragment: config

Other Index Fibres
^^^^^^^^^^^^^^^^^^

See :ref:`how-to-recipes-access-index-fibres` for more information.
