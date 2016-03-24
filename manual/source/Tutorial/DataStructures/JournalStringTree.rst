.. sidebar:: ToC

    .. contents::

.. _tutorial-datastructures-journaledstringtree:

Journaled String Tree
=====================

Learning Objective
  This tutorial introduces you to the new data structure Journaled String Tree.
  You will learn how to use this data structure and how to exploit it for an efficient analysis while searching multiple sequences simultaneously.

Difficulty
  Average

Duration
  45 min

Prerequisites
  :ref:`tutorial-datastructures-sequences`, :ref:`tutorial-datastructures-sequences-string-sets`

A typical task in bioinformatics is to find patterns in biological sequences e.g. transcription factors, or to examine different biological traits and the effects of modifications on such traits.
With the current advances in sequencing technologies, sequences of whole populations have been made available.
But the time for searching in all these sequences is proportional to the number of sequences searched.
That's why it is important to find novel strategies to cope with the deluge of sequencing data.
Since, many biological problems often involve the analysis of sequences of the same species, one effective strategy would be to exploit the similarities of such sequences.

For this special purpose we will introduce you to two data structures that are designed to improve these algorithmic tasks.
The first one is the :dox:`JournaledString` and the second is the :dox:`JournaledStringTree`.

Journaled String
----------------

The :dox:`JournaledString` data structure behaves like a normal :dox:`String` in SeqAn, except that it is composed of two data structures.

#. The first data structure is a :dox:`Holder` which stores a sequence.
#. The second data structure stores modifications that are made to this particular sequence using a **journal** (see `Journaling Filesystems <http://en.wikipedia.org/wiki/Journaling_file_system>`_ for more information).
   This journal contains a list of deletions and insertions.
   The inserted characters are stored in an additional **insertion buffer**.

The advantage of this data structure lies in representing a String as a "patch" to another String.
The journaled data structure can be modified without loosing the original context.
We want to show you how to work with these data structures so you can build your own algorithms based on this.

First, we show you how to work with the Journaled String so you can learn the basic principles.
To get access to the Journaled String implementation you have to include the ``<seqan/sequence_journaled.h>`` header file.
Note that you will need the ``<seqan/stream.h>`` too in order to print the sequences.

.. includefrags:: demos/tutorial/journaled_set/example_journal_string_basic.cpp
   :fragment: main

In the next step we define the Journaled String type.
A Journaled String is a specialization of the String class and is defined as ``String<TValue, Journaled<THostSpec, TJournalSpec, TBufferSpec> >``.
The specialization takes two parameters: (1) ``TValue`` defines the alphabet type used for the Journaled String and (2) ``Journaled<>`` selects the Journaled String specialization of the String class.

``Journaled<>`` is further specialized with

* ``THostSpec`` selects the specialization of the underlying host sequence (``Alloc<>`` for [dox:AllocString Alloc String),
* ``TJournaleSpec`` selects the used implementation to manage the journaled differences (here: ``SortedArray``), and
* ``TBufferSpec`` selects the used specialization for the internally managed insertion buffer (here: ``Alloc<>`` as well).

In our scenario we use a ``char`` alphabet and [dox:AllocString Alloc String for the host string and the insertion buffer.
Additionally, we use a ``Sorted Array`` as the model to manage the recorded differences.

We use the metafunction :dox:`HostedConcept#Host` to get the type of the underlying host string used for the Journaled String.

.. includefrags:: demos/tutorial/journaled_set/example_journal_string_basic.cpp
   :fragment: typedef

Now we can define the variables holding data structures.
First, we construct our host sequence and after that we construct the Journaled String.
Then, we set the host sequence using the function :dox:`JournaledString#setHost`.
Afterwards, we examine the data structure in more detail and print the host sequence the constructed journaled sequence and the nodes of it.

.. includefrags:: demos/tutorial/journaled_set/example_journal_string_basic.cpp
   :fragment: init

.. tip::

    The Journal

    A node in the Journaled String represents either a part of the host sequence or a part of the insertion buffer.
    The type of a node is distinguished by the member variable **segmentSource** and can be of value ``SOURCE_ORIGINAL`` to refere to a part in the host or ``SOURCE_PATCH`` to refere to a part in the insertion buffer.
    A node further consists of three variables which specify the **virtual position**, the **physical position** and the **length** of this part.
    The **virtual position** gives the relative position of the Journaled String after all modifications before this position have been "virtually" applied.
    The **physical position** gives the absolute position where this part of the journal maps to either the host sequence or the insertion buffer.

This is followed by modifying our Journaled String.
We insert the string ``"modified"`` at position ``7`` and delete the suffix ``"sequence"`` at position ``19``.
Note that position ``19`` refers to the string after the insertion of ``"modified"`` at position ``7``.
Again we print the host, the journaled sequence and the nodes that represent the modifications to see how our changes affect the host and the journaled sequence.

.. includefrags:: demos/tutorial/journaled_set/example_journal_string_basic.cpp
   :fragment: modification

All of this is followed by calling :dox:`JournaledString#flatten` on our journeld string.
This call applies all journaled changes to the host sequence.
Again we print the sequences to see the effects.

.. includefrags:: demos/tutorial/journaled_set/example_journal_string_basic.cpp
   :fragment: flatten

Here is the output of our small program.

.. includefrags:: demos/tutorial/journaled_set/example_journal_string_basic.cpp.stdout


.. important::

   Be careful when using the :dox:`JournaledString#flatten` function as it modifies the underlying host sequence.
   This might affect other journaled sequences that share the same host sequence.
   This becomes important especially when working with Journaled Sets where a whole set of sequences is journaled to the same reference.

Simultaneously Searching Multiple Sequences
-------------------------------------------

Now, we come to a simple example to demonstrate the use of the :dox:`JournaledStringTree Journaled String Tree` (JST).
As you could imagine, the JST internally uses a set of :dox:`JournaledString Journaled Strings` to buffer the sequences, while requiring only a low memory footprint.

In this article, we are going to create a small JST, which we will use to search for a pattern using the :dox:`HorspoolPattern Horspool` algorithm.

Let's just start with the include headers.
In order to make the JST implementation visible to our source code we need to include the header ``include <seqan/journaled_string_tree.h>``.

.. includefrags:: demos/tutorial/journaled_string_tree/journaled_string_tree_base.cpp
    :fragment: include

In the next step, we are going to define the type of the JST.

.. includefrags:: demos/tutorial/journaled_string_tree/journaled_string_tree_base.cpp
    :fragment: typedef

The only information that is required is the type of the reference sequence used for the underlying sequence.
We also defined the pattern type and  a traverser type, which we will explain soon.

Now, we are ready to initialize the JST.
To construct a JST, we need to know the reference sequence and how many sequences should be represented by the JST.
In our case we assume 10 sequences.
The JST supports insertion or deletion of :dox:`DeltaTypeTags delta events`.
A delta event is a tuple consisting of four parameters: The reference position, the value, the coverage and the delta type.
The reference position determines the position within the reference sequence, where this event occurs.
The value represents the actual modification applied to the sequences, that are determined by the coverage.
The type of the value depends on the delta type.

.. tip::
    The internal types, e.g. the types of the different delta events, of the :dox:`JournaledStringTree JST` can be overloaded with a second optional traits object.
    If no trait object is given :dox:`DefaultJstConfig` is taken as default.
    See the API documentation for more information.

The following listing creates a JST and inserts some delta events into the object:

.. includefrags:: demos/tutorial/journaled_string_tree/journaled_string_tree_base.cpp
    :fragment: init

After creating the JST, we can now prepare the search.
To do so, we first define a needle that we want to search.
Second, we need to instantiate a traverser object.
A traverser represents the current state of the traversal over the JST.
It is comparable to an iterator, but it is not lightweight, as it uses a state stack to implement the traversal over the JST.
The traverser is initialized with two arguments: The instance of the JST and the context length, which is in our case the length of the needle.

Here is the listing:

.. includefrags:: demos/tutorial/journaled_string_tree/journaled_string_tree_base.cpp
    :fragment: prepare_search

In line 4 and 5 in the listing above we initialize the pattern with the needle and then create an ``JstExtension`` object.
This ``JstExtension`` is needed to extend the :dox:`Pattern Pattern class` of SeqAn with some auxiliary functions necessary for the JST based search.
The only thing required, is that ``pat`` is fully initialized when passing it to ``ext``.

The last preparation step we need before invoking the search algorithm is to create a functor that is called, whenever the search algorithm finds a match.
In our scenario we simply want to print the sequences and the positions where the hit occurs.
Therefor we create a simple ``MatchPrinter`` functor:

.. includefrags:: demos/tutorial/journaled_string_tree/journaled_string_tree_base.cpp
    :fragment: match_printer

This match printer, holds a reference to the actual traverser.
So we can call the ``position`` function on the traverser, when the function-call-operator is invoked by the search algorithm.

Now we can invoke the search using the ``find`` interface:

.. includefrags:: demos/tutorial/journaled_string_tree/journaled_string_tree_base.cpp
    :fragment: search

And finally the output:

.. includefrags:: demos/tutorial/journaled_string_tree/journaled_string_tree_base.cpp.stdout

The following list gives an overview of the available search algorithms:

Horspool
  Exact online search using :dox:`HorspoolPattern` as base pattern class.

ShiftAnd
  Exact online search using :dox:`ShiftAndPattern` as base pattern class.

ShiftOr
  Exact online search using :dox:`ShiftOrPattern` as base pattern class.

MyersUkkonen
  Approximate online search using :dox:`MyersUkkonen` as base pattern class.

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective

     Use the code from above and find all patterns of the needle ``CCTCCA`` with up to 2 errors.

   Hints
     .. container:: foldable

        When searching with errors, the context size needs to be updated accordingly.

   Solution
     .. container:: foldable

        Since we are trying to find the needle approximatively, we need to use the ``Myers' bitvector`` algorithm.
        Here is the entire solution:

        .. includefrags:: demos/tutorial/journaled_string_tree/journaled_string_tree_base_assignment1.cpp

        And here is the output:

        .. includefrags:: demos/tutorial/journaled_string_tree/journaled_string_tree_base_assignment1.cpp.stdout

