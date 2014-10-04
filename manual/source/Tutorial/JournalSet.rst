.. sidebar:: ToC

   .. contents::


.. _tutorial-data-journaling:

Journaled Set
=============

Learning Objective
  This tutorial introduces you to the new data structures Journaled Set and Journaled String.
  You will learn how to use them and how to exploit these data structures for an efficient analysis, while implementing a native online search.

Difficulty
  Advanced

Duration
  2 h

Prerequisites
  :ref:`tutorial-sequences`, :ref:`tutorial-string-sets`, :ref:`tutorial-iterators`

A typical task in bioinformatics is to find patterns in biological sequences e.g. transcription factors, or to examine different biological traits and the effects of modifications on such traits.
With the current advances in sequencing technologies, sequences of whole populations have been made available.
But the time for searching in all these sequences is proportional to the number of sequences searched.
That's why it is important to find novel strategies to cope with the deluge of sequencing data.
Since, many biological problems often involve the analysis of sequences of the same species, one effective strategy would be to exploit the similarities of such sequences.

For this special purpose we provide two data structures that are designed to improve the algorithmic tasks.
The first one is the :dox:`JournaledString` and the second is the :dox:`JournaledSet`.

In this tutorial, we will introduce you to both data structures and implement a simple online search step by step.

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
For this reason we keep the applicational background simple and implement an native online-search algorithm by which we examine different properties of the data structures.

Before we start implementing our online search algorithm, we show you how to work with the Journaled String to learn the basic principles.
To get access to the Journaled String implementation you have to include the ``<seqan/sequence_journaled.h>`` header file.
Note that you will need the ``<seqan/file.h>`` too in order to print the sequences.

.. includefrags:: extras/demos/tutorial/data_journaling/example_journal_string_basic.cpp
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

.. includefrags:: extras/demos/tutorial/data_journaling/example_journal_string_basic.cpp
   :fragment: typedef

Now we can define the variables holding data structures.
First, we construct our host sequence and after that we construct the Journaled String.
Then, we set the host sequence using the function :dox:`JournaledString#setHost`.
Afterwards, we examine the data structure in more detail and print the host sequence the constructed journaled sequence and the nodes of it.

.. includefrags:: extras/demos/tutorial/data_journaling/example_journal_string_basic.cpp
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

.. includefrags:: extras/demos/tutorial/data_journaling/example_journal_string_basic.cpp
   :fragment: modification

All of this is followed by calling :dox:`JournaledString#flatten` on our journeld string.
This call applies all journaled changes to the host sequence.
Again we print the sequences to see the effects.

.. includefrags:: extras/demos/tutorial/data_journaling/example_journal_string_basic.cpp
   :fragment: flatten

Here is the output of our small program.

.. code-block:: console

   After creating the Journaled String:
   Host: thisisahostsequence
   Journal: thisisahostsequence
   Nodes: JournalEntries({segmentSource=1, virtualPosition=0, physicalPosition=0, length=19})

   After modifying the Journaled String:
   Host: thisisahostsequence
   Journal: thisisamodifiedhost
   Nodes: JournalEntries({segmentSource=1, virtualPosition=0, physicalPosition=0, length=7}, {segmentSource=2, virtualPosition=7, physicalPosition=0, length=8}, {segmentSource=1, virtualPosition=15, physicalPosition=7, length=4})

   After flatten the Journaled String:
   Host: thisisamodifiedhost
   Journal: thisisamodifiedhost
   Nodes: JournalEntries({segmentSource=1, virtualPosition=0, physicalPosition=0, length=19})

.. important::

   Be careful when using the :dox:`JournaledString#flatten` function as it modifies the underlying host sequence.
   This might affect other journaled sequences that share the same host sequence.
   This becomes important especially when working with Journal Sets where a whole set of sequences is journaled to the same reference.

Journaled Set
-------------

The :dox:`JournaledSet` is a specialization of the :dox:`StringSet` which can be used exactly as such but also provides some additional functions optimized to work with :dox:`JournaledString JournaledStrings`.
The general interface is equal to the interface of the StringSet.
But it also provides some interfaces specialized for the use of Journaled Strings.
One of these interfaces is the :dox:`JournaledSet#join` function which journales a contained Journaled String to the previously set global reference.
The following code snippet demonstrates the usage of the Journal Set and how to join a sequence to the previously set reference sequence.

As usual we include the necessary headers.
We need the header ``<seqan/journal_set.h>`` to get access to the Journal Set.
Then we define a type for journaled sequences.
After that we define our Journal Set.
The Journal Set is a specialization of the :dox:`OwnerStringSet Owner` concept of StringSets and is defined as ``StringSet<TJournalString, Owner<JournaledSet> >``.

.. includefrags:: extras/demos/tutorial/data_journaling/example_join.cpp
   :fragment: main

In the subsequent steps we want to set a reference sequence to the Journal Set and add some sequences to it.
We can set the reference sequence by using the function :dox:`JournaledSet#setHost`.
This function stores only a pointer to the given sequence.
In some cases it might be necessary to copy the reference sequence instead.
For this purpose you can use the function :dox:`JournaledSet#createHost`.

.. includefrags:: extras/demos/tutorial/data_journaling/example_join.cpp
   :fragment: init

Just adding sequences to the Journal Set does not automatically journal them to the global reference sequence of the set.
One can explicitly trigger this process using the function :dox:`JournaledSet#join`.
This function takes as parameters the Journal Set and the position of the contained Journaled String which is to be journaled to the reference sequence.
Thus, the programmer is free in the decision which sequence has to be journaled and which not.
Furthermore, we can use the :dox:`JoinConfig` object to specify the method that shall be used for the journaling process.

.. includefrags:: extras/demos/tutorial/data_journaling/example_join.cpp
   :fragment: join

.. tip::

    Configuration of the Join Methods

    The :dox:`JoinConfig` object differentiates two methods in general and each method further differs in the used strategy.
    The two methods are the :dox:`GlobalAlign` and the :dox:`GlobalChain` method.
    They differ in the approach of computing the alignment which is necessary to construct the journal.
    The first method uses a global alignment algorithm while the second one uses an anchor approach in which first exact seeds are found using a q-gram index and after that the optimal chain between the identified anchors is computed.
    For each method the user can specify a different strategy.
    The first strategy is triggered by using :dox:`JoinStrategiesTags JournaledManhatten`.
    This means for the the GlobalAlign method, that the complete sequence is inserted and the complete reference is deleted, while for the GlobalChain methods this means that the gaps between the anchors are connected through the Manhatten distance.
    The second strategy is specified using the :dox:`JoinStrategiesTags JournaledCompact` tag. It computes the most compact form of a journal by menas of memory requirements.

Here is the output of the program.

.. code-block:: console

   Reference: DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE
   Journaled Sequence 0: DPKKPRGKMVNSPPAFFVQTSREEHKKKHPDASVFSKKCSERWKTMSAKEKGKFEDMAKARYEREMKTTYIPKGETYIPPKGE
   Journaled Sequence 1: DPHHPPKPRGKMVNSPPAFFVQTSREEHKPDASVFSKKCSERRMPNHHTMSAKEKGKFEDMAKARYEREMKTTYIPKGETYIPPKGE
   Journaled Sequence 2: DPKKPRGKMSSYAFFVQTSREEHKKKHPKKCDEFSKKCSERWKTMSAKEKGKFEDARYEREMKTYIPPKGE

Implementing an Online-Search
-----------------------------

Now we have all foundations laid down to implement the online-search algorithm.
Let us begin with the first assignment where we read in some sequences and use the currently learned things about the Journal Set.

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Review, Application

   Objective
     Download the fasta file :download:`sequences.fasta <sequences.fasta>` which contains some DNA sequences.
     Write a method called ``loadAndJoin`` that gets a Journal Set and a stream file pointing to the downloaded fasta file.
     The method reads in the sequences one after another using SeqAn's :dox:`RecordReader`.
     The first read sequences is set to the reference sequence.
     All following sequences are first appended to the StringSet and afterwards joined to the StringSet using a global alignment strategy and the most compact form.

   Hints
     .. container:: foldable

        You can start using the following code snippet. Replace the path of the iostream such that it points to your path and fill in the missing parts ``A``, ``B`` and ``C`` in the function ``loadAndJoin`` (Altogether, you will need 4 lines of code.).

       .. includefrags:: extras/demos/tutorial/data_journaling/example_online_search_assignment1_hint.cpp
          :fragment: main

   Solution
     .. container:: foldable

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment1.cpp
           :fragment: main

Now we have loaded and journaled our sequences and we use the minimal possible memory requirements for our sequences.
Let's continue and implement the exact online-search on the Journal Set.
For this purpose we write a function called ``searchPattern`` which takes a StringSet of ``String<int>`` which we use to store each hit for each sequence in the Journal Set including the reference sequence.
First we have to check whether the reference sequence is set.
If not we abort the search since we cannot guarantee a correct search when the reference is not set.
We also have to clear our ``hitSet`` to ensure there remain no phantom hits from previous searches.
Then we resize it to the number of contained Journaled Strings plus an additional space for the global reference sequence.

.. includefrags:: extras/demos/tutorial/data_journaling/example_online_search.cpp
   :fragment: searchPatternPart1

Before we can search for the pattern in the Journaled Strings we have to find all occurrences in the reference sequence.
Therefore we call the function ``findPatternInReference`` which takes a ``String<int>`` which we use to store the hits, the global reference and the pattern.

.. includefrags:: extras/demos/tutorial/data_journaling/example_online_search.cpp
   :fragment: searchPatternPart2

After that we implement the body to search for occurrences in the Journaled Strings.
Therefore we use ``for``-loop to iterate over all contained sequences and call for each sequence the function ``findPatternInJournalString``.
The function gets as parameters the corresponding ``String<int>`` from the ``hitSet``, the journaled sequence the pattern and the identified hits in the reference sequence.

.. includefrags:: extras/demos/tutorial/data_journaling/example_online_search.cpp
   :fragment: searchPatternPart3

So far our program won't compile. We have to first implement the both functions ``findPatternInReference`` and ``findPatternInJournalString``.

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Implement the function ``findPatternInReference`` using a brute force pattern search algorithm.
     Store all found hits in the passed ``hits`` variable.
     Print all occurrences in the end of the ``main`` function.

   Hints
     .. container:: foldable
        The following code snippet provides you with the backbone for the search algorithm. Fill in the missing parts ``[A]``, ``[B]``, ``[C]`` and ``[D]``.

        .. includefrags:: extras/demos/tutorial/data_journaling/example_online_search_assignment2_hint.cpp
           :fragment: findPatternInReferenceHint

   Solution
     .. container:: foldable

        Here is the solution for this function.
        Click **more...** below, to see a complete solution of everything we have done so far.

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment2.cpp
           :fragment: findPatternInReference

        .. container:: foldable

           Include the necessary headers.

           .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment2.cpp
              :fragment: include

           Implementation of the `findPatternInReference` function.

           .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment2.cpp
              :fragment: findPatternInReference

           Implementation of the `searchPattern` function. Note that we haven't implemented the function `findPatternInJournalString` yet.

           .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment2.cpp
              :fragment: searchPattern

           Implementation of the `loadAndJoin` function.

           .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment2.cpp
              :fragment: loadAndJoin

           Implementation of the `main` function.

           .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment2.cpp
              :fragment: main

           Printing the hits of the reference sequence.

           .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment2.cpp
              :fragment: printResult

           And here is the result.

           .. code-block:: console

              Search for: GTGGT:
              Hit in reference  at 311: GTGGT	644: GTGGT

We know can search for all occurrences of a pattern in the reference sequence.
Now we can try to find all occurrences in the journaled sequences.
Therefore we implement the function ``findPatternInJournalString``.
Our function gets the variable ``hitsTarget`` which stores the hits found in the JournaledString.
It gets the search text and the pattern and finally the hits detected in the reference sequence.
Instead of searching each position in the Journaled String, we only search in areas that are new to the search.
This involves all inserted parts and all parts where the pattern crosses a border to another node.
So instead of iterating over each position we iterate over the nodes of the Journaled String.
To do so we have to determine the type of the data structure that is used by the Journaled String to manage the nodes.
We can use the metafunction :dox:`JournaledString#JournalType` for this task.
Afterwards we define an Iterator over the so called ``TJournalEntries`` data structure.

Again we check first whether the pattern fits into our sequence.

.. includefrags:: extras/demos/tutorial/data_journaling/example_online_search.cpp
   :fragment: findPatternInJournalStringPart1

We then iterate over all nodes beginning from the first until we have reached the node in which the pattern reaches the end of the Journaled Sequence.
The function findInJournalEntries helps us to find the corresponding node.
We increment the position of the iterator by one such that it points behind the last element which is included by the search.

.. includefrags:: extras/demos/tutorial/data_journaling/example_online_search.cpp
   :fragment: findPatternInJournalStringPart2

Now we search in each node until we have reached the end.
For each node we first check the type of the journaled operation.
If we are in an "original" node we call the function ``_findInOriginalNode``.
If we are in a "patch" node we call the function ``_findInPatchNode`` and in the end we call the function ``_searchAtBorder`` which is called for each node type and scans all possible hits at the border of a node.

.. includefrags:: extras/demos/tutorial/data_journaling/example_online_search.cpp
   :fragment: findPatternInJournalStringPart3

Let us begin with the implementation of the function ``_findInOriginalNode``.
In this function we exploit the journaling concept such that we can speed up the search algorithm from :math:`\mathcal{O}(p \cdot n)` to :math:`\mathcal{O}(\log_2(k))`, where :math:`p` is the length of the pattern, :math:`n` is the length of the search text, and ``k`` is the number of hits identified in the reference sequence.
We need at most :math:`\log_2(k)` comparisons to find the first hit which occurred in the reference sequence that also occurs in the current original node.

Assignment 3
^^^^^^^^^^^^

.. container:: assignment

   Type
     Transfer

   Objective
     Implement the function ``_findInOriginalNode``, which identifies all shared hits between the current ``original`` node and the corresponding part in the reference sequence.
     Note you do not need to scan all positions again.
     In the end print all occurrences to the console.

   Hints
     .. container:: foldable
       The following code snippet provides you with the backbone for this function.
       Fill in the missing parts ``[A]``, ``[B]``, ``[C]``, ``[D]`` and ``[E]``.

       Use the STL function `std::upper_bound <http://www.cplusplus.com/reference/algorithm/upper_bound/>`_ to conduct a binary search to find the first possible hit from the reference that is also represented by the current node.

       .. includefrags:: extras/demos/tutorial/data_journaling/example_online_search_assignment3_hint.cpp
          :fragment: findInOriginalNode

   Solution
     .. container:: foldable

       Here is the solution to this function.
       Click **more...** below, to see a complete solution of everything we have done so far.

       .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
          :fragment: findInOriginalNode

       .. container:: foldable

          Include the necessary headers.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
             :fragment: include

          We know implement the method to search for hits in an original node.
          We only need to check if the current node covers a region of the reference in which we've found a hit.
          We use the function `std::upper_bound <http://www.cplusplus.com/reference/algorithm/upper_bound/>`_ to find the first element that is greater than the current physical position.
          Since, we've found an upper bound we have to check additionally if there exists a previous hit that lies directly on the physical begin position of our current node.
          We then include all hits that fit into this current node until we have found the first position where the pattern would cross the border of this node or we have reached the end of the reference hit set.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
            :fragment: findInOriginalNode

          Implementing the backbone to search for a pattern in the reference string.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
            :fragment: findPatternInJournalString

          Implementing the search within the reference sequence.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
            :fragment: findPatternInReference

          Implementing the backbone of the search.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
            :fragment: searchPattern

          Implement the `laodAndJoin` method.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
            :fragment: loadAndJoin

          Implementing the main method.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
            :fragment: main

          Printing the hits of the reference sequence.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
            :fragment: printResultReference

          Printing the hits of the journaled sequences.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment3.cpp
            :fragment: printResultJournalSequence

          And here is the result.

          .. code-block:: console

             Search for: GTGGT:
             Hit in reference  at 311: GTGGT	644: GTGGT
             Hit in sequence 0 at 312: GTGGT
             Hit in sequence 1 at 308: GTGGT
             Hit in sequence 2 at 311: GTGGT
             Hit in sequence 3 at 327: GTGGT
             Hit in sequence 4 at 317: GTGGT
             Hit in sequence 5 at 320: GTGGT

We are almost at the end of our online-search algorithm.
Let's now implement the method ``_findInPatchNode``.
We basically had this already implemented when we wrote the search function for the reference. Let's recall this part together.

First we write the body of our function and define now an Iterator over the Journaled String.

.. includefrags:: extras/demos/tutorial/data_journaling/example_online_search.cpp
   :fragment: findInPatchNodePart1

We know specify the range in which we are searching for the pattern.
This range starts at the current physical position of the insertion buffer and ends at the last position of this node where the pattern completely fits.

.. includefrags:: extras/demos/tutorial/data_journaling/example_online_search.cpp
   :fragment: findInPatchNodePart2

We need to use a second temporary iterator which is used to compare the current value with the pattern.
If all positions matches then we report a match at the current virtual position.

.. includefrags:: extras/demos/tutorial/data_journaling/example_online_search.cpp
   :fragment: findInPatchNodePart3

To ensure that we are not missing any hits we also have to scan the regions where the pattern is leaving the current node.
You can solve this problem in the next assignment.

Assignment 4
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     Implement the last function ``_searchAtBorder``, which identifies all hits that cross the border of the current node to the next node.

   Hints
     .. container:: foldable

       The following code snippet provides you with the backbone for this function.
       Fill in the missing parts ``[A]``, ``[B]``, ``[C]`` and ``[D]``.

       .. includefrags:: extras/demos/tutorial/data_journaling/example_online_search_assignment4_hint.cpp
          :fragment: searchAtBorder

   Solution
     .. container:: foldable

        Here is the solution to this function.
        Click **more...** below, to see a complete solution of everything we have done so far.

       .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
          :fragment: searchAtBorder

       .. container:: foldable

          Include the necessary headers.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: include

          Search at the border the current node for the pattern.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: searchAtBorder

          Search for the pattern in the insertion region covered by the current node.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: findInPatchNode

          Check if hit was reported for this region in the reference sequence.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: findInOriginalNode

          Implementing the backbone of the search for the Journaled String.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: findPatternInJournalStringPart1

          Implementing the search for the reference sequence.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: findPatternInReference

          The backbone of the search method.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: searchPatternPart1

          Loading and joining the sequences.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: loadAndJoin

          Implementing the main function.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: main

          Reporting the identified hits.

          .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_assignment4.cpp
            :fragment: printResult

          And here is the result.

          .. code-block:: console

             Search for: GTGGT:
             Hit in reference  at 311: GTGGT	644: GTGGT
             Hit in sequence 0 at 151: GTGGT	312: GTGGT
             Hit in sequence 1 at 308: GTGGT
             Hit in sequence 2 at 311: GTGGT	507: GTGGT
             Hit in sequence 3 at 327: GTGGT
             Hit in sequence 4 at 307: GTGGT	312: GTGGT	317: GTGGT
             Hit in sequence 5 at 0: GTGGT	320: GTGGT	986: GTGGT

Congratulations!
You have just implemented a cool online-search which is speed up by exploiting the parallelism given by the data set.
And here is the final result.

.. code-block:: console

   Search for: GTGGT:
   Hit in reference  at 311: GTGGT 644: GTGGT
   Hit in sequence 0 at 151: GTGGT 312: GTGGT
   Hit in sequence 1 at 308: GTGGT
   Hit in sequence 2 at 311: GTGGT 507: GTGGT
   Hit in sequence 3 at 327: GTGGT
   Hit in sequence 4 at 307: GTGGT 312: GTGGT  317: GTGGT
   Hit in sequence 5 at 0: GTGGT   320: GTGGT  986: GTGGT

Assignment 5
""""""""""""

.. container:: assignment

   Type
     Transfer

   Objective
     Try to replace the brute force versions using using SeqAn's :dox:`Finder` and :dox:`Pattern` concept.
     You can find additional material to this topic in the :ref:`tutorial-pattern-matching` Tutorial.

   Solution
     .. container:: foldable

        Now we want to replace the brute force methods with some cool pattern matching algorithms.
        Therefore we include the header `<seqan/finder.h>`.

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_finder.cpp
          :fragment: include

        Now we can use the :dox:`Finder` interface of SeqAn.
        One cool thing of the usage of the Finder class is that we don't have to check for the borders anymore.
        This will do the Finder for us.
        We only have to specify the correct infix over which the Finder should iterate to find the pattern.
        We first compute the positions that enclose the search region.
        Afterwards, we get an infix for this region and pass it to the Finder's constructor.
        We also have to define the :dox:`Pattern` object which gets the pattern we are searching for.
        Then we can simply call the function :dox:`Finder#find` until we there is no more match.
        Be careful when storing the position that the Finder is returning.
        We have to recompute the correct virtual position since we used an infix of the original search text.

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_finder.cpp
          :fragment: searchAtBorder

        So the biggest change is done.
        We simply repeat the changes from above and watch to get the correct virtual position.

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_finder.cpp
          :fragment: findInPatchNodePart1

        Of course we don't need to change anything for the original nodes.

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_finder.cpp
          :fragment: findInOriginalNode

        Also the function `findPatternInJournalString` remains the same.

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_finder.cpp
          :fragment: findPatternInJournalString

        We will switch to the Finder concept for the function `findPatternInReference` too.
        This is done quickly, since we have the basis already laid down in the previous functions.

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_finder.cpp
          :fragment: findPatternInReference

        From here on, we don't have to change anything.

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_finder.cpp
          :fragment: searchPatternPart1

        We write the same main body ...

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_finder.cpp
          :fragment: laodAndJoin

        and finally print the results.

        .. includefrags:: extras/demos/tutorial/data_journaling/solution_online_search_finder.cpp
          :fragment: main

        And here is the result using the Finder and Pattern concept of SeqAn.

        .. code-block:: console

           Search for: GTGGT:
           Hit in reference  at 311: GTGGT	644: GTGGT
           Hit in sequence 0 at 151: GTGGT	312: GTGGT
           Hit in sequence 1 at 308: GTGGT
           Hit in sequence 2 at 311: GTGGT	507: GTGGT
           Hit in sequence 3 at 327: GTGGT
           Hit in sequence 4 at 307: GTGGT	312: GTGGT	317: GTGGT
           Hit in sequence 5 at 0: GTGGT	320: GTGGT	986: GTGGT
