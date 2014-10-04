.. sidebar:: ToC

   .. contents::


.. _tutorial-sequences-in-depth:

Sequences In-Depth
------------------

Learning Objective
  You will learn in detail how to optimize the usage of sequences dependent on your needs.

Difficulty
  Advanced

Duration
  20 min

Prerequisites
  :ref:`tutorial-sequences`

Sequences, particularly :dox:`String Strings`, are fundamental in SeqAn.
You learned already how to use the default implementation of strings and how to easily work with them.
In the most cases the default string specialization is well suited as well as the default behavior for capicity changes.
Nevertheless, sometimes you might want to change the default behavior for efficiency reasons and adjust it to your specific needs.

String Specializations
~~~~~~~~~~~~~~~~~~~~~~

In this section you will learn about the different string specializations and when to use them.

The user can specify the kind of string that should be used in an optional second template argument of :dox:`String`.

.. code-block:: cpp

   String<Dna>           dnaSeq1; // The default string implementation: Alloc
   String<Dna, Alloc<> > dnaSeq2; // The same as above

In most cases, the implementation :dox:`AllocString Alloc String` (the default when using a ``String<T>``) is the best choice.
Exceptions are when you want to process extremely large strings that are a bit larger than the available memory (consider :dox:`AllocString Alloc String`) or much larger so most of them are stored on the hard disk and only parts of them are loaded in main memory (consider :dox:`ExternalString External String`).

The following list describes in detail the different specializations:

Specialization :dox:`AllocString Alloc String`
  * **Description**
    Expandable string that is stored on the heap.
  * **Applications**
    The default string implementation that can be used for general purposes.
  * **Limitations**
    Changing the :dox:`SequenceConcept#capacity` can be very costly since all values must be copied.

Specialization :dox:`ArrayString Array String`
  * **Description**
    Fast but non-expandable string.Fast storing of fixed-size sequences.
  * **Limitations**
    :dox:`SequenceConcept#capacity Capacity` must already be known at compile time. Not suitable for storing large sequences.

Specialization :dox:`BlockString Block String`
  * **Description**
    String that stores its sequence characters in blocks.
  * **Applications**
    The :dox:`SequenceConcept#capacity` of the string can quickly be increased. Good choice for growing strings or stacks.
  * **Limitations**
    Iteration and random access to values is slightly slower than for :dox:`AllocString Alloc String`.

Specialization :dox:`PackedString Packed String`
  * **Description**
    A string that stores as many values in one machine word as possible.
  * **Applications**
    Suitable for storing large strings in memory.
  * **Limitations**
    Slower than other in-memory strings.

Specialization :dox:`ExternalString External String`
  * **Description**
    String that is stored in secondary memory.
  * **Applications**
    Suitable for storing very large strings (>2GB). Parts of the string are automatically loaded from secondary memory on demand.
  * **LimitationsApplications**
    Slower than other string classes.

Specialization :dox:`CStyleString CStyle String`
  * **Description**
    Allows adaption of strings to C-style strings.
  * **Applications**
    Used for transforming other String classes into C-style strings (i.e. null terminated char arrays). Useful for calling functions of C-libraries.
  * **Limitations**
    Only sensible if value type is ``char`` or ``wchar_t``.

.. code-block:: cpp

   // String with maximum length 100.
   String<char, Array<100> > myArrayString;
   // String that takes only 2 bits per nucleotide.
   String<Dna, Packed<> > myPackedString;

Overflow Strategies
~~~~~~~~~~~~~~~~~~~

The following section will describe how you can improve capacity changes for your sequences.

Each sequence object has a capacity, i.e. the reserved space for this object.
The capacity can be set explicitly by functions such as :dox:`String#reserve` or :dox:`SequenceConcept#resize`.
It can also bet set implicitly by functions like :dox:`ContainerConcept#append`, :dox:`AssignableConcept#assign`, :dox:`SequenceConcept#insert` or :dox:`SequenceConcept#replace`, if the operation's result exceeds the length of the target sequence.

If the current capacity of a sequence is exceeded by chaning the length, we say that the sequence overflows.
There are several overflow strategies that determine what actually happens when a string should be expanded beyond its capacity.
The user can specify this for a function call by additionally handing over a tag.
If no overflow strategy is specified, a default overflow strategy is selected depending on the type of the sequence.

The following overflow strategies exist:

:dox:`OverflowStrategyTags#Exact`
  Expand the sequence exactly as far as needed. The capacity is only changed if the current capacity is not large enough.

:dox:`OverflowStrategyTags#Generous`
  Whenever the capacity is exceeded, the new capacity is chosen somewhat larger than currently needed.
  This way, the number of capacity changes islimited in a way that resizing the sequence only takes amortized constant time.

:dox:`OverflowStrategyTags#Limit`
  Instead of changing the capacity, the contents are limited to current capacity.
  All values that exceed the capacity are lost.

:dox:`OverflowStrategyTags#Insist`
  No capacity check is performed, so the user has to ensure that the container's capacity is large enough.

The next example illustrates how the different strategies could be used:

.. code-block:: cpp

   String<Dna> dnaSeq;
   // Sets the capacity of dnaSeq to 5.
   resize(dnaSeq, 4, Exact());
   // Only "TATA" is assigned to dnaSeq, since dnaSeq is limited to 4.
   assign(str, "TATAGGGG", Limit());
   std::cout << dnaSeq << std::endl;
   // Use the default expansion strategy.
   append(dnaSeq, "GCGCGC");
   std::cout << dnaSeq << std::endl;

.. code-block:: console

   TATA
   TATAGCGCGC

Workshop Assignment 1
^^^^^^^^^^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Build a string of Dna (default specialization) and use the function ``appendValue`` to append a million times the nucleotide 'A'.
     Do it both using the overflow strategy ``Exact`` and ``Generous``.
     Measure the time for the two different strategies.

   Solution
      Click **more...** to see the solution.

      .. container:: foldable

         .. includefrags:: core/demos/tutorial/sequences_in_depth/assignment_exact_generous_solution.cpp
