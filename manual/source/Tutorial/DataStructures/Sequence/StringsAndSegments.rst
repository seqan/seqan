.. sidebar:: ToC

    .. contents::

.. _tutorial-datastructures-sequences-strings-and-segments:

Strings and Segments
====================

Learning Objective
  You will learn about the SeqAn sequence concept and its main class :dox:`String` as well as the class :dox:`Segment`.
  After completing this tutorial, you will be able to use important functionalities of sequences in SeqAn and you will be ready to continue with the more specific tutorials, e.g. :ref:`tutorial-datastructures-alignment`, or :ref:`tutorial-algorithms-alignment-pairwise-sequence-alignment`.

Difficulty
  Very basic

Duration
  45 min

Prerequisites
  Basic C or C++ knowledge, the :ref:`tutorial-getting-started-first-steps-in-seqan` tutorial helps.

Sequences are the core concept of SeqAn.
A sequence is a container that stores an ordered list of values.
In SeqAn, there are three kinds of sequences: Strings, Sequence Adaptions and Segments.

The :dox:`String` class is one of the most fundamental classes in SeqAn.
It is designed as a generic data structure that can be instantiated for all kinds of values, both simple (e.g. ``char``, :dox:`Dna`, :dox:`AminoAcid`) and non-simple value types (e.g. :dox:`Tuple`, :dox:`String`).
With sequence adaptions, SeqAn offers an interface for accessing data types that are not part of SeqAn, namely standard library strings and c-style char arrays.
Thus those built-in types can be handled in a similar way as SeqAn strings, for example with the :dox:`ContainerConcept#length` function.
:dox:`Segment Segments` are contiguous subsequences that represent parts of other sequences.

This tutorial will deal with the SeqAn sequence classes :dox:`String` and :dox:`Segment`.

Strings
-------

In this section, we will have a detailed look at the SeqAn class :dox:`String`.
You will learn how to build and expand strings as well as how to compare and convert them.

Defining Strings
^^^^^^^^^^^^^^^^

Let's first have a look at a simple example on how to define a :dox:`String`.
The type of the contained value is specified by the first template argument, e.g. ``char`` or ``int``.

.. includefrags:: demos/tutorial/sequences/base.cpp
    :fragment: string_example

To fill the string with contents, we can simply assign a string literal to the created variable:

.. includefrags:: demos/tutorial/sequences/base.cpp
    :fragment: simple_string_construction

Any type that provides a default constructor, a copy constructor and an assignment operator can be used as the alphabet / contained type of a :dox:`String`.
This includes the C++ `POD types <https://isocpp.org/wiki/faq/intrinsic-types#pod-types>`_, e.g. ``char``, ``int``, ``double`` etc., or even more complex types complex types, such as :dox:`String Strings`.

.. includefrags:: demos/tutorial/sequences/base.cpp
    :fragment: string_of_strings_example

.. hint::

   Nested Sequences (aka "Strings of Strings")

   A collection of sequences can either be stored in a sequence of sequences, for example in a ``String< String<char> >``, or in a :dox:`StringSet`.
   The latter one allows for more auxiliary functionalities to improve the efficiency of working with large sequence collections.
   You can learn more about it in the tutorial :ref:`tutorial-datastructures-sequences-string-sets`.

SeqAn also provides the following types that are useful in bioinformatics: :dox:`AminoAcid`, :dox:`Dna`, :dox:`Dna5`, :dox:`DnaQ`, :dox:`Dna5Q`, :dox:`Finite`, :dox:`Iupac`, :dox:`Rna`, :dox:`Rna5`.
You can find detailed information in the tutorial :ref:`tutorial-datastructures-sequences-alphabets`.

.. includefrags:: demos/tutorial/sequences/base.cpp
    :fragment: special_types_example

For commonly used string parameterizations, SeqAn has a range of shortcuts implemented, e.g. :dox:`DnaString`, :dox:`RnaString` and :dox:`Peptide`.

.. includefrags:: demos/tutorial/sequences/base.cpp
    :fragment: shortcuts_example

Working with Strings
^^^^^^^^^^^^^^^^^^^^

The SeqAn String implementation provides the common C++ operators that you know already from the `vector <http://en.cppreference.com/w/cpp/container/vector>`_ class of the STL.
For example:

.. includefrags:: demos/tutorial/sequences/example_functionality1.cpp
    :fragment: main

.. includefrags:: demos/tutorial/sequences/example_functionality1.cpp.stdout

Each sequence object has a capacity, i.e. the maximum length of a sequence that can be stored in this object.
While some sequence types have a fixed capacity, the capacity of other sequence classes like :dox:`AllocString Alloc String` or ``std::basic_string`` can be changed at runtime.
The capacity can be set explicitly by functions such as :dox:`String#reserve` or :dox:`StringConcept#resize`.
It can also be set implicitly by functions like :dox:`StringConcept#append` or :dox:`StringConcept#replace`, if the operation's result exceeds the length of the target string.

In the following example we create a :dox:`String` of :dox:`Dna5String`. We first set the new length of the container with :dox:`StringConcept#resize` to two elements.
After assigning two elements we append one more element with :dox:`StringConcept#appendValue`.
In the last step the capacity is implicitly changed.

.. includefrags:: demos/tutorial/sequences/example_functionality2.cpp
    :fragment: main

Using the function :dox:`ContainerConcept#length`, we can now get the length of our strings, e.g.:

.. includefrags:: demos/tutorial/sequences/example_functionality2.cpp
    :fragment: print

.. includefrags:: demos/tutorial/sequences/example_functionality2.cpp.stdout

To empty a :dox:`String`, the function :dox:`StringConcept#clear` resets the object.

.. includefrags:: demos/tutorial/sequences/example_functionality2.cpp
    :fragment: clear

SeqAn offers a range of other functions for the work with the :dox:`String` class, e.g. :dox:`AssignableConcept#assign`, :dox:`RandomAccessContainerConcept#assignValue`, :dox:`ContainerConcept#empty`, etc.
The full list of functions you can find in the documentation :dox:`String`.

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     In the following assignment, you will write a small function that builds the reverse complement of a given string.
     Copy the code below and add the following functionalities:

     #. Use the ``resize`` function to ``resize`` the ``revComplGenome`` variable.
     #. Using the ``getRevCompl`` function, get the reverse complement for every nucleotide ``genome`` and store it in reverse order ``revComplGenome``.
     #. Print out the original genome and the reverse complement.

        .. includefrags:: demos/tutorial/sequences/assignment_1_solution.cpp
           :fragment: top

        .. code-block:: cpp

           // Your code snippet here


        .. includefrags:: demos/tutorial/sequences/assignment_1_solution.cpp
           :fragment: bottom

   Hints
     Remember that the last element in ``genome`` is stored at position ``length(genome) - 1``.

   Solution
     Click *more...* to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/sequences/assignment_1_solution.cpp
            :fragment: full

        Your output should look like this:

        .. includefrags:: demos/tutorial/sequences/assignment_1_solution.cpp.stdout

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     In this assignment, you will do some simple string building tasks, and write a simple alignment of the given reads and chromosomes.
     Use the given code template to solve these subtasks:

   #. Assume we have mapped the reads to the positions 7, 100, 172, and 272 in 'chr1'.
      Store these positions in another string 'alignPosList'.
   #. Build another String bsChr1 as a copy of chr1, and exchange every 'C' with a 'T', as in a bisulfite treated genome.
   #. Print alignments of the reads and chr1 (or bschr1) using the function ``printAlign`` and the string ``alignPosList``.

    .. includefrags:: demos/tutorial/sequences/assignment_2_solution.cpp
          :fragment: one

    .. code-block:: cpp

        // Your code snippet here for 1.+2.

    .. includefrags:: demos/tutorial/sequences/assignment_2_solution.cpp
          :fragment: two

    .. code-block:: cpp

        // Your code snippet here for 3.

    .. includefrags:: demos/tutorial/sequences/assignment_2_solution.cpp
          :fragment: three

    .. code-block:: cpp

        // Your code snippet here for 3.

    .. includefrags:: demos/tutorial/sequences/assignment_2_solution.cpp
          :fragment: four

   Hints
     You have to create a copy of the fragment in chr1 (bsChr1) that is aligned to the read.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/sequences/assignment_2_solution.cpp
          :fragment: full

        .. includefrags:: demos/tutorial/sequences/assignment_2_solution.cpp.stdout

Comparisons
^^^^^^^^^^^

Two sequences can be lexicographically **compared** using standard operators such as ``<`` or ``>=``.

.. includefrags:: demos/tutorial/sequences/example_comparisons.cpp
    :fragment: main

.. includefrags:: demos/tutorial/sequences/example_comparisons.cpp.stdout

Each comparison involves a scan of the two sequences for searching the first mismatch between the strings.
This could be costly if the two sequences share a long common prefix.
Suppose we want to branch in a program depending on whether ``a < b``, ``a == b``, or ``a > b``.

.. includefrags:: demos/tutorial/sequences/example_comparisons.cpp
    :fragment: first

In this case, although only one scan would be enough to decide what case is to be applied, each operator ``>`` and ``<`` performs a new comparison.
SeqAn offers the class :dox:`Lexical` to avoid unnecessary sequence scans.
Lexicals can store the result of a comparison, for example:

.. includefrags:: demos/tutorial/sequences/example_comparisons.cpp
    :fragment: second

Conversions
^^^^^^^^^^^

A sequence of type A values can be converted into a sequence of type B values, if A can be converted into B.
SeqAn offers different conversion alternatives.

**Copy conversion.**
The source sequence is copied into the target sequence.
This can be done by assignment (``operator=``) or using the function :dox:`AssignableConcept#assign`.

.. includefrags:: demos/tutorial/sequences/example_conversions_copy.cpp
    :fragment: main

.. includefrags:: demos/tutorial/sequences/example_conversions_copy.cpp.stdout

**Move conversion.**
If the source sequence is not needed any more after the conversion, it is always advisable to use :dox:`AssignableConcept#move` instead of :dox:`AssignableConcept#assign`.
The function :dox:`AssignableConcept#move` does not make a copy but can reuse the source sequence storage.
In some cases, :dox:`AssignableConcept#move` can also perform an in-place conversion.

.. includefrags:: demos/tutorial/sequences/example_conversions_move.cpp
    :fragment: main

.. includefrags:: demos/tutorial/sequences/example_conversions_move.cpp.stdout

Assignment 3
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     In this assignment you will sort nucleotides.
     Copy the code below. Adjust the code such that all nucleotides, which are lexicographically smaller than a Dna5 ``'G'`` are stored in a list ``lesser``, while all nucleotides which are greater, should be stored in a list ``greater``.
     Print out the final lists.

     .. includefrags:: demos/tutorial/sequences/assignment_3.cpp

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/sequences/assignment_3_solution.cpp

        .. includefrags:: demos/tutorial/sequences/assignment_3_solution.cpp.stdout

Assignment 4
""""""""""""

.. container:: assignment

   Type
     Transfer

   Objective
     In this task you will compare whole sequences.
     Reuse the code from above. Instead of a ``String<Dna5>`` we will now deal with a ``String<Dna5String>``.
     Build a string which contains the Dna5Strings "ATATANGCGT", "AAGCATGANT" and "TGAAANTGAC".
     Now check for all elements of the container, if they are lexicographically smaller or bigger than the  given subject sequence "GATGCATGAT" and append them to a appropriate list.
     Print out the final lists.

   Hints
     Try to avoid unnecessary sequence scans.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/sequences/assignment_4_solution.cpp

        .. includefrags:: demos/tutorial/sequences/assignment_4_solution.cpp.stdout

.. _tutorial-datastructures-sequences-strings-and-segments-iterators:

Iteration
^^^^^^^^^

Very often you will be required to iterate over your string to either retrieve what's stored in the string or to write something at a specific position.
For this purpose SeqAn provides Iterators for all container types.
The metafunction :dox:`ContainerConcept#Iterator` can be used to determine the appropriate iterator type for a given a container.

An iterator always points to one value of the container.
The operator :dox:`IteratorAssociatedTypesConcept#operator*` can be used to access this value by reference.
Functions like :dox:`InputIteratorConcept#operator++(prefix)` or :dox:`BidirectionalIteratorConcept#operator--(prefix)` can be used to move the iterator to other values within the container.

The functions :dox:`ContainerConcept#begin` and :dox:`ContainerConcept#end`, applied to a container, return iterators to the begin and to the end of the container.
Note that similar to C++ standard library iterators, the iterator returned by :dox:`ContainerConcept#end` does not point to the last value of the container but to the position behind the last one.
If the container is empty then ``end() == begin()``.

The following code prints out a sequence and demonstrates how to iterate over a string.

.. includefrags:: demos/tutorial/iterators/base.cpp
    :fragment: use-case

.. includefrags:: demos/tutorial/iterators/base.cpp.stdout
    :fragment: use-case


Different Iterator Types
""""""""""""""""""""""""

Some containers offer several kinds of iterators, which can be selected by an optional template parameter of the Iterator class.
For example, the tag :dox:`ContainerIteratorTags#Standard` can be used to get an iterator type that resembles the C++ standard random access iterator.
For containers there is also a second variant available, the so called :dox:`ContainerIteratorTags#Rooted` iterator.
The rooted iterator knows its container by pointing back to it.
This gives us a nice interface to access member functions of the underlying container while operating on a rooted iterator.
The construction of an iterator in SeqAn, e.g. for a :dox:`DnaString Dna String`, could look like the following:

.. includefrags:: demos/tutorial/iterators/base.cpp
    :fragment: construction

.. tip::

   The default iterator implementation is :dox:`ContainerIteratorTags#Standard`.
   Rooted iterators offer some convenience interfaces for the user.
   They offer additional functions like :dox:`RootedIteratorConcept#container` for determining the container on which the iterator works, and they simplify the interface for other functions like :dox:`RootedIteratorConcept#atEnd`.
   Moreover, rooted iterators may change the container’s length or capacity, which makes it possible to implement a more intuitive variant of a remove algorithm.

   While rooted iterators can usually be converted into standard iterators, it is not always possible to convert standard iterators back into rooted iterators, since standard iterators may lack the information about the container they work on.
   Therefore, many functions that return iterators like :dox:`ContainerConcept#begin` or :dox:`ContainerConcept#end` return rooted iterators instead of standard iterators; this way, they can be used to set both rooted and standard iterator variables.
   Alternatively it is possible to specify the returned iterator type explicitly by passing the iterator kind as a tag argument, e.g. ``begin(str, Standard())``.

Assignment 5
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     Copy the code below, which replaces all N's of a given :dox:`String` with A's.
     Adjust the code to use iterators to traverse the container.
     Use the :dox:`ContainerIteratorTags#Standard` iterator.

     .. includefrags:: demos/tutorial/iterators/assignment_1.cpp

    Solution

      Click **more...** to see the solution.

      .. container:: foldable

         .. includefrags:: demos/tutorial/iterators/assignment_1_solution.cpp

Assignment 6
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Use the code from above and change the :dox:`ContainerIteratorTags#Standard` to a :dox:`ContainerIteratorTags#Rooted` iterator.
     Try to shorten the code wherever possible.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/iterators/assignment_2_solution.cpp

String Allocation Strategies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each sequence object has a capacity, i.e. the reserved space for this object.
The capacity can be set explicitly by functions such as :dox:`String#reserve` or :dox:`StringConcept#resize`.
It can also bet set implicitly by functions like :dox:`ContainerConcept#append`, :dox:`AssignableConcept#assign`, :dox:`StringConcept#insert` or :dox:`StringConcept#replace`, if the operation's result exceeds the length of the target sequence.

If the current capacity of a sequence is exceeded by chaining the length, we say that the sequence overflows.
There are several overflow strategies that determine what actually happens when a string should be expanded beyond its capacity.
The user can specify this for a function call by additionally handing over a tag.
If no overflow strategy is specified, a default overflow strategy is selected depending on the type of the sequence.

The following overflow strategies exist:

:dox:`OverflowStrategyTags#Exact`
  Expand the sequence exactly as far as needed. The capacity is only changed if the current capacity is not large enough.

:dox:`OverflowStrategyTags#Generous`
  Whenever the capacity is exceeded, the new capacity is chosen somewhat larger than currently needed.
  This way, the number of capacity changes is limited in a way that resizing the sequence only takes amortized constant time.

:dox:`OverflowStrategyTags#Limit`
  Instead of changing the capacity, the contents are limited to current capacity.
  All values that exceed the capacity are lost.

:dox:`OverflowStrategyTags#Insist`
  No capacity check is performed, so the user has to ensure that the container's capacity is large enough.

The next example illustrates how the different strategies could be used:

.. includefrags:: demos/tutorial/sequences_in_depth/example_overflow.cpp
   :fragment: example

.. includefrags:: demos/tutorial/sequences_in_depth/example_overflow.cpp.stdout

Assignment 7
""""""""""""

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

         .. includefrags:: demos/tutorial/sequences_in_depth/assignment_exact_generous_solution.cpp

String Specializations
^^^^^^^^^^^^^^^^^^^^^^

The user can specify the kind of string that should be used in an optional second template argument of :dox:`String`.
The default string implementation is :dox:`AllocString Alloc String`.

.. includefrags:: demos/tutorial/sequences/base.cpp
    :fragment: sdefault_type

In most cases, the implementation :dox:`AllocString Alloc String` (the default when using a ``String<T>``) is the best choice.
Exceptions are when you want to process extremely large strings that are a bit larger than the available memory (consider :dox:`AllocString Alloc String`) or much larger so most of them are stored on the hard disk and only parts of them are loaded in main memory (consider :dox:`ExternalString External String`).
The following list describes in detail the different specializations:

Specialization :dox:`AllocString Alloc String`
  * **Description**
    Expandable string that is stored on the heap.
  * **Applications**
    The default string implementation that can be used for general purposes.
  * **Limitations**
    Changing the :dox:`StringConcept#capacity` can be very costly since all values must be copied.

Specialization :dox:`ArrayString Array String`
  * **Description**
    Fast but non-expandable string. Fast storing of fixed-size sequences.
  * **Limitations**
    :dox:`StringConcept#capacity Capacity` must already be known at compile time. Not suitable for storing large sequences.

Specialization :dox:`BlockString Block String`
  * **Description**
    String that stores its sequence characters in blocks.
  * **Applications**
    The :dox:`StringConcept#capacity` of the string can quickly be increased. Good choice for growing strings or stacks.
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

Specialization :dox:`JournaledString Journaled String`
  * **Description**
    String that stores differences to an underlying text rather than applying them directly.
  * **Applications**
    Suitable for efficiently storing similar strings, if their differences to an underlying reference sequence are known.
  * **LimitationsApplications**
    Slower than other string classes, due to logarithmic penalty for random accesses.

Specialization :dox:`CStyleString CStyle String`
  * **Description**
    Allows adaption of strings to C-style strings.
  * **Applications**
    Used for transforming other String classes into C-style strings (i.e. null terminated char arrays). Useful for calling functions of C-libraries.
  * **Limitations**
    Only sensible if value type is ``char`` or ``wchar_t``.

.. includefrags:: demos/tutorial/sequences_in_depth/base.cpp
      :fragment: type_examples

.. includefrags:: demos/tutorial/sequences/base.cpp
    :fragment: external_string_spec

.. tip::

   String Simplify Memory Management

   One advantage of using Strings is that the user does not need to reserve memory manually with **new** and does not need **delete** to free memory.
   Instead, those operations are automatically handled by the :dox:`String` class.

   .. includefrags:: demos/tutorial/sequences/base.cpp
        :fragment: initialization_example

Segments
--------

The following section will introduce you into the :dox:`Segment` class of SeqAn.

:dox:`Segment Segments` are contiguous subsequences that represent parts of other sequences.
Therefore, their functionality is similar to the :dox:`String` functionality.
In SeqAn, there are three kinds of segments: :dox:`InfixSegment`, :dox:`PrefixSegment`, and :dox:`SuffixSegment`.
The metafunctions :dox:`SegmentableConcept#Infix`, :dox:`SegmentableConcept#Prefix`, and :dox:`SegmentableConcept#Suffix`, respectively, return the appropriate segment data type for a given sequence type.

For prefixes, we use the function :dox:`SegmentableConcept#prefix` to build the prefix.
The first parameter is the sequence we build the prefix from, the second the **excluding** end position.
For :dox:`SegmentableConcept#infix`\ es, we have to provide both the including start and the excluding end position.
For :dox:`SegmentableConcept#suffix`\ es, the second parameter of the function denotes the including starting position of the suffix:

.. includefrags:: demos/tutorial/sequences/example_segments.cpp
    :fragment: main

.. includefrags:: demos/tutorial/sequences/example_segments.cpp.stdout


Segments store a pointer on the underlying sequence object, the *host*, and an start and/or end position, depending on the type of segment.
The segment is *not* a copy of the sequence segment.

.. warning::

   Please note that it is not possible anymore to change the underlying sequence by changing the segment.
   If you want to change the host sequence, you have to explicitly modify this.
   If you want to modify only the segment, you have to explicitly make a copy of the string.

Assignment 8
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     In this task you will use a segment to pass over an infix of a given sequence to a function without copying the corresponding fragment.
     Use the code given below.
     Lets assume that we have given a ``genome`` and a ``read`` sequence as well as the begin position of a given alignment.
     In the main function a fragment of the Dna5String ``genome`` is copied and passed together with the Dna5String ``read`` to a ``print`` function.
     Adjust the code to use an infix of the genome, instead of copying the corresponding fragment.

     .. includefrags:: demos/tutorial/sequences/assignment_5_solution.cpp
          :fragment: top


     .. includefrags:: demos/tutorial/sequences/base.cpp
          :fragment: assignment5_code_to_change

     .. includefrags:: demos/tutorial/sequences/assignment_5_solution.cpp
          :fragment: bottom

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/sequences/assignment_5_solution.cpp
            :fragment: full

        .. includefrags:: demos/tutorial/sequences/assignment_5_solution.cpp.stdout

Assignment 9
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Take the solution from the workshop assignment above and change it to use Segments for building the genome fragment.

   Hints
     Note that because ``printAlign`` uses templates, you don't have to change the function even though the type of ``genomeFragment`` is different.

   Solution
    Click **more...** to see the solution.

    .. container:: foldable

       .. includefrags:: demos/tutorial/sequences/assignment_6_solution.cpp

       .. includefrags:: demos/tutorial/sequences/assignment_6_solution.cpp.stdout
