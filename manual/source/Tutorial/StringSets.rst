.. sidebar:: ToC

   .. contents::


.. _tutorial-string-sets:

String Sets
-----------

Learning Objective
  You will learn the advantages of StringSets and how to work with them.

Difficulty
  Basic

Duration
  15 min

Prerequisites
  :ref:`tutorial-sequences`

A set of sequences can either be stored in a sequence of sequences, for example in a ``String<String<char> >``, or in a :dox:`StringSet`.
This tutorial will introduce you to the SeqAn class :dox:`StringSet`, its background and how to use it.

Background
~~~~~~~~~~

One advantage of using :dox:`StringSet` is that it supports the function :dox:`StringSet#concat` that returns a *concatenator* of all sequences in the string set.
A *concatenator* is an object that represents the concatenation of a set of strings.
This way, it is possible to build up index data structures for multiple sequences by using the same construction methods as for single sequences.

There are two kinds of :dox:`StringSet` specializations in SeqAn: :dox:`OwnerStringSet Owner StringSet`, the default specialisation, and :dox:`DependentStringSet Dependent StringSet`; see the list below for details.
:dox:`OwnerStringSet Owner StringSets` actually store the sequences, whereas :dox:`DependentStringSet Dependent StringSets` just refer to sequences that are stored outside of the string set.

.. code-block:: cpp

   StringSet<DnaString>               ownerSet;
   StringSet<DnaString, Owner<> >     ownerSet2;      // same as above
   StringSet<DnaString, Dependent<> > dependentSet;

The specialization :dox:`ConcatDirectStringSet ConcatDirecet StringSet` already stores the sequences in a concatenation.
The concatenators for all other specializations of :dox:`StringSet` are **virtual** sequences, that means their interface **simulates** a concatenation of the sequences, but they do not literally concatenate the sequences into a single sequence.
Hence, the sequences do not need to be copied when a concatenator is created.

One string can be an element of several :dox:`DependentStringSet Dependent StringSets`.
Typical tasks are, e.g., to find a specific string in a string set, or to test whether the strings in two string sets are the same.
Therefore a mechanism to identify the strings in the string set is needed, and, for performance reasons, this identification should not involve string comparisons.
SeqAn solves this problem by introducing *ids*, which are by default ``unsigned int`` values.

The following list lists the different :dox:`StringSet` specializations:

Specialization ``Owner<ConcatDirect>``
  The sequences are stored as parts of a long string.
  Since the sequences are already concatenated, :dox:`StringSet#concat` just needs to return this string.
  The string set also stores lengths and starting positions of the strings.
  Inserting new strings into the set or removing strings from the set is more expensive than for the default :dox:`OwnerStringSet` specialization, since this involves moving all subsequent sequences in memory.

Specialization ``Dependent<Tight>``
  This specialization stores sequence pointers consecutively in an array.
  Another array stores an id value for each sequence.
  That means that accessing given an id needs a search through the id array.

Specialization ``Dependent<Generous>``
  The sequence pointers are stored in an array at the position of their ids.
  If a specific id is not present, the array stores a zero at this position.
  The advantage of this specialization is that accessing the sequence given its id is very fast.
  On the other hand, accessing a sequence given its position ``i`` can be expensive, since this means we have to find the *i*-th non-zero value in the array of sequence pointers.
  The space requirements of a string set object depends on the largest id rather than the number of sequences stored in the set.
  This could be inefficient for string sets that store a small subset out of a large number of sequences.

Building String Sets
~~~~~~~~~~~~~~~~~~~~

Use the function :dox:`StringConcept#appendValue` to append strings to string sets.

.. code-block:: cpp

   StringSet<DnaString> stringSet;
   DnaString str0 = "TATA";
   DnaString str1 = "CGCG";
   appendValue(stringSet, str0);
   appendValue(stringSet, str1);

Functionality
~~~~~~~~~~~~~

This section will give you a short overview of the functionality of the class :dox:`StringSet`.

There are two ways for accessing the sequences in a string set: (1) the function :dox:`RandomAccessContainerConcept#value` returns a reference to the sequence at a specific *position* within the sequence of sequences, and (2) :dox:`StringSet#valueById` accesses a sequence given its *id*.
We can retrieve the *id* of a sequence in a :dox:`StringSet` with the function :dox:`StringSet#positionToId`.

.. code-block:: cpp

   // (1) Access by position
   std::cout << "Owner: " << '\n';
   std::cout << "Position 0: " << value(stringSet, 0) << '\n';

   // Get the corresponding ids
   unsigned id0 = positionToId(stringSet, 0);
   unsigned id1 = positionToId(stringSet, 1);

   // (2) Access by id
   std::cout << "Id 0:  " << valueById(stringSet, id0) << '\n';

.. code-block:: console

   Owner:
   Position 0: TATA
   Id       0: TATA

In the case of :dox:`OwnerStringSet Owner StringSets`, id and position of a string are always the same, but for :dox:`DependentStringSet Dependent StringSets`, the ids can differ from the positions.
For example, if a :dox:`DependentStringSet Dependent StringSet` is used to represent subsets of strings that are stored in :dox:`OwnerStringSet Owner StringSets`, one can use the position of the string within the :dox:`OwnerStringSet Owner StringSet` as id of the strings.
With the function :dox:`StringSet#assignValueById`, we can add the string with a given id from the source string set to the target string set.

.. code-block:: cpp

   // Let's create a string set of type dependent to represent strings,
   // which are stored in the StringSet of type Owner
   StringSet<DnaString, Dependent<Tight> > depSet;
   // We assign the first two strings of the owner string set to the dependent StringSet,
   // but in a reverse order
   assignValueById(depSet, stringSet, id1);
   assignValueById(depSet, stringSet, id0);

   std::cout << "Dependent: " << '\n';
   // (1) Access by position
   std::cout << "Pos 0: " << value(depSet, 0) << '\n';
   // (2) Access by id
   std::cout << "Id 0:  " << valueById(depSet, id0) << '\n';

.. code-block:: console

   Dependent:
   Position 0: CGCG
   Id       0: TATA

With the function :dox:`StringSet#positionToId` we can show that, in this case, the position and the id of a string are different.


.. code-block:: cpp

   std::cout << "Position 0: Id " << positionToId(depSet, 0) << '\n';
   std::cout << "Position 1: Id " << positionToId(depSet, 1) << '\n';

.. code-block:: console

   Position 0: Id 1
   Position 1: Id 0

Iterating over String Sets
~~~~~~~~~~~~~~~~~~~~~~~~~~

As well as for other containers, SeqAn has implemented iterators for :dox:`StringSet StringSets`.
The generall usage of iterators is described in the tutorial :ref:`tutorial-iterators`.
The following example illustrates, how to iterate over the :dox:`StringSet`.

.. code-block:: cpp

   typedef Iterator<StringSet<DnaString> >::Type TStringSetIterator;
   for (TStringSetIterator it = begin(stringSet); it != end(stringSet); ++it)
   {
       std::cout << *it << '\n';
   }

.. code-block:: console

   TATA
   CGCG

If we want to iterate over the contained :dox:`String Strings` as well, as if the :dox:`StringSet` would be one sequence, we can use the function :dox:`StringSet#concat` to get the concatenation of all sequences.
Therefore we first use the metafunction :dox:`StringSet#Concatenator` to receive the type of the concatenation.
Then, we can simply build an iterator for this type and iterate over the concatenation of all strings.

.. code-block:: cpp

    typedef Concatenator<StringSet<DnaString> >::Type TConcat;
    TConcat concatSet = concat(stringSet);

    Iterator<TConcat>::Type it = begin(concatSet);
    Iterator<TConcat>::Type itEnd = end(concatSet);
    for (; it != itEnd; goNext(it))
    {
        std::cout << getValue(it) << " ";
    }
    std::cout << '\n';

.. code-block:: console

   T A T A C G C G

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Build a string set with default specialization and which contains the strings ``"AAA"``, ``"CCC"``, ``"GGG"`` and ``"TTT"``.
     After that print the length of the string set and use a simple for-loop to print all elements of the strings set.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/string_sets/assignment_1_solution.cpp

Assignment 2
^^^^^^^^^^^^

.. container:: assignment

    Type
      Application

    Objective
      In this task you will test, whether a :dox:`DependentStringSet Dependent StringSet` contains a string without comparing the actual sequences.
      Use the given code frame below and adjust it in the following way:

      #. Build a :dox:`OwnerStringSet Owner StringSet` to store the given strings.
      #. Get the corresponding ids for each position and store them.
      #. Build a :dox:`DependentStringSet` and assign the strings of the owner string set from position 0,1 and 3 by their id to it.
      #. Write a function ``isElement`` which takes a ``StringSet<Dependent<> >`` and a ``Id`` as arguments and checks whether a string set contains a string with a given id.
      #. Check if the string set contains the string of position ``3`` and ``2`` and print the result.

      .. code-block:: cpp

         #include <iostream>
         #include <seqan/sequence.h>
         #include <seqan/file.h>

         using namespace seqan;


         int main()
         {
             // Build strings
             DnaString str0 = "TATA";
             DnaString str1 = "CGCG";
             DnaString str2 = "TTAAGGCC";
             DnaString str3 = "ATGC";
             DnaString str4 = "AGTGTCA";

             // Your code

             return 0;
         }

   Hints
     You can use the SeqAn functions :dox:`StringSet#positionToId` and :dox:`StringSet#assignValueById`.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/string_sets/assignment_2_solution.cpp

Workshop Assignment 4
^^^^^^^^^^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     In this assignment, we pick up the example from the workshop assignments from the sequences and iterators tutorials.
     Take the last solution and change the code to build and use StringSets.

     #. Build a StringSet of readList. Reuse the Rooted iterator above.
     #. Iterate over the StringSet and print out the values.

     .. includefrags:: demos/tutorial/string_sets/assignment_3b_workshop_solution.cpp

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/string_sets/assignment_4_workshop_solution.cpp
