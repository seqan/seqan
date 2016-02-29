.. sidebar:: ToC

    .. contents::

.. _tutorial-datastructures-modifiers:

Modifiers
=========

Learning Objective
  In this tutorial you will learn how to modify the elements of a container without copying them using SeqAn modifiers.
  You will learn about the different specializations and how to work with them.

Difficulty
  Basic

Duration
  20 min

Prerequisites
  :ref:`tutorial-getting-started-first-steps-in-seqan`, :ref:`tutorial-datastructures-sequences`

Overview
--------

Modifiers give a different view to other classes.
They can be used to change the elements of a container without touching them.
For example, someone gave you an algorithm that works on two arbitrary :dox:`String` objects, but you want to use it for the special pair of a string and its reverse (left-to-right mirror).
The classical approach would be to make a copy of the one string, where all elements are mirrored from left to right and call the algorithm with both strings.
With modifiers, e.g. a :dox:`ModifiedString`, you can create the reverse in :math:`\mathcal{O}(1)` extra memory without copying the original string.
This can be handy if the original sequence is large.

Modifiers implement a certain concept (e.g. :dox:`ContainerConcept`, :dox:`RandomAccessIteratorConcept Iterator`, ...) or class interface (:dox:`String`, ...) and thus can be used as such.
The mirror modifier is already part of SeqAn and implements the class interface of :dox:`String` and can be used in every algorithm that works on strings.

The Modified String
-------------------

The :dox:`ModifiedString ModifiedString` is a modifier that implements the :dox:`String` interface and thus can be used like a :dox:`String`.
It has two template parameters.
The first one specifies a sequence type (e.g. :dox:`String`, :dox:`Segment`, ...) and the second one specifies the modifiers behavior.
That can be :dox:`ModReverseString` for mirroring a string left to right or :dox:`ModViewModifiedString` for applying a function to every single character (like 'C'->'G', 'A'->'T', ...).

ModReverse
^^^^^^^^^^

We begin with the specialization :dox:`ModReverseString` from the example above.
Now we have a given string:

.. includefrags:: demos/tutorial/modifiers/modreverse.cpp
   :fragment: main

and want to get the reverse.
So we need a :dox:`ModifiedString` specialized with ``String<char>`` and :dox:`ModReverseString`.
We create the modifier and link it with ``myString``:

.. includefrags:: demos/tutorial/modifiers/modreverse.cpp
   :fragment: modifier

The result is:

.. includefrags:: demos/tutorial/modifiers/modreverse.cpp
   :fragment: output1

.. includefrags:: demos/tutorial/modifiers/modreverse.cpp.stdout
   :fragment: output1

To verify that we didn't copy ``myString``, we replace an infix of the original string and see that, as a side effect, the modified string has also changed:

.. includefrags:: demos/tutorial/modifiers/modreverse.cpp
   :fragment: output2

.. includefrags:: demos/tutorial/modifiers/modreverse.cpp.stdout
   :fragment: output2

ModView
^^^^^^^

Another specialization of the :dox:`ModifiedString` is the :dox:`ModViewModifiedString` modifier.
Assume we need all characters of ``myString`` to be in upper case without copying ``myString``.
In SeqAn you first create a functor (a STL unary function) which converts a character to its upper-case character.

.. includefrags:: demos/tutorial/modifiers/modview.cpp
   :fragment: functor

and then create a :dox:`ModifiedString` specialized with ``ModView<MyFunctor>``:

.. includefrags:: demos/tutorial/modifiers/modview.cpp
   :fragment: mod_str

The result is:

.. includefrags:: demos/tutorial/modifiers/modview.cpp
   :fragment: output

.. includefrags:: demos/tutorial/modifiers/modview.cpp.stdout

The upper-case functor and some other predefined functors are part of SeqAn (in ``seqan/modifier/modifier_functors.h``) already.
The following functors can be used as an argument of :dox:`ModViewModifiedString`:

``FunctorUpcase<TValue>``
  Converts each character of type ``TValue`` to its upper-case character

``FunctorLowcase<TValue>``
  Converts each character to type ``TValue`` to its lower-case character

``FunctorComplement<Dna>``
  Converts each nucleotide to its complementary nucleotide

``FunctorComplement<Dna5>``
  The same for the :dox:`Dna5` alphabet

``FunctorConvert<TInValue,TOutValue>``
  Converts the type of each character from ``TInValue`` to ``TOutValue``

So instead of defining your own functor we could have used a predefined one:

.. includefrags:: demos/tutorial/modifiers/modview.cpp
   :fragment: predefined

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     In this assignment you will create a modifier using your own functor.
     Assume you have given two Dna sequences as strings as given in the code example below.
     Let's assume you know that in one of your Dna sequences a few 'C' nucleotides are converted into 'T' nucleotides, but you still want to compare the sequence.
     Extend the code example as follows:

     #. Write a functor which converts all 'C' nucleotides to 'T' nucleotides.
     #. Define a :dox:`ModifiedString` with the specialization :dox:`ModViewModifiedString` using this functor.
     #. Now you can modify both sequences to compare them, treating all 'Cs' as 'Ts'.
        Print the results.

    .. includefrags:: demos/tutorial/modifiers/assignment1.cpp

   Solution
      .. container:: foldable

         .. includefrags:: demos/tutorial/modifiers/assignment1_solution.cpp

         .. includefrags:: demos/tutorial/modifiers/assignment1_solution.cpp.stdout

^^^^^^^^^

For some commonly used modifiers you can use the following shortcuts:

+-----------------------------------+---------------------------------------------------------------------------------+
| Shortcut                          | Substitution                                                                    |
+===================================+=================================================================================+
| ``ModComplementDna``              | ``ModView<FunctorComplement<Dna> >``                                            |
+-----------------------------------+---------------------------------------------------------------------------------+
| ``ModComplementDna5``             | ``ModView<FunctorComplement<Dna5> >``                                           |
+-----------------------------------+---------------------------------------------------------------------------------+
| ``DnaStringComplement``           | ``ModifiedString<DnaString, ModComplementDna>``                                 |
+-----------------------------------+---------------------------------------------------------------------------------+
| ``Dna5StringComplement``          | ``ModifiedString<Dna5String, ModComplementDna5>``                               |
+-----------------------------------+---------------------------------------------------------------------------------+
| ``DnaStringReverse``              | ``ModifiedString<DnaString, ModReverse>``                                       |
+-----------------------------------+---------------------------------------------------------------------------------+
| ``Dna5StringReverse``             | ``ModifiedString<Dna5String, ModReverse>``                                      |
+-----------------------------------+---------------------------------------------------------------------------------+
| ``DnaStringReverseComplement``    | ``ModifiedString<ModifiedString<DnaString, ModComplementDna>, ModReverse>``     |
+-----------------------------------+---------------------------------------------------------------------------------+
| ``Dna5StringReverseComplement``   | ``ModifiedString<ModifiedString<Dna5String, ModComplementDna5>, ModReverse>``   |
+-----------------------------------+---------------------------------------------------------------------------------+

The Modified Iterator
---------------------

We have seen how a :dox:`ModifiedString` can be used to modify strings without touching or copying original data.
The same can be done with iterators.
The :dox:`ModifiedIterator` implements the :dox:`RandomAccessIteratorConcept Iterator` concept and thus can be used in every algorithm or data structure that expects an iterator.
In fact, we have already used the :dox:`ModifiedIterator` unknowingly in the examples above, as in our cases the :dox:`ModifiedString` returns a corresponding :dox:`ModifiedIterator` via the :dox:`ContainerConcept#Iterator` meta-function.
The main work is done in the :dox:`ModifiedIterator`, whereas the :dox:`ModifiedString` only overloads the :dox:`ContainerConcept#begin` and :dox:`ContainerConcept#end`.
Normally, you are going to use the :dox:`ModifiedString` and maybe the result of its :dox:`ContainerConcept#Iterator` meta-function instead of a :dox:`ModifiedIterator` directly.

Nested Modifiers
----------------

As modifiers implement a certain concept and depend on classes of this concept, two modifiers can be chained to create a new modifier.
We have seen how the :dox:`ModifiedString` specialized with :dox:`ModReverseString` and :dox:`ModViewModifiedString` can be used.
Now we want to combine them to create a modifier for the reverse complement of a :dox:`DnaString` We begin with the original string:

.. includefrags:: demos/tutorial/modifiers/nested.cpp
   :fragment: string

Then we define the modifier that complements a :dox:`DnaString`:

.. includefrags:: demos/tutorial/modifiers/nested.cpp
   :fragment: complement

This modifier now should be reversed from left to right:

.. includefrags:: demos/tutorial/modifiers/nested.cpp
   :fragment: reverse

The original string can be given to the constructor.

.. includefrags:: demos/tutorial/modifiers/nested.cpp
   :fragment: constructor

The result is:

.. includefrags:: demos/tutorial/modifiers/nested.cpp
   :fragment: output

.. includefrags:: demos/tutorial/modifiers/nested.cpp.stdout
    :fragment: output


Using a predefined shortcut, the whole example could be reduced to:

.. includefrags:: demos/tutorial/modifiers/nested.cpp
   :fragment: alternative
