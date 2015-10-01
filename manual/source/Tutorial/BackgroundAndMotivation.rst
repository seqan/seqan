.. sidebar:: ToC

   .. contents::


.. _tutorial-background-and-motivation:

Background and Motivation
-------------------------

Learning Objective
  You will learn about the design goals and fundamental ideas used in the SeqAn library.
  Also, you will see how the SeqAn library obtains genericity while still retaining high performance.

Difficulty
  Very basic

Duration
  1h

Prerequisites
  Basic C or C++ knowledge

This tutorial is meant to be the first chapter you read in the SeqAn Tutorial.
It will give you an overview about the design goals, design decisions, and explain the motivation for these decisions.
The next chapter :ref:`First Examples <tutorial-first-steps-in-seqan>` will flesh out the most important points of this chapter with code examples of everyday SeqAn use.

Library Design Aims
~~~~~~~~~~~~~~~~~~~

The following lists some library design aims of the SeqAn library.
Note that they are contradicting.
The focus is on efficiency but small trade-offs are allowed to improve consistency and ease of use.

#. **Efficiency**.
   The focus of SeqAn is to provide a library of efficient and reuseable algorithmic components for biological sequence analysis.
   Algorithms should have good practical implementations with low overhead, even at the cost of being harder to use.
#. **Consistency**.
   Be consistent wherever possible, even at slight costs of efficiency.
#. **Ease of use**.
   The library should be easily usable wherever possible, even at slight costs of efficiency.
#. **Reuseability and Generosity**.
   The algorithms in SeqAn should be reuseable and generic, even at small costs of efficiency.

Modern C++ (C++98)
~~~~~~~~~~~~~~~~~~

C++ is sometimes described as a language that most people know only 20% of but everyone knows a different 20%.
This section gives an overview over some C++98 idioms we use.
This might be no news if you are a seasoned C++ programmer who is apt at using the STL and Boost libraries.
However, programmers coming from C and Java might find them interesting.

References
  References are alternatives to pointers in C++ to construct value aliases.
  Also see `Wikipedia on C++ references <http://en.wikipedia.org/wiki/Reference_(C%2B%2B)>`_.

Templates
  C++ allows you to perform `generic programming <http://en.wikipedia.org/wiki/Generic_programming>`_ using templates.
  While similar to generics in Java (C++ templates are more than a decade older), C++ templates are designed to write zero-overhead abstractions that can be written to be as efficient as hand-written code while retaining a high level of abstraction.
  See `cplusplus.com on C++ templates <http://www.cplusplus.com/doc/tutorial/templates/>`_.
  Note that there is no way to restrict the type that can be used in templates, there is no mechanism such as Java's ``? extends T`` in C++.
  Using an incompatible type leads to compiler errors because some operator or function could not be found.

Memory Management / No Pointers
  Instead of using raw pointers and memory, memory management should be done using containers.
  The STL provides containers such as `std::vector <http://www.cplusplus.com/reference/stl/vector/>`_ and SeqAn offers :dox:`String`.

Memory Management in SeqAn
~~~~~~~~~~~~~~~~~~~~~~~~~~

C++ allows to allocate complex objects on the stack (in contrast to Java where objects are always constructed on the heap).
The objects are constructed when the code execution enters the scope/block they are defined in and freed when the block is left.
Allocation of resources (e.g. memory) happens on construction and deallocation happens when the current block is left.
This is best explained in an example.

.. code-block:: cpp

   #include <seqan/sequence.h>

   int main(int argc, char const ** argv)
   {
       seqan::String<char> programName = argv[0];
       if (argc > 1)
       {
           seqan::String<char> firstArg = argv[1];
           if (argc > 2)
               return 1;
       }
       return 0;
   };

``seqan::String<char>`` is a class (actually an instantiation of the class template :dox:`String`) that allows to store strings of ``char`` values, similar to ``std::vector<char>`` or ``std::string``.

When the variable ``programName`` is allocated, the constructor of the ``String<char>`` class is called.
It allocates sufficient memory to store the value of ``argv[0]`` and then copies over the values from this string.
The variable exists until the current block is left.
Since it is defined in the ``main()`` function, this can only happen in the last line of ``main()`` at the ``return 0``.
When the variable goes out of scope, its value is deconstructed and all allocated memory is freed.

If an argument was given to the program, the block in the ``if`` clause is entered.
When this happens, the variable ``firstArg`` is constructed, memory is allocated and the value of ``argv[1]`` is copied into the buffer.
When the block is left, the variable is deconstructed and all memory is deallocated.

Note that all memory is released when the ``main()`` function is left, regardless whether it is left in the ``return 0`` or the ``return 1``.
Corresponding code in C would be (arguably) more messy, either requiring ``goto`` or multiple ``free()`` calls, one before either ``return``.

Motivation for Template Programming
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section, we will give a short rationale why C++ with heavy use of template programming was used for SeqAn.

Any sequence analysis will have sequence data structures and algorithms on sequences at its heart.
Even when only considering DNA and amino acid alphabets, there are various variants for alphabets that one has to consider.
Otherwise, important applications in bioinformatics cannot be covered:

* 4-character DNA,
* 5-character DNA with ``N``,
* 15-character IUPAC, and
* 27-character amino acids.

A simple implementation could simply store such strings as ASCII characters.
However, there are some implementation tricks that can lead to great reduction of memory usage (e.g. encoding eight 4-character DNA characters in one byte) or running time (fast lookup tables for characters or q-grams) for small alphabets.
Thus, simply using a ``std::string`` would come at high costs to efficiency.

Given that in the last 10-15 years, Java and C# have gained popularity, one could think about an object oriented solution: strings could simply be arrays of ``Character`` objects.
Using polymorphism (e.g. overwriting of functions in subclasses), one could then write generic and reuseable algorithms.
For example, the Java 2 platform defines the sort function for all objects implementing a ``Comparable`` interface.
Note that such an implementation would have to rely on `virtual functions <http://en.wikipedia.org/wiki/Virtual_function>`_ of some sort.
However, as we will see in the section OOP vs. Template Subclassing, **this comes at a high performance cost, being in conflict with efficiency**.
For a sequence library, we could implement functions that map values from an alphabet to an ordinal value between ``0`` and ``S - 1`` where ``S`` is the number of elements in the alphabet.

Generic programming offers one way out: C++ templates allow to define template classes, e.g. the STL's ``std::vector<T>`` or SeqAn's :dox:`String`.
Here, instead of creating a string class around an array of ``char`` values (or objects), we can leave the type of the array's elements open.
We can then introduce different types, e.g. ``Dna`` or ``Dna5`` for 4- and 5-character DNA alphabets.

Algorithms can be implemented using templated functions and the template types are fixed at compile time.
Thus, the compiler does not have to use virtual function tables and other "crutches", less indirection is involved, and more code can be inlined and aggressively optimized.
When written appropriately, such algorithms can also work on different string implementations! Also, when defining our own alphabet types, we can directly influence how their abstractions (and APIs) work.

Thus, C++ allows us to implement (1) a generic and reuseable library with (2) high level abstractions (and thus ease of use) that still allows the compiler to employ aggressive optimization and thus achieves (3) efficiency.
With the words of the C++ inventor `Bjarne Stroustrup <http://www.artima.com/intv/abstreffi.html>`_:

   A high level of abstraction is good, not just in C++, but in general.
   We want to deal with problems at the level we are thinking about those problems.
   When we do that, we have no gap between the way we understand problems and the way we implement their solutions.
   We can understand the next guy's code. We don't have to be the compiler.

OOP vs. Generic Programming
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In SeqAn, we use a technique called `template subclassing <tutorial-template-subclassing>`_ which is based on generic programming.
This technique provides `polymorphism <http://en.wikipedia.org/wiki/Polymorphism_in_object-oriented_programming>`_ into C++ programs at **compile time** using templates.
Such static polymorphism is different from **runtime polymorphism** which is supported in C++ using subclassing and virtual functions.
It comes at the cost of some additional typing but has the advantage that the compiler can inline all function calls and thus achieve better performance.
An example will be given in `the section "From OOP to SeqAn" in the First Steps Tutorial <tutorial-first-steps-in-seqan>`_.

The important point is that in contrast to runtime polymorphism such static polymorphism allows the compiler to inline functions.
