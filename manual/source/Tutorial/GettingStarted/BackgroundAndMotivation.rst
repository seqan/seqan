.. sidebar:: ToC

    .. contents::

.. _tutorial-getting-started-background-and-motivation:

Background and Motivation
=========================

Learning Objective
  You will learn about the design goals and fundamental ideas used in the SeqAn library.
  Also, you will see how the SeqAn library can be generic while still retaining high performance.

Difficulty
  Very basic

Duration
  Take the time you need!

Prerequisites
  Basic C or C++ knowledge

Hi, we are glad you made it here.
You being here, and reading these lines means you are eager to learn more about SeqAn and this is the right place to start.
In this tutorial, we will give you an overview about the design goals, design decisions of the SeqAn library, and explain the motivation for these decisions.
The next chapter :ref:`First Steps <tutorial-getting-started-first-steps-in-seqan>` will flesh out the most important points of this chapter with code examples of everyday SeqAn use.

Library Design Aims
-------------------

The following lists some library design aims of the SeqAn library.
Note that they are contradicting.
The focus is on efficiency but small trade-offs are allowed to improve consistency and ease of use.

#. **Efficiency**.
   The focus of SeqAn is to provide a library of efficient and reusable algorithmic components for biological sequence analysis.
   Algorithms should have good practical implementations with low overhead, even at the cost of being harder to use.
#. **Consistency**.
   Be consistent wherever possible, even at slight costs of efficiency.
#. **Ease of use**.
   The library should be easily usable wherever possible, even at slight costs of efficiency.
#. **Reusability and Generosity**.
   The algorithms in SeqAn should be reusable and generic, even at small costs of efficiency.

Modern C++
----------

C++ is sometimes described as a language that most people know only 20% of but everyone knows a different 20%.
This section gives an overview over some C++ idioms we use.
This might be no news if you are a seasoned C++ programmer who is apt at using the STL and Boost libraries.
However, programmers coming from C and Java might find them interesting (We still encourage to read the `C++ FAQ <https://isocpp.org/faq>`_ if you are new to C++).

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
  Object oriented programming is another key programming paradigm made available with C++ (Compared to C).
  This means, that instead of using raw pointers to allocated chunks of memory, memory management should be done using containers.
  The STL provides containers such as `std::vector <http://www.cplusplus.com/reference/stl/vector/>`_ and SeqAn offers :dox:`String`.

Memory Management in SeqAn
--------------------------

C++ allows to allocate complex objects on the stack (in contrast to Java where objects are always constructed on the heap).
The objects are constructed when the code execution enters the scope/block they are defined in and freed when the block is left.
Allocation of resources (e.g. memory) happens on construction and deallocation happens when the current block is left.
This is best explained in an example.

.. includefrags:: demos/tutorial/background_and_motivation/example.cpp

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
-----------------------------------

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
Using polymorphism (e.g. overwriting of functions in subclasses), one could then write generic and reusable algorithms.
For example, the Java 2 platform defines the sort function for all objects implementing a ``Comparable`` interface.
Note that such an implementation would have to rely on `virtual functions <http://en.wikipedia.org/wiki/Virtual_function>`_ of some sort.
However, as we will see in the section OOP vs. Generic Progamming, **this comes at a high performance cost, being in conflict with efficiency**.
For a sequence library, we could implement functions that map values from an alphabet to an ordinal value between ``0`` and ``S - 1`` where ``S`` is the number of elements in the alphabet.

Generic programming offers one way out: C++ templates allow to define template classes, e.g. the STL's ``std::vector<T>`` or SeqAn's :dox:`String`.
Here, instead of creating a string class around an array of ``char`` values (or objects), we can leave the type of the array's elements open.
We can then introduce different types, e.g. ``Dna`` or ``Dna5`` for 4- and 5-character DNA alphabets.

Algorithms can be implemented using template functions and the template types are fixed at compile time.
Thus, the compiler does not have to use virtual function tables and other "crutches", less indirection is involved, and more code can be inlined and aggressively optimized.
When written appropriately, such algorithms can also work on different string implementations! Also, when defining our own alphabet types, we can directly influence how their abstractions (and APIs) work.

Thus, C++ allows us to implement (1) a generic and reusable library with (2) high level abstractions (and thus ease of use) that still allows the compiler to employ aggressive optimization and thus achieves (3) efficiency.
With the words of the C++ inventor `Bjarne Stroustrup <http://www.artima.com/intv/abstreffi.html>`_:

   A high level of abstraction is good, not just in C++, but in general.
   We want to deal with problems at the level we are thinking about those problems.
   When we do that, we have no gap between the way we understand problems and the way we implement their solutions.
   We can understand the next guy's code. We don't have to be the compiler.

OOP vs. Generic Programming
---------------------------

In SeqAn, we use a technique called :ref:`template subclassing <template-subclassing>` which is based on generic programming.
This technique provides `polymorphism <http://en.wikipedia.org/wiki/Polymorphism_in_object-oriented_programming>`_ into C++ programs at **compile time** using templates.
Such static polymorphism is different from **runtime polymorphism** which is supported in C++ using subclassing and virtual functions.
It comes at the cost of some additional typing but has the advantage that the compiler can inline all function calls and thus achieve better performance.
An example will be given in :ref:`the section "From OOP to SeqAn" in the First Example Tutorial <oop-to-seqan>`.

.. todo::
    We need a little code example here.

The important point is that in contrast to runtime polymorphism such static polymorphism allows the compiler to inline functions, which has huge effect on the overall performance of the program.
Which as you recall correctly from above, is the main objective of the SeqAn library :)!

Global Function Interface
-------------------------

As we already stated, using template subclassing to achieve OOP like behavior in a more efficient way comes with a certain drawback.
Subclassed objects are seen by the compiler as singular instances of a specific type.
That means a subclassed object does not inherit the member or member functions of the alleged base class.
In order to reduce the overhead of reimplementing the same member functions for every subclassed object, we use global interface functions.

You might already have seen global function interfaces while working with the STL.
With the new C++11 standard the STL now provides some global interface functions, e.g., the `begin <http://en.cppreference.com/w/cpp/iterator/begin>`_ or `end <http://en.cppreference.com/w/cpp/iterator/end>`_ interface.

The rationale behind this is the following observation.
Global interface functions allow us to implement a general functionality that is used for all subclassed objects of this template class (assuming the accessed member variables exists in all subclassed objects as in the base template class, otherwise the compiler will complain).
If the behavior for any subclassed object changes, the corresponding global function will be reimplemented for this special type covering the desired functionality.
Due to template deduction the compiler already chooses the correct function and inlines the kernel if possible, which very likely improves the performance of the program.
By this design, we can avoid code duplication, and by that increasing maintainability and reducing subtle errors due to less copy-and-paste code.

So, while most C++ developers, who are familiar with the STL and have a strong background in OO programming, are used to the typical dot notation, in SeqAn you have to get used to global function interfaces instead.
But, cheer up! You will adapt to this very quickly. Promised!

Meta-Programming
----------------

Generic algorithms usually have to know certain types that correspond to their arguments.
An algorithm on containers may need to know which type of values are stored in the string, or what kind of iterator we need to access it.
The usual way in the STL is to define the value type of a class like ``vector`` as a *member typedef* of this class, so it can be retrieved by ``vector::value_type``.

Unfortunately member typedef declarations have the same disadvantages as any members: Since they are specified by the class definition, they cannot be changed or added to the class without changing the code of the class, and it is not possible in C++ to define members for built-in types.
What we need therefore is a mechanism that returns an output type (e.g. the value type) given an input type (e.g. the string) and doing so does not rely on members of the input type, but instead uses some kind of global interface.

Such tasks can be performed by **metafunctions**, also known as **type traits**.
A metafunction is a construct to map some types or constants to other entities like types, constants, functions, or objects at compile time.
The name metafunction comes from fact that they can be regarded as part of a meta-programming language that is evaluated during compilation.

In SeqAn we use class templates to implement metafunctions in C++.
Generic algorithms usually have to know certain types that correspond to their arguments: An algorithm on strings may need to know which type of characters are stored in the string, or what kind of iterator can be used to browse it.
SeqAn uses Metafunctions (also known as "traits") for that purpose.

Looking at an Example
^^^^^^^^^^^^^^^^^^^^^

Assuming that we define a string of amino acids:

.. includefrags:: demos/tutorial/metafunctions/base.cpp
    :fragment: amino

Now lets define a function that exchanges the first two values in a string:

.. includefrags:: demos/tutorial/metafunctions/base.cpp
    :fragment: func_exchange1

Since this function only works for instances of :dox:`String String<`:dox:`AminoAcid AminoAcid>`, we could try to make it more general by making a template out of it.

.. includefrags:: demos/tutorial/metafunctions/base.cpp
    :fragment: func_exchange2

Now the function works for all sequence types ``T`` that store ``AminoAcid`` objects, but it will fail for other value types as soon as the variable temp cannot store ``str[0]`` anymore.
To overcome this problem, we must redefine ``temp`` in a way that it can store a value of the correct type.
The question is: "Given a arbitrary type ``T``, what is the value type of ``T``?"

The metafunction :dox:`ContainerConcept#Value` answers this question: "The value type of ``T`` is given by ``Value<T>::Type``."

Hence, the final version of our function ``exchangeFirstValues`` reads as follows:

.. includefrags:: demos/tutorial/metafunctions/base.cpp
    :fragment: func_exchange3

We can view ``Value`` as a kind of "function" that takes ``T`` as an argument (in angle brackets) and returns the required value type of ``T``.
In fact, ``Value`` is not implemented as a C++ function, but as a class template.
This class template is specialized for each sequence type ``T`` in a way that the ``typedef Type`` provides the value type of ``T``.
Unfortunately, the current C++ language standard does not allow to write simply "``Value<T> temp``;", so we must select the return value by appending "``::Type``".
The leading "``typename``" becomes necessary since ``Value<T>::Type`` is a type that depends on a template parameter of the surrounding function template.

And now?
--------

Wow, this was quite some information to digest, wasn't it?
We suggest you take a break!
Get some fresh air!
Grab something to drink or to eat!
Let the information settle down.

Do you think you've got everything?
Well, if not don't worry!
Follow the :ref:`First Steps <tutorial-getting-started-first-steps-in-seqan>` tutorial which will cover the topics discussed above.
This gives you the chance to apply the recently discussed paradigms to an actual (uhm, simplistic) use case.
But it will help you to better understand the way data structures and algorithms are implemented in SeqAn.

We recommend you to also read the :ref:`Argument Parser Tutorial <tutorial-getting-started-parsing-command-line-arguments>`.
This tutorial will teach you how to easily add command line arguments for your program and how to generate a help page for the options.
Or you go back to the :ref:`main page <manual-main-tutorials>` and stroll through the other tutorials.
You are now ready to dive deeper into SeqAn.
Enjoy!
