.. sidebar:: ToC

   .. contents::


.. _tutorial-basics:

Basics
------

Alphabets
~~~~~~~~~

In SeqAn, alphabets are value types that can take a limited number of values and which hence can be mapped to a range of natural numbers.
We can retrieve the number of different values of an alphabet, the alphabet size, by the metafunction :dox:`FiniteOrderedAlphabetConcept#ValueSize`.
Another useful metafunction called :dox:`AlphabetConcept#BitsPerValue` can be used to determine the number of bits needed to store a value of a given alphabet.
The order of a character in the alphabet (i.e. its corresponding natural number) can be retrieved by calling the function :dox:`FiniteOrderedAlphabetConcept#ordValue`.
In SeqAn, several standard alphabets are already predefined, for example :dox:`Dna` :dox:`Dna5`, :dox:`Rna`, :dox:`Rna5`, :dox:`Iupac`, :dox:`AminoAcid`, ....

Let's start with a simple task. We want to write a program that outputs all letters of the predefined :dox:`AminoAcid` alphabet.
First we include the corresponding header files and specify that we are using the namespace ``seqan``.

.. includefrags:: core/demos/tutorial/basics/show_alphabets.cpp
   :fragment: includes

Next, we will define a template function ``template<typename TAlphabet> void showAllLettersOfMyAlphabet(TAlphabet const&)`` which will iterate over all the characters of the alphabet and outputs them.
For this, we need to determine the alphabet size using the metafunction :dox:`FiniteOrderedAlphabetConcept#ValueSize ValueSize<TAlphabet>::VALUE`.
Then we increment a counter from 0 to the alphabet size minus one and output the counter as well as the corresponding letter of the alphabet using a conversion constructor.

.. includefrags:: core/demos/tutorial/basics/show_alphabets.cpp
   :fragment: showAllLettersOfMyAlphabet

In the main program we simply call the above function using a number of alphabets that are predefined in SeqAn.

.. includefrags:: core/demos/tutorial/basics/show_alphabets.cpp
   :fragment: main

This program produces the following output:

.. code-block:: console

     darwin10.0 : ./show_alphabets
    0,A  1,R  2,N  3,D  4,C  5,Q  6,E  7,G  8,H  9,I  10,L  11,K  12,M  13,F  14,P  15,S  16,T  17,W  18,Y  19,V  20,B  21,Z  22,X  23,*
    0,A  1,C  2,G  3,T
    0,A  1,C  2,G  3,T  4,N

Iterators
~~~~~~~~~

An iterator is an object that is used to browse through the values of a container.
The metafunction :dox:`ContainerConcept#Iterator` can be used to determine an appropriate iterator type given a container.
Some containers offer several kinds of iterators, which can be selected by an optional argument of Iterator.
For example, the tag :dox:`ContainerIteratorTags#Standard` can be used to get an iterator type that resembles the C++ standard random access iterator.
The more elaborated :dox:`RootedIteratorConcept Rooted\ Iterator`, i.e., an iterator that knows its container, can be selected by specifying the :dox:`ContainerIteratorTags#Rooted` tag.

Rooted iterators offer some convenience for the user: They offer additional functions like :dox:`RootedIteratorConcept#container` for determining the container on which the iterator works, and they simplify the interface for other functions like :dox:`RootedIteratorConcept#atEnd`.
Moreover, rooted iterators may change the container’s length or capacity, which makes it possible to implement a more intuitive variant of a remove algorithm.

While rooted iterators can usually be converted into standard iterators, it is not always possible to convert standard iterators back into rooted iterators, since standard iterators may lack the information about the container they work on.
Therefore, many functions that return iterators like :dox:`ContainerConcept#begin` or :dox:`ContainerConcept#end` return rooted iterators instead of standard iterators; this way, they can be used to set both rooted and standard iterator variables.
Alternatively it is possible to specify the returned iterator type explicitly by passing the iterator kind as a tag argument.

The following code piece shows examples for creating Iterators for :dox:`ContainerConcept Containers`.
If no iterator kind is specified, the metafunction :dox:`ContainerConcept#Iterator` assumes :dox:`ContainerIteratorTags#Standard` and the function :dox:`ContainerConcept#begin` assumes :dox:`ContainerIteratorTags#Rooted`.
Both ``it1`` and ``it2`` are standard iterators, whereas ``it3`` and ``it4`` are rooted iterators.

.. code-block:: cpp

    String<char> str = "ACME";
    Iterator<String<char> >::Type it1 = begin(str); //a standard iterator
    Iterator<String<char>, Standard>::Type it2 = begin(str);  //same as above
    Iterator<String<char>, Rooted>::Type it3 = begin(str);  //a rooted iterator
    Iterator<String<char>, Rooted>::Type it4 = begin(str, Rooted());  //same as above

.. comment

    An iterator is stable if it stays valid even if its container is expanded, otherwise it is unstable. For example, the standard iterator of :dox:`AllocString` – which is a simple pointer to a value in the string – is unstable, since during the expansion of an Alloc String, all values are moved to new memory addresses.
    A typical implementation of stable iterators for strings stores the position instead of a pointer to the current value.
    The :dox:`Iterator` metafunction called with the [seqan:"Tag.Iterator Spec" Stable] tag returns a type for stable iterators.

    Stable tag does not appear in Doku. Clarify with Andreas.

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Transfer

   Objective

      Write a program which does the following:
        #. Create an amino acid string of the following sequence: "MQDRVKRPMNAFIVWSRDQRRKMALEN".
        #. Iterate through the sequence and replace all ‘R’ with ‘A’.
        #. Create a second string where you count the number of occurrences of each amino acid.
        #. Iterate through the latter string and output the frequency table.

   Hints

      After a few hours browsing through the demos you should be able to solve this.

   Solution

      .. container:: foldable

         In this assignment we practice the use of alphabets, iterators and metafunctions in SeqAn. We start by including the seqan basic header and enter the namespace ``seqan`` to avoid writing it as a prefix (as we do with the namespace ``std`` in this example).
         In the ``main`` function we first define a a type ``TAmincoAcidString`` which is a ``String<AminoAcid>`` (Note the SeqAn naming conventions).
         Then we define a variable ``sourceSeq`` of this type and initialize it with a string constant.

         .. comment

            Add link to naming conventions

         .. includefrags:: core/demos/tutorial/basics/strings.cpp
            :fragment: create-string

         Then we define an iterator type using the SeqAn metafunction :dox:`ContainerConcept#Iterator`.
         Using the correct iterator we iterate over our amino acid string using the SeqAn functions :dox:`ContainerConcept#begin`, :dox:`ContainerConcept#end`, and :dox:`InputIteratorConcept#goNext`.
         In the body of the while loop we use the SeqAn function :dox:`IteratorAssociatedTypesConcept#value` to access the value the iterator is pointing to.
         Note that this function returns a reference which allows us to replace the occurrence of all ``R``'s with ``A``'s.
         So at this point we have solved parts a) and b) of the assignment.

         .. includefrags:: core/demos/tutorial/basics/strings.cpp
            :fragment: iterate-and-replace

         In the next part of the code we want to count, how often a specific letter of the alphabet occurs in the string.
         To obtain the size type of the used alphabet we call the SeqAn metafunction :dox:`ContainerConcept#Size Size` and define a :dox:`String` of that type to hold the counters.
         The :dox:`String` has here basically the same functionality as a STL ``vector``.
         Since alphabets are mapped to a contiguous interval of the natural numbers, we can initialize the counter up to the size of the alphabet which we obtain by a call to the SeqAn metafunction :dox:`ContainerConcept#ValueSize ValueSize`.
         We then iterate over the amino acid string and increment the counter for the corresponding letter of the alphabet.
         In order to know the corresponding natural number of an alphabet letter, we use the SeqAn function :dox:`FiniteOrderedAlphabetConcept#ordValue`.
         Note the use of the :dox:`IteratorAssociatedTypesConcept#value` function.
         In this example one could also use the ``operator[]`` to write ``counter[ordValue(value(it))]++``.

         .. includefrags:: core/demos/tutorial/basics/strings.cpp
            :fragment: count-occurrences

         Finally we iterate through the counter String and output the i-th aminoacid (by calling a constructor with the letter's ordinal value) ad its frequency.

         .. includefrags:: core/demos/tutorial/basics/strings.cpp
            :fragment: frequency-table

         The result looks like this:

         .. code-block:: console

             $darwin10.0 : basics//strings
             M,Q,D,A,V,K,A,P,M,N,A,F,I,V,W,S,A,D,Q,A,A,K,M,A,L,E,N,
             A:7
             R:0
             N:2
             D:2
             C:0
             Q:2
             E:1
             G:0
             H:0
             I:1
             L:1
             K:2
             M:3
             F:1
             P:1
             S:1
             T:0
             W:1
             Y:0
             V:2
             B:0
             Z:0
             X:0
             *:0


Memory Allocation
~~~~~~~~~~~~~~~~~

Controlling memory allocation is one of the big advantages of C++ compared to other programming languages as for example Java.
Depending on the size of objects and the pattern they are allocated during the program execution, certain memory allocation strategies have advantages compared to others.
SeqAn supports a variety of memory allocation strategies.

The two functions :dox:`Allocator#allocate` and :dox:`Allocator#deallocate` are used in SeqAn to allocate and deallocate dynamic memory.
Both functions take an allocator as an argument.
An :dox:`Allocator` is an object that is responsible for allocated memory.
The default implementations of :dox:`Allocator#allocate` and :dox:`Allocator#deallocate` completely ignore the allocator but simply call the basic operators ``new`` and ``delete``.
Although in principle every kind of object can be used as allocator, typically the object that stores the pointer to the allocated memory is used as allocator.
For example, if memory is allocated for an :dox:`AllocString Alloc String`, this string itself acts as allocator.
A memory block should be deallocated using the same allocator object as it was allocated for.
The following allocators are available in SeqAn and support the :dox:`Allocator#clear` function.
This function deallocates at once all memory blocks that were previously
allocated.

| :dox:`SimpleAllocator Simple Allocator`
|    General purpose allocator.
| :dox:`SinglePoolAllocator Single Pool Allocator`
|    Allocator that pools memory blocks of specific size. Blocks of different sizes are not pooled.
| :dox:`MultiPoolAllocator Multi Pool Allocator`
|    Allocator that pools memory blocks. Only blocks up to a certain size are pooled. The user can specify the size limit in a template argument.

The function :dox:`Allocator#allocate` has an optional argument to specify the intended allocator usage for the requested memory.
The user can
thereby specialize :dox:`Allocator#allocate` for different allocator applications.
For example, the tag :dox:`AllocatorUsageTags#TagAllocateTemp` specifies that the memory will only be used temporarily, whereas :dox:`AllocatorUsageTags#TagAllocateStorage` indicates that the memory will be used in the long run for storing values of a container.

SeqAn also offers more complex allocators which support the function :dox:`Allocator#clear`.
The library predefines some allocator specializations for different uses (see above).
Most of these allocators are pool allocators.
A pool allocator implements its own memory management.
It reserves storage for multiple memory blocks at a time and recycles deallocated blocks.
This reduces the number of expensive ``new`` and ``delete`` calls and speeds up the allocation and deallocation.

Assignment 2
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective

      Write a program which compares the runtimes of the :dox:`SimpleAllocator Simple Allocator` and the :dox:`MultiPoolAllocator Multi Pool Allocator` for pool sizes (10,100,1000) for allocating and deallocating memory.

   Hint

      .. container:: foldable

         For timing the allocation you can use :dox:`sysTime`.

   Solution

      .. container:: foldable

         We start in this assignment by including the ``basic.h`` SeqAn header and defining two different allocators, one :dox:`MultiPoolAllocator Multi Pool Allocator` and one :dox:`SimpleAllocator Simple Allocator`.

         .. includefrags:: core/demos/tutorial/basics/allocator.cpp
            :fragment: definitions

         Given these fixed allocators we allocate now various size blocks, namely of size 10, 100, and 1000.
         We repeat the allocation a number of times and then clear the allocated memory.
         For each of the block sizes we output the system time needed to allocate and clear the memory.

         .. includefrags:: core/demos/tutorial/basics/allocator.cpp
            :fragment: time-measurements

         Running this program results in the following output which shows the advantage of the :dox:`MultiPoolAllocator Multi Pool Allocator`:

         .. code-block:: console

            $ darwin10.0 : cd ~/seqan/projects/library/demos/tutorial
            $ darwin10.0 : ./basics/allocator
            Allocating and clearing 100000 times blocks of size 10 with MultiPool Allocator took 0.00200295
            Allocating and clearing 100000 times blocks of size 10 with Standard Allocator took 0.0451179
            Allocating and clearing 100000 times blocks of size 100 with MultiPool Allocator took 0.0599239
            Allocating and clearing 100000 times blocks of size 100 with Standard Allocator took 0.127033
            Allocating and clearing 100000 times blocks of size 1000 with MultiPool Allocator took 0.368732
            Allocating and clearing 100000 times blocks of size 1000 with Standard Allocator took 0.560434
