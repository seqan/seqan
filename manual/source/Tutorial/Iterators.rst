.. sidebar:: ToC

   .. contents::


.. _tutorial-iterators:

Iterators
---------

Learning Objective
  You will learn how to use iterators to traverse containers in SeqAn.
  After this tutorial, you will be ready to continue with the tutorials about iterating on more complex structures, e.g. :ref:`tutorial-index-iterators`.

Difficulty
  Basic

Duration
  30 min

Prerequsites
  :ref:`tutorial-sequences`

Iterators are objects that can be used to browse through the values of containers such as :dox:`String Strings` or :dox:`Segment Segments`.
SeqAn also offers a range of iterators to traverse efficiently more complex data structures, e.g. :dox:`Graph Graphs`, whose specific usage will be explained in the corresponding tutorials.
This tutorial will introduce you into the basic concept of iterators using :dox:`String` iterators as illustration.

Defining Iterators
~~~~~~~~~~~~~~~~~~

This section will show you how to define different kinds of iterators.

The metafunction :dox:`ContainerConcept#Iterator` can be used to determine the appropriate iterator type for a given a container.
Some containers offer several kinds of iterators, which can be selected by an optional argument of Iterator.
For example, the tag :dox:`ContainerIteratorTags#Standard` can be used to get an iterator type that resembles the C++ standard random access iterator.
The more elaborated :dox:`ContainerIteratorTags#Rooted` iterator, i.e., an iterator that knows its container, can be selected by specifying the :dox:`ContainerIteratorTags#Rooted` tag.
The construction of an iterator in SeqAn, e.g. for a :dox:`DnaString Dna String`, could look like the following:

.. code-block:: cpp

   Iterator<DnaString>::Type           it1;  // A standard iterator
   Iterator<DnaString, Standard>::Type it2;  // Same as above
   Iterator<DnaString, Rooted>::Type   it3;  // A rooted iterator

.. tip::

   The default iterator implementation is :dox:`ContainerIteratorTags#Standard`.
   Rooted iterators offer some convenience for the user.
   They offer additional functions like :dox:`RootedIteratorConcept#container` for determining the container on which the iterator works, and they simplify the interface for other functions like :dox:`RootedIteratorConcept#atEnd`.
   Moreover, rooted iterators may change the container’s length or capacity, which makes it possible to implement a more intuitive variant of a remove algorithm.

   While rooted iterators can usually be converted into standard iterators, it is not always possible to convert standard iterators back into rooted iterators, since standard iterators may lack the information about the container they work on.
   Therefore, many functions that return iterators like :dox:`ContainerConcept#begin` or :dox:`ContainerConcept#end` return rooted iterators instead of standard iterators; this way, they can be used to set both rooted and standard iterator variables.
   Alternatively it is possible to specify the returned iterator type explicitly by passing the iterator kind as a tag argument, e.g. ``begin(str, Standard())``.

Traversing Containers
~~~~~~~~~~~~~~~~~~~~~

In this section you will learn how to iterate over a container using the basic functionality of iterators.

An iterator always points to one value of the container.
The function :dox:`RandomAccessContainerConcept#value`, which is equivalent to the ``operator*``, can be used to access this value by reference.
In contrast :dox:`RandomAccessContainerConcept#getValue` return a copy of the value.
Functions like :dox:`InputIteratorConcept#goNext` or :dox:`BidirectionalIteratorConcept#goPrevious`, which are equivalent to ``operator++`` and ``operator--`` respectively, can be used to move the iterator to other values within the container.

The functions :dox:`ContainerConcept#begin` and :dox:`ContainerConcept#end`, applied to a container, return iterators to the begin and to the end of the container.
Note that similar to C++ standard library iterators, the iterator returned by :dox:`ContainerConcept#end` does not point to the last value of the container but to the position behind the last one.
If the container is empty then ``end() == begin()``.

The following code prints out a sequence and demonstrates how to iterate over a string.

.. code-block:: cpp

   DnaString genome = "ACGTACGTACGT";
   typedef Iterator<DnaString>::Type TIterator;
   for (TIterator it = begin(genome); it != end(genome); goNext(it))
   {
       std::cout << value(it);
   }

.. code-block:: console

    ACGTACGTACGT

A Working Example
~~~~~~~~~~~~~~~~~

Let us now clarify the usage of iterators with a working example.
The following program demonstrates the usage of iterators.

.. includefrags:: core/demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: includes

The metafunction :dox:`ContainerConcept#Iterator` returns the iterator type for a given container type.
In this case the default implementation :dox:`ContainerIteratorTags#Standard` is used.

.. includefrags:: core/demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: metafunctions

We can use iterators to iterate over the elements of a container, e.g.  to print the elements.

.. includefrags:: core/demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: iterators

Instead of comparing the two iterators ``it`` and ``itEnd``, we could also use the function :dox:`RootedIteratorConcept#atEnd` to check whether we reached the end of the container.

.. includefrags:: core/demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: standard-iterators

Next we will use :dox:`RootedIteratorConcept Rooted Iterators`.
Since :dox:`RootedIteratorConcept Rooted Iterators` know their container, the functions :dox:`RootedRandomAccessIteratorConcept#goBegin` and :dox:`RootedIteratorConcept#atEnd` do not need to get the container as an argument.
The following example prints for each element of the :dox:`Dna5String Dna5 String` ``genome`` its complement:

.. includefrags:: core/demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: rooted-iterators

Some iterators support iteration in reverse order with :dox:`BidirectionalIteratorConcept#goPrevious` as you can see in the next example.
Note that :dox:`BidirectionalIteratorConcept#goPrevious` is called before the value of ``it2`` is accessed.
Remember that the end position of a container is always the position behind the last item in the container.

.. includefrags:: core/demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: iterator-reverse

:dox:`RandomAccessContainerConcept#assignValue` can be used to change the value of an iterator.

.. includefrags:: core/demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: assign-value

The output of the program is as follows.

.. code-block:: console

   TATANNNGCGCG
   TATANNNGCGCG
   ATATNNNCGCGC
   GCGCGNNNATAT
   NATANNNGCGCG

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Copy the code below, which replaces all N's of a given :dox:`String` with A's.
     Adjust the code to use iterators to traverse the container.
     Use the :dox:`ContainerIteratorTags#Standard` iterator.

     .. code-block::cpp

        #include <iostream>
        #include <seqan/sequence.h>
        #include <seqan/file.h>

        using namespace seqan;

        int main()
        {
            Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
            for (unsigned i = 0; i < length(genome); ++i){
                if (genome[i] == 'N')
                    genome[i] = 'A';
            }
            std::cout << "Modified genome: " << genome << std::endl;
            return 0;
        }

    Solution

      Click **more...** to see the solution.

      .. container:: foldable

         .. includefrags:: core/demos/tutorial/iterators/iterators_assignment_1_solution.cpp

Assignment 2
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     Use the code from above and change the :dox:`ContainerIteratorTags#Standard` to a :dox:`ContainerIteratorTags#Rooted` iterator.
     Try to shorten the code wherever possible.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: core/demos/tutorial/iterators/iterators_assignment_2_solution.cpp

Workshop Assignment 3
^^^^^^^^^^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     In this assignment, we pick up the example from the workshop assignments from the sequences tutorial.
     Take the last solution and change the code to use Iterators.
     First, use Standard Iterators to do this.

     .. code-block:: cpp

        #include <iostream>
        #include <seqan/sequence.h>
        #include <seqan/file.h>

        using namespace seqan;
        // Function to print simple alignment between two sequences with the same length
        template <typename TText1, typename TText2>
        void printAlign(TText1 const & genomeFragment, TText2 const & read)
        {
                std::cout <<  "Alignment " << std::endl;
                std::cout << "  genome : " << genomeFragment << std::endl;
                std::cout << "  read   : " << read << std::endl;
        }

        int main(int, char const **)
        {
            // Build reads and genomes
            DnaString chr1 = "TATAATATTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGC"
                             "TCTGCGATATATCGCGCTAGATGTGCAGCTCGATCGAATGCACGTGT"
                             "GTGCGATCGATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTAT"
                             "CGGACGATCATATTAGCGGTCTAGCATTTAG";

            // Build List containing all reads
            typedef String<DnaString> DnaList;
            DnaList readList;
            resize(readList, 4);
            readList[0] = "TTGCTATCGCGATATCGCTAGCTAGCTACGGATTATGCGCTCTGCGATATATCGCGCT";
            readList[1] = "TCGATTAGCGTCGATCATCGATCTATATTAGCGCGCGGTATCGGACGATCATATTAGCGGTCTAGCATT";
            readList[2] = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGTTCATGTGCGCTGAAGCACACATGCACA";
            readList[3] = "CGTGCACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGTACTGCTGCTGACA";

            // Append a second chromosome sequence fragment to chr1
            DnaString chr2 = "AGCCTGCGTACGTTGCAGTGCGTGCGTAGACTGTTGCAAGCCGGGGGT"
                             "TCATGTGCGCTGAAGCACACATGCACACGTCTCTGTGTTCCGACGTGT"
                             "GTCACGTGCACTGCTGACGTCGTGGTTGTCACATCGTCGTGCGTGCGT"
                             "ACTGCTGCTGACACATGCTGCTG";
            append(chr1, chr2);

            // Print readlist
            std::cout << " \n Read list: " << std::endl;
            for(unsigned i = 0; i < length(readList); ++i)
                std::cout << readList[i] << std::endl;

            // Assume we have mapped the 4 reads to chr1 (and chr2) and now have the mapping start positions (no gaps).
            // Store the start position in a String: 7, 100, 172, 272
            String<unsigned> alignPosList;
            resize(alignPosList, 4);
            alignPosList[0] = 7;
            alignPosList[1] = 100;
            alignPosList[2] = 172;
            alignPosList[3] = 272;

            // Print alignments using Segment
            std::cout << " \n Print alignment using Segment: " << std::endl;
            for(unsigned i = 0; i < length(readList); ++i)
            {
                // Begin and end position of a given alignment between the read and the genome
                unsigned beginPosition = alignPosList[i];
                unsigned endPosition = beginPosition + length(readList[i]);
                // Build infix
                Infix<DnaString>::Type genomeFragment = infix(chr1, beginPosition, endPosition);
                // Call of our function to print the simple alignment
                printAlign(genomeFragment, readList[i]);
            }

            // Iterators :)
            // Print alignments using Iterators: Do the same as above, but use Iterators to iterate over the read list.
            // First, use Standard Iterators: Build two iterators it and itEnd to traverse readList.

            std::cout << " \n Print alignment using Standard Iterators: " << std::endl;

            return 1;
        }

   Solution
     Click **more...** to see the solution

     .. container:: foldable

        .. includefrags:: core/demos/tutorial/iterators/iterators_assignment_3_workshop_solution.cpp

Workshop Assignment 4
^^^^^^^^^^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Now, use rooted iterators in the example from Workshop ASsignment 3.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: core/demos/tutorial/iterators/iterators_assignment_4_workshop_solution.cpp
