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

.. includefrags:: demos/tutorial/iterators/base.cpp
    :fragment: construction

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

.. includefrags:: demos/tutorial/iterators/base.cpp
    :fragment: use-case

.. includefrags:: demos/tutorial/iterators/base.cpp.stdout
    :fragment: use-case

A Working Example
~~~~~~~~~~~~~~~~~

Let us now clarify the usage of iterators with a working example.
The following program demonstrates the usage of iterators.

.. includefrags:: demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: includes

The metafunction :dox:`ContainerConcept#Iterator` returns the iterator type for a given container type.
In this case the default implementation :dox:`ContainerIteratorTags#Standard` is used.

.. includefrags:: demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: metafunctions

We can use iterators to iterate over the elements of a container, e.g.  to print the elements.

.. includefrags:: demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: iterators

Instead of comparing the two iterators ``it`` and ``itEnd``, we could also use the function :dox:`RootedIteratorConcept#atEnd` to check whether we reached the end of the container.

.. includefrags:: demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: standard-iterators

Next we will use :dox:`RootedIteratorConcept Rooted Iterators`.
Since :dox:`RootedIteratorConcept Rooted Iterators` know their container, the functions :dox:`RootedRandomAccessIteratorConcept#goBegin` and :dox:`RootedIteratorConcept#atEnd` do not need to get the container as an argument.
The following example prints for each element of the :dox:`Dna5String Dna5 String` ``genome`` its complement:

.. includefrags:: demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: rooted-iterators

Some iterators support iteration in reverse order with :dox:`BidirectionalIteratorConcept#goPrevious` as you can see in the next example.
Note that :dox:`BidirectionalIteratorConcept#goPrevious` is called before the value of ``it2`` is accessed.
Remember that the end position of a container is always the position behind the last item in the container.

.. includefrags:: demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: iterators-reverse

:dox:`RandomAccessContainerConcept#assignValue` can be used to change the value of an iterator.

.. includefrags:: demos/tutorial/iterators/sequence_iterator_demo.cpp
   :fragment: assign-value

The output of the program is as follows.

.. includefrags:: demos/tutorial/iterators/sequence_iterator_demo.cpp.stdout

Assignment 1
^^^^^^^^^^^^

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

        .. includefrags:: demos/tutorial/iterators/assignment_2_solution.cpp

Workshop Assignment 3
^^^^^^^^^^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     In this assignment, we pick up the example from the workshop assignments from the sequences tutorial.
     Take the last solution and change the code to use Iterators.
     First, use Standard Iterators to do this.

     .. includefrags:: demos/tutorial/iterators/assignment_3_workshop.cpp

   Solution
     Click **more...** to see the solution

     .. container:: foldable

        .. includefrags:: demos/tutorial/iterators/assignment_3_workshop_solution.cpp

Workshop Assignment 4
^^^^^^^^^^^^^^^^^^^^^

.. container:: assignment

   Type
     Review

   Objective
     Now, use rooted iterators in the example from Workshop Assignment 3.

   Solution
     Click **more...** to see the solution.

     .. container:: foldable

        .. includefrags:: demos/tutorial/iterators/assignment_4_workshop_solution.cpp
