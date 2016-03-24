.. sidebar:: ToC

    .. contents::

.. _tutorial-how-to-recipes-metafunctions:

Metafunctions
=============

Type Metafunctions
------------------

For example, the metafunction :dox:`ContainerConcept#Iterator` is a type metafunction, i.e. it is used to determine a type.
Type metafunctions have the form:

``typename TypeMetaFunc<T1, T2, ..., TN>::Type``

``TypeMetaFunc``
  The name of the metafunction

``T1, T2, ..., TN``
  Arguments (types or constants)

``Type``
  The resulting type

The keyword ``typename`` must be stated if one of the arguments ``T1, T2, ..., TN`` is or uses a template parameter.
For example the following piece of code uses the metafunction ``Iterator`` to determine an iterator type for a string class:

.. includefrags:: demos/tutorial/metafunctions/base.cpp
    :fragment: iterator

.. includefrags:: demos/tutorial/metafunctions/base.cpp.stdout
    :fragment: iterator

Value Metafunctions
-------------------

Metafunctions can also be used to determine constant values at compile time.
The general form of value metafunctions is:

``VALUE_META_FUNC<T1, T2, ..., TN>::VALUE``

``VALUE_META_FUNC``
  The name of the metafunction

``T1, T2, ..., TN``
  Arguments (types or constants)

``VALUE``
  The resulting constant value

For example the following function prints the length of a fixed sized string using the value metafunction :dox:`LENGTH`:

.. includefrags:: demos/tutorial/metafunctions/base.cpp
    :fragment: length

.. includefrags:: demos/tutorial/metafunctions/base.cpp.stdout
    :fragment: length

.. important::
      Different uses of "*Value*":

      Please note that :dox:`ContainerConcept#Value` (``Value<TSomeType>::Type``) is a **Type Metafunction**, because it returns a Type (e.g. of values in a container) and not a value.

Assignment 1
""""""""""""

.. container:: assignment

   Objective
     Write a generic function ``checkContainerForDna(T & container)`` that prints out a message if the value inside the container is of the type :dox:`Dna`. The type ``T`` of the container should be specified as a template argument. Test you function with some examples.

   Hint
      * Use the **Type Metafunction** :dox:`ContainerConcept#Value` to access the (alphabet-)type of the elements in the container.
      * Use the **Value Metafunction** :dox:`IsSameType` to check for type equality.

   Solution
     .. container:: foldable

        Your program should look something like this:

        .. includefrags:: demos/tutorial/metafunctions/assignment1_solution.cpp

        Note: Because the Value Metafunction ``IsSameType<>`` is evaluated at compile time, the part of the if-statement code that does not apply won't even appear in the compiled code. This can be an improvement to the runtime of your code.

        The output is the following:

        .. includefrags:: demos/tutorial/metafunctions/assignment1_solution.cpp.stdout

