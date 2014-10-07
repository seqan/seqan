.. sidebar:: ToC

   .. contents::


.. _tutorial-basic-techniques:

Basic Techniques
----------------

.. todo::

   Here should come an easy introduction into STL code. Has low priority

Generic Programming
~~~~~~~~~~~~~~~~~~~

SeqAn adopts **generic programming**, a paradigm that was proven to be an efficient design strategy in the C++ standard.
The standard template library (STL) as part of the C++ standard is a prototypical example for generic programming.
Generic programming designs algorithms and data structures in a way that they work on all types that meet a minimal set of requirements.
An example for a generic data structure in the STL is the class ``vector``: It is a container for storing objects of a type ``T`` that are assignable, which means that we can assign one instance ``s`` of ``T`` to another instance ``t`` of ``T``, i.e. the code ``T t = s`` is valid.

This kind of requirement to the interface of a type ``T`` is called a ``concept``, and we say that a type ``T`` *implements* a concept, if it fulfills all requirements stated by that concept; for example the concept assignable is implemented by all built-in types and every class that has both a copy assignment operator and a copy constructor.

Generic programming has two implications:

#. Data structures and algorithms work on *all* types ``T`` that implement the relevant concept, i.e. relevant is not the type ``T`` itself but its interface, and
#. this concept is *minimal* in the sense that it contains only those requirements that are essential for the data structure or algorithm to work on ``T``.

This way data structures and algorithms can be applied to as many types as possible, and hence generic programming promotes the generality of the library.

Generic data types and algorithms can be implemented in C++ using templates.
A class template parameterizes a class with a list of types or constants.
For example, a declaration for the class ``vector`` could be:

.. code-block:: cpp

   template <typename T> class vector;

where ``T`` stands for the *value type*, i.e. the type of the values that will be stored in ``vector``.
The template is generic, it can be applied to any type ``T``.
For example, a vector for storing ``int`` values is instantiated by:

.. code-block:: cpp

   vector<int> my_vector;

That is we use ``int`` as template argument for ``T``, and the result of the instantiation is an object ``my_vector`` of the complete type ``vector<int>``.
The compiler employs the same template, i.e. the same piece of code, for different template argument types.
The compilation succeeds if the applied template argument type supports all uses of the parameter ``T`` within the template code, so the C++ template instantiation process implies the minimality of the concepts.

Further reading
^^^^^^^^^^^^^^^

For more information about generic programming and the STL we recommend reading:

* Vandervoorde, Josuttis: C++ Templates - The complete guide, Addison-Wesley
