.. sidebar:: ToC

   .. contents::


.. _tutorial-template-subclassing:

Generic Programming
===================

A generic algorithm that is applicable to a type ``T`` needs not to be optimal for that type.
The algorithm find in the standard library (ISO/IEC 1998, section 25.3.1.1) for example performs a sequential linear time search and is therefore capable of finding a given value in any standard compliant container.
However, the container map was designed to support a faster logarithmic time search, so the algorithm find – though applicable – is not optimal for searching in map.
This shows that sometimes a special algorithm could be faster than a generic algorithm.
Hence, in order to achieve better performance, SeqAn supports refinements of algorithms.
A special version is only useful if it really allows a speedup in some cases, and only in this case it will actually be implemented.
Therefore we assume that for a given case always the most special applicable variant is the best, where we have to assure that there is always a definite most special candidate according to the C++ function overload resolution rules (ISO/IEC 1998, sections 13.3 and 14.5.8).
We can write ``find(obj)`` for any container type ``obj``, and this invokes the most suitable implementation of ``find`` depending on the type of ``obj``.
We call this approach **template subclassing**.

The technique of template subclassing may be summarized as follows:

*  The data types are realized as default implementation or specialization of class templates, e.g., ``Class``, which have at least one template parameter ``TSpec``.
*  Refinements of ``Class`` are specified by using in ``TSpec`` a tag class, e.g., ``Subclass``, that means they are implemented as class template specializations ``Class<Subclass>``.
*  Whenever further refinements may be possible, we declare the tag classes as class templates with at least one template parameter ``TSpec``, in which more tag classes can be used.
   For example we may implement a class template specialization ``Class<Subclass<Subsubclass<...> > >``.
   This way, we can reach arbitrary levels of specialization.
*  Algorithms can be implemented for each level of specialization.
   If multiple implementations for different levels of specialization exist, then the C++ function overload resolution selects the most special from all applicable variants.

Example: Generic q-gram hashing
-------------------------------

In many applications in bioinformatics so called q-grams are used.
A short string of length q can be interpreted as a number to the base of the cardinality of the alphabet.
So for example for the Dna alphabet ``cgta=0*1+3*4+2*16+1*64=108``.
q-grams can be gapped or ungapped.
In the gapped case they are often called :dox:`Shape shapes` and are simply an ordered list of integers.
The number of integers in the list is called the ``size`` of a shape whereas the largest element -1 is called the ``span``.

The following code sniplet shows a generic algorithm for computing all hash values for a shape.
The function ``span`` applied to the shape s = ⟨s1, . . . , sq⟩ returns sq − 1.

.. code-block:: cpp

   template <typename TShape, typename TString> void hashAll(TShape & shape, TString & str)
      typedef typename Iterator<TString>::Type TIterator;
      TIterator it = begin(str);
      TIterator it_end = end(str) - span(shape);
      while (it != it_end) {
         unsigned int hash_value = hash(shape, it);
        /* do some things with the hash value */ ++it;
      }

Each shape has to know the alphabet :math:`\Sum`, so we specify this value type in the first template parameter of ``Shape``.
The actual specialization is selected in the second template parameter ``TSpec``.

.. code-block:: cpp

    template <typename TValue, typename TSpec = SimpleShape> class Shape;

the default is :dox:`SimpleShape` which is simply an ungapped shape storing merely the length of the shape from which it can deduce its span and size.

.. code-block:: cpp

   template <typename TValue> class Shape< TValue, SimpleShape >
   {
      public:
         unsigned int span;
   };

If we know q at compile time, then we can specify it in a template parameter and define span as a static member:

.. code-block:: cpp

   template <unsigned int q = 0> struct UngappedShape<q>;

   template <typename TValue, unsigned int q> class Shape< TValue, UngappedShape<q> >
   {
      public:
         static unsigned int const span = q;
   };

The question is now, whether we can speed up the above ``hashAll`` functions for specializations of the class ``shape`` like ungapped shapes.
A little thinking yields a positive answer to that question.
Indeed, for ungapped shapes, we can incrementally compute the next hash value form a given hashvalue in constant time using the formula ``hash(a_{i+1}...a{_i+q})=hash(a_{i}...a_{i+q−1})|Σ|−a_{i}|Σ|^q +a_{i+q``}, that means when shifting the shape along the sequence, we have to subtract the effect of the leftmost letter and add the effect of the rightmost letter, all scaled with the corresponding factor. All digits in between are *shifted* by multiplying them with the alphabet size.
Obviously this allows for a much more efficient implementation of the ``hashAll`` function. This functionality can be encoded in the following function :dox:`Shape#hashNext`.

.. code-block:: cpp

   template <typename TValue, unsigned int q, typename TIterator>
   inline unsigned int
   hashNext(Shape< TValue, UngappedShape<q> > const & shape, TIterator it, unsigned int prev)
   {
      unsigned int val = prev * ValueSize<TValue>::VALUE - *it * shape.fac
                                  + *(it + shape.span);
      return val;
    // shape.fac stores |Σ|^q
   }

SeqAn aims at not using virtual functions for introducing polymorphism.
Instead the concept is called ``template subclassing``.
Hence we can now define a specialized ``hashAll`` function for all ungapped shapes as follows.

.. code-block:: cpp

   template <typename TValue, unsigned int q, typename TString>
   void hashAll(Shape< TValue, UngappedShape<q> > & shape, TString & str)
      typedef typename Iterator<TString>::Type TIterator;
      TIterator it = begin(str); TIterator it_end = end(str) - span(shape);
      unsigned int hash_value = hash(shape, it);
      /* do some things with the hash value */

      while (++it != it_end) {
         unsigned int hash_value = hashNext(shape, it, hash_value);
         /* do some things with the hash value */
      }
   }

Thats pretty much it.
The C++ resolution mechanisms will ensure that whenever you use an ungapped shape in your code, the more efficient ``hashAll`` function above will be compiled.
Note that this decision is made at *compile time* as opposed to the virtual function mechanism which is invoked at *run time*.

Further reading
---------------

For more information about generic programming and the STL we recommend reading.

* Vandervoorde, Josuttis: C++ Templates - The complete guide, Addison-Wesley

Template Subclassing Demo
-------------------------

Here is an example of template subclassing.

.. includefrags:: core/demos/tutorial/generic_programming/template_subclassing.cpp
