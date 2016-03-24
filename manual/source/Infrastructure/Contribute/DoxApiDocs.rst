.. sidebar:: ToC

    .. contents::

.. _infra-contribute-dox:

API Documentation System (dox)
==============================

Since the 1.4.1 release, SeqAn uses a new documentation system. The
syntax is similar to `Doxygen <http://doxygen.sf.net>`_ but slightly
different. The main differences are (1) not identifying functions by
their signatures but only by their names, (2) adding the idea of
metafunctions, (3) adding the idea of interface functions and (4) an
extension to SeqAn-specific things like documenting concepts.

General Documentation Structure
-------------------------------

Dox comments are placed in C-style comments with an exclamation mark (see below).
The first dox tag should be placed on the next line, each line should begin with a correctly indented star.
The first line only contains the slash-star-exclamation-mark and the last line only contains the star-slash.

.. code-block:: cpp

   /*!
    * @fn myFunction
    * @signature void myFunction()
    */

The documentation and the code are independent. Each item to be documented (adaption, class, concept, enum, function, group, macro, metafunction, page, tag, typedef, variable) has to be epxlicitely given (see tags below). The available top level tags are [#adaption @adaption], [#class @class], [#concept @concept], [#defgroup @defgroup], [#enum @enum], [#fn @fn], [#macro @macro], [#metafunction @mfn], [#page @page], [#tag @tag], [#typedef @typedef], and [#variable @var].

Each top-level tag creates a documentation entry. For example, the
following defines a class ``Klass`` with two global interface functions
``f1`` and ``f2`` for this class:

.. code-block:: cpp

    /*!
     * @class Klass
     * @fn Klass#f1
     * @fn Klass#f2
     */

Member functions are given using ``::``, the same as in the C++
language:

.. code-block:: cpp

    /*!
     * @class Klass
     * @fn Klass::memberFunc
     */

Global interface functions are global functions that belong to the
interface of a type. Similar, interface metafunctions are metafunctions
that belong to the interface of a type. Their fully qualified name for
dox consists of the type name, followed by a hash ``#`` and the
function/metafunction name:

.. code-block:: cpp

    /*!
     * @class Klass
     * @fn Klass#interfaceFunc
     * @mfn Klass#InterfaceMetaFunc
     */

Below the top-level tags, come the second-level tags. The first kind of
second-level tags defines properties of an entry. Important such
second-level entries are ``@brief``, ``@signature``, ``@see``,
``@param``, ``@tparam``, ``@return``. You can also write text for the
description of your entity and use tags such as ``@section``,
``@subsection``, ``@snippet``, ``@code`` to format the description. You
can use HTML tags for formatting the documentation.

Example:

.. code-block:: cpp

    /*!
     * @class Align
     * @brief Store a tabular alignment.
     *
     * @signature template <typename TSource, typename TSpec>
     *            class Align;
     *
     * @tparam TSource The type of the underlying sequence.
     * @tparam TSpec   Tag for selecting the specialization of the Align class.
     *
     * The <tt>Align</tt> class provides a tabular alignment of sequences with the
     * same type.  The sequences are given with <tt>TSource</tt>.  An <tt>Align</tt>
     * object will use a <a href="seqan:Gaps">Gaps</a> object for each sequence.
     * The specialization of the <a href="seqan:Gaps">Gaps</a> object can be selected
     * using the <tt>TSpec</tt> template parameter.
     *
     * @see Gaps
     * @see globalAlignment
     */

Images are included using ``<img src="${PATH}">`` where ``${PATH}`` is
relative to the source image directory.

Tag Documentation
-----------------

Below, we differentiate between **names** and **labels**.

**Names** are used to identify documentation items and must follow
extended C++ identifier rules. A sub name consists of only alphanumeric
characters and the underscore is allowed, must not start with a number.
Sub names can be glued together with ``::`` for class members and ``#``
for interface functions. In contrast, **labels** are used for the
display to the user. For example, the alloc string has the name
``AllocString`` but the label "Alloc String", the constructor of
``AllocString`` has name ``AllocString::String``, and its length
function has name ``AllocString#length``.

@adaption
^^^^^^^^^

**Signature** ``@adaption AdaptionName [Adaption Label]``

Top-level tag.

Defines an adaption with the given name and an optional label.

An adaption is a collection of global interface functions and
metafunctions that adapt a type outside the SeqAn library to a concept
in the SeqAn library. For example, the STL ``std::string`` class can be
adapted to the interface of the ``StringConcept`` concept.

.. code-block:: cpp

    /*!
     * @adaption StdStringToStringConcept std::string to Sequence concept
     * @brief The <tt>std::string</tt> class is adapted to the Sequence concept.
     */

@aka
^^^^

**Signature** ``@aka OtherName``

Second-level entry.

Assigns an alias name for a function, metafunction, class, concept, or
enum. The list of aliases will be printed for each code entry. Also, the
aliases will be incorporated into search results.

.. code-block:: cpp

    /*!
     * @class InfixSegment
     * @brief Represents a part of a string.
     *
     * @aka substring
     */

    template <typename TSequence>
    class InfixSegment<TSequence, Infix>;

@brief
^^^^^^

**Signature** ``@brief Brief description.``

Second-level tag.

Defines the brief description of the top-level entry it belongs to. You
can use HTML in the description.

.. code-block:: cpp

    /*!
     * @fn f
     * @brief A minimal function.
     * @signature void f();
     */

    void f();

@class
^^^^^^

**Signature** ``@class ClassName [Class Label]``

Top-level tag.

Defines a class with the given name ``ClassName`` and an optional label.

.. code-block:: cpp

    /*!
     * @class AllocString Alloc String
     * @extends String
     * @brief Implementation of the String class using dynamically allocated array.
     *
     * @signature template <typename TAlphabet, typename TSpec>
     * class String<TAlphabet, Alloc<TSpec> >;
     * @tparam TAlphabet Type of the alphabet (the string's value).
     * @tparam TSpec     Tag for the further specialization.
     */

    template <typename TAlphabet, typename TSpec>
    class String<TAlphabet, Alloc<TSpec> >
    {
        // ...
    };

@code
^^^^^

**Signature** ``@code{.ext} ... @endcode``

Second-level tag.

Provides the means to display code blocks in the documentation. The
extension ``.ext`` is used for identifying the type (use ``.cpp`` for
C++ code) and selecting the appropriate highlighting.

.. code-block:: cpp

    /*!
     * @fn f
     * @brief Minimal function.
     * @signature void f();
     *
     * @code{.cpp}
     * int main()
     * {
     *     f();  // Call function.
     *     return 0;
     * }
     * @endcode
     */

    void f();

Note that you can use the extension value ``.console`` to see console output.

.. code-block:: cpp

   /*!
    * @fn f
    * @brief Some function
    *
    * @section Examples
    *
    * @include demos/module/demo_f.cpp
    *
    * The output is as follows:
    *
    * @code{.console}
    * This is some output of the program.
    * @endcode
     */

@concept
^^^^^^^^

**Signature** ``@concept ConceptName [Concept Label]``

Top-level tag.

Creates a documentation entry for a concept with the given name and an
optional label. All concept names should have the suffix ``Concept``.
Use the fake keyword ``concept`` in the ``@signature``.

A concept is the C++ equivalent to interfaces known in other classes.
C++ provides no real way for concepts so at the moment they are a formal
construct used in the documentation.

.. code-block:: cpp

    /*!
     * @concept StringConcept Sequence
     * @signature concept StringConcept;
     * @extends ContainerConcept
     * @brief Concept for sequence types.
     */

@defgroup
^^^^^^^^^

**Signature** ``@defgroup GroupName [Group Label]``

Top-level tag.

Creates a documentation entry for a group with a given name and an
optional label. Groups are for rough grouping of global functions and/or
tags.

You can put types and functions into a group similar to making global
interface functions and metafunctions part of the interface of a class
or concept.

.. code-block:: cpp

    /*!
     * @defgroup FastxIO FASTA/FASTQ I/O
     * @brief Functionality for FASTA and FASTQ I/O.
     *
     * @fn FastxIO#readRecord
     * @brief Read one record from FASTA/FASTQ files.
     *
     * @fn FastxIO#writeRecord
     * @brief Write one record to FASTA/FASTQ files.
     *
     * @fn FastxIO#readBatch
     * @brief Read multiple records from FASTA/FASTQ file, limit to a given count.
     *
     * @fn FastxIO#writeBatch
     * @brief Write multiple records to FASTA/FASTQ file, limit to a given count.
     *
     * @fn FastxIO#readAll
     * @brief Read all records from a FASTA/FASTQ file.
     *
     * @fn FastxIO#writeAll
     * @brief Write all records to a FASTA/FASTQ file.
     */

@deprecated
^^^^^^^^^^^

**Signature** ``@deprecated message``

Second-level entry.

Marks a given function, metafunction, class, concept, or enum as
deprecated. A deprecation message will be generated in the API
documentation.

.. code-block:: cpp

    /*!
     * @fn f
     * @deprecated Use @link g @endlink instead.
     * @brief Minimal function.
     */

    void f();

@enum
^^^^^

**Signature** ``@enum EnumName [Enum Label]``

Top-level entry.

Provides documentation for an enum with given name and optional label.

.. code-block:: cpp

    /*!
     * @enum MyEnum
     * @brief An enum.
     *
     * @val MyEnum VALUE1
     * @brief VALUE1 value of enum MyEnum.
     *
     * @val MyEnum VALUE2
     * @brief VALUE2 value of enum MyEnum.
     */

    enum MyEnum
    {
      VALUE1,
      VALUE2
    };

@extends
^^^^^^^^

**Signature** ``@extends OtherName``

Gives a parent class for a given class or a parent concept for a given
concept.

.. code-block:: cpp

    /*!
     * @concept OneConcept
     *
     * @concept TwoConcept
     * @extends OneConept
     *
     * @class MyClass
     *
     * @class OtherClass
     * @extends MyClass
     */

@fn
^^^

**Signature** ``@fn FunctionName [Function Label]``

Top-level entry.

Documents a function (global, global interface, or member) with given
name and label. The type of the function is given by its name.

.. code-block:: cpp

    /*!
     * @fn globalAlignment
     * @brief Pairwise, DP-based global alignment.
     */

@headerfile
^^^^^^^^^^^

**Signature** ``@headerfile path``

Second-level entry.

Gives the required ``#include`` path for a code entry.

**Note:** Use angular brackets as below for SeqAn includes.

.. code-block:: cpp

    /*!
     * @fn f
     * @brief A minimal function.
     * @headerfile <seqan/module.h>
     */

@implements
^^^^^^^^^^^

**Signature** ``@implements ConceptName``

Second-level entry.

Marks a class to implement a given concept.

.. code-block:: cpp

    /*!
     * @concept MyConcept
     *
     * @class ClassName
     * @implements MyConcept
     */

@include
^^^^^^^^

**Signature** ``@include path/to/file``

Second-level entry.

Includes a C++ source file as an example. See [#snippet @snippet] for
including fragments.

.. code-block:: cpp

    /*!
     * @fn f
     * @brief Minimal function.
     *
     * The following example shows the usage of the function.
     * @include demos/use_f.cpp
     */

@internal
^^^^^^^^^

**Signature** ``@internal [ignored comment``

Second-level entry.

Marks a given function, metafunction, class, concept, or enum as
internal. You can also provide a comment that is ignored/not used in the
output.

.. code-block:: cpp

    /*!
     * @fn f
     * @internal
     * @brief Minimal function.
     */

    void f();

@link
^^^^^

**Signature** ``@link TargetName target label``

In-text tag.

Provides tag to link to a documentation entry with a given label.

The difference to [#see @see] is that ``@link .. @endlink`` is used
inline in text whereas ``@see`` is a second-level tag and adds a ``see``
property to the documented top-level entry. Use ``@link`` to link to
entries within the documentation and the HTML ``<a>`` tag to link to
external resources.

.. code-block:: cpp

    /*!
     * @fn f
     * @brief Minimal function.
     *
     * The function is mostly useful with the @link String string class@endlink.
     */

@macro
^^^^^^

**Signature** ``@macro MacroName [Macro Label]``

Top-level tag.

Documents a macro.

.. code-block:: cpp

    /*!
     * @macro MY_MACRO
     * @brief Multiply two values.
     *
     * @signature #define MY_MACRO(i, j) ...
     * @param i A value for i.
     * @param j A value for j.
     * @return The product of i and j: (i * j)
     */

    #define MY_MACRO(i, j) (i * j)

@mfn
^^^^

**Signature** ``@mfn MetafunctionName [Metafunction Label]``

Top-level tag.

Documents a metafunction.

.. code-block:: cpp

    /*!
     * @mfn Identity
     * @brief Identity function for types.
     *
     * @signature Identity<T>::Type
     * @tparam T The type to pass in.
     * @returns The type T.
     */

    template <typename T>
    struct Identity
    {
        typedef T Type;
    };

@note
^^^^^

**Signature** ``@note message``

Second-level entry.

Adds an informative note to a function, metafunction, class, concept,
enum, or group.

.. code-block:: cpp

    /*!
     * @fn f
     * @note Very useful if used together with @link g @endlink.
     * @brief Minimal function.
     */

    void f();

@page
^^^^^

**Signature** ``@page PageName [Page Title]``

Top-level entry.

Creates a documentation page.

.. code-block:: cpp

    /*!
     * @page SomePage Page Title
     *
     * A very simple page
     *
     * @section Section
     *
     * A section!
     *
     * @subsection Subsection
     *
     * A subsection!
     */

@param
^^^^^^

**Signature** ``@param Name Label``

Second-level entry.

Documents a value (and non-type) parameter from a function or member
function.

.. code-block:: cpp

    /*!
     * @fn square
     * @brief Compute the square of an <tt>int</tt> value.
     *
     * @signature int square(x);
     * @param x The value to compute square of (type <tt>int</tt>).
     * @return int The square of <tt>x</tt>.
     */

    int square(int x);

@return
^^^^^^^

**Signature** ``@return Type Label``

Defines the return value for a function or metafunction.

Also see the example for [#param @param].

When documenting functions and the result type is the result of a
metafunction then use a ``TXyz`` return type in ``@return`` and document
``TXyz`` in the text of ``@return`` as follows:

.. code-block:: cpp

    /*!
     * @fn lengthSquare
     * @brief Compute the square of the length of a container.
     *
     * @signature TSize square(c);
     *
     * @param c The container to compute the squared length of.
     * @return TSize squared length of <tt>c</tt>.  <tt>TSize</tt> is the size type of <tt>c</tt>.
     */

    template <typename TContainer>
    typename Size<TContainer>::Type lengthSquare(TContainer const & c);

@throw
^^^^^^

**Signature** ``@return Exception Label``

Adds a note on a function or macro throwing an exception.

.. code-block:: cpp

    /*!
     * @fn myFunction
     * @brief Writes things to a file.
     * @signature void myFunction(char const * filename);
     *
     * @param[in] filename File to write to.
     *
     * @throw std::runtime_error If something goes wrong.
     */
    void myFunction(char const * filename);

@datarace
^^^^^^^^^

**Signature** ``@datarace Description``

Describes possible data races for functions and macros.
If this value is not specified it defaults to ``Thread safety unknown!``

.. code-block:: cpp

    /*!
     * @fn myFunction
     * @brief Writes things to a file.
     * @signature void myFunction(char const * filename);
     *
     * @param[in] filename File to write to.
     *
     * @datarace This function is not thread safe and concurrent writes to the file might invalidate the output.
     */
    void myFunction(char const * filename);


@section
^^^^^^^^

**Signature** ``@section Title``

Second-level entry.

Adds a section to the documentation of an entry.

See the example for [#page @page].

@see
^^^^

**Signature** ``@see EntryName``

Second-level entry.

Adds "see also" link to a documentation entry.

.. code-block:: cpp

    /*!
     * @fn f
     * @brief A simple function.
     *
     * Here is a snippet:
     *
     * @snippet demos/use_f.cpp Simple Function
     */

And here is the file with the snippet.

.. code-block:: cpp

    /* Some code */

    int main(int argc, char const ** argv)
    {
    //![Simple Function]
        return 0;
    //![Simple Function]
    }

    /* Some more code */

@tag
^^^^

**Signature** ``@tag TagName``

Top-level entry.

Documents a tag. Mostly, you would group tags in a group using [#defgroup
@defgroup].

.. code-block:: cpp

    /*!
     * @defgroup MyTagGroup My Tag Group
     *
     * @tag MyTagGroup#TagName
     * @tag MyTagGroup#MyOtherTagName
     */

@tparam
^^^^^^^

**Signature** ``@tparam TArg``

Second-level entry.

Documents a template parameter of a metafunction or class template.

.. code-block:: cpp

    /*!
     * @mfn MetaFunc
     * @signature MetaFunc<T1, T2>::Type
     *
     * @tparam T1 First type.
     * @tparam T2 Second type.
     */

@typedef
^^^^^^^^

**Signature** ``@typedef TypedefName``

Top-level entry.

Documents a typedef.

.. code-block:: cpp

    /*!
     * @typedef CharString
     * @brief An AllocString of character.
     *
     * @signature typedef String<char, Alloc<> > CharString;
     */

@var
^^^^

**Signature** ``@var VariableType VariableName``

Top-level entry. Document a global variable or member variable.

.. code-block:: cpp

    /*!
     * @class MyClass
     *
     * @var int MyClass::iVar
     */

    class MyClass
    {
    public:
        int iVar;
    };

@val
^^^^

**Signature** ``@val EnumType EnumValueName``

Top-level entry.
Documents an enum value.

.. code-block:: cpp

    /*!
     * @enum EnumName
     * @brief My enum.
     * @signature enum EnumName;
     *
     * @val EnumName::VALUE1;
     * @brief The first enum value.
     *
     * @val EnumName::VALUE2;
     * @brief The second enum value.
     */

    enum MyEnum
    {
        VALUE1,
        VALUE2
    };

@warning
^^^^^^^^

**Signature** ``@warning message``

Second-level entry.

Adds a warning to a function, metafunction, class, concept, enum, or group.

.. code-block:: cpp

    /*!
     * @fn f
     * @note Using this function can lead to memory leaks.  Try to use @link g @endlink instead.
     * @brief Minimal function.
     */

    void f();

Best Practice
-------------

This section describes the best practice when writing documentation.

Clarifying Links
^^^^^^^^^^^^^^^^

Our usability research indicates that some functionality is confusing
(e.g. see #1050) but cannot be removed. One example is the function
``reserve()`` which can be used to *increase* the *capacity* of a
container whereas the function ``resize()`` allows to change the *size*
of a container, *increasing or decreasing* its size.

The documentation of such functions should contain a clarifying text and
a link to the other function.

.. code-block:: cpp

    /*!
     * @fn Sequence#reserve
     *
     * Can be used to increase the <b>capacity</b> of a sequence.
     *
     * Note that you can only modify the capacity of the sequence.  If you want to modify the
     * <b>length</b> of the sequence then you have to use @link Sequence#resize @endlink.
     */

Documentation Location
^^^^^^^^^^^^^^^^^^^^^^

**Add the documentation where it belongs.** For example, when
documenting a class with multiple member functions, put the dox comments
for the class before the class, the documentation of the member
functions in front of the member functions. For another example, if you
have to define multiple signatures for a global interface function or
metafunctions, put the documentation before the first function.

.. code-block:: cpp

    /*!
     * @class Klass
     * @brief A class.
     */
    class Klass
    {
    public:
        /*!
         * @var int Klass::x
         * @brief The internal value.
         */
        int x;

        /*!
         * @fn Klass::Klass
         * @brief The constructor.
         *
         * @signature Klass::Klass()
         * @signature Klass::Klass(i)
         * @param i The initial value for the member <tt>x</tt> (type <tt>int</tt>).
         */
        Klass() : x(0)
        {}

        Klass(int x) : x(0)
        {}

        /*!
         * @fn Klass::f
         * @brief Increment member <tt>x</tt>
         * @signature void Klass::f()
         */
        void f()
        {
            ++x;
        }
    };

Signatures
^^^^^^^^^^

Always document the return type of a function. If it is the result of a
metafunction or otherwise depends on the input type, use ``TResult`` or
so and document it with ``@return``.

HTML Subset
-----------

You can use inline HTML to format your description and also for creating
links.

*  Links into the documentation can be generated using ``<a>`` if the scheme in ``href`` is ``seqan:``: ``<a href="seqan:AllocString">the alloc string</a>.``
*  Use ``<i>`` for italic/emphasized text.
*  Use ``<b>`` for bold text.
*  Use ``<tt>`` for typewriter text.

Tag Ordering
^^^^^^^^^^^^

**Keep consistent ordering of second-level tags.** The following order
should be used, i.e. if several of the following tags appear, they
should appear in the order below.

#. ``@internal``
#. ``@deprecated``
#. ``@warning``
#. ``@note``
#. ``@brief``
#. ``@extends``
#. ``@implements``
#. ``@signature``
#. ``@param``
#. ``@tparam``
#. ``@return``
#. ``@headerfile``
#. The documentation body with the following tags in any order (as fit for the documentation text) and possibly interleaved with text:
   ``@code``, ``@snippet``, ``@include``, ``@section``, ``@subsection``.
#. ``@see``
#. ``@aka``

Documenting Concepts
^^^^^^^^^^^^^^^^^^^^

All concepts should have the suffix ``Concept``.

Use the pseudo keyword ``concept`` in the ``@signature``.

Use the following template:

.. code-block:: cpp

    /*!
     * @concept MyConcept
     * @brief The concept title.
     *
     * @signature concept MyConcept;
     *
     * The concept description possibly using include, snippet, and <b><i>formatting</i></b> etc.
     */

Documenting Classes
^^^^^^^^^^^^^^^^^^^

Use the following template:

.. code-block:: cpp

    /*!
     * @class AllocString Alloc String
     * @brief A string storing its elements on dynamically heap-allocated arrays.
     *
     * @signature template <typename TAlphabet, typename TSpec>
     * class AllocString<TAlphabet, Alloc<TSpec> >;
     * @tparam TAlphabet The alphabet/value type to use.
     * @tparam TSpec    The tag to use for further specialization.
     *
     * The class description possibly using include, snippet, and <b><i>formatting</i></b> etc.
     */

Documenting Functions
^^^^^^^^^^^^^^^^^^^^^

Use the following template:

.. code-block:: cpp

    /*!
     * @fn globalAlignment
     * @brief Global DP-based pairwise alignment.
     *
     * @signature TScore globalAlignment(align, scoringScheme);
     * @signature TScore globalAlignment(align, scoringScheme, lowerBand, upperBand);
     * @param align Align object to store the result in. Must have length 2 and be filled with sequences.
     * @param scoringScheme Score object to use for scoring.
     * @param lowerBand The lower band of the alignment (<tt>int</tt>).
     * @param upperBAnd The upper band of the alignment (<tt>int</tt>).
     * @return TScore The alignment score of type <tt>Value<TScore>::Type</tt> where <tt>TScore</tt> is the type of <tt>scoringScheme</tt>.
     *
     * The function description possibly using include, snippet, and <b><i>formatting</i></b> etc.
     */

Documenting Metafunctions
^^^^^^^^^^^^^^^^^^^^^^^^^

Use the following template:

.. code-block:: cpp

    /*!
     * @mfn Size
     * @brief Return size type of another type.
     *
     * @signature Size<T>::Type
     * @tparam T The type to query for its size type.
     * @return TSize The size type to use for T.
     *
     * The class description possibly using include, snippet, and <b><i>formatting</i></b> etc.
     */

Documenting Enums
^^^^^^^^^^^^^^^^^

.. code-block:: cpp

    /*!
     * @enum EnumName
     * @brief My enum.
     * @signature enum EnumName;
     *
     * @var EnumName::VALUE
     * @summary The enum's first value.
     *
     * @var EnumName::VALUE2
     * @summary The enum's second value.
     */

Difference to Doxygen
---------------------

If you already know Doxygen, the following major differences apply.

* The documentation is more independent of the actual code.
  Doxygen creates a documentation entry for all functions that are present in the code and allows the additional documentation, e.g. using ``@fn`` for adding functions.
  With the SeqAn dox system, you have to explicitely use a top level tag for adding documentation items.
* Documentation entries are not identified by their signature but by their name.
* We allow the definition of interface functions and metafunctions (e.g. ``@fn Klass#func`` and ``@mfn Klass#Func``) in addition to member functions (``@fn Klass::func``).
* We do not allow tags with backslashes but consistently use at signs (``@``).
