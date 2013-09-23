// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Default implementations for sequences.
//
// TODO(holtgrew): There appears to be some overlap with string_base.h. Maybe it is a good idea to move everything related to strings to string_base.h and remove all default-container behaviour.
// TODO(holtgrew): Each value is a container by itself. This is highly undesirable since it introduces the set-of-sets problem and when users confuse atomic values with containers, bugs are hard to find. This feature should be removed.
// ==========================================================================

// TODO(holtgrew): These functions have (documentation wise) mostly gone into Container, Sequence Concepts and String class.  This is where they belong.

#ifndef SEQAN_SEQUENCE_SEQUENCE_INTERFACE_H_
#define SEQAN_SEQUENCE_SEQUENCE_INTERFACE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @defgroup OverflowStrategyTags
 * @brief The strategy for resizing containers.
 * 
 * @section Remarks
 * 
 * Changing the capacity of a container can invalidate the iterators of this container.
 * 
 * If no overflow tag is specified, most operations use the default overflow strategy given by @link
 * DefaultOverflowImplicit @endlink or @link DefaultOverflowExplicit @endlink, depending on the kind of operation.
 * 
 * @see computeGenerousCapacity
 * @see DefaultOverflowImplicit
 * @see DefaultOverflowExplicit
 * 
 * @tag OverflowStrategyTags#Limit
 * @headerfile <seqan/sequence.h> 
 * @brief Limit the contents to current capacity.
 * 
 * @typedef Tag<TagLimit_> Limit;
 *
 * @section Remarks
 * 
 * All entries that exceed the capacity are lost.
 * @headerfile <seqan/sequence.h> 
 * @tag OverflowStrategyTags#Generous
 * 
 * @brief Expand if needed, get precautionary extra space.
 * 
 * @typedef Tag<TagGenerous_> Generous;
 *
 * @section Remarks
 * 
 * Whenever the capacity has to be increased, the new capacity is choosen somewhat large than actually needed.  This
 * strategy limits the number of capacity changes, so that resizing takes armotized constant time.  Use this strategy if
 * the total amount of storage is unkown at first.
 * 
 * The new capacity is computed by @link computeGenerousCapacity @endlink. By default, it is guaranteed not to exceed
 * about  tree halfs of the space that is used to store the data.  The user can overload @link computeGenerousCapacity
 * @endlink in order to change this behavior.
 * 
 * @tag OverflowStrategyTags#Exact
 * @headerfile <seqan/sequence.h> 
 * @brief Expand as far as needed.
 * 
 * @typedef Tag<TagExact_> Exact;
 *
 * @section Remarks
 * 
 * The capacity is only changed if the current capacity is not large enough.  If the capacity can only be expanded up to
 * a certain ammount, it will be increased as far as possible  and the contents are limited to the new capacity.
 * 
 * Note that the capacity will never be shrinked.  Use @link ContainerConcept#shrinkToFit @endlink to resize the
 * capacity down to the current length.
 * 
 * @tag OverflowStrategyTags#Insist
 * @headerfile <seqan/sequence.h> 
 * @brief No capacity check.
 *
 * @typedef Tag<TagInsist_> Insist;
 * 
 * @section Remarks
 * 
 * The user has to ensure that the container's capacity is large enough.
 */

/**
.Tag.Overflow Strategy:
..cat:Sequences
..summary:The strategy for resizing containers.
..tag.Insist:No capacity check.
...remarks:The user has to ensure that the container's capacity is large enough.
..tag.Limit:Limit the contents to current capacity.
...remarks: All entries that exceed the capacity are lost.
..tag.Exact:Expand as far as needed.
...remarks: The capacity is only changed if the current capacity is not large enough.
 If the capacity can only be expanded up to a certain ammount, it will be increased as far as possible
 and the contents are limited to the new capacity.
...remarks:Note that the capacity will never be shrinked.
 Use @Function.shrinkToFit@ to resize the capacity down to the current length.
..tag.Generous:Expand if needed, get precautionary extra space.
...remarks:Whenever the capacity has to be increased, the new capacity is choosen somewhat large than actually needed.
 This strategy limits the number of capacity changes, so that resizing takes armotized constant time.
 Use this strategy if the total amount of storage is unkown at first.
...remarks:The new capacity is computed by @Function.computeGenerousCapacity@.
By default, it is guaranteed not to exceed about
 tree halfs of the space that is used to store the data.
 The user can overload @Function.computeGenerousCapacity@ in order to change this behavior.
..remarks:Changing the capacity of a container can invalidate the iterators of this container.
..remarks:If no overflow tag is specified, most operations use the default overflow strategy given by @Metafunction.DefaultOverflowImplicit@
or @Metafunction.DefaultOverflowExplicit@, depending on the kind of operation.
..include:seqan/sequence.h
*/
struct TagInsist_;
typedef Tag<TagInsist_> Insist;
typedef Tag<TagInsist_> Tight;  // TODO(holtgrew): Necessary?

struct TagLimit_;
typedef Tag<TagLimit_> Limit;

struct TagGenerous_;
typedef Tag<TagGenerous_> Generous;

struct TagExact_;
typedef Tag<TagExact_> Exact;

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// --------------------------------------------------------------------------

/*!
 * @mfn DefaultOverflowImplicit
 * @headerfile <seqan/sequence.h>
 * @brief The default overflow strategy for implicit resize.
 *
 * @signature DefaultOverflowImplicit<T>::Type;
 *
 * @tparam T The type to get the default overflow tag for.
 *
 * @return Type The default overflow tag.  The default implementation returns <tt>Generous</tt>.
 *
 * @section Remarks
 *
 * This function is used for functions that cause an implicit change of a container's size, like e.g. assign, append,
 * and replace.
 */

/**
.Metafunction.DefaultOverflowImplicit:
..hidefromindex
..class:Class.String
..summary:The default overflow strategy for implicit resize.
..signature:DefaultOverflowImplicit<T>::Type
..param.T:Type for which the overflow strategy is determined.
...type:Class.String
..returns.param.Type:Expansion tag for type of $T$.
..remarks:This function is used for functions that cause an implicit change of a container's size, like
e.g. @Function.assign@, @Function.append@, and @Function.replace@.
..see:Tag.Overflow Strategy
..include:seqan/sequence.h
*/
template <typename T>
struct DefaultOverflowImplicit
{
    typedef Generous Type;
};

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowExplicit
// --------------------------------------------------------------------------

/*!
 * @mfn DefaultOverflowExplicit
 * @headerfile <seqan/sequence.h>
 * @brief The default overflow strategy for explicit resize.
 *
 * @signature DefaultOverflowExplicit<T>::Type;
 *
 * @tparam T The type to determine overflow strategy.
 *
 * @return Type The resulting expantion tag for <tt>T</tt>.
 *
 * @section Remarks
 *
 * This function is used for functions that change a container's size explicit, like e.g. resize.
 */

/**
.Metafunction.DefaultOverflowExplicit:
..hidefromindex
..class:Class.String
..summary:The default overflow strategy for explicit resize.
..signature:DefaultOverflowExplicit<T>::Type
..param.T:Type for which the overflow strategy is determined.
...type:Class.String
..returns.param.Type:Expansion tag for type of $T$.
..remarks:This function is used for functions that change a container's size explicit, like e.g. @Function.resize@.
..see:Tag.Overflow Strategy
..see:Metafunction.DefaultOverflowImplicit
..include:seqan/sequence.h
*/
template <typename T>
struct DefaultOverflowExplicit
{
    typedef Generous Type;
};

// --------------------------------------------------------------------------
// Metafunction IsContiguous
// --------------------------------------------------------------------------

/*!
 * @mfn IsContiguous
 * @headerfile <seqan/sequence.h>
 * @brief Determines whether a container stores its elements contiguously in memory.
 *
 * @signature IsContiguous<T>::Type;
 * @signature IsContiguous<T>::VALUE;
 *
 * @tparam T The type that is tested for being a string.
 *
 * @return Type  Either <tt>True</tt> or <tt>False</tt>, depending on whether <tt>T</tt> is stored contiguously.
 * @return VALUE Either <tt>true</tt> or <tt>false</tt>, depending on whether <tt>T</tt> is stored contiguously.
 *
 * @section Remarks
 *
 * A sequence container is "contiguous", if its elements are stored in a single contiguous array.  Examples for
 * contiguous sequences are AllocString or char arrays.
 *
 * If an object <tt>obj</tt> is a contiguous sequence, then <tt>begin(obj)</tt> can be converted to a pointer to the
 * first element of the content array.
 */

/**
.Metafunction.IsContiguous:
..cat:Sequences
..summary:Determines whether a container stores its elements in a contiguous array.
..signature:IsContiguous<T>::VALUE
..param.T:Type that is tested for being a string.
..returns.param.VALUE:$true$ if $T$ is a string, $false$ otherwise.
..remarks:Definition: A sequence container is "contiguous", if its elements
    are stored in a single contiguous array.
    Examples for contiguous sequences are @Spec.Alloc String@ or @Adaption.char array@.
..remarks:If an object $obj$ is a contiguous sequence, then $begin(obj)$ can be
    converted to a pointer to the first element of the content array.
..include:seqan/sequence.h
*/
template <typename T>
struct IsContiguous
{
    typedef False Type;
    enum { VALUE = false };
};

template <typename T>
struct IsContiguous<T const>
    : public IsContiguous<T> {};

// --------------------------------------------------------------------------
// Metafunction IsSequence
// --------------------------------------------------------------------------

// TODO(holtgrew): Deprecate in favour of Is<SequenceConcept>?

/*!
 * @mfn IsSequence
 * @headerfile <seqan/sequence.h>
 * @brief Determines whether a type is a sequence.
 *
 * @signature IsSequence<T>::Type;
 * @signature IsSequence<T>::VALUE;
 *
 * @tparam T The type to query.
 *
 * @return Type  <tt>True</tt> if <tt>T</tt> is a sequence and <tt>False</tt> otherwise.
 * @return VALUE <tt>true</tt> if <tt>T</tt> is a sequence and <tt>false</tt> otherwise.
 *
 * @section Remarks
 *
 * For example, String and Segment as <tt>T</tt> return true.
 */

/**
.Metafunction.IsSequence:
..cat:Sequences
..summary:Determines whether a container stores its elements in sequential order.
..signature:IsSequence<T>::VALUE
..param.T:Type that is tested for being a sequence.
..returns.param.VALUE:$true$ if $T$ is a sequence, $false$ otherwise.
..remarks:For example @Class.String@ and @Class.Segment@ return $true$.
..include:seqan/sequence.h
*/
template <typename T>
struct IsSequence
{
    typedef False Type;
    enum { VALUE = false };
};

template <typename T>
struct IsSequence<T const>
    : public IsSequence<T> {};

// --------------------------------------------------------------------------
// Metafunction AllowsFastRandomAccess
// --------------------------------------------------------------------------

/*!
 * @mfn AllowsFastRandomAccess
 * @headerfile <seqan/sequence.h>
 * @brief Determines whether a sequence efficiently supports random access.
 *
 * @signature AllowsFastRandomAccess<T>::Type;
 * @signature AllowsFastRandomAccess<T>::VALUE;
 *
 * @tparam T The type to query.
 *
 * @return Type  <tt>True</tt> if <tt>T</tt> allows for fast random access and <tt>False</tt> otherwise.
 * @return VALUE <tt>true</tt> if <tt>T</tt> allows for fast random access and <tt>false</tt> otherwise.
 *
 * @section Remarks
 *
 * For example, String and std::vector allow for fast random access, while std::list does not.
 */

/**
.Metafunction.AllowsFastRandomAccess:
..cat:Sequences
..summary:Determines whether a sequence efficiently supports random access.
..signature:AllowsFastRandomAccess<T>::VALUE
..param.T:Type that is tested for fast random access.
..returns.param.VALUE:$true$ if $T$ supports fast random access, $false$ otherwise.
..remarks:For example @Spec.Alloc String@, @Class.Segment@, and @Spec.Block String@ return $true$.
..include:seqan/sequence.h
*/
template <typename T>
struct AllowsFastRandomAccess
{
    typedef True Type;
    enum { VALUE = true };
};

template <typename T>
struct AllowsFastRandomAccess<T const>
    : public AllowsFastRandomAccess<T> {};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function getObjectId()
// --------------------------------------------------------------------------

/*!
 * @fn getObjectId
 * @header <seqan/sequence.h>
 * @brief A value that identifies the underlying sequence.
 *
 * @signature TVoidPtr getObjectId(object);
 *
 * @param object The object for which to determine the id.
 *
 * @return TVoidPtr a <tt>void const *</tt> value identying the object.
 *
 * @section Remarks
 *
 * Two sequences should have the same id, if they share the same resource, e.g. the same memory buffer.
 *
 * The exact semantic of the returned id can vary for different classes.  Typically, the id of a string is a <tt>void
 * const *</tt> to the end of the string.
 *
 * @section Examples
 *
 * @code{.cpp}
 * String<char> str = "hallo seqan";
 * bool b1 = (getObjectId(str) == getObjectId(infix(str, 3, 7));   //true
 * bool b2 = (getObjectId(str) == getObjectId(String<char>(str))); //false
 * bool b3 = (getObjectId(str) == getObjectId(toCString(str)));
 * @endcode
 *
 * In this example, <tt>b1</tt> is <tt>true</tt., since the segment object returned by <tt>infix()</tt> is just a filter
 * and uses the buffer of it's host object str.
 *
 * <tt>String&lt;char&gt;(str)</tt> constructs a temporary copy of <tt>str</tt>, so these two strings have different id values.
 *
 * The result of the last comparison depends on the implementation of <tt>toCString</tt> and cannot be predicted at
 * compile time.
 */

/**
.Function.getObjectId:
..cat:Miscellaneous
..class:Class.String
..summary:A value that identifies the underlying sequence.
..signature:void const * getObjectId(object)
..param.object:The object for which the id will be determined.
..returns:The id of $sequence$.
..remarks.text:Two sequences should have the same id, if they share the same resource, e.g. the same memory buffer.
..remarks.text:The exact semantic of the returned id can vary for different classes.
Typically, the id of a string is a $void const *$ to the end of the string.
..remarks.note:The id of a single character need not to be the id of its container.
..example.code:String<char> str = "hallo seqan";
bool b1 = (getObjectId(str) == getObjectId(infix(str, 3, 7));   //true
bool b2 = (getObjectId(str) == getObjectId(String<char>(str))); //false
bool b3 = (getObjectId(str) == getObjectId(toCString(str)));
..example.text:In this example, $b1$ is $true$, since the segment object returned by $infix()$
is just a filter and uses the buffer of it's host object $str$.
..example.text:$String<char>(str)$ constructs a temporary copy of $str$, so these two
strings have different id values.
..example.text:The result of the last comparison depends on the implementation of $toCString$
and cannot be predicted at compile time.
..include:seqan/sequence.h
*/
template <typename T>
inline void const *
getObjectId(T const & me)
{
    SEQAN_CHECKPOINT;
    return end(me, Standard());
}

// --------------------------------------------------------------------------
// Function shareResources()
// --------------------------------------------------------------------------

/*!
 * @fn shareResources
 * @headerfile <seqan/sequence.h>
 * @brief Determines whether two sequences share the same resource.
 * 
 * @signature bool shareResources(s1, s2);
 *
 * @param[in] s1 First sequence.
 * @param[in] s2 Second sequence.
 *
 * @return bool <tt>true</tt> if the two sequences share resources and <tt>false</tt> if not.
 */

/**
.Function.shareResources:
..cat:Miscellaneous
..class:Class.String
..summary:Determines whether two sequences share the same resource.
..signature:bool shareResources(sequence1, sequence2)
..param.sequence1, sequence2:Two sequences.
..returns:$false$ if it can be guaranteed that $sequence1$ and $sequence2$ can be modified without changing each other, $true$ otherwise.
..remarks:Non-sequences are interpreted as sequences of size 1.
..remarks:Note that this function may not work properly for argument types that are not listed here.
..include:seqan/sequence.h
*/

template <typename T1, typename T2>
inline bool
shareResources(T1 const & obj1,
               T2 const & obj2)
{
    SEQAN_CHECKPOINT;
    return getObjectId(obj1) == getObjectId(obj2);
}

// --------------------------------------------------------------------------
// Function _beginDefault()
// --------------------------------------------------------------------------

//* ???Anti Default Sequences
// TODO(holtgrew): Evil -- each value is a container of length 1.
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, Standard>::Type
_beginDefault(T & me,
               Standard)
{
    SEQAN_CHECKPOINT;
    return & me;
}
// TODO(holtgrew): Evil -- each value is a container of length 1.
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Standard>::Type
_beginDefault(T const & me,
               Standard)
{
    SEQAN_CHECKPOINT;
    return & me;
}
//*/

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, Rooted>::Type
_beginDefault(T & me,
               Rooted)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<T, Rooted>::Type TIterator;
    return TIterator(me, begin(me, Standard()));
}
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Rooted>::Type
_beginDefault(T const & me,
               Rooted)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<T const, Rooted>::Type TIterator;
    return TIterator(me, begin(me, Standard()));
}

// --------------------------------------------------------------------------
// Function begin()
// --------------------------------------------------------------------------

/*!
 * @fn begin
 * @headerfile <seqan/sequences.h>
 * @brief Returns the begin iterator of a container.
 *
 * @signature TIterator begin(object[, tag]);
 *
 * @param[in] object The container to get the begin iterator for.
 * @param[in] tag    The tag to use for picking the resulting type.
 *
 * @return TIterator An iterator to the first item in <tt>object</tt>. The type is the result of
 *                   <tt>Iterator<TContainer, TTag>::Type</tt> where <tt>TContainer</tt> is the type of
 *                   <tt>object</tt> and <tt>TTag</tt> is the type of <tt>tag</tt>.
 *
 * @section Remarks
 *
 * If the container is empty, the resulting iterator does not point to a valid element in the container.
 */

/**
.Function.begin:
..cat:Iteration
..cat:Containers
..class:Class.String
..class:Concept.ContainerConcept
..summary:The begin of a container.
..signature:Iterator begin(object [, tag])
..param.object:A container.
...type:Class.String
...concept:Concept.ContainerConcept
..param.tag:An @Tag.Iterator Spec.iterator spec@ tag that specifies the kind of the iterator returned. (optional)
...default:Given by @Metafunction.DefaultGetIteratorSpec@.
..returns:An iterator to the first item in $object$. The type is the result of $Iterator<TContainer, TTag>::Type$ where $TContainer$ is the type of $object$ and $TTag$ is the type of $tag$.
...metafunction:Metafunction.Iterator
..remarks.text:If the container does not contain any items at all, the function may return 0.
..see:Function.end
..see:Metafunction.Iterator
..include:seqan/sequence.h
*/
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type
begin(T & me)
{
    SEQAN_CHECKPOINT;
    return begin(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type
begin(T const & me)
{
    SEQAN_CHECKPOINT;
    return begin(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

//folgende forward Deklaration wurde wegen Phaenomene bei VC++ 2003 hinzugenommen
//implemented in string_pointer.h
template <typename TValue>
SEQAN_HOST_DEVICE inline typename Iterator<TValue const *, Standard>::Type
begin(TValue const * me,
      Standard);

template <typename T, typename TSpec>
SEQAN_HOST_DEVICE inline typename Iterator<T, Tag<TSpec> const>::Type
begin(T & me,
      Tag<TSpec> const tag_)
{
    SEQAN_CHECKPOINT;
    return _beginDefault(me, tag_);
}
template <typename T, typename TSpec>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Tag<TSpec> const>::Type
begin(T const & me,
      Tag<TSpec> const tag_)
{
    SEQAN_CHECKPOINT;
    return _beginDefault(me, tag_);
}

/*
template <typename TValue>
inline typename Iterator<TValue *, Standard>::Type
begin(TValue * me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return me;
}

//folgende Version wurde wegen eines seltsamen Phaenomens bei VC++ hinzugenommen
template <typename TValue>
inline typename Iterator<TValue const *, Standard>::Type
begin(TValue const * me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return me;
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue *, Standard>::Type
begin(TValue * me,
      Tag<TSpec> const tag_)
//    Standard)
{
    SEQAN_CHECKPOINT;
    return me;
}

//folgende Version wurde wegen eines seltsamen Phaenomens bei VC++ hinzugenommen
template <typename TValue, typename TSpec>
inline typename Iterator<TValue const *, Standard>::Type
begin(TValue const * me,
      Tag<TSpec> const tag_)
//    Standard)
{
    SEQAN_CHECKPOINT;
    return me;
}
*/

// --------------------------------------------------------------------------
// Function beginPosition()
// --------------------------------------------------------------------------

/*!
 * @fn beginPosition
 * @headerfile <seqan/sequence.h>
 * @brief Begin position of an object in a host.
 *
 * @signature TPosition beginPosition(object);
 *
 * @param[in] object An object.
 *
 * @return TPosition The position of the first item in <tt>host(object)</tt> that belongs to <tt>object</tt>.
 *                   The type TPosition is the position type of <tt>object</tt>.
 *
 * @section Remarks
 *
 * For most classes, this function returns 0.  Exceptions are Segments.
 *
 * @section Examples
 * CharString str = "ABCDEF";
 * std::cout &lt;&lt; beginPosition(str) &lt;&lt; std::endl;
 * 
 * Infix&lt;CharString&gt;::Type myInfix = infix(str, 1, 5);
 * std::cout &lt;&lt; beginPosition(myInfix) &lt;&lt; std::endl;
 * @code{.cpp}
 * @endcode
 */

/**
.Function.beginPosition:
..cat:Containers
..summary:Begin position of object in host.
..signature:Position beginPosition(object)
..class:Class.String
..class:Concept.ContainerConcept
..param.object:An object.
...type:Class.String
...concept:Concept.ContainerConcept
..returns:The position of the first item in $host(object)$ that belongs of $object$.
...metafunction:Metafunction.Position
..remarks
...text:For most classes $beginPosition$ always returns 0. Exceptions are e.g. @Spec.InfixSegment@ and @Spec.SuffixSegment@.
..see:Function.begin
..include:seqan/sequence.h
..example.code:
CharString str = "ABCDEF";
std::cout << beginPosition(str) << std::endl;

Infix<CharString >::Type myInfix = infix(str, 1, 5);
std::cout << beginPosition(myInfix) << std::endl;
*/
template <typename T>
inline typename Position<T>::Type
beginPosition(T &)
{
    SEQAN_CHECKPOINT;
    return 0;
}

template <typename T>
inline typename Position<T>::Type
beginPosition(T const &)
{
    SEQAN_CHECKPOINT;
    return 0;
}

// --------------------------------------------------------------------------
// Function _endDefault()
// --------------------------------------------------------------------------

//* ???Anti Default Sequences
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, Standard>::Type
_endDefault(T & me,
             Standard)
{
    SEQAN_CHECKPOINT;
    return (& me) + 1;
}
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Standard>::Type
_endDefault(T const & me,
             Standard)
{
    SEQAN_CHECKPOINT;
    return (& me) + 1;
}
//*/

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, Rooted>::Type
_endDefault(T & me,
             Rooted)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<T, Rooted>::Type TIterator;
    return TIterator(me, end(me, Standard()));
}
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Rooted>::Type
_endDefault(T const & me,
             Rooted)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<T const, Rooted>::Type TIterator;
    return TIterator(me, end(me, Standard()));
}

// --------------------------------------------------------------------------
// Function end()
// --------------------------------------------------------------------------

/*!
 * @fn end
 * @headerfile <seqan/sequence.h>
 * @brief Return iterator to the end of a container.
 *
 * @signature TIterator end(object[, tag]);
 *
 * @param[in] object The container to get the end iterator for.
 * @param[in] tag    The tag to use for picking the resulting type.
 *
 * @return TIterator An iterator to the first item in <tt>object</tt>. The type is the result of
 *                   <tt>Iterator<TContainer, TTag>::Type</tt> where <tt>TContainer</tt> is the type of
 *                   <tt>object</tt> and <tt>TTag</tt> is the type of <tt>tag</tt>.
 */

/**
.Function.end:
..cat:Iteration
..cat:Containers
..class:Class.String
..concept:Concept.ContainerConcept
..summary:The end of a container.
..signature:Iterator end(object [, tag])
..param.object:A container.
...type:Class.String
...concept:Concept.ContainerConcept
..param.tag:An @Tag.Iterator Spec.iterator spec@ tag that specifies the kind of the iterator returned. (optional)
...default:Given by @Metafunction.DefaultGetIteratorSpec@.
..returns:An iterator that points behind the last item in $object$. The type is the result of $Iterator<TContainer, TTag>::Type$ where $TContainer$ is the type of $object$ and $TTag$ is the type of $tag$.
...metafunction:Metafunction.Iterator
..remarks.text:If the container does not contain any items at all, the function may return 0.
..see:Function.begin
..see:Metafunction.Iterator
..include:seqan/sequence.h
*/
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type
end(T & me)
{
    SEQAN_CHECKPOINT;
    return end(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type
end(T const & me)
{
    SEQAN_CHECKPOINT;
    return end(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

template <typename T, typename TSpec>
SEQAN_HOST_DEVICE inline typename Iterator<T, Tag<TSpec> const>::Type
end(T & me,
    Tag<TSpec> const tag_)
{
    SEQAN_CHECKPOINT;
    return _endDefault(me, tag_);
}

template <typename T, typename TSpec>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Tag<TSpec> const>::Type
end(T const & me,
    Tag<TSpec> const tag_)
{
    SEQAN_CHECKPOINT;
    return _endDefault(me, tag_);
}

// --------------------------------------------------------------------------
// Function endPosition()
// --------------------------------------------------------------------------

/*!
 * @fn endPosition
 * @headerfile <seqan/sequence.h>
 * @brief End position of an object in its host.
 *
 * @signature TPosition endPosition(object);
 *
 * @param[in] object The object to query for its end position.
 *
 * @return TPosition The position behind the last item in <tt>host(object)</tt> that belongs to <tt>object</tt>.
 */

/**
.Function.endPosition:
..cat:Containers
..class:Class.String
..concept:Concept.ContainerConcept
..summary:End position of object in host.
..signature:Position endPosition(object)
..param.object:An object.
...concept:Concept.ContainerConcept
...type:Class.String
..returns:The position behind the last item in $host(object)$ that belongs of $object$.
...metafunction:Metafunction.Position
..see:Function.end
..see:Function.beginPosition
..include:seqan/sequence.h
*/
template <typename T>
inline typename Position<T>::Type
endPosition(T & me)
{
    SEQAN_CHECKPOINT;
    return length(me);
}

template <typename T>
inline typename Position<T>::Type
endPosition(T const & me)
{
    SEQAN_CHECKPOINT;
    return length(me);
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

/*!
 * @fn value
 * @headerfile <seqan/sequence.h>
 * @brief Reference to the value.
 *
 * @signature TReference value(container, position);
 *
 * @param[in] container A container of values.
 * @param[in] position  The position to query the value for.
 *
 * @return TReference A reference or proxy to the value.
 */

/**
.Function.value:
..cat:Iteration
..cat:Containers
..concept:Concept.ContainerConcept
..summary:Reference to the value.
..signature:Reference value(container, position)
..param.container:A container of values.
..param.position:A position in $container$ on which the value should be accessed.
..returns:A reference or proxy to the value.
...metafunction:Metafunction.Reference
..include:seqan/sequence.h
*/
//* ???Anti Default Sequences
template <typename T, typename TPos>
SEQAN_HOST_DEVICE inline typename Reference<T>::Type
value(T & me,
      TPos /*pos*/)
{
    SEQAN_CHECKPOINT;
    return me;
}

template <typename T, typename TPos>
SEQAN_HOST_DEVICE inline typename Reference<T const>::Type
value(T const & me,
      TPos /*pos*/)
{
    SEQAN_CHECKPOINT;
    return me;
}
//*/

// --------------------------------------------------------------------------
// Function getValue()
// --------------------------------------------------------------------------

/*!
 * @fn getValue
 * @headerfile <seqan/sequence.h>
 * @brief Access a value.
 *
 * @signature TGetValue getValue(container, pos);
 *
 * @param[in] container The container to query for a value.
 * @param[in] pos       The position to query the value at.
 *
 * @return TGetValue The get-value in <tt>container</tt> at position <tt>pos</tt>.
 *
 * @section Remarks
 *
 * The get value is always convertible into the Value type of the container.  It can either be a temporary copy of the
 * given item or a const reference into the container.
 *
 * If <tt>pos</tt> is outside of <tt>container</tt> then the behaviour of the function is undefined.
 */

/**
.Function.getValue:
..summary:Access to the value.
..cat:Containers
..cat:Content Manipulation
..concept:Concept.ContainerConcept
..signature:GetValue getValue(container, pos)
..param.container:A container.
...concept:Concept.ContainerConcept
..param.pos:The position of an item in $object$.
...remarks:$pos$ should be convertible to $Position<T>::Type$ for $container$-type $T$.
..returns:The item at position $pos$ in $container$.
This can either be a reference to the item or a temporary copy of the item.
...metafunction:Metafunction.GetValue
..remarks:
...text:If $pos$ is out of range, then the behavior of the function is undefined.
..see:Metafunction.GetValue
..see:Metafunction.Position
..see:Function.value
..include:seqan/sequence.h
*/
template <typename T, typename TPos>
inline typename GetValue<T>::Type
getValue(T & me,
         TPos pos)
{
    SEQAN_CHECKPOINT;
    return (typename GetValue<T>::Type) value(me, pos);
}

template <typename T, typename TPos>
inline typename GetValue<T const>::Type
getValue(T const & me,
         TPos pos)
{
    SEQAN_CHECKPOINT;
    return value(me, pos);
}

// --------------------------------------------------------------------------
// Function front()
// --------------------------------------------------------------------------

/*!
 * @fn front
 * @headerfile <seqan/sequence.h>
 * @brief Return the first item of a container.
 *
 * @signature TReference front(container);
 *
 * @param[in] container The container to query.
 *
 * @return TReference Reference to the first item in <tt>container</tt>.
 * 
 * @section Remarks
 *
 * This function is equivaeltn to <tt>value(container, 0)</tt>.
 */

/**
.Function.front:
..cat:Containers
..concept:Concept.ContainerConcept
..summary:The first item in container.
..signature:Iterator front(container)
..param.container:A container.
...concept:Concept.ContainerConcept
..returns:A @Metafunction.Reference.reference@ of the first item in $container$.
...metafunction:Metafunction.Reference
..remarks:This function is equivalent to $value(me, 0)$.
..see:Function.value
..see:Function.begin
..include:seqan/sequence.h
*/

template <typename T>
inline typename Reference<T>::Type
front(T & me)
{
    SEQAN_CHECKPOINT;
    return value(me, 0);
}
template <typename T>
inline typename Reference<T const>::Type
front(T const & me)
{
    SEQAN_CHECKPOINT;
    return value(me, 0);
}

// --------------------------------------------------------------------------
// Function back()
// --------------------------------------------------------------------------

/*!
 * @fn back
 * @headerfile <seqan/sequence.h>
 * @brief Return the last item of a container.
 *
 * @signature TReference back(container);
 *
 * @param[in] container The container to query.
 *
 * @return TReference Reference to the last item in <tt>container</tt>.
 * 
 * @section Remarks
 *
 * This function is equivaeltn to <tt>value(container, length(container) - 1)</tt>.
 */

/**
.Function.back:
..cat:Containers
..concept:Concept.ContainerConcept
..summary:The last item in container.
..signature:Iterator back(container)
..param.container:A container.
...concept:Concept.ContainerConcept
..returns:A @Metafunction.Reference.reference@ of the last item in $container$.
...metafunction:Metafunction.Reference
..remarks:This function is equivalent to $value(me, length(me) - 1)$.
..see:Function.value
..see:Function.end
..see:Function.front
..include:seqan/sequence.h
*/

template <typename T>
SEQAN_HOST_DEVICE inline typename Reference<T const>::Type
back(T const & me)
{
    SEQAN_CHECKPOINT;
    return value(me, length(me) - 1);
}

template <typename T>
SEQAN_HOST_DEVICE inline typename Reference<T>::Type
back(T & me)
{
    SEQAN_CHECKPOINT;
    return value(me, length(me) - 1);
}

// --------------------------------------------------------------------------
// Function iter()
// --------------------------------------------------------------------------

/*!
 * @fn iter
 * @headerfile <seqan/sequence.h>
 * @brief Iterator to the item at the given position in a container.
 *
 * @signature TIterator iter(object, pos[, tag]);
 *
 * @param[in] object The container to get the iterator for.
 * @param[in] pos    The position to get the iterator for.
 * @param[in] tag    The tag to pick the type of the iterator.
 *
 * @return TIterator The resulting iterator.  If <tt>TTag</tt> is the type of <tt>tag</tt> and <tt>TContainer</tt> the
 *                   type of <tt>object</tt> then TIterator is of the type <tt>Iterator&lt;TContainer,
 *                   TTag&gt;::Type</tt>.
 *
 * @section Remarks
 *
 * If <tt>pos</tt> is out of range then the iterator is invalid.
 */

/**
.Function.iter:
..cat:Containers
..class:Class.String
..summary:Iterator to item at given position.
..signature:Iterator iter(object, pos [, tag])
..param.object:A container.
...type:Class.String
..param.pos:The position of an item in $object$.
...metafunction:Metafunction.Position
..param.tag:An @Tag.Iterator Spec.iterator spec@ tag that specifies the kind of the iterator returned. (optional)
...default:Given by @Metafunction.DefaultGetIteratorSpec@.
..returns:An iterator to the item at position $pos$ in $object$. The type is the result of $Iterator<TContainer, TTag>::Type$ where $TContainer$ is the type of $object$ and $TTag$ is the type of $tag$.
...metafunction:Metafunction.Iterator
..remarks:
...text:If $pos$ is out of range, then the behavior of the function is undefined.
..see:Function.value
..see:Metafunction.Iterator
..see:Metafunction.Position
..include:seqan/sequence.h
 */
template <typename T, typename TPos>
inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type
iter(T & me,
     TPos pos)
{
    SEQAN_CHECKPOINT;
    return iter(me, pos, typename DefaultGetIteratorSpec<T>::Type());
}

template <typename T, typename TPos>
inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type
iter(T const & me,
     TPos pos)
{
    SEQAN_CHECKPOINT;
    return iter(me, pos, typename DefaultGetIteratorSpec<T>::Type());
}

template <typename T, typename TPos, typename TTag>
inline typename Iterator<T, Tag<TTag> const>::Type
iter(T & me,
     TPos pos,
     Tag<TTag> const tag_)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LEQ_MSG(pos, static_cast<TPos>(length(me)), "Trying to get an iterator behind a container through iter().");
    return begin(me, tag_) + pos;
}

template <typename T, typename TPos, typename TTag>
inline typename Iterator<T const, Tag<TTag> const>::Type
iter(T const & me,
     TPos pos,
     Tag<TTag> const tag_)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LEQ_MSG(pos, static_cast<TPos>(length(me)), "Trying to get an iterator behind a container through iter().");
    return begin(me, tag_) + pos;
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

/*!
 * @fn assignValue
 * @headerfile <seqan/sequence.h>
 * @brief Assign a value of a container at a given position.
 *
 * @signature void assignValue(container, pos, value);
 *
 * @param[in,out] container The container to manipulate.
 * @param[in]     pos       The position of the item in the container to manipulate.
 * @param[in]     value     The value to assign to <tt>container[pos]</tt>.
 */

/**
.Function.assignValue:
..cat:Content Manipulation
..concept:Concept.ContainerConcept
..signature:assignValue(container, pos, value)
..param.container:A container.
...concept:Concept.ContainerConcept
..param.pos:Position of the item in $container$ to that $value$ is assigned.
..remarks:If $object$ is a container (that is $pos$ is not specified),
    the whole content of $object$ is replaced by $value$.
..remarks.text:
    If $value$ is not used again after calling this function,
    then consider to use @Function.moveValue@ that could be faster in some cases instead.
..include:seqan/sequence.h
*/

template <typename T, typename TValue, typename TPos>
inline void
assignValue(T & me,
            TPos pos,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    assign(value(me, pos), _value);
}

// --------------------------------------------------------------------------
// Function moveValue()
// --------------------------------------------------------------------------

/*!
 * @fn moveValue
 * @headerfile <seqan/sequence.h>
 * @brief Move a value of a container to a given position.
 *
 * @signature void assignValue(container, pos, value);
 *
 * @param[in,out] container The container to manipulate.
 * @param[in]     pos       The position of the item in the container to manipulate.
 * @param[in,out] value     The value to move to <tt>container[pos]</tt>.
 */

/**
.Function.moveValue:
..cat:Content Manipulation
..concept:Concept.ContainerConcept
..signature:moveValue(container, pos, value)
..param.object:
...concept:Concept.ContainerConcept
..param.container:A container.
...concept:Concept.ContainerConcept
..param.pos:Position of the item in $container$ to that $value$ is moved to.
..remarks:If $object$ is a container (that is $pos$ is not specified),
the whole content of $object$ is replaced by $value$.
..remarks.text:
    This function possibly clears $value$.
    If $value$ should be used further, consider to use @Function.assignValue@ instead.
..include:seqan/sequence.h
*/

template <typename T, typename TValue, typename TPos>
inline void
moveValue(T & me,
          TPos pos,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    move(value(me, pos), _value);
}

template <typename T, typename TValue, typename TPos>
inline void
moveValue(T const & me,
          TPos pos,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    move(value(me, pos), _value);
}

// --------------------------------------------------------------------------
// Function length()
// --------------------------------------------------------------------------

/*!
 * @fn length
 * @headerfile <seqan/sequence.h>
 * @brief Return the number of items in a container.
 *
 * @signature TSize length(object);
 *
 * @param[in] object The container to query.
 *
 * @return TSize The resulting size of the container.
 *
 * @section Remarks
 *
 * The length of a sequence can never exceed its capacity.
 *
 * @see capacity
 */

/**
.Function.length:
..cat:Containers
..concept:Concept.ContainerConcept
..summary:The number of items/characters.
..signature:Size length(object)
..param.object:A container.
...concept:Concept.ContainerConcept
..returns:The number of items/characters in $object$.
...metafunction:Metafunction.Size
..remarks.text:The length of a sequence can never exceed it's capacity.
..see:Function.capacity
..include:seqan/sequence.h
*/

//* ???Anti Default Sequences
template <typename T>
inline typename Size<T>::Type
length(T const & /*me*/)
{
    SEQAN_CHECKPOINT;
    return 1;
}
//*/

// --------------------------------------------------------------------------
// Function capacity()
// --------------------------------------------------------------------------

/*!
 * @fn capacity
 * @headerfile <seqan/sequence.h>
 * @brief Return the capacity of a container.
 *
 * @signature TSize capacity(object);
 *
 * @param[in] object The container to query for its capacity.
 *
 * @return TSize Returns the size of the container.
 *
 * @section Remarks
 *
 * The size of a sequence can never exceed its capacity but some container support resizing of the capacity.  Some
 * functions do that implicitely if they are called with a suitable @link OverflowStrategyTags tag @endlink.  The
 * function reserve can be used to change the capacity explicitely.
 */

/**
.Function.capacity:
..cat:Containers
..concept:Concept.ContainerConcept
..summary:The maximal length.
..signature:Size capacity(object)
..param.object:A container.
...remarks: If $object$ cannot be converted to one of these types, the function returns 1.
..returns:The maximal number of items/characters that can be stored in $object$.
...metafunction:Metafunction.Size
..remarks.text:The size of a sequence can never exceed it's capacity, but some containers support
resizing of the capacity.
Some functions do that implicitely if they are called with a suitable @Tag.Overflow Strategy.overflow strategy@.
The function @Function.reserve@ can be used to change the capacity explicitely.
..include:seqan/sequence.h
*/
template <typename T>
inline SEQAN_HOST_DEVICE typename Size<T const>::Type
capacity(T const & me)
{
    SEQAN_CHECKPOINT;
    return length(me);
}

// --------------------------------------------------------------------------
// Function empty()
// --------------------------------------------------------------------------

/*!
 * @fn empty
 * @headerfile <seqan/sequence.h>
 * @brief Test a container for being empty.
 *
 * @signature bool empty(object);
 *
 * @param[in] object The container to query for being empty.
 *
 * @return bool <tt>true</tt> or <tt>false</tt>, depending on whether <tt>object</tt> is empty or not.
 */

/**
.Function.empty:
..cat:Containers
..concept:Concept.ContainerConcept
..summary:Test a container for being empty.
..signature:bool empty(object)
..param.object:A container.
..returns:$true$ if $object$ contains no elements, otherwise $false$.
..remarks.text:$empty(x)$ is guaranteed to be at least as fast as $length(me) == 0$,
but can be significantly faster in some cases.
..see:Function.length
..include:seqan/sequence.h
*/
template <typename T>
SEQAN_HOST_DEVICE inline bool
empty(T const & me)
{
    SEQAN_CHECKPOINT;
    return (length(me) == 0);
}

// --------------------------------------------------------------------------
// Function _computeSizeForCapacity()
// --------------------------------------------------------------------------

// note: for value types of size 1 or 2,
// an extra position for the termination character is allocated.
// This speeds up a conversion to a c style string (see Spec.CStyle String)
// note that this extra position is necessary not only for char and wchar_t,
// but also for all other value types of size 1 and 2 to make the application
// of the funciton move for in-place alphabet conversion.


template <typename T, typename TSize>
inline TSize
_computeSizeForCapacity(T const & /*me*/,
                      TSize capacity)
{
    SEQAN_CHECKPOINT;
    typedef typename Value<T>::Type TValue;
    if (sizeof(TValue) <= 2) return capacity + 1;
    else return capacity;
}

// --------------------------------------------------------------------------
// Function computeGenerousCapacity()
// --------------------------------------------------------------------------

/*!
 * @fn computeGenerousCapacity
 * @headerfile <seqan/sequence.h>
 * @brief Capacity for generous expansion.
 *
 * @signature TSize computeGenerousCapacity(container, capacity);
 *
 * @param[in,out] container The container to compute the generous capacity for.
 * @param[in      capacity  The minimal required capacity.
 *
 * @return TSize A value larger than <tt>capacity</tt> that should be used as the new capacity for <tt>container</tt>
 *               when it is expanded using the <tt>Generous</tt> overflow strategy.
 */

// TODO(holtgrew): This is a helper and should conceptually not be in the "interface" header.
/**
.Function.computeGenerousCapacity:
..hidefromindex
..cat:Containers
..concept:Concept.ContainerConcept
..summary:Capacity for generous expansion.
..signature:Size computeGenerousCapacity(container, capacity)
..param.container:A container that should be expanded.
..param.capacity:Minimal capacity needed.
..returns:A value larger than $capacity$ that should be used as new capacity for $container$
when it is expanded using the @Tag.Overflow Strategy."Generous" overflow strategy@.
...metafunction:Metafunction.Size
..see:Tag.Overflow Strategy
..include:seqan/sequence.h
*/
template <typename T, typename TSize>
inline TSize
computeGenerousCapacity(T const & /*me*/,
                         TSize capacity)
{
    SEQAN_CHECKPOINT;
    if (capacity < 32) return 32;       // returned value is implicitly >= capacity + 1
    return capacity + (capacity >> 1);
}

// --------------------------------------------------------------------------
// Function _storageUpdated()
// --------------------------------------------------------------------------

/*
template <typename T>
inline void
_storageUpdated(T & me,
                void const *)
{
}

template <typename T>
inline void
_storageUpdated(T & me)
{
    _storageUpdated_(me, (T *) 0);
}

template <typename T>
inline void
_storageUpdated(T const & me)
{
    _storageUpdated_(me, (T const *) 0);
}
*/

// --------------------------------------------------------------------------
// Function assign()
// --------------------------------------------------------------------------

template<typename TTarget, typename TSource>
inline void
assign(TTarget & target,
       TSource & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
assign(TTarget const & target,
       TSource & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TSource>
inline void
assign(TTarget & target,
       TSource const & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
assign(TTarget const & target,
       TSource const & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

// --------------------------------------------------------------------------
// Function append()
// --------------------------------------------------------------------------

/*!
 * @fn append
 * @headerfile <seqan/sequence.h>
 * @brief Concatenate a container to another.
 *
 * @signature void append(target, source);
 *
 * @param[in,out] target The container to append <tt>source</tt> to.
 * @param[in]     source This container will be appended to <tt>source</tt>.
 */

/**
.Function.append
..summary:Concatenate two containers.
..description:This function allows to append the contents of one container to another.
...notes:Appending a value/entry to a container can be done using @Function.appendValue@.
..cat:Content Manipulation
..concept:Concept.ContainerConcept
..signature:append(target, source [, limit] [,resize_tag])
..param.target: A container $source$ is append to.
..param.source: A container that is append to $target$.
...remarks:The function does not modify this container.
..param.limit: The maximal length of $target$ after the operation. (optional)
..param.resize_tag: Specifies the strategy that is applied if $target$ has not enough capacity to store the complete content. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowImplicit@ of the $target$ type.
..remarks:The result of this operation is stored in $target$.
..see:Function.assign
..include:seqan/sequence.h
*/

template<typename TTarget, typename TSource>
inline void
append(TTarget & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    append(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget const & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    append(target, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    append(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget const & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    append(target, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget & target,
       TSource & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget const & target,
       TSource & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget & target,
       TSource const & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget const & target,
       TSource const & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

// --------------------------------------------------------------------------
// Function appendValue()
// --------------------------------------------------------------------------

/*!
 * @fn appendValue
 * @headerfile <seqan/sequence.h>
 * @brief Append a value to a container.
 *
 * @signature void appendValue(target, val[, tag]);
 *
 * @param[in,out] target The container to append <tt>val</tt> to.
 * @param[in]     val    The value to append to <tt>target</tt>.
 * @param[in]     tag    The resize tag to use.  Defaults to What DefaultOverflowImplicit returns for the type of
 *                       <tt>target</tt>.
 */

/**
.Function.appendValue:
..signature:appendValue(target, value [, resize_tag])
..description:This function allows to append an element/value to a container.
...notes:Appending the elements of a container to another container can be done using @Function.append@.
..cat:Content Manipulation
..concept:Concept.ContainerConcept
..summary:Appends a value to a container.
..param.target:A container.
..param.value:Value that is appended to $target$.
..param.resize_tag:
..param.resize_tag: Specifies the strategy that is applied if $target$ has not enough capacity to store the complete content. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowImplicit@ of the $target$ type.
..include:seqan/sequence.h
*/

template <typename T, typename TValue>
inline void
appendValue(T & me,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    appendValue(me, _value, typename DefaultOverflowImplicit<T>::Type());
}

template <typename T, typename TValue>
inline void
appendValue(T const & me,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    appendValue(me, _value, typename DefaultOverflowImplicit<T const>::Type());
}

// --------------------------------------------------------------------------
// Function insert()
// --------------------------------------------------------------------------

/*!
 * @fn insert
 * @headerfile <seqan/sequence.h>
 * @brief Inserts a sequence into a container.
 *
 * @signature void insert(target, pos, insertSeq[, tag]);
 *
 * @param[in,out] target    The container to insert elements into.
 * @param[in]     pos       The position to start inserting at.
 * @param[in]     insertSeq The sequence to insert at start.
 * @param[[in]    tag       The resize tag, defaults to what <tt>OverflowStrategyImplicit</tt> returns.
 */

/**
.Function.insert:
..cat:Content Manipulation
..summary:Inserts a sequence into a container.
..signature:insert(target, pos, insertSeq [, resize_tag])
..concept:Concept.ContainerConcept
..param.target:The container
..param.pos:Position within $target$ at which $insertSeq$ is to be inserted.
..param.insertSeq:Sequence that will be inserted into $target$.
..param.resize_tag:Strategy that is applied if $target$ has not enough capacity to store the complete content.
...type:Tag.Overflow Strategy
..see:Function.insertValue
..see:Function.append
..include:seqan/sequence.h
*/

template <typename T, typename TPosition, typename TSeq, typename TExpand>
inline void
insert(T & me,
       TPosition pos,
       TSeq const & insertSeq,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    replace(me, pos, pos, insertSeq, Tag<TExpand>());
}

template <typename T, typename TPosition, typename TSeq, typename TExpand>
inline void
insert(T const & me,
       TPosition pos,
       TSeq const & insertSeq,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    replace(me, pos, pos, insertSeq, Tag<TExpand>());
}

template <typename T, typename TPosition, typename TSeq>
inline void
insert(T & me,
       TPosition pos,
       TSeq const & insertSeq)
{
    SEQAN_CHECKPOINT;
    replace(me, pos, pos, insertSeq, typename DefaultOverflowImplicit<T>::Type());
}

template <typename T, typename TPosition, typename TSeq>
inline void
insert(T const & me,
       TPosition pos,
       TSeq const & insertSeq)
{
    SEQAN_CHECKPOINT;
    replace(me, pos, pos, insertSeq, typename DefaultOverflowImplicit<T const>::Type());
}

// --------------------------------------------------------------------------
// Function insertValue()
// --------------------------------------------------------------------------

/*!
 * @fn insertValue
 * @headerfile <seqan/sequence.h>
 * @brief Inserts an element into a container.
 *
 * @signature void insertValue(target, pos, val[, tag]);
 *
 * @param[in,out] target    The container to insert element into.
 * @param[in]     pos       The position to insert at.
 * @param[in]     val       The value to insert at start.
 * @param[[in]    tag       The resize tag, defaults to what <tt>OverflowStrategyImplicit</tt> returns.
 */

/**
.Function.insertValue:
..cat:Content Manipulation
..summary:Inserts a single value into a container.
..signature:insertValue(target, pos, value [, resize_tag])
..concept:Concept.ContainerConcept
..param.target:The container
..param.pos:Position within $target$ at which $value$ is to be inserted.
..param.value:Value that will be inserted into $target$.
..param.resize_tag:Strategy that is applied if $target$ has not enough capacity to store the complete content.
...type:Tag.Overflow Strategy
..see:Function.insert
..see:Function.appendValue
..include:seqan/sequence.h
*/

template <typename T, typename TPosition, typename TValue>
inline void
insertValue(T & me,
            TPosition pos,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    insertValue(me, pos, _value, typename DefaultOverflowImplicit<T>::Type());
}

template <typename T, typename TPosition, typename TValue>
inline void
insertValue(T const & me,
            TPosition pos,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    insertValue(me, pos, _value, typename DefaultOverflowImplicit<T const>::Type());
}

// --------------------------------------------------------------------------
// Function replace()
// --------------------------------------------------------------------------

/*!
 * @fn replace
 * @headerfile <seqan/sequence.h>
 * @brief Replaces a part of a container with another container.
 *
 * @signature void replace(target, posBegin, posEnd, source[, limit][, resizeTag]);
 *
 * @param[in,out] target    The container to modify.
 * @param[in]     posBegin  The begin position of the range to replace.
 * @param[in]     posEnd    The end position of the range to replace.
 * @param[in]     source    The source sequence to replace <tt>[posBegin, posEnd)</tt> with.
 * @param[in]     limit     Largest size of <tt>target</tt> after the operation.
 * @param[in]     resizeTag Specify the resizing behaviour.  Defaults to what <tt>DefaultOverflowImplicit</tt>
 *                          returns.
 */

/**
.Function.replace:
..summary:Replaces a part of a container with another container.
..cat:Content Manipulation
..concept:Concept.ContainerConcept
..signature:replace(target, pos_begin, pos_end, source [, limit] [,resize_tag])
..param.target: A container that is modified.
..param.pos_begin: Begin of replaced area.
...text:The first position in $target$ of the area that is replaced by $source$.
..param.pos_end: End of replaced area.
...text:The position behind the last position in $target$ of the area that is replaced by $source$.
..param.source: A container that is inserted into $target$.
...remarks:The function does not modify this container.
..param.limit: The maximal length of $target$ after the operation. (optional)
..param.resize_tag: Specifies the strategy that is applied if $target$ has not enough capacity to store the complete content. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowImplicit@ of the $target$ type.
..see:Function.assign
..see:Function.append
..remarks.text:Some compilers have difficulties if $pos_begin$ and $pos_end$ are both 0, since 0 can be
both a position or an iterator. The workaround is to convert at least one of these arguments
explicite to the position or to the interator type.
..include:seqan/sequence.h
*/
template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource & source)
{
    replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget const & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource & source)
{
    replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source)
{
    replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget const & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source)
{
    replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource & source,
        typename Size<TTarget>::Type limit)
{
    replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget const & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource & source,
        typename Size<TTarget>::Type limit)
{
    replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source,
        typename Size<TTarget>::Type limit)
{
    replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget const & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source,
        typename Size<TTarget>::Type limit)
{
    replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

// --------------------------------------------------------------------------
// Function _capacityReturned()
// --------------------------------------------------------------------------

// TODO(holtgrew): Is this internal or a helper?

template <typename T, typename TSize, typename TExpand>
inline typename Size<T>::Type
_capacityReturned(T & me,
                  TSize,
                  Tag<TExpand>)
{
    return capacity(me);
}

template <typename T, typename TSize>
inline typename Size<T>::Type
_capacityReturned(T &,
                  TSize new_capacity,
                  Insist const & )
{
    return new_capacity;
}

// --------------------------------------------------------------------------
// Function reserve()
// --------------------------------------------------------------------------

/*!
 * @fn reserve
 * @headerfile <seqan/sequence.h>
 * @fn String#reserve
 * @brief Increases the capacity.
 * 
 * @signature TSize reserve(object, new_capacity[, tag]);
 * 
 * @param[in,out] object      A container.
 * @param[in]     newCapacity The new capacity <tt>object</tt> will get.
 * @param[in]     tag         Specifies the strategy that is applied for changing the capacity.
 * 
 * @return TSize The amount of the requested capacity that was available.  That is the function returns the minimum of
 *               <tt>newCapacity</tt> and <tt>capacity(me)</tt>.
 * 
 * This function allows to increase the capacity but not the length of a container.
 * 
 * Use @link resize @endlink if you want to change the size of a container.
 * 
 * @section Remarks
 * 
 * For std::basic_string STL Adaptions, <tt>reserve</tt> is only guaranteed to have the specified behaviour with
 * <tt>Insist</tt> and <tt>Generous</tt>.
 * 
 * For std::vector STL Adaptions, <tt>reserve</tt> is only guaranteed to have the specified behaviour with
 * <tt>Insist</tt> and <tt>Generous</tt>.
 * 
 * At the end of the operation, <tt>capacity(me)</tt> can be larger than <tt>new_capacity</tt>.  If
 * <tt>new_capacity</tt> is smaller than <tt>capacity(me)</tt> at the beginning of the operation, the operation need not
 * to change the capacity at all.
 * 
 * This operation does not changes the content of <tt>object</tt>.
 * 
 * This operation may invalidate iterators of <tt>object</tt>.
 */

/**
.Function.reserve:
..cat:Containers
..summary:Increases the capacity.
..description:This function allows to increase the capacity but not the length of a container.
...note:Use @Function.resize@ if you want to change the size of a container.
..concept:Concept.ContainerConcept
..signature:Size reserve(object, new_capacity [, resize_tag])
..param.object: A container.
..param.new_capacity: The new capacity $object$ will get.
..param.resize_tag: Specifies the strategy that is applied for changing the capacity. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowExplicit@.
..returns:The amount of the requested capacity that was available.
That is the function returns the minimum of $new_capacity$ and $capacity(me)$.
...metafunction:Metafunction.Size
..remarks:At the end of the operation, $capacity(me)$ can be larger than $new_capacity$.
If $new_capacity$ is smaller than $capacity(me)$ at the beginning of the operation,
the operation need not to change the capacity at all.
..remarks:This operation does not changes the content of $object$.
...note:This operation may invalidate iterators of $object$.
..see:Function.capacity
..include:seqan/sequence.h
*/

template <typename T, typename TSize, typename TExpand>
inline typename Size<T>::Type
reserve(T & me,
        TSize const & new_capacity,
        Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    return _capacityReturned(me, new_capacity, tag);
}

template <typename T, typename TSize>
inline typename Size<T>::Type
reserve(T & me,
        TSize const & new_capacity)
{
    SEQAN_CHECKPOINT;
    return reserve(me, new_capacity, typename DefaultOverflowExplicit<T>::Type());
}

// --------------------------------------------------------------------------
// Function resize()
// --------------------------------------------------------------------------

/*!
 * @fn resize
 * @headerfile <seqan/sequence.h>
 * @brief Resizes a container.
 *
 * If the new length exceeds the old length then the new elements are filled with copies of <tt>x</tt>.
 *
 * @signature TSize resize(object, newLength[, x[, tag]]);
 *
 * @param[in,out] object  The container to resize.
 * @param[in]     newSize The new container size.
 * @param[in]     x       The value to use as the prototype when increasing the size.
 * @param[in]     tag     The strategy to apply if the capacity of <tt>object</tt> is less than <tt>newLength</tt..
 *
 * @return TSize The new length <tt>length(object)</tt>.
 */

/**
.Function.resize:
..cat:Containers
..concept:Concept.ContainerConcept
..summary:Resizes a container. If the new length exceeds the old length the new elements are filled with copies of $value$.
..signature:Size resize(object, newLength [,value [,resizeTag]])
..param.object: A container.
...type:Class.String
..param.newLength: The new length $object$ will get.
..param.value: Value that is copied if new items are created in $object$.
...remarks:If the $value$ argument is omitted, the items are not initialized if @Metafunction.IsSimple@ returns `False`.
..param.resizeTag: Specifies the strategy that is applied if the capacity of $object$ is less than $newLength$. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowExplicit@.
..returns:The new length $length(object)$.
...metafunction:Metafunction.Size
..remarks:This function can be used both for expanding and for shrinking $object$.
..see:Function.length
..see:Function.reserve
..include:seqan/sequence.h
*/

template <typename T, typename TSize>
inline typename Size<T>::Type
resize(T & me,
       TSize new_length)
{
    SEQAN_CHECKPOINT;
    return resize(me, new_length, typename DefaultOverflowExplicit<T>::Type());
}

template <typename T, typename TSize, typename TValue>
inline typename Size<T>::Type
resize(T & me,
       TSize new_length,
       TValue const & val)
{
    SEQAN_CHECKPOINT;
    return resize(me, new_length, val, typename DefaultOverflowExplicit<T>::Type());
}

// --------------------------------------------------------------------------
// Function resizeSpace()
// --------------------------------------------------------------------------

// TODO(holtgrew): Deprecated!

/*!
 * @fn String#resizeSpace
 * @headerfile <seqan/sequence.h>
 * @brief Makes free space in container
 * 
 * @signature TSize resizeSpace(object, size, posBegin, posEnd [, limit][, resizeTag]);
 * 
 * @param[in,out] object    The container. Types: String
 * @param[in]     size      Number of characters that should be freed.
 * @param[in]     posEnd    Position behind the last item in <tt>object</tt> that is to be destroyed.  If
 *                          <tt>posEnd == posBegin</tt>, no item in <tt>object</tt> will be destroyed.
 * @param[in]     posBegin  Position of the first item in <tt>object</tt> that is to be destroyed.
 * @param[in]     limit     Maximal length <tt>object</tt> can get after this operation. (optional)
 * @param[in]     resizeTag Strategy that is applied if <tt>object</tt> has not enough capacity to store the
 *                          complete content. (optional)
 * 
 * @return TSize The number of free characters.Depeding on resizeTag, this could be <tt>size</tt> or less than
 *               <tt>size</tt> if <tt>object</tt> has not enough <tt>capacity</tt>.
 */

/**
.Function.resizeSpace:
..cat:Containers
..concept:Concept.ContainerConcept
..summary:Makes free space in container
..signature:Size resizeSpace(object, size, pos_begin, pos_end [, limit] [, resize_tag])
..param.object:The container.
...type:Class.String
..param.size:Number of characters that should be freed.
..param.pos_begin:Position of the first item in $object$ that is to be destroyed.
..param.pos_end:Position behind the last item in $object$ that is to be destroyed.
...remarks:If $pos_end == pos_begin$, no item in $object$ will be destroyed.
..param.limit:Maximal length $object$ can get after this operation. (optional)
..param.resize_tag:Strategy that is applied if $object$ has not enough capacity to store the complete content. (optional)
...metafunction:Metafunction.DefaultOverflowExplicit
..returns:The number of free characters.
...metafunction:Metafunction.Size
...remarks:Depeding on the @Tag.Overflow Strategy.overflow strategy@ specified by $resize_tag$,
this could be $size$ or less than $size$ if $object$ has not enough @Function.capacity@.
..include:seqan/sequence.h
*/

template<typename T, typename TSize, typename TBeginPosition, typename TEndPosition>
inline TSize
resizeSpace(T & me,
            TSize size,
            TBeginPosition pos_begin,
            TEndPosition pos_end)
{
    SEQAN_CHECKPOINT;
    return resizeSpace(me, size, pos_begin, pos_end, typename DefaultOverflowExplicit<T>::Type());
}

template<typename T, typename TSize, typename TBeginPosition, typename TEndPosition, typename TLimit>
inline TSize
resizeSpace(T & me,
            TSize size,
            TBeginPosition pos_begin,
            TEndPosition pos_end,
            TLimit limit)
{
    SEQAN_CHECKPOINT;
    return resizeSpace(me, size, pos_begin, pos_end, limit, typename DefaultOverflowExplicit<T>::Type());
}

// --------------------------------------------------------------------------
// Function erase()
// --------------------------------------------------------------------------

/*!
 * @fn erase
 * @headerfile <seqan/sequence.h>
 * @brief Erases a part of a container.
 *
 * @signature void erase(object, pos[, posEnd]);
 *
 * @param[in,out] object The container to erase elements of.
 * @param[in]     pos    The first position to erase.
 * @param[in]     posEnd The last position to erase, defaults to <tt>pos + 1</tt>.
 */

/**
.Function.erase:
..summary:Erases a part of a container
..cat:Containers
..concept:Concept.ContainerConcept
..signature:erase(object, pos [, pos_end])
..param.object:The container.
...type:Class.String
..param.pos:Position of the first item in $object$ that is to be destroyed.
..param.pos_end:Position behind the last item in $object$ that is to be destroyed. (optional)
...default:$pos + 1$
...remarks:If $pos_end$ is omitted, only one element in $object$ at position $pos$ is destroyed.
..remarks:$erase(object, pos, pos_end)$ is semantically the same as @Function.resizeSpace.resizeSpace(object, 0, pos, pos_end)@.
..see:Function.eraseBack
..include:seqan/sequence.h
*/

template<typename T, typename TBeginPosition, typename TEndPosition>
inline void
erase(T & me,
      TBeginPosition pos,
      TEndPosition pos_end)
{
    SEQAN_CHECKPOINT;
    resizeSpace(me, 0, pos, pos_end);
}

template<typename T, typename TPosition>
inline void
erase(T & me,
      TPosition pos)
{
    SEQAN_CHECKPOINT;
    resizeSpace(me, 0, pos, pos + 1);
}

// For segments, we also have to define the version for const-containers.

template<typename T, typename TBeginPosition, typename TEndPosition>
inline void
erase(T const & me,
      TBeginPosition pos,
      TEndPosition pos_end)
{
    SEQAN_CHECKPOINT;
    resizeSpace(me, 0, pos, pos_end);
}

template<typename T, typename TPosition>
inline void
erase(T const & me,
      TPosition pos)
{
    SEQAN_CHECKPOINT;
    resizeSpace(me, 0, pos, pos + 1);
}

// --------------------------------------------------------------------------
// Function eraseBack()
// --------------------------------------------------------------------------

/*!
 * @fn eraseBack
 * @headerfile <seqan/sequence.h>
 * @brief Deletes the last item of a container and reduces its size by 1.  The container must have a size greater than
 *        or equal to 1.
 *
 * @signature void eraseBack(object);
 *
 * @param[in,out] object The container to modify.
 *
 * @section Remarks
 *
 * This is equivalent to <tt>erase(object, length(object) - 1)</tt>.
 */

/**
.Function.eraseBack:
..summary:Deletes the last item of a container and reduces its size by 1.  The container must have a size greater than or equal to 1.
..cat:Containers
..concept:Concept.ContainerConcept
..signature:eraseBack(object)
..param.object:The container.
...type:Class.String
..remarks:$erase(object)$ is semantically the same as @Function.erase.erase(me, length(me) - 1)@.
..see:Function.erase
..include:seqan/sequence.h
*/

template <typename T>
inline void eraseBack(T & me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GT_MSG(length(me), 0u, "String must have more than 0 characters in eraseBack()!");
	resize(me, length(me) - 1);
}

// --------------------------------------------------------------------------
// Function shrinkToFit()
// --------------------------------------------------------------------------

/*!
 * @fn ContainerConcept#shrinkToFit
 * @headerfile <seqan/sequence.h>
 * @brief Resizes container to minimum capacity.
 *
 * @signature void shrinkToFit(cont);
 *
 * @param[in] cont The container to shrink.
 */

/**
.Function.shrinkToFit:
..cat:Containers
..concept:Concept.ContainerConcept
..summary:Resizes container to minimum capacity
..signature:shrinkToFit(object)
..param.object: A container.
...type:Concept.ContainerConcept
..remarks
...text:$shrinkToFit(object)$ is equivalent to $reserve(object, length(object), Exact())$.
..see:Function.capacity
..see:Function.length
..see:Function.reserve
..include:seqan/sequence.h
*/

template <typename T>
inline void
shrinkToFit(T & me)
{
    SEQAN_CHECKPOINT;

//  following line has no effect as in SeqAn it is not yet possible
//  to reduce the memory consumption of a string with resize/reserve
//
//  reserve(me, length(me), Exact());

    T tmp;
    assign(tmp, me, Exact());
    swap(me, tmp);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_SEQUENCE_INTERFACE_H_
