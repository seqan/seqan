// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Manual forwards for the sequence module.
// ==========================================================================

#ifndef SEQAN_HEADER_SEQUENCE_FORWARDS_H
#define SEQAN_HEADER_SEQUENCE_FORWARDS_H

#if !defined(_MSC_VER)

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////

namespace seqan
{

// ==========================================================================
// Adaption Forwards
// ==========================================================================

// TODO(holtgrew): I wonder whether everything below will still be necessary after auto-sequence feature removal? See note below.

// NOTE(holtgrew): My guess / understanding why we need forwards here.
//
// The problem with needing forwards here appears to be that we have default
// implementations of metafunctions Reference<>, Length<> etc.
//
// Consider the setup for a template function A() using getValue(): A tries to
// use getValue().  If there was no default implementation of Reference<> then
// the function getValue() would not be instantiated at this point.  Since
// there is such a default implementation, however, getValue() gets instantiated,
// returns a Reference<std::string> == (std::string &) and this is where the
// compiler balks.
//
// I think that instantiation would get deferred in the case of Reference<> not
// being defined, but I am not sure.  I need to do more research about this.

// --------------------------------------------------------------------------
// Forwards from sequence_interface.h.
// --------------------------------------------------------------------------

template <typename T> struct AllowsFastRandomAccess;
template <typename T> struct DefaultOverflowExplicit;
template <typename T> struct DefaultOverflowImplicit;
// template <typename T, class Enable = void> struct IsContiguous;
template <typename T> struct IsSequence;
struct TagExact_;
struct TagGenerous_;
struct TagInsist_;
struct TagLimit_;
typedef Tag<TagExact_> Exact;
typedef Tag<TagGenerous_> Generous;
typedef Tag<TagInsist_> Insist;
typedef Tag<TagLimit_> Limit;
typedef Tag<TagInsist_> Tight;
template <typename T> inline typename Iterator<T, Standard>::Type _beginDefault(T & me, Standard);
template <typename T> inline typename Iterator<T const, Standard>::Type _beginDefault(T const & me, Standard);
template <typename T> inline typename Iterator<T, Rooted>::Type _beginDefault(T & me, Rooted);
template <typename T> inline typename Iterator<T const, Rooted>::Type _beginDefault(T const & me, Rooted);
template <typename T, typename TSize> inline TSize _computeSizeForCapacity(T const & , TSize capacity);
template <typename T> inline typename Iterator<T, Standard>::Type _endDefault(T & me, Standard);
template <typename T> inline typename Iterator<T const, Standard>::Type _endDefault(T const & me, Standard);
template <typename T> inline typename Iterator<T, Rooted>::Type _endDefault(T & me, Rooted);
template <typename T> inline typename Iterator<T const, Rooted>::Type _endDefault(T const & me, Rooted);
template <typename TContainer> inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, void) assign(TContainer && me, typename RemoveReference<TContainer>::Type source);
template <typename TContainer, typename TSource> inline SEQAN_FUNC_ENABLE_IF(And<Not<IsSameType<typename RemoveReference<TContainer>::Type, TSource> >, Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> > >, void) assign(TContainer && me, TSource const & source);
template <typename TContainer, typename TSource> inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, void) assign(TContainer && me, TSource const & source, typename Size<TSource>::Type limit);
template<typename TTarget, typename TSource> inline SEQAN_FUNC_DISABLE_IF(Is<StlContainerConcept<typename RemoveReference<TTarget>::Type> >, void) assign(TTarget && target, TSource && source, typename Size<TTarget>::Type const limit);
template <typename T> inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type begin(T & me);
template <typename T> inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type begin(T const & me);
template <typename T, typename TSpec> inline SEQAN_FUNC_DISABLE_IF(Is<StlContainerConcept<typename RemoveReference<T>::Type> >, typename Iterator<T, Tag<TSpec> const>::Type) begin(T & me, Tag<TSpec> const tag_);
template <typename T, typename TSpec> inline SEQAN_FUNC_DISABLE_IF(Is<StlContainerConcept<typename RemoveReference<T>::Type> >, typename Iterator<T const, Tag<TSpec> const>::Type) begin(T const & me, Tag<TSpec> const tag_);
template <typename T> inline typename Position<T>::Type beginPosition(T &);
template <typename T> inline typename Position<T>::Type beginPosition(T const &);
template <typename T, typename TSize> inline TSize computeGenerousCapacity(T const & , TSize capacity);
template <typename T> inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type end(T & me);
template <typename T> inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type end(T const & me);
template <typename T, typename TSpec> inline typename Iterator<T, Tag<TSpec> const>::Type end(T & me, Tag<TSpec> const tag_);
template <typename T, typename TSpec> inline typename Iterator<T const, Tag<TSpec> const>::Type end(T const & me, Tag<TSpec> const tag_);
template <typename T> inline typename Position<T>::Type endPosition(T & me);
template <typename T> inline typename Position<T>::Type endPosition(T const & me);
template <typename T> inline SEQAN_FUNC_DISABLE_IF(Is<StlContainerConcept<typename RemoveReference<T>::Type> >, void const *) getObjectId(T const & me);
template <typename TContainer> inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<typename RemoveReference<TContainer>::Type> >, void const *) getObjectId(TContainer && me);
template <typename T, typename TPos> inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type iter(T & me, TPos const pos);
template <typename T, typename TPos> inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type iter(T const & me, TPos const pos);
template <typename T, typename TPos, typename TTag> inline typename Iterator<T, Tag<TTag> const>::Type iter(T & me, TPos const pos, Tag<TTag> const &);
template <typename T, typename TPos, typename TTag> inline typename Iterator<T const, Tag<TTag> const>::Type iter(T const & me, TPos const pos, Tag<TTag> const &);
template <typename T> inline SEQAN_FUNC_DISABLE_IF(Is<StlContainerConcept<typename RemoveReference<T>::Type> >, typename Size<T>::Type) length(T const & /*me*/);
template <typename TChar, typename TAlloc> inline typename Size<std::forward_list<TChar, TAlloc> >::Type length(std::forward_list<TChar, TAlloc> const & me);
template <typename T, typename TValue, typename TPos> inline void moveValue(T & me, TPos pos, TValue const & _value);
template <typename T, typename TValue, typename TPos> inline void moveValue(T const & me, TPos pos, TValue const & _value);
template <typename TContainer> inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, typename Size<TContainer>::Type) length(TContainer const & me);
template <typename T, typename TSize, typename TBeginPosition, typename TEndPosition> inline TSize resizeSpace(T & me, TSize size, TBeginPosition pos_begin, TEndPosition pos_end);
template <typename T, typename TSize, typename TBeginPosition, typename TEndPosition, typename TLimit> inline TSize resizeSpace(T & me, TSize size, TBeginPosition pos_begin, TEndPosition pos_end, TLimit limit);
template <typename T1, typename T2> inline bool shareResources(T1 const & obj1, T2 const & obj2);
template <typename T> inline void shrinkToFit(T & me);
template <typename T, typename TPos> inline SEQAN_FUNC_DISABLE_IF(Is<StlContainerConcept<typename RemoveReference<T>::Type> >, typename Reference<T>::Type) value(T & me, TPos /*pos*/);
template <typename T, typename TPos> inline SEQAN_FUNC_DISABLE_IF(Is<StlContainerConcept<typename RemoveReference<T>::Type> >, typename Reference<T const>::Type) value(T const & me, TPos /*pos*/);
template <typename TContainer, typename TPos> inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, HasSubscriptOperator<TContainer> >, typename Reference<TContainer>::Type) value(TContainer & me, TPos const pos);
template <typename TContainer, typename TPos> inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, HasSubscriptOperator<TContainer> >, typename Reference<TContainer const>::Type) value(TContainer const & me, TPos const pos);
template <typename TContainer, typename TPos> inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, HasSubscriptOperator<TContainer> >, typename Value<TContainer>::Type) value(TContainer && me, TPos const pos);
template <typename TContainer, typename TPos> inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, Not<HasSubscriptOperator<TContainer> > >, typename Reference<TContainer>::Type) value(TContainer & me, TPos const pos);
template <typename TContainer, typename TPos> inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, Not<HasSubscriptOperator<TContainer> > >, typename Reference<TContainer const>::Type) value(TContainer const & me, TPos const pos);
template <typename TContainer, typename TPos> inline SEQAN_FUNC_ENABLE_IF(And<Is<StlContainerConcept<TContainer> >, Not<HasSubscriptOperator<TContainer> > >, typename Value<TContainer>::Type) value(TContainer && me, TPos const pos);

#if !(defined(STDLIB_VS) || __cplusplus > 201402L)
template <typename TContainer> inline SEQAN_FUNC_ENABLE_IF(Is<StlContainerConcept<TContainer> >, bool) empty(TContainer const & me);
template <typename T> inline SEQAN_FUNC_DISABLE_IF(Is<StlContainerConcept<typename RemoveReference<T>::Type> >, bool) empty(T const & me);
#endif

// --------------------------------------------------------------------------
// Forwards For arrays and pointers.
// --------------------------------------------------------------------------

template <typename TValue> struct DefaultOverflowExplicit;
template <typename TValue> struct DefaultOverflowImplicit;
template <typename TValue, typename TExpand> inline size_t _clearSpace(TValue * me, size_t size, Tag<TExpand>);
template <typename TValue, typename TExpand> inline size_t _clearSpace(TValue * me, size_t size, size_t limit, Tag<TExpand>);
template <typename TValue, typename TPosition, typename TExpand> inline size_t _clearSpace(TValue * me, size_t size, TPosition pos_begin, TPosition pos_end, Tag<TExpand>);
template <typename TValue, typename TPosition, typename TExpand> inline size_t _clearSpace(TValue * me, size_t size, TPosition pos_begin, TPosition pos_end, size_t limit, Tag<TExpand>);
template <typename TValue> inline void _setLength(TValue * me, size_t new_length);
template <typename TTargetValue, typename TSource, typename TExpand> inline void append(TTargetValue * target, TSource const & source, Tag<TExpand>);
template <typename TTargetValue, typename TSource, typename TExpand> inline void append(TTargetValue * target, TSource const & source, size_t limit, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void append(TTargetValue * target, TSourceValue const * source, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void append(TTargetValue * target, TSourceValue const * source, size_t limit, Tag<TExpand>);
template <typename TTargetValue, typename TSource> inline typename EnableIf<IsCharType<TTargetValue> >::Type assign(TTargetValue * target, TSource & source);
template <typename TTargetValue, typename TSource> inline typename EnableIf<IsCharType<TTargetValue> >::Type assign(TTargetValue * target, TSource const & source);
template <typename TTargetValue, typename TSource, typename TExpand> inline void assign(TTargetValue * target, TSource const & source, Tag<TExpand>);
template <typename TTargetValue, typename TSource, typename TExpand> inline void assign(TTargetValue * target, TSource const & source, size_t limit, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void assign(TTargetValue * target, TSourceValue const * source, Tag<TExpand>);
template <typename TTargetValue, typename TSourceValue, typename TExpand> inline void assign(TTargetValue * target, TSourceValue const * source, size_t limit, Tag<TExpand>);
template <typename TValue, typename TPos> inline void assignValue(TValue * me, TPos pos, TValue const & _value);
template <typename TValue> inline bool atEnd(TValue * pos);
template <typename TValue> inline bool atEnd(TValue const * pos, TValue const * );
template <typename T> inline typename Iterator<T *, typename DefaultGetIteratorSpec<T>::Type>::Type begin(T * me);
template <typename TValue> inline typename Iterator<TValue *, Standard>::Type begin(TValue * me, Standard);
template <typename TValue> inline typename Iterator<TValue const *, Standard>::Type begin(TValue const * me, Standard);
template <typename TValue, typename TSpec> inline typename Iterator<TValue *, Tag<TSpec> const>::Type begin(TValue * me, Tag<TSpec> const);
template <typename TValue, typename TSpec> inline typename Iterator<TValue const *, Tag<TSpec> const>::Type begin(TValue const * me, Tag<TSpec> const);
template <typename TValue> inline void clear(TValue * me);
template <typename TValue> inline bool empty(TValue * me);
template <typename TValue> inline typename Iterator<TValue *, Standard>::Type end(TValue * me, Standard);
template <typename TValue> inline typename Iterator<TValue const *, Standard>::Type end(TValue const * me, Standard);
template <typename TValue, typename TSpec> inline typename Iterator<TValue *, Tag<TSpec> const>::Type end(TValue * me, Tag<TSpec> const tag_);
template <typename TValue, typename TSpec> inline typename Iterator<TValue const *, Tag<TSpec> const>::Type end(TValue const * me, Tag<TSpec> const tag_);
template <typename TLeftValue, typename TRight > inline bool isEqual(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight> inline bool isGreater(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight> inline bool isGreaterOrEqual(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight> inline bool isLess(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight> inline bool isLessOrEqual(TLeftValue * left, TRight const & right);
template <typename TLeftValue, typename TRight > inline bool isNotEqual(TLeftValue * left, TRight const & right);
template <typename TValue> inline size_t length(TValue * me);
template <typename TValue> inline size_t length(TValue const * me);
inline size_t length(char * me);
inline size_t length(char const * me);
template <typename TTargetValue, typename TSource> inline void move(TTargetValue * & target, TSource & source);
template <typename TTargetValue, typename TSource> inline void move(TTargetValue * & target, TSource const & source);
template <typename TValue, typename TPos> inline void moveValue(TValue * me, TPos pos, TValue const & _value);
template <typename TValue, typename TSize, typename TExpand> inline size_t resize( TValue * me, TSize new_length, Tag<TExpand>);
template <typename TValue, typename TSize, typename TExpand> inline size_t resize( TValue * me, TSize new_length, TValue const & val, Tag<TExpand>);
template <typename TValue, typename TPos> inline TValue & value(TValue * me, TPos pos);
template <typename TValue, typename TPos> inline TValue const & value(TValue const * me, TPos pos);

}  // namespace seqan

#endif  // #if !defined(_MSC_VER)

#endif

