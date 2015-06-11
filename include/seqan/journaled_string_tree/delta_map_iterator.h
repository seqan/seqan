// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements the iterator over the delta map.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ITERATOR_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ITERATOR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <typename TType, typename TOtherType> struct IsConstructible;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class DeltaMapIterator
 * @implements BidirectionalIteratorConcept
 * @extends Iter
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Bidirectional iterator over a @link DeltaMap @endlink.
 *
 * @signature template <typename TDeltaMap>
 *            class Iter<TDeltaMap, DeltaMapIteratorSpec>;
 *
 * @tparam TDeltaMap The delta map to iterate.
 */

template <typename TDeltaMap>
class Iter<TDeltaMap, DeltaMapIteratorSpec>
{
public:
    typedef typename RemoveConst<TDeltaMap>::Type TNonConstDeltaMap_;

    typedef typename Member<TDeltaMap, DeltaMapEntriesMember>::Type TMapEntries_;
    typedef typename Iterator<TMapEntries_, Standard>::Type TMapIterator;

    TDeltaMap*   _mapPtr;
    TMapIterator _mapIter;

    //Default C'tor.
    Iter() : _mapPtr(NULL), _mapIter()
    {}

    // Copy C'tor.
    template <typename TDeltaMapOther>
    Iter(Iter<TDeltaMapOther, DeltaMapIteratorSpec> const & other,
         SEQAN_CTOR_ENABLE_IF(IsConstructible<TDeltaMap, TDeltaMapOther>)) :
        _mapPtr(other._mapPtr),
        _mapIter(other._mapIter)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Custom C'tor.
    template <typename TPos>
    Iter(TDeltaMap & container, TPos pos) : _mapPtr(&container), _mapIter()
    {
        _mapIter = begin(_mapPtr->_entries, Standard()) + pos;
    }

    // Assignment Operator.
    template <typename TDeltaMapOther>
    SEQAN_FUNC_DISABLE_IF(IsConstructible<TDeltaMap, TDeltaMapOther>, Iter<TDeltaMap, DeltaMapIteratorSpec> &)
    operator=(Iter<TDeltaMapOther, DeltaMapIteratorSpec> const & other)
    {
        _mapPtr = other._mapPtr;
        _mapIter = other._mapIter;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
struct Difference<Iter<TDeltaMap, DeltaMapIteratorSpec> >
{
    typedef Iter<TDeltaMap, DeltaMapIteratorSpec> TIter_;
    typedef typename TIter_::TMapIterator TMapIterator_;
    typedef typename Difference<TMapIterator_>::Type Type;
};

template <typename TDeltaMap>
struct Difference<Iter<TDeltaMap, DeltaMapIteratorSpec> const > :
    Difference<Iter<TDeltaMap, DeltaMapIteratorSpec> >{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename Container<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type &
container(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    return *iter._mapPtr;
}

template <typename TDeltaMap>
inline typename Container<Iter<TDeltaMap, DeltaMapIteratorSpec> const>::Type &
container(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    return *iter._mapPtr;
}

// ----------------------------------------------------------------------------
// Function deltaValue()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMapIterator#deltaValue
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Returns the delta value associated with the current iterator position.
 *
 * @signature TDeltaValue deltaValue(it, tag);
 * @param[in]   it  The iterator to query the SNP event for.
 * @param[in]   tag The tag used to select the requested delta type.
 *
 * @return TDeltaValue A reference to the delta at the current iterator position of type @link DeltaMap#DeltaValue @endlink.
 */

template <typename TDeltaMap, typename TTag>
inline typename DeltaValue<typename Container<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type, TTag>::Type &
deltaValue(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter, TTag const & tag)
{
    return deltaValue(container(iter)._deltaStore, getDeltaRecord(value(iter)).i2, tag);
}

template <typename TDeltaMap, typename TTag>
inline typename DeltaValue<typename Container<Iter<TDeltaMap, DeltaMapIteratorSpec> const>::Type, TTag>::Type &
deltaValue(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter, TTag const & tag)
{
    return deltaValue(container(iter)._deltaStore, getDeltaRecord(value(iter)).i2, tag);
}

// ----------------------------------------------------------------------------
// Function operator*()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename Reference<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type
operator*(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    return *iter._mapIter;
}

template <typename TDeltaMap>
inline typename Reference<Iter<TDeltaMap, DeltaMapIteratorSpec> const>::Type
operator*(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    return *iter._mapIter;
}

// ----------------------------------------------------------------------------
// Function operator++()                                               [prefix]
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline Iter<TDeltaMap, DeltaMapIteratorSpec> &
operator++(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    ++iter._mapIter;
    return iter;
}

// ----------------------------------------------------------------------------
// Function operator++()                                              [postfix]
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline Iter<TDeltaMap, DeltaMapIteratorSpec>
operator++(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter,  int /*postfix*/)
{
    Iter<TDeltaMap, DeltaMapIteratorSpec> temp(iter);
    ++iter;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator+=()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSize>
inline Iter<TDeltaMap, DeltaMapIteratorSpec> &
operator+=(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter,  TSize const & len)
{
    iter._mapIter += len;
    return iter;
}

// ----------------------------------------------------------------------------
// Function operator+()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSize>
inline Iter<TDeltaMap, DeltaMapIteratorSpec>
operator+(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter,  TSize const & len)
{
    Iter<TDeltaMap, DeltaMapIteratorSpec> temp(iter);
    temp += len;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator--()                                              [postfix]
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline Iter<TDeltaMap, DeltaMapIteratorSpec> &
operator--(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    --iter._mapIter;
    return iter;
}

// ----------------------------------------------------------------------------
// Function operator--()                                               [prefix]
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline Iter<TDeltaMap, DeltaMapIteratorSpec>
operator--(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter,  int /*postfix*/)
{
    Iter<TDeltaMap, DeltaMapIteratorSpec> temp(iter);
    --iter;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSize>
inline Iter<TDeltaMap, DeltaMapIteratorSpec> &
operator-=(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter,  TSize const & len)
{
    iter._mapIter -= len;
    return iter;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSize>
inline Iter<TDeltaMap, DeltaMapIteratorSpec>
operator-(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter,  TSize const & len)
{
    Iter<TDeltaMap, DeltaMapIteratorSpec> temp(iter);
    temp -= len;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename Difference<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type
operator-(Iter<TDeltaMap, DeltaMapIteratorSpec> const & lhs,
          Iter<TDeltaMap, DeltaMapIteratorSpec> const & rhs)
{
    return lhs._mapIter - rhs._mapIter;
}

template <typename TDeltaMap>
inline typename Difference<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type
operator-(Iter<TDeltaMap, DeltaMapIteratorSpec> const & lhs,
          typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & rhs)
{
    return lhs._mapIter - rhs._mapIter;
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline bool
operator==(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           Iter<TDeltaMap, DeltaMapIteratorSpec> const & b)
{
    return a._mapIter == b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator==(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & b)
{
    return a._mapIter == b._mapIter;
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline bool
operator!=(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           Iter<TDeltaMap, DeltaMapIteratorSpec> const & b)
{
    return !(a == b);
}

template <typename TDeltaMap>
inline bool
operator!=(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & b)
{
    return !(a == b);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline bool
operator<(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           Iter<TDeltaMap, DeltaMapIteratorSpec> const & b)
{
    return a._mapIter < b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator<(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
          typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & b)
{
    return a._mapIter < b._mapIter;
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline bool
operator<=(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           Iter<TDeltaMap, DeltaMapIteratorSpec> const & b)
{
    return a._mapIter <= b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator<=(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & b)
{
    return a._mapIter <= b._mapIter;
}

// ----------------------------------------------------------------------------
// Function operator>()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline bool
operator>(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           Iter<TDeltaMap, DeltaMapIteratorSpec> const & b)
{
    return a._mapIter > b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator>(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & b)
{
    return a._mapIter > b._mapIter;
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline bool
operator>=(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           Iter<TDeltaMap, DeltaMapIteratorSpec> const & b)
{
    return a._mapIter >= b._mapIter;
}

template <typename TDeltaMap>
inline bool
operator>=(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & b)
{
    return a._mapIter >= b._mapIter;
}

//template <typename TDeltaMapLhs, typename TDeltaMapRhs>
//inline void
//swap(Iter<TDeltaMapLhs, DeltaMapIteratorSpec> & lhs,
//     Iter<TDeltaMapRhs, DeltaMapIteratorSpec> & rhs)
//{
//    swap(lhs._mapPtr, rhs._mapPtr);
//    swap(lhs._mapIter, rhs._mapIter);
//}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ITERATOR_H_
