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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements the iterator over the delta map.
// ==========================================================================
#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ITERATOR_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ITERATOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

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
    typedef typename GetMapValueString_<TDeltaMap>::Type TDeltaMapKeys;
    typedef typename Iterator<TDeltaMapKeys, Standard>::Type TMapIterator;

    typedef typename GetDeltaStore_<TDeltaMap>::Type TDeltaStore_;
    typedef typename Iterator<TDeltaStore_, Standard>::Type TDeltaStoreIter;

    typedef typename GetDeltaCoverageStore_<TDeltaMap>::Type TDeltaCoverageStore_;
    typedef typename Iterator<TDeltaCoverageStore_, Standard>::Type TCoverageIter;

    TMapIterator    _mapIter;
    TDeltaStoreIter _deltaStoreIter;
    TCoverageIter   _deltaCoverageIter;

    Iter() : _mapIter(), _deltaStoreIter(), _deltaCoverageIter()
    {}

    // Copy Constructor
    Iter(Iter<TDeltaMap, DeltaMapIteratorSpec> const & other) : _mapIter(other._mapIter),
                                                                 _deltaStoreIter(other._deltaStoreIter),
                                                                 _deltaCoverageIter(other._deltaCoverageIter)
    {}

    Iter(typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & other) :
                    _mapIter(other._mapIter),
                    _deltaStoreIter(other._deltaStoreIter),
                   _deltaCoverageIter(other._deltaCoverageIter)
    {}

    // Assignment Operator
    Iter<TDeltaMap, DeltaMapIteratorSpec> &
    operator=(Iter<TDeltaMap, DeltaMapIteratorSpec> const & other)
    {
        if (this != &other)
        {
            _mapIter = other._mapIter;
            _deltaStoreIter = other._deltaStoreIter;
            _deltaCoverageIter = other._deltaCoverageIter;
        }
        return *this;
    }

    Iter<TDeltaMap, DeltaMapIteratorSpec> &
    operator=(typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & other)
    {
        _mapIter = other._mapIter;
        _deltaStoreIter = other._deltaStoreIter;
        _deltaCoverageIter = other._deltaCoverageIter;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
struct Value<Iter<TDeltaMap, DeltaMapIteratorSpec> >
{
    typedef typename Value<TDeltaMap>::Type Type;
};

template <typename TDeltaMap>
struct Value<Iter<TDeltaMap, DeltaMapIteratorSpec> const >
{
    typedef typename Value<TDeltaMap const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
struct Reference<Iter<TDeltaMap, DeltaMapIteratorSpec> >
{
    typedef typename Reference<TDeltaMap>::Type Type;
};

template <typename TDeltaMap>
struct Reference<Iter<TDeltaMap, DeltaMapIteratorSpec> const >
{
    typedef typename Reference<TDeltaMap const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
struct GetValue<Iter<TDeltaMap, DeltaMapIteratorSpec> > :
    Value<Iter<TDeltaMap, DeltaMapIteratorSpec> >{};

template <typename TDeltaMap>
struct GetValue<Iter<TDeltaMap, DeltaMapIteratorSpec> const> :
    Value<Iter<TDeltaMap, DeltaMapIteratorSpec> const>{};

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

template <typename TDeltaMap>
inline void
_initBegin(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter,
           TDeltaMap & map)
{
    iter._mapIter = begin(map._keys, Standard());
    iter._deltaStoreIter = begin(map._deltaStore, Standard());
    iter._deltaCoverageIter = begin(map._deltaCoverageStore, Standard());
}

template <typename TDeltaMap>
inline void
_initEnd(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter,
           TDeltaMap & map)
{
    iter._mapIter = end(map._keys, Standard());
    iter._deltaStoreIter = end(map._deltaStore, Standard());
    iter._deltaCoverageIter = end(map._deltaCoverageStore, Standard());
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename Value<TDeltaMap>::Type const &
value(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    return value(iter._mapIter);
}

template <typename TDeltaMap>
inline typename Reference<TDeltaMap const>::Type
value(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    return value(iter._mapIter);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename GetValue<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const &
getValue(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    return getValue(iter._mapIter);
}

template <typename TDeltaMap>
inline typename GetValue<Iter<TDeltaMap, DeltaMapIteratorSpec> const>::Type
getValue(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    return getValue(iter._mapIter);
}

// ----------------------------------------------------------------------------
// Function deltaType()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMapIterator#deltaType
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Returns the id of the delta event the current iterator points to.
 *
 * @signature TId deltaType(it)
 * @param[in]   it  The iterator to query the delta event key for.
 *
 * @return TId The id for the current delta event of type <tt>DeltaType::TValue</tt>.
 */

template <typename TDeltaMap>
inline DeltaType::TValue
deltaType(Iter<TDeltaMap, DeltaMapIteratorSpec> const& iter)
{
    return deltaType(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaPosition()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline DeltaType::TValue
deltaPosition(Iter<TDeltaMap, DeltaMapIteratorSpec> const& iter)
{
    return deltaPosition(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaSnp()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMapIterator#deltaSnp
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Returns the SNP associated with the current iterator position.
 *
 * @signature TSnp deltaSnp(it)
 * @param[in]   it  The iterator to query the SNP event for.
 *
 * @return TSnp A reference to the SNP at the current iterator position of type @link DeltaMap#DeltaValue @endlink.
 *
 * @see DeltaMapIterator#deltaDel
 * @see DeltaMapIterator#deltaIns
 * @see DeltaMapIterator#deltaIndel
 */

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    SEQAN_ASSERT_EQ(deltaType(iter), DeltaType::DELTA_TYPE_SNP);
    return deltaSnp(iter._deltaStoreIter);
}

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    SEQAN_ASSERT_EQ(deltaType(iter), DeltaType::DELTA_TYPE_SNP);
    return deltaSnp(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaDel()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMapIterator#deltaDel
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Returns the deletion associated with the current iterator position.
 *
 * @signature TDel deltaDel(it)
 * @param[in]   it  The iterator to query the deletion event for.
 *
 * @return TDel A reference to the deletion at the current iterator position of type @link DeltaMap#DeltaValue @endlink.
 *
 * @see DeltaMapIterator#deltaSnp
 * @see DeltaMapIterator#deltaIns
 * @see DeltaMapIterator#deltaIndel
 */

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    SEQAN_ASSERT_EQ(deltaType(iter), DeltaType::DELTA_TYPE_DEL);
    return deltaDel(iter._deltaStoreIter);
}

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    SEQAN_ASSERT_EQ(deltaType(iter), DeltaType::DELTA_TYPE_DEL);
    return deltaDel(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaIns()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMapIterator#deltaIns
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Returns the insertion associated with the current iterator position.
 *
 * @signature TIns deltaIns(it)
 * @param[in]   it  The iterator to query the insertion event for.
 *
 * @return TIns A reference to the insertion at the current iterator position of type @link DeltaMap#DeltaValue @endlink.
 *
 * @see DeltaMapIterator#deltaSnp
 * @see DeltaMapIterator#deltaDel
 * @see DeltaMapIterator#deltaIndel
 */

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    SEQAN_ASSERT_EQ(deltaType(iter), DeltaType::DELTA_TYPE_INS);
    return deltaIns(iter._deltaStoreIter);
}

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    SEQAN_ASSERT_EQ(deltaType(iter), DeltaType::DELTA_TYPE_INS);
    return deltaIns(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaIndel()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMapIterator#deltaIndel
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Returns the replacement associated with the current iterator position.
 *
 * @signature TIndel deltaIndel(it)
 * @param[in]   it  The iterator to query the replacement event for.
 *
 * @return TIndel A reference to the replacement at the current iterator position of type @link DeltaMap#DeltaValue @endlink.
 *
 * @see DeltaMapIterator#deltaSnp
 * @see DeltaMapIterator#deltaDel
 * @see DeltaMapIterator#deltaIns
 */

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    SEQAN_ASSERT_EQ(deltaType(iter), DeltaType::DELTA_TYPE_INDEL);
    return deltaIndel(iter._deltaStoreIter);
}

template <typename TDeltaMap>
inline typename DeltaValue<TDeltaMap const, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    SEQAN_ASSERT_EQ(deltaType(iter), DeltaType::DELTA_TYPE_INDEL);
    return deltaIndel(iter._deltaStoreIter);
}

// ----------------------------------------------------------------------------
// Function deltaCoverage()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMapIterator#deltaCoverage
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Returns the coverage associated with the current iterator position.
 *
 * @signature TCov deltaCoverage(it)
 * @param[in]   it  The iterator to query the coverage for.
 *
 * @return TCov A reference to the coverage at the current iterator position of type @link DeltaMap#DeltaCoverage @endlink.
 */

template <typename TDeltaMap>
inline typename DeltaCoverage<TDeltaMap>::Type &
deltaCoverage(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    return value(iter._deltaCoverageIter);
}

template <typename TDeltaMap>
inline typename DeltaCoverage<TDeltaMap>::Type &
deltaCoverage(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    return value(iter._deltaCoverageIter);
}

// ----------------------------------------------------------------------------
// Function operator*()
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline typename Value<TDeltaMap>::Type const &
operator*(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    return value(iter);
}

template <typename TDeltaMap>
inline typename Value<TDeltaMap>::Type const &
operator*(Iter<TDeltaMap, DeltaMapIteratorSpec> const & iter)
{
    return value(iter);
}

// ----------------------------------------------------------------------------
// Function operator++()                                               [prefix]
// ----------------------------------------------------------------------------

template <typename TDeltaMap>
inline Iter<TDeltaMap, DeltaMapIteratorSpec> &
operator++(Iter<TDeltaMap, DeltaMapIteratorSpec> & iter)
{
    ++iter._mapIter;
    ++iter._deltaStoreIter;
    ++iter._deltaCoverageIter;
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
    iter._deltaStoreIter += len;
    iter._deltaCoverageIter += len;
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
    --iter._deltaStoreIter;
    --iter._deltaCoverageIter;
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
    iter._deltaStoreIter -= len;
    iter._deltaCoverageIter -= len;
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
    if (a._mapIter != b._mapIter)
        return false;
    return true;
}

template <typename TDeltaMap>
inline bool
operator==(Iter<TDeltaMap, DeltaMapIteratorSpec> const & a,
           typename IterComplementConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type const & b)
{
    typedef typename IterMakeConst<Iter<TDeltaMap, DeltaMapIteratorSpec> >::Type TConstIter;
    return static_cast<TConstIter>(a) == static_cast<TConstIter>(b);
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

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_ITERATOR_H_
