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
// Implements a data structure to store deltas efficiently.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_STORE_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_STORE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Struct DeltaType
// ----------------------------------------------------------------------------

/*!
 * @enum DeltaType
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Keys for specifying the delta type to be accessed.
 *
 * @val DeltaType DELTA_TYPE_SNP
 * @brief Id to denote SNP events.
 *
 * @val DeltaType DELTA_TYPE_DEL
 * @brief Id to denote deletion events.
 *
 * @val DeltaType DELTA_TYPE_INS
 * @brief Id to denote insertion events.
 *
 * @val DeltaType DElTA_TYPE_INDEL
 * @brief Id to denote replacement events.
 */

struct DeltaType
{
    typedef size_t TValue;

    static const TValue MASK_DELTA_TYPE;
    static const TValue MASK_DELTA_POSITION;

    static const TValue DELTA_TYPE_SNP;
    static const TValue DELTA_TYPE_DEL;
    static const TValue DELTA_TYPE_INS;
    static const TValue DELTA_TYPE_INDEL;
};

// We make: 00 = SNP
//          01 = DEL
//          10 = INS
//          11 = INS_DEL -> INS_SNP can be replaced by INS_DEL.
const size_t DeltaType::MASK_DELTA_TYPE = 3ull << (BitsPerValue<size_t>::VALUE - 2);
const size_t DeltaType::MASK_DELTA_POSITION = ~MASK_DELTA_TYPE;
const size_t DeltaType::DELTA_TYPE_SNP = 0ull;
const size_t DeltaType::DELTA_TYPE_DEL = 1ull << (BitsPerValue<size_t>::VALUE - 2);
const size_t DeltaType::DELTA_TYPE_INS = 2ull << (BitsPerValue<size_t>::VALUE - 2);
const size_t DeltaType::DELTA_TYPE_INDEL = 3ull << (BitsPerValue<size_t>::VALUE - 2);

template <typename T>
struct GetDataMap_;

// ----------------------------------------------------------------------------
// Class DeltaStore
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
class DeltaStore
{
public:
    typedef String<TAlphabet> TSnpData;
    typedef String<TSize>   TDelData;
    typedef String<String<TAlphabet> > TInsData;
    typedef Pair<TSize, String<TAlphabet> > TInsDelDataValue;
    typedef String<TInsDelDataValue> TInsDelData;
    typedef typename GetDataMap_<DeltaStore>::Type TMap;

    // TODO(rmaerker): Elaborate on these ideas!
    // Idea a) Use ConcatStringSet for insertions. Use as global insertion buffer for all journal sequences.
    // Idea b) Use bit encoding for DNA alphabet.
    // Idea c) Instead of insertion buffer, we append the inserted strings to the reference and only use original nodes.

    TMap _varDataMap;

    TSnpData    _snpData;
    TDelData    _delData;
    TInsData    _insData;
    TInsDelData _indelData;

    DeltaStore()
    {}
};

struct SpecDeltaStoreIterator_;

template <typename TDeltaStore>
class Iter<TDeltaStore,  SpecDeltaStoreIterator_>
{
public:
    typedef typename GetDataMap_<TDeltaStore>::Type TMap_;
    typedef typename Iterator<TMap_, Standard>::Type TMapIterator;

    TDeltaStore* _containerPtr;
    TMapIterator _dataIter;

    Iter() : _containerPtr(NULL), _dataIter()
    {}

    // Copy Constructor
    Iter(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & other) : _containerPtr(other._containerPtr),
                                                                     _dataIter(other._dataIter)
    {}

    Iter(typename IterComplementConst<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type const & other) :
                     _containerPtr(other._containerPtr),
                     _dataIter(other._dataIter)
    {}

    // Assignment Operator
    Iter<TDeltaStore, SpecDeltaStoreIterator_> &
    operator=(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & other)
    {
        if (this != &other)
        {
            _containerPtr = other._containerPtr;
            _dataIter = other._dataIter;
        }
        return *this;
    }

    Iter<TDeltaStore, SpecDeltaStoreIterator_> &
    operator=(typename IterComplementConst<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type const & other)
    {
        _containerPtr = other._containerPtr;
        _dataIter = other._dataIter;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct Value<DeltaStore<TSize, TAlphabet> >
{
    typedef typename DeltaType::TValue Type;
};

template <typename TSize, typename TAlphabet>
struct Value<DeltaStore<TSize, TAlphabet> const>
{
    typedef typename DeltaType::TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct Reference<DeltaStore<TSize, TAlphabet> >
{
    typedef typename DeltaType::TValue & Type;
};

template <typename TSize, typename TAlphabet>
struct Reference<DeltaStore<TSize, TAlphabet> const>
{
    typedef typename DeltaType::TValue const & Type;
};


// ----------------------------------------------------------------------------
// Metafunction GetDataMap_
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct GetDataMap_<DeltaStore<TSize, TAlphabet> >
{
    typedef typename DeltaType::TValue TValue_;
    typedef String<TValue_> Type;
};

template <typename TSize, typename TAlphabet>
struct GetDataMap_<DeltaStore<TSize, TAlphabet> const>
{
    typedef typename DeltaType::TValue TValue_;
    typedef String<TValue_> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
struct Difference<Iter<TDeltaStore, SpecDeltaStoreIterator_> >
{
    typedef Iter<TDeltaStore, SpecDeltaStoreIterator_> TIter_;
    typedef typename TIter_::TMapIterator TMapIterator_;
    typedef typename Difference<TMapIterator_>::Type Type;
};

template <typename TDeltaStore>
struct Difference<Iter<TDeltaStore, SpecDeltaStoreIterator_> const> :
    Difference<Iter<TDeltaStore, SpecDeltaStoreIterator_> >{};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_SNP]
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename DeltaType::TValue>
struct DeltaValue;

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_SNP>
{
    typedef TAlphabet Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_SNP>
{
    typedef TAlphabet const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_DEL]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_DEL>
{
    typedef TSize Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_DEL>
{
    typedef TSize const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                     [DELTA_TYPE_INS]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INS>
{
    typedef String<TAlphabet> Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_INS>
{
    typedef String<TAlphabet> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue                                   [DELTA_TYPE_INDEL]
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INDEL>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INS>::Type TIns_;
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_DEL>::Type TDel_;
    typedef Pair<TDel_, TIns_> Type;
};

template <typename TSize, typename TAlphabet>
struct DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_INDEL>
{
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INS>::Type TIns_;
    typedef typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_DEL>::Type TDel_;
    typedef Pair<TDel_, TIns_> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
struct Iterator<DeltaStore<TSize, TAlphabet>, Standard>
{
    typedef Iter<DeltaStore<TSize, TAlphabet>, SpecDeltaStoreIterator_> Type;
};

template <typename TAlphabet, typename TSize>
struct Iterator<DeltaStore<TSize, TAlphabet> const, Standard>
{
    typedef Iter<DeltaStore<TSize, TAlphabet> const, SpecDeltaStoreIterator_> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
inline typename Iterator<DeltaStore<TSize, TAlphabet>, Standard>::Type
begin(DeltaStore<TSize, TAlphabet> & deltaStore, Standard const /*tag*/)
{
    typedef typename Iterator<DeltaStore<TSize, TAlphabet>, Standard>::Type TIterator;
    TIterator tmp;
    tmp._dataIter = begin(deltaStore._varDataMap, Standard());
    tmp._containerPtr = &deltaStore;
    return tmp;
}

template <typename TSize, typename TAlphabet>
inline typename Iterator<DeltaStore<TSize, TAlphabet> const, Standard>::Type
begin(DeltaStore<TSize, TAlphabet> const & deltaStore, Standard const /*tag*/)
{
    typedef typename Iterator<DeltaStore<TSize, TAlphabet> const, Standard>::Type TIterator;
    TIterator tmp;
    tmp._dataIter = begin(deltaStore._varDataMap, Standard());
    tmp._containerPtr = &deltaStore;
    return tmp;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet>
inline typename Iterator<DeltaStore<TSize, TAlphabet>, Standard>::Type
end(DeltaStore<TSize, TAlphabet> & deltaStore, Standard const /*tag*/)
{
    typedef typename Iterator<DeltaStore<TSize, TAlphabet>, Standard>::Type TIterator;
    TIterator tmp;
    tmp._dataIter = end(deltaStore._varDataMap, Standard());
    tmp._containerPtr = &deltaStore;
    return tmp;
}

template <typename TSize, typename TAlphabet>
inline typename Iterator<DeltaStore<TSize, TAlphabet> const, Standard>::Type
end(DeltaStore<TSize, TAlphabet> const & deltaStore, Standard const /*tag*/)
{
    typedef typename Iterator<DeltaStore<TSize, TAlphabet> const, Standard>::Type TIterator;
    TIterator tmp;
    tmp._dataIter = end(deltaStore._varDataMap, Standard());
    tmp._containerPtr = &deltaStore;
    return tmp;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline typename Reference<TDeltaStore>::Type
value(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter)
{
    return value(iter._dataIter);
}

template <typename TDeltaStore>
inline typename Reference<TDeltaStore const>::Type
value(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & iter)
{
    return value(iter._dataIter);
}

// ----------------------------------------------------------------------------
// Function operator++()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline Iter<TDeltaStore, SpecDeltaStoreIterator_> &
operator++(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter)
{
    ++iter._dataIter;
    return iter;
}

template <typename TDeltaStore>
inline Iter<TDeltaStore, SpecDeltaStoreIterator_>
operator++(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter,  int /*postfix*/)
{
    Iter<TDeltaStore, SpecDeltaStoreIterator_> temp(iter);
    ++iter;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator+=()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TSize>
inline Iter<TDeltaStore, SpecDeltaStoreIterator_> &
operator+=(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter,  TSize const & len)
{
    iter._dataIter += len;
    return iter;
}

// ----------------------------------------------------------------------------
// Function operator+()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TSize>
inline Iter<TDeltaStore, SpecDeltaStoreIterator_>
operator+(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & iter,  TSize const & len)
{
    Iter<TDeltaStore, SpecDeltaStoreIterator_> temp(iter);
    temp += len;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator--()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline Iter<TDeltaStore, SpecDeltaStoreIterator_> &
operator--(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter)
{
    --iter._dataIter;
    return iter;
}

template <typename TDeltaStore>
inline Iter<TDeltaStore, SpecDeltaStoreIterator_>
operator--(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter,  int /*postfix*/)
{
    Iter<TDeltaStore, SpecDeltaStoreIterator_> temp(iter);
    --iter;
    return temp;
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TSize>
inline Iter<TDeltaStore, SpecDeltaStoreIterator_> &
operator-=(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter,  TSize const & len)
{
    iter._dataIter -= len;
    return iter;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TDeltaStore, typename TSize>
inline Iter<TDeltaStore, SpecDeltaStoreIterator_>
operator-(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & iter,  TSize const & len)
{
    Iter<TDeltaStore, SpecDeltaStoreIterator_> temp(iter);
    temp -= len;
    return temp;
}

template <typename TDeltaStore>
inline typename Difference<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type
operator-(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & lhs,
          Iter<TDeltaStore, SpecDeltaStoreIterator_> const & rhs)
{
    return lhs._dataIter - rhs._dataIter;
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline bool
operator==(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           Iter<TDeltaStore, SpecDeltaStoreIterator_> const & b)
{
    if (a._dataIter != b._dataIter)
        return false;
    if (a._containerPtr != b._containerPtr)
        return false;
    return true;
}

template <typename TDeltaStore>
inline bool
operator==(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type const & b)
{
    typedef typename IterMakeConst<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type TConstIter;
    return static_cast<TConstIter>(a) == static_cast<TConstIter>(b);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline bool
operator!=(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           Iter<TDeltaStore, SpecDeltaStoreIterator_> const & b)
{
    return !(a == b);
}

template <typename TDeltaStore>
inline bool
operator!=(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type const & b)
{
    return !(a == b);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline bool
operator<(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           Iter<TDeltaStore, SpecDeltaStoreIterator_> const & b)
{
    return a._dataIter < b._dataIter;
}

template <typename TDeltaStore>
inline bool
operator<(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
          typename IterComplementConst<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type const & b)
{
    return a._dataIter < b._dataIter;
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline bool
operator<=(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           Iter<TDeltaStore, SpecDeltaStoreIterator_> const & b)
{
    return a._dataIter <= b._dataIter;
}

template <typename TDeltaStore>
inline bool
operator<=(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type const & b)
{
    return a._dataIter <= b._dataIter;
}

// ----------------------------------------------------------------------------
// Function operator>()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline bool
operator>(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           Iter<TDeltaStore, SpecDeltaStoreIterator_> const & b)
{
    return a._dataIter > b._dataIter;
}

template <typename TDeltaStore>
inline bool
operator>(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type const & b)
{
    return a._dataIter > b._dataIter;
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TDeltaStore>
inline bool
operator>=(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           Iter<TDeltaStore, SpecDeltaStoreIterator_> const & b)
{
    return a._dataIter >= b._dataIter;
}

template <typename TDeltaStore>
inline bool
operator>=(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & a,
           typename IterComplementConst<Iter<TDeltaStore, SpecDeltaStoreIterator_> >::Type const & b)
{
    return a._dataIter >= b._dataIter;
}

// ----------------------------------------------------------------------------
// Function addSnpDelta()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline void
addSnpDelta(DeltaStore<TSize, TAlphabet> & variantData,
            TAlphabet const & delta,
            TPosition const & insPos)
{
    typedef typename DeltaType::TValue TValue SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(length(variantData._snpData), (static_cast<TValue>(1) << (BitsPerValue<TValue>::VALUE - 2)));

    appendValue(variantData._snpData, delta);
    insertValue(variantData._varDataMap, insPos, ((length(variantData._snpData) - 1) | DeltaType::DELTA_TYPE_SNP));
}

template <typename TSize, typename TAlphabet>
inline void
addSnpDelta(DeltaStore<TSize, TAlphabet> & variantData,
            TAlphabet const & delta)
{
    addSnpDelta(variantData, delta, length(variantData._varDataMap));
}

// ----------------------------------------------------------------------------
// Function addInsDelta()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TStringSpec, typename TPosition>
inline void
addInsDelta(DeltaStore<TSize, TAlphabet> & variantData,
            String<TAlphabet, TStringSpec> const & delta,
            TPosition const & insPos)
{
    typedef typename DeltaType::TValue TValue SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(length(variantData._insData), (static_cast<TValue>(1) << (BitsPerValue<TValue>::VALUE - 2)));

    appendValue(variantData._insData, delta);
    insertValue(variantData._varDataMap, insPos, ((length(variantData._insData) - 1) | DeltaType::DELTA_TYPE_INS));
}

template <typename TSize, typename TAlphabet, typename TStringSpec>
inline void
addInsDelta(DeltaStore<TSize, TAlphabet> & variantData,
            String<TAlphabet, TStringSpec> const & delta)
{
    addInsDelta(variantData, delta, length(variantData._varDataMap));
}

// ----------------------------------------------------------------------------
// Function addDelDelta()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TSize2, typename TPosition>
inline void
addDelDelta(DeltaStore<TSize, TAlphabet> & variantData,
            TSize2 const & delta,
            TPosition const & insPos)
{
    typedef typename DeltaType::TValue TValue SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(length(variantData._delData), (static_cast<TValue>(1) << (BitsPerValue<TValue>::VALUE - 2)));

    appendValue(variantData._delData, delta);
    insertValue(variantData._varDataMap, insPos, ((length(variantData._delData) - 1) | DeltaType::DELTA_TYPE_DEL));
}

template <typename TSize, typename TAlphabet, typename TSize2>
inline void
addDelDelta(DeltaStore<TSize, TAlphabet> & variantData,
            TSize2 const & delta)
{
    addDelDelta(variantData, delta, length(variantData._varDataMap));
}

// ----------------------------------------------------------------------------
// Function addIndelDelta()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TStringSpec, typename TSize2, typename TPosition>
inline void
addIndelDelta(DeltaStore<TSize, TAlphabet> & variantData,
              Pair<TSize2, String<TAlphabet, TStringSpec> > const & delta,
              TPosition const & insPos)
{
    typedef typename DeltaType::TValue TValue SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT(length(variantData._indelData), (static_cast<TValue>(1) << (BitsPerValue<TValue>::VALUE - 2)));

    appendValue(variantData._indelData, delta);
    insertValue(variantData._varDataMap, insPos, ((length(variantData._indelData) - 1) | DeltaType::DELTA_TYPE_INDEL));
}

template <typename TSize, typename TAlphabet, typename TStringSpec, typename TSize2>
inline void
addIndelDelta(DeltaStore<TSize, TAlphabet> & variantData,
              Pair<TSize2, String<TAlphabet, TStringSpec> > const & delta)
{
    addIndelDelta(variantData, delta, length(variantData._varDataMap));
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename Reference<DeltaStore<TSize, TAlphabet> >::Type
value(DeltaStore<TSize, TAlphabet> & deltaStore,
      TPosition const & pos)
{
    return value(deltaStore._varDataMap, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename Reference<DeltaStore<TSize, TAlphabet> const>::Type
value(DeltaStore<TSize, TAlphabet> const & deltaStore,
      TPosition const & pos)
{
    return value(deltaStore._varDataMap, pos);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline void
clear(DeltaStore<TSize, TAlphabet> & deltaStore)
{
    clear(deltaStore._varDataMap);
    clear(deltaStore._delData);
    clear(deltaStore._indelData);
    clear(deltaStore._insData);
    clear(deltaStore._snpData);
}

// ----------------------------------------------------------------------------
// Function deltaType()
// ----------------------------------------------------------------------------

inline DeltaType::TValue
deltaType(DeltaType::TValue const& val)
{
    return val & DeltaType::MASK_DELTA_TYPE;
}

template <typename TDeltaStore>
inline DeltaType::TValue
deltaType(Iter<TDeltaStore, SpecDeltaStoreIterator_> const& iter)
{
    return value(iter) & DeltaType::MASK_DELTA_TYPE;
}

// ----------------------------------------------------------------------------
// Function deltaPosition()
// ----------------------------------------------------------------------------

inline DeltaType::TValue
deltaPosition(DeltaType::TValue const& val)
{
    return val & DeltaType::MASK_DELTA_POSITION;
}

template <typename TDeltaStore>
inline DeltaType::TValue
deltaPosition(Iter<TDeltaStore, SpecDeltaStoreIterator_> const& iter)
{
    return value(iter) & DeltaType::MASK_DELTA_POSITION;
}

// ----------------------------------------------------------------------------
// Function deltaSnp()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(DeltaStore<TSize, TAlphabet> & store, TPosition const & pos)
{
    return value(store._snpData, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(DeltaStore<TSize, TAlphabet> const & store, TPosition const & pos)
{
    return value(store._snpData, pos);
}

template <typename TDeltaStore>
inline typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter)
{
    return value(iter._containerPtr->_snpData, deltaPosition(iter));
}

template <typename TDeltaStore>
inline typename DeltaValue<TDeltaStore const, DeltaType::DELTA_TYPE_SNP>::Type &
deltaSnp(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & iter)
{
    return value(iter._containerPtr->_snpData, deltaPosition(iter));
}

// ----------------------------------------------------------------------------
// Function deltaDel()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(DeltaStore<TSize, TAlphabet> & store, TPosition const & pos)
{
    return value(store._delData, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(DeltaStore<TSize, TAlphabet> const & store, TPosition const & pos)
{
    return value(store._delData, pos);
}

template <typename TDeltaStore>
inline typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter)
{
    return value(iter._containerPtr->_delData, deltaPosition(iter));
}

template <typename TDeltaStore>
inline typename DeltaValue<TDeltaStore const, DeltaType::DELTA_TYPE_DEL>::Type &
deltaDel(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & iter)
{
    return value(iter._containerPtr->_delData, deltaPosition(iter));
}

// ----------------------------------------------------------------------------
// Function deltaIns()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(DeltaStore<TSize, TAlphabet> & store, TPosition const & pos)
{
    return value(store._insData, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(DeltaStore<TSize, TAlphabet> const & store, TPosition const & pos)
{
    return value(store._insData, pos);
}

template <typename TDeltaStore>
inline typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter)
{
    return value(iter._containerPtr->_insData, deltaPosition(iter));
}

template <typename TDeltaStore>
inline typename DeltaValue<TDeltaStore const, DeltaType::DELTA_TYPE_INS>::Type &
deltaIns(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & iter)
{
    return value(iter._containerPtr->_insData, deltaPosition(iter));
}

// ----------------------------------------------------------------------------
// Function deltaIndel()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet>, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(DeltaStore<TSize, TAlphabet> & store, TPosition const & pos)
{
    return value(store._indelData, pos);
}

template <typename TSize, typename TAlphabet, typename TPosition>
inline typename DeltaValue<DeltaStore<TSize, TAlphabet> const, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(DeltaStore<TSize, TAlphabet> const & store, TPosition const & pos)
{
    return value(store._indelData, pos);
}

template <typename TDeltaStore>
inline typename DeltaValue<TDeltaStore, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(Iter<TDeltaStore, SpecDeltaStoreIterator_> & iter)
{
    return value(iter._containerPtr->_indelData, deltaPosition(iter));
}

template <typename TDeltaStore>
inline typename DeltaValue<TDeltaStore const, DeltaType::DELTA_TYPE_INDEL>::Type &
deltaIndel(Iter<TDeltaStore, SpecDeltaStoreIterator_> const & iter)
{
    return value(iter._containerPtr->_indelData, deltaPosition(iter));
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_STORE_H_
