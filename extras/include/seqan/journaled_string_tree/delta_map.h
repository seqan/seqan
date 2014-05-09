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
// Implements the variant map. This is a facade that contains the delta
// store and the delta coverage store. It maps to each reference position
// the corresponding delta and coverage.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct SpecDeltaMapIterator_;
typedef Tag<SpecDeltaMapIterator_> DeltaMapIteratorSpec;

template <typename TMap>
struct GetMapValueString_;

template <typename TMap>
struct GetDeltaStore_;

template <typename TMap>
struct GetDeltaCoverageStore_;

/*!
 * @class DeltaMap
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Stores delta information and maps them to a common coordinate system.
 *
 * @signature template <typename TValue, typename TAlphabet>
 *            class DeltaMap
 * @tparam TValue The value type of the keys.
 * @tparam TAlphabet The alphabet typed used to for SNPs and insertions.
 *
 * This map stores biological delta events, i.e. replacements, insertions and deletions, for multiple sequences
 * based to a common reference sequence. A bitvector is used to denote the coverage of a delta.
 * Note, the keys must be comparable using the less-than operator.
 *
 * The delta events are stored in a multi-container fashion. To access a delta event at any given iterator position
 * of the delta map the delta type must be known beforehand. The function @link DeltaMapIterator#deltaType @endlink can
 * be used to access the id of the corresponding delta event. Based on the type the correct position must be called
 * to get access to the stored value.
 *
 * Deletions are stored as the length of the deletion. Insertions are stored in an concatenated string set internally
 * and replacements are stored as a pair of a deletion length as first parameter and an insertion string as second.
 * Thus variable replacements, where the deletion length can differ from the length of the insertion.
 */

template <typename TValue, typename TAlphabet, typename TSpec = Default>
class DeltaMap
{
public:

    typedef typename GetMapValueString_<DeltaMap>::Type TKeys;
    typedef typename GetDeltaStore_<DeltaMap>::Type TDeltaStore_;
    typedef typename GetDeltaCoverageStore_<DeltaMap>::Type TCoverageStore;

    TKeys            _keys;  // Key string containing the positions within the reference.
    TDeltaStore_     _deltaStore;
    TCoverageStore   _deltaCoverageStore;

    DeltaMap() : _keys(), _deltaStore(), _deltaCoverageStore()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetMapValueString_
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetMapValueString_<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef String<TValue> Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetMapValueString_<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef String<TValue> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetDeltaStore_
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetDeltaStore_<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename Size<TDeltaMap_>::Type TSize_;
    typedef DeltaStore<TSize_, TAlphabet> Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetDeltaStore_<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename Size<TDeltaMap_>::Type TSize_;
    typedef DeltaStore<TSize_, TAlphabet> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetDeltaCoverageStore_
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetDeltaCoverageStore_<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef DeltaCoverageStore<> Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct GetDeltaCoverageStore_<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef DeltaCoverageStore<> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

/*!
 * @mfn DeltaMap#Value
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns type of the reference coordinate system.
 *
 * @signature Value<TDeltaMap>::Type
 * @tparam TDeltaMap The type to query the value type for.
 * @return TValue The value type to use for <tt>TDeltaMap</tt>.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
struct Value<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Value<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

/*!
 * @mfn DeltaMap#Reference
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns reference type of the reference coordinate system.
 *
 * @signature Reference<TDeltaMap>::Type
 * @tparam TDeltaMap The type to query the reference type for.
 * @return TReference The reference type to use for <tt>TDeltaMap</tt>.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
struct Reference<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef TValue & Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Reference<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef TValue const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/*!
 * @mfn DeltaMap#Size
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns size type for a delta map.
 *
 * @signature Size<TDeltaMap>::Type
 * @tparam TDeltaMap The type to query the size type for.
 * @return TSize The size type to use for <tt>TDeltaMap</tt>.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
struct Size<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename GetMapValueString_<TDeltaMap_>::Type TKeys_;

    typedef typename Size<TKeys_>::Type Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Size<DeltaMap<TValue, TAlphabet, TSpec> const > :
    Size<DeltaMap<TValue, TAlphabet, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

/*!
 * @mfn DeltaMap#Position
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns position type for a delta map.
 *
 * @signature Position<TDeltaMap>::Type
 * @tparam TDeltaMap The type to query the position type for.
 * @return TPosition The position type to use for <tt>TDeltaMap</tt>.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
struct Position<DeltaMap<TValue, TAlphabet, TSpec> >  :
    Size<DeltaMap<TValue, TAlphabet, TSpec> >{};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Position<DeltaMap<TValue, TAlphabet, TSpec> const >  :
    Size<DeltaMap<TValue, TAlphabet, TSpec> const>{};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

/*!
 * @mfn DeltaMap#Iterator
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns iterator type for a delta map.
 *
 * @signature Iterator<TDeltaMap, Standard>::Type
 * @tparam TDeltaMap The type to query the iterator type for.
 * @return TIterator The iterator type to use for <tt>TDeltaMap</tt>. @link DeltaMapIterator @endlink.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
struct Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_, DeltaMapIteratorSpec> Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_ const, DeltaMapIteratorSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultGetIteratorSpec
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
struct DefaultGetIteratorSpec<DeltaMap< TValue, TAlphabet, TSpec> >
{
    typedef DeltaMapIteratorSpec Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct DefaultGetIteratorSpec<DeltaMap< TValue, TAlphabet, TSpec> const> :
    DefaultGetIteratorSpec<DeltaMap< TValue, TAlphabet, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue
// ----------------------------------------------------------------------------

/*!
 * @mfn DeltaMap#DeltaValue
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns value type of a specific delta.
 *
 * @signature DeltaValue<TDeltaMap, ID>::Type
 * @tparam TDeltaMap The type to query the iterator type for.
 * @tparam ID        Id to specify the delta type. One of @link DeltaType @endlink.
 *
 * The delta map stores four different delta events: SNPs, insertions, deletions and variable replacements.
 * This metafunction returns the correct type for the specified event.
 */

template <typename TValue, typename TAlphabet, typename TSpec, typename DeltaType::TValue TYPE>
struct DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, TYPE>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename TDeltaMap_::TDeltaStore_ TDeltaStore_;
    typedef typename DeltaValue<TDeltaStore_, TYPE>::Type Type;
};

template <typename TValue, typename TAlphabet, typename TSpec, typename DeltaType::TValue TYPE>
struct DeltaValue<DeltaMap<TValue, TAlphabet, TSpec> const, TYPE>
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap_;
    typedef typename TDeltaMap_::TDeltaStore_ TDeltaStore_;
    typedef typename DeltaValue<TDeltaStore_ const, TYPE>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaCoverage
// ----------------------------------------------------------------------------

/*!
 * @mfn DeltaMap#DeltaCoverage
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns coverage type of a delta map.
 *
 * @signature DeltaCoverage<TDeltaMap, KEY>::Type
 * @tparam TDeltaMap The type to query the iterator type for.
 */

template <typename TMap>
struct DeltaCoverage;

template <typename TValue, typename TAlphabet, typename TSpec>
struct DeltaCoverage<DeltaMap<TValue, TAlphabet, TSpec> >
{
    typedef typename Value<DeltaCoverageStore<> >::Type Type;
};

template <typename TValue, typename TAlphabet, typename TSpec>
struct DeltaCoverage<DeltaMap<TValue, TAlphabet, TSpec> const>
{
    typedef typename Value<DeltaCoverageStore<> >::Type const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#find
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Locates the given key within the map.
 *
 * @signature TIterator find(deltaMap, key)
 *
 * @param[in] deltaMap The delta map to search for the key.
 * @param[in] key   The key to be searched.
 *
 * @return TIterator An @link DeltaMap#Iterator @endlink pointing to the first mapped entry that compares equal to the key.
 *  If the key is not contained @link DeltaMap#end @endlink is returned.
 *
 * @note The runtime is logarithmic in the size of the map.
 */

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
find(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
     TKey const & key)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;

    SEQAN_ASSERT(!empty(deltaMap));

    TMapIterator it = std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), key);
    return (*it == key) ? it : end(deltaMap, Standard());
}

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
find(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap,
     TKey const & key)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> const TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;

    SEQAN_ASSERT(!empty(deltaMap));

    TMapIterator it = std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), key);
    return (*it == key) ? it : end(deltaMap, Standard());
}

// ----------------------------------------------------------------------------
// Function _insert()                                                     [DEL]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey, typename TPosition>
inline void
_insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
       TKey const & key,
       TPosition const & insPos,
       typename DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, DeltaType::DELTA_TYPE_DEL>::Type const & delta)
{
    SEQAN_ASSERT_LEQ(insPos, length(deltaMap));

    insertValue(deltaMap._keys, insPos, key);  // TODO(rmaerker): Adapt to binary search later.
    addDelDelta(deltaMap._deltaStore, delta, insPos);
}

// ----------------------------------------------------------------------------
// Function _insert()                                                     [SNP]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey, typename TPosition>
inline void
_insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
        TKey const & key,
        TPosition const & insPos,
        typename DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, DeltaType::DELTA_TYPE_SNP>::Type const & delta)
{
    SEQAN_ASSERT_LEQ(insPos, length(deltaMap));

    insertValue(deltaMap._keys, insPos, key);  // TODO(rmaerker): Adapt to binary search later.
    addSnpDelta(deltaMap._deltaStore, delta, insPos);
}

// ----------------------------------------------------------------------------
// Function _insert()                                                     [INS]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey, typename TPosition>
inline void
_insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
        TKey const & key,
        TPosition const & insPos,
        typename DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, DeltaType::DELTA_TYPE_INS>::Type const & delta)
{
    SEQAN_ASSERT_LEQ(insPos, length(deltaMap));

    insertValue(deltaMap._keys, insPos, key);
    addInsDelta(deltaMap._deltaStore, delta, insPos);
}

// ----------------------------------------------------------------------------
// Function _insert()                                                   [INDEL]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey, typename TPosition>
inline void
_insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
        TKey const & key,
        TPosition const & insPos,
        typename DeltaValue<DeltaMap<TValue, TAlphabet, TSpec>, DeltaType::DELTA_TYPE_INDEL>::Type const & delta)
{
    SEQAN_ASSERT_LEQ(insPos, length(deltaMap));

    insertValue(deltaMap._keys, insPos, key);
    addIndelDelta(deltaMap._deltaStore, delta, insPos);
}

// ----------------------------------------------------------------------------
// Function insert()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#insert
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Inserts a new delta event.
 *
 * @signature TIterator insert(deltaMap, key, delta, cov)
 *
 * @param[in,out] deltaMap  The map to insert the new delta.
 * @param[in]     key       The key the new event maps to.
 * @param[in]     delta     The delta event of type @link DeltaMap#DeltaValue @endlink.
 * @param[in]     cov       The coverage of this delta event.
 *
 * @return TIterator An iterator of type @link DeltaMap#Iterator @endlink pointing to the inserted element.
 *
 * @note At the moment the map is implemented as a vector, thus the runtime is linear in worst case.
 */

template <typename TValue, typename TAlphabet, typename TSpec, typename TKey,
          typename TDelta>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
       TKey const & key,
       TDelta const & delta,
       typename DeltaCoverage<DeltaMap<TValue, TAlphabet, TSpec> >::Type const & deltaCoverage)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Position<TDeltaMap>::Type TPosition;

    TPosition insPos = std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), key) -
                       begin(deltaMap, Standard());
    addCoverage(deltaMap._deltaCoverageStore, deltaCoverage, insPos);
    _insert(deltaMap, key, insPos, delta);
    return iter(deltaMap, insPos, Standard());
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#begin
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns an iterator pointing to the beginning of the map.
 *
 * @signature TIterator begin(deltaMap, tag)
 *
 * @param[in] deltaMap  The map to get the iterator for.
 * @param[in] tag       The iterator tag. Of type @link ContainerIteratorTags @endlink.
 *
 * @return TIterator An iterator of type @link DeltaMap#Iterator @endlink pointing to the beginning of the map.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
begin(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, Standard const & /*tag*/)
{
    typedef typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type TIterator;
    TIterator tmp;
    _initBegin(tmp, deltaMap);
    return tmp;
}

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
begin(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap, Standard const &/*tag*/)
{
    typedef typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type TIterator;
    TIterator tmp;
    _initBegin(tmp, deltaMap);
    return tmp;
}

// ----------------------------------------------------------------------------
// Function iter()
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
iter(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
     TPos const & pos,
     Standard const & /*tag*/)
{
    return begin(deltaMap, Standard()) + pos;
}

template <typename TValue, typename TAlphabet, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
iter(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap,
     TPos const & pos,
     Standard const &/*tag*/)
{
    return begin(deltaMap, Standard()) + pos;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#end
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns an iterator pointing to the end of the map.
 *
 * @signature TIterator end(deltaMap, tag)
 *
 * @param[in] deltaMap  The map to get the iterator for.
 * @param[in] tag       The iterator tag. Of type @link ContainerIteratorTags @endlink.
 *
 * @return TIterator An iterator of type @link DeltaMap#Iterator @endlink pointing to the end of the map.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
end(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, Standard const &/*tag*/)
{
    typedef typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type TIterator;
    TIterator tmp;
    _initEnd(tmp, deltaMap);
    return tmp;
}

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
end(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap, Standard const &/*tag*/)
{
    typedef typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type TIterator;
    TIterator tmp;
    _initEnd(tmp, deltaMap);
    return tmp;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#clear
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Clears the delta map.
 *
 * @signature clear(deltaMap)
 *
 * @param[in,out] deltaMap  The map to be cleared.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
inline void
clear(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap)
{
    clear(deltaMap._keys);
    clear(deltaMap._deltaStore);
    clear(deltaMap._deltaCoverageStore);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#empty
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Checks if the delta map is empty.
 *
 * @signature bool empty(deltaMap)
 *
 * @param[in] deltaMap  The map to be checked for.
 *
 * @return bool <tt>true</tt> if empty, otherwise <tt>false</tt>
 */

template <typename TValue, typename TAlphabet, typename TSpec>
inline bool
empty(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap)
{
    return empty(deltaMap._keys);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#length
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the number of mapped delta events.
 *
 * @signature TSize length(deltaMap)
 *
 * @param[in] deltaMap  The map to get the length for.
 *
 * @return TSize The number of delta events stored in the map.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Size<DeltaMap<TValue, TAlphabet, TSpec> >::Type
length(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap)
{
    return length(deltaMap._keys);
}

// ----------------------------------------------------------------------------
// Function coverageSize()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#coverageSize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the number of sequences covered per delta event.
 *
 * @signature TSize coverageSize(deltaMap)
 *
 * @param[in] deltaMap  The map to get the coverage size for.
 *
 * @return TSize The number of sequences covering a delta event.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Size<DeltaMap<TValue, TAlphabet, TSpec> >::Type
coverageSize(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap)
{
    return coverageSize(deltaMap._deltaCoverageStore);
}

// ----------------------------------------------------------------------------
// Function setCoverageSize()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#setCoverageSize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Sets the number of sequences covering a delta event.
 *
 * @signature  setCoverageSize(deltaMap, size)
 *
 * @param[in,out] deltaMap  The map to set the coverage size for.
 * @param[in]      size      The new coverage size.
 *
 * This function sets the coverage size globally to all delta events contained in the map.
 */

template <typename TValue, typename TAlphabet, typename TSpec, typename TSize>
inline void
setCoverageSize(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, TSize const & size)
{
    return setCoverageSize(deltaMap._deltaCoverageStore, size);
}

// TODO(rmaerker): Implement erase
}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_H_
