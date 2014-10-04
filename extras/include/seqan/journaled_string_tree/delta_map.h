// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

template <typename TMap>
struct GetMapValueString_;

template <typename TMap>
struct GetDeltaStore_;

template <typename TMap>
struct GetDeltaCoverageStore_;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag DeltaMapIteratorSpec
// ----------------------------------------------------------------------------

struct SpecDeltaMapIterator_;
typedef Tag<SpecDeltaMapIterator_> DeltaMapIteratorSpec;

// ----------------------------------------------------------------------------
// Struct DeltaMapConfig
// ----------------------------------------------------------------------------

//template <typename TDelSize = __uint32, typename TAlphabet = Dna, typename TSpec = Default>
//struct DeltaMapConfig
//{
//    typedef TDelSize  TDeletionValue;
//    typedef TAlphabet TSnpValue;
//};

// ----------------------------------------------------------------------------
// Tag DeltaMapEntriesMember
// ----------------------------------------------------------------------------

struct DeltaMapEntriesMember_;
typedef Tag<DeltaMapEntriesMember_> DeltaMapEntriesMember;

// ----------------------------------------------------------------------------
// Tag DeltaStoreMember
// ----------------------------------------------------------------------------

struct DeltaMapStoreMember_;
typedef Tag<DeltaMapStoreMember_> DeltaMapStoreMember;

// ----------------------------------------------------------------------------
// Class DeltaMap
// ----------------------------------------------------------------------------

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

    typedef typename Member<DeltaMap, DeltaMapEntriesMember>::Type TDeltaEntries;
    typedef typename Member<DeltaMap, DeltaMapStoreMember>::Type TDeltaStore;
    typedef typename DeltaCoverage<DeltaMap>::Type TCoverage_;
    typedef typename Size<TCoverage_>::Type TCoverageSize;

    TCoverageSize  _coverageSize;
    TDeltaEntries  _entries;
    TDeltaStore    _deltaStore;

    DeltaMap() : _coverageSize(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Member                                  [DeltaMapEntriesMember]
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct Member<DeltaMap<TRefPos, TAlphabet, TSpec>, DeltaMapEntriesMember>
{
    typedef typename Value<DeltaMap<TRefPos, TAlphabet, TSpec> >::Type TValue_;
    typedef String<TValue_> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member                                    [DeltaMapStoreMember]
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct Member<DeltaMap<TRefPos, TAlphabet, TSpec>, DeltaMapStoreMember>
{
    typedef DeltaStore<unsigned, TAlphabet> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member
// ----------------------------------------------------------------------------

// Const version.
template <typename TRefPos, typename TAlphabet, typename TSpec, typename TTag>
struct Member<DeltaMap<TRefPos, TAlphabet, TSpec> const, TTag>
{
    typedef typename Member<DeltaMap<unsigned, TAlphabet, TSpec>, TTag>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// TODO(rmaerker): Enable config struct.
template <typename TRefPos, typename TAlphabet, typename TSpec>
struct Value<DeltaMap<TRefPos, TAlphabet, TSpec> >
{
    typedef typename Member<DeltaMap<TRefPos, TAlphabet, TSpec>, DeltaMapStoreMember>::Type TDeltaStore_;
    typedef typename Size<TDeltaStore_>::Type TSize_;
    typedef DeltaMapEntry<TRefPos, TSize_> Type;
};

template <typename TSize, typename TAlphabet, typename TSpec>
struct Value<DeltaMap<TSize, TAlphabet, TSpec> const>
{
    typedef typename Value<DeltaMap<TSize, TAlphabet, TSpec> >::Type const Type;
};

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

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct Size<DeltaMap<TRefPos, TAlphabet, TSpec> >
{
    typedef typename Member<DeltaMap<TRefPos, TAlphabet, TSpec>, DeltaMapEntriesMember>::Type TEntries;
    typedef typename Size<TEntries>::Type Type;
};

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct Size<DeltaMap<TRefPos, TAlphabet, TSpec> const > :
    Size<DeltaMap<TRefPos, TAlphabet, TSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Iterator                                             [Standard]
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

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct Iterator<DeltaMap<TRefPos, TAlphabet, TSpec>, Standard>
{
    typedef DeltaMap<TRefPos, TAlphabet, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_, DeltaMapIteratorSpec> Type;
};

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct Iterator<DeltaMap<TRefPos, TAlphabet, TSpec> const, Standard>
{
    typedef DeltaMap<TRefPos, TAlphabet, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_ const, DeltaMapIteratorSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator                                               [Rooted]
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct Iterator<DeltaMap<TRefPos, TAlphabet, TSpec>, Rooted>
{
    typedef DeltaMap<TRefPos, TAlphabet, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_, DeltaMapIteratorSpec> Type;
};

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct Iterator<DeltaMap<TRefPos, TAlphabet, TSpec> const, Rooted>
{
    typedef DeltaMap<TRefPos, TAlphabet, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_ const, DeltaMapIteratorSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultGetIteratorSpec
// ----------------------------------------------------------------------------

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct DefaultGetIteratorSpec<DeltaMap<TRefPos, TAlphabet, TSpec> >
{
    typedef DeltaMapIteratorSpec Type;
};

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct DefaultGetIteratorSpec<DeltaMap<TRefPos, TAlphabet, TSpec> const> :
    DefaultGetIteratorSpec<DeltaMap<TRefPos, TAlphabet, TSpec> >{};

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

template <typename TRefPos, typename TAlphabet, typename TSpec, typename TDeltaType>
struct DeltaValue<DeltaMap<TRefPos, TAlphabet, TSpec>, TDeltaType>
{
    typedef DeltaMap<TRefPos, TAlphabet, TSpec> TDeltaMap_;
    typedef typename Member<TDeltaMap_, DeltaMapStoreMember>::Type TDeltaStore_;
    typedef typename DeltaValue<TDeltaStore_, TDeltaType>::Type Type;
};

template <typename TRefPos, typename TAlphabet, typename TSpec, typename TDeltaType>
struct DeltaValue<DeltaMap<TRefPos, TAlphabet, TSpec> const, TDeltaType>
{
    typedef DeltaMap<TRefPos, TAlphabet, TSpec> TDeltaMap_;
    typedef typename Member<TDeltaMap_, DeltaMapStoreMember>::Type TDeltaStore_;
    typedef typename DeltaValue<TDeltaStore_ const, TDeltaType>::Type Type;
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

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct DeltaCoverage<DeltaMap<TRefPos, TAlphabet, TSpec> >
{
    typedef typename Value<DeltaMap<TRefPos, TAlphabet, TSpec> >::Type TEntry_;
    typedef typename DeltaCoverage<TEntry_>::Type Type;
};

template <typename TRefPos, typename TAlphabet, typename TSpec>
struct DeltaCoverage<DeltaMap<TRefPos, TAlphabet, TSpec> const>
{
    typedef typename Value<DeltaMap<TRefPos, TAlphabet, TSpec> >::Type TEntry_;
    typedef typename DeltaCoverage<TEntry_ const>::Type Type;
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

template <typename TRefPos, typename TAlphabet, typename TSpec, typename TPosition>
inline typename Iterator<DeltaMap<TRefPos, TAlphabet, TSpec>, Standard>::Type
find(DeltaMap<TRefPos, TAlphabet, TSpec> & deltaMap, TPosition refPosition)
{
    typedef DeltaMap<TRefPos, TAlphabet, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    SEQAN_ASSERT(!empty(deltaMap));

    TEntry tmp;
    tmp.deltaPosition = refPosition;
    return std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), tmp,
                            DeltaMapEntryCompareLessByDeltaPosition_());
}

template <typename TRefPos, typename TAlphabet, typename TSpec, typename TPosition>
inline typename Iterator<DeltaMap<TRefPos, TAlphabet, TSpec> const, Standard>::Type
find(DeltaMap<TRefPos, TAlphabet, TSpec> const & deltaMap, TPosition refPosition)
{
    typedef DeltaMap<TRefPos, TAlphabet, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    SEQAN_ASSERT(!empty(deltaMap));

    TEntry tmp;
    tmp.deltaPosition = refPosition;
    return std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), tmp,
                            DeltaMapEntryCompareLessByDeltaPosition_());
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
 * @param[in,out] deltaMap  The map to insert the new delta operation for.
 * @param[in]     refPos    The reference position to which the delta operation maps to.
 * @param[in]     deltaVal  The value of the delta operation.
 * @param[in]     cov       The coverage of the delta operation.
 * @param[in]     tag       A specifier to select the correct delta type. One of @link DeltaTypeTags @endlink.
 *
 * @return TIterator An iterator of type @link DeltaMap#Iterator @endlink pointing to the inserted element.
 *
 * @remark The map is implemented as a vector and the insertion time is linear in worst case.
 */

template <typename TValue, typename TAlphabet, typename TSpec, typename TDeltaPos, typename TDeltaValue,
          typename TCoverage, typename TTag>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
insert(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
       TDeltaPos deltaPos,
       TDeltaValue const & deltaValue,
       TCoverage const & coverage,
       TTag const & deltaType)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;
    typedef typename Member<TDeltaMap, DeltaMapStoreMember>::Type TDeltaStore;
    typedef typename Position<TDeltaStore>::Type TStorePos;
    typedef typename Value<TDeltaMap>::Type TEntry;
    typedef typename DeltaRecord<TEntry>::Type TDeltaRecord;

    // First update the coverage size of all delta events if necessary.
    if (getCoverageSize(deltaMap) < length(coverage))
        setCoverageSize(deltaMap, length(coverage));

    TStorePos storePos = addDeltaValue(deltaMap._deltaStore, deltaValue, deltaType);

    if (empty(deltaMap))
    {
        appendValue(deltaMap._entries, TEntry(deltaPos, TDeltaRecord(selectDeltaType(deltaType), storePos), coverage));
        resize(deltaCoverage(begin(deltaMap, Standard())), getCoverageSize(deltaMap), false, Exact());
        return begin(deltaMap, Standard());
    }
    // First we need to search for the insert position.
    TMapIterator mapIt = find(deltaMap, deltaPos);
    insertValue(deltaMap._entries, mapIt - begin(deltaMap, Standard()),
                TEntry(deltaPos, TDeltaRecord(selectDeltaType(deltaType), storePos), coverage));
    resize(deltaCoverage(mapIt), getCoverageSize(deltaMap), false, Exact());
    return mapIt;
}

// ----------------------------------------------------------------------------
// Function begin()                                                  [Standard]
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
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type(deltaMap, 0);
}

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
begin(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap, Standard const &/*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type(deltaMap, 0);
}

// ----------------------------------------------------------------------------
// Function begin()                                                    [Rooted]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Rooted>::Type
begin(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, Rooted const & /*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Rooted>::Type(deltaMap, 0);
}

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Rooted>::Type
begin(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap, Rooted const &/*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Rooted>::Type(deltaMap, 0);
}

// ----------------------------------------------------------------------------
// Function iter()                                                   [Standard]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type
iter(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
     TPos const & pos,
     Standard const & /*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type(deltaMap, pos);
}

template <typename TValue, typename TAlphabet, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
iter(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap,
     TPos const & pos,
     Standard const &/*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type(deltaMap, pos);
}

// ----------------------------------------------------------------------------
// Function iter()                                                     [Rooted]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Rooted>::Type
iter(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, TPos pos, Rooted const & /*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Rooted>::Type(deltaMap, pos);
}

template <typename TValue, typename TAlphabet, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Rooted>::Type
iter(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap, TPos pos, Rooted const & /*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Rooted>::Type(deltaMap, pos);
}

// ----------------------------------------------------------------------------
// Function end()                                                    [Standard]
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
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Standard>::Type(deltaMap, length(deltaMap));
}

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type
end(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap, Standard const &/*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Standard>::Type(deltaMap, length(deltaMap));
}

// ----------------------------------------------------------------------------
// Function end()                                                      [Rooted]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Rooted>::Type
end(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, Rooted const &/*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec>, Rooted>::Type(deltaMap, length(deltaMap));
}

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Rooted>::Type
end(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap, Rooted const &/*tag*/)
{
    return typename Iterator<DeltaMap<TValue, TAlphabet, TSpec> const, Rooted>::Type(deltaMap, length(deltaMap));
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
    deltaMap._coverageSize = 0;
    clear(deltaMap._entries);
    clear(deltaMap._deltaStore);
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
    return empty(deltaMap._entries);
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
    return length(deltaMap._entries);
}

// ----------------------------------------------------------------------------
// Function getCoverageSize()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#getCoverageSize
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the number of sequences covered per delta event.
 *
 * @signature TSize getCoverageSize(deltaMap)
 *
 * @param[in] deltaMap  The map to get the coverage size for.
 *
 * @return TSize The number of sequences covering a delta event.
 */

template <typename TValue, typename TAlphabet, typename TSpec>
inline typename Size<typename DeltaCoverage<DeltaMap<TValue, TAlphabet, TSpec> >::Type>::Type
getCoverageSize(DeltaMap<TValue, TAlphabet, TSpec> const & deltaMap)
{
    return deltaMap._coverageSize;
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
 * @param[in]     size      The new coverage size.
 *
 * This function sets the coverage size globally to all delta events contained in the map.
 */

template <typename TValue, typename TAlphabet, typename TSpec, typename TSize>
inline void
setCoverageSize(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap, TSize const & newSize)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;

    deltaMap._coverageSize = newSize;
    for (TMapIterator it = begin(deltaMap, Standard()); it != end(deltaMap, Standard()); ++it)
        resize(deltaCoverage(it), newSize, Exact());
}

// TODO(rmaerker): Implement erase
}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_H_
