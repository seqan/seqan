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
// Implements the delta map to efficiently store delta entries ordered by
// their position within the base sequence in ascending order.
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

// ----------------------------------------------------------------------------
// Tag DeltaMapIteratorSpec
// ----------------------------------------------------------------------------

struct SpecDeltaMapIterator_;
typedef Tag<SpecDeltaMapIterator_> DeltaMapIteratorSpec;

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


// TODO(rrahn): Add demo.
/*!
 * @class DeltaMap
 *
 * @headerfile <seqan/journaled_string_tree.h>
 *
 * @brief Stores delta information and maps them to a common coordinate system.
 *
 * @signature template <typename TConfig>
 *            class DeltaMap
 * @tparam TConfig   A config type to set the types for the different delta values.
 *
 * This map stores delta events, i.e. replacements, insertions and deletions, for multiple sequences
 * based on a common reference sequence. A bitvector is used to store the coverage of a delta.
 * The types of the correspondinf delta values must be set with the <tt>TConfig<\tt> object, which can be any
 * object which defines the following types: 
 * 
 * <tt>TDeltaPos<\tt>: The value type used to store the position of the delta within the reference.
 * <tt>TDeltaSnp<\tt>: The value type used to store SNPs.
 * <tt>TDeltaDel<\tt>: The value type used to store deletions.
 * <tt>TDeltaIns<\tt>: The value type used to store insertions.
 * <tt>TDeltaSV<\tt>:  The value type used to store structural variants.
 *
 * The delta values are stored in a multi-container. To access a delta value at any given iterator position
 * of the delta map the delta type (see @link DeltaTypeTags @endlink) must be known.
 * The function @link DeltaMapIterator#deltaType @endlink can be used to access the id of the corresponding delta event. 
 * Given the delta type the function @link DeltaMapIterator#deltaValue @endlink can be used to access the corresponding 
 * value.
 */

template <typename TConfig, typename TSpec = Default>
class DeltaMap
{
public:

    typedef typename Member<DeltaMap, DeltaMapEntriesMember>::Type TDeltaEntries;
    typedef typename Member<DeltaMap, DeltaMapStoreMember>::Type TDeltaStore;
    typedef typename DeltaCoverage<DeltaMap>::Type TCoverage_;
    typedef typename Size<TCoverage_>::Type TCoverageSize;

    TDeltaEntries  _entries;
    TDeltaStore    _deltaStore;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Member                                  [DeltaMapEntriesMember]
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec>
struct Member<DeltaMap<TConfig, TSpec>, DeltaMapEntriesMember>
{
    typedef typename Value<DeltaMap<TConfig, TSpec> >::Type TValue_;
    typedef String<TValue_> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member                                    [DeltaMapStoreMember]
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec>
struct Member<DeltaMap<TConfig, TSpec>, DeltaMapStoreMember>
{
    typedef typename TConfig::TSnpValue TSnpValue_;
    typedef typename TConfig::TInsValue TInsValue_;
    typedef typename TConfig::TDelValue TDelValue_;
    typedef typename TConfig::TSVValue TSVValue_;
    typedef DeltaStore<TSnpValue_, TDelValue_, TInsValue_, TSVValue_> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member
// ----------------------------------------------------------------------------

// Const version.
template <typename TConfig, typename TSpec, typename TTag>
struct Member<DeltaMap<TConfig, TSpec> const, TTag>
{
    typedef typename Member<DeltaMap<TConfig, TSpec>, TTag>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec>
struct Value<DeltaMap<TConfig, TSpec> >
{
    typedef typename Member<DeltaMap<TConfig, TSpec>, DeltaMapStoreMember>::Type TDeltaStore_;
    typedef typename Size<TDeltaStore_>::Type TSize_;
    typedef DeltaMapEntry<typename TConfig::TDeltaPos, TSize_> Type;
};

template <typename TConfig, typename TSpec>
struct Value<DeltaMap<TConfig, TSpec> const>
{
    typedef typename Value<DeltaMap<TConfig, TSpec> >::Type const Type;
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

template <typename TConfig, typename TSpec>
struct Size<DeltaMap<TConfig, TSpec> >
{
    typedef typename Member<DeltaMap<TConfig, TSpec>, DeltaMapEntriesMember>::Type TEntries;
    typedef typename Size<TEntries>::Type Type;
};

template <typename TConfig, typename TSpec>
struct Size<DeltaMap<TConfig, TSpec> const > :
    Size<DeltaMap<TConfig, TSpec> >{};

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

template <typename TConfig, typename TSpec>
struct Iterator<DeltaMap<TConfig, TSpec>, Standard>
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_, DeltaMapIteratorSpec> Type;
};

template <typename TConfig, typename TSpec>
struct Iterator<DeltaMap<TConfig, TSpec> const, Standard>
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_ const, DeltaMapIteratorSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator                                               [Rooted]
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec>
struct Iterator<DeltaMap<TConfig, TSpec>, Rooted>
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_, DeltaMapIteratorSpec> Type;
};

template <typename TConfig, typename TSpec>
struct Iterator<DeltaMap<TConfig, TSpec> const, Rooted>
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_ const, DeltaMapIteratorSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaValue
// ----------------------------------------------------------------------------

/*!
 * @mfn DeltaMap#DeltaValue
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns value type for a specific delta.
 *
 * @signature DeltaValue<TDeltaMap, TType>::Type
 * @tparam TDeltaMap The type of the delta map.
 * @tparam TType     The type of the delta value. One of @link DeltaTypeTags @endlink.
 *
 * The delta map stores four different delta events: SNPs, insertions, deletions and variable replacements.
 * This metafunction returns the correct type for the specified event.
 */

template <typename TConfig, typename TSpec, typename TDeltaType>
struct DeltaValue<DeltaMap<TConfig, TSpec>, TDeltaType>
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap_;
    typedef typename Member<TDeltaMap_, DeltaMapStoreMember>::Type TDeltaStore_;
    typedef typename DeltaValue<TDeltaStore_, TDeltaType>::Type Type;
};

template <typename TConfig, typename TSpec, typename TDeltaType>
struct DeltaValue<DeltaMap<TConfig, TSpec> const, TDeltaType>
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap_;
    typedef typename Member<TDeltaMap_, DeltaMapStoreMember>::Type TDeltaStore_;
    typedef typename DeltaValue<TDeltaStore_ const, TDeltaType>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction DeltaCoverage
// ----------------------------------------------------------------------------

/*!
 * @mfn DeltaMap#DeltaCoverage
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns coverage type for a delta map.
 *
 * @signature DeltaCoverage<TDeltaMap>::Type
 * @tparam TDeltaMap The type of the delta map.
 */

template <typename TConfig, typename TSpec>
struct DeltaCoverage<DeltaMap<TConfig, TSpec> >
{
    typedef typename Value<DeltaMap<TConfig, TSpec> >::Type TEntry_;
    typedef typename DeltaCoverage<TEntry_>::Type Type;
};

template <typename TConfig, typename TSpec>
struct DeltaCoverage<DeltaMap<TConfig, TSpec> const>
{
    typedef typename Value<DeltaMap<TConfig, TSpec> >::Type TEntry_;
    typedef typename DeltaCoverage<TEntry_ const>::Type Type;
};


// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::lbWrapper
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec, typename TEntry, typename TFunctor>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
lbWrapper(DeltaMap<TConfig, TSpec> const & deltaMap, TEntry const & entry, TFunctor const & f)
{
    return std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), entry, f);
}

template <typename TConfig, typename TSpec, typename TEntry, typename TFunctor>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
lbWrapper(DeltaMap<TConfig, TSpec> & deltaMap, TEntry const & entry, TFunctor const & f)
{
    return std::lower_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), entry, f);
}

// ----------------------------------------------------------------------------
// Function impl::lowerBound
// ----------------------------------------------------------------------------

// Only searches for the position.
template <typename TConfig, typename TSpec, typename TPosition>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
lowerBound(DeltaMap<TConfig, TSpec> const & deltaMap, TPosition refPosition)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    TEntry entry;
    entry.deltaPosition = refPosition;
    return lbWrapper(deltaMap, entry, DeltaMapEntryPosLessThanComparator_());
}

template <typename TConfig, typename TSpec, typename TPosition>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
lowerBound(DeltaMap<TConfig, TSpec> & deltaMap, TPosition refPosition)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    TEntry entry;
    entry.deltaPosition = refPosition;
    return lbWrapper(deltaMap, entry, DeltaMapEntryPosLessThanComparator_());
}

// Searches for the position and the delta type.
template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
lowerBound(DeltaMap<TConfig, TSpec> const & deltaMap, TPosition refPosition, TDeltaType /*deltaType*/)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    TEntry entry;
    entry.deltaPosition = refPosition;
    entry.deltaRecord.i1 = selectDeltaType(TDeltaType());
    return lbWrapper(deltaMap, entry, DeltaMapEntryPosAndTypeLessThanComparator_());
}

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
lowerBound(DeltaMap<TConfig, TSpec> & deltaMap, TPosition refPosition, TDeltaType /*deltaType*/)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    TEntry entry;
    entry.deltaPosition = refPosition;
    entry.deltaRecord.i1 = selectDeltaType(TDeltaType());
    return lbWrapper(deltaMap, entry, DeltaMapEntryPosAndTypeLessThanComparator_());
}

template <typename TDeltaMap, typename TDeltaPos, typename TTag>
inline bool
checkNoDuplicate(TDeltaMap const & map, TDeltaPos pos, TTag /*deltaType*/)
{
    typedef typename Value<TDeltaMap>::Type TEntry;
    typedef typename DeltaPosition<TEntry>::Type TEntryPos;

    auto it = lowerBound(map, pos, TTag());
    if (it != end(map, Standard()))
        if (getDeltaPosition(*it) == static_cast<TEntryPos>(pos))
            return getDeltaRecord(*it).i1 != selectDeltaType(TTag());
    return true;
}

// ----------------------------------------------------------------------------
// Function impl::insert()
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TDeltaPos, typename TDeltaValue, typename TCoverage, typename TTag>
inline void
insert(Iter<TDeltaMap, DeltaMapIteratorSpec> const & mapIt,
       TDeltaPos deltaPos,
       TDeltaValue const & deltaValue,
       TCoverage const & coverage,
       TTag const & deltaType)
{
    typedef typename Value<TDeltaMap>::Type TEntry;
    typedef typename DeltaRecord<TEntry>::Type TDeltaRecord;

    SEQAN_ASSERT(checkNoDuplicate(*mapIt._mapPtr, deltaPos, deltaType));     // Check valid insert position.

    insertValue(mapIt._mapPtr->_entries, mapIt - begin(*mapIt._mapPtr, Standard()),
                TEntry(deltaPos, TDeltaRecord(selectDeltaType(deltaType),
                                              addDeltaValue(mapIt._mapPtr->_deltaStore, deltaValue, deltaType)), coverage));
}

}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#find
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Finds the element specified by the given delta position and delta type.
 *
 * @signature TIterator find(deltaMap, pos, type)
 *
 * @param[in] deltaMap  The delta map that is searched for the element.
 * @param[in] pos       The delta position to be searched for.
 * @param[in] type      The type of the delta operation. Must be of type @link DeltaTypeTags @endlink.
 *
 * @return TIterator An @link DeltaMap#Iterator @endlink pointing to the corresponding element.
 *  If the key is not contained @link DeltaMap#end @endlink is returned.
 *
 * @remark The runtime is logarithmic in the size of the map.
 */

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
find(DeltaMap<TConfig, TSpec> const & deltaMap, TPosition refPosition, TDeltaType /*deltaType*/)
{
    auto it = impl::lowerBound(deltaMap, refPosition, TDeltaType());
    if (getDeltaPosition(*it) == refPosition && getDeltaRecord(*it).i1 == selectDeltaType(TDeltaType()))
        return it;
    return end(deltaMap, Standard());
}

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
find(DeltaMap<TConfig, TSpec> & deltaMap, TPosition refPosition, TDeltaType /*deltaType*/)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;
    typedef typename DeltaPosition<TEntry>::Type TDeltaPos;
    auto it = impl::lowerBound(deltaMap, refPosition, TDeltaType());
    if (getDeltaPosition(*it) == static_cast<TDeltaPos>(refPosition) &&
        getDeltaRecord(*it).i1 == selectDeltaType(TDeltaType()))
        return it;
    return end(deltaMap, Standard());
}

// ----------------------------------------------------------------------------
// Function insert()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#insert
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Inserts a new delta entry.
 *
 * @signature bool insert(deltaMap, pos, val, cov, type);
 *
 * @param[in,out] deltaMap  The map to insert the new delta operation. Of type @link DeltaMap @endlink.
 * @param[in]     pos       The position of the inserted delta entry.
 * @param[in]     deltaVal  The value of the delta operation.
 * @param[in]     cov       The coverage of the delta operation.
 * @param[in]     type      A specifier to select the correct delta type. One of @link DeltaTypeTags @endlink.
 *
 * @return bool <tt>false<\tt> if an entry with the same <tt>pos<\tt> and <tt>type<\tt> already exists, <tt>true<\tt> otherwise.
 *
 * @remark The map is implemented as a vector and the insertion time is linear in worst case.
 */

template <typename TConfig, typename TSpec, typename TDeltaPos, typename TDeltaValue,
          typename TCoverage, typename TTag>
inline bool
insert(DeltaMap<TConfig, TSpec> & deltaMap,
       TDeltaPos deltaPos,
       TDeltaValue const & deltaValue,
       TCoverage const & coverage,
       TTag const & deltaType)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;
    typedef typename DeltaRecord<TEntry>::Type TDeltaRecord;
    typedef typename DeltaPosition<TEntry>::Type TEntryPos;

    if (SEQAN_UNLIKELY(empty(deltaMap)))
    {
        appendValue(deltaMap._entries,
                    TEntry(deltaPos, TDeltaRecord(selectDeltaType(deltaType),
                                                  addDeltaValue(deltaMap._deltaStore, deltaValue, deltaType)), coverage));
        return true;
    }

    auto it = impl::lowerBound(deltaMap, deltaPos, deltaType);
    if (SEQAN_UNLIKELY(it != end(deltaMap, Standard()) &&
                       getDeltaPosition(*it) == static_cast<TEntryPos>(deltaPos) &&
                       getDeltaRecord(*it).i1 == selectDeltaType(deltaType)))
        return false;
    impl::insert(it, deltaPos, deltaValue, coverage, deltaType);
    return true;
}

// ----------------------------------------------------------------------------
// Function erase()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#erase
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Erases an existing delta entry.
 *
 * @signature bool erase(deltaMap, pos, type);
 *
 * @param[in,out] deltaMap  The map to erase the delta from. Of type @link DeltaMap @endlink.
 * @param[in]     pos       The position of the targeted delta entry.
 * @param[in]     type      The type of the targeted delta entry. One of @link DeltaTypeTags @endlink.
 *
 * @return bool <tt>false<\tt> if such an entry does not exist, <tt>true<\tt> otherwise.
 *
 * @remark The map is implemented as a vector and the insertion time is linear in worst case.
 */

template <typename TConfig, typename TSpec, typename TDeltaPos, typename TDeltaType>
inline bool
erase(DeltaMap<TConfig, TSpec> & deltaMap,
      TDeltaPos deltaPos,
      TDeltaType /*deltaType*/)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Member<TDeltaMap, DeltaMapStoreMember>::Type TStore;
    typedef typename Value<TDeltaMap>::Type TEntry;
    typedef typename Size<TStore>::Type TStoreSize;

    if (SEQAN_UNLIKELY(empty(deltaMap)))  // Check for empty deltaMap.
        return false;

    auto it = find(deltaMap, deltaPos, TDeltaType());
    if (SEQAN_UNLIKELY(it == end(deltaMap, Standard())))
        return false;  // Element does not exists.

    // 1. Erase the delta record from the corresponding delta store.
    TStoreSize storePos = getDeltaRecord(value(it)).i2;
    SEQAN_ASSERT_LT(storePos, length(getDeltaStore(deltaMap._deltaStore, TDeltaType())));
    eraseDeltaValue(deltaMap._deltaStore, storePos, TDeltaType());
    // 2. Erase corresponding entry.
    erase(deltaMap._entries, it - begin(deltaMap, Standard()));
    // 3. Update record position for all entries
    forEach(deltaMap, [storePos](TEntry & entry)
    {
        if (getDeltaRecord(entry).i1 == selectDeltaType(TDeltaType()) && getDeltaRecord(entry).i2 > storePos)
            --getDeltaRecord(entry).i2;  // Decrease record position by one.
    });
    return true;
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

template <typename TConfig, typename TSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
begin(DeltaMap<TConfig, TSpec> & deltaMap, Standard const & /*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type(deltaMap, 0);
}

template <typename TConfig, typename TSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
begin(DeltaMap<TConfig, TSpec> const & deltaMap, Standard const &/*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type(deltaMap, 0);
}

// ----------------------------------------------------------------------------
// Function begin()                                                    [Rooted]
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Rooted>::Type
begin(DeltaMap<TConfig, TSpec> & deltaMap, Rooted const & /*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec>, Rooted>::Type(deltaMap, 0);
}

template <typename TConfig, typename TSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Rooted>::Type
begin(DeltaMap<TConfig, TSpec> const & deltaMap, Rooted const &/*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec> const, Rooted>::Type(deltaMap, 0);
}

// ----------------------------------------------------------------------------
// Function iter()                                                   [Standard]
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
iter(DeltaMap<TConfig, TSpec> & deltaMap,
     TPos const & pos,
     Standard const & /*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type(deltaMap, pos);
}

template <typename TConfig, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
iter(DeltaMap<TConfig, TSpec> const & deltaMap,
     TPos const & pos,
     Standard const &/*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type(deltaMap, pos);
}

// ----------------------------------------------------------------------------
// Function iter()                                                     [Rooted]
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Rooted>::Type
iter(DeltaMap<TConfig, TSpec> & deltaMap, TPos pos, Rooted const & /*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec>, Rooted>::Type(deltaMap, pos);
}

template <typename TConfig, typename TSpec, typename TPos>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Rooted>::Type
iter(DeltaMap<TConfig, TSpec> const & deltaMap, TPos pos, Rooted const & /*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec> const, Rooted>::Type(deltaMap, pos);
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

template <typename TConfig, typename TSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
end(DeltaMap<TConfig, TSpec> & deltaMap, Standard const &/*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type(deltaMap, size(deltaMap));
}

template <typename TConfig, typename TSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
end(DeltaMap<TConfig, TSpec> const & deltaMap, Standard const &/*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type(deltaMap, size(deltaMap));
}

// ----------------------------------------------------------------------------
// Function end()                                                      [Rooted]
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Rooted>::Type
end(DeltaMap<TConfig, TSpec> & deltaMap, Rooted const &/*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec>, Rooted>::Type(deltaMap, size(deltaMap));
}

template <typename TConfig, typename TSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Rooted>::Type
end(DeltaMap<TConfig, TSpec> const & deltaMap, Rooted const &/*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec> const, Rooted>::Type(deltaMap, size(deltaMap));
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

template <typename TConfig, typename TSpec>
inline void
clear(DeltaMap<TConfig, TSpec> & deltaMap)
{
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

template <typename TConfig, typename TSpec>
inline bool
empty(DeltaMap<TConfig, TSpec> const & deltaMap)
{
    return empty(deltaMap._entries);
}

// ----------------------------------------------------------------------------
// Function size()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#size
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the number of mapped delta events.
 *
 * @signature TSize size(deltaMap)
 *
 * @param[in] deltaMap  The map to get the size for.
 *
 * @return TSize The number of delta events stored in the map.
 */

template <typename TConfig, typename TSpec>
inline typename Size<DeltaMap<TConfig, TSpec> >::Type
size(DeltaMap<TConfig, TSpec> const & deltaMap)
{
    return length(deltaMap._entries);
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_H_
