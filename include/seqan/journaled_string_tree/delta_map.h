// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
 * The types of the corresponding delta values must be set with the <tt>TConfig<\tt> object, which can be any
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
 *
 * The delta map implements the interfaces of the <b>AssociativeContainerConcept<\b> and is a multi-map.
 */

template <typename TConfig, typename TSpec = Default>
class DeltaMap
{
public:

    typedef typename Member<DeltaMap, DeltaMapEntriesMember>::Type TDeltaEntries;
    typedef typename Member<DeltaMap, DeltaMapStoreMember>::Type TDeltaStore;

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
    typedef impl::DeltaStore<TSnpValue_, TDelValue_, TInsValue_, TSVValue_> Type;
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

/*!
 * @mfn DeltaMap#Value
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns value type for the delta map.
 *
 * @signature Value<TDeltaMap>::Type
 * @tparam TDeltaMap The type to query the value type for.
 * @return TValue The value type to use for <tt>TDeltaMap</tt>.
 */

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

//template <typename TConfig, typename TSpec>
//struct Iterator<DeltaMap<TConfig, TSpec>, Standard>
//{
//    typedef DeltaMap<TConfig, TSpec> TDeltaMap_;
//    typedef Iter<TDeltaMap_, DeltaMapIteratorSpec> Type;
//};
//
//template <typename TConfig, typename TSpec>
//struct Iterator<DeltaMap<TConfig, TSpec> const, Standard>
//{
//    typedef DeltaMap<TConfig, TSpec> TDeltaMap_;
//    typedef Iter<TDeltaMap_ const, DeltaMapIteratorSpec> Type;
//};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec, typename TIteratorSpec>
struct Iterator<DeltaMap<TConfig, TSpec>, Tag<TIteratorSpec> const>
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap_;
    typedef Iter<TDeltaMap_, DeltaMapIteratorSpec> Type;
};

template <typename TConfig, typename TSpec, typename TIteratorSpec>
struct Iterator<DeltaMap<TConfig, TSpec> const, Tag<TIteratorSpec> const>
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
 * This metafunction returns the correct type for the specified event given the delta type tag.
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
// Function impl::ubWrapper
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec, typename TEntry, typename TFunctor>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
ubWrapper(DeltaMap<TConfig, TSpec> const & deltaMap, TEntry const & entry, TFunctor const & f)
{
    return std::upper_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), entry, f);
}

template <typename TConfig, typename TSpec, typename TEntry, typename TFunctor>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
ubWrapper(DeltaMap<TConfig, TSpec> & deltaMap, TEntry const & entry, TFunctor const & f)
{
    return std::upper_bound(begin(deltaMap, Standard()), end(deltaMap, Standard()), entry, f);
}

// ----------------------------------------------------------------------------
// Function impl::lowerBound();
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
lowerBound(DeltaMap<TConfig, TSpec> const & deltaMap,
           TPosition const refPosition,
           DeltaEndType const endType,
           TDeltaType /*deltaType*/)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    TEntry entry;
    entry.deltaPosition = refPosition;
    entry.deltaRecord.i1 = selectDeltaType(TDeltaType());
    entry.deltaTypeEnd = endType;
    return impl::lbWrapper(deltaMap, entry, DeltaMapEntryPosAndTypeLessThanComparator_());
}

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
lowerBound(DeltaMap<TConfig, TSpec> & deltaMap,
           TPosition const refPosition,
           DeltaEndType const endType,
           TDeltaType /*deltaType*/)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    TEntry entry;
    entry.deltaPosition = refPosition;
    entry.deltaRecord.i1 = selectDeltaType(TDeltaType());
    entry.deltaTypeEnd = endType;
    return impl::lbWrapper(deltaMap, entry, DeltaMapEntryPosAndTypeLessThanComparator_());
}

// ----------------------------------------------------------------------------
// Function impl::upperBound();
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
upperBound(DeltaMap<TConfig, TSpec> const & deltaMap,
           TPosition const refPosition,
           DeltaEndType const endType,
           TDeltaType /*deltaType*/)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    TEntry entry;
    entry.deltaPosition = refPosition;
    entry.deltaRecord.i1 = selectDeltaType(TDeltaType());
    entry.deltaTypeEnd = endType;
    return impl::ubWrapper(deltaMap, entry, DeltaMapEntryPosAndTypeLessThanComparator_());
}

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
upperBound(DeltaMap<TConfig, TSpec> & deltaMap,
           TPosition refPosition,
           DeltaEndType const endType,
           TDeltaType /*deltaType*/)
{
    typedef DeltaMap<TConfig, TSpec> TDeltaMap;
    typedef typename Value<TDeltaMap>::Type TEntry;

    TEntry entry;
    entry.deltaPosition = refPosition;
    entry.deltaRecord.i1 = selectDeltaType(TDeltaType());
    entry.deltaTypeEnd = endType;
    return impl::ubWrapper(deltaMap, entry, DeltaMapEntryPosAndTypeLessThanComparator_());
}

// ----------------------------------------------------------------------------
// Function impl::find();
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
find(DeltaMap<TConfig, TSpec> const & deltaMap,
     TPosition const refPosition,
     DeltaEndType const endType,
     TDeltaType /*deltaType*/)
{
    auto it = lowerBound(deltaMap, refPosition, endType, TDeltaType());
    if (it != end(deltaMap, Standard()) &&
        static_cast<TPosition>(getDeltaPosition(*it)) == refPosition &&
        getDeltaRecord(*it).i1 == selectDeltaType(TDeltaType()))
        return it;
    return end(deltaMap, Standard());
}

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
find(DeltaMap<TConfig, TSpec> & deltaMap,
     TPosition const refPosition,
     DeltaEndType const endType,
     TDeltaType /*deltaType*/)
{
    auto it = lowerBound(deltaMap, refPosition, endType, TDeltaType());
    if (it != end(deltaMap, Standard()) &&
        static_cast<TPosition>(getDeltaPosition(*it)) == refPosition &&
        getDeltaRecord(*it).i1 == selectDeltaType(TDeltaType()))
        return it;
    return end(deltaMap, Standard());
}

// ----------------------------------------------------------------------------
// Function impl::checkNoDuplicate()
// ----------------------------------------------------------------------------

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

    insertValue(mapIt._mapPtr->_entries, mapIt - begin(*mapIt._mapPtr, Standard()),
                TEntry(deltaPos, TDeltaRecord(selectDeltaType(deltaType),
                                              addDeltaValue(mapIt._mapPtr->_deltaStore, deltaValue, deltaType)),
                                              coverage,
                                              DeltaEndType::IS_BOTH));
}

template <typename TDeltaMap, typename TDeltaPos, typename TDeltaValue, typename TCoverage>
inline void
insert(Iter<TDeltaMap, DeltaMapIteratorSpec> mapIt,
       TDeltaPos deltaPos,
       TDeltaValue const & deltaValue,
       TCoverage const & coverage,
       DeltaTypeDel const & /*deltaType*/)
{
    typedef typename Value<TDeltaMap>::Type TEntry;
    typedef typename DeltaRecord<TEntry>::Type TDeltaRecord;

    DeltaEndType endType = DeltaEndType::IS_BOTH;
    if (deltaValue > 1)
        endType = DeltaEndType::IS_LEFT;

    auto storePos = addDeltaValue(mapIt._mapPtr->_deltaStore, deltaValue, DeltaTypeDel());
    insertValue(mapIt._mapPtr->_entries, mapIt - begin(*mapIt._mapPtr, Standard()),
                TEntry(deltaPos, TDeltaRecord(DELTA_TYPE_DEL, storePos), coverage, endType));

    if (endType == DeltaEndType::IS_LEFT)
    {
        deltaPos += deltaValue - 1;
        mapIt = lowerBound(*mapIt._mapPtr, deltaPos, DeltaTypeDel());
        insertValue(mapIt._mapPtr->_entries, mapIt - begin(*mapIt._mapPtr, Standard()),
                    TEntry(deltaPos, TDeltaRecord(DELTA_TYPE_DEL, storePos), coverage, DeltaEndType::IS_RIGHT));
    }
}

template <typename TDeltaMap, typename TDeltaPos, typename TDeltaValue, typename TCoverage>
inline void
insert(Iter<TDeltaMap, DeltaMapIteratorSpec> mapIt,
       TDeltaPos deltaPos,
       TDeltaValue const & deltaValue,
       TCoverage const & coverage,
       DeltaTypeSV const & /*deltaType*/)
{
    typedef typename Value<TDeltaMap>::Type TEntry;
    typedef typename DeltaRecord<TEntry>::Type TDeltaRecord;

    DeltaEndType endType = DeltaEndType::IS_BOTH;
    if (deltaValue.i1 > 1)
        endType = DeltaEndType::IS_LEFT;

    auto storePos = addDeltaValue(mapIt._mapPtr->_deltaStore, deltaValue, DeltaTypeSV());
    insertValue(mapIt._mapPtr->_entries, mapIt - begin(*mapIt._mapPtr, Standard()),
                TEntry(deltaPos, TDeltaRecord(DELTA_TYPE_SV, storePos), coverage, endType));

    if (endType == DeltaEndType::IS_LEFT)
    {
        deltaPos += deltaValue.i1 - 1;
        mapIt = lowerBound(*mapIt._mapPtr, deltaPos, DeltaTypeSV());
        insertValue(mapIt._mapPtr->_entries, mapIt - begin(*mapIt._mapPtr, Standard()),
                    TEntry(deltaPos, TDeltaRecord(DELTA_TYPE_SV, storePos), coverage, DeltaEndType::IS_RIGHT));
    }
}

}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function lowerBound()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#lowerBound
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Finds the first element that compares not less than the specified key.
 *
 * @signature TIterator lowerBound(deltaMap, pos, type)
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

// Searches for the position and the delta type.

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type
lowerBound(DeltaMap<TConfig, TSpec> const & deltaMap,
           TPosition const refPosition,
           TDeltaType /*deltaType*/)
{
    return impl::lowerBound(deltaMap, refPosition, DeltaEndType::IS_LEFT, TDeltaType());
}

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
lowerBound(DeltaMap<TConfig, TSpec> & deltaMap,
           TPosition const refPosition,
           TDeltaType /*deltaType*/)
{
    return impl::lowerBound(deltaMap, refPosition, DeltaEndType::IS_LEFT, TDeltaType());
}

// ----------------------------------------------------------------------------
// Function upperBound()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#upperBound
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Finds the first element that compares not less or equal to the specified key.
 *
 * @signature TIterator upperBound(deltaMap, pos, type)
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
upperBound(DeltaMap<TConfig, TSpec> const & deltaMap,
           TPosition const refPosition,
           TDeltaType /*deltaType*/)
{
    return impl::upperBound(deltaMap, refPosition, DeltaEndType::IS_BOTH, TDeltaType());
}

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
upperBound(DeltaMap<TConfig, TSpec> & deltaMap,
           TPosition const refPosition,
           TDeltaType /*deltaType*/)
{
    return impl::upperBound(deltaMap, refPosition, DeltaEndType::IS_BOTH, TDeltaType());
}

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
find(DeltaMap<TConfig, TSpec> const & deltaMap,
     TPosition const refPosition,
     TDeltaType /*deltaType*/)
{
    return impl::find(deltaMap, refPosition, DeltaEndType::IS_LEFT, TDeltaType());
}

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Standard>::Type
find(DeltaMap<TConfig, TSpec> & deltaMap,
     TPosition const refPosition,
     TDeltaType /*deltaType*/)
{
    return impl::find(deltaMap, refPosition, DeltaEndType::IS_LEFT, TDeltaType());
}

// ----------------------------------------------------------------------------
// Function count()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#count
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Counts the number of elements that compare equal to the specified key.
 *
 * @signature TSize count(deltaMap, pos, type)
 *
 * @param[in] deltaMap  The delta map that is searched for the element.
 * @param[in] pos       The delta position to be searched for.
 * @param[in] type      The type of the delta operation. Must be of type @link DeltaTypeTags @endlink.
 *
 * @return TSize the number of elements with the specified key. Of type @link DeltaMap#Size @endlink.
 *
 * @remark The runtime is logarithmic in the size of the map.
 */

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline typename Size<DeltaMap<TConfig, TSpec> >::Type
count(DeltaMap<TConfig, TSpec> const & deltaMap,
      TPosition refPosition,
      Tag<TDeltaType> /*deltaType*/)
{
    auto itB = lowerBound(deltaMap, refPosition, Tag<TDeltaType>());
    auto count = 0;
    while (itB != end(deltaMap, Standard()) &&
           static_cast<TPosition>(getDeltaPosition(*itB)) == refPosition &&
           getDeltaRecord(*itB).i1 == selectDeltaType(Tag<TDeltaType>()))
    {
        ++count;
        ++itB;
    }
    return count;
}

// ----------------------------------------------------------------------------
// Function equalRange()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#equalRange
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief Returns the range over all elements comparing equal to the specified key.
 *
 * @signature Pair<TIterator> equalRange(deltaMap, pos, type)
 *
 * @param[in] deltaMap  The delta map that is searched for the element.
 * @param[in] pos       The delta position to be searched for.
 * @param[in] type      The type of the delta operation. Must be of type @link DeltaTypeTags @endlink.
 *
 * @return Pair<TIterator> A @link Pair @endlink of iterator types @link DeltaMap#Iterator @endlink. The first value points
 *  to the first element that compares not less than the specified key or to the @link DeltaMap#end @endlink if such an elment could not be found.
 * The second value points to the first element that does not compare less than or equal to the specified key or to the @link DeltaMap#end @endlink if such an elment could not be found.
 *
 * @remark The runtime is logarithmic in the size of the map.
 */

template <typename TConfig, typename TSpec, typename TPosition, typename TDeltaType>
inline Pair<typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type>
equalRange(DeltaMap<TConfig, TSpec> const & deltaMap, TPosition refPosition, TDeltaType /*deltaType*/)
{
    auto rBeg = lowerBound(deltaMap, refPosition, TDeltaType());
    auto rEnd = upperBound(deltaMap, refPosition, TDeltaType());
    return Pair<typename Iterator<DeltaMap<TConfig, TSpec> const, Standard>::Type>(rBeg, rEnd);
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
 *
 * @remark The map is implemented as a vector and the insertion time is linear in worst case.
 */

// TODO(rrahn): Change to return iterator as specified for insert of associative containers of the STL.
template <typename TConfig, typename TSpec, typename TDeltaPos, typename TDeltaValue,
          typename TCoverage, typename TDeltaType>
inline void
insert(DeltaMap<TConfig, TSpec> & deltaMap,
       TDeltaPos deltaPos,
       TDeltaValue const & value,
       TCoverage const & coverage,
       TDeltaType const & /*deltaType*/)
{
    if (SEQAN_UNLIKELY(empty(deltaMap)))
        reserve(deltaMap._entries, 1);

    auto it = upperBound(deltaMap, deltaPos, TDeltaType());
    impl::insert(it, deltaPos, value, coverage, TDeltaType());
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

// TODO(rrahn): Change to return iterator as specified for insert of associative containers of the STL.
// This is the key based method.
// Also implement the iterator base method.
// TODO(rrahn): Put position and delta type into a key type.
template <typename TConfig, typename TSpec, typename TDeltaPos, typename TDeltaType>
inline typename Size<DeltaMap<TConfig, TSpec> >::Type
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

    auto itRange = equalRange(deltaMap, deltaPos, TDeltaType());
    auto count = itRange.i2 - itRange.i1;
    if (itRange.i1 == end(deltaMap, Standard()))
        return count;  // If lower bound is at end, there is no element with the given key.

    for (auto it = itRange.i1; it != itRange.i2; ++it)
    {
        // Find potential right end of this delta.
        SEQAN_IF_CONSTEXPR (!IsSameType<TDeltaType, DeltaTypeIns>::VALUE)
        {
            auto itR = impl::find(deltaMap,
                                  deltaPos + deletionSize(deltaMap._deltaStore, getDeltaRecord(*it).i2, TDeltaType()) - 1,
                                  DeltaEndType::IS_RIGHT,
                                  TDeltaType());

            SEQAN_ASSERT_LEQ(position(it, deltaMap), position(itR, deltaMap));
            if (it != itR)
            {
                erase(deltaMap._entries, itR - begin(deltaMap, Standard()));
                ++count;  // Increase count number by removed entry.
            }
        }

        // Erase the delta record from the corresponding delta store.
        TStoreSize storePos = getDeltaRecord(*it).i2;
        SEQAN_ASSERT_LT(storePos, length(getDeltaStore(deltaMap._deltaStore, TDeltaType())));
        eraseDeltaValue(deltaMap._deltaStore, storePos, TDeltaType());
        // Update record position for all entries.
        forEach(deltaMap, [storePos](TEntry & entry)
        {
            if (getDeltaRecord(entry).i1 == selectDeltaType(TDeltaType()) && getDeltaRecord(entry).i2 > storePos)
                --getDeltaRecord(entry).i2;  // Decrease record position by one.
        });
    }

    // Erase all entries in the range of equal values.
    erase(deltaMap._entries, position(itRange.i1, deltaMap), position(itRange.i2, deltaMap));

    return count;
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

template <typename TConfig, typename TSpec, typename TIteratorSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Tag<TIteratorSpec> const>::Type
begin(DeltaMap<TConfig, TSpec> & deltaMap, Tag<TIteratorSpec> const /*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec>, Tag<TIteratorSpec> const>::Type(deltaMap, 0);
}

template <typename TConfig, typename TSpec, typename TIteratorSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Tag<TIteratorSpec> const>::Type
begin(DeltaMap<TConfig, TSpec> const & deltaMap, Tag<TIteratorSpec> const/*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec> const, Tag<TIteratorSpec> const>::Type(deltaMap, 0);
}

// ----------------------------------------------------------------------------
// Function iter()
// ----------------------------------------------------------------------------

template <typename TConfig, typename TSpec, typename TIteratorSpec, typename TPos>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Tag<TIteratorSpec> const>::Type
iter(DeltaMap<TConfig, TSpec> & deltaMap,
     TPos const & pos,
     Tag<TIteratorSpec> const /*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec>, Tag<TIteratorSpec> const>::Type(deltaMap, pos);
}

template <typename TConfig, typename TSpec, typename TIteratorSpec, typename TPos>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Tag<TIteratorSpec> const>::Type
iter(DeltaMap<TConfig, TSpec> const & deltaMap,
     TPos const & pos,
     Tag<TIteratorSpec> const &/*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec> const, Tag<TIteratorSpec> const>::Type(deltaMap, pos);
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

template <typename TConfig, typename TSpec, typename TIteratorSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec>, Tag<TIteratorSpec> const>::Type
end(DeltaMap<TConfig, TSpec> & deltaMap, Tag<TIteratorSpec> const /*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec>, Tag<TIteratorSpec> const>::Type(deltaMap, size(deltaMap));
}

template <typename TConfig, typename TSpec, typename TIteratorSpec>
inline typename Iterator<DeltaMap<TConfig, TSpec> const, Tag<TIteratorSpec> const>::Type
end(DeltaMap<TConfig, TSpec> const & deltaMap, Tag<TIteratorSpec> const /*tag*/)
{
    return typename Iterator<DeltaMap<TConfig, TSpec> const, Tag<TIteratorSpec> const>::Type(deltaMap, size(deltaMap));
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
inline auto
size(DeltaMap<TConfig, TSpec> const & deltaMap) -> decltype(length(deltaMap._entries))
{
    return length(deltaMap._entries);
}

// ----------------------------------------------------------------------------
// Function maxSize()
// ----------------------------------------------------------------------------

/*!
 * @fn DeltaMap#maxSize
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
constexpr typename Size<DeltaMap<TConfig, TSpec> >::Type
maxSize(DeltaMap<TConfig, TSpec> const & /*deltaMap*/)
{
    return std::numeric_limits<typename Size<DeltaMap<TConfig, TSpec> >::Type>::max();
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_DELTA_MAP_H_
