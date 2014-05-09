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
// Implements some utility functions for data parallel processing of journal
// strings.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_UTIL_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_UTIL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct MergePointSyncResize_;
typedef Tag<MergePointSyncResize_> MergePointSyncResize;

struct MergePointSyncReadOnly_;
typedef Tag<MergePointSyncReadOnly_> MergePointSyncReadOnly;

// ----------------------------------------------------------------------------
// Class MergePointMap_
// ----------------------------------------------------------------------------

template <typename TVariantMap>
class MergePointMap_
{
public:
    typedef typename Value<MergePointMap_>::Type TValue;
    typedef typename DeltaCoverage<TVariantMap>::Type TDeltaCoverage;
    typedef String<TValue> TMergePoints;

    TVariantMap*    _varMapPtr;
    mutable TDeltaCoverage _mergeCoverage;
    mutable TMergePoints    _mergePoints;  // Stores the reference position of the merge point

    MergePointMap_() : _varMapPtr(NULL), _mergeCoverage(), _mergePoints()
    {}

    MergePointMap_(TVariantMap & map) : _varMapPtr(&map), _mergeCoverage(), _mergePoints()
    {
        resize(_mergeCoverage, coverageSize(map), false, Exact());
    }

};

// ----------------------------------------------------------------------------
// Class MergePointComparator_
// ----------------------------------------------------------------------------

template <typename TMergePoint>
struct MergePointComparator_
{
    inline bool operator()(TMergePoint const & lhs, TMergePoint const & rhs)
    {
        return !(lhs.i1 < rhs.i1);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TVariantMap>
struct Value<MergePointMap_<TVariantMap> >
{
    typedef typename Position<TVariantMap>::Type TPosition_;
    typedef Pair<TPosition_, TPosition_> Type;
};

template <typename TVariantMap>
struct Value<MergePointMap_<TVariantMap> const>
{
    typedef typename Position<TVariantMap>::Type TPosition_;
    typedef Pair<TPosition_, TPosition_> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TVariantMap>
struct Reference<MergePointMap_<TVariantMap> >
{
    typedef typename Value<MergePointMap_<TVariantMap> >::Type & Type;
};

template <typename TVariantMap>
struct Reference<MergePointMap_<TVariantMap> const>
{
    typedef typename Value<MergePointMap_<TVariantMap> const>::Type & Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TVariantMap>
inline void clear(MergePointMap_<TVariantMap> & mergePointStore)
{
    clear(mergePointStore._mergePoints);
}

// ----------------------------------------------------------------------------
// Function push()
// ----------------------------------------------------------------------------

template <typename TVariantMap, typename TPos, typename TIter>
inline void push(MergePointMap_<TVariantMap> & mergePointStore,
                 TPos const & pos,
                 TIter const & branchNodeIt)
{
    typedef MergePointMap_<TVariantMap> TMergePointStore;
    typedef typename TMergePointStore::TMergePoints TMergePoints;
    typedef typename Iterator<TMergePoints, Rooted>::Type TIterator;
    typedef typename Value<TMergePointStore>::Type TValue;

    TValue tmp(pos, branchNodeIt - begin(*mergePointStore._varMapPtr, Standard()));
    if (empty(mergePointStore._mergePoints))
    {
        insertValue(mergePointStore._mergePoints, 0, tmp);
        return;
    }

    TIterator it = std::upper_bound(begin(mergePointStore._mergePoints, Rooted()),
                                    end(mergePointStore._mergePoints, Rooted()),
                                    tmp,
                                    MergePointComparator_<TValue>());
    insertValue(mergePointStore._mergePoints, position(it), tmp);
    // Record the new merge point character.
    mergePointStore._mergeCoverage |= deltaCoverage(branchNodeIt);
}

template <typename TVariantMap>
inline void pop(MergePointMap_<TVariantMap> & mergePointStore)
{
    eraseBack(mergePointStore._mergePoints);
}

template <typename TVariantMap>
inline typename Reference<MergePointMap_<TVariantMap> >::Type
topMergePoint(MergePointMap_<TVariantMap> & mergePointStore)
{
    return back(mergePointStore._mergePoints);
}

template <typename TVariantMap>
inline typename Reference<MergePointMap_<TVariantMap> const>::Type
topMergePoint(MergePointMap_<TVariantMap> const & mergePointStore)
{
    return back(mergePointStore._mergePoints);
}

template <typename TVariantMap>
inline typename DeltaCoverage<TVariantMap>::Type &
topMergeCoverage(MergePointMap_<TVariantMap> & mergePointStore)
{
    return deltaCoverage(iter(*mergePointStore._varMapPtr, back(mergePointStore._mergePoints).i2));
}

template <typename TVariantMap>
inline typename DeltaCoverage<TVariantMap>::Type const &
topMergeCoverage(MergePointMap_<TVariantMap> const & mergePointStore)
{
    return deltaCoverage(iter(*mergePointStore._varMapPtr, back(mergePointStore._mergePoints).i2));
}


template <typename TVariantMap, typename TPosition>
inline typename DeltaCoverage<TVariantMap>::Type &
getMergeCoverage(MergePointMap_<TVariantMap> & mergePointStore,
                 TPosition const & pos)
{
    return deltaCoverage(iter(*mergePointStore._varMapPtr, mergePointStore._mergePoints[pos].i2, Standard()));
}

template <typename TVariantMap, typename TPosition>
inline typename DeltaCoverage<TVariantMap>::Type const &
getMergeCoverage(MergePointMap_<TVariantMap> const & mergePointStore,
                 TPosition const & pos)
{
    return deltaCoverage(iter(*mergePointStore._varMapPtr, mergePointStore._mergePoints[pos].i2, Standard()));
}


template <typename TVariantMap, typename TPosition>
inline bool
_updateMergePoints(MergePointMap_<TVariantMap> & mergePointStack,
                   TPosition const & pos)
{
    typedef MergePointMap_<TVariantMap> TMergePointStack;
    typedef typename Value<TMergePointStack>::Type TMergePoint;
    typedef typename TMergePointStack::TMergePoints const TMergePoints;
    typedef typename Iterator<TMergePoints const, Standard>::Type TMergePointIt;

    if (length(mergePointStack._mergePoints) == 1u)
        return false;

    TMergePoint tmp(pos, 0);
    TMergePointIt itEnd = std::upper_bound(begin(mergePointStack._mergePoints, Standard()),
                                           end(mergePointStack._mergePoints, Standard()), tmp,
                                           MergePointComparator_<TMergePoint>());
    --itEnd;
    if (position(itEnd, mergePointStack._mergePoints) + 1 != length(mergePointStack._mergePoints))
    {
        erase(mergePointStack._mergePoints, position(itEnd, mergePointStack._mergePoints) + 1,
              length(mergePointStack._mergePoints));
        arrayFill(begin(mergePointStack._mergeCoverage, Standard()), end(mergePointStack._mergeCoverage, Standard()),
                  false);
        for (TMergePointIt it = begin(mergePointStack._mergePoints, Standard()) + 1;
             it != end(mergePointStack._mergePoints, Standard()); ++it)
            mergePointStack._mergeCoverage |= deltaCoverage(iter(*mergePointStack._varMapPtr, it->i2, Standard()));
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function _syncToMergePoint()
// ----------------------------------------------------------------------------

template <typename TCoverage, typename TMergePointStack, typename TPosition>
inline unsigned
_syncToMergePoint(TCoverage & target,
                  TMergePointStack const & mergePointStack,
                  TPosition const & pos,
                  MergePointSyncReadOnly const & /*tag*/)
{
    typedef typename Value<TMergePointStack>::Type TMergePoint;
    typedef typename TMergePointStack::TMergePoints const TMergePoints;
    typedef typename Iterator<TMergePoints const, Rooted>::Type TMergePointIt;

    if (length(mergePointStack._mergePoints) == 1u)
        return 1u;

    TMergePoint tmp(pos, 0);
    TMergePointIt itEnd = std::upper_bound(begin(mergePointStack._mergePoints, Rooted()),
                                           end(mergePointStack._mergePoints, Rooted()), tmp,
                                           MergePointComparator_<TMergePoint>());
    --itEnd;
    TMergePointIt it = end(mergePointStack._mergePoints, Rooted()) -1;

    SEQAN_ASSERT(mergePointStack._varMapPtr != NULL);  // Check if the pointer is initialized.
    for (; it != itEnd; --it)
        transform(target, target, mappedCoverage(*mergePointStack._varMapPtr, it->i2), FunctorBitwiseOr());

    return position(itEnd) + 1;
}

template <typename TCoverage, typename TMergePointStack, typename TPosition>
inline unsigned
_syncToMergePoint(TCoverage & target,
                  TMergePointStack & mergePointStack,
                  TPosition const & pos,
                  MergePointSyncResize const & /*tag*/)
{
    unsigned newLength = _syncToMergePoint(target, mergePointStack, pos, MergePointSyncReadOnly());
    SEQAN_ASSERT_LEQ(newLength, length(mergePointStack._mergePoints));

    if (newLength != length(mergePointStack._mergePoints))
        erase(mergePointStack._mergePoints, newLength, length(mergePointStack._mergePoints));
    return newLength;
}

// ----------------------------------------------------------------------------
// Function _mapVirtualToVirtual()
// ----------------------------------------------------------------------------

// Rebases the target iterator to the source iterator's position based on the common host.
template <typename TIter, typename TIterSrc, typename TBranchNodeIt, typename TDeltaMap, typename TProxyId>
inline void
_mapVirtualToVirtual(TIter & target,
                     TIterSrc const & source,
                     TBranchNodeIt const & branchNodeIt,
                     TDeltaMap const & variantStore,
                     TProxyId const & proxyId)
{
    typedef typename Position<TIter>::Type TPosition;
    // Check if both journals point to the same reference.
    SEQAN_ASSERT_EQ(&host(*target._journalStringPtr), &host(*source._journalStringPtr));

    if (source._journalEntriesIterator->segmentSource == SOURCE_PATCH)
    {
        SEQAN_ASSERT_GEQ(*branchNodeIt, _physicalOriginPosition(source));
        if (_physicalOriginPosition(source) == 0)
        {
            setPosition(target, *branchNodeIt + _localEntryPosition(source));
            return;
        }

        unsigned hostPos = _physicalOriginPosition(source) - 1;  // In the patch node -> the first position that is original to the left plus 1.
        unsigned mappedVirtPos = hostToVirtualPosition(*target._journalStringPtr, hostPos);
        setPosition(target, mappedVirtPos);

        SEQAN_ASSERT_EQ(target._journalEntriesIterator->segmentSource, SOURCE_ORIGINAL);  // TODO(rmaerker): Check end condition!

        if (_physicalOriginPosition(target) >= *branchNodeIt)
        {
            while(!atBegin(target._journalEntriesIterator, target._journalStringPtr->_journalEntries) &&
                (--target)._journalEntriesIterator->segmentSource == SOURCE_PATCH);

            if (atBegin(target._journalEntriesIterator, target._journalStringPtr->_journalEntries) &&
                target._journalEntriesIterator->segmentSource == SOURCE_PATCH)
            {
                setPosition(target, static_cast<TPosition>(0));
                hostPos = 0;
            }
            else
            {
                SEQAN_ASSERT_EQ(target._journalEntriesIterator->segmentSource, SOURCE_ORIGINAL);
                _updateSegmentIteratorsLeft(target);
                hostPos = _physicalOriginPosition(target);
            }
        }

        TBranchNodeIt tmpIt = const_cast<TBranchNodeIt&>(branchNodeIt);  // We are at this position.
        unsigned virtOffset = target._journalEntriesIterator->length - _localEntryPosition(target) - 1;
        while (!atBegin(tmpIt, variantStore) && *(--tmpIt) > hostPos)
        {
            if (deltaCoverage(tmpIt)[proxyId] != true)  // Irrelevant variant.
                continue;

            // If between hostPos and breakpoint are any other insertion or deletion, then we keep track of this virtual offset.
            if (deltaType(tmpIt) == DeltaType::DELTA_TYPE_INS)
                virtOffset += length(deltaIns(tmpIt));
            else if (deltaType(tmpIt) == DeltaType::DELTA_TYPE_SNP)
                ++virtOffset;
            else if (deltaType(tmpIt) == DeltaType::DELTA_TYPE_INDEL)
                virtOffset += length(deltaIndel(tmpIt).i2);
        }
        target += (1 + virtOffset + _localEntryPosition(source));
    }
    else
        setPosition(target, hostToVirtualPosition(*target._journalStringPtr, _physicalOriginPosition(source)));
}

// ----------------------------------------------------------------------------
// Function _mapHostToVirtual()
// ----------------------------------------------------------------------------

template <typename TIterator, typename TValue, typename THostSpec, typename TBuffSpec, typename TDeltaMap,
          typename TProxyId, typename THostPos>
inline void
_mapHostToVirtual(TIterator & resultIt,
                  String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > & js,
                  TDeltaMap & variantStore,
                  TProxyId const & proxyId,
                  THostPos const & hostPos)
{
    typedef String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > TJournalString;
    typedef typename JournalType<TJournalString>::Type TJournalEntries;
    typedef typename Value<TJournalEntries>::Type TCargo;
    typedef typename Iterator<TJournalEntries>::Type TEntriesIterator;
    typedef typename Position<TCargo>::Type TCargoPos;
    typedef typename Size<TCargo>::Type TCargoSize;
    typedef JournalEntryLtByPhysicalOriginPos<TCargoPos, TCargoSize> TComp;

    typedef typename Iterator<TDeltaMap, Standard>::Type TVarIterator;

    // We need to set the iterator to the correct position within the proxy sequence given the host pos.
    TJournalEntries & journalEntries = _journalEntries(js);

    if (empty(journalEntries._journalNodes))
    {
        resultIt = end(js);  // Put the iterator into a valid state.
        return;
    }

    resultIt = begin(js);

    TCargo refCargo;
    refCargo.physicalOriginPosition = hostPos;
    TEntriesIterator it = std::lower_bound(begin(journalEntries._journalNodes, Standard()),
                                           end(journalEntries._journalNodes, Standard()), refCargo, TComp());

    // This is now the first position whose var is equal or greater to the host pos.
    // Since this is either a position that is deleted
    // or a position after the insertion made -> Even for a SNP
    // We have to go backwards.
    if (it != begin(journalEntries, Standard()))
        --it;

    while (it != begin(journalEntries, Standard()) && it->segmentSource == SOURCE_PATCH)
        --it;

    if (it->segmentSource == SOURCE_PATCH)  // The iterator has to be at the beginning.
    {
        TVarIterator itVar = begin(variantStore, Standard());
        SEQAN_ASSERT_LEQ(*itVar, static_cast<unsigned const>(hostPos));

        unsigned virtualOffset = 0;
        // Now we move to the right until we find the node that we are looking for and reconstruct the offset of the virtual positions.
        while(*itVar != static_cast<unsigned const>(hostPos) && !atEnd(itVar, variantStore))
        {
            if (deltaCoverage(itVar)[proxyId] != true)  // irrelevant variant.
            {
                ++itVar;
                continue;
            }
//            TMappedDelta deltaKey = mappedDelta(variantStore, position(itVar));
            if (deltaType(itVar) == DeltaType::DELTA_TYPE_INS)
                virtualOffset += length(deltaIns(itVar));
            else if (deltaType(itVar) == DeltaType::DELTA_TYPE_SNP)
                ++virtualOffset;
            else if (deltaType(itVar) == DeltaType::DELTA_TYPE_INDEL)
                virtualOffset += length(deltaIndel(itVar).i2);
            ++itVar;
        }
        resultIt += virtualOffset;  // Set the iterator to the beginning of the variant.
        return;
    }

    SEQAN_ASSERT_EQ(it->segmentSource, SOURCE_ORIGINAL);

    // We assume that the operation begins here!
    resultIt._journalEntriesIterator = it;
    if (it->physicalOriginPosition + it->length > static_cast<unsigned const>(hostPos))
    {
        _updateSegmentIterators(resultIt);
        if (it->physicalOriginPosition < hostPos)
            resultIt += hostPos - it->physicalOriginPosition;
        return;
    }

    _updateSegmentIteratorsLeft(resultIt);  // Set the iterator to the end of the current original node.
    if (_physicalPosition(resultIt) + 1 == static_cast<unsigned const>(hostPos))
    {
        ++resultIt;
        return;
    }

    // TODO(rmaerker): Can remove the binary Search here!
    // Find the first node that is left or equal to the current physical position!
    TVarIterator itVar = std::upper_bound(begin(variantStore, Standard()),
                                          end(variantStore, Standard()),
                                          _physicalPosition(resultIt));

    SEQAN_ASSERT_LEQ(*itVar, static_cast<unsigned const>(hostPos));

    unsigned virtualOffset = 0;
    // Now we move to the right until we find the node that we are looking for and reconstruct the offset of the virtual positions.
    while(*itVar != static_cast<unsigned const>(hostPos) && !atEnd(itVar, variantStore))
    {
        if (deltaCoverage(itVar)[proxyId] != true)  // irrelevant variant.
        {
            ++itVar;
            continue;
        }

//        TMappedDelta deltaKey = mappedDelta(variantStore, position(itVar));
        if (deltaType(itVar) == DeltaType::DELTA_TYPE_INS)
            virtualOffset += length(deltaIns(itVar));
        else if (deltaType(itVar) == DeltaType::DELTA_TYPE_SNP)
            ++virtualOffset;
        else if (deltaType(itVar) == DeltaType::DELTA_TYPE_INDEL)
            virtualOffset += length(deltaIndel(itVar).i2);
        ++itVar;
    }
    resultIt += virtualOffset + 1;  // Set the iterator to the beginning of the variant.
}

// ----------------------------------------------------------------------------
// Function _testEqual()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostSpec>
inline bool
_testEqual(String<TValue, Packed<THostSpec> > const & lhs,
           String<TValue, Packed<THostSpec> > const & rhs)
{
    typedef String<TValue, Packed<THostSpec> > TPackedString;
    typedef typename Host<TPackedString>::Type TPackedHost;
    typedef typename Iterator<TPackedHost const, Standard>::Type TConstPackedHostIterator;
    typedef PackedTraits_<TPackedString> TPackedTraits;

    TConstPackedHostIterator itLOperand = begin(host(lhs), Standard());
    TConstPackedHostIterator itROperand = begin(host(rhs), Standard());
    TConstPackedHostIterator itEndLOperand = end(host(lhs), Standard()) - 1;
    TConstPackedHostIterator itEndROperand = end(host(rhs), Standard()) - 1;

    while(itLOperand != itEndLOperand && itROperand != itEndROperand)
    {
        if (*itLOperand != *itROperand)
            return false;
        ++itLOperand;
        ++itROperand;
    }

    return (*itLOperand >> (TPackedTraits::VALUES_PER_HOST_VALUE - (length(lhs) % TPackedTraits::VALUES_PER_HOST_VALUE)))
        == (*itROperand >> (TPackedTraits::VALUES_PER_HOST_VALUE - (length(rhs) % TPackedTraits::VALUES_PER_HOST_VALUE)));
}

}

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_UTIL_H_
