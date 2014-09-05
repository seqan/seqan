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
// Implements methods to adapt a journaled set to a variant store and
// vice versa.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_ADAPTER_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_ADAPTER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <typename T> struct GetStringSet;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class CompareType_
// ----------------------------------------------------------------------------

// TODO(rrahn): Remove after changing adaptTo method for journaled string set to delta map.
template <typename TSize, typename TAlphabet>
struct CompareType_
{
//    typedef typename TValue TDeltaType;


    DeltaType _deltaType;
    TSize      _del;
    String<TAlphabet> _ins;

    CompareType_(DeltaType const & deltaType, TSize const & del) : _deltaType(deltaType), _del(del)
    {}

    CompareType_(DeltaType const & deltaType, String<TAlphabet> const & ins) : _deltaType(deltaType), _del(0), _ins(ins)
    {}

    CompareType_(DeltaType const & deltaType, TAlphabet const & snp) : _deltaType(deltaType), _del(0)
    {
        appendValue(_ins, snp);
    }

    CompareType_(DeltaType const & deltaType, TSize const & del, String<TAlphabet> const & ins) : _deltaType(deltaType), _del(del), _ins(ins)
    {}

    // operator=

    CompareType_ & operator=(CompareType_ const & other)
    {
        if (this != &other)
        {
            _deltaType = other._deltaType;
            _del = other._del;
            _ins = other._ins;
        }
        return *this;
    }
};

// ----------------------------------------------------------------------------
// Class MapKeyCompareLessFunctor_
// ----------------------------------------------------------------------------

struct MapKeyCompareLessFunctor_
{
    template <typename TValue>
    inline bool
    operator()(TValue const & lhs, TValue const & rhs) const
    {
        if (lhs.i1 != rhs.i1)
            return lhs.i1 < rhs.i1;
        if (lhs.i2._deltaType != rhs.i2._deltaType)
            return lhs.i2._deltaType < rhs.i2._deltaType;
        if (lhs.i2._del != rhs.i2._del)
            return lhs.i2._del < rhs.i2._del;
        return isLess(lhs.i2._ins, rhs.i2._ins);
    }
};

// ----------------------------------------------------------------------------
// Class JournalDeltaContext_
// ----------------------------------------------------------------------------

template <typename TJst, typename TIterator>
struct JournalDeltaContext_
{
    typedef typename Container<TJst>::Type TDeltaMap_;
    typedef typename Iterator<TDeltaMap_, Standard>::Type TMapIter_;
    
    TJst &          jst;
    TIterator       deltaIterator;
    TMapIter_       deltaIteratorBegin;
    unsigned        bitVecBegin;
    unsigned        bitVecEnd;

    JournalDeltaContext_(TJst & _jst,
                         TIterator & _iter,
                         unsigned _begin,
                         unsigned _end) :
        jst(_jst),
        deltaIterator(_iter),
        deltaIteratorBegin(),
        bitVecBegin(_begin),
        bitVecEnd(_end)
    {
        deltaIteratorBegin = begin(container(jst), Standard());
    }

    template <typename TTag>
    inline bool operator()(TTag const & tag)
    {
        typedef typename Container<TIterator>::Type TContainer;
        typedef typename DeltaCoverage<TContainer>::Type TBitVector;
        typedef typename Iterator<TBitVector, Standard>::Type TBitVecIter;

        if (isDeltaType(deltaType(deltaIterator), tag))
        {
            TBitVecIter itVecBegin = begin(deltaCoverage(deltaIterator), Standard());
            TBitVecIter itVec = itVecBegin + bitVecBegin;
            TBitVecIter itVecEnd = begin(deltaCoverage(deltaIterator), Standard()) + bitVecEnd;
            for (;itVec != itVecEnd; ++itVec)
            {
                SEQAN_ASSERT_NOT(empty(host(value(stringSet(jst), itVec - itVecBegin))));

                if (getValue(itVec))
                    journalDelta(jst._journalSet[itVec - itVecBegin], deltaPosition(deltaIterator), deltaValue(deltaIterator, tag), tag);
            }
            return true;
        }
        return false;
    }
};

// ----------------------------------------------------------------------------
// Class  BufferJournaledStringsEndHelper_
// ----------------------------------------------------------------------------

template <typename TJst, typename TContextSize>
struct BufferJournaledStringsEndHelper_
{
    typedef typename GetStringSet<TJst>::Type TStringSet;
    typedef typename Position<typename Host<TJst>::Type>::Type THostPos;
    typedef typename Size<TStringSet>::Type TSetSize;

    typedef typename Container<TJst>::Type TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;

    TJst &       jst;
    TContextSize contextSize;
    TSetSize     jsPos;
    THostPos     localContextEndPos;
    TMapIterator mapIter;

    BufferJournaledStringsEndHelper_(TJst & _jst, TContextSize _contextSize) : jst(_jst), contextSize(_contextSize)
    {}

    template <typename TTag>
    inline bool operator()(TTag const & deltaTypeTag)
    {
        if (isDeltaType(deltaType(mapIter), deltaTypeTag))
        {
            if (deltaCoverage(mapIter)[jsPos])
            {
                journalDelta(jst._journalSet[jsPos], deltaPosition(mapIter), deltaValue(mapIter, TTag()), TTag());
                if (IsSameType<TTag, DeltaTypeDel>::VALUE)
                    localContextEndPos += deltaValue(mapIter, DeltaTypeDel());
                if (IsSameType<TTag, DeltaTypeIns>::VALUE)
                    localContextEndPos -= _min(contextSize,
                                               static_cast<TContextSize>(length(deltaValue(mapIter, DeltaTypeIns()))));
                if (IsSameType<TTag, DeltaTypeSV>::VALUE)
                {
                    localContextEndPos += deltaValue(mapIter, DeltaTypeSV()).i1;
                    localContextEndPos -= _min(contextSize,
                                               static_cast<TContextSize>(length(deltaValue(mapIter, DeltaTypeSV()).i2)));
                }
            }
            return true;
        }
        return false;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function journalDelta()                                       [DeltaTypeSnp]
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TSnp>
inline void
journalDelta(TTarget & target,
             TPos refPos,
             TSnp const & snp,
             DeltaTypeSnp const & /*tag*/)
{
    typedef typename JournalType<TTarget>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries, Standard>::Type TEntriesIterator;
    typedef typename Position<TJournalEntries>::Type TEntryPos;

    TEntriesIterator entryIt = end(_journalEntries(target), Standard()) -1;
    SEQAN_ASSERT_EQ(entryIt->segmentSource, SOURCE_ORIGINAL);
    SEQAN_ASSERT_GEQ(refPos, entryIt->physicalOriginPosition);
    SEQAN_ASSERT_GT(entryIt->physicalOriginPosition + entryIt->length, refPos);

    TEntryPos virtPos = entryIt->virtualPosition + (refPos - entryIt->physicalOriginPosition);
    appendValue(target._insertionBuffer, snp);
    _doRecordInsertion(target._journalEntries, entryIt, virtPos, length(target._insertionBuffer) - 1, 1u);
    entryIt = end(_journalEntries(target), Standard()) -1;
    ++virtPos;
    _doRecordErase(target._journalEntries, entryIt, virtPos, virtPos + 1);
}

// ----------------------------------------------------------------------------
// Function journalDelta()                                       [DeltaTypeDel]
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TSize>
inline void
journalDelta(TTarget & target,
             TPos refPos,
             TSize const & delLength,
             DeltaTypeDel const & /*tag*/)
{
    typedef typename JournalType<TTarget>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries, Standard>::Type TEntryIterator;
    typedef typename Position<TJournalEntries>::Type TEntryPos;

    TEntryIterator entryIt = end(target._journalEntries, Standard()) - 1;
    SEQAN_ASSERT_EQ(entryIt->segmentSource, SOURCE_ORIGINAL);
    SEQAN_ASSERT_GEQ(refPos, entryIt->physicalOriginPosition);
    SEQAN_ASSERT_GT(entryIt->physicalOriginPosition + entryIt->length, refPos);

    target._length -= delLength;
    TEntryPos virtPos = entryIt->virtualPosition + (refPos - entryIt->physicalOriginPosition);
    _doRecordErase(target._journalEntries, entryIt, virtPos, virtPos + delLength);
    if (length(target._journalEntries) == 0)
        clear(target._insertionBuffer);
}

// ----------------------------------------------------------------------------
// Function journalDelta()                                       [DeltaTypeIns]
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TInsertion>
inline void
journalDelta(TTarget & target,
             TPos refPos,
             TInsertion const & insSeq,
             DeltaTypeIns const & /*tag*/)
{
    typedef typename JournalType<TTarget>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries, Standard>::Type TEntryIterator;
    typedef typename Position<TJournalEntries>::Type TEntryPos;

    TEntryIterator entryIt = end(target._journalEntries, Standard()) - 1;
    SEQAN_ASSERT_EQ(entryIt->segmentSource, SOURCE_ORIGINAL);
    SEQAN_ASSERT_GEQ(refPos, entryIt->physicalOriginPosition);
    SEQAN_ASSERT_GT(entryIt->physicalOriginPosition + entryIt->length, refPos);

    target._length += length(insSeq);
    TEntryPos virtPos = entryIt->virtualPosition + (refPos - entryIt->physicalOriginPosition);
    TEntryPos physPos = length(target._insertionBuffer);
    append(target._insertionBuffer, insSeq);
    _doRecordInsertion(target._journalEntries, entryIt, virtPos, physPos, length(insSeq));
}

// ----------------------------------------------------------------------------
// Function journalDelta()                                        [DeltaTypeSV]
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TSV>
inline void
journalDelta(TTarget & target,
             TPos refPos,
             TSV const & sv,
             DeltaTypeSV const & /*tag*/)
{
    typedef typename JournalType<TTarget>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries, Standard>::Type TEntryIterator;
    typedef typename Position<TJournalEntries>::Type TEntryPos;

    TEntryIterator entryIt = end(target._journalEntries, Standard()) - 1;
    SEQAN_ASSERT_EQ(entryIt->segmentSource, SOURCE_ORIGINAL);
    SEQAN_ASSERT_GEQ(refPos, entryIt->physicalOriginPosition);
    SEQAN_ASSERT_GT(entryIt->physicalOriginPosition + entryIt->length, refPos);

    target._length -= sv.i1;
    TEntryPos virtPos = entryIt->virtualPosition + (refPos - entryIt->physicalOriginPosition);
    _doRecordErase(target._journalEntries, entryIt, virtPos, virtPos + sv.i1);

    entryIt = end(target._journalEntries, Standard()) - 1;
    target._length += length(sv.i2);
    TEntryPos physPos = length(target._insertionBuffer);
    append(target._insertionBuffer, sv.i2);
    _doRecordInsertion(target._journalEntries, entryIt, virtPos, physPos, length(sv.i2));
}

// ----------------------------------------------------------------------------
// Function _getInsertion()
// ----------------------------------------------------------------------------

template <typename TJournalIt, typename TJournalSequence>
inline typename InsertionBuffer<TJournalSequence>::Type
_getInsertion(TJournalIt const & it, TJournalSequence const & journalSeq)
{
    return infix(journalSeq._insertionBuffer, it->physicalPosition, it->physicalPosition + it->length);
}

// TODO(rrahn): Consider revision. Yet there is no real scenario where this is actually needed.
template <typename TValue, typename TAlphabet, typename TSpec, typename TJournalSequence>
inline void
adaptTo(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
        StringSet<TJournalSequence, Owner<JournaledSet> > const & journalSet)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;

    typedef typename GetDeltaCoverageStore_<TDeltaMap>::Type TDeltaCoverageStore;
    typedef typename Value<TDeltaCoverageStore>::Type TBitString;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TDeltaAlphabet;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeDel>::Type TDeltaDel;

    typedef StringSet<TJournalSequence, Owner<JournaledSet> > const TJournalSet;
    typedef typename Position<TJournalSet>::Type TSetPosition;
    typedef typename Position<TJournalSequence const>::Type TPosition;
    typedef typename Size<TJournalSequence const>::Type TSize;
    typedef typename InsertionBuffer<TJournalSequence>::Type TInsertionBuff;
    typedef typename JournalType<TJournalSequence const>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries, Standard>::Type TEntiresIterator;

    typedef CompareType_<TDeltaDel, TDeltaAlphabet> TCompareType;
    typedef Pair<TPosition, TCompareType> TMapKey_;
    typedef String<TSetPosition> TMapValue_;
    typedef std::map<TMapKey_, TMapValue_, MapKeyCompareLessFunctor_> TMap;
    typedef typename TMap::value_type TMapInsertType;
    typedef typename TMap::iterator TMapIterator;
    typedef std::pair<TMapIterator,bool> TMapInsertResult;

    TMap _map;
    for (unsigned seqId = 0; seqId < length(journalSet); ++seqId)
    {
        TJournalSequence const & journalSeq = journalSet[seqId];
        TEntiresIterator it = begin(_journalEntries(journalSeq), Standard());
        TEntiresIterator itEnd = end(_journalEntries(journalSeq), Standard());
        TMapValue_ tmpValue;
        appendValue(tmpValue, seqId);

        TPosition lastPhysPos = 0;
        TPosition currPhysBeginPos = 0;
        for (; it != itEnd; ++it)
        {
            String<TEntiresIterator> tmpInsertionEntries;

            for (; it != itEnd && it->segmentSource == SOURCE_PATCH; ++it)
                appendValue(tmpInsertionEntries, it);

            if (it == itEnd)  // Up to the end only insertions detected.
            {
                TInsertionBuff insBuff;
                for (unsigned j = 0; j < length(tmpInsertionEntries); ++j)
                    append(insBuff, _getInsertion(tmpInsertionEntries[j], journalSeq));
                TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos, TCompareType(DELTA_TYPE_INS, insBuff)), tmpValue);
                TMapInsertResult res = _map.insert(_mapInsValue);
                if (!res.second)
                    appendValue(res.first->second, seqId);
                continue;
            }

            if (it->physicalOriginPosition > lastPhysPos)  // We detected a deletion.
            {
                // What if we found multiple insertions in between.
                TSize delSize = it->physicalPosition - lastPhysPos;
                unsigned j = 0;
                for (; (j < length(tmpInsertionEntries)) && (delSize > 0); ++j)
                {
                    TEntiresIterator tmpInsIt = tmpInsertionEntries[j];
                    unsigned k = 0;
                    for (; (k < tmpInsIt->length) && (delSize > 0); ++k, --delSize, ++currPhysBeginPos)
                    {
                        TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos, TCompareType(DELTA_TYPE_SNP, journalSeq._insertionBuffer[tmpInsIt->physicalPosition + k])), tmpValue);
                        TMapInsertResult res = _map.insert(_mapInsValue);
                        if (!res.second)
                            appendValue(res.first->second, seqId);

                    }
                    if (k < tmpInsIt->length)  // Deletion was smaller than current insertion, so we need to add the remaining insertion characters.
                    {
                        TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos,
                                     TCompareType(DELTA_TYPE_INS,
                                                  infix(journalSeq._insertionBuffer,
                                                        tmpInsIt->physicalPosition + k,
                                                        tmpInsIt->physicalPosition + tmpInsIt->length))), tmpValue);
                        TMapInsertResult res = _map.insert(_mapInsValue);
                        if (!res.second)
                            appendValue(res.first->second, seqId);
                    }

                }
                if (delSize > 0)  // Case 1: Add the trailing deletion.
                {
                    TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos, TCompareType(DELTA_TYPE_DEL, delSize)), tmpValue);
                    TMapInsertResult res = _map.insert(_mapInsValue);
                    if (!res.second)
                        appendValue(res.first->second, seqId);
                }
                else if (j < length(tmpInsertionEntries))  // Case 2: Add remaining untouched insertions.
                {
                    TInsertionBuff insBuff;
                    for (unsigned k = 0; k < length(tmpInsertionEntries); ++k)
                        append(insBuff, _getInsertion(tmpInsertionEntries[k], journalSeq));
                    TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos,TCompareType(DELTA_TYPE_INS, insBuff)), tmpValue);
                    TMapInsertResult res = _map.insert(_mapInsValue);
                    if (!res.second)  // The element already exists.
                        appendValue(res.first->second, seqId);  // Appends seq id to existing insertion.
                }
            }
            else if (!empty(tmpInsertionEntries))  // We have do deal with an insertion.
            {
                TInsertionBuff insBuff;
                for (unsigned j = 0; j < length(tmpInsertionEntries); ++j)
                    append(insBuff, _getInsertion(tmpInsertionEntries[j], journalSeq));
                TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos,TCompareType(DELTA_TYPE_INS, insBuff)), tmpValue);
                TMapInsertResult res = _map.insert(_mapInsValue);
                if (!res.second)
                    appendValue(res.first->second, seqId);

            }
            // Update the physical positions.
            lastPhysPos = it->physicalOriginPosition + it->length;
            currPhysBeginPos = lastPhysPos;
        }
    }

    //NOTE(rmaerker): We could adapt the delta map to use a binary balanced tree in the background.
    // Fill deltaMap.
    setCoverageSize(deltaMap, length(journalSet));
    TMapIterator mapIt = _map.begin();
    TMapIterator mapItEnd = _map.end();

    for (; mapIt != mapItEnd; ++mapIt)
    {
        // Transform coverage.
        TBitString cov;
        resize(cov, getCoverageSize(deltaMap), false, Exact());
        for (unsigned i = 0; i < length(mapIt->second); ++i)
            cov[mapIt->second[i]] = true;

        switch(mapIt->first.i2._deltaType)
        {
            case DELTA_TYPE_SNP:
                insert(deltaMap, mapIt->first.i1, mapIt->first.i2._ins[0], cov);
                break;
            case DELTA_TYPE_DEL:
                insert(deltaMap, mapIt->first.i1, mapIt->first.i2._del, cov);
                break;
            case DELTA_TYPE_INS:
                insert(deltaMap, mapIt->first.i1, mapIt->first.i2._ins, cov);
                break;
            // TODO(rmaerker): Add Case for INDEL
        }
    }
}

}

#endif // EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_ADAPTER_H_
