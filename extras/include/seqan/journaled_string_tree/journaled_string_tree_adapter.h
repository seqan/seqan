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

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSize, typename TAlphabet>
struct CompareType_
{
    typedef typename DeltaType::TValue TDeltaType;

    TDeltaType _deltaType;
    TSize      _del;
    String<TAlphabet> _ins;

    CompareType_(TDeltaType const & deltaType, TSize const & del) : _deltaType(deltaType), _del(del)
    {}

    CompareType_(TDeltaType const & deltaType, String<TAlphabet> const & ins) : _deltaType(deltaType), _del(0), _ins(ins)
    {}

    CompareType_(TDeltaType const & deltaType, TAlphabet const & snp) : _deltaType(deltaType), _del(0)
    {
        appendValue(_ins, snp);
    }

    CompareType_(TDeltaType const & deltaType, TSize const & del, String<TAlphabet> const & ins) : _deltaType(deltaType), _del(del), _ins(ins)
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

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _journalSnp()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TSnp>
inline void
_journalSnp(TTarget & target,
           TPos refPos,
           TSnp const & snp)
{
    typedef typename JournalType<TTarget>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries, Standard>::Type TEntriesIterator;
    typedef typename Position<TJournalEntries>::Type TEntryPos;

    TEntriesIterator entryIt = end(_journalEntries(target), Standard()) -1;
    SEQAN_ASSERT_EQ(entryIt->segmentSource, SOURCE_ORIGINAL);
    SEQAN_ASSERT_GEQ(refPos, entryIt->physicalOriginPosition);
    SEQAN_ASSERT_GT(entryIt->physicalOriginPosition + entryIt->length, refPos);

    TEntryPos virtPos = entryIt->virtualPosition + (refPos - entryIt->physicalOriginPosition);
    _doRecordErase(target._journalEntries, entryIt, virtPos, virtPos + 1);
    entryIt = end(_journalEntries(target), Standard()) -1;
    TEntryPos physPos = length(target._insertionBuffer);
    appendValue(target._insertionBuffer, snp);
    _doRecordInsertion(target._journalEntries, entryIt, virtPos, physPos, 1u);
}


// ----------------------------------------------------------------------------
// Function _journalSnp()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TCoverage, typename TSnp>
inline void
_journalSnp(TTarget & target,
           TPos refPos,
           TCoverage const & coverage,
           TSnp const & snp)
{
    typedef typename Value<TTarget>::Type TJournalSeq;
    typedef typename Iterator<TCoverage const>::Type TCoverageIterator;

    TCoverageIterator itBegin = begin(coverage);
    TCoverageIterator it = begin(coverage);
    TCoverageIterator itEnd = end(coverage);

    for (; it != itEnd; ++it)
    {
        if (getValue(it) == true)
        {
            TJournalSeq & journal = value(target, it - itBegin);
            assignValue(journal, hostToVirtualPosition(journal, refPos), snp);
        }
    }
}

// ----------------------------------------------------------------------------
// Function _journalDel()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TSize>
inline void
_journalDel(TTarget & target,
           TPos refPos,
           TSize const & delLength)
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
// Function _journalDel()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TCoverage, typename TSize>
inline void
_journalDel(TTarget & target,
           TPos refPos,
           TCoverage const & coverage,
           TSize const & delLength)
{
    typedef typename Value<TTarget>::Type TJournalSeq;
    typedef typename Iterator<TCoverage const>::Type TCoverageIterator;

    TCoverageIterator itBegin = begin(coverage);
    TCoverageIterator it = begin(coverage);
    TCoverageIterator itEnd = end(coverage);

    for (; it != itEnd; ++it)
    {
        if (getValue(it) == true)
        {
            TJournalSeq & journal = value(target, it - itBegin);
            unsigned virtPos = hostToVirtualPosition(journal, refPos);
            erase(journal, virtPos, virtPos + delLength);
        }
    }
}

// ----------------------------------------------------------------------------
// Function _journalIns()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TInsertion>
inline void
_journalIns(TTarget & target,
           TPos refPos,
           TInsertion const & insSeq)
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
// Function _journalIns()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TCoverage, typename TInsertion>
inline void
_journalIns(TTarget & target,
           TPos refPos,
           TCoverage const & coverage,
           TInsertion const & ins)
{
    typedef typename Value<TTarget>::Type TJournalSeq;
    typedef typename Iterator<TCoverage const>::Type TCoverageIterator;

    TCoverageIterator itBegin = begin(coverage);
    TCoverageIterator it = begin(coverage);
    TCoverageIterator itEnd = end(coverage);

    for (; it != itEnd; ++it)
    {
        if (getValue(it) == true)
        {
            TJournalSeq & journal = value(target, it - itBegin);
            insert(journal, hostToVirtualPosition(journal, refPos), ins);
        }
    }
}

// ----------------------------------------------------------------------------
// Function _journalIns()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TIndel>
inline void
_journalIndel(TTarget & target,
              TPos refPos,
              TIndel const & indel)
{
    typedef typename JournalType<TTarget>::Type TJournalEntries;
    typedef typename Iterator<TJournalEntries, Standard>::Type TEntryIterator;
    typedef typename Position<TJournalEntries>::Type TEntryPos;

    TEntryIterator entryIt = end(target._journalEntries, Standard()) - 1;
    SEQAN_ASSERT_EQ(entryIt->segmentSource, SOURCE_ORIGINAL);
    SEQAN_ASSERT_GEQ(refPos, entryIt->physicalOriginPosition);
    SEQAN_ASSERT_GT(entryIt->physicalOriginPosition + entryIt->length, refPos);

    target._length -= indel.i1;
    TEntryPos virtPos = entryIt->virtualPosition + (refPos - entryIt->physicalOriginPosition);
    _doRecordErase(target._journalEntries, entryIt, virtPos, virtPos + indel.i1);

    entryIt = end(target._journalEntries, Standard()) - 1;
    SEQAN_ASSERT_EQ(entryIt->segmentSource, SOURCE_ORIGINAL);
    SEQAN_ASSERT_EQ(refPos + indel.i1, entryIt->physicalOriginPosition);

    target._length += length(indel.i2);
    TEntryPos physPos = length(target._insertionBuffer);
    append(target._insertionBuffer, indel.i2);
    _doRecordInsertion(target._journalEntries, entryIt, virtPos, physPos, length(indel.i2));
}

// ----------------------------------------------------------------------------
// Function _journalIns()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TPos, typename TCoverage, typename TIndel>
inline void
_journalIndel(TTarget & target,
              TPos refPos,
              TCoverage const & coverage,
              TIndel const & ins)
{
    typedef typename Value<TTarget>::Type TJournalSeq;
    typedef typename Iterator<TCoverage const>::Type TCoverageIterator;

    TCoverageIterator itBegin = begin(coverage);
    TCoverageIterator it = begin(coverage);
    TCoverageIterator itEnd = end(coverage);

    for (; it != itEnd; ++it)
    {
        if (getValue(it) == true)
        {
            TJournalSeq & journal = value(target, it - itBegin);
            insert(journal, hostToVirtualPosition(journal, refPos), ins);
        }
    }
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

// ----------------------------------------------------------------------------
// Function _transformJournalCoverage()
// ----------------------------------------------------------------------------

template <typename TSize, typename TAlphabet, typename TSpec, typename TBitVector>
inline void
_transformJournalCoverage(DeltaMap<TSize, TAlphabet, TSpec> & variantMap,
                          TBitVector const & bitVec)
{
    typedef typename Iterator<TBitVector const, Standard>::Type TBitVectorIter;

    TBitVectorIter itBegin = begin(bitVec);
    TBitVectorIter itEnd = end(bitVec);

    // Variant two: We search the reference sequence and then all diffs.
    Splitter<TBitVectorIter> vecSplitter(itBegin, itEnd, Parallel());
    SEQAN_OMP_PRAGMA(parallel for)
    for (unsigned job = 0; job < length(vecSplitter); ++job)
    {
        // TODO(rmaerker): Remove debug code.
        //        printf("Thread %ld of %ld threads.\n", job, length(refSplitter));
        //        THostView tmpHostView;
        //        if (job == 0)
        //        {
        //            tmpHostView._begin = refSplitter[job];
        //            tmpHostView._end = refSplitter[job+1];
        //        }
        //        else
        //        {
        //            tmpHostView._begin = refSplitter[job] - overlapSize;
        //            tmpHostView._end = refSplitter[job+1];
        //        }
        for (TBitVectorIter parallelIt = vecSplitter[job]; parallelIt != vecSplitter[job + 1]; ++parallelIt)
            appendValue(value(deltaCoverageStore(variantMap)._coverageData, parallelIt - itBegin), *parallelIt);
    }
}

/*!
 * @fn adaptTo
 * @brief Transforms delta map into journal string set.
 *
 * @signature adaptTo(target, src, tag)
 */

template <typename TJournalSequence, typename TValue, typename TAlphabet, typename TSpec,
          typename TBlockBegin, typename TBlockEnd, typename TParallelTag>
inline void
adaptTo(StringSet<TJournalSequence, Owner<JournaledSet> > & journalSet,
        DeltaMap<TValue, TAlphabet, TSpec> & variantMap,
        TBlockBegin blockBegin,
        TBlockEnd blockEnd,
        Tag<TParallelTag> parallelTag = Serial())
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type TMapIterator;

    typedef typename GetDeltaCoverageStore_<TDeltaMap>::Type TDeltaCoverageStore;
    typedef typename Value<TDeltaCoverageStore>::Type TBitVec;
    typedef typename Iterator<TBitVec, Standard>::Type TBitVecIter;

    typedef StringSet<TJournalSequence, Owner<JournaledSet> > TJournalSet;
    typedef typename Iterator<TJournalSet, Standard>::Type TJournalSetIter;
    typedef typename Size<TJournalSet>::Type TSize;

    // Ensure the size of the journal set is sufficient.
    resize(journalSet, coverageSize(variantMap), Exact());

    // Temporary string keeping track of the last journaled operation.
    String<TSize> _lastVisitedNodes;
    resize(_lastVisitedNodes, length(journalSet), 0, Exact());

    Splitter<TJournalSetIter> jSetSplitter(begin(journalSet, Standard()), end(journalSet, Standard()), parallelTag);
    SEQAN_OMP_PRAGMA(parallel for)
    for (unsigned jobId = 0; jobId < length(jSetSplitter); ++jobId)
    {
//        printf("Thread: %i of %i\n", jobId, omp_get_num_threads());
        unsigned jobBegin = jSetSplitter[jobId] - begin(journalSet, Standard());
        unsigned jobEnd = jSetSplitter[jobId + 1] - begin(journalSet, Standard());

        // First we set the reference sequence for all journal strings.
        for (TJournalSetIter it = jSetSplitter[jobId]; it != jSetSplitter[jobId + 1]; ++it)
            setHost(*it, host(journalSet));

//        printf("Thread %i: jobBegin %i - jobEnd %i\n", jobId, jobBegin, jobEnd);

        TMapIterator itMapBegin = begin(variantMap, Standard());
        TMapIterator itMap = itMapBegin + blockBegin;
        TMapIterator itMapEnd = begin(variantMap, Standard()) + blockEnd;
//        TCoverageIterator it = begin(deltaCoverageStore(variantMap), Standard()) + blockBegin;

        for (; itMap != itMapEnd; ++itMap)
        {
//            TMappedDelta varKey = mappedDelta(variantMap, itMap - itMapBegin);

            TBitVecIter itVecBegin = begin(deltaCoverage(itMap), Standard());
            TBitVecIter itVec = itVecBegin + jobBegin;
            TBitVecIter itVecEnd = begin(deltaCoverage(itMap), Standard()) + jobEnd;
            for (;itVec != itVecEnd; ++itVec)
            {
//                printf("Thread %i: vecPos %i\n", jobId, (int) (itVec - itVecBegin));
                SEQAN_ASSERT_NOT(empty(host(journalSet[itVec - itVecBegin])));

                if (!(*itVec))
                    continue;

//                _lastVisitedNodes[itVec - itVecBegin] = itMap - itMapBegin;
                switch(deltaType(itMap))
                {
                    case DeltaType::DELTA_TYPE_SNP:
                        _journalSnp(journalSet[itVec - itVecBegin], *itMap, deltaSnp(itMap));
                        break;
                    case DeltaType::DELTA_TYPE_DEL:
                        _journalDel(journalSet[itVec - itVecBegin], *itMap, deltaDel(itMap));
                        break;
                    case DeltaType::DELTA_TYPE_INS:
                        _journalIns(journalSet[itVec - itVecBegin], *itMap, deltaIns(itMap));
                        break;
                        // TODO(rmaerker): Add Case for INDEL
                }
            }
        }

    }

    // Now we have journaled all variants in the current block.

//    {
//        unsigned jobBegin = 0; //jSetSplitter[jobId] - begin(journalSet, Standard());
//        unsigned jobEnd = length(journalSet);//jSetSplitter[jobId + 1] - begin(journalSet, Standard());
//        // Second we run through all deltas and journal them. In case of parallelization we do it separately for each block of journal strings separately.
//        TMapIterator itMapBegin = begin(variantMap, Standard());
//        TMapIterator itMap = itMapBegin;
//        TMapIterator itMapEnd = end(variantMap, Standard());
//        TCoverageIterator it = begin(deltaCoverageStore(variantMap), Standard());
//
//        for (; itMap != itMapEnd; ++it, ++itMap)
//        {
//            TMappedDelta varKey = mappedDelta(variantMap, itMap - itMapBegin);
//
//            TBitVecIter itVecBegin = begin(*it, Standard()) + jobBegin;
//            TBitVecIter itVec = itVecBegin;
//            TBitVecIter itVecEnd = begin(*it, Standard()) + jobEnd;
//            for (;itVec != itVecEnd; ++itVec)
//            {
//                SEQAN_ASSERT_NOT(empty(host(journalSet[itVec - itVecBegin])));
//                if (!(*itVec))
//                    continue;
//
//                switch(deltaType(varKey))
//                {
//                    case DeltaType::DELTA_TYPE_SNP:
//                        _journalSnp(journalSet[itVec - itVecBegin], *itMap, deltaSnp(variantMap, deltaPosition(varKey)));
//                        break;
//                    case DeltaType::DELTA_TYPE_DEL:
//                        _journalDel(journalSet[itVec - itVecBegin], *itMap, deltaDel(variantMap, deltaPosition(varKey)));
//                        break;
//                    case DeltaType::DELTA_TYPE_INS:
//                        _journalIns(journalSet[itVec - itVecBegin], *itMap, deltaIns(variantMap, deltaPosition(varKey)));
//                        break;
//                        // TODO(rmaerker): Add Case for INDEL
//                }
//            }
//        }
//    }

}

template <typename TValue, typename TAlphabet, typename TSpec, typename TJournalSequence>
inline void
adaptTo(DeltaMap<TValue, TAlphabet, TSpec> & deltaMap,
        StringSet<TJournalSequence, Owner<JournaledSet> > const & journalSet)
{
    typedef DeltaMap<TValue, TAlphabet, TSpec> TDeltaMap;

    typedef typename GetDeltaCoverageStore_<TDeltaMap>::Type TDeltaCoverageStore;
    typedef typename Value<TDeltaCoverageStore>::Type TBitString;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_SNP>::Type TDeltaAlphabet;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type TDeltaDel;

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
                TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos, TCompareType(DeltaType::DELTA_TYPE_INS, insBuff)), tmpValue);
                TMapInsertResult res = _map.insert(_mapInsValue);
                if (!res.second)
                    appendValue(res.first->second, seqId);
                continue;
            }

            if (it->physicalOriginPosition > lastPhysPos)  // We detected a deletion.
            {
                // What if we found multiple insertions in between.
                TSize delSize = it->physicalPosition - lastPhysPos;
                register unsigned j = 0;
                for (; (j < length(tmpInsertionEntries)) && (delSize > 0); ++j)
                {
                    TEntiresIterator tmpInsIt = tmpInsertionEntries[j];
                    register unsigned k = 0;
                    for (; (k < tmpInsIt->length) && (delSize > 0); ++k, --delSize, ++currPhysBeginPos)
                    {
                        TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos, TCompareType(DeltaType::DELTA_TYPE_SNP, journalSeq._insertionBuffer[tmpInsIt->physicalPosition + k])), tmpValue);
                        TMapInsertResult res = _map.insert(_mapInsValue);
                        if (!res.second)
                            appendValue(res.first->second, seqId);

                    }
                    if (k < tmpInsIt->length)  // Deletion was smaller than current insertion, so we need to add the remaining insertion characters.
                    {
                        TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos,
                                     TCompareType(DeltaType::DELTA_TYPE_INS,
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
                    TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos, TCompareType(DeltaType::DELTA_TYPE_DEL, delSize)), tmpValue);
                    TMapInsertResult res = _map.insert(_mapInsValue);
                    if (!res.second)
                        appendValue(res.first->second, seqId);
                }
                else if (j < length(tmpInsertionEntries))  // Case 2: Add remaining untouched insertions.
                {
                    TInsertionBuff insBuff;
                    for (unsigned k = 0; k < length(tmpInsertionEntries); ++k)
                        append(insBuff, _getInsertion(tmpInsertionEntries[k], journalSeq));
                    TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos,TCompareType(DeltaType::DELTA_TYPE_INS, insBuff)), tmpValue);
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
                TMapInsertType _mapInsValue = TMapInsertType(TMapKey_(currPhysBeginPos,TCompareType(DeltaType::DELTA_TYPE_INS, insBuff)), tmpValue);
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
        resize(cov, coverageSize(deltaMap), false, Exact());
        for (unsigned i = 0; i < length(mapIt->second); ++i)
            cov[mapIt->second[i]] = true;

        switch(mapIt->first.i2._deltaType)
        {
            case DeltaType::DELTA_TYPE_SNP:
                insert(deltaMap, mapIt->first.i1, mapIt->first.i2._ins[0], cov);
                break;
            case DeltaType::DELTA_TYPE_DEL:
                insert(deltaMap, mapIt->first.i1, mapIt->first.i2._del, cov);
                break;
            case DeltaType::DELTA_TYPE_INS:
                insert(deltaMap, mapIt->first.i1, mapIt->first.i2._ins, cov);
                break;
            // TODO(rmaerker): Add Case for INDEL
        }
    }

}

}

#endif // EXTRAS_INCLUDE_SEQAN_DATA_PARALLEL_DATA_PARALLEL_ADAPTER_H_
