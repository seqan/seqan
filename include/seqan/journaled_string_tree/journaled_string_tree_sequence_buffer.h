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
// Jst traversal buffer used to construct sequence contents.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_SEQUENCE_BUFFER_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_SEQUENCE_BUFFER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Member Tags
// ----------------------------------------------------------------------------

struct JstBufferSetMember_;
typedef Tag<JstBufferSetMember_> JstBufferSetMember;

struct JstBufferDeltaMapMember_;
typedef Tag<JstBufferDeltaMapMember_> JstBufferDeltaMapMember;

// ----------------------------------------------------------------------------
// Class JstBuffer_
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TDirection = ForwardTraversal>
class JstBuffer_
{
public:

    // __ Types ___________________________________________________________________

    typedef typename Member<JstBuffer_, JstBufferSetMember>::Type           TJournaledSet;
    typedef typename Member<TJournaledStringTree, JstSourceMember>::Type    TSource;
    typedef typename Iterator<TSource, Rooted>::Type                        TSourceIterator;

    typedef typename Member<JstBuffer_, JstBufferDeltaMapMember>::Type      TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type                    TDeltaIterator;

    typedef typename Size<TSource>::Type                                    TSourceSize;
    typedef String<TSourceSize>                                             TStringPositions;

    // __ Members _________________________________________________________________

    bool                    _isSynchronized;

    TDeltaMap*              _deltaMapPtr;

    TSourceIterator         _sourceBegin;       // Iterator to the clipped begin position of the source.
    TSourceIterator         _sourceEnd;         // Iterator to the clipped end position of the source.

    TDeltaIterator          _deltaRangeBegin;
    TDeltaIterator          _deltaRangeEnd;

    TStringPositions        _startPositions; // The begin positions of the strings within the chunk.
    TJournaledSet           _journaledSet;   // The actual sequences over the current chunk.

    // __ Constructor _____________________________________________________________

    JstBuffer_()
    {}

    JstBuffer_(TDeltaMap & map) : _deltaMapPtr(&map)
    {}

    JstBuffer_(TDeltaMap & map, TSourceIterator const & srcBeg, TSourceIterator const & srcEnd) :
        _deltaMapPtr(&map),
        _sourceBegin(srcBeg),
        _sourceEnd(srcEnd)
    {
        sync(*this);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TSpec>
struct Size<JstBuffer_<TJournaledStringTree, TSpec> >
{
    typedef typename Source<TJournaledStringTree>::Type TSource_;
    typedef typename Size<TSource_>::Type               Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member                               [JstSeqBufferJournaledSet]
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TSpec>
struct Member<JstBuffer_<TJournaledStringTree, TSpec>, JstBufferSetMember>
{
    typedef typename Member<TJournaledStringTree, JstSourceMember>::Type TSource_;
    typedef StringSet<TSource_, Owner<JournaledSet> >   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member                                [JstBufferDeltaMapMember]
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
struct Member<JstBuffer_<TJst, TSpec>, JstBufferDeltaMapMember>
{
    typedef typename Member<TJst, JstDeltaMapMember>::Type Type;
};

template <typename TJst, typename TSpec>
struct Member<JstBuffer_<TJst, TSpec> const, JstBufferDeltaMapMember>
{
    typedef typename Member<TJst, JstDeltaMapMember>::Type const Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::journalDelta()                                 [DeltaTypeSnp]
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
// Function impl::journalDelta()                                 [DeltaTypeDel]
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
// Function impl::journalDelta()                                 [DeltaTypeIns]
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
    SEQAN_ASSERT_GEQ(entryIt->physicalOriginPosition + entryIt->length, refPos);  // Special case for the insertion where it is valid to insert behind the sequence.

    target._length += length(insSeq);
    TEntryPos virtPos = entryIt->virtualPosition + (refPos - entryIt->physicalOriginPosition);
    TEntryPos physPos = length(target._insertionBuffer);
    append(target._insertionBuffer, insSeq);
    _doRecordInsertion(target._journalEntries, entryIt, virtPos, physPos, length(insSeq));
}

// ----------------------------------------------------------------------------
// Function impl::journalDelta()                                  [DeltaTypeSV]
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
// Function impl::create()
// ----------------------------------------------------------------------------

template <typename TMapIter, typename TJSIter>
struct JournaledStringCreateFunctor
{
    TMapIter mapIt;
    TJSIter  setIt;

    template <typename TTag>
    inline void
    operator()(TTag const & /*deltaType*/)
    {
        if ((*mapIt).deltaTypeEnd != DeltaEndType::IS_RIGHT)
            impl::journalDelta(*setIt, getDeltaPosition(*mapIt), deltaValue(mapIt, TTag()), TTag());
    }
};

template <typename TJst, typename TSpec>
inline void
create(JstBuffer_<TJst, TSpec> & buffer)
{
    typedef typename Member<JstBuffer_<TJst, TSpec>, JstBufferSetMember>::Type      TJSet;
    typedef typename Value<TJSet>::Type                                             TJString;
    typedef typename Iterator<TJSet, Standard>::Type                                TJSetIter;
    typedef typename Member<JstBuffer_<TJst, TSpec>, JstBufferDeltaMapMember>::Type TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type                            TDeltaMapIter;

    SEQAN_ASSERT(!empty(buffer._startPositions));

    clear(buffer._journaledSet);  // Guarentee empty set.

    // Initialize the journaled set.
    resize(buffer._journaledSet, length(buffer._startPositions), Exact());

    forEach(buffer._journaledSet, [&buffer](TJString &jStr)
    {
        setHost(jStr, host(container(buffer._sourceBegin)));
    }, Parallel());

    // Construct the journaled string context

    Splitter<TJSetIter> jSetSplitter(begin(buffer._journaledSet, Standard()), end(buffer._journaledSet, Standard()),
                                     Parallel());
    SEQAN_OMP_PRAGMA(parallel for)
    for (auto jobId = 0; jobId < static_cast<int>(length(jSetSplitter)); ++jobId)
    {
        impl::JournaledStringCreateFunctor<TDeltaMapIter, TJSetIter> f;
        f.mapIt = buffer._deltaRangeBegin;
        for (; f.mapIt != buffer._deltaRangeEnd; ++f.mapIt)
        {
            f.setIt = jSetSplitter[jobId];
            if (SEQAN_UNLIKELY((*f.mapIt).deltaTypeEnd == DeltaEndType::IS_RIGHT))
                continue;

            auto covIt = begin(getDeltaCoverage(*f.mapIt), Standard()) +
                (f.setIt - begin(buffer._journaledSet, Standard()));
            for (; f.setIt != jSetSplitter[jobId + 1]; ++f.setIt, ++covIt)
            {
                DeltaTypeSelector deltaSelector;
                if (*covIt)  // If the current sequence covers the current delta.
                    applyOnDelta(f, getDeltaType(*f.mapIt), deltaSelector);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function impl::synchronize()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
synchronize(JstBuffer_<TJst, TSpec> & buffer)
{
    typedef typename Size<JstBuffer_<TJst, TSpec> >::Type                           TSize;

    SEQAN_ASSERT(buffer._sourceBegin < buffer._sourceEnd);

    // TODO(rrahn): Implement block wise computation.
    // Discover overlapping deletions at begin and end.
    // 0123456789012345678901234...
    // xxxxx------------xxxxxxxx...
    //           [------xxxxxxxx...[  // Beginning of the buffer lies within a deletion of a sequence.
    //
    // When detecting the deletion check if end point lies behind source begin position.
    // If false: ignore deletion.
    // If true: Add deletion point to the mapExtension -> Maybe we can even add it as a new branching point.
        // In this case we add a deletion right from the beginning.
        // We just communicate through this object anyway. -> Then the algorithm does not have to change.
    // If this is at the end: we cut the deletion -> But we have one at the same position which is just longer.
    // But you would cut after the deletion anyway.
    // We could only represent infixes over the journaled strings after they have been constructed?
    // So we have a fixed end position and begin position.

    buffer._deltaRangeBegin = begin(*buffer._deltaMapPtr, Standard());
    buffer._deltaRangeEnd = end(*buffer._deltaMapPtr, Standard());

// TODO(rrahn): Enable this when allowing chunking.
//    buffer._deltaRangeBegin = std::lower_bound(begin(buffer._deltaMap, Standard()),
//                                               end(buffer._deltaMap, Standard()),
//                                               position(buffer._sourceBegin), DeltaExtensionCompareLessPos_());
//    buffer._deltaRangeEnd = std::lower_bound(begin(buffer._deltaMap, Standard()), end(buffer._deltaMap, Standard()),
//                                             position(buffer._sourceEnd), DeltaExtensionCompareLessPos_());

    // Stream from the beginning to the expected range to get the begin positions of the current segment.
    auto mapIt = begin(*buffer._deltaMapPtr, Standard());
    resize(buffer._startPositions, length(getDeltaCoverage(*mapIt)), position(buffer._sourceBegin), Exact());

    if (position(buffer._sourceBegin) == 0)  // Special case: We also set the begin position to 0 even if there is an insertion at the 0th position.
        return true;

    for (; mapIt != buffer._deltaRangeBegin; ++mapIt)
    {
        auto net = netSize(mapIt);
        auto covBegin = begin(getDeltaCoverage(*mapIt), Standard());
        for (auto covIt = covBegin; covIt != end(getDeltaCoverage(*mapIt), Standard()); ++covIt)
        {
            if (*covIt)
            {
                if (SEQAN_UNLIKELY(static_cast<TSize>(std::abs(net)) > buffer._startPositions[covIt - covBegin] &&
                                   (net < 0)))
                    buffer._startPositions[covIt - covBegin] = 0;  // In case the entire prefix of this sequence is deleted.
                else
                    buffer._startPositions[covIt - covBegin] += net;
            }
        }
    }
    buffer._isSynchronized = true;
    return true;
}

}  // namespace impl

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function sourceBegin()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename JstBuffer_<TJst, TSpec>::TSourceIterator
sourceBegin(JstBuffer_<TJst, TSpec> & buffer)
{
    return buffer._sourceBegin;
}

// ----------------------------------------------------------------------------
// Function setSourceBegin()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TSourceIter>
inline void
setSourceBegin(JstBuffer_<TJst, TSpec> & buffer,
               TSourceIter const & begin)
{
    buffer._isSynchronized = false;
    buffer._sourceBegin = begin;
}

// ----------------------------------------------------------------------------
// Function sourceEnd()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename JstBuffer_<TJst, TSpec>::TSourceIterator
sourceEnd(JstBuffer_<TJst, TSpec> & buffer)
{
    return buffer._sourceEnd;
}

// ----------------------------------------------------------------------------
// Function setSourceEnd()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec, typename TSourceIter>
inline void
setSourceEnd(JstBuffer_<TJst, TSpec> & buffer,
             TSourceIter const & end)
{
    buffer._isSynchronized = false;
    buffer._sourceEnd = end;
}

// ----------------------------------------------------------------------------
// Function setDeltaMap()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
setDeltaMap(JstBuffer_<TJst, TSpec> & buffer,
            typename Member<JstBuffer_<TJst, TSpec>, JstBufferDeltaMapMember>::Type & map)
{
    if (buffer._deltaMapPtr != &map)
    {
        markModified(buffer);
        buffer._deltaMapPtr = &map;
    }
}

// ----------------------------------------------------------------------------
// Function sync()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
sync(JstBuffer_<TJst, TSpec> & buffer)
{
    if (isSynchronized(buffer))
        return true;
    return impl::synchronize(buffer);
}

// ----------------------------------------------------------------------------
// Function create()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
create(JstBuffer_<TJst, TSpec> & buffer)
{
    if (!isSynchronized(buffer))
        if (!sync(buffer))
            return false;
    impl::create(buffer);
    return true;
}

// ----------------------------------------------------------------------------
// Function isSynchronized()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
isSynchronized(JstBuffer_<TJst, TSpec> const & buffer)
{
    return buffer._isSynchronized;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
clear(JstBuffer_<TJst, TSpec> & buffer)
{
    buffer._isSynchronized = false;
    buffer._deltaMapPtr = nullptr;
    clear(buffer._startPositions);
}

// ----------------------------------------------------------------------------
// Function markModified()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
markModified(JstBuffer_<TJst, TSpec> & buffer)
{
    buffer._isSynchronized = false;
}

// ----------------------------------------------------------------------------
// Function init();
// ----------------------------------------------------------------------------

// TODO(rrahn): Fix constness issue.
template <typename TJst, typename TSpec>
inline void
init(JstBuffer_<TJst, TSpec> & me,
     TJst & jst)
{
    setDeltaMap(me, impl::member(jst, JstDeltaMapMember()));
    setSourceBegin(me, begin(impl::member(jst, JstSourceMember()), Standard()));
    setSourceEnd(me, end(impl::member(jst, JstSourceMember()), Standard()));
    sync(me);
    create(me);
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_SEQUENCE_BUFFER_H_
