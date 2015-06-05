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
// Basic defintions and forwards used globally for this module.
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

// ----------------------------------------------------------------------------
// Class JstBuffer
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TDirection = ForwardTraversal>
class JstBuffer
{
public:

    // __ Types ___________________________________________________________________

    typedef typename Member<JstBuffer, JstBufferSetMember>::Type    TJournaledSet;
    typedef typename Source<TJournaledStringTree>::Type             TSource;
    typedef typename Iterator<TSource const, Rooted>::Type          TSourceIterator;

    typedef typename RemoveConst<TJournaledStringTree>::Type const  TConstJst;
    typedef typename Host<TConstJst>::Type                          TConstDeltaMap;
    typedef typename Iterator<TConstDeltaMap, Standard>::Type       TDeltaIterator;

    typedef typename Size<TSource>::Type                            TSourceSize;
    typedef String<TSourceSize>                                     TStringPositions;

    // __ Members _________________________________________________________________

    TConstJst *             _jstPtr;         // Pointer to the underlying host.
    bool                    _isSynchronized;

    TDeltaIterator          _deltaRangeBegin;
    TDeltaIterator          _deltaRangeEnd;

    TStringPositions        _startPositions; // The begin positions of the strings within the chunk.
    TJournaledSet           _journaledSet;   // The actual sequences over the current chunk.

    // __ Constructor _____________________________________________________________

    JstBuffer() : _jstPtr(nullptr)
    {}

    JstBuffer(TJournaledStringTree const & jst)
    {
        setSource(*this, jst);
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
struct Size<JstBuffer<TJournaledStringTree, TSpec> >
{
    typedef typename Source<TJournaledStringTree>::Type TSource_;
    typedef typename Size<TSource_>::Type               Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member                               [JstSeqBufferJournaledSet]
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TSpec>
struct Member<JstBuffer<TJournaledStringTree, TSpec>, JstBufferSetMember>
{
    typedef typename Source<TJournaledStringTree>::Type TSource_;
    typedef StringSet<TSource_, Owner<JournaledSet> >   Type;
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
    SEQAN_ASSERT_GT(entryIt->physicalOriginPosition + entryIt->length, refPos);

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
        impl::journalDelta(*mapIt, getDeltaPosition(*mapIt), deltaValue(mapIt, TTag()), TTag());
    }
};

template <typename TJst, typename TSpec>
inline void
create(JstBuffer<TJst, TSpec> & buffer)
{
    typedef typename Member<JstBuffer<TJst, TSpec>, JstBufferSetMember>::Type   TJSet;
    typedef typename Value<TJSet>::Type                                         TJString;
    typedef typename Iterator<TJSet, Standard>::Type                            TJSetIter;
    typedef typename Container<TJst>::Type                                      TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type                        TMapIter;

    SEQAN_ASSERT(buffer._jstPtr != nullptr);  // This must be not nullptr.
    SEQAN_ASSERT(empty(buffer._startPositions));

    clear(buffer._journaledSet);  // Guarentee empty set.

    // Initialize the journaled set.
    resize(buffer._journaledSet, length(buffer._startPositions), Exact());

    forEach(begin(buffer._journaledSet, Standard()), end(buffer._journaledSet, Standard()),
    [&buffer](TJSetIter & it)
    {
        setHost(*it, host(source(buffer._jstPtr)));
    }, Parallel());

    // Construct the journaled string context

    Splitter<TJSetIter> jSetSplitter(begin(buffer._journalSet, Standard()), end(buffer._journalSet, Standard()),
                                     Parallel());
    SEQAN_OMP_PRAGMA(parallel for)
    for (auto jobId = 0; jobId < static_cast<int>(length(jSetSplitter)); ++jobId)
    {
        impl::JournaledStringCreateFunctor<TMapIter, TJSetIter> f;
        f.mapIt = buffer._deltaRangeBegin;
        for (; f.mapIt != buffer._deltaRangeEnd; ++f.mapIt)
        {
            f.setIt = jSetSplitter[jobId];
            auto covIt = begin(deltaCoverage(*f.mapIt), Standard()) + (f.setIt - begin(buffer._journalSet, Standard()));
            for (; f.setIt != jSetSplitter[jobId + 1]; ++f.setIt, ++covIt)
            {
                DeltaTypeSelector deltaSelector;
                if (*covIt)  // If the current sequence covers the current delta.
                    applyOnDelta(f, getDeltaType(*f.mapIt), deltaSelector);
            }
        }
    }

    // Now we can take into account the end positions.
}

// ----------------------------------------------------------------------------
// Functor impl::NetSizeExtractor
// ----------------------------------------------------------------------------

template <typename TSignedSize, typename TMapIter>
struct NetSizeExtractor
{
    TMapIter    mapIt;
    TSignedSize val;

    template <typename TTag>
    inline void operator()(TTag const &)
    {
        val = netSize(container(mapIt)._deltaStore, getStorePosition(*mapIt), TTag());
    }
};

// ----------------------------------------------------------------------------
// Function impl::synchronize()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
synchronize(JstBuffer<TJst, TSpec> & buffer)
{
    typedef typename Size<JstBuffer<TJst, TSpec> >::Type                TSize;
    typedef typename Container<TJst>::Type                              TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type                TMapIter;
    typedef typename MakeSigned<typename Size<TDeltaMap>::Type>::Type   TSignedSize;

    SEQAN_ASSERT(buffer.sourceRange.begin < buffer.sourceRange.end);

    //
    impl::NetSizeExtractor<TSignedSize, TMapIter> f;
    f.mapIt = begin(*buffer._deltaMapPtr, Standard());
    resize(buffer._startPositions, length(getDeltaCoverage(*f.mapIt)), sourceBeginPosition(*buffer._mapPtr), Exact());

    if (position(buffer.sourceRange.begin) == 0)  // Special case: We also set the begin position to 0 even if there is an insertion at the 0th position.
        return;

    for (; f.mapIt != buffer.deltaRange.begin; ++f.mapIt)
    {
        DeltaTypeSelector selector;
        applyOnDelta(f, getDeltaType(*f.mapIt), selector);
        auto covBegin = begin(getDeltaCoverage(*f.mapIt), Standard());
        for (auto covIt = covBegin; covIt != end(getDeltaCoverage(*f.mapIt), Standard()); ++covIt)
        {
            if (*covIt)
            {
                if (SEQAN_UNLIKELY(static_cast<TSize>(std::abs(f.val)) > buffer._startPositions[covIt - covBegin] &&
                                   (f.val < 0)))
                    buffer._startPositions[covIt - covBegin] = 0;  // In case the entire prefix of this sequence is deleted.
                else
                    buffer._startPositions[covIt - covBegin] += f.val;
            }
        }
    }
    buffer._isSynchronized = true;
}

// ----------------------------------------------------------------------------
// Function impl::sourceBegin()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename JstBuffer<TJst, TSpec>::TSourceIterator
sourceBegin(JstBuffer<TJst, TSpec> & buffer)
{
    return sourceBegin(*buffer._jstPtr);
}

// ----------------------------------------------------------------------------
// Function impl::sourceEnd()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline typename JstBuffer<TJst, TSpec>::TSourceIterator
sourceEnd(JstBuffer<TJst, TSpec> & buffer)
{
    return sourceEnd(*buffer._jstPtr);
}

}  // namespace impl

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
setHost(JstBuffer<TJst, TSpec> & buffer,
        TJst & jst)
{
    if (SEQAN_UNLIKELY(empty(jst)))
        return;

    buffer._isSynchronized = false;
    buffer._jstPtr = &jst;
    buffer._deltaRangeBegin = impl::lowerBound(host(jst), sourceBeginPosition(jst));
    buffer._deltaRangeEnd   = impl::lowerBound(host(jst), sourceEndPosition(jst));
}

// ----------------------------------------------------------------------------
// Function create()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
create(JstBuffer<TJst, TSpec> & buffer)
{
    if (!isSynchronized(buffer))
        return false;
    return impl::create(buffer);
}

// ----------------------------------------------------------------------------
// Function sync()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
sync(JstBuffer<TJst, TSpec> & buffer)
{
    if (isSynchronized(buffer))
        return;
    impl::synchronize(buffer);
}

// ----------------------------------------------------------------------------
// Function sync()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
unsync(JstBuffer<TJst, TSpec> & buffer)
{
    buffer._isSynchronized = false;
}

// ----------------------------------------------------------------------------
// Function isSynchronized()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline bool
isSynchronized(JstBuffer<TJst, TSpec> const & buffer)
{
    return buffer._isSynchronized;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline void
clear(JstBuffer<TJst, TSpec> & buffer)
{
    buffer._jstPtr = nullptr;
    buffer._isSynchronized = false;
    clear(buffer._startPositions);
}

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

template <typename TJst, typename TSpec>
inline TJst const &
host(JstBuffer<TJst, TSpec> & buffer)
{
    return *buffer._jstPtr;
}

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_SEQUENCE_BUFFER_H_
