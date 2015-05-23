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
// Tag Traversal Tags
// ----------------------------------------------------------------------------

struct ForwardTraversal_;
typedef Tag<ForwardTraversal_> ForwardTraversal;

//struct BackwardTraversal_;
//typedef Tag<BackwardTraversal_> BackwardTraversal;


// ----------------------------------------------------------------------------
// Tag Buffer Tags
// ----------------------------------------------------------------------------

struct StaticBuffer_;
typedef Tag<StaticBuffer_> StaticBuffer;

template <unsigned BUFFER_SIZE = 10000u>
struct DynamicBuffer;


// ----------------------------------------------------------------------------
// Class JstBuffer
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TDirection = ForwardTraversal, typename TBufferSpec = StaticBuffer>
class JstBuffer
{
public:

    // __ Types ___________________________________________________________________

    typedef typename Member<JstBuffer, JstBufferSetMember>::Type    TJournaledSet;
    typedef typename Value<JournaledSet>::Type                      TJournaledSeq;
    typedef typename Iterator<TJournaledSeq, Standard>::Type        TJournaledSeqIt;
    typedef typename Size<TJournaledSeq>::Type                      TSize;

    typedef typename Container<typename RemoveConst<TJournaledStringTree>::Type >::Type  TDeltaMap;
    typedef typename Iterator<TDeltaMap const, Standard>::Type      TDeltaIterator;

    typedef typename Position<TJournaledSet>::Type                  TJSetPos;
    typedef String<TJSetPos>                                        TStringPositions;

    // __ Members _________________________________________________________________

    TDeltaMap const *       _deltaMapPtr;
    TSize                   _chunkSize;      // Size of the chunk.
    bool                    _isSynchronized;

    Range<TDeltaIterator>   deltaRange;     // The range over the deltas that need to be buffered.
    Range<TJournaledSeqIt>  sourceRange;    // Range represented by the Journaled-String-Tree.

    TJournaledSeqIt         _chunkBegin;
    TJournaledSeqIt         _chunkEnd;

    TStringPositions        _startPositions; // The begin positions of the strings within the chunk.
    TJournaledSet           _journaledSet;   // The actual sequences over the current chunk.

    // __ Constructor _____________________________________________________________

    JstBuffer() : _deltaMapPtr(nullptr)
    {}

    JstBuffer(TDeltaMap const & deltaMap) : _deltaMapPtr(&deltaMap)
    {}

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Member                               [JstSeqBufferJournaledSet]
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TDirection>
struct Member<JstBuffer<TJournaledStringTree, TDirection>, JstBufferSetMember>
{
    typedef typename Source<TJournaledStringTree>::Type                         Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{
    
// ----------------------------------------------------------------------------
// Function impl::getChunkSize()
// ----------------------------------------------------------------------------

template <typename TJst, typename TDir, typename TSpec>
inline constexpr typename Size<JstBuffer<TJst, TDir, TSpec> >::Type
getChunkSize(JstBuffer<TJst, TDir, TSpec> const & /*buffer*/)
{
    return MaxValue<typename Size<JstBuffer<TJst, TDir, TSpec> >::Type>::VALUE;
}

template <typename TJst, typename TDir, unsigned BUFFER_SIZE>
inline constexpr typename Size<JstBuffer<TJst, TDir, DynamicBuffer<BUFFER_SIZE> > >::Type
getChunkSize(JstBuffer<TJst, TDir, DynamicBuffer<BUFFER_SIZE> > const & /*buffer*/)
{
    return BUFFER_SIZE;
}

template <typename TSignedSize, typename TDeltaStore>
struct NetSizeExtractor
{
    TDeltaStore &   _store;
    TSignedSize     val;

    NetSizeExtractor(TDeltaStore & store) : _store(store), val(0)
    {}

    template <typename TEntry, typename TTag>
    inline void operator()(TEntry const & entry, TTag const &)
    {
        val = netSize(_store, getStorePos(entry), TTag());
    }
};

// ----------------------------------------------------------------------------
// Function impl::synchronize()
// ----------------------------------------------------------------------------

template <typename TJst, typename TDir, typename TSpec>
inline void
synchronize(JstBuffer<TJst, TDir, TSpec> & buffer)
{
    typedef typename Size<JstBuffer<TJst, TDir, TSpec> >::Type          TSize;
    typedef typename Container<TJst>::Type                              TDeltaMap;
    typedef typename Member<TDeltaMap, DeltaMapStoreMember>::Type       TDeltaStore;
    typedef typename MakeSigned<typename Size<TDeltaStore>::Type>::Type TSignedSize;

    SEQAN_ASSERT(buffer.sourceRange.begin < buffer.sourceRange.end);

    buffer._chunkSize = std::min(impl::getChunkSize(buffer), static_cast<TSize>(buffer.sourceRange.end - buffer.sourceRange.begin));
    // Now we need to stream to the position.

    auto deltaIt = begin(*buffer._deltaMapPtr, Standard());
    resize(buffer._startPositions, length(getDeltaCoverage(*deltaIt)), position(buffer.sourceRange.begin), Exact());

    if (position(buffer.sourceRange.begin) == 0)  // Special case: We also set the begin position to 0 even if there is an insertion at the 0th position.
        return;

    impl::NetSizeExtractor<TSignedSize, TDeltaStore const> f(buffer._deltaMapPtr->_deltaStore);
    for (; deltaIt != buffer.deltaRange.begin; ++deltaIt)
    {
        DeltaTypeSelector selector;
        applyOnDelta(f, *deltaIt, selector);
        auto covBegin = begin(getDeltaCoverage(*deltaIt), Standard());
        for (auto covIt = covBegin; covIt != end(getDeltaCoverage(*deltaIt), Standard()); ++covIt)
        {
            if (*covIt)
            {
                if (SEQAN_UNLIKELY(std::abs(f.val) > buffer._startPositions[covIt - covBegin]))
                    buffer._startPositions[covIt - covBegin] = 0;
                else
                    buffer._startPositions[covIt - covBegin] += f.val;
            }
        }
    }
    buffer._isSynchronized = true;
}

}  // namespace impl

// ============================================================================
// Public Functions
// ============================================================================

// TODO(rrahn): sourceRange();
// TODO(rrahn): Implement setSourceBegin() -> Stream to correct buffer position.
// TODO(rrahn): Implement setSourceEnd()  -> Stream to correct buffer position.
// TODO(rrahn): Implement setSourceRange()    -> stream to correct buffer position.

// TODO(rrahn): chunk() -> Returns buffered journaled Set -> Read-Only.
// TODO(rrahn): setChunkSize();
// TODO(rrahn): chunkSize();

// TODO(rrahn): advanceChunk();

// ----------------------------------------------------------------------------
// Function setSourceRange()
// ----------------------------------------------------------------------------

template <typename TJst, typename TDir, typename TSpec>
inline void
setDeltaMap(JstBuffer<TJst, TDir, TSpec> & buffer,
            typename Container<TJst>::Type const & container)
{
    clear(buffer);
    buffer._deltaMapPtr = &container;
    if (empty(container))
        return;
    buffer.deltaRange.begin = begin(*buffer._deltaMapPtr, Standard());
    buffer.deltaRange.end   = end(*buffer._deltaMapPtr, Standard());
}

// ----------------------------------------------------------------------------
// Function setSourceRange()
// ----------------------------------------------------------------------------

template <typename TJst, typename TDir, typename TSpec,
          typename TIter>
inline void
setSourceRange(JstBuffer<TJst, TDir, TSpec> & buffer, TIter const & srcBegin, TIter const & srcEnd)
{
    SEQAN_ASSERT(buffer._deltaMapPtr != nullptr);
    SEQAN_ASSERT(!empty(*buffer._deltaMapPtr));

    buffer._isSynchronized = false;

    buffer.deltaRange.begin = impl::lowerBound(*buffer._deltaMapPtr, position(srcBegin));
    buffer.deltaRange.end   = impl::lowerBound(*buffer._deltaMapPtr, position(srcEnd));
}

// ----------------------------------------------------------------------------
// Function setSourceRange()
// ----------------------------------------------------------------------------

template <typename TJst, typename TDir, typename TSpec>
inline void
sync(JstBuffer<TJst, TDir, TSpec> & buffer)
{
    if (buffer._isSynchronized)
        return;
    impl::synchronize(buffer);
}

template <typename TJst, typename TDir, typename TSpec>
inline void
clear(JstBuffer<TJst, TDir, TSpec> & buffer)
{
    buffer._deltaMapPtr = nullptr;
    buffer._isSynchronized = false;
}


}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_SEQUENCE_BUFFER_H_
