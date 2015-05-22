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
// Tag Direction Tags
// ----------------------------------------------------------------------------

struct TraverseForward_;
typedef Tag<TraverseForward_> TraverseForward;

struct TraverseBackward_;
typedef Tag<TraverseBackward_> TraverseBackward;

struct TraverseBidirectional_;
typedef Tag<TraverseBidirectional_> TraverseBidirectional;

// ----------------------------------------------------------------------------
// Tag Member Tags
// ----------------------------------------------------------------------------

struct JstSeqBufferJournaledSet_;
typedef Tag<JstSeqBufferJournaledSet_> JstSeqBufferJournaledSet;

// ----------------------------------------------------------------------------
// Tag Buffering Tags
// ----------------------------------------------------------------------------

struct StaticBuffer_;
typedef Tag<StaticBuffer_> StaticBuffer;

template <unsigned BLOCK_SIZE = 10000>
struct DynamicBuffer;

// ----------------------------------------------------------------------------
// Class JstBuffer
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TDirection>
class JstBuffer
{
public:

    // __ Types ___________________________________________________________________

    typedef typename RemoveConst<TJournaledStringTree>::Type                    THost;
    typedef typename Member<JstBuffer, JstSeqBufferJournaledSet>::Type          TJournaledSet;
    typedef typename Value<JournaledSet>::Type                                  TJournaledSeq;
    typedef typename Iterator<TJournaledSeq, Standard>::Type                    TJournaledSeqIt;
    typedef typename Size<TJournaledSeq>::Type                                  TSize;

    typedef typename Container<TJournaledStringTree>::Type                      TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type                        TDeltaIterator;

    typedef typename Position<TJournaledSet>::Type                              TJSetPos;
    typedef String<TJSetPos>                                                    TStringPositions;

    // __ Members _________________________________________________________________

    THost const *           hostPtr;        // Pointer to the underlying host.
    TSize                   chunkSize;


    Range<TJournaledSeqIt>  sourceRange;    // Range to be traversed over the Journaled-String-Tree.
    Range<TJournaledSeqIt>  bufferedChunk;  // The chunk currently buffered.
    Range<TDeltaIterator>   deltaRange;     // The range of the delta map currentyl buffered.

    TStringPositions        startPositions; // The begin positions of the strings within the chunk.
    // TStringPositions        endPositions;   // The end positions of the strings within the chunk.

    TJournaledSeq           sourceJournaled;  // The original source sequence represented as JournaledString.
    TJournaledSet           journaledSet;   // The actual sequences over the current chunk.

    // __ Constructor _____________________________________________________________

    JstBuffer() : hostPtr(nullptr)
    {}

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Member                               [JstSeqBufferJournaledSet]
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TDirection>
struct Member<JstBuffer<TJournaledStringTree, TDirection>, JstSeqBufferJournaledSet>
{
    typedef typename Source<TJournaledStringTree>::Type                         Type;
};

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TJst, typename TDir>
struct Host<JstBuffer<TJst, TDir> >
{
    typedef TJst Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{
    
// TODO(rrahn): Implement streamTo()

// ----------------------------------------------------------------------------
// Function impl::setNewHost()
// ----------------------------------------------------------------------------

template <typename TJst, typename TDir>
inline void
setNewHost(JstBuffer<TJst, TDir> & buffer,
           TJst const & newHost)
{
    buffer.hostPtr = &newHost;
    clear(buffer.startPositions);
    resize(buffer.startPositions, depth(newHost), 0, Exact());
    clear(buffer.journaledSet);
    resize(buffer.journaledSet);
    setHost(buffer.journaledSet, baseSequence(newHost));
    setHost(buffer.sourceJournaled, host(buffer.journaledSet));
    assignRange(buffer.sourceRange, begin(buffer.sourceJournaled, Standard()), end(buffer.sequenceJournaled, Standard()));
    buffer.bufferedChunk = buffer.sourceRange;
    assignRange(buffer.deltaRange, begin(container(newHost), Standard()), end(container(newHost), Standard()));
    buffer.chunkSize = -1;
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

// TODO(rrahn): setHost();

template <typename TJst, typename TDir>
inline void
setHost(JstBuffer<TJst, TDir> & buffer,
        TJst const & host)
{
    clear(buffer);
    impl::setNewHost(host);
}


template <typename TJst, typename TDir>
inline typename Host<JstBuffer<TJst, TDir> const>::Type &
host(JstBuffer<TJst, TDir> const & buffer)
{
    return *buffer.hostPtr;
}

// TODO(rrahn): clear();


}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_SEQUENCE_BUFFER_H_
