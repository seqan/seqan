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
// Class JstSequenceBuffer
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TDirection = TraverseForward, typename TSpec = StaticBuffer>
class JstSequenceBuffer
{

};

template <typename TJournaledStringTree>
class JstSequenceBuffer<TJournaledStringTree, TraverseForward>
{
public:

    typedef typename Member<JstSequenceBuffer, JstSeqBufferJournaledSet>::Type  TJournaledSet;
    typedef typename Value<JournaledSet>::Type                                  TJournaledSeq;
    typedef typename Iterator<TJournaledSeq, Standard>::Type                    TJournaledSeqIt;

    typedef typename Container<TJournaledStringTree>::Type                      TDeltaMap;
    typedef typename Iterator<TDeltaMap, Standard>::Type                        TDeltaIterator;

    typedef typename Position<TJournaledSet>::Type                              TJSetPos;
    typedef String<TJSetPos>                                                    TStringPositions;

    TJournaledSeq           source;
    Range<TJournaledSeqIt>  sourceRange;    // Range to be traversed over the Journaled-String-Tree.
    Range<TJournaledSeqIt>  bufferedChunk;  // The chunk currently buffered.
    Range<TDeltaIterator>   deltaRange;
    TJournaledSet           journaledSet;   // The actual sequences over the current chunk.
    TStringPositions        startPositions; // The begin positions of the strings within the chunk.
    TStringPositions        endPositions;   // The end positions of the strings within the chunk.

    JstSequenceBuffer()
    {}

};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Member                               [JstSeqBufferJournaledSet]
// ----------------------------------------------------------------------------

template <typename TJournaledStringTree, typename TDirection>
struct Member<JstSequenceBuffer<TJournaledStringTree, TDirection>, JstSeqBufferJournaledSet>
{
    typedef typename Member<TJournaledStringTree, JstBaseSequenceMember>::Type  TBaseSeq_;
    typedef typename Value<TBaseSeq_>::Type                                     TBaseSeqVal_;
    typedef String<TBaseSeqVal_, Journaled<Alloc<>, SortedArray,Alloc<> > >     TJournaledString_;
    typedef StringSet<TJournaledString_, Owner<Journaled> >                     Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{
    
// TODO(rrahn): Implement streamTo()

}  // namespace impl

// ============================================================================
// Public Functions
// ============================================================================

// TODO(rrahn): Implement sync
// TODO(rrahn): Implement setBeginPosition()
// TODO(rrahn): Implement setEndPosition()

// TODO(rrahn): window();
// TODO(rrahn): sourceRange();

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_SEQUENCE_BUFFER_H_
