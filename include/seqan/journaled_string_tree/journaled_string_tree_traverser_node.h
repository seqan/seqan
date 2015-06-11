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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_NODE_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_NODE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


template <typename TJst>
class JstTraversalNode
{
public:
    typedef typename Source<TJst>::Type                                 TSequence;
    typedef typename Iterator<TSequence, Standard>::Type                TSeqIterator;
    typedef typename Position<TSequence>::Type                          TPosition;
    typedef typename Member<TJst, JstBufferMember>::Type                TJstBuffer;
    typedef typename Member<TJstBuffer, JstBufferExtensionMap>::Type    TExtMap;
    typedef typename Iterator<TExtMap, Standard>::Type                  TDeltaIterator;
    typedef typename Host<TJst>::Type                                   TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type                     TCoverage;
    typedef typename Size<TJst>::Type                                   TSize;

    TPosition           mappedSrcEndPos;
    TSize               remainingSize;
    TDeltaIterator      curDelta;  // current delta iterator.
    TDeltaIterator      nextDelta;
    TSeqIterator        begEdgeIt;    // Current iterator of this edge.
    TSeqIterator        curEdgeIt;    // Current iterator of this edge.
    TSeqIterator        endEdgeIt;    // End of this edge.
    TCoverage           coverage;  // Coverage of this node.
    bool                isBase;
    bool                fromBase;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TJstLhs, typename TJstRhs>
inline void
swap(JstTraversalNode<TJstLhs> & lhs,
     JstTraversalNode<TJstRhs> & rhs)
{
    std::swap(lhs.mappedSrcEndPos, rhs.mappedSrcEndPos);
    std::swap(lhs.remainingSize, rhs.remainingSize);
    std::swap(lhs.curDelta, rhs.curDelta);
    std::swap(lhs.nextDelta, rhs.nextDelta);
    std::swap(lhs.begEdgeIt, rhs.begEdgeIt);
    std::swap(lhs.curEdgeIt, rhs.curEdgeIt);
    std::swap(lhs.endEdgeIt, rhs.endEdgeIt);
    swap(lhs.coverage, rhs.coverage);
    std::swap(lhs.isBase, rhs.isBase);
    std::swap(lhs.fromBase, rhs.fromBase);
}

}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSER_NODE_H_
