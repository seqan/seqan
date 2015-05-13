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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSAL_ENTRY_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSAL_ENTRY_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BaseJstTraversalEntry
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSequence>
class BaseJstTraversalEntry
{
public:

    typedef typename Iterator<TDeltaMap, Standard>::Type    TDeltaIter;
    typedef typename Iterator<TSequence, Standard>::Type    TSeqIterator;
    typedef typename Position<TSeqIterator>::Type           TSeqPosition;
    typedef typename DeltaCoverage<TDeltaMap>::Type         TCoverage;

    // ____________ Delta State ________________________

    TSeqPosition    bpNextVirtual; // Virtual position of next branch point within current journaled string.
    TDeltaIter      bp;        // Initial branch point from the source sequence. Only needed within branch.
    TDeltaIter      bpNext;    // Next branch point.

    // ____________ Sequence State  ____________________

    TSeqIterator    cur;  // Current iterator position.
    TSeqIterator    end;  // End of branch/source sequence.
    TSeqIterator    begBp;  // Iterator to beginning of the branchPoint.

    // ____________ Coverage ___________________________

    TCoverage coverage;             // The current coverage.
};

// ----------------------------------------------------------------------------
// Class JstTraversalEntry
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSequence, typename TState = Nothing>
class JstTraversalEntry : public BaseJstTraversalEntry<TDeltaMap, TSequence>
{
public:
    
    TState state;           // External state set by the algorithm.
    
    JstTraversalEntry()
    {}
    
    // TODO(rrahn): Check initializer list for construction.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Private Functions
// ============================================================================

// ============================================================================
// Public Functions
// ============================================================================

template <typename TDeltaMap, typename TSequence>
inline void
swap(BaseJstTraversalEntry<TDeltaMap, TSequence> & lhs,
     BaseJstTraversalEntry<TDeltaMap, TSequence> & rhs)
{
    swap(lhs.bpNextVirtual, rhs.bpNextVirtual);
    swap(lhs.bp, rhs.bp);
    swap(lhs.bpNext, rhs.bpNext);
    swap(lhs.cur, rhs.cur);
    swap(lhs.end, rhs.end);
    swap(lhs.begBp, rhs.begBp);
    swap(lhs.coverage, rhs.coverage);
}

template <typename TDeltaMap, typename TSequence, typename TState>
inline void
swap(JstTraversalEntry<TDeltaMap, TSequence, TState> & lhs,
     JstTraversalEntry<TDeltaMap, TSequence, TState> & rhs)
{
    typedef BaseJstTraversalEntry<TDeltaMap, TSequence> TBase;
    swap(static_cast<TBase&>(lhs), static_cast<TBase&>(rhs));
    swap(lhs.state, rhs.state);
}

// TODO (rmaerker): Move to separate header containing stl stack adaptions.

// ----------------------------------------------------------------------------
// Function pop()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline void pop(TContainer & stack)  // calls either pop_back or resize(length(cont) - 1);
{
    eraseBack(stack);
}

// ----------------------------------------------------------------------------
// Function top()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline typename Reference<TContainer>::Type
top(TContainer & stack)
{
    return back(stack);
}

template <typename TContainer>
inline typename Reference<TContainer const>::Type
top(TContainer const & stack)
{
    return back(stack);
}

// ----------------------------------------------------------------------------
// Function push()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline void
push(TContainer & stack,
     typename Value<TContainer>::Type SEQAN_FORWARD_CARG entry)
{
    typedef typename Value<TContainer>::Type TValue;
    appendValue(stack, SEQAN_FORWARD(TValue, entry));  // Enable perfect forwarding to allow move semantics of rvalue references.
}

// TODO (rmaerker): Implement size()
// TODO (rmaerker): Implement swap()

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSAL_ENTRY_H_
