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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct OnHit_;
typedef Tag<OnHit_> OnHit;

template <typename TState, typename TStack = String<TState, Block<> > >
struct StackEventObserver  // implements ObserverConcept
{
    TState &    state;
    TStack      stack;

    ExternalStateStack(TState & obj) : state(obj)
    {
        appendValue(stack, state);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TState, typename TStack>
inline void update(StackEventObserver<TState, TStack> & observer,
                   PushEvent const & /*tag*/)
{
    appendValue(observer.stack, observer.state);
}

template <typename TState, typename TStack>
inline void update(StackEventObserver<TState, TStack> & observer,
                   PopEvent const & /*tag*/)
{
    eraseBack(observer.stack);
}

// ----------------------------------------------------------------------------
// Function find
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TObserver, typename TState>
inline void
find(TraverserImpl<TContainer, JstTraverserSpec<TSpec>, TObserver> & input,
     TState & state)
{
    typedef StackStateObserver<TState>                                              TInternalObserver;
    typedef TraverserImpl<TContainer, JstTraverserSpec<TSpec>, TInternalObserver>   TInternalTraverser;

    TInternalObserver observer(state);

    TInternalTraverser traverser(input);  // Original traverser.
    addObserver(traverser, observer);

    while (!atEnd(traverser))
    {
        Pair<TSize, bool> result = back(observer.stack)(stringContext(traverser));
        if (result.i2)
            notify(input, OnHit());
        goNext(traverser, result.i1);
    }
}

// my algorithm:
// Jst ...
// trav = traverser(jst, observer);
// trav = traverser(jst, void);
// Stack<TAlgo> stack;
// trav = find(trav, stack);

}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSAL_H_
