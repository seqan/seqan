// ==========================================================================
//                             traverser
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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSER_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSER_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TData, typename TSpec = void>
class Traverser{};

// We need to get the state.
// We need to set the algorithm to the function.
// How can we do this efficently?
template <typename TReference, typename TConfig, typename TJstSpec, typename TSpec>
class Traverser<JournaledStringTree<TReference, TConfig, TJstSpec>, TSpec>
{
public:
    typedef JournaledStringTree<TReference, TConfig, TJstSpec>  TData;
    typedef TraverserContext<TData>                             TContext;

    // ----------------------------------------------------------------------------
    // Member variables.
    // ----------------------------------------------------------------------------

    TContext    _contextOwner;          // Uses own context.
    TContext*   _contextDependent;      // Uses dependent context.

    itRef        // Iterator of window within ref.
    itBranch     // Iterator of window within branch.

    windowSize   // Size of the underlying window.

    void * externalLink;  // Links to external source.



    // ----------------------------------------------------------------------------
    // Member functions.
    // ----------------------------------------------------------------------------

    Traverser() : _contextOwner(), _contextDependent(nullptr)
    {}

    Traverser(TData const & data) : _contextOwner(data), _contextDependent(&_contextOwner)
    {}

    Traverser(TContext & other) :
        _contextOwner(),
        _contextDependent(&_other)
    {}


};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

}

// ============================================================================
// Public Functions
// ============================================================================

// TODO(rrahn): Impl me: goNext(me, lambda)
// TODO(rrahn): Impl me: goPrevious(me, lambda)
// TODO(rrahn): Impl me: context(me)
// TODO(rrahn): Impl me: setWindowSize(me, size)
// TODO(rrahn): Impl me: windowSize(me)
// TODO(rrahn): Impl me: windowRange(me) -> Range of iterator, refIt or not?
// TODO(rrahn): Impl me: position(me)
// TODO(rrahn): Impl me: attach(me, src) -> add pointer to external algo, but how to get the state?
// TODO(rrahn): Impl me: detach(me, src)

}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_TRAVERSER_H_
