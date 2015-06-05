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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSOR_OBSERVABLE_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSOR_OBSERVABLE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TObserver>
struct ImplementObservable;

template <typename TJst, typename TObserver>
class TraversorImpl<TJst, JstTraversalSpec<ImplementObservable<TObserver> > > : public Observable<TObserver>
{
public:

    typedef TraversorImpl<TJst, JstTraversalSpec<void> >    THost;
    typedef typename THost::TTraversalNode                  TTraversalNode;
    typedef typename THost::TTraversalStack                 TTraversalStack;
    typedef typename Container<TraversorImpl>::Type         TContainer;
    typedef typename Size<TraversorImpl>::Type              TSize;

    TContainer*     _contPtr;
    TSize           _contextSize;
    TTraversalNode  _tmp;
    TTraversalStack _stack;

    // __ Constructors ________________________________________________________

    TraversorImpl() : _contPtr(nullptr), _contextSize(1)
    {}

    TraversorImpl(TJst & jst) : _contPtr(nullptr), _contextSize(1), _stack()
    {
        _init(*this, jst);
    }

    // Copy Constructor!
    template <typename TOtherJst, typename TSpec>
    TraversorImpl(TraversorImpl<TOtherJst, JstTraversalSpec<TSpec> > const & other,
                  SEQAN_CTOR_ENABLE_IF(IsConstructible<TJst, TOtherJst>)) :
    _contPtr(other._contPtr),
    _contextSize(other._contextSize),
    _tmp(other._tmp),
    _stack(other._stack)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // __ Member Functions ____________________________________________________

    template <typename TOtherJst, typename TSpec>
    inline SEQAN_FUNC_ENABLE_IF(IsConstructible<TJst, TOtherJst>, TraversorImpl &)
    operator=(TraversorImpl<TOtherJst, JstTraversalSpec<TSpec> > const & other)
    {
        _contPtr = other._contPtr;
        _contextSize = other._contextSize;
        _tmp = other._tmp;
        _stack = other._stack;
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TJst, typename TSpec, typename TObserver>
struct MakeObservable<TraversorImpl<TJst, JstTraversalSpec<TSpec> >, TObserver>
{
    typedef TraversorImpl<TJst, JstTraversalSpec<ImplementObservable<TObserver> > > Type;
};

template <typename TJst, typename TObserver>
struct MakeObservable<TraversorImpl<TJst, JstTraversalSpec<ImplementObservable<TObserver> > >, TObserver>
{
    typedef TraversorImpl<TJst, JstTraversalSpec<ImplementObservable<TObserver> > > Type;
};

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

template <typename TJst, typename TObserver, typename TTraversalNode>
inline void
pushNode(TraversorImpl<TJst, JstTraversalSpec<ImplementObservable<TObserver> > > & traversor,
               TTraversalNode SEQAN_FORWARD_CARG node)
{
    appendValue(traversor._stack, SEQAN_FORWARD(TTraversalNode, node));
    notify(traversor, PushEvent());
}

template <typename TJst, typename TObserver, typename TTraversalNode>
inline void
popNode(TraversorImpl<TJst, JstTraversalSpec<ImplementObservable<TObserver> > > & traversor)
{
    eraseBack(traversor._stack);
    notify(traversor, PopEvent());
}

}

// ============================================================================
// Functions
// ============================================================================

}

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_TRAVERSOR_OBSERVABLE_H_
