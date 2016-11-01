// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_BASE_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_BASE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetPatternState
// ----------------------------------------------------------------------------

template <typename TExtension>
struct GetPatternState
{
    using Type = Nothing;
};

    // ----------------------------------------------------------------------------
    // Metafunction ProxySelectionMethod
    // ----------------------------------------------------------------------------

template <typename TAlgorithm>
struct ProxySelectionMethod
{
    using Type = SelectValidProxy;
};

template <typename TAlgorithm>
struct ProxySelectionMethod<TAlgorithm const> :
    ProxySelectionMethod<TAlgorithm>{};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TPattern>
class JstExtension;

// ----------------------------------------------------------------------------
// Group ContextPostion
// ----------------------------------------------------------------------------

struct ContextBegin_;
using ContextBegin = Tag<ContextBegin_>;

struct ContextEnd_;
using ContextEnd = Tag<ContextEnd_>;

struct ContextRange_;
using ContextRange = Tag<ContextRange_>;

// ----------------------------------------------------------------------------
// Class JstExtensionBase
// ----------------------------------------------------------------------------

template <typename TExtension, typename TCxtPosition = ContextRange>
class JstExtensionBase
{
public:
    using TState = typename GetPatternState<TExtension>::Type;

    TExtension&         _derived;
    mutable TState      _state;

    JstExtensionBase(TExtension & ext) : _derived(ext)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction HasState
// ----------------------------------------------------------------------------

template <typename TObject>
struct HasState;

template <typename TPattern>
struct HasState<JstExtension<TPattern> > :
    If<IsSameType<typename GetPatternState<JstExtension<TPattern> >::Type, Nothing>, False, True>::Type{};

template <typename TExtension, typename TCxtPosition>
struct HasState<JstExtensionBase<TExtension, TCxtPosition> > :
    public HasState<TExtension>{};

// ----------------------------------------------------------------------------
// Metafunction ObservedValue
// ----------------------------------------------------------------------------

template <typename TPattern>
struct ObservedValue<JstExtension<TPattern> > :
    public GetPatternState<JstExtension<TPattern> >{};

template <typename TPattern>
struct ObservedValue<JstExtension<TPattern> const> :
    public GetPatternState<JstExtension<TPattern> const>{};

template <typename TExtension, typename TCxtPosition>
struct ObservedValue<JstExtensionBase<TExtension, TCxtPosition> > :
    public ObservedValue<TExtension>{};

template <typename TExtension, typename TCxtPosition>
struct ObservedValue<JstExtensionBase<TExtension, TCxtPosition> const> :
    public ObservedValue<TExtension const>{};

// ============================================================================
// Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function impl::run()
// ----------------------------------------------------------------------------

template <typename TExtension, typename TTraverser>
inline auto
run(TExtension & extension, TTraverser const & traverser, ContextBegin /*tag*/) ->
    decltype(run(extension, contextBegin(traverser)))
{
    return run(extension, contextBegin(traverser));
}

template <typename TExtension, typename TTraverser>
inline auto
run(TExtension & extension, TTraverser const & traverser, ContextEnd /*tag*/) ->
    decltype(run(extension, contextEnd(traverser)))
{
    return run(extension, contextEnd(traverser));
}

template <typename TExtension, typename TTraverser>
inline auto
run(TExtension & extension, TTraverser const & traverser, ContextRange /*tag*/) ->
    decltype(run(extension, contextBegin(traverser), contextEnd(traverser)))
{
    return run(extension, contextBegin(traverser), contextEnd(traverser));
}

}  // namespace impl

// Returns the state.

// ----------------------------------------------------------------------------
// Function state();
// ----------------------------------------------------------------------------

template <typename TExtension, typename TCxtPosition>
inline typename GetPatternState<TExtension>::Type &
state(JstExtensionBase<TExtension, TCxtPosition> & extBase)
{
    return extBase._state;
}

template <typename TExtension, typename TCxtPosition>
inline typename GetPatternState<TExtension const>::Type &
state(JstExtensionBase<TExtension, TCxtPosition> const & extBase)
{
    return extBase._state;
}

// ----------------------------------------------------------------------------
// Function setState();
// ----------------------------------------------------------------------------

template <typename TExtension, typename TCxtPosition>
inline void
setState(JstExtensionBase<TExtension, TCxtPosition> & extBase,
         typename GetPatternState<TExtension>::Type && state)
{
    extBase._state = state;
}

// ----------------------------------------------------------------------------
// Function getObservedValue();
// ----------------------------------------------------------------------------

template <typename TExtension, typename TCxtPosition>
inline auto
getObservedValue(JstExtensionBase<TExtension, TCxtPosition> & extBase) -> decltype(state(extBase))
{
    return state(extBase);
}

template <typename TExtension, typename TCxtPosition>
inline auto
getObservedValue(JstExtensionBase<TExtension, TCxtPosition> const & extBase) -> decltype(state(extBase))
{
    return state(extBase);
}

// ----------------------------------------------------------------------------
// Function setObservedValue();
// ----------------------------------------------------------------------------

template <typename TExtension, typename TCxtPosition, typename TValue>
inline void
setObservedValue(JstExtensionBase<TExtension, TCxtPosition> & extBase,
                 TValue && val)
{
    setState(extBase, std::forward<TValue>(val));
}

// ----------------------------------------------------------------------------
// Function run()
// ----------------------------------------------------------------------------

template <typename TExtension, typename TCxtPosition,
          typename TTraverser,
          typename TDelegate>
inline auto
run(JstExtensionBase<TExtension, TCxtPosition> & extension,
    TTraverser const & traverser,
    TDelegate && delegate)
#if !defined(COMPILER_WINTEL)
// the intel compiler on windows fails with this decltype, but can auto infer
// the return type itself (possible as of c++14).
    -> decltype(impl::run(extension._derived, traverser, TCxtPosition()).first)
#endif
{
    auto res = impl::run(extension._derived, traverser, TCxtPosition());
    if (res.second)
        delegate();
    return res.first;
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec,
          typename TAlgorithm,
          typename TDelegate,
          typename TObserver>
inline void
find(TraverserImpl<TContainer, JstTraversalSpec<TSpec> > & traverser,
     TAlgorithm & algorithm,
     TDelegate && delegate,
     TObserver & observer)
{
    using TProxySelector = typename ProxySelectionMethod<TAlgorithm>::Type;
    init(traverser, observer, TProxySelector());

    while (!atEnd(traverser))
    {
#if defined(JST_FIND_DEBUG)
        _fillTestSet(traverser);
#endif
        auto steps = run(algorithm, traverser, delegate);
        advance(traverser, steps, observer, TProxySelector());
    }
}

template <typename TContainer, typename TSpec,
          typename TAlgorithm,
          typename TDelegate>
inline void
find(TraverserImpl<TContainer, JstTraversalSpec<TSpec> > & traverser,
     TAlgorithm & algorithm,
     TDelegate && delegate)
{
    if (HasState<TAlgorithm>::VALUE)
    {
        StackObserver<TAlgorithm> algoObserver(algorithm);
        auto observer = makeObserverList(algoObserver);
        find(traverser, algorithm, delegate, observer);
    }
    else
    {  // Disabled observer.
        auto observer = makeObserverList();
        find(traverser, algorithm, delegate, observer);
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_BASE_H_
