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

template <typename TPattern>
struct GetPatternState
{
    using Type = Nothing;
};

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

template <typename TPattern, typename THasState, typename TCxtPosition = ContextRange>
class JstExtensionBase;

template <typename TPattern, typename TCxtPosition>
class JstExtensionBase<TPattern, False, TCxtPosition>
{
public:
    TPattern & _derived;

    JstExtensionBase(TPattern & pattern) : _derived(pattern)
    {}
};

template <typename TExtension, typename TCxtPosition>
class JstExtensionBase<TExtension, True, TCxtPosition>
{
public:
    using TState = typename GetPatternState<TExtension>::Type;

    TExtension& _derived;
    TState      _state;

    JstExtensionBase(TExtension & ext) : _derived(ext)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TObject>
struct HasState : False{};

template <typename TExtension, typename THasState, typename TCxtPosition>
struct HasState<JstExtensionBase<TExtension, THasState, TCxtPosition> > : THasState{};

// ============================================================================
// Functions
// ============================================================================

namespace impl
{
template <typename TExtension, typename TTraverser>
inline auto
run(TExtension & extension, TTraverser & traverser, ContextBegin /*tag*/) ->
    decltype(run(extension, contextBegin(traverser)))
{
    return run(extension, contextBegin(traverser));
}

template <typename TExtension, typename TTraverser>
inline auto
run(TExtension & extension, TTraverser & traverser, ContextEnd /*tag*/) ->
    decltype(run(extension, contextEnd(traverser)))
{
    return run(extension, contextEnd(traverser));
}

template <typename TExtension, typename TTraverser>
inline auto
run(TExtension & extension, TTraverser & traverser, ContextRange /*tag*/) ->
    decltype(run(extension, contextBegin(traverser), contextEnd(traverser)))
{
    return run(extension, contextBegin(traverser), contextEnd(traverser));
}

}  // namespace impl

// Returns the state.
template <typename TExtension, typename THasState, typename TCxtPosition>
inline typename GetPatternState<TExtension>::Type &
state(JstExtensionBase<TExtension, THasState, TCxtPosition> & extBase)
{
    return extBase._state;
}

template <typename TExtension, typename THasState, typename TCxtPosition>
inline typename GetPatternState<TExtension const>::Type &
state(JstExtensionBase<TExtension, THasState, TCxtPosition> const & extBase)
{
    return extBase._state;
}

// ----------------------------------------------------------------------------
// Function find
// ----------------------------------------------------------------------------

template <typename TExtension, typename THasState, typename TCxtPosition,
          typename TTraverser,
          typename TDelegate>
inline auto
run(JstExtensionBase<TExtension, THasState, TCxtPosition> & extension,
    TTraverser & traverser,
    TDelegate && delegate) -> decltype(impl::run(extension._derived, traverser, TCxtPosition()).first)
{
    auto res = impl::run(extension._derived, traverser, TCxtPosition());
    if (res.second)
        delegate();
    return res.first;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_JST_EXTENSION_BASE_H_
