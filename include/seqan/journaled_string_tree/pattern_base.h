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

#ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_PATTERN_BASE_H_
#define INCLUDE_SEQAN_JOURNALED_STRING_TREE_PATTERN_BASE_H_

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

template <typename TNeedle, typename TSpec>
class Pattern2;

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
// Class PatternBase
// ----------------------------------------------------------------------------

template <typename TPattern, typename THasState, typename TCxtPosition = ContextRange>
class PatternBase;

template <typename TPattern, typename TCxtPosition>
class PatternBase<TPattern, False, TCxtPosition>
{
public:
    TPattern & _derived;

    PatternBase(TPattern & pattern) : _derived(pattern)
    {}
};


template <typename TPattern, typename TCxtPosition>
class PatternBase<TPattern, True, TCxtPosition>
{
public:
    using TState = typename GetPatternState<TPattern>::Type;

    TPattern & _derived;
    TState _state;

    PatternBase(TPattern & pattern) : _derived(pattern)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TPattern, typename TState, typename TCxtPosition>
struct HasState : public TState{};

/*!
 * @mfn Pattern#Container
 * @brief Returns the needle type of the pattern.
 *
 * @signature Container<TPattern>::Type;
 *
 * @tparam TPattern The pattern to query for its needle type.
 *
 * @return Type The needle type.
 */

template <typename TNeedle, typename TSpec>
struct Container< Pattern2<TNeedle, TSpec> > {
    typedef TNeedle Type;
};

template <typename TNeedle, typename TSpec>
struct Container< Pattern2<TNeedle, TSpec> const > {
    typedef TNeedle const Type;
};

/*!
 * @mfn Pattern#Host
 * @brief Returns the host type of the pattern.
 *
 * @signature Host<TPattern>::Type;
 *
 * @tparam TPattern The pattern to query for its host type.
 *
 * @return Type The host type.
 */

template <typename TNeedle, typename TSpec>
struct Host< Pattern2<TNeedle, TSpec> > {
    typedef TNeedle Type;
};

template <typename TNeedle, typename TSpec>
struct Host< Pattern2<TNeedle, TSpec> const > {
    typedef TNeedle const Type;
};

template <typename TNeedle, typename TSpec>
struct Needle< Pattern2<TNeedle, TSpec> > : Host<Pattern2<TNeedle, TSpec> >
{};

template <typename TNeedle, typename TSpec>
struct Needle< Pattern2<TNeedle, TSpec> const > : Host<Pattern2<TNeedle, TSpec> const>
{};


/*!
 * @mfn Pattern#Value
 * @brief Returns the value type of the underlying pattern.
 *
 * @signature Value<TPattern>::Type;
 *
 * @tparam TPattern The Pattern to query.
 *
 * @return Type The value type.
 */

template <typename TPattern, typename TSpec>
struct Value< Pattern2<TPattern, TSpec> > {
    typedef typename Value<TPattern>::Type Type;
};

/*!
 * @mfn Pattern#Position
 * @brief Returns the position type of the underlying pattern.
 *
 * @signature Position<TPattern>::Type;
 *
 * @tparam TPattern The Pattern to query.
 *
 * @return Type The position type.
 */

template <typename TPattern, typename TSpec>
struct Position< Pattern2<TPattern, TSpec> > {
    typedef typename Position<TPattern>::Type Type;
};

/*!
 * @mfn Pattern#Difference
 * @brief Returns the difference type of the underlying pattern.
 *
 * @signature Difference<TPattern>::Type;
 *
 * @tparam TPattern The Pattern to query.
 *
 * @return Type The difference type.
 */

template <typename TPattern, typename TSpec>
struct Difference< Pattern2<TPattern, TSpec> > {
    typedef typename Difference<TPattern>::Type Type;
};

/*!
 * @mfn Pattern#Size
 * @brief Returns the size type of the underlying pattern.
 *
 * @signature Size<TPattern>::Type;
 *
 * @tparam TPattern The Pattern to query.
 *
 * @return Type The size type.
 */

template <typename TPattern, typename TSpec>
struct Size< Pattern2<TPattern, TSpec> > {
    typedef typename Size<TPattern>::Type Type;
};


// ============================================================================
// Functions
// ============================================================================

namespace impl
{
template <typename TPattern, typename TTraverser>
inline auto
run(TPattern & pattern, TTraverser & traverser, ContextBegin /*tag*/) ->
    decltype(run(pattern, contextBegin(traverser)))
{
    return run(pattern, contextBegin(traverser));
}

template <typename TPattern, typename TTraverser>
inline auto
run(TPattern & pattern, TTraverser & traverser, ContextEnd /*tag*/) ->
    decltype(run(pattern, contextEnd(traverser)))
{
    return run(pattern, contextEnd(traverser));
}

template <typename TPattern, typename TTraverser>
inline auto
run(TPattern & pattern, TTraverser & traverser, ContextRange /*tag*/) ->
    decltype(run(pattern, contextBegin(traverser), contextEnd(traverser)))
{
    return run(pattern, contextBegin(traverser), contextEnd(traverser));
}

}  // namespace impl

template <typename TNeedle, typename TSpec>
inline Holder<TNeedle> &
_dataHost(Pattern2<TNeedle, TSpec> & me)
{
    return me.data_host;
}

template <typename TNeedle, typename TSpec>
inline Holder<TNeedle> const &
_dataHost(Pattern2<TNeedle, TSpec> const & me)
{
    return me.data_host;
}

template <typename TNeedle, typename TSpec, typename TNeedle2>
inline void
setHost(Pattern2<TNeedle, TSpec> & me,
        TNeedle2 && ndl)
{
    SEQAN_ASSERT(!empty(ndl));
    setValue(_dataHost(me), std::forward<TNeedle2>(ndl));
    _reinitPattern(me);
}

template < typename TNeedle, typename TSpec >
inline typename Needle< Pattern2<TNeedle, TSpec> >::Type &
needle(Pattern2<TNeedle, TSpec> & obj)
{
    return host(obj);
}

template < typename TNeedle, typename TSpec >
inline typename Needle< Pattern2<TNeedle, TSpec> const>::Type &
needle(Pattern2<TNeedle, TSpec> const & obj)
{
    return host(obj);
}

// Returns the state.
template <typename TPattern, typename THasState, typename TCxtPosition>
inline typename GetPatternState<TPattern>::Type &
state(PatternBase<TPattern, THasState, TCxtPosition> & pattern)
{
    return pattern._state;
}

template <typename TPattern, typename THasState, typename TCxtPosition>
inline typename GetPatternState<TPattern const>::Type &
state(PatternBase<TPattern, THasState, TCxtPosition> const & pattern)
{
    return pattern._state;
}

template <typename TPattern, typename THasState, typename TCxtPosition,
          typename TTraverser,
          typename TDelegate>
inline auto
find2(PatternBase<TPattern, THasState, TCxtPosition> & pattern,
     TTraverser & traverser,
     TDelegate & delegate) -> decltype(impl::run(pattern._derived, traverser, TCxtPosition()).first)
{
    auto res = impl::run(pattern._derived, traverser, TCxtPosition());
    if (res.second)
        delegate();
    return res.first;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_STRING_TREE_PATTERN_BASE_H_
