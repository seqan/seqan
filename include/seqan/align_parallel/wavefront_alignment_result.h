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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_ALIGNMENT_RESULT_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_ALIGNMENT_RESULT_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TTraits>
struct WavefrontAlignmentResult
{
    // ----------------------------------------------------------------------------
    // Member Types.

    using TState = std::pair<typename TTraits::TScoreValue, typename TTraits::THostPosition>;

    // ----------------------------------------------------------------------------
    // Member Variables

    TState  mMaxState{minValue<typename TTraits::TScoreValue>(), typename TTraits::THostPosition{}};
    size_t  mTileCol{0};
    size_t  mTileRow{0};

    // ----------------------------------------------------------------------------
    // Constructors.

    // Note: Although, this could be an aggregate type, the icpc-17 crashes,
    // when compiling without the defaulted constructor.
    WavefrontAlignmentResult() = default;

    // ----------------------------------------------------------------------------
    // Member Functions.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

namespace impl
{

template <typename TIntermediate,
          typename TState>
inline void
updateMax(TIntermediate & me,
          TState const & state,
          size_t const tileCol,
          size_t const tileRow)
{
    if (state.first > me.mMaxState.first)
    {
        me.mMaxState = state;
        me.mTileCol = tileCol;
        me.mTileRow = tileRow;
    }
}
}  // namespace impl

template <typename ...TArgs>
inline void
updateMax(WavefrontAlignmentResult<TArgs...> & me,
          typename WavefrontAlignmentResult<TArgs...>::TState const & state,
          size_t const tileCol,
          size_t const tileRow)
{
    impl::updateMax(me, state, tileCol, tileRow);
}

template <typename ...TArgs>
inline void
updateMax(WavefrontAlignmentResult<TArgs...> & lhs,
          WavefrontAlignmentResult<TArgs...> const & rhs)
{
    impl::updateMax(lhs, rhs.mMaxState, rhs.mTileCol, rhs.mTileRow);
}

template <typename ...TArgs>
inline void
clear(WavefrontAlignmentResult<TArgs...> & me)
{
    WavefrontAlignmentResult<TArgs...> tmp;
    swap(me, tmp);
}

template <typename ...TArgs>
inline void
clear(WavefrontAlignmentResult<TArgs...> && me)
{
    WavefrontAlignmentResult<TArgs...> tmp;
    swap(me, tmp);
}

template <typename ...TArgs>
inline typename WavefrontAlignmentResult<TArgs...>::TState const &
value(WavefrontAlignmentResult<TArgs...> const & me)
{
    return me.mMaxState;
}

template <typename ...TArgs>
inline void
swap(WavefrontAlignmentResult<TArgs...> & lhs,
     WavefrontAlignmentResult<TArgs...> & rhs)
{
    // TODO (rrahn): report issue with Intel
    WavefrontAlignmentResult<TArgs...> tmp = std::move(lhs);
    lhs = std::move(rhs);
    rhs = std::move(tmp);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_WAVEFRONT_ALIGNMENT_RESULT_H_
