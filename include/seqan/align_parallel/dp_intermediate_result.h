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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_INTERMEDIATE_RESULT_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_INTERMEDIATE_RESULT_H_

namespace seqan
{
namespace impl
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TDPScout>
struct IndermediateResultTraits
{
    using TScore = typename std::decay<decltype(maxScore(TDPScout{}))>::type;
    using TPos   = typename std::decay<decltype(maxHostPosition(TDPScout{}))>::type;

    using TState      = std::pair<TScore, TPos>;
    using TComparator = std::less<TState>;
};

template <typename TDPScout,
          typename TTraits = IndermediateResultTraits<TDPScout>>
struct IntermediateResult
{
    // Typedefs
    // ----------------------------------------------------------------------------

    // Member Variables
    // ----------------------------------------------------------------------------

    TState  mState{};       // Requires: Default-Constructible!, Copy-Constructible!, Move-Constructible+.
    size_t  mTileCol{0};
    size_t  mTileRow{0};

    // Constructors
    // ----------------------------------------------------------------------------

    // Member Functions.
    // ----------------------------------------------------------------------------
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TDPScout, typename TTraits>
struct Traits<IntermediateResult<TDPScout, TTraits>>
{
    using Type = TTraits;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TScout, typename TTraits, typename TState>
void update(IntermediateResult<TScout, TTraits> & me,
            TState && state,
            size_t tileCol,
            size_t tileRow)
{
    using TComp = typename TTraits::TComparator;
    if (TComp{}(me.mState, std::forward<TSingleMaxState>(state)))
    {
        std::swap(me.mState, state);
        me.mTileCol = tileCol;
        me.mTileRow = tileRow;
    }
}

template <typename TScout, typename TTraits, typename TSingleMaxState>
void reset(IntermediateResult<TScout, TTraits> & me)
{
    using std::swap;

    using TState = typename IntermediateResult<TScout>::TState;

    swap(me.mState, TState{});
    mTileCol = 0;
    mTileRow = 0;
}

template <typename TScout, typename TTraits, typename TSingleMaxState>
auto const & value(IntermediateResult<TScout, TTraits> const & me)
{
    return me.mState;
}

template <typename TAlgorithm>
bool void isTrackingEnabled(bool const inLastColumn, bool const inLastRow)
{
    if (inLastColumn && inLastRow)
        return true;

}

}  // namespace impl
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_INTERMEDIATE_RESULT_H_
