// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_SETTINGS_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_SETTINGS_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Translates global function interface into setting struct.
template <typename TScoringScheme_, typename TDPTraits = DPTraits::GlobalLinear>
struct DPSettings
{
    using TTraits        = TDPTraits;
    using TScoringScheme = TScoringScheme_;
    using TBandConfig    = DPBandConfig<typename TDPTraits::TBandType>;

    TScoringScheme  scoringScheme;
    TBandConfig     bandScheme;

    DPSettings() = default;

    explicit DPSettings(TScoringScheme score) : scoringScheme(std::move(score))
    {}
};

#ifdef SEQAN_SIMD_ENABLED
// Simd version of DP settings.
template <typename TDPSettings, typename TOffsetSpec = False>
struct SimdDPSettings : public TDPSettings
{
    //-------------------------------------------------------------------------
    // Member Types.

    using TTraits = typename TDPSettings::TTraits;
    using TScoringScheme = typename TDPSettings::TScoringScheme;
    using TScoreValue = typename Value<TScoringScheme>::Type;
    using TScoreValueSimd = typename SimdVector<
                                        std::conditional_t<std::is_same<TOffsetSpec, BlockOffsetOptimization>::value,
                                                           int16_t,
                                                           TScoreValue>>::Type;
    using TSimdScoringScheme = Score<TScoreValueSimd, ScoreSimdWrapper<TScoringScheme>>;

    //-------------------------------------------------------------------------
    // Members.

    TSimdScoringScheme  simdScoringScheme;

    //-------------------------------------------------------------------------
    // Constructor.

    SimdDPSettings() = default;

    explicit SimdDPSettings(TScoringScheme score) :
        TDPSettings(std::move(score)),
        simdScoringScheme(score)
    {}
};
#endif  // SEQAN_SIMD_ENABLED
// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_DP_SETTINGS_H_
