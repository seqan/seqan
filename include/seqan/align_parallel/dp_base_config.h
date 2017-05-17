// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BASE_CONFIG_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BASE_CONFIG_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Parallel Modes:

// SEQ_NOSIMD_NOTILING  -> for loop over all sequences/or sequence and call single alignment.
// SEQ_NOSIMD_TILING    -> NaN!
// SEQ_SIMD_NOTILING    -> for loop over all sequences in batches and call on single thread -> do not spawn new thread.
// SEQ_SIMD_TILING      -> NaN!

// For latter options, we want to set the thread number. => default is std::thread::hardware_concurrency();
// PAR_NOSIMD_NOTILING  -> parallel for loop over all sequences
// PAR_NOSIMD_TILING    -> alignment_scheduler (QueueSpec<unlimited,limited>, TBBQueue?, STDQueue?)+ thread pool.  // Set the limit on alignment scheduler, set the block size
// PAR_SIMD_NOTILING    -> parallel for loop over all sequences plus process in vector batches.
// PAR_SIMD_TILING      -> alignment_scheduler (QueueSpec<unlimited,limited>, TBBQueue?, STDQueue?)+ thread pool + dp task pool.  // Set the limit on alignment scheduler, set the block size

// Use interface as: DPSettings<TScore, TBand, Traits<>>
template <typename TScoringScheme_, typename TDPTraits = DPTraits::GlobalLinear>
struct DPSettings
{
    using TTraits        = TDPTraits;
    using TScoringScheme = TScoringScheme_;
    using TBandConfig    = DPBandConfig<typename TDPTraits::TBandType>;

    TScoringScheme  mScoringScheme;
    TBandConfig     mBandScheme;

    DPSettings() = default;

    explicit DPSettings(TScoringScheme score) : mScoringScheme(std::move(score))
    {}
};

#ifdef SEQAN_SIMD_ENABLED
// Aggregate type.
template <typename TDPSettings, typename TOffsetSpec = False>
struct SimdDPSettings : public TDPSettings
{
    //-------------------------------------------------------------------------
    // Member Types.

    using TTraits = typename TDPSettings::TTraits;
    using TScoringScheme = typename TDPSettings::TScoringScheme;
    using TScoreValue = typename Value<TScoringScheme>::Type;
    using TScoreValueSimd = typename SimdVector<std::conditional_t<TOffsetSpec::VALUE, int16_t, TScoreValue>>::Type;
    using TSimdScoringScheme = Score<TScoreValueSimd, ScoreSimdWrapper<TScoringScheme>>;

    //-------------------------------------------------------------------------
    // Members.

    TSimdScoringScheme  mSimdScoringScheme;

    //-------------------------------------------------------------------------
    // Constructor.

    SimdDPSettings() = default;

    explicit SimdDPSettings(TScoringScheme score) :
        TDPSettings(std::move(score)),
        mSimdScoringScheme(score)
    {}
};
#endif
// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BASE_CONFIG_H_
