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

template <typename TSpec = Default>
struct DPBasicTraits
{
    // The algorithm to choose.
    using TAlgorithm    = typename DefaultDPAlgorithm<TSpec>::Type;
    // The Gaps to choos
    using TGap          = typename DefaultDPGaps<TSpec>::Type;
    // The Band to choose.
    using TBand         = typename DefaultDPBand<TSpec>::Type;
    // The traceback.
    using TTraceback    = typename DefaultDPTraceback<TSpec>::Type;
    // The output to choose.
    using TFormat       = typename DefaultDPFormat<TSpec>::Type;
};

// We have to fill in the policy type.
template <typename TSpec = Default>
struct DPParallelTraits
{
    using TVectorizationPolicy  = typename DefaultDPVectorization<TSpec>::Type;     // Sequential-default, Vectorial<>
    using TParallelPolicy       = typename DefaultDPParallelPolicy<TSpec>::Type;    // Sequential-default, Parallel<Std>, Parallel<Tbb>, Parallel<Omp>?
    using TSchedulingPolicy     = typename DefaultDPSchedulingPolicy<TSpec>::Type;  // Static<>-default, Dynamic<TilingPolicy>  -> is ignored if ParallelPolicy is Sequential!
};

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

template <typename TScore,
          typename TDPBasicTraits = DPBasicTraits<>,
          typename TDPParallelTraits = DPParallelTraits<>>
struct DPConfig : public DPBandConfig<typename TDPBasicTraits::TBandSpec>,
                  public typename impl::DPParallelConfigSelector<typename TDPParallelTraits::TParallelPolicy, typename TDPParallelTraits::TSchedulingPolicy>::Type
{
    // Typedefs.
    // ----------------------------------------------------------------------------
    
    using TWrappedScore = typename impl::DPScoreSelector<typename TScore, typename TDPTraits::TVectorizationPolicy>::Type;
    using TCache        = typename impl::DPCacheSelector<typename TScore, typename TDPTraits::TGap, typename TDPTraits::TVectorizationPolicy>::Type;
    using TScout        = typename impl::DPScoutSelector<typename TDPBasicTraits::TAlgorithm, typename TDPParallelTraits::TParallelPolicy, typename TDPParallelTraits::TVectorizationPolicy>::Type;
    
    // Member Variables.
    // ----------------------------------------------------------------------------
    
    // Required dp policies for the computation.
    TWrappedScore   mScore;
    
    // Presumably private and not intended for user interaction but needed for injecting special behavior to customize the dp code by library developers.
    TCache          _mCache;
    TScout          _mScout;  // We might as well only want the state of the scout?
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename ...Ts>
inline auto const &
getScore(DPConfig<Ts...> const & config)
{
    return config.mScore;
}

template <typename TScore, typename ...Ts>
inline void
setScore(DPConfig<TScore, Ts...> & config,
         TScore const & score)
{
    config.mScore = score;  // Need to implement operator= in score simd wrapper.
}
    
}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BASE_CONFIG_H_