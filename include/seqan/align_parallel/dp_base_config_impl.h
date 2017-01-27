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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BASE_CONFIG_IMPL_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BASE_CONFIG_IMPL_H_

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

// Used to optimize out unused config mode.
struct EmptyClass
{};

struct DPParallelConfigBase {
    
    size_t mNumThreads = std::thread::hardware_concurrency();

    void setThreads(size_t const threads)
    {
        mNumThreads = threads;
    }
    
    auto getThreads() const
    {
        return mNumThreads;
    }
};

template <typename TSpec>
class DPParallelConfig;

template <typename TSpec>
class DPParallelConfig<StaticScheduling<TSpec>> : public DPParallelConfigBase
{
    // Can be further specialized later on.
    
};

template <typename TSpec>
class DPParallelConfig<DynamicScheduling<TSpec>> : public DPParallelConfigBase
{
public:
    size_t mTileSize = 100;
    
    void setTileSize(size_t const tileSize)
    {
        mTileSize = tileSize;
    }
    
    auto getTileSize() const
    {
        return mTileSize;
    }
};

template <>
class DPParallelConfig<DynamicScheduling<Limit>> : public DPParallelConfigBase
{
public:
    size_t mTileSize                   = 100;
    size_t mParallelAlignmentInstances = std::thread::hardware_concurrency() << 1;  // twice as threads available..
    
    // Getter and Setter.
    void setTileSize(size_t const tileSize)
    {
        mTileSize = tileSize;
    }
    
    auto getTileSize() const
    {
        return mTileSize;
    }
    
    void setParallelInstances(size_t const num)
    {
        mParallelAlignmentInstances = num;
    }
    
    auto getParallelInstances() const
    {
        return mParallelAlignmentInstances;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// MF: DPCacheSelector
template <typename TScore, typename TGapSpec, typename TSimdFlag>
struct DPCacheSelector
{
    using TScore_           = typename Value<TScore>::Type;
    using TCellValue_       = typename If<IsSimdMode<TSimdFlag>, SimdVector<TScore_>, TScore_>::Type;
    using TDPCell_          = DPCell_<TCellValue_, TGapSpec>;
    using TTraceValue_      = typename TraceBitMap_<TCellValue_>::Type;
    using TAllocatorSpec_   = typename If<IsSimdMode<TSimdFlag>, OverAligned, void>::Type;
                            // Rename to cache.
    using Type              = DPContext<TDPCell_, TTraceValue_,
                                        String<TDPCell_, Alloc<TAllocatorSpec_>>,
                                        String<TTraceValue_, Alloc<TAllocatorSpec_>>>;
};

//MF: DPParallelConfigSelector
template <typename TParallelPolicy, typename TSchedulingPolicy>
struct DPParallelConfigSelector
{
    using Type = typename If<IsSequential<TParallelPolicy>,                // Cond: is_sequential
                             EmptyClass,                                   // then: no config
                             DPParallelConfig<TSchedulingPolicy> >::Type;  // else: parallel config subclassed by scheduling policy.
};

//MF: DPScoreSelector
template <typename TScore, typename TVectorizationPolicy>
struct DPScoreSelector
{
    using Type = typename If<IsSequential<TVectorizationPolicy>,                                                // cond: is sequential
                             TScore,                                                                            // then: unity
                             Score<SimdVector<typename Value<TScore>::Type>, ScoreSimdWrapper<TScore>>>::Type;  // else: score wrapper.
};

// ============================================================================
// Functions
// ============================================================================
    
}  // namespace impl
}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_DP_BASE_CONFIG_IMPL_H_