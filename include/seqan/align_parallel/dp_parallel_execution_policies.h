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
// Policies used for parallel alignment computation.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_PARALLEL_EXECUTION_PLOCIES_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_PARALLEL_EXECUTION_PLOCIES_H_

namespace seqan
{
// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSpec = void>
struct ParallelAlignWavefront;

template <typename TSpec, typename TVectorizationSpec>
struct ExecutionPolicy<ParallelAlignWavefront<TSpec>, TVectorizationSpec> :
    public ExecutionPolicy<Parallel, TVectorizationSpec>
{
    size_t mBlockSize{100};
    size_t mParallelInstances{numThreads << 1};
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TSpec, typename TVectorizationSpec>
inline auto
blockSize(ExecutionPolicy<ParallelAlignWavefront<TSpec>, TVectorizationSpec> const & p)
{
    return p.mBlockSize;
}

template <typename TSpec, typename TVectorizationSpec>
inline void
setBlockSize(ExecutionPolicy<ParallelAlignWavefront<TSpec>, TVectorizationSpec> & p,
             size_t const bs)
{
    p.mBlockSize = bs;
}

template <typename TSpec, typename TVectorizationSpec>
inline auto
parallelInstances(ExecutionPolicy<ParallelAlignWavefront<TSpec>, TVectorizationSpec> const & p)
{
    return p.mParallelInstances;
}

template <typename TSpec, typename TVectorized>
inline void
setParallelInstances(ExecutionPolicy<ParallelAlignWavefront<TSpec>, TVectorizationSpec> & p,
                     size_t const pi)
{
    p.mParallelInstances = pi;
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_DP_PARALLEL_EXECUTION_PLOCIES_H_
