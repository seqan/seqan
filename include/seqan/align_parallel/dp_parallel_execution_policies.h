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

/*!
 * @tag BlockOffsetOptimization
 * @brief Optimization for vectorized wave-front execution model.
 * @headerfile <seqan/align_parallel.h>
 * @see WavefrontExecutionPolicy
 */
struct BlockOffsetOptimization_;
using BlockOffsetOptimization = Tag<BlockOffsetOptimization_>;

/*!
 * @class WavefrontExecutionPolicy
 * @headerfile <seqan/align_parallel.h>
 * @extends ExecutionPolicy
 * @brief Policy to select runtime execution mode for algorithms.
 * @signature template<typename TWaveSpec, typename TVectorizationMode>
 *            struct ExecutionPolicy<WavefrontAlignment<TWaveSpec>, TVectorizationMode>;
 * @tparam TWaveSpec Type specializing the wave-front threading model.
 *         Can be <tt>void</tt> (default) or @link BlockOffsetOptimization @endlink.
 * @tparam TVectorizationMode Type specifying the vectorization model.
 *         Can be @link ParallelismTags#Vectorial @endlink or @link ParallelismTags#Serial @endlink (default).
 *
 * Special execution policy for computing sequence alignments with wave-front parallelization strategy.
 * In the wave-front execution the DP matrix is partitioned into blocks which can be executed
 * in parallel along the minor diagonal of the DP matrix.
 * The execution policy can be further specialized if used in combination with the @link ParallelismTags#Vectorial @endlink
 * execution mode (see @link WavefrontExecutionPolicy @endlink).
 *
 * @section Vectorization
 *
 * In the vectorization mode, the blocks are gathered into SIMD registers.
 * The @link BlockOffsetOptimization @endlink can be used to always ensure that <tt>sizeof(SIMD) / 2</tt> many blocks
 * can be packed into one SIMD register.
 * This requires, that the available instruction set supports 16 bit packed SIMD operations (e.g. SSE4, AVX2)
 * and the score value type (@link Score @endlink) is bigger then 16 bit.
 * In the default mode, the optimization is disabled and the number of packed alignment blocks is solely determined by
 * the score value type passed to the algorithm as a parameter (e.g. see @link globalAlignmentScore @endlink).
 */
 template <typename TSpec = void>
 struct WavefrontAlignment;

template <typename TSpec, typename TVectorizationSpec>
struct ExecutionPolicy<WavefrontAlignment<TSpec>, TVectorizationSpec> :
    public ExecutionPolicy<Parallel, TVectorizationSpec>
{
    /*!
     *@var size_t WavefrontExecutionPolicy::blockSize
     * @brief The size of the blocks to use. Defaults to 100.
     */
    size_t blockSize{100};
    /*!
     * @var size_t WavefrontExecutionPolicy::parallelAlignments
     * @brief Number of alignments scheduled concurrently. Defaults to <tt>std::thread::hardware_concurrency()</tt>.
     */
    size_t parallelAlignments{std::thread::hardware_concurrency()};
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/*!
 * @fn WavefrontExecutionPolicy#blockSize
 * @brief Getter for the current block size.
 * @signature size_t blockSize(exec);
 * @param[in] exec The wave-front execution policy to query.
 */
template <typename TSpec, typename TVectorizationSpec>
inline auto
blockSize(ExecutionPolicy<WavefrontAlignment<TSpec>, TVectorizationSpec> const & p)
{
    return p.blockSize;
}

/*!
 * @fn WavefrontExecutionPolicy#setBlockSize
 * @brief Setter for the current block size.
 * @signature void setBlockSize(exec, bs);
 * @param[in,out] exec The wave-front execution policy to update.
 * @param[in] bs The new block size to set. Must be a positive integral number greater or equal than 5.
 */
template <typename TSpec, typename TVectorizationSpec>
inline void
setBlockSize(ExecutionPolicy<WavefrontAlignment<TSpec>, TVectorizationSpec> & p,
             size_t const bs)
{
    SEQAN_ASSERT_GEQ(bs, static_cast<size_t>(5));
    p.blockSize = bs;
}

/*!
 * @fn WavefrontExecutionPolicy#parallelAlignments
 * @brief Getter for the current number of alignments executed in parallel.
 * @signature void parallelAlignments(exec);
 * @param[in] exec The wave-front execution policy to update.
 */
template <typename TSpec, typename TVectorizationSpec>
inline auto
parallelAlignments(ExecutionPolicy<WavefrontAlignment<TSpec>, TVectorizationSpec> const & p)
{
    return p.parallelAlignments;
}

/*!
 * @fn WavefrontExecutionPolicy#setParallelAlignments
 * @brief Setter for the current number of alignments executed in parallel.
 * @signature void setParallelAlignments(exec, pa);
 * @param[in,out] exec The wave-front execution policy to update.
 * @param[in] pa The number of alignments to execute in parallel. Must be a positive integral number greater than 0.
 */
template <typename TSpec, typename TVectorizationSpec>
inline void
setParallelAlignments(ExecutionPolicy<WavefrontAlignment<TSpec>, TVectorizationSpec> & p,
                      size_t const pi)
{
    SEQAN_ASSERT_GT(pi, static_cast<size_t>(0));
    p.parallelAlignments = pi;
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_PARALLEL_DP_PARALLEL_EXECUTION_PLOCIES_H_
