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
// Author: David Weese <david.weese@fu-berlin.de>
//         Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Tags to select serial/parallel algorithms.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_TAGS_H_
#define SEQAN_PARALLEL_PARALLEL_TAGS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Parallel
// ----------------------------------------------------------------------------

/*!
 * @defgroup ParallelismTags Parallelism Tags
 * @brief Tags for enabling/disabling parallelism.
 *
 * @tag ParallelismTags#Parallel
 * @headerfile <seqan/parallel.h>
 * @brief Tag to select the parallel implementation of an algorithm.
 *
 * @tag ParallelismTags#Serial
 * @headerfile <seqan/parallel.h>
 * @brief Tag to select the serial implementation of an algorithm.
 */

//struct Parallel_;
//typedef Tag<Parallel_> Parallel;

// ----------------------------------------------------------------------------
// Tag ExecutionPolicy
// ----------------------------------------------------------------------------

// Dynamic execution policy.
template <typename TThreadModel = Serial, typename TVectorSpec = Serial>
struct ExecutionPolicy
{
    // Member variables.
    size_t numThreads = 1;
};

// ----------------------------------------------------------------------------
// Tag ExecutionPolicy  [Sequential]
// ----------------------------------------------------------------------------

constexpr ExecutionPolicy<> seq{};

// ----------------------------------------------------------------------------
// Tag Vectorial
// ----------------------------------------------------------------------------

struct Vectorial_;
using Vectorial = Tag<Vectorial_>;
ExecutionPolicy<Serial, Vectorial> vec{};

// ----------------------------------------------------------------------------
// Tag Parallel Flags depending on parallel option.
// ----------------------------------------------------------------------------

#if defined(SEQAN_TBB)
struct ParallelTbb_;
using ParallelTbb = Tag<ParallelTbb_>;
using Parallel = ParallelTbb;
ExecutionPolicy<ParallelTbb> par{std::thread::hardware_concurrency()};
ExecutionPolicy<ParallelTbb, Vectorial> parVec{std::thread::hardware_concurrency()};
#elsif defined(_OPENMP)
struct ParallelOmp_;
using ParallelOmp = Tag<ParallelOmp_>;
using Parallel = ParallelOmp;
ExecutionPolicy<ParallelOmp> par{std::thread::hardware_concurrency()};
ExecutionPolicy<ParallelOmp, Vectorial> parVec{std::thread::hardware_concurrency()};
#else
struct ParallelStd_;
using ParallelStd = Tag<ParallelStd_>;
using Parallel = ParallelStd;
ExecutionPolicy<ParallelStd> par{std::thread::hardware_concurrency()};
ExecutionPolicy<ParallelStd, Vectorial> parVec{std::thread::hardware_concurrency()};
#endif

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction IsVectorized
// ----------------------------------------------------------------------------

template <typename T>
struct IsVectorized : False
{};

template <>
struct IsVectorized<Vectorial> : True
{};

template <typename TPar, typename TVec>
struct IsVectorized<ExecutionPolicy<TPar, TVec> > :
    public IsVectorized<TVec>
{};

// ----------------------------------------------------------------------------
// Metafunction IsParallel
// ----------------------------------------------------------------------------

template <typename T>
struct IsParallel : False
{};

template <>
struct IsParallel<Parallel> : True
{};

template <typename TPar, typename TVec>
struct IsParallel<ExecutionPolicy<TPar, TVec> > :
    public IsParallel<TPar>
{};

// ----------------------------------------------------------------------------
// Metafunction IsExecutionPolicy
// ----------------------------------------------------------------------------

template <typename T>
struct IsExecutionPolicy : False
{};

template <typename TPar, typename TVec>
struct IsExecutionPolicy<ExecutionPolicy<TPar, TVec> > : True
{};

// ----------------------------------------------------------------------------
// Metafunction DefaultParallelSpec
// ----------------------------------------------------------------------------

template <typename TObject>
struct DefaultParallelSpec
{
    typedef Parallel Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TParallelSpec, typename TVectorizationSpec>
inline auto
numThreads(ExecutionPolicy<TParallelSpec, TVectorizationSpec> const & p)
{
    return p.numThreads;
}

template <typename TParallelSpec, typename TVectorizationSpec>
inline void
setNumThreads(ExecutionPolicy<TParallelSpec, TVectorizationSpec> & p,
              size_t const nt)
{
    p.numThreads = nt;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_TAGS_H_
