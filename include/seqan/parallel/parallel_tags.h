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

template <typename T>
struct IsExecutionPolicy;

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

struct Parallel_;
typedef Tag<Parallel_> Parallel;

// ----------------------------------------------------------------------------
// Tag ExecutionPolicy
// ----------------------------------------------------------------------------

// Dynamic execution policy.
template <typename TParallelSpec, typename TVectorSpec = Default>
struct ExecutionPolicy{};

// ----------------------------------------------------------------------------
// Tag SerialExecutionPolicy
// ----------------------------------------------------------------------------

constexpr ExecutionPolicy<Serial> ser{};

// ----------------------------------------------------------------------------
// Tag VectorExecutionPolicy
// ----------------------------------------------------------------------------

struct VectorExecutionPolicy_;
using VectorExecutionPolicy = Tag<VectorExecutionPolicy_>;
constexpr ExecutionPolicy<Serial, VectorExecutionPolicy> vec{};

// ----------------------------------------------------------------------------
// Tag ParallelExecutionPolicyTbb
// ----------------------------------------------------------------------------

#if defined(SEQAN_TBB)
struct ParallelExecutionPolicyTbb_;
using ParallelExecutionPolicyTbb = Tag<ParallelExecutionPolicyTbb_>;
constexpr ExecutionPolicy<ParallelExecutionPolicyTbb> parTbb{};
constexpr ExecutionPolicy<ParallelExecutionPolicyTbb, VectorExecutionPolicy> parTbbVec{};
#endif

// ----------------------------------------------------------------------------
// Tag ParallelExecutionPolicyOmp
// ----------------------------------------------------------------------------

#if defined(_OPENMP)
struct ParallelExecutionPolicyOmp_;
using ParallelExecutionPolicyOmp = Tag<ParallelExecutionPolicyOmp_>;
constexpr ExecutionPolicy<ParallelExecutionPolicyOmp> parOmp{};
constexpr ExecutionPolicy<ParallelExecutionPolicyOmp, VectorExecutionPolicy> parOmpVec{};
#endif

// ----------------------------------------------------------------------------
// Tag ParallelExecutionPolicyNative
// ----------------------------------------------------------------------------

struct ParallelExecutionPolicyNative_;
using ParallelExecutionPolicyNative = Tag<ParallelExecutionPolicyNative_>;
constexpr ExecutionPolicy<ParallelExecutionPolicyNative> parNative{};
constexpr ExecutionPolicy<ParallelExecutionPolicyNative, VectorExecutionPolicy> parNativeVec{};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction IsVectorExecutionPolicy
// ----------------------------------------------------------------------------

template <typename T>
struct IsVectorExecutionPolicy : False
{};

template <>
struct IsVectorExecutionPolicy<VectorExecutionPolicy> : True
{};

template <typename TPar, typename TVec>
struct IsVectorExecutionPolicy<ExecutionPolicy<TPar, TVec> > : IsVectorExecutionPolicy<TVec>
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

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_TAGS_H_
