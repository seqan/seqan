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
 *
 * @tag ParallelismTags#Vectorial
 * @headerfile <seqan/parallel.h>
 * @brief Tag to select the vectorized implementation of an algorithm.
 */

struct Parallel_;
typedef Tag<Parallel_> Parallel;

// ----------------------------------------------------------------------------
// Tag Vectorial
// ----------------------------------------------------------------------------

struct Vectorial_;
using Vectorial = Tag<Vectorial_>;

// ----------------------------------------------------------------------------
// Class ExecutionPolicy
// ----------------------------------------------------------------------------

/*!
 * @class ExecutionPolicy
 * @headerfile <seqan/parallel.h>
 * @brief Policy to select runtime execution mode for algorithms.
 * @signature template<typename TThreadingMode, typename TVectorizationMode>
 *            struct ExecutionPolicy;
 * @tparam TThreadingMode Type specifying the threading model.
 *         Can be @link ParallelismTags#Parallel @endlink or @link ParallelismTags#Serial @endlink (default).
 * @tparam TVectorizationMode Type specifying the vectorization model.
 *         Can be @link ParallelismTags#Vectorial @endlink or @link ParallelismTags#Serial @endlink (default).
 *
 * The <tt>ExecutionPolicy</tt> class is used to select different execution models for certain algorithms.
 * Depending on the specialization of the template parameters 4 different modes are possible: sequential, parallel,
 * vectorized and parallel+vectorized. The number of threads for the parallel execution modes can be configured via a
 * member variable.
 */
// Dynamic execution policy.
template <typename TThreadModel = Serial, typename TVectorSpec = Serial>
struct ExecutionPolicy
{
    // Member variables.
    /*!
     * @var size_t ExecutionPolicy::numThreads
     * @brief The number of threads to use for the parallel execution. Defaults to 1.
     */
    size_t _numThreads = 1;
};

// ----------------------------------------------------------------------------
// Tag ExecutionPolicy  [Sequential]
// ----------------------------------------------------------------------------

/*!
 * @typedef Sequential
 * @brief A typedef for the sequential @link ExecutionPolicy @endlink version.
 * @headerfile <seqan/parallel.h>
 *
 * @signature using Sequential = ExecutionPolicy<Serial, Serial>;
 */
using Sequential = ExecutionPolicy<>;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction IsVectorized
// ----------------------------------------------------------------------------

/*!
 * @mfn ExecutionPolicy#IsVectorized
 * @brief Evaluates to @link LogicalValuesTags#True @endlink if the execution policy is
 *        specialized with @link ParallelismTags#Vectorial @endlink.
 * @headerfile <seqan/parallel.h>
 *
 * @signature IsVectorized<TExecPolicy>::VALUE
 * @tparam TTExecPolicy The @link ExecutionPolicy @endlink to check.
 * @return bool <tt>true</tt> if @link ExecutionPolicy @endlink is @link ParallelismTags#Vectorial @endlink,
 *              otherwise <tt>false</tt>.
 */
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

/*!
 * @mfn ExecutionPolicy#IsParallel
 * @brief Evaluates to @link LogicalValuesTags#True @endlink if the execution policy is
 *        specialized with @link ParallelismTags#Parallel @endlink.
 * @headerfile <seqan/parallel.h>
 *
 * @signature IsParallel<TExecPolicy>::VALUE
 * @tparam TTExecPolicy The @link ExecutionPolicy @endlink to check.
 * @return bool <tt>true</tt> if @link ExecutionPolicy @endlink is @link ParallelismTags#Parallel @endlink,
 *              otherwise <tt>false</tt>.
 */
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

/*!
 * @mfn ExecutionPolicy#IsExecutionPolicy
 * @brief Checks if a given type is an @link ExecutionPolicy @endlink.
 * @headerfile <seqan/parallel.h>
 *
 * @signature IsExecutionPolicy<T>::VALUE
 * @tparam T The type to check.
 * @return bool <tt>true</tt> if the <tt>T</tt> is an @link ExecutionPolicy @endlink type,
 *              otherwise <tt>false</tt>.
 */
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

/*!
 * @fn ExecutionPolicy#numThreads
 * @brief Getter function for the thread number.
 * @headerfile <seqan/parallel.h>
 *
 * @signature auto numThreads(exec);
 * @param[in] exec The @link ExecutionPolicy @endlink to get the number of threads for.
 * @return auto The number of threads.
 */
template <typename TParallelSpec, typename TVectorizationSpec>
inline auto
numThreads(ExecutionPolicy<TParallelSpec, TVectorizationSpec> const & p)
{
    return p._numThreads;
}

/*!
 * @fn ExecutionPolicy#setNumThreads
 * @brief Setter function for the thread number.
 * @headerfile <seqan/parallel.h>
 *
 * @signature void setNumThreads(exec, n);
 * @param[in,out] exec The @link ExecutionPolicy @endlink to set the number of threads for.
 * @param[in] n The number of threads used for the parallel execution.
 */
template <typename TParallelSpec, typename TVectorizationSpec>
inline void
setNumThreads(ExecutionPolicy<TParallelSpec, TVectorizationSpec> & p,
              size_t const nt)
{
    p._numThreads = nt;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_TAGS_H_
