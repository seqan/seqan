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
// Author: Manuel Holtgrewe <manuel.holtgrew@fu-berlin.de>
// ==========================================================================
// Tests for the module parallel.
// ==========================================================================

#include <seqan/basic.h>     // Testing infrastructure.
#include <seqan/stream.h>      // Required to print strings in tests.
#include <seqan/parallel.h>  // Header under test.

#if defined(_OPENMP)
#include <omp.h>
#endif  // #if defined(_OPENMP)

#include "test_parallel_atomic_primitives.h"
#include "test_parallel_atomic_misc.h"
#include "test_parallel_splitting.h"
#include "test_parallel_algorithms.h"
#include "test_parallel_queue.h"
#include "test_parallel_thread_pool.h"
#include "test_parallel_enumerable_thread_local.h"

SEQAN_BEGIN_TESTSUITE(test_parallel) {
#if defined(_OPENMP)
    // Set number of threads to >=2 so there actually is parallelism.
    if (omp_get_max_threads() < 2)
        omp_set_num_threads(2);
    std::cout << "PARALLELISM ENABLED (CTEST_FULL_OUTPUT)" << std::endl;
#endif  // #if defined(_OPENMP)

    // TODO(holtgrew): Re-enable tests on LLVM when bug 9041 is fixed.
    // LLVM has problems with atomic operation builtins, re-enable when
    // this problem is fixed. See http://llvm.org/bugs/show_bug.cgi?id=9041
    //
    // There is a problem with compare-and-swap on MinGW, too.
#if !defined(__llvm__)
    // Tests for atomic primitives.
    SEQAN_CALL_TEST(test_parallel_atomic_inc);
    SEQAN_CALL_TEST(test_parallel_atomic_dec);
    SEQAN_CALL_TEST(test_parallel_atomic_add);
    SEQAN_CALL_TEST(test_parallel_atomic_or);
    SEQAN_CALL_TEST(test_parallel_atomic_xor);
    SEQAN_CALL_TEST(test_parallel_atomic_cas);

    // Tests for misc simpmle atomic operations.
    SEQAN_CALL_TEST(test_parallel_atomic_min);
    SEQAN_CALL_TEST(test_parallel_atomic_max);
#endif  // #if !defined(__llvm__)

    SEQAN_CALL_TEST(test_parallel_splitter_equidistant);
    SEQAN_CALL_TEST(test_parallel_splitting_compute_splitters);
    SEQAN_CALL_TEST(test_parallel_sum);
    SEQAN_CALL_TEST(test_parallel_partial_sum);

    // Tests for parallel queue.
    SEQAN_CALL_TEST(test_parallel_queue_simple);
    SEQAN_CALL_TEST(test_parallel_queue_resize);
    SEQAN_CALL_TEST(test_parallel_queue_non_pod);

    if (std::thread::hardware_concurrency() >= 2u)
    {
        SEQAN_CALL_TEST(test_parallel_queue_spsc_fixedsize);
        SEQAN_CALL_TEST(test_parallel_queue_spsc_dynamicsize);
        SEQAN_CALL_TEST(test_parallel_queue_spmc_fixedsize);
        SEQAN_CALL_TEST(test_parallel_queue_spmc_dynamicsize);
        SEQAN_CALL_TEST(test_parallel_queue_mpsc_fixedsize);
//        SEQAN_CALL_TEST(test_parallel_queue_mpsc_dynamicsize);
//        SEQAN_CALL_TEST(test_parallel_queue_mpmc_fixedsize);
//        SEQAN_CALL_TEST(test_parallel_queue_mpmc_dynamicsize);
    }

    // -----------------------------------------------------------------------
    // Test thread pool.
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_parallel_thread_pool_construct);
    SEQAN_CALL_TEST(test_parallel_thread_pool_spawn);
    SEQAN_CALL_TEST(test_parallel_thread_pool_join);
    SEQAN_CALL_TEST(test_parallel_thread_pool_destruct);

    // -----------------------------------------------------------------------
    // Test Enumerable Thread Specific.
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_parallel_enumerable_thread_local_construct);
    SEQAN_CALL_TEST(test_parallel_enumerable_thread_local_local);
    SEQAN_CALL_TEST(test_parallel_enumerable_thread_local_enumerate);
    SEQAN_CALL_TEST(test_parallel_enumerable_thread_local_combine_unary);
    SEQAN_CALL_TEST(test_parallel_enumerable_thread_local_combine_binary);
}
SEQAN_END_TESTSUITE
