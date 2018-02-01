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

#include <seqan/basic.h>

#include "test_align_wavefront_task_scheduler.h"
#include "test_align_wavefront_alignment_scheduler.h"
#include "test_align_wavefront_intermediate_dp_result.h"
#include "test_align_wavefront_alignment_thread_local.h"

SEQAN_BEGIN_TESTSUITE(test_align_parallel_data_structures)
{
    // -----------------------------------------------------------------------
    // Test wavefront task scheduler.
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_parallel_wavefront_task_scheduler_construct);
    SEQAN_CALL_TEST(test_align_parallel_wavefront_task_scheduler_async);

    // -----------------------------------------------------------------------
    // Test wavefront alignment scheduler.
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_parallel_wavefront_alignment_scheduler_construct);
    SEQAN_CALL_TEST(test_align_parallel_wavefront_alignment_scheduler_async);
    SEQAN_CALL_TEST(test_align_parallel_wavefront_alignment_scheduler_async_with_exception);

    // -----------------------------------------------------------------------
    // Test WavefrontAlignmentResult
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_parallel_intermediate_dp_result_construct);
    SEQAN_CALL_TEST(test_align_parallel_intermediate_dp_result_update_max);
    SEQAN_CALL_TEST(test_align_parallel_intermediate_dp_result_clear);

    // -----------------------------------------------------------------------
    // Test WavefrontAlignmentThreadLocalStorage
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_parallel_wavefront_alignment_thread_local_storage_construt);
    SEQAN_CALL_TEST(test_align_parallel_wavefront_alignment_thread_local_storage_cache);
    SEQAN_CALL_TEST(test_align_parallel_wavefront_alignment_thread_local_storage_intermediate);
}
SEQAN_END_TESTSUITE
