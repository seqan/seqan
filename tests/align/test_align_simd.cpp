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
// Author: René Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/stream.h>

#include "test_align_simd.h"

SEQAN_BEGIN_TESTSUITE(test_align_simd)
{

#if defined(__SSE3__) || defined(__AVX2__)
    // Global alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_global_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_global_affine);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_global_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_align_global_affine_banded);

    // Overlap alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_overlap_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_overlap_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_overlap_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_overlap_affine);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_overlap_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_align_overlap_affine_banded);

    // Semi-global alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_semi_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_semi_global_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_align_gaps_semi_global_affine);
    SEQAN_CALL_TEST(test_alignment_algorithms_score_semi_global_affine);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_semi_global_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_align_semi_global_affine_banded);

    // Local alignment.
    SEQAN_CALL_TEST(test_alignment_algorithms_align_local_linear);
    SEQAN_CALL_TEST(test_alignment_algorithms_align_local_affine);

    SEQAN_CALL_TEST(test_alignment_algorithms_align_local_linear_banded);
    SEQAN_CALL_TEST(test_alignment_algorithms_align_local_affine_banded);
#endif  // defined(__SSE3__) || defined(__AVX2__)

}
SEQAN_END_TESTSUITE
