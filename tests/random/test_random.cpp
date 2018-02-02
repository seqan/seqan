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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for SeqAn's module random.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/random.h>  // The module under test.

#include "test_random_basic.h"
#include "test_random_beta.h"
#include "test_random_util.h"


SEQAN_BEGIN_TESTSUITE(test_random) {
    // Call Tests.
    SEQAN_CALL_TEST(test_default_rng);

    SEQAN_CALL_TEST(test_random_beta_constructors);
    SEQAN_CALL_TEST(test_random_beta_alpha);
    SEQAN_CALL_TEST(test_random_beta_beta);
    SEQAN_CALL_TEST(test_random_beta_min);
    SEQAN_CALL_TEST(test_random_beta_max);
    SEQAN_CALL_TEST(test_random_beta_param);
    SEQAN_CALL_TEST(test_random_beta_set_param);
    SEQAN_CALL_TEST(test_random_beta_write);
    SEQAN_CALL_TEST(test_random_beta_read);
    SEQAN_CALL_TEST(test_random_beta_functor);

    SEQAN_CALL_TEST(test_random_shuffle);
    SEQAN_CALL_TEST(test_random_cvt_beta_param);
    SEQAN_CALL_TEST(test_random_cvt_lognormal_param);
}
SEQAN_END_TESTSUITE
