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
// Tests for the shuffle() function.
// ==========================================================================

#ifndef TEST_RANDOM_TEST_RANDOM_UTIL_H_
#define TEST_RANDOM_TEST_RANDOM_UTIL_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/random.h>

SEQAN_DEFINE_TEST(test_random_shuffle)
{
    using namespace seqan;

    std::mt19937 mt(0);
    CharString container = "Hello!";
    CharString const before = container;

    shuffle(container, mt);
    SEQAN_ASSERT_NEQ(before, container);
}

SEQAN_DEFINE_TEST(test_random_cvt_beta_param)
{
    using namespace seqan;

    using TParam = BetaDistribution<double>::param_type;
    TParam p = cvtBetaDistParam(0.3, 0.2);
    SEQAN_ASSERT_IN_DELTA(p.alpha(), 1.275, 0.01);
    SEQAN_ASSERT_IN_DELTA(p.beta(), 2.975, 0.01);
}

SEQAN_DEFINE_TEST(test_random_cvt_lognormal_param)
{
    using namespace seqan;

    using TParam = std::lognormal_distribution<double>::param_type;

    TParam p = cvtLogNormalDistParam(1.0, 1.0);
    SEQAN_ASSERT_IN_DELTA(p.m(), -0.346574, 0.01);
    SEQAN_ASSERT_IN_DELTA(p.s(), 0.832555, 0.01);
}

#endif  // TEST_RANDOM_TEST_RANDOM_UTIL_H_
