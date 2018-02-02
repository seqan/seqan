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
// Tests for misc/misc_accumulators.h.
// ==========================================================================

#ifndef SEQAN_TESTS_MISC_TEST_MISC_ACCUMULATORS_H_
#define SEQAN_TESTS_MISC_TEST_MISC_ACCUMULATORS_H_

#include <seqan/misc/accumulators.h>

#include <seqan/basic.h>

SEQAN_DEFINE_TEST(test_misc_accumulators_average_accumulator_int_average)
{
    using namespace seqan;

    Accumulator<int, AccuAverage> acc;

    SEQAN_ASSERT_EQ(sum(acc), 0);

    push(acc, 1);
    push(acc, 9);
    push(acc, 5);

    SEQAN_ASSERT_IN_DELTA(average(acc), 5.0, 0.001);
}

SEQAN_DEFINE_TEST(test_misc_accumulators_average_accumulator_int_count)
{
    using namespace seqan;

    Accumulator<int, AccuAverage> acc;

    SEQAN_ASSERT_EQ(sum(acc), 0);

    push(acc, 1);
    push(acc, 9);
    push(acc, 5);

    SEQAN_ASSERT_EQ(count(acc), 3u);
}

SEQAN_DEFINE_TEST(test_misc_accumulators_average_accumulator_int_sum)
{
    using namespace seqan;

    Accumulator<int, AccuAverage> acc;

    SEQAN_ASSERT_EQ(sum(acc), 0);

    push(acc, 1);
    push(acc, 9);
    push(acc, 5);

    SEQAN_ASSERT_EQ(sum(acc), 15);
}

SEQAN_DEFINE_TEST(test_misc_accumulators_average_accumulator_int_clear)
{
    using namespace seqan;

    Accumulator<int, AccuAverage> acc;

    SEQAN_ASSERT_EQ(count(acc), 0u);

    push(acc, 1);

    SEQAN_ASSERT_EQ(count(acc), 1u);

    clear(acc);

    SEQAN_ASSERT_EQ(count(acc), 0u);
}

SEQAN_DEFINE_TEST(test_misc_accumulators_average_accumulator_double_average)
{
    using namespace seqan;

    Accumulator<double, AccuAverage> acc;

    SEQAN_ASSERT_EQ(sum(acc), 0);

    push(acc, 1);
    push(acc, 10);
    push(acc, 5);

    SEQAN_ASSERT_IN_DELTA(average(acc), 5.3333, 0.001);
}

SEQAN_DEFINE_TEST(test_misc_accumulators_average_accumulator_double_count)
{
    using namespace seqan;

    Accumulator<double, AccuAverage> acc;

    SEQAN_ASSERT_EQ(sum(acc), 0);

    push(acc, 1);
    push(acc, 9);
    push(acc, 5);

    SEQAN_ASSERT_EQ(count(acc), 3u);
}

SEQAN_DEFINE_TEST(test_misc_accumulators_average_accumulator_double_sum)
{
    using namespace seqan;

    Accumulator<double, AccuAverage> acc;

    SEQAN_ASSERT_EQ(sum(acc), 0);

    push(acc, 1.1);
    push(acc, 9.2);
    push(acc, 5.3);

    SEQAN_ASSERT_IN_DELTA(sum(acc), 15.6, 0.001);
}

SEQAN_DEFINE_TEST(test_misc_accumulators_average_accumulator_double_clear)
{
    using namespace seqan;

    Accumulator<double, AccuAverage> acc;

    SEQAN_ASSERT_EQ(count(acc), 0u);

    push(acc, 1);

    SEQAN_ASSERT_EQ(count(acc), 1u);

    clear(acc);

    SEQAN_ASSERT_EQ(count(acc), 0u);
}

#endif  // #ifndef SEQAN_TESTS_MISC_TEST_MISC_ACCUMULATORS_H_
