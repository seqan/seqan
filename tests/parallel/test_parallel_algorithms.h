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
// ==========================================================================
// Tests for parallel algorithms.
// ==========================================================================

#ifndef TEST_PARALLEL_TEST_PARALLEL_ALGORITHMS_H_
#define TEST_PARALLEL_TEST_PARALLEL_ALGORITHMS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>

template <typename T1, typename T2>
void compare(T1 &t1, T2 &t2)
{
    SEQAN_ASSERT_EQ_MSG(length(t1), length(t2), "Sequence lengths differ!");
    bool equal = true;
    for(unsigned i = 0; i < length(t1); ++i)
    {
        if (t1[i] != t2[i])
        {
            std::cout << i << '\t' << t1[i] << " != " << t2[i] << std::endl;
            equal = false;
        }
    }
    SEQAN_ASSERT_MSG(equal, "Sequences differ!");
}

SEQAN_DEFINE_TEST(test_parallel_sum)
{
    seqan::String<int> ints;
    appendValue(ints, 4);
    appendValue(ints, 1);
    appendValue(ints, 5); // 10
    appendValue(ints, 3);
    appendValue(ints, 2);
    appendValue(ints, 7); // 22
    appendValue(ints, 4);
    appendValue(ints, 1);
    appendValue(ints, 4); // 31
    appendValue(ints, 5);
    appendValue(ints, 2);
    appendValue(ints, 9); // 47

    SEQAN_ASSERT_EQ(sum(ints, seqan::Serial()), 47);
    SEQAN_ASSERT_EQ(sum(ints, seqan::Parallel()), 47);
}

SEQAN_DEFINE_TEST(test_parallel_partial_sum)
{
    seqan::String<int> ints;
    appendValue(ints, 4);
    appendValue(ints, 1);
    appendValue(ints, 5); // 10
    appendValue(ints, 3);
    appendValue(ints, 2);
    appendValue(ints, 7); // 22
    appendValue(ints, 4);
    appendValue(ints, 1);
    appendValue(ints, 4); // 31
    appendValue(ints, 5);
    appendValue(ints, 2);
    appendValue(ints, 9); // 47

    seqan::String<int> sum1, sum2;
    SEQAN_ASSERT_EQ(partialSum(sum1, ints, seqan::Serial()), 47);
    SEQAN_ASSERT_EQ(partialSum(sum2, ints, seqan::Parallel()), 47);

    SEQAN_ASSERT_EQ(sum1[0], 4);
    SEQAN_ASSERT_EQ(sum1[1], 5);
    SEQAN_ASSERT_EQ(sum1[2], 10);
    SEQAN_ASSERT_EQ(sum1[3], 13);
    SEQAN_ASSERT_EQ(sum1[4], 15);
    SEQAN_ASSERT_EQ(sum1[5], 22);
    SEQAN_ASSERT_EQ(sum1[6], 26);
    SEQAN_ASSERT_EQ(sum1[7], 27);
    SEQAN_ASSERT_EQ(sum1[8], 31);
    SEQAN_ASSERT_EQ(sum1[9], 36);
    SEQAN_ASSERT_EQ(sum1[10], 38);
    SEQAN_ASSERT_EQ(sum1[11], 47);

    compare(sum1, sum2);
    SEQAN_ASSERT_EQ(sum1, sum2);

    SEQAN_ASSERT_EQ(partialSum(ints, ints, seqan::Parallel()), 47);

    compare(sum1, ints);
    SEQAN_ASSERT_EQ(sum1, ints);
}

#endif  // TEST_PARALLEL_TEST_PARALLEL_ALGORITHMS_H_
