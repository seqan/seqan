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
// Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_PLATFORM_TEST_PLATFORM_RANGE_BASED_FOR_LOOPS_H_
#define TESTS_PLATFORM_TEST_PLATFORM_RANGE_BASED_FOR_LOOPS_H_

#include <string>
#include <vector>

struct String
{
    std::string me;
};

inline char *
begin(String & str)
{
    return &str.me[0];
}

inline char *
end(String & str)
{
    return &str.me[0] + str.me.length();
}

void foobar1()
{
    String hello{"world"};

    std::vector<char> hits;
    for (auto c: hello)
    {
        hits.emplace_back(c);
    }

    SEQAN_ASSERT_EQ(hits.size(), 5u);
    SEQAN_ASSERT_EQ(hits[0], 'w');
    SEQAN_ASSERT_EQ(hits[1], 'o');
    SEQAN_ASSERT_EQ(hits[2], 'r');
    SEQAN_ASSERT_EQ(hits[3], 'l');
    SEQAN_ASSERT_EQ(hits[4], 'd');
}

template <bool speed_up>
void foobar2()
{
    String hello{"world"};

    std::vector<char> hits;
    for (auto c: hello)
    {
        hits.emplace_back(c);
    }

    SEQAN_ASSERT_EQ(hits.size(), 5u);
    SEQAN_ASSERT_EQ(hits[0], 'w');
    SEQAN_ASSERT_EQ(hits[1], 'o');
    SEQAN_ASSERT_EQ(hits[2], 'r');
    SEQAN_ASSERT_EQ(hits[3], 'l');
    SEQAN_ASSERT_EQ(hits[4], 'd');
}

SEQAN_DEFINE_TEST(test_platform_range_based_for_loops1) {
    foobar1();
}

SEQAN_DEFINE_TEST(test_platform_range_based_for_loops2) {
    #if defined(__INTEL_COMPILER)
        #if __INTEL_COMPILER < 1600 || (__INTEL_COMPILER == 1600 && __INTEL_COMPILER_UPDATE <= 3)
            SEQAN_SKIP_TEST;
        #endif
    #endif

    foobar2<true>();
}

#endif
