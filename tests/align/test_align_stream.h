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
// Writing of alignment data structures to streams.
// ==========================================================================

#ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_STREAM_H_
#define SEQAN_TESTS_ALIGN_TEST_ALIGN_STREAM_H_

#include <sstream>

#include <seqan/align.h>

// Test function write() on Array Gaps.

SEQAN_DEFINE_TEST(test_align_stream_gaps_write)
{
    seqan::Gaps<seqan::Dna5String> gaps;
    assignSource(gaps, "CGATTTAT");
    insertGaps(gaps, 8, 2);
    insertGaps(gaps, 5, 1);
    insertGaps(gaps, 0, 2);

    std::stringstream ss;
    ss << gaps;

    std::stringstream expected;
    expected << "--CGATT-TAT--";
    SEQAN_ASSERT_EQ(expected.str(), ss.str());
}

// Test operator<<() on Array Gaps.

SEQAN_DEFINE_TEST(test_align_stream_gaps_stream)
{
    seqan::Gaps<seqan::Dna5String> gaps;
    assignSource(gaps, "CGATTTAT");
    insertGaps(gaps, 8, 2);
    insertGaps(gaps, 5, 1);
    insertGaps(gaps, 0, 2);

    std::stringstream ss;
    ss << gaps;

    std::stringstream expected;
    expected << "--CGATT-TAT--";
    SEQAN_ASSERT_EQ(expected.str(), ss.str());
}

// Test write() on Align.

SEQAN_DEFINE_TEST(test_align_stream_align_write)
{
    seqan::Align<seqan::Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "CGATTAGTG");
    assignSource(row(align, 1), "ATTAGTA");
    globalAlignment(align, seqan::EditDistanceScore());

    std::stringstream ss;
    ss << align;

    std::stringstream expected;
    expected << "      0     .     \n"
             << "        CGATTAGTG\n"
             << "          |||||| \n"
             << "        --ATTAGTA\n\n\n";
    SEQAN_ASSERT_EQ(expected.str(), ss.str());
}

// Test operator<<() on Align.

SEQAN_DEFINE_TEST(test_align_stream_align_stream)
{
    seqan::Align<seqan::Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "CGATTAGTG");
    assignSource(row(align, 1), "ATTAGTA");
    globalAlignment(align, seqan::EditDistanceScore());

    std::stringstream ss;
    ss << align;

    std::stringstream expected;
    expected << "      0     .     \n"
             << "        CGATTAGTG\n"
             << "          |||||| \n"
             << "        --ATTAGTA\n\n\n";
    SEQAN_ASSERT_EQ(expected.str(), ss.str());
}

#endif  // SEQAN_TESTS_ALIGN_TEST_ALIGN_STREAM_H_
