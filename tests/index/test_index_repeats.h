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

// TODO(holtgrew): More tests.

#ifndef SEQAN_TESTS_INDEX_TEST_INDEX_REPEATS_H_
#define SEQAN_TESTS_INDEX_TEST_INDEX_REPEATS_H_

#include <seqan/index.h>

// Repeats with period 1.  Test that all locations of Ns are reported
// although the length might be smaller.

SEQAN_DEFINE_TEST(test_index_repeats_period_1_ignore_minlen_for_ns)
{
    using namespace seqan;

    typedef Repeat<unsigned, unsigned> TRepeat;

    String<TRepeat> repeats;
    Dna5String text = "CCCAGATNCGTGATNNNCCCACA";

    findRepeats(repeats, text, 3);  // TODO(holtgrew): repeat len thresh is 1-off

    SEQAN_ASSERT_EQ(length(repeats), 2u);

    SEQAN_ASSERT_EQ(repeats[0].beginPosition, 7u);
    SEQAN_ASSERT_EQ(repeats[0].endPosition, 8u);
    SEQAN_ASSERT_EQ(repeats[0].period, 1u);

    SEQAN_ASSERT_EQ(repeats[1].beginPosition, 14u);
    SEQAN_ASSERT_EQ(repeats[1].endPosition, 17u);
    SEQAN_ASSERT_EQ(repeats[1].period, 1u);
}

// Repeats with period 1.  Test with many ns such that the field of ns spans over threads.

SEQAN_DEFINE_TEST(test_index_repeats_period_1_many_ns)
{
    using namespace seqan;

    typedef Repeat<unsigned, unsigned> TRepeat;

    String<TRepeat> repeats;
    // The following text has a length of 320 characters such that we can also search for repeats in parallel (condition
    // is length > 2 * minLength * num-threads).  It is written such that one repeat spans the border between two threads
    // for four threads and one thread gets no repeats.
    Dna5String text =
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" \
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCGCGCGCGCGAAGCGCGCGCGCGCGCGCGCGC" \
            "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCAAGCGCGCGCGCGCGCGCGCAAGCGCGGCGCGCGCGCGCGCGCGCGC" \
            "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCCGGCGCGCGCGCGCAAA";

    findRepeats(repeats, text, 2);  // TODO(holtgrew): repeat len thresh is 1-off

    SEQAN_ASSERT_EQ(length(repeats), 2u);

    SEQAN_ASSERT_EQ(repeats[0].beginPosition, 0u);
    SEQAN_ASSERT_EQ(repeats[0].endPosition, 127u);
    SEQAN_ASSERT_EQ(repeats[0].period, 1u);

    SEQAN_ASSERT_EQ(repeats[1].beginPosition, 317u);
    SEQAN_ASSERT_EQ(repeats[1].endPosition, 320u);
    SEQAN_ASSERT_EQ(repeats[1].period, 1u);
}

// Repeats with period 1.  Test case with no Ns, here the repeat
// length is always heeded.

SEQAN_DEFINE_TEST(test_index_repeats_period_1_no_ns)
{
    using namespace seqan;

    typedef Repeat<unsigned, unsigned> TRepeat;

    String<TRepeat> repeats;
    // The following text has a length of 320 characters such that we can also search for repeats in parallel (condition
    // is length > 2 * minLength * num-threads).  It is written such that one repeat spans the border between two threads
    // for four threads and one thread gets no repeats.
    Dna5String text =
            "AAAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCA" \
            "AACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGAAGCGCGCGCGCGCGCGCGCGC" \
            "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCAAGCGCGCGCGCGCGCGCGCAAGCGCGGCGCGCGCGCGCGCGCGCGC" \
            "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCCGGCGCGCGCGCGCAAA";

    findRepeats(repeats, text, 2);  // TODO(holtgrew): repeat len thresh is 1-off

    SEQAN_ASSERT_EQ(length(repeats), 3u);

    SEQAN_ASSERT_EQ(repeats[0].beginPosition, 0u);
    SEQAN_ASSERT_EQ(repeats[0].endPosition, 3u);
    SEQAN_ASSERT_EQ(repeats[0].period, 1u);

    SEQAN_ASSERT_EQ(repeats[1].beginPosition, 79u);
    SEQAN_ASSERT_EQ(repeats[1].endPosition, 82u);
    SEQAN_ASSERT_EQ(repeats[1].period, 1u);

    SEQAN_ASSERT_EQ(repeats[2].beginPosition, 317u);
    SEQAN_ASSERT_EQ(repeats[2].endPosition, 320u);
    SEQAN_ASSERT_EQ(repeats[2].period, 1u);
}

#endif  // #ifndef SEQAN_TESTS_INDEX_TEST_INDEX_REPEATS_H_
