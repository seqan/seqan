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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Tests for splitting.
// ==========================================================================

#ifndef TEST_PARALLEL_TEST_PARALLEL_SPLITTING_H_
#define TEST_PARALLEL_TEST_PARALLEL_SPLITTING_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>

SEQAN_DEFINE_TEST(test_parallel_splitter_equidistant)
{
    using namespace seqan;

    // Simple case.
    {
        Splitter<unsigned, Equidistant> splitters(0, 10, 2);
        SEQAN_ASSERT_EQ(length(splitters), 2u);
        SEQAN_ASSERT_EQ(splitters[0], 0u);
        SEQAN_ASSERT_EQ(splitters[1], 5u);
        SEQAN_ASSERT_EQ(splitters[2], 10u);
    }

    // Simple case with offset.
    {
        Splitter<unsigned, Equidistant> splitters(3, 13, 2);
        SEQAN_ASSERT_EQ(length(splitters), 2u);
        SEQAN_ASSERT_EQ(splitters[0], 3u);
        SEQAN_ASSERT_EQ(splitters[1], 8u);
        SEQAN_ASSERT_EQ(splitters[2], 13u);
    }

    // One chunk.
    {
        Splitter<unsigned, Equidistant> splitters(0, 10, 1);
        SEQAN_ASSERT_EQ(length(splitters), 1u);
        SEQAN_ASSERT_EQ(splitters[0], 0u);
        SEQAN_ASSERT_EQ(splitters[1], 10u);
    }

    // Case with empty chunks.
    {
        Splitter<unsigned, Equidistant> splitters(0, 3, 5);
        SEQAN_ASSERT_EQ(length(splitters), 5u);
        SEQAN_ASSERT_EQ(splitters[0], 0u);
        SEQAN_ASSERT_EQ(splitters[1], 1u);
        SEQAN_ASSERT_EQ(splitters[2], 2u);
        SEQAN_ASSERT_EQ(splitters[3], 3u);
        SEQAN_ASSERT_EQ(splitters[4], 3u);
        SEQAN_ASSERT_EQ(splitters[5], 3u);
    }
}

SEQAN_DEFINE_TEST(test_parallel_splitting_compute_splitters)
{
    using namespace seqan;

    // Simple case.
    {
        String<unsigned> splitters;
        computeSplitters(splitters, 10, 2);
        SEQAN_ASSERT_EQ(length(splitters), 3u);
        SEQAN_ASSERT_EQ(splitters[0], 0u);
        SEQAN_ASSERT_EQ(splitters[1], 5u);
        SEQAN_ASSERT_EQ(splitters[2], 10u);
    }

    // One chunk.
    {
        String<unsigned> splitters;
        computeSplitters(splitters, 10, 1);
        SEQAN_ASSERT_EQ(length(splitters), 2u);
        SEQAN_ASSERT_EQ(splitters[0], 0u);
        SEQAN_ASSERT_EQ(splitters[1], 10u);
    }

    // Case with empty chunks.
    {
        String<unsigned> splitters;
        computeSplitters(splitters, 3, 5);
        SEQAN_ASSERT_EQ(length(splitters), 6u);
        SEQAN_ASSERT_EQ(splitters[0], 0u);
        SEQAN_ASSERT_EQ(splitters[1], 1u);
        SEQAN_ASSERT_EQ(splitters[2], 2u);
        SEQAN_ASSERT_EQ(splitters[3], 3u);
        SEQAN_ASSERT_EQ(splitters[4], 3u);
        SEQAN_ASSERT_EQ(splitters[5], 3u);
    }
}

#endif  // TEST_PARALLEL_TEST_PARALLEL_SPLITTING_H_
