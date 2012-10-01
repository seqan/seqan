// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Test the specialization Unordered SeedSet.
// ==========================================================================

#ifndef TEST_SEEDS_TEST_SEEDS_SEED_SET_UNORDERED_H_
#define TEST_SEEDS_TEST_SEEDS_SEED_SET_UNORDERED_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

template <typename TSeedSpec>
void testSeedsSeedSetConstructors(TSeedSpec const &)

// Test container functions for specialization unordered.
SEQAN_DEFINE_TEST(test_seeds_seed_set_container_functions_unordered)
{
    using namespace seqan;

    { // Construct with begin/end in both dimensions.
        // Define Seed type and declare a variable.
        typedef Seed<Simple> TSeed;
        TSeed s(1, 2, 3, 5);

        // Check values from construction.
        SEQAN_ASSERT_EQ(1u, getBeginDim0(s));
        SEQAN_ASSERT_EQ(2u, getBeginDim1(s));
        SEQAN_ASSERT_EQ(3u, getEndDim0(s));
        SEQAN_ASSERT_EQ(5u, getEndDim1(s));
        SEQAN_ASSERT_EQ(1, getLowerDiagonal(s));
        SEQAN_ASSERT_EQ(2, getUpperDiagonal(s));
        SEQAN_ASSERT_EQ(1, getStartDiagonal(s));
        SEQAN_ASSERT_EQ(2, getEndDiagonal(s));
    }
    { // Construct from ChainedSeed object.
        typedef Seed<ChainedSeed> TSeed2;
        TSeed2 s2(1, 2, 3);
        typedef Seed<Simple> TSeed;
        TSeed s(s2);

        // Check values from construction.
        SEQAN_ASSERT_EQ(1u, getBeginDim0(s));
        SEQAN_ASSERT_EQ(4u, getEndDim0(s));
        SEQAN_ASSERT_EQ(2u, getBeginDim1(s));
        SEQAN_ASSERT_EQ(5u, getEndDim1(s));
        SEQAN_ASSERT_EQ(1, getLowerDiagonal(s));
        SEQAN_ASSERT_EQ(1, getUpperDiagonal(s));
        SEQAN_ASSERT_EQ(1, getStartDiagonal(s));
        SEQAN_ASSERT_EQ(1, getEndDiagonal(s));
    }
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_SET_UNORDERED_H_
