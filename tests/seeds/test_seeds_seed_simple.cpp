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
// Test the specialization Simple Seed.
// ==========================================================================

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/stream.h>   // Required to print strings in tests.

#include <seqan/seeds.h>  // Include module under test.

// Test specific construtors.
SEQAN_DEFINE_TEST(test_seeds_seed_simple_constructors)
{
    using namespace seqan;

    { // Construct with begin/end in both dimensions.
        // Define Seed type and declare a variable.
        typedef Seed<Simple> TSeed;
        TSeed s(1, 2, 3, 5);

        // Check values from construction.
        SEQAN_ASSERT_EQ(1u, beginPositionH(s));
        SEQAN_ASSERT_EQ(2u, beginPositionV(s));
        SEQAN_ASSERT_EQ(3u, endPositionH(s));
        SEQAN_ASSERT_EQ(5u, endPositionV(s));
        SEQAN_ASSERT_EQ(-2, lowerDiagonal(s));
        SEQAN_ASSERT_EQ(-1, upperDiagonal(s));
        SEQAN_ASSERT_EQ(-1, beginDiagonal(s));
        SEQAN_ASSERT_EQ(-2, endDiagonal(s));
    }
    { // Construct from ChainedSeed object.
        typedef Seed<ChainedSeed> TSeed2;
        TSeed2 s2(1, 2, 3);
        typedef Seed<Simple> TSeed;
        TSeed s(s2);

        // Check values from construction.
        SEQAN_ASSERT_EQ(1u, beginPositionH(s));
        SEQAN_ASSERT_EQ(4u, endPositionH(s));
        SEQAN_ASSERT_EQ(2u, beginPositionV(s));
        SEQAN_ASSERT_EQ(5u, endPositionV(s));
        SEQAN_ASSERT_EQ(-1, lowerDiagonal(s));
        SEQAN_ASSERT_EQ(-1, upperDiagonal(s));
        SEQAN_ASSERT_EQ(-1, beginDiagonal(s));
        SEQAN_ASSERT_EQ(-1, endDiagonal(s));
    }
}

// Test setters that are specific to Simple Seeds.
SEQAN_DEFINE_TEST(test_seeds_seed_simple_setters)
{
    using namespace seqan;

    // Define Seed type and declare a variable.
    typedef Seed<Simple> TSeed;
    TSeed s(1, 2, 3, 5);

    // Run setters.
    setBeginPositionH(s, 2);
    setBeginPositionV(s, 2);
    setEndPositionH(s, 4);
    setEndPositionV(s, 4);
    SEQAN_ASSERT_EQ(2u, beginPositionH(s));
    SEQAN_ASSERT_EQ(2u, beginPositionV(s));
    SEQAN_ASSERT_EQ(4u, endPositionH(s));
    SEQAN_ASSERT_EQ(4u, endPositionV(s));
    SEQAN_ASSERT_EQ(-2, lowerDiagonal(s));
    SEQAN_ASSERT_EQ(-1, upperDiagonal(s));
    SEQAN_ASSERT_EQ(0, beginDiagonal(s));
    SEQAN_ASSERT_EQ(0, endDiagonal(s));
}

SEQAN_BEGIN_TESTSUITE(test_seeds_seed_simple)
{
    SEQAN_CALL_TEST(test_seeds_seed_simple_constructors);
    SEQAN_CALL_TEST(test_seeds_seed_simple_setters);
}
SEQAN_END_TESTSUITE
