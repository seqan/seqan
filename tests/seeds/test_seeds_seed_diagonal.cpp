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
// Test the header seeds_seed_diagonal.h containing the SeedDiagonal class.
// ==========================================================================

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/stream.h>   // Required to print strings in tests.

#include <seqan/seeds.h>  // Include module under test.

// Test constructors of the SeedDiagonal class.  This also tests that
// all member variables are there.
SEQAN_DEFINE_TEST(test_seeds_seed_diagonal_constructors)
{
    using namespace seqan;

    typedef SeedDiagonal<int, int> TSeedDiagonal;

    // Default constructor.
    {
        TSeedDiagonal sd;
        SEQAN_ASSERT_EQ(0, sd.beginPositionH);
        SEQAN_ASSERT_EQ(0, sd.beginPositionV);
        SEQAN_ASSERT_EQ(0, sd.length);
    }
    // Only other constructor has all properties.
    {
        TSeedDiagonal sd(1, 2, 3);
        SEQAN_ASSERT_EQ(1, sd.beginPositionH);
        SEQAN_ASSERT_EQ(2, sd.beginPositionV);
        SEQAN_ASSERT_EQ(3, sd.length);
    }
}


// Test metafunctions of the SeedDiagonal class.
SEQAN_DEFINE_TEST(test_seeds_seed_diagonal_metafunctions)
{
    using namespace seqan;

    // Test the parametrization expected to be used most.
    {
        typedef SeedDiagonal<size_t, size_t> TSeedDiagonal;
        typedef Position<TSeedDiagonal>::Type TPosition;
        bool b1 = IsSameType<size_t, TPosition>::VALUE;
        SEQAN_ASSERT(b1);
        typedef Size<TSeedDiagonal>::Type TSize;
        bool b2 = IsSameType<size_t, TSize>::VALUE;
        SEQAN_ASSERT(b2);
    }
    // Test another parametrization.
    {
        typedef SeedDiagonal<double, int> TSeedDiagonal;
        typedef Position<TSeedDiagonal>::Type TPosition;
        bool b1 = IsSameType<double, TPosition>::VALUE;
        SEQAN_ASSERT(b1);
        typedef Size<TSeedDiagonal>::Type TSize;
        bool b2 = IsSameType<int, TSize>::VALUE;
        SEQAN_ASSERT(b2);
    }
}

SEQAN_BEGIN_TESTSUITE(test_seeds_seed_diagonal)
{
    // Tests for seed diagonals.
    SEQAN_CALL_TEST(test_seeds_seed_diagonal_constructors);
    SEQAN_CALL_TEST(test_seeds_seed_diagonal_metafunctions);
}
SEQAN_END_TESTSUITE
