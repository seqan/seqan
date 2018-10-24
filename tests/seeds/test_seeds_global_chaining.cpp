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
// Tests for the header seeds_global_chaining.h.
// ==========================================================================

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/stream.h>   // Required to print strings in tests.

#include <seqan/seeds.h>  // Include module under test.

// -----------------------------------------------------------------------------
// Test global chaining weighting the seeds by their length only (Simple Seed)
// -----------------------------------------------------------------------------

using namespace seqan;

typedef SeedSet<Seed<Simple>, Unordered> TSeedSet;
typedef Value<TSeedSet>::Type TSeed;
typedef String<TSeed> TSeedChain;

SEQAN_DEFINE_TEST(test_seeds_global_chaining_one_seed)
{
    TSeedSet seedSet;
    addSeed(seedSet, TSeed(1, 2, 3), Single());

    TSeedChain result;
    chainSeedsGlobally(result, seedSet, SparseChaining());

    SEQAN_ASSERT_EQ(1u, length(result));
    SEQAN_ASSERT_EQ(TSeed(1, 2, 3), front(result));
}

SEQAN_DEFINE_TEST(test_seeds_global_chaining_two_seeds)
{
    TSeedSet seedSet;
    addSeed(seedSet, TSeed(1, 2, 3), Single());
    addSeed(seedSet, TSeed(4, 5, 6), Single());

    TSeedChain result;
    chainSeedsGlobally(result, seedSet, SparseChaining());

    SEQAN_ASSERT_EQ(2u, length(result));
    SEQAN_ASSERT_EQ(TSeed(1, 2, 3), result[0]);
    SEQAN_ASSERT_EQ(TSeed(4, 5, 6), result[1]);
}


// Test with two seeds, only first one is part of the chain.
SEQAN_DEFINE_TEST(test_seeds_global_chaining_first_seed)
{
    TSeedSet seedSet;
    addSeed(seedSet, TSeed(1, 2, 3), Single());
    addSeed(seedSet, TSeed(2, 1, 2), Single());

    TSeedChain result;
    chainSeedsGlobally(result, seedSet, SparseChaining());

    SEQAN_ASSERT_EQ(1u, length(result));
    SEQAN_ASSERT_EQ(TSeed(1, 2, 3), front(result));
}

// A bit larger example.
SEQAN_DEFINE_TEST(test_seeds_global_chaining_four_seeds)
{
    TSeedSet seedSet;
    addSeed(seedSet, TSeed(0, 0, 2), Single());
    addSeed(seedSet, TSeed(3, 5, 2), Single());
    addSeed(seedSet, TSeed(4, 2, 3), Single());
    addSeed(seedSet, TSeed(9, 9, 2), Single());

    TSeedChain result;
    chainSeedsGlobally(result, seedSet, SparseChaining());

    SEQAN_ASSERT_EQ(3u, length(result));
    SEQAN_ASSERT_EQ(TSeed(0, 0, 2), result[0]);
    SEQAN_ASSERT_EQ(TSeed(4, 2, 3), result[1]);
    SEQAN_ASSERT_EQ(TSeed(9, 9, 2), result[2]);
}

// Another non-trivial example that caused a problem earlier.  The
// seeds are all overlapping here.
SEQAN_DEFINE_TEST(test_seeds_global_chaining_overlapping_seeds)
{
    TSeedSet seedSet;
    addSeed(seedSet, TSeed(0, 93, 281, 342), Single());
    addSeed(seedSet, TSeed(3, 237, 127, 364), Single());
    addSeed(seedSet, TSeed(3, 284, 86, 368), Single());
    addSeed(seedSet, TSeed(5, 146, 239, 374), Single());

    TSeedChain result;
    chainSeedsGlobally(result, seedSet, SparseChaining());

    SEQAN_ASSERT_EQ(1u, length(result));
    SEQAN_ASSERT_EQ(TSeed(0, 93, 281, 342), result[0]);
}

// Issue #2082
SEQAN_DEFINE_TEST(test_seeds_global_chaining_issue_2082)
{
    {
        TSeedSet seedSet;

        addSeed(seedSet, TSeed(0, 0, 3), Single());
        addSeed(seedSet, TSeed(2, 3, 2), Single());

        TSeedChain result;
        chainSeedsGlobally(result, seedSet, SparseChaining());

        SEQAN_ASSERT_EQ(1u, length(result));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 3), result[0]);
    }

    {
        TSeedSet seedSet;

        addSeed(seedSet, TSeed(0, 0, 100), Single());
        addSeed(seedSet, TSeed(95, 95, 10), Single());

        TSeedChain result;
        chainSeedsGlobally(result, seedSet, SparseChaining());

        SEQAN_ASSERT_EQ(1u, length(result));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 100), result[0]);
    }
}

// Issue #2318
SEQAN_DEFINE_TEST(test_seeds_global_chaining_issue_2318)
{
    TSeedSet simple;
    addSeed(simple, TSeed(1,  5, 10, 15), Single());
    addSeed(simple, TSeed(5,  9, 40, 34), Single());
    addSeed(simple, TSeed(40, 10, 45, 17), Single());

    TSeedChain result;
    chainSeedsGlobally(result, simple, SparseChaining());

    SEQAN_ASSERT_EQ(1u, length(result));
    SEQAN_ASSERT_EQ(TSeed(5,  9, 40, 34), result[0]);
}

// Another test case, without deleting seeds (from list L) during the chaining algorithm
SEQAN_DEFINE_TEST(test_seeds_global_chaining_no_deleting_intermediate_solution)
{
    TSeedSet simple;
    addSeed(simple, TSeed(27, 21, 37, 24), Single());
    addSeed(simple, TSeed(25, 14, 33, 19), Single());
    addSeed(simple, TSeed(15, 7, 21, 11), Single());
    addSeed(simple, TSeed(13, 13, 22, 14), Single());
    addSeed(simple, TSeed(2,  8, 9, 12), Single());
    addSeed(simple, TSeed(1,  1, 7, 5), Single());

    TSeedChain result;
    chainSeedsGlobally(result, simple, SparseChaining());

    SEQAN_ASSERT_EQ(3u, length(result));
    SEQAN_ASSERT_EQ(TSeed(2,  8, 9, 12), result[0]);
    SEQAN_ASSERT_EQ(TSeed(13, 13, 22, 14), result[1]);
    SEQAN_ASSERT_EQ(TSeed(27, 21, 37, 24), result[2]);
}

// Another test case, with deleting seeds (from list L) during the chaining algorithm
SEQAN_DEFINE_TEST(test_seeds_global_chaining_deleting_intermediate_solution)
{
    TSeedSet simple;
    addSeed(simple, TSeed(2, 16, 10, 21), Single());
    addSeed(simple, TSeed(7, 12, 19, 19), Single());
    addSeed(simple, TSeed(12, 8, 28, 17), Single());
    addSeed(simple, TSeed(22, 3, 39, 12), Single());

    TSeedChain result;
    chainSeedsGlobally(result, simple, SparseChaining());

    SEQAN_ASSERT_EQ(1u, length(result));
    SEQAN_ASSERT_EQ(TSeed(22, 3, 39, 12), result[0]);
}

// Another test case, with bordering rectangles (each seed touches its neighbours)
SEQAN_DEFINE_TEST(test_seeds_global_chaining_bordering_rectangles)
{
    TSeedSet simple;
    addSeed(simple, TSeed(1, 1, 3, 3), Single());
    addSeed(simple, TSeed(3, 1, 6, 3), Single());
    addSeed(simple, TSeed(6, 1, 10, 3), Single());
    addSeed(simple, TSeed(1, 3, 3, 6), Single());
    addSeed(simple, TSeed(3, 3, 6, 6), Single());
    addSeed(simple, TSeed(6, 3, 10, 6), Single());
    addSeed(simple, TSeed(1, 6, 3, 10), Single());
    addSeed(simple, TSeed(3, 6, 6, 10), Single());
    addSeed(simple, TSeed(6, 6, 10, 10), Single());

    TSeedChain result;
    chainSeedsGlobally(result, simple, SparseChaining());

    SEQAN_ASSERT_EQ(3u, length(result));
    SEQAN_ASSERT_EQ(TSeed(1, 1, 3, 3), result[0]);
    SEQAN_ASSERT_EQ(TSeed(3, 3, 6, 6), result[1]);
    SEQAN_ASSERT_EQ(TSeed(6, 6, 10, 10), result[2]);
}

// Another test case, with rectangles on the same height (y axis)
SEQAN_DEFINE_TEST(test_seeds_global_chaining_same_height)
{
    {
        TSeedSet simple;
        addSeed(simple, TSeed(1, 1, 3, 3), Single());

        addSeed(simple, TSeed(3, 6, 5,  8), Single());
        addSeed(simple, TSeed(5, 5, 8,  8), Single());
        addSeed(simple, TSeed(8, 3, 13, 8), Single());

        addSeed(simple, TSeed(13, 8, 15, 10), Single());

        TSeedChain result;
        chainSeedsGlobally(result, simple, SparseChaining());

        SEQAN_ASSERT_EQ(3u, length(result));
        SEQAN_ASSERT_EQ(TSeed(1, 1, 3, 3), result[0]);
        SEQAN_ASSERT_EQ(TSeed(8, 3, 13, 8), result[1]);
        SEQAN_ASSERT_EQ(TSeed(13, 8, 15, 10), result[2]);
    }
    {
        TSeedSet simple;
        addSeed(simple, TSeed(1, 1, 3, 3), Single());

        addSeed(simple, TSeed(3,  3, 8,  8), Single());
        addSeed(simple, TSeed(8,  6, 10, 8), Single());
        addSeed(simple, TSeed(10, 5, 13, 8), Single());

        addSeed(simple, TSeed(13, 8, 15, 10), Single());

        TSeedChain result;
        chainSeedsGlobally(result, simple, SparseChaining());

        SEQAN_ASSERT_EQ(3u, length(result));
        SEQAN_ASSERT_EQ(TSeed(1,  1, 3,  3), result[0]);
        SEQAN_ASSERT_EQ(TSeed(3,  3, 8,  8), result[1]);
        SEQAN_ASSERT_EQ(TSeed(13, 8, 15, 10), result[2]);
    }
    {
        TSeedSet simple;
        addSeed(simple, TSeed(1, 1, 3, 3), Single());

        addSeed(simple, TSeed(3,  5, 6,  8), Single());
        addSeed(simple, TSeed(6,  3, 11, 8), Single());
        addSeed(simple, TSeed(11, 6, 13, 8), Single());

        addSeed(simple, TSeed(13, 8, 15, 10), Single());

        TSeedChain result;
        chainSeedsGlobally(result, simple, SparseChaining());

        SEQAN_ASSERT_EQ(3u, length(result));
        SEQAN_ASSERT_EQ(TSeed(1,  1, 3,  3), result[0]);
        SEQAN_ASSERT_EQ(TSeed(6,  3, 11, 8), result[1]);
        SEQAN_ASSERT_EQ(TSeed(13, 8, 15, 10), result[2]);
    }
}

SEQAN_BEGIN_TESTSUITE(test_seeds_global_chaining)
{
    // Test global chaining of seeds.
    SEQAN_CALL_TEST(test_seeds_global_chaining_one_seed);
    SEQAN_CALL_TEST(test_seeds_global_chaining_two_seeds);
    SEQAN_CALL_TEST(test_seeds_global_chaining_first_seed);
    SEQAN_CALL_TEST(test_seeds_global_chaining_four_seeds);
    SEQAN_CALL_TEST(test_seeds_global_chaining_overlapping_seeds);
    SEQAN_CALL_TEST(test_seeds_global_chaining_issue_2082);
    SEQAN_CALL_TEST(test_seeds_global_chaining_issue_2082);
    SEQAN_CALL_TEST(test_seeds_global_chaining_issue_2318);
    SEQAN_CALL_TEST(test_seeds_global_chaining_no_deleting_intermediate_solution);
    SEQAN_CALL_TEST(test_seeds_global_chaining_deleting_intermediate_solution);
    SEQAN_CALL_TEST(test_seeds_global_chaining_bordering_rectangles);
    SEQAN_CALL_TEST(test_seeds_global_chaining_same_height);
}

SEQAN_END_TESTSUITE
