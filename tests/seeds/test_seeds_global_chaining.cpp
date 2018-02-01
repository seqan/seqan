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


// Test global chaining weighting the seeds by their length only.
SEQAN_DEFINE_TEST(test_seeds_global_chaining_sparse_length)
{
    using namespace seqan;

    typedef SeedSet<Seed<Simple>, Unordered> TSeedSet;
    typedef Value<TSeedSet>::Type TSeed;
    typedef String<TSeed> TSeedChain;

    // Test with one seed.
    {
        TSeedSet seedSet;
        addSeed(seedSet, TSeed(1, 2, 3), Single());

        TSeedChain result;
        chainSeedsGlobally(result, seedSet, SparseChaining());

        SEQAN_ASSERT_EQ(1u, length(result));
        SEQAN_ASSERT_EQ(TSeed(1, 2, 3), front(result));
    }
    // Test with two seeds, both are part of the chain.
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

    { // Issue #2082
       TSeedSet seedSet;

       addSeed(seedSet, TSeed(0, 0, 3), Single());
       addSeed(seedSet, TSeed(2, 3, 2), Single());

       TSeedChain result;
       chainSeedsGlobally(result, seedSet, SparseChaining());

       SEQAN_ASSERT_EQ(1u, length(result));
       SEQAN_ASSERT_EQ(TSeed(0, 0, 3), result[0]);
   }

   { // Issue #2082
        TSeedSet seedSet;

        addSeed(seedSet, TSeed(0, 0, 100), Single());
        addSeed(seedSet, TSeed(95, 95, 10), Single());

        TSeedChain result;
        chainSeedsGlobally(result, seedSet, SparseChaining());

        SEQAN_ASSERT_EQ(1u, length(result));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 100), result[0]);
    }
}

SEQAN_BEGIN_TESTSUITE(test_seeds_global_chaining)
{
    // Test global chaining of seeds.
    SEQAN_CALL_TEST(test_seeds_global_chaining_sparse_length);
}
SEQAN_END_TESTSUITE
