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
// Test the specialization Chained Seed.
// ==========================================================================

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/stream.h>   // Required to print strings in tests.

#include <seqan/seeds.h>  // Include module under test.

struct TestSmallSeedConfig
{
    typedef unsigned TPosition;
    typedef unsigned TSize;
    typedef int TDiagonal;
    typedef seqan::True THasScore;
    typedef int TScoreValue;
};

// Test assignment of chained seeds.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_assign)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;

    TSeed s(0, 0, 3);

    TSeed s2 = s;
    s2 = s;
    assign(s2, s);
}

// Test the metafunctions of the ChainedSeed specialization.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_metafunctions)
{
    using namespace seqan;

    // Test with the default configuration.
    {
        typedef Seed<ChainedSeed> TSeed;
        typedef Value<TSeed>::Type TSeedDiagonal;
        bool b;
        b = IsSameType<size_t, Position<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<size_t, Size<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<size_t, Position<TSeed>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<size_t, Size<TSeed>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<MakeSigned_<size_t>::Type, Diagonal<TSeed>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<int, SeedScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Test with other specialization.
    {
        typedef Seed<ChainedSeed, TestSmallSeedConfig> TSeed;
        typedef Value<TSeed>::Type TSeedDiagonal;
        bool b;
        b = IsSameType<unsigned, Position<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<unsigned, Size<TSeedDiagonal>::Type>::VALUE;
        SEQAN_ASSERT(b);

        b = IsSameType<unsigned, Position<TSeed>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<unsigned, Size<TSeed>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<int, Diagonal<TSeed>::Type>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<int, SeedScore<TSeed>::Type>::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Test the front() and back() functions for Chained Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_front_back)
{
    using namespace seqan;
    
    typedef Seed<ChainedSeed> TSeed;
    typedef Value<TSeed>::Type TSeedDiagonal;

    {
        TSeed s(1, 3, 4);
        SEQAN_ASSERT_EQ(TSeedDiagonal(1, 3, 4), front(s));
        SEQAN_ASSERT_EQ(TSeedDiagonal(1, 3, 4), back(s));
    }
    {
        TSeed const cs(1, 3, 4);
        SEQAN_ASSERT_EQ(TSeedDiagonal(1, 3, 4), front(cs));
        SEQAN_ASSERT_EQ(TSeedDiagonal(1, 3, 4), back(cs));
    }
}

// Test the appendDiagonal() function for Chained Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_append_diagonal)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<TSeed>::Type TSeedDiagonal;
    TSeed s(1, 3, 4);

    SEQAN_ASSERT_EQ(1u, beginPositionH(s));
    SEQAN_ASSERT_EQ(3u, beginPositionV(s));
    SEQAN_ASSERT_EQ(5u, endPositionH(s));
    SEQAN_ASSERT_EQ(7u, endPositionV(s));
    SEQAN_ASSERT_EQ(-2, beginDiagonal(s));
    SEQAN_ASSERT_EQ(-2, endDiagonal(s));
    SEQAN_ASSERT_EQ(1u, length(s));

    appendDiagonal(s, TSeedDiagonal(5, 7, 3));

    SEQAN_ASSERT_EQ(1u, beginPositionH(s));
    SEQAN_ASSERT_EQ(3u, beginPositionV(s));
    SEQAN_ASSERT_EQ(8u, endPositionH(s));
    SEQAN_ASSERT_EQ(10u, endPositionV(s));
    SEQAN_ASSERT_EQ(-2, beginDiagonal(s));
    SEQAN_ASSERT_EQ(-2, endDiagonal(s));
    SEQAN_ASSERT_EQ(2u, length(s));
}

// Test the truncateDiagonal() function for Chained Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_truncate_diagonals)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<TSeed>::Type TSeedDiagonal;
    TSeed s(1, 3, 4);

    appendDiagonal(s, TSeedDiagonal(5, 7, 3));
    SEQAN_ASSERT_EQ(2u, length(s));
    appendDiagonal(s, TSeedDiagonal(10, 10, 3));
    SEQAN_ASSERT_EQ(3u, length(s));
    typedef Iterator<TSeed, Standard>::Type TIterator;
    TIterator it = begin(s);
    ++it;
    truncateDiagonals(s, it);
    SEQAN_ASSERT_EQ(1u, length(s));
    TSeedDiagonal const diag = *begin(s, Standard());
    SEQAN_ASSERT(diag == TSeedDiagonal(1, 3, 4));
}

// Test the begin/end functions for chained seeds.
SEQAN_DEFINE_TEST(test_seeds_seed_chained_iterators)
{
    using namespace seqan;

    typedef Seed<ChainedSeed> TSeed;
    typedef Value<TSeed>::Type TSeedDiagonal;

    TSeed s(1, 2, 3);
    appendDiagonal(s, TSeedDiagonal(4, 5, 3));

    {  // non-const seed
        typedef Iterator<TSeed, Standard>::Type TIterator;
        TIterator it = begin(s);
        SEQAN_ASSERT_EQ(1u, it->beginPositionH);
        SEQAN_ASSERT_EQ(2u, it->beginPositionV);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_EQ(4u, it->beginPositionH);
        SEQAN_ASSERT_EQ(5u, it->beginPositionV);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT(it == end(s));
    }
    {  // const seed
        TSeed const & cs = s;
        typedef Iterator<TSeed const, Standard>::Type TIterator;
        TIterator it = begin(cs);
        SEQAN_ASSERT_EQ(1u, it->beginPositionH);
        SEQAN_ASSERT_EQ(2u, it->beginPositionV);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT_EQ(4u, it->beginPositionH);
        SEQAN_ASSERT_EQ(5u, it->beginPositionV);
        SEQAN_ASSERT_EQ(3u, it->length);
        ++it;
        SEQAN_ASSERT(it == end(cs));
    }
}

SEQAN_BEGIN_TESTSUITE(test_seeds_seed_chained)
{
    SEQAN_CALL_TEST(test_seeds_seed_chained_assign);
    SEQAN_CALL_TEST(test_seeds_seed_chained_metafunctions);
    SEQAN_CALL_TEST(test_seeds_seed_chained_append_diagonal);
    SEQAN_CALL_TEST(test_seeds_seed_chained_truncate_diagonals);
    SEQAN_CALL_TEST(test_seeds_seed_chained_iterators);
    SEQAN_CALL_TEST(test_seeds_seed_chained_front_back);
}
SEQAN_END_TESTSUITE
