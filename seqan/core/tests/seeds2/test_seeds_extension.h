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
// Tests for the header seeds_extension.h.
// ==========================================================================

#ifndef TEST_SEEDS_TEST_SEEDS_EXTENSION_H_
#define TEST_SEEDS_TEST_SEEDS_EXTENSION_H_

template <typename TSeedSpec>
void
testSeedsExtensionMatchExtension(TSeedSpec const &)
{
    using namespace seqan;

    typedef Seed<TSeedSpec> TSeed;

    // Test extension to the left only.  Extension possible by 2.
    {
        DnaString query = "AAACCCCCCAA";
        DnaString database = "CCCCCCCCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_LEFT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 4), seed);
    }
    // Test extension to the left only.  No extension is possible.
    {
        DnaString query = "AAAAACCCCAA";
        DnaString database = "CCCCCCCCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_LEFT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
    // Test extension to the right only.  Extension possible by 2.
    {
        DnaString query = "AAAAAAAAACC";
        DnaString database = "AAAAAAAAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_RIGHT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 4), seed);
    }
    // Test extension to the right only.  No extension is possible.
    {
        DnaString query = "AAAAAAAAAA";
        DnaString database = "AAAAAACC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_RIGHT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
    // Test extension in both directions.  Extension possible by 2 in both directions.
    {
        DnaString query = "AAACCCCCCC";
        DnaString database = "GGCCCCCCAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 6), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the left.
    {
        DnaString query = "AAACCCCCCC";
        DnaString database = "GGCCCCAAAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 4), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the right.
    {
        DnaString query = "AAAAACCCCC";
        DnaString database = "GGCCCCCCAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 4), seed);
    }
    // Test extension in both directions.  No Extension possible in either direction.
    {
        DnaString query = "AAAAAAAAAAA";
        DnaString database = "CCCCAACCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
}


template <typename TSeedSpec>
void
testSeedsExtensionUnGappedXDropExtension(TSeedSpec const &)
{
    using namespace seqan;

    typedef Seed<TSeedSpec> TSeed;

    // Test extension to the left only.  Extension possible by 2 over a gap of length 1.
    {
        DnaString query = "AAACACCCCAA";
        DnaString database = "CCCCCCCCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_LEFT, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 4), seed);
    }
    // Test extension to the left only.  No extension is possible.
    {
        DnaString query = "AAAAACCACAA";
        DnaString database = "CCCCCCCCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_LEFT, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
    // Test extension to the right only.  Extension possible by 2 over a gap of length 1.
    {
        DnaString query = "AAAAAAACACC";
        DnaString database = "AAAAAAAAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_RIGHT, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 4), seed);
    }
    // Test extension to the right only.  No extension is possible.
    {
        DnaString query = "AAAAAAAAAA";
        DnaString database = "AAAAAACC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_RIGHT, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
    // Test extension in both directions.  Extension possible by 2 in both directions over a gap of length 1 on each side.
    {
        DnaString query = "AAACCCCCCC";
        DnaString database = "GGCCCCCCAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 6), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the left over a gap of length 1.
    {
        DnaString query = "AAACCCCCCC";
        DnaString database = "GGCCCCAAAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(3, 2, 4), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the right over a gap of length 1.
    {
        DnaString query = "AAAAACCCCC";
        DnaString database = "GGCCCCCCAA";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 4), seed);
    }
    // Test extension in both directions.  No Extension possible in either direction.
    {
        DnaString query = "AAAAAAAAAAA";
        DnaString database = "CCCCAACCC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_BOTH, Score<int, EditDistance>(), 2, UnGappedXDrop());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 2), seed);
    }
}


template <typename TSeedSpec>
void
testSeedsExtensionGappedXDropExtension(TSeedSpec const &)
{
    using namespace seqan;

    typedef Seed<TSeedSpec> TSeed;
	typedef Score<int, Simple> TScoringScheme;

    // Each of the following blocks contains a test case for the gapped X-drop
    // seed extension.
    
    { // Test 1
        DnaString query =    "AAACCCTTTGGGTTTTT";
        DnaString database = "AACCCCTTTGGTGAAAAA";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(4, 4, 3);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(getBeginDim0(seed), 0u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 13u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 0u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 14u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), 2);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), -2);
    }
    { // Test 2
        DnaString query =    "aaaacgatcgatgc";
        DnaString database = "ttttcgatcgatgcttttt";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(4, 4, 9);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(getBeginDim0(seed), 3u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 14u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 3u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 15u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), 1);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), 0);
    }
    { // Test 3
        DnaString query =    "aaaacgatcgatgc";
        DnaString database = "ttttcgatcgatgcttttt";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(5, 5, 7);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(getBeginDim0(seed), 3u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 14u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 3u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 15u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), 1);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), -1);
    }
    { // Test 4
        DnaString query =    "aaaacgatcgatgc";
        DnaString database = "ttttcgatcgatgcttttt";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(6, 6, 6);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 0, GappedXDrop());

        SEQAN_ASSERT_EQ(getBeginDim0(seed), 4u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 14u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 4u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 14u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), 0);
    }
    { // Test 5
        DnaString query =        "cgatcgatgcaaaaaaaaa";
        DnaString database = "ttttcgatcgatgc";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(1, 5, 5);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(getBeginDim0(seed), 0u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 11u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 3u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 14u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), 5);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), 3);
    }
    { // Test 6
        DnaString query = "aaaaaaaaacgatcgatgcaaaaaaaaa";
        DnaString database =       "cgatcgatgccaact";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(11,2,4);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(getBeginDim0(seed), 8u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 23u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 0u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 14u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), -7);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), -10);
    }
    { // Test 7
        DnaString query =      "ttttagtgacgttttaaaaaa";
        DnaString database = "ccccagctgatcgtttgcccccc";
        TScoringScheme scoringScheme(1,-5,-5);
        TSeed seed(9,11,5);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 7, GappedXDrop());
        
        SEQAN_ASSERT_EQ(getBeginDim0(seed), 4u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 15u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 4u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 17u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), 2);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), 0);
    }
    {
        DnaString query =    "aaaaaattttgcagtgatttt";
        DnaString database = "ccccccgtttgctagtcgacccc";
        TScoringScheme scoringScheme(1, -5, -5);
        TSeed seed(7, 7, 5);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 7, GappedXDrop());
        
        SEQAN_ASSERT_EQ(getBeginDim0(seed), 6u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 17u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 6u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 19u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), 2);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), 0);
    }
    {
        DnaString query =  "ccccagctgatcgtttgcccccc";
        DnaString database = "ttttagtgacgttttaaaaaa";
        TScoringScheme scoringScheme(1, -5, -5);
        TSeed seed(11, 9, 5);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 7, GappedXDrop());

        SEQAN_ASSERT_EQ(getBeginDim0(seed), 4u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 17u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 4u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 15u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), -2);
    }
    { // Test 10
        DnaString query =  "CTGACAGCGAGAAACAGTAACCAGCTAGCCT";
        DnaString database =    "AGCAGCAAACAGTAAGCCAGCAGCCT";
        TScoringScheme scoringScheme(1, -9, -9);
        TSeed seed(26, 21, 5);

        extendSeed(seed, query, database, EXTEND_BOTH, scoringScheme, 45, GappedXDrop());

        SEQAN_ASSERT_EQ(getBeginDim0(seed), 1u);
        SEQAN_ASSERT_EQ(getBeginDim1(seed), 0u);
        SEQAN_ASSERT_EQ(getEndDim0(seed), 31u);
        SEQAN_ASSERT_EQ(getEndDim1(seed), 26u);
        SEQAN_ASSERT_EQ(getUpperDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(getLowerDiagonal(seed), -9);
    }
}


// Test the seed extension algorithm with match extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_match_extension_simple)
{
    using namespace seqan;
    testSeedsExtensionMatchExtension(Simple());
}

// Test the seed extension algorithm with ungapped x-drop extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_ungapped_xdrop_extension_simple)
{
    using namespace seqan;
    testSeedsExtensionUnGappedXDropExtension(Simple());
}

// Test the seed extension algorithm with gapped x-drop extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_gapped_xdrop_extension_simple)
{
    using namespace seqan;
    testSeedsExtensionGappedXDropExtension(Simple());
}

// Test the seed extension algorithm with ungapped x-drp extension for chained seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_match_extension_chained)
{
    using namespace seqan;
    testSeedsExtensionMatchExtension(ChainedSeed());
}

// Test the seed extension algorithm with ungapped x-drop extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_ungapped_xdrop_extension_chained)
{
    using namespace seqan;
    testSeedsExtensionUnGappedXDropExtension(ChainedSeed());
}

// Test the seed extension algorithm with gapped x-drop extension for simple seeds.
SEQAN_DEFINE_TEST(test_seeds_extension_gapped_xdrop_extension_chained)
{
    using namespace seqan;
    testSeedsExtensionGappedXDropExtension(ChainedSeed());
}

#endif  // TEST_SEEDS_TEST_SEEDS_EXTENSION_H_
