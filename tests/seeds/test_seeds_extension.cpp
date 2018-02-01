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
// Tests for the header seeds_extension.h.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/stream.h>  // for printing seqan::String<>

#include <seqan/seeds.h>

template <typename TSeedSpec>
void
testSeedsExtensionMatchExtension(TSeedSpec const &)
{
    using namespace seqan;

    typedef Seed<TSeedSpec> TSeed;

    // Test extension to the left only.  Extension possible by 2.
    {
        DnaString database = "CCCCCCCCC";
        DnaString query = "AAACCCCCCAA";
        TSeed seed(4, 5, 2);
        extendSeed(seed, database, query, EXTEND_LEFT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(2, 3, 4), seed);
    }
    // Test extension to the left only.  No extension is possible.
    {
        DnaString database = "CCCCCCCCC";
        DnaString query = "AAAAACCCCAA";
        TSeed seed(4, 5, 2);
        extendSeed(seed, database, query, EXTEND_LEFT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(4, 5, 2), seed);
    }
    // Test extension to the right only.  Extension possible by 2.
    {
        DnaString database = "AAAAAAAAA";
        DnaString query = "AAAAAAAAACC";
        TSeed seed(5, 4, 2);
        extendSeed(seed, query, database, EXTEND_RIGHT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(5, 4, 4), seed);
    }
    // Test extension to the right only.  No extension is possible.
    {
        DnaString database = "AAAAAACC";
        DnaString query = "AAAAAAAAAA";
        TSeed seed(4, 5, 2);
        extendSeed(seed, database, query, EXTEND_RIGHT, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(4, 5, 2), seed);
    }
    // Test extension in both directions.  Extension possible by 2 in both directions.
    {
        DnaString database = "GGCCCCCCAA";
        DnaString query = "AAACCCCCCC";
        TSeed seed(4, 5, 2);
        extendSeed(seed, database, query, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(2, 3, 6), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the left.
    {
        DnaString database = "GGCCCCAAAA";
        DnaString query = "AAACCCCCCC";
        TSeed seed(4, 5, 2);
        extendSeed(seed, database, query, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(2, 3, 4), seed);
    }
    // Test extension in both directions.  Extension possible by 2 to the right.
    {
        DnaString database = "GGCCCCCCAA";
        DnaString query = "AAAAACCCCC";
        TSeed seed(4, 5, 2);
        extendSeed(seed, database, query, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(4, 5, 4), seed);
    }
    // Test extension in both directions.  No Extension possible in either direction.
    {
        DnaString database = "CCCCAACCC";
        DnaString query = "AAAAAAAAAAA";
        TSeed seed(4, 5, 2);
        extendSeed(seed, database, query, EXTEND_BOTH, MatchExtend());
        SEQAN_ASSERT_EQ(TSeed(4, 5, 2), seed);
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
        DnaString database = "AACCCCTTTGGTGAAAAA";
        DnaString query =    "AAACCCTTTGGGTTTTT";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(4, 4, 3);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 13u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 14u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 2);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -2);
    }
    { // Test 2
        DnaString database = "ttttcgatcgatgcttttt";
        DnaString query =    "aaaacgatcgatgc";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(4, 4, 9);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 3u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 14u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 3u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 15u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 1);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), 0);
    }
    { // Test 3
        DnaString database = "ttttcgatcgatgcttttt";
        DnaString query =    "aaaacgatcgatgc";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(5, 5, 7);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 3u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 14u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 3u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 15u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 1);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -1);
    }
    { // Test 4
        DnaString database = "ttttcgatcgatgcttttt";
        DnaString query =    "aaaacgatcgatgc";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(6, 6, 6);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 0, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 4u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 14u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 4u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 14u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), 0);
    }
    { // Test 5
        DnaString database = "ttttcgatcgatgc";
        DnaString query =        "cgatcgatgcaaaaaaaaa";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(5, 1, 5);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 11u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 3u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 14u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 5);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), 3);
    }
    { // Test 6
        DnaString database =       "cgatcgatgccaact";
        DnaString query = "aaaaaaaaacgatcgatgcaaaaaaaaa";
        TScoringScheme scoringScheme(2, -1, -1);
        TSeed seed(2, 11, 4);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 8u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 23u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 14u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), -7);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -10);
    }
    { // Test 7
        DnaString database = "ccccagctgatcgtttgcccccc";
        DnaString query =      "ttttagtgacgttttaaaaaa";
        TScoringScheme scoringScheme(1,-5,-5);
        TSeed seed(11, 9, 5);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 7, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 4u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 15u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 4u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 17u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 2);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), 0);
    }
    {
        DnaString database = "ccccccgtttgctagtcgacccc";
        DnaString query =    "aaaaaattttgcagtgatttt";
        TScoringScheme scoringScheme(1, -5, -5);
        TSeed seed(7, 7, 5);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 7, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 6u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 17u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 6u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 19u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 2);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), 0);
    }
    {
        DnaString database = "ttttagtgacgttttaaaaaa";
        DnaString query =  "ccccagctgatcgtttgcccccc";
        TScoringScheme scoringScheme(1, -5, -5);
        TSeed seed(9, 11, 5);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 7, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 4u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 17u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 4u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 15u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -2);
    }
    { // Test 10
        DnaString database =    "AGCAGCAAACAGTAAGCCAGCAGCCT";
        DnaString query =  "CTGACAGCGAGAAACAGTAACCAGCTAGCCT";
        TScoringScheme scoringScheme(1, -9, -9);
        TSeed seed(21, 26, 5);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme, 45, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 1u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 31u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 26u);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -9);
    }

    Blosum62 scoringScheme2;
    setScoreGap(scoringScheme2, -9);

    { // Test 11 (simple)
        Peptide database = "CSNNHKMMMMAAGGW";
        Peptide query =    "CSNNHKMMMMAAGGW";
        TSeed seed(6, 6, 3);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme2, 1, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionV(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 15u);
        SEQAN_ASSERT_EQ(beginPositionH(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 15u);
        SEQAN_ASSERT_EQ(beginDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(endDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), 0);
    }
    { // Test 12 (gaps)
        Peptide database = "CSNNHK""MMMMAAGGW";
        Peptide query =    "CSNNH" "MMMMAAGGW";
        TSeed seed(6, 5, 3);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme2, 10, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionH(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 15u);
        SEQAN_ASSERT_EQ(beginPositionV(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 14u);
        SEQAN_ASSERT_EQ(beginDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(endDiagonal(seed), 1);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 2);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -1);
    }
    { // Test 13 (good mismatches)
        Peptide database = "CSNNHKMMMMAAGGW";
        Peptide query =    "CSNNHRMLMLAAGGW";
        TSeed seed(6, 6, 3);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme2, 10, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionH(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 15u);
        SEQAN_ASSERT_EQ(beginPositionV(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 15u);
        SEQAN_ASSERT_EQ(beginDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(endDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 1);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -1);
    }
    { // Test 14 (bad mismatches)
        Peptide database = "CSNNHKMMMMAAGGW";
        Peptide query =    "WWNNHKMMMMAAGLV";
        TSeed seed(6, 6, 3);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme2, 10, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionH(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 15u);
        SEQAN_ASSERT_EQ(beginPositionV(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 15u);
        SEQAN_ASSERT_EQ(beginDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(endDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 1);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -1);
    }
    { // Test 15 (mixed mismatches)
        Peptide database = "CSNNH""KMMMMAAGGW";
        Peptide query =    "WSNNH" "MLMWAAGLV";
        TSeed seed(6, 5, 3);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme2, 10, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionH(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 15u);
        SEQAN_ASSERT_EQ(beginPositionV(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 14u);
        SEQAN_ASSERT_EQ(beginDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(endDiagonal(seed), 1);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 2);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -1);
    }
    { // Test 16
        Peptide database = "CSNNH""KMMMMA"  "GWW";
        Peptide query =    "WSNNH" "MLMWA""ARGWW";
        TSeed seed(6, 5, 3);

        extendSeed(seed, database, query, EXTEND_BOTH, scoringScheme2, 10, GappedXDrop());

        SEQAN_ASSERT_EQ(beginPositionH(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionH(seed), 14u);
        SEQAN_ASSERT_EQ(beginPositionV(seed), 0u);
        SEQAN_ASSERT_EQ(endPositionV(seed), 13u);
        SEQAN_ASSERT_EQ(beginDiagonal(seed), 0);
        SEQAN_ASSERT_EQ(endDiagonal(seed), 1);
        SEQAN_ASSERT_EQ(upperDiagonal(seed), 2);
        SEQAN_ASSERT_EQ(lowerDiagonal(seed), -1);
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

SEQAN_BEGIN_TESTSUITE(test_seeds_extension)
{
    // Tests for seed extension algorithms
    SEQAN_CALL_TEST(test_seeds_extension_match_extension_simple);
    SEQAN_CALL_TEST(test_seeds_extension_ungapped_xdrop_extension_simple);
    SEQAN_CALL_TEST(test_seeds_extension_gapped_xdrop_extension_simple);
    SEQAN_CALL_TEST(test_seeds_extension_match_extension_chained);
    SEQAN_CALL_TEST(test_seeds_extension_ungapped_xdrop_extension_chained);

    // Disabled the test for now.  Extension function contains a
    // force-failure assertion and instruction show to implement this.
    // See http://trac.mi.fu-berlin.de/seqan/ticket/344 for details.
    // 
    // TODO(holtgrew): Implement this.
    // 
    // SEQAN_CALL_TEST(test_seeds_extension_gapped_xdrop_extension_chained);
}
SEQAN_END_TESTSUITE
