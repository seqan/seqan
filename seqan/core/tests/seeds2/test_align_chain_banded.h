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
// Author: Carsten Kemena <carsten.kemena@crg.es>
// ==========================================================================
// Tests for the banded chain alignment.  Contains code that is based on
// the test code by Carsten Kemena.
// ==========================================================================

#ifndef TEST_SEEDS_TEST_ALIGN_CHAIN_BANDED_H_
#define TEST_SEEDS_TEST_ALIGN_CHAIN_BANDED_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test overlap computation functions.
SEQAN_DEFINE_TEST(test_align_chain_banded_compute_upper_left_overlap)
{
    using namespace seqan;

    typedef Seed<Simple> TSimpleSeed;
    typedef Score<int, Simple> TScoringScheme;

    CharString sequence0 = "unimportant";
    CharString sequence1 = "also unimportant";
    AlignmentChain_<CharString, TScoringScheme, NeedlemanWunsch> alignmentChain(2, TScoringScheme(), sequence0, sequence1);

    // Test with a seed where start/end = lower/upper diagonal
    {
        TSimpleSeed seed(2, 1, 4, 5);
        // Diagonals: start/lower = -1, end/upper = 1;

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeUpperLeftOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(2u, overlap0);
        SEQAN_ASSERT_EQ(4u, overlap1);
    }
    // Test with a seed where start/end = upper/lower diagonal
    {
        TSimpleSeed seed(1, 2, 5, 4);
        // Diagonals: start/upper = 1, end/lower = -1;

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeUpperLeftOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(4u, overlap0);
        SEQAN_ASSERT_EQ(2u, overlap1);
    }
    // Test with a seed where start = end = lower = upper diagonal
    {
        TSimpleSeed seed(1, 1, 5, 5);
        // Diagonals: all are = 0

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeUpperLeftOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(3u, overlap0);
        SEQAN_ASSERT_EQ(3u, overlap1);
    }
    // Test with a seed where {start, end} != {lower, upper diagonal}.
    {
        TSimpleSeed seed(1, 2, 5, 4);
        // Diagonals: start = 1, end = -1;
        setUpperDiagonal(seed, 3);
        setLowerDiagonal(seed, -2);

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeUpperLeftOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(4u, overlap0);
        SEQAN_ASSERT_EQ(2u, overlap1);
    }
}


// Test overlap computation functions.
SEQAN_DEFINE_TEST(test_align_chain_banded_compute_lower_right_overlap)
{
    using namespace seqan;

    typedef Seed<Simple> TSimpleSeed;
    typedef Score<int, Simple> TScoringScheme;

    CharString sequence0 = "unimportant";
    CharString sequence1 = "also unimportant";
    AlignmentChain_<CharString, TScoringScheme, NeedlemanWunsch> alignmentChain(2, TScoringScheme(), sequence0, sequence1);

    // Test with a seed where start/end = lower/upper diagonal
    {
        TSimpleSeed seed(2, 1, 4, 5);
        // Diagonals: start/lower = -1, end/upper = 1;

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeLowerRightOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(2u, overlap0);
        SEQAN_ASSERT_EQ(4u, overlap1);
    }
    // Test with a seed where start/end = upper/lower diagonal
    {
        TSimpleSeed seed(1, 2, 5, 4);
        // Diagonals: start/upper = 1, end/lower = -1;

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeLowerRightOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(4u, overlap0);
        SEQAN_ASSERT_EQ(2u, overlap1);
    }
    // Test with a seed where start = end = lower = upper diagonal
    {
        TSimpleSeed seed(1, 1, 5, 5);
        // Diagonals: all are = 0

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeLowerRightOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(3u, overlap0);
        SEQAN_ASSERT_EQ(3u, overlap1);
    }
    // Test with a seed where {start, end} != {lower, upper diagonal}.
    {
        TSimpleSeed seed(1, 2, 5, 4);
        // Diagonals: start = 1, end = -1;
        setUpperDiagonal(seed, 3);
        setLowerDiagonal(seed, -2);

        // Compute overlaps.
        unsigned overlap0, overlap1;
        _computeLowerRightOverlap(overlap0, overlap1, seed, alignmentChain);
        SEQAN_ASSERT_EQ(4u, overlap0);
        SEQAN_ASSERT_EQ(2u, overlap1);
    }
}

// Test banded alignment algorithm around a chain.  Linear gap cost
// case.
SEQAN_DEFINE_TEST(test_align_chain_banded_align_linear)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    // Test alignment with chain of length 1 with various bandwidths
    // to check that too large bandwidths do not do any harm.
    for (unsigned k = 0; k < 6; ++k)
    {
        CharString sequence0 = "NNNAAANNN";
        CharString sequence1 =  "NCAAACNN";  // TODO(holtgrew): Switch back again.
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(3, 2, 6, 5));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        int result = bandedChainAlignment(alignment, seedChain, k, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ_MSG(result, 8, "k == %u", k);
    }

    // Test on whole strings.
    {
        // Resulting alignment should be something like this (seeds
        // are marked by < and >).
        //
        //     > <    > <    >  <
        //   GGCGATNNNCAT--GGCACA
        //   --CGA-ATCCATCCCACACA
        CharString sequence0 = "GGCGATNNNCATGGCACA";
        CharString sequence1 = "CGAATCCATCCCACACA";  // TODO(holtgrew): Switch back again.
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(2, 0, 6, 5));
        appendValue(seedChain, TSeed(9, 6, 12, 9));
        appendValue(seedChain, TSeed(14, 11, 16, 17));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        int result = bandedChainAlignment(alignment, seedChain, 2, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        //
        // Note that leading and trailing gaps are not stored, so the
        // start and end gaps do not appear in the strings.
        SEQAN_ASSERT(row(alignment, 0) == "GGCGATNNNCAT--GGCACA"/*--*/);
        SEQAN_ASSERT(row(alignment, 1) == /*--*/"CGA-ATCCATCCCACACA");
    }
    // Test on infixes.
    {
        // The test data is the same as above but with a T and a TT prepended.
        CharString sequence0 = "TGGCGATNNNCATGGCACA";
        CharString sequence1 = "TTCGAATCCATCCCACACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed > seedChain;
        appendValue(seedChain, TSeed(2, 0, 6, 5));
        appendValue(seedChain, TSeed(9, 6, 12, 9));
        appendValue(seedChain, TSeed(14, 11, 16, 17));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        // TODO(holtgrew): This infix assignment for alignments is a bit creepy, maybe one of the too many shortcuts?
        assignSource(row(alignment, 0), sequence0);
        setClippedBeginPosition(row(alignment, 0), 1);
        setClippedEndPosition(row(alignment, 0), length(sequence0));
        assignSource(row(alignment, 1), sequence1);
        setClippedBeginPosition(row(alignment, 1), 2);
        setClippedEndPosition(row(alignment, 1), length(sequence1));

        int result = bandedChainAlignment(alignment, seedChain, 1, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        //
        // Note that leading and trailing gaps are not stored, so the
        // start and end gaps do not appear in the strings.
        SEQAN_ASSERT(row(alignment, 0) == "GGCGATNNNCAT--GGCACA"/*--*/);
        SEQAN_ASSERT(row(alignment, 1) == /*--*/"CGA-ATCCATCCCACACA");
    }
}


// Test banded alignment algorithm around a chain.  Affine gap cost
// case.
SEQAN_DEFINE_TEST(test_align_chain_banded_align_affine)
{
    using namespace seqan;
    typedef Seed<Simple> TSeed;

    // Test alignment with chain of length 1 with various bandwidths
    // to check that too large bandwidths do not do any harm.
    for (unsigned k = 0; k < 6; ++k)
    {
        CharString sequence0 = "NNNAAANNN";
        CharString sequence1 =  "NCAAACNN";  // TODO(holtgrew): Switch back again.
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(3, 2, 6, 5));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        int result = _bandedChainAlignment(alignment, seedChain, k, scoringScheme, AlignConfig<false, false, false, false>(), Gotoh());
        SEQAN_ASSERT_EQ_MSG(result, 8, "k == %u", k);

        // Compare alignment rows.
        //
        // Note that leading and trailing gaps are not stored, so the
        // start and end gaps do not appear in the strings.
        // SEQAN_ASSERT_MSG(row(alignment, 0) == "NNNAAANNN", "k == %u", k);
        // SEQAN_ASSERT_MSG(row(alignment, 1) == /*-*/"NCAAACNN", "k == %u", k);
    }

    // Test on whole strings.
    {
        CharString query = "ACGTCCTCGTACACCGTCTTAA";
        CharString database = "TACGATCCACACCGCGTCT";
        Score<int, Simple> scoringScheme(3, -2, -1, -3);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(1, 2, 6, 8));
        appendValue(seedChain, TSeed(12, 10, 19, 19));
        
        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), query);
        assignSource(row(alignment, 1), database);
        
        int result = bandedChainAlignment(alignment, seedChain, 2, scoringScheme, AlignConfig<false, false, false, false>());
        SEQAN_ASSERT_EQ(result, 24);

        // Compare alignment rows.
        //
        // Note that leading and trailing gaps are not stored, so the
        // start and end gaps do not appear in the strings.
        SEQAN_ASSERT(row(alignment, 0) == "ACG-TCCTCGTACACCG--TCTTAA");
        SEQAN_ASSERT(row(alignment, 1) == "TACGATCC----ACACCGCGTCT");
    }
    // Test on infixes.
    {
        CharString query = "TACGTCCTCGTACACCGTCTTAA";
        CharString database = "TTTACGATCCACACCGCGTCT";
        Score<int, Simple> scoringScheme(3, -2, -1, -3);
        
        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(1, 2, 6, 8));
        appendValue(seedChain, TSeed(12, 10, 19, 19));
        
        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), query);
        setClippedBeginPosition(row(alignment, 0), 1);
        setClippedEndPosition(row(alignment, 0), length(query));
        assignSource(row(alignment, 1), database);
        setClippedBeginPosition(row(alignment, 1), 2);
        setClippedEndPosition(row(alignment, 1), length(database));
        
        int result = _bandedChainAlignment(alignment, seedChain, 2, scoringScheme, AlignConfig<false, false, false, false>(), Gotoh());
        SEQAN_ASSERT_EQ(result, 24);

        // Compare alignment rows.
        //
        // Note that leading and trailing gaps are not stored, so the
        // start and end gaps do not appear in the strings.
        SEQAN_ASSERT(row(alignment, 0) == "ACG-TCCTCGTACACCG--TCTTAA");
        SEQAN_ASSERT(row(alignment, 1) == "TACGATCC----ACACCGCGTCT");
    }
    // Test on whole strings -- linear gap costs.
    {
        // Resulting alignment should be something like this (seeds
        // are marked by < and >).
        //
        //     > <    > <    >  <
        //   GGCGATNNNCAT--GGCACA
        //   --CGA-ATCCATCCCACACA
        CharString sequence0 = "GGCGATNNNCATGGCACA";
        CharString sequence1 = "CGAATCCATCCCACACA";  // TODO(holtgrew): Switch back again.
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(2, 0, 6, 5));
        appendValue(seedChain, TSeed(9, 6, 12, 9));
        appendValue(seedChain, TSeed(14, 11, 16, 17));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        int result = _bandedChainAlignment(alignment, seedChain, 2, scoringScheme, AlignConfig<false, false, false, false>(), Gotoh());
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        //
        // Note that leading and trailing gaps are not stored, so the
        // start and end gaps do not appear in the strings.
        SEQAN_ASSERT(row(alignment, 0) == "GGCGA-TNNNCATGGCACA"/*--*/);
        SEQAN_ASSERT(row(alignment, 1) == /*--*/"CGAATC--CATCCCACACA");
    }
    // Test on infixes -- linear gap costs.
    {
        // The test data is the same as above but with a T and a TT prepended.
        CharString sequence0 = "TGGCGATNNNCATGGCACA";
        CharString sequence1 = "TTCGAATCCATCCCACACA";
        Score<int, Simple> scoringScheme(2, -1, -2);

        String<TSeed > seedChain;
        appendValue(seedChain, TSeed(2, 0, 6, 5));
        appendValue(seedChain, TSeed(9, 6, 12, 9));
        appendValue(seedChain, TSeed(14, 11, 16, 17));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        // TODO(holtgrew): This infix assignment for alignments is a bit creepy, maybe one of the too many shortcuts?
        assignSource(row(alignment, 0), sequence0);
        setClippedBeginPosition(row(alignment, 0), 1);
        setClippedEndPosition(row(alignment, 0), length(sequence0));
        assignSource(row(alignment, 1), sequence1);
        setClippedBeginPosition(row(alignment, 1), 2);
        setClippedEndPosition(row(alignment, 1), length(sequence1));

        int result = _bandedChainAlignment(alignment, seedChain, 2, scoringScheme, AlignConfig<false, false, false, false>(), Gotoh());
        SEQAN_ASSERT_EQ(result, 5);

        // Compare alignment rows.
        //
        // Note that leading and trailing gaps are not stored, so the
        // start and end gaps do not appear in the strings.
        SEQAN_ASSERT(row(alignment, 0) == "GGCGA-TNNNCATGGCACA"/*--*/);
        SEQAN_ASSERT(row(alignment, 1) == /*--*/"CGAATC--CATCCCACACA");
    }
}

#endif  // TEST_SEEDS_TEST_ALIGN_CHAIN_BANDED_H_
