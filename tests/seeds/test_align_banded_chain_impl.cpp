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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Tests the different locations of seeds within the seed chain and the
// special conditions at the beginning and the end of the global grid.
// ==========================================================================

#include <sstream>

#include <seqan/basic.h>
#include <seqan/stream.h>  // for printing seqan::String<>

#include <seqan/seeds.h>

template <typename TGapCosts>
void testBandedChainAlignmentEmptyChain(TGapCosts const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<Simple>, Unordered> TSeedSet;
    AlignConfig<true, true, true, true> alignConfig;

    CharString seqH = "ACGATCGATCGACTGACT";
    CharString seqV = "ACGATCGATCGACTGACT";

    TSeedSet seedSet;

    Align<CharString> align;
    resize(rows(align), 2);
    assignSource(row(align,0), seqH);
    assignSource(row(align,1), seqV);

    Score<int,Simple> scoreScheme(5, -3, -5);

    if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        scoreScheme = Score<int,Simple>(5, -3, -1, -5);

    int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);
    SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());
    SEQAN_ASSERT_EQ(row(align,0), seqH);
    SEQAN_ASSERT_EQ(row(align,1), seqV);
}

template <typename TGapCosts>
void testBandedChainAlignmentOneSeed(TGapCosts const &)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    typedef SeedSet<TSeed, Unordered> TSeedSet;
    AlignConfig<true, true, true, true> alignConfig;

    // Tested overlapping anchor without shift which overlaps in all directions.
    {
                        //012345678901
        DnaString seqH = "AAAAAAAAAAA";
        DnaString seqV = "AAAAAAAAAAA";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(0,0,11,11), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);
        SEQAN_ASSERT_EQ(score, 55);
        SEQAN_ASSERT_EQ(row(align,0), seqH);
        SEQAN_ASSERT_EQ(row(align,1), seqV);
    }

    // Tested overlapping anchor without shift which overlaps in all directions.
    {
                       // 0123456789012345678901
        DnaString seqH = "CCCCCAAAAAAAAAAATTTTT";
        DnaString seqV = "GGGGGAAAAAAAAAAAGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(5,5,16,16), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);

        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 37);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCC-----AAAAAAAAAAA-----TTTTT");
            SEQAN_ASSERT_EQ(row(align,1), "-----GGGGGAAAAAAAAAAAGGGGG-----");
        }
        if (IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 25);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCCAAAAAAAAAAATTTTT");
            SEQAN_ASSERT_EQ(row(align,1), "GGGGGAAAAAAAAAAAGGGGG");
        }
    }

    // Tested overlapping anchor without shift which overlaps with beginning of horizontal sequence and end of vertical sequence.
    {                  // 0         1         2         3
                       // 01234567890123456789012345678901
        DnaString seqH = "CCCCCGTGTGAAAAAAAAAAATCGATTTTTT";
        DnaString seqV = "GGGGGCGTGCAAAAAAAAAAAGTAGCGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(10,10,21,21), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);

        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 56);
            SEQAN_ASSERT_EQ(row(align,0), "-----CCCCCGTGTG-AAAAAAAAAAA-T--CG----ATTTTTT");
            SEQAN_ASSERT_EQ(row(align,1), "GGGGGC------GTGCAAAAAAAAAAAGTAGCGGGGG-------");
        }
        if (IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 43);
            SEQAN_ASSERT_EQ(row(align,0), "-CCCCCGTGTGAAAAAAAAAAA-T--CGATTTTTT");
            SEQAN_ASSERT_EQ(row(align,1), "GGGGGCGTG-CAAAAAAAAAAAGTAGCGGGGG---");
        }
    }

    // Tested overlapping anchor without shift which overlaps with beginning of vertical sequence and end of horizontal sequence.
    {                  // 0         1         2         3
                       // 01234567890123456789012345678901
        DnaString seqH = "CCCCCAAAAAAAAAAATCGATTTTTT";
        DnaString seqV = "GGGGGCGTGCAAAAAAAAAAAGTAGCGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(5,10,16,21), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);

        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 52);
            SEQAN_ASSERT_EQ(row(align,0), "-----CCCCCAAAAAAAAAAA-T--CG----ATTTTTT");
            SEQAN_ASSERT_EQ(row(align,1), "GGGGGCGTGCAAAAAAAAAAAGTAGCGGGGG-------");
        }
        if (IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 44);
            SEQAN_ASSERT_EQ(row(align,0), "-----CCCCCAAAAAAAAAAA-T--CGATTTTTT");
            SEQAN_ASSERT_EQ(row(align,1), "GGGGGCGTGCAAAAAAAAAAAGTAGCGGGGG---");
        }
    }

    // Tested none overlapping anchor without shift.
    {                  // 0         1         2         3
                       // 01234567890123456789012345678901
        DnaString seqH = "CCCCCGTGTGAAAAAAAAAAATCGAT";
        DnaString seqV = "CGTGCAAAAAAAAAAAGTAGCGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(10,5,21,16), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);

        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 66);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCCGTGTGAAAAAAAAAAATCGAT--------");
            SEQAN_ASSERT_EQ(row(align,1), "----CGTG-CAAAAAAAAAAA--G-TAGCGGGGG");
        }
        if (IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 63);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCCGTGTGAAAAAAAAAAA-TCGAT----");
            SEQAN_ASSERT_EQ(row(align,1), "----CGTG-CAAAAAAAAAAAGTAGCGGGGG");
        }
    }

    // Tested none overlapping anchor with band shift in vertical direction.
    {                  // 0         1         2         3
                       // 01234567890123456789012345678901
        DnaString seqH = "CCCCCGTGTGAAAAAAAAAAATCGAT";
        DnaString seqV = "CGTGCAAAAAAAAAAAGTAGCGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(12,5,19,16), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);

        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 66);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCCGTGTGAAAAAAAAAAATCGAT--------");
            SEQAN_ASSERT_EQ(row(align,1), "----CGTG-CAAAAAAAAAAA--G-TAGCGGGGG");
        }
        if (IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 63);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCCGTGTGAAAAAAAAAAA-TCGAT----");
            SEQAN_ASSERT_EQ(row(align,1), "----CGTG-CAAAAAAAAAAAGTAGCGGGGG");
        }
    }

    // Tested none overlapping anchor with band shift in horizontal direction.
    {                  // 0         1         2         3
                       // 01234567890123456789012345678901
        DnaString seqH = "CCCCCGTGTGAAAAAAAAAAATCGAT";
        DnaString seqV = "CGTGCAAAAAAAAAAAGTAGCGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(5,12,16,14), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);

        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 64);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCCGTGTGAAAAAAAAAA--ATCGAT----");
            SEQAN_ASSERT_EQ(row(align,1), "----CGTGCAAAAAAAAAAAGTAGCG--GGGG");
        }
        if (IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 62);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCCGTGTGAAAAAAAAAAATCGAT----");
            SEQAN_ASSERT_EQ(row(align,1), "----CGTGCAAAAAAAAAAAGTAGCGGGGG");
        }
    }

    // Tested one shifted anchor which exceeds the beginning of the horizontal sequence and the end of the vertical sequence.
    {                  // 0         1         2         3
                       // 01234567890123456789012345678901
        DnaString seqH = "CCCCCGTGTGAAAAAAAAAAATCGAT";
        DnaString seqV = "CGTGCAAAAAAAAAAAGTAGCGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(12,5,19,16), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 8);

        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 66);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCCGTGTGAAAAAAAAAAATCGAT--------");
            SEQAN_ASSERT_EQ(row(align,1), "----CGTG-CAAAAAAAAAAA--G-TAGCGGGGG");
        }
        if (IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 63);
            SEQAN_ASSERT_EQ(row(align,0), "CCCCCGTGTGAAAAAAAAAAA-TCGAT----");
            SEQAN_ASSERT_EQ(row(align,1), "----CGTG-CAAAAAAAAAAAGTAGCGGGGG");
        }
    }

    // Tested shifted anchor which exceeds the beginning of the vertical sequence and the end of the horizontal sequence.
    {                  // 0         1         2         3
                       // 01234567890123456789012345678901
        DnaString seqH = "CGTGCAAAAAAAAAAAGTAGCGGGGG";
        DnaString seqV = "CCCCCGTGTGAAAAAAAAAAATCGAT";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(5,12,16,19), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 8);

        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 66);
            SEQAN_ASSERT_EQ(row(align,0), "----CGTG-CAAAAAAAAAAA--G-TAGCGGGGG");
            SEQAN_ASSERT_EQ(row(align,1), "CCCCCGTGTGAAAAAAAAAAATCGAT--------");
        }
        if (IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 63);
            SEQAN_ASSERT_EQ(row(align,0), "----CGTG-CAAAAAAAAAAAGTAGCGGGGG");
            SEQAN_ASSERT_EQ(row(align,1), "CCCCCGTGTGAAAAAAAAAAA-TCGAT----");
        }
    }
}

template <typename TGapCosts>
void testBandedChainAlignmentTwoSeeds(TGapCosts const &)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    typedef SeedSet<TSeed, Unordered> TSeedSet;
    AlignConfig<true, true, true, true> alignConfig;
    
    // Test two anchors that does not overlap.
    {
                        //012345678901234567890123456789
        DnaString seqH = "AAAAAAAAAAATTTTTTTTTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCC";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(0,2,11,11), Single());
        addSeed(seedSet, TSeed(21,21,26,29), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 65);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAATTTTTTTTTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCC");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 67);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAATTTTTTTTTT----------CCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAA----------GGGGGGGGGGCCCCCCCC");
        }
    }

    // Test two anchors one position apart in vertical direction.
    {
                        //012345678901234567890123456789
        DnaString seqH = "AAAAAAAAAAATTTTTTTTTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAACCCCCCCCCC";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(2,0,10,11), Single());
        addSeed(seedSet, TSeed(21,11,28,20), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 49);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAATTTTTTTTTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAA--------CCCCCCCCCC");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 81);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAATTTTTTTTTTCCCCCCCC--");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAA----------CCCCCCCCCC");
        }
    }

    // Test two anchors one position apart in horizontal direction.
    {
                        //012345678901234567890123456789
        DnaString seqH = "AAAAAAAAAAACCCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAATTTTTTTTTTCCCCCCCC";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(0,2,11,10), Single());
        addSeed(seedSet, TSeed(11,21,20,28), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 49);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAA--------CCCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAATTTTTTTTTTCCCCCCCC");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 81);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAA----------CCCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAATTTTTTTTTTCCCCCCCC--");
        }
    }

    // Test two anchors one position apart in horizontal & vertical direction.
    {
                        //012345678901234567890123456789
        DnaString seqH = "AAAAAAAAAAATTTTTTTTTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCC";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(2,6,17,15), Single());
        addSeed(seedSet, TSeed(17,15,29,29), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 5);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 65);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAATTTTTTTTTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCC");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 55);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAATTTTTTTTTTCC----------CCCCCC--");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAA------------GGGGGGGGGGCCCCCCCC");
        }
    }

    // Test two anchors one position apart in horizontal & vertical direction.
    {
                        //012345678901234567890123456789
        DnaString seqH = "AAAAAAAAAAATTTTTTTTTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCC";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(2,6,17,15), Single());
        addSeed(seedSet, TSeed(17,15,29,29), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 7);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 65);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAATTTTTTTTTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCC");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 67);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAATTTTTTTTTT----------CCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAA----------GGGGGGGGGGCCCCCCCC");
        }
    }

    // Test two anchors one position apart in horizontal & vertical direction.
    {
                        //012345678901234567890123456789
        DnaString seqH = "AAAAAAAAAAATTTTTTTTTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCC";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        addSeed(seedSet, TSeed(2,6,17,15), Single());
        addSeed(seedSet, TSeed(17,15,29,29), Single());

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 1);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 65);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAATTTTTTTTTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCC");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 61);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAA----TTTTTTTT----TTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGG--------GGGGGGCCCCCCCC");
        }
    }
}

template <typename TGapCosts>
void testBandedChainAlignmentThreeSeeds(TGapCosts const &)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    typedef StringSet<TSeed> TSeedSet;
    AlignConfig<true, true, true, true> alignConfig;

    // Test two anchors one position apart in horizontal & vertical direction.
    {                   //0         1         2         3         4         5         6
                        //0123456789012345678901234567890123456789012345678901234567890
        DnaString seqH = "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        appendValue(seedSet, TSeed(10,11,21,21));
        appendValue(seedSet, TSeed(21,29,32,40));
        appendValue(seedSet, TSeed(38,40,40,48));

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 1);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 5);
            SEQAN_ASSERT_EQ(row(align,0), "-TTTTTTTTTTGGGGGGGGGG-------GAAAAAAAAAAAGGGGGGTTCCCCCCCC----");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 65);
            SEQAN_ASSERT_EQ(row(align,0), "-----------TTTTTTTTTTGGGGGGGGGG-------GAAAAAAAAAAAGGGGGGTTCCCCCCCC------------------");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAA----------GGGGGGGGGGCCCCCCCCAAAAAAAAAAA------TT--------TTTTTTGGGGGGGGGGGG");
        }

        clearGaps(row(align, 0));
        clearGaps(row(align, 1));
        score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 15);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 109);
            SEQAN_ASSERT_EQ(row(align,0), "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGG--GGTTCCCCCCCC-------------------------------");
            SEQAN_ASSERT_EQ(row(align,1), "---------------------AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 113);
            SEQAN_ASSERT_EQ(row(align,0), "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGG--GGTTCCCCCCCC-------------------------------");
            SEQAN_ASSERT_EQ(row(align,1), "---------------------AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
        }

        clearGaps(row(align, 0));
        clearGaps(row(align, 1));
        score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 6);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 3);
            SEQAN_ASSERT_EQ(row(align,0), "-TTTTTTTTTTGGGGGGGGGG-------GAAAAAAAAAAA-------GGGGGGTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG---");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 85);
            SEQAN_ASSERT_EQ(row(align,0), "-----------TTTTTTTTTTGGGGGGGGGG-------GAAAAAAAAAAA--------GGGGGG------TTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAA----------GGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG----------");
        }
    }


    // Test two anchors one position apart in horizontal & vertical direction.
    {                   //0         1         2         3         4         5         6
                        //0123456789012345678901234567890123456789012345678901234567890
        DnaString seqH = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG";
        DnaString seqV = "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        appendValue(seedSet, TSeed(11,10,21,21));
        appendValue(seedSet, TSeed(29,21,40,32));
        appendValue(seedSet, TSeed(40,38,48,40));

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 1);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 5);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
            SEQAN_ASSERT_EQ(row(align,1), "-TTTTTTTTTTGGGGGGGGGG-------GAAAAAAAAAAAGGGGGGTTCCCCCCCC----");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 65);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAA----------GGGGGGGGGGCCCCCCCCAAAAAAAAAAA------TT--------TTTTTTGGGGGGGGGGGG");
            SEQAN_ASSERT_EQ(row(align,1), "-----------TTTTTTTTTTGGGGGGGGGG-------GAAAAAAAAAAAGGGGGGTTCCCCCCCC------------------");
        }

        clearGaps(row(align, 0));
        clearGaps(row(align, 1));
        score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 15);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 109);
            SEQAN_ASSERT_EQ(row(align,0), "---------------------AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
            SEQAN_ASSERT_EQ(row(align,1), "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGG--GGTTCCCCCCCC-------------------------------");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 113);
            SEQAN_ASSERT_EQ(row(align,0), "---------------------AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
            SEQAN_ASSERT_EQ(row(align,1), "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGG--GGTTCCCCCCCC-------------------------------");
        }

        clearGaps(row(align, 0));
        clearGaps(row(align, 1));
        score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 6);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 3);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG---");
            SEQAN_ASSERT_EQ(row(align,1), "-TTTTTTTTTTGGGGGGGGGG-------GAAAAAAAAAAA-------GGGGGGTTCCCCCCCC");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 85);
            SEQAN_ASSERT_EQ(row(align,0), "AAAAAAAAAAA----------GGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG----------");
            SEQAN_ASSERT_EQ(row(align,1), "-----------TTTTTTTTTTGGGGGGGGGG-------GAAAAAAAAAAA--------GGGGGG------TTCCCCCCCC");
        }
    }

    // Test two anchors one position apart in horizontal & vertical direction.
    {                   //0         1         2         3         4         5         6
                        //0123456789012345678901234567890123456789012345678901234567890
        DnaString seqH = "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        appendValue(seedSet, TSeed(10,11,21,21));
        appendValue(seedSet, TSeed(21,29,32,40));
        appendValue(seedSet, TSeed(38,40,40,48));

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 16);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 109);
            SEQAN_ASSERT_EQ(row(align,0), "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGG--GGGTTCCCCCCCC-------------------------------");
            SEQAN_ASSERT_EQ(row(align,1), "---------------------AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 113);
            SEQAN_ASSERT_EQ(row(align,0), "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGG--GGGTTCCCCCCCC-------------------------------");
            SEQAN_ASSERT_EQ(row(align,1), "---------------------AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
        }
    }
}

template <typename TGapCosts>
void testBandedChainAlignmentRareCases(TGapCosts const &)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    typedef StringSet<TSeed> TSeedSet;
    AlignConfig<false, true, false, true> alignConfig;

    // Test one anchor of minimal size with huge band extension.
    {                   //0         1         2         3         4         5         6
                        //0123456789012345678901234567890123456789012345678901234567890
        DnaString seqH = "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        appendValue(seedSet, TSeed(25,25,26,26));

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 50);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 89);
            SEQAN_ASSERT_EQ(row(align,0), "---------------------------------------TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG---------------------------");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 90);
            SEQAN_ASSERT_EQ(row(align,0), "-----------TTTTTTTTTTGGGGGGGGGGG--------AAAAAAAAAAA--------------GGGGGGTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAA-----------GGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG----------");
        }
    }

    // Test three anchors where the first two are covered by the third.
    {                   //0         1         2         3         4         5         6
                        //0123456789012345678901234567890123456789012345678901234567890
        DnaString seqH = "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        appendValue(seedSet, TSeed(0,0,2,2));
        appendValue(seedSet, TSeed(10,11,21,21));
        appendValue(seedSet, TSeed(25,25,26,26));

        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 50);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 89);
            SEQAN_ASSERT_EQ(row(align,0), "---------------------------------------TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG---------------------------");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 90);
            SEQAN_ASSERT_EQ(row(align,0), "-----------TTTTTTTTTTGGGGGGGGGGG--------AAAAAAAAAAA--------------GGGGGGTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAA-----------GGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG----------");
        }
    }

    // Test five anchors where the first two are covered by the third anchor and the fith anchor is covered by the fourth anchor.
    {                   //0         1         2         3         4         5         6
                        //0123456789012345678901234567890123456789012345678901234567890
        DnaString seqH = "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        appendValue(seedSet, TSeed(0,0,2,2));
        appendValue(seedSet, TSeed(10,11,21,21));

        appendValue(seedSet, TSeed(22,22,23,23));

        appendValue(seedSet, TSeed(30,29,37,40));
        appendValue(seedSet, TSeed(38,40,40,48));


        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 22);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, 82);
            SEQAN_ASSERT_EQ(row(align,0), "----------------------------------------TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAAT--TTTTTTTGGGGGGGGGGGG--------------------------");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, 53);
            SEQAN_ASSERT_EQ(row(align,0), "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAA--GGGGGGTTCCCCCCCC-------------------------------");
            SEQAN_ASSERT_EQ(row(align,1), "---------------------AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
        }
    }

    // Test five tiny anchors with minimal band extension of 1.
    {                   //0         1         2         3         4         5         6
                        //0123456789012345678901234567890123456789012345678901234567890
        DnaString seqH = "TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAAAGGGGGGTTCCCCCCCC";
        DnaString seqV = "AAAAAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG";

        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align,0), seqH);
        assignSource(row(align,1), seqV);

        TSeedSet seedSet;

        appendValue(seedSet, TSeed(1,1,2,2));
        appendValue(seedSet, TSeed(8,7,9,8));
        appendValue(seedSet, TSeed(22,22,23,23));
        appendValue(seedSet, TSeed(30,30,31,31));
        appendValue(seedSet, TSeed(35,40,36,41));


        Score<int,Simple> scoreScheme(5, -3, -5);
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
            scoreScheme = Score<int,Simple>(5, -3, -1, -5);

        int score = bandedChainAlignment(align, seedSet, scoreScheme, alignConfig, 1);
        if(IsSameType<TGapCosts, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ(score, -73);
            SEQAN_ASSERT_EQ(row(align,0), "--TTTTTTTTTTGGGGGGGGGGGAAAAAAAAAA--AGGGGG-----GTT----CCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAA-AAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG");
        }
        else
        {
            SEQAN_ASSERT_EQ(score, -15);
            SEQAN_ASSERT_EQ(row(align,0), "--TTTTTTTT------TTGGGGGGGGGGGAAAAAAAAAA------AGGGG---------------GGTTCCCCCCCC");
            SEQAN_ASSERT_EQ(row(align,1), "AAA-------AAAAAAAAGGGGGGGGGGCCCCCCCCAAAAAAAAAAATTTTTTTTGGGGGGGGGGGG----------");
        }
    }
}

template <typename TGapSpec>
void testBandedChainAlignmentBandExtension(TGapSpec)
{
    using namespace seqan;

    typedef Seed<Simple> TSeed;
    // Test alignment with chain of length 1 with various bandwidths
    // to check that too large bandwidths do not do any harm.
    for (unsigned k = 1; k < 6; ++k)
    {
        CharString sequence0 = "NNNAAANNN";
        CharString sequence1 =  "NCAAACNN";

        Score<int, Simple> scoringScheme(2, -1, -2);
        if (IsSameType<TGapSpec, AffineGaps>::VALUE)
            scoringScheme = Score<int, Simple>(2, -1, -1, -3);

        String<TSeed> seedChain;
        appendValue(seedChain, TSeed(3,2,6,5));

        Align<CharString, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        int result = bandedChainAlignment(alignment, seedChain, scoringScheme, AlignConfig<false, false, false, false>(), k);
        if (IsSameType<TGapSpec, LinearGaps>::VALUE)
        {
            SEQAN_ASSERT_EQ_MSG(result, 8, "k == %u", k);
            SEQAN_ASSERT_EQ(row(alignment,0), "NNNAAANNN");
            SEQAN_ASSERT_EQ(row(alignment,1), "-NCAAACNN");
        }
        else
        {
            SEQAN_ASSERT_EQ_MSG(result, 7, "k == %u", k);
            SEQAN_ASSERT_EQ(row(alignment,0), "NNNAAANNN");
            SEQAN_ASSERT_EQ_MSG(row(alignment,1), "-NCAAACNN", "k == %u", k);
        }

    }
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_issue_1020)
{
    using namespace seqan;

    DnaString query = "ATCTCTCTCAACAAAACAACGAGGAGGAGTGAAAAGAGAGAGAT";
    DnaString ref   = "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT";

    typedef Seed<Simple> TSeed;
    String<TSeed> seedChain;
    appendValue(seedChain, TSeed( 0,  0, 14));
    appendValue(seedChain, TSeed(30, 31, 14));
    Score<int, Simple> scoringScheme(2, -1, -2);

    Align<DnaString, ArrayGaps> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), query);
    assignSource(row(align, 1), ref);

    int res = bandedChainAlignment(align, seedChain, scoringScheme, 14);
    SEQAN_ASSERT_EQ(res, 80);
    SEQAN_ASSERT_EQ(row(align, 0), "ATCTCTCTCAACAA-AACAAC-GAGGAGGAGTGAAAAGAGAGAGAT");
    SEQAN_ASSERT_EQ(row(align, 1), "ATCTCTCTCAACAACAACAACGGAGGAGGAG-GAAAAGAGAGAGAT");
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_empty_set_linear)
{
    testBandedChainAlignmentEmptyChain(seqan::LinearGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_empty_set_affine)
{
    testBandedChainAlignmentEmptyChain(seqan::AffineGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_one_seed_linear)
{
    testBandedChainAlignmentOneSeed(seqan::LinearGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_one_seed_affine)
{
    testBandedChainAlignmentOneSeed(seqan::AffineGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_two_seeds_linear)
{
    testBandedChainAlignmentTwoSeeds(seqan::LinearGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_two_seeds_affine)
{
    testBandedChainAlignmentTwoSeeds(seqan::AffineGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_three_seeds_linear)
{
    testBandedChainAlignmentThreeSeeds(seqan::LinearGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_three_seeds_affine)
{
    testBandedChainAlignmentThreeSeeds(seqan::AffineGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_special_seeds_linear)
{
    testBandedChainAlignmentRareCases(seqan::LinearGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_special_seeds_affine)
{
    testBandedChainAlignmentRareCases(seqan::AffineGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_band_extensions_linear)
{
    testBandedChainAlignmentBandExtension(seqan::LinearGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_alignment_band_extensions_affine)
{
    testBandedChainAlignmentBandExtension(seqan::AffineGaps());
}

SEQAN_DEFINE_TEST(test_banded_chain_score_overflow_detection)
{
    using namespace seqan;
    Score<int8_t> score(2, -5, -2);
    SEQAN_ASSERT(!_checkScoreOverflow(127, score));
    SEQAN_ASSERT(_checkScoreOverflow(50, score));
    SEQAN_ASSERT(_checkScoreOverflow(24, score));
    SEQAN_ASSERT(_checkScoreOverflow(15, score));

    Score<int32_t> score2(2, -5, -2);
    SEQAN_ASSERT(_checkScoreOverflow(655536, score2));
    SEQAN_ASSERT(_checkScoreOverflow(127, score2));
    SEQAN_ASSERT(_checkScoreOverflow(50, score2));
    SEQAN_ASSERT(_checkScoreOverflow(24, score2));
    SEQAN_ASSERT(_checkScoreOverflow(15, score2));
}

SEQAN_BEGIN_TESTSUITE(test_banded_chain_impl)
{
    SEQAN_CALL_TEST(test_banded_chain_alignment_empty_set_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_empty_set_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_one_seed_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_one_seed_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_two_seeds_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_two_seeds_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_three_seeds_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_three_seeds_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_special_seeds_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_special_seeds_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_band_extensions_linear);
    SEQAN_CALL_TEST(test_banded_chain_alignment_band_extensions_affine);
    SEQAN_CALL_TEST(test_banded_chain_alignment_issue_1020);
    SEQAN_CALL_TEST(test_banded_chain_score_overflow_detection);
}
SEQAN_END_TESTSUITE
