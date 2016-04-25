// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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

#ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_H_
#define TESTS_ALIGN_TEST_ALIGN_SIMD_H_

#include <seqan/basic.h>
#include <seqan/align.h>

// These are interface tests.

// How can we reduce redundant test instances?

// Global, Overlap, Semi-Global
// Local, affine
// Align, Graph, Gaps, Fragments, Score

template <typename TAlignTraits>

// Problem is clearly that the result is different.
SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_global_linear)
{
    using namespace seqan;

    struct AlignTrait
    {
        
    }
    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));

        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 6);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
            SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
        }
    }

    // Alignment with both leading gaps - Simd version
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 5);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
            SEQAN_ASSERT_EQ(ssV.str(), "--AAA---G--TT");
        }
    }

    // Alignment with both leading and trailing gaps in different rows - Simd version
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }
        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 5);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTT-----GGG");
            SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
        }
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_global_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet stringsSimdH, stringsSimdV;

        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 6);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet stringsSimdH, stringsSimdV;

        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 5);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet stringsSimdH, stringsSimdV;

        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 5);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_global_affine)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Alignment with gaps in horizontal sequence - Simd version
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 2);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
            SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
        }
    }

    // Alignment with gaps vertical sequence - Simd version
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 1);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
            SEQAN_ASSERT_EQ(ssV.str(), "AAA-----GTT--");
        }
    }

    // Alignment with gaps in different rows - Simd version
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 1);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
            SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
        }
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_global_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;
    AlignConfig<> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 2);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 1);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 1);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_overlap_linear)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 6);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
            SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
        }
    }

    // Alignment with both leading and trailing gaps in one row - Simd version
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 9);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
            SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
        }
    }

    // Alignment with both leading and trailing gaps in different rows - Simd version
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 13);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
            SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
        }
    }

    // Alignment with no overlap - Simd version
    {
        Dna5String strH = "GGGGGGGGG";
        Dna5String strV = "CCCCCCCCC";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 0);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "---------GGGGGGGGG");
            SEQAN_ASSERT_EQ(ssV.str(), "CCCCCCCCC---------");
        }
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_overlap_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 6);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 9);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }
        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 13);
    }

    // Alignment with no overlap.
    {
        DnaString strH = "CCCCCCCCC";
        Dna5String strV = "GGTG";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 0);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_overlap_affine)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 4);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "----ATGT");
            SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT--");
        }
    }

    // Alignment with both leading and trailing gaps in one row - Simd version
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 7);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
            SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
        }
    }

    // Alignment with both leading and trailing gaps in different rows - Simd version
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 13);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
            SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
        }
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_overlap_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;

    AlignConfig<true, true, true, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet stringsSimdH, stringsSimdV;

        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 4);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 7);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 13);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_semi_global_linear)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 6);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
            SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
        }
    }

    // Alignment with both leading and trailing gaps in one row - Simd version
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 9);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
            SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
        }
    }

    // Alignment with both leading and trailing gaps in different rows - Simd version
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 8);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTT-----GGG");
            SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
        }
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_semi_global_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;
    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 6);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 9);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 8);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_semi_global_affine)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 2);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
            SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
        }
    }

    // Alignment with both leading and trailing gaps in one row - Simd version
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 7);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
            SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
        }
    }

    // Alignment with both leading and trailing gaps in different rows - Simd version
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        StringSet<Align<Dna5String> > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<Dna5String> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), strH);
            assignSource(row(align, 1), strV);
            appendValue(alignments, align);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignment(alignments, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], 6);
            std::stringstream ssH, ssV;
            ssH << row(alignments[i], 0);
            ssV << row(alignments[i], 1);
            SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
            SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
        }
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_semi_global_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 2);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 7);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet stringsSimdH, stringsSimdV;
        for(unsigned i = 0; i < 34; ++i)
        {
            appendValue(stringsSimdH, strH);
            appendValue(stringsSimdV, strV);
        }

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        String<int> scores = globalAlignmentScore(stringsSimdH, stringsSimdV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));
        for(size_t i = 0; i < 34; ++i)
            SEQAN_ASSERT_EQ(scores[i], 6);
    }
}

#endif  // #ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_H_
