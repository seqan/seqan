// ==========================================================================
//                     test_alignment_algorithms_global.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_GLOBAL_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_GLOBAL_H_

#include <sstream>

#include <seqan/basic.h>
#include <seqan/align.h>

// ==========================================================================
// Long Interfaces
// ==========================================================================

// ----------------------------------------------------------------------------
// Global alignments.
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_global_linear)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Alignment with both leading gaps.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 5);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---G--TT");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 5);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTT-----GGG");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_gaps_global_linear)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());
        SEQAN_ASSERT_EQ(score, 5);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---G--TT");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 5);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTT-----GGG");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_gaps_global_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .  \n"
                   << "        AT-G-T\n"
                   << "        || | |\n"
                   << "        ATAGAT\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
        // std::cerr << align << "\n";
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 5);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    \n"
                   << "        AAAAAGGGGTTTT\n"
                   << "          |||   |  ||\n"
                   << "        --AAA---G--TT\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 5);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .     \n"
                   << "        AAAAAATTTTT-----GGG\n"
                   << "              |||||     |||\n"
                   << "        ---TTTTTTTTGGGGGGGG\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_gaps_global_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 5, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 3, 1));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 0, 2));
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 5);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 11, 1, 4, 2));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 8, 1, 3, 1));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 2, 1, 0, 3));
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 5);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 11, 1, 13, 3));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 0, 8));
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

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 6);
        SEQAN_ASSERT_EQ(score2, 6);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 5);
        SEQAN_ASSERT_EQ(score2, 5);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 5);
        SEQAN_ASSERT_EQ(score2, 5);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_global_affine)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Alignment with gaps in horizontal sequence.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Alignment with gaps vertical sequence.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Align<Dna5String> align2;
        resize(rows(align2), 2);
        assignSource(row(align2, 0), strH);
        assignSource(row(align2, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);

        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());
        SEQAN_ASSERT_EQ(score, 1);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-----GTT--");
    }

    // Alignment with gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 1);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_gaps_global_affine)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 1);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-----GTT--");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);

        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 1);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_gaps_global_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .  \n"
                   << "        AT-G-T\n"
                   << "        || | |\n"
                   << "        ATAGAT\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 1);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    \n"
                   << "        AAAAAGGGGTTTT\n"
                   << "        |||     |||  \n"
                   << "        AAA-----GTT--\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 1);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .     \n"
                   << "        AAAAAATTTTTGGG-----\n"
                   << "              ||||||||     \n"
                   << "        ---TTTTTTTTGGGGGGGG\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_gaps_global_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;
    AlignConfig<> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 2);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 5, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 3, 1));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 0, 2));
        // std::cerr << align << "\n";
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 1);
        SEQAN_ASSERT_EQ(length(fragments), 2u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 8, 1, 3, 3));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 0, 1, 0, 3));
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 1);
        SEQAN_ASSERT_EQ(length(fragments), 1u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 0, 11));
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

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score1, 2);
        SEQAN_ASSERT_EQ(score2, 2);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score1, 1);
        SEQAN_ASSERT_EQ(score2, 1);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score1, 1);
        SEQAN_ASSERT_EQ(score2, 1);
    }
}

// ----------------------------------------------------------------------------
// Overlap alignments.
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_overlap_linear)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 9);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 13);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }

    // Alignment with no overlap.
    {
        Dna5String strH = "GGGGGGGGG";
        Dna5String strV = "CCCCCCCCC";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 0);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "---------GGGGGGGGG");
        SEQAN_ASSERT_EQ(ssV.str(), "CCCCCCCCC---------");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_gaps_overlap_linear)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";
        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 9);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 13);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }

    // Alignment with no overlap.
    {
        DnaString strH = "GGGGGGGGGGGGG";
        Dna5String strV = "CCCCCCCCCC";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 0);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "----------GGGGGGGGGGGGG");
        SEQAN_ASSERT_EQ(ssV.str(), "CCCCCCCCCC-------------");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_gaps_overlap_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
    AlignConfig<true, true, true, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .  \n"
                   << "        AT-G-T\n"
                   << "        || | |\n"
                   << "        ATAGAT\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 9);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    \n"
                   << "        AAAAAGGGGTTTT\n"
                   << "          |||   |||  \n"
                   << "        --AAA---GTT--\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 13);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .     \n"
                   << "        AAAAAATTTTTGGG-----\n"
                   << "              ||||||||     \n"
                   << "        ---TTTTTTTTGGGGGGGG\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with no ovrlap.
    {
        Dna5String strH = "CCCCCCCCCC";
        Dna5String strV = "AAATTATTATATTATA";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 0);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :    .  \n"
                   << "        ----------------CCCCCCCCCC\n"
                   << "                                  \n"
                   << "        AAATTATTATATTATA----------\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_gaps_overlap_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 5, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 3, 1));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 0, 2));
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 9);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 8, 1, 3, 3));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 0, 3));
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 13);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 0, 11));
    }

    // Alignment with no overlap.
    {
        DnaString strH = "CCCCCCCCC";
        Dna5String strV = "GGTG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 0);

        SEQAN_ASSERT_EQ(empty(fragments), true);
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

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 6);
        SEQAN_ASSERT_EQ(score2, 6);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 9);
        SEQAN_ASSERT_EQ(score2, 9);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 13);
        SEQAN_ASSERT_EQ(score2, 13);
    }

    // Alignment with no overlap.
    {
        DnaString strH = "CCCCCCCCC";
        Dna5String strV = "GGTG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 0);
        SEQAN_ASSERT_EQ(score2, 0);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_overlap_affine)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 4);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "----ATGT");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT--");
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 7);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 13);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_gaps_overlap_affine)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 4);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "----ATGT");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT--");
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 7);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 13);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_gaps_overlap_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 4);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    \n"
                   << "        ----ATGT\n"
                   << "            ||  \n"
                   << "        ATAGAT--\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 7);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    \n"
                   << "        AAAAAGGGGTTTT\n"
                   << "          |||   |||  \n"
                   << "        --AAA---GTT--\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 13);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .     \n"
                   << "        AAAAAATTTTTGGG-----\n"
                   << "              ||||||||     \n"
                   << "        ---TTTTTTTTGGGGGGGG\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_gaps_overlap_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;
    AlignConfig<true, true, true, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 4);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 0, 1, 4, 2));
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 7);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 8, 1, 3, 3));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 0, 3));
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 13);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 0, 11));
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

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score1, 4);
        SEQAN_ASSERT_EQ(score2, 4);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score1, 7);
        SEQAN_ASSERT_EQ(score2, 7);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score1, 13);
        SEQAN_ASSERT_EQ(score2, 13);
    }
}

// ----------------------------------------------------------------------------
// Semi-Global alignments.
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_semi_global_linear)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
        // std::cerr << align << "\n";
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 9);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 8);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTT-----GGG");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_gaps_semi_global_linear)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 9);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 8);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTT-----GGG");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_gaps_semi_global_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .  \n"
                   << "        AT-G-T\n"
                   << "        || | |\n"
                   << "        ATAGAT\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 9);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    \n"
                   << "        AAAAAGGGGTTTT\n"
                   << "          |||   |||  \n"
                   << "        --AAA---GTT--\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 8);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .     \n"
                   << "        AAAAAATTTTT-----GGG\n"
                   << "              |||||     |||\n"
                   << "        ---TTTTTTTTGGGGGGGG\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_gaps_semi_global_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 5, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 3, 1));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 0, 2));
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 9);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 8, 1, 3, 3));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 0, 3));
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 8);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 11, 1, 13, 3));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 3, 1, 0, 8));
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

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 6);
        SEQAN_ASSERT_EQ(score2, 6);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 9);
        SEQAN_ASSERT_EQ(score2, 9);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score1, 8);
        SEQAN_ASSERT_EQ(score2, 8);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_semi_global_affine)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
        // std::cerr << align << "\n";
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 7);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_gaps_semi_global_affine)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 7);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "--AAA---GTT--");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_gaps_semi_global_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .  \n"
                   << "        AT-G-T\n"
                   << "        || | |\n"
                   << "        ATAGAT\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 7);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    \n"
                   << "        AAAAAGGGGTTTT\n"
                   << "          |||   |||  \n"
                   << "        --AAA---GTT--\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strSet;
        appendValue(strSet, strH);
        appendValue(strSet, strV);

        TAlignmentGraph alignmentGraph(strSet);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .     \n"
                   << "        AAAAAATTTTTGGG-----\n"
                   << "              ||||||||     \n"
                   << "        ---TTTTTTTTGGGGGGGG\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_gaps_semi_global_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps.
    {
        DnaString strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 2);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 5, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 3, 1));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 0, 2));
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 7);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 8, 1, 3, 3));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 0, 3));
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 6);
        SEQAN_ASSERT_EQ(length(fragments), 1u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 0, 11));
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

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score1, 2);
        SEQAN_ASSERT_EQ(score2, 2);
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        DnaString strH = "AAAAAGGGGTTTT";
        Dna5String strV = "AAAGTT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score1, 7);
        SEQAN_ASSERT_EQ(score2, 7);
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTGGG";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);
        int score1 = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
        int score2 = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score1, 6);
        SEQAN_ASSERT_EQ(score2, 6);
    }
}

// ==========================================================================
// Shorter Interfaces
// ==========================================================================

// We only test the resulting scores and compare with the result from the
// longer interfaces (those are already tested above).

SEQAN_DEFINE_TEST(test_align_global_alignment_shorter_interfaces_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Gaps<Dna5String, ArrayGaps> TGaps;
    typedef Align<Dna5String> TAlign;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;


    AlignConfig<> alignConfig;
    String<TTraceSegment> traceSegments;

    Dna5String strH = "ATGT";
    Dna5String strV = "ATAGAT";

    TStringSet strings;
    appendValue(strings, strH);
    appendValue(strings, strV);

    TAlign refAlign(strings);

    TGaps refGapsH(strH);
    TGaps refGapsV(strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    // Compute reference results from long interfaces.

    int refRes = globalAlignment(refAlign, scoringScheme, alignConfig, NeedlemanWunsch());

    int res1 = globalAlignment(refGapsH, refGapsV, scoringScheme, alignConfig, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res1, refRes);

    TFragmentString refFragments;
    int res2 = globalAlignment(refFragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res2, refRes);

    TAlignmentGraph refAlignmentGraph(strings);
    int res3 = globalAlignment(refAlignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res3, refRes);

    // Compute with the variant that does not have an AlignConfig parameter.
    {
        TAlign align(strings);
        res1 = globalAlignment(align, scoringScheme, NeedlemanWunsch());
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(align == refAlign);

        TGaps gapsH(strH), gapsV(strV);
        res1 = globalAlignment(gapsH, gapsV, scoringScheme, NeedlemanWunsch());
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(gapsH == refGapsH);
        SEQAN_ASSERT(gapsV == refGapsV);

        TFragmentString fragments;
        res1 = globalAlignment(fragments, strings, scoringScheme, NeedlemanWunsch());
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(fragments == refFragments);

        TAlignmentGraph alignmentGraph(strings);
        res1 = globalAlignment(alignmentGraph, scoringScheme, NeedlemanWunsch());
        SEQAN_ASSERT_EQ(res1, refRes);
        std::stringstream ss1, ss2;
        ss1 << alignmentGraph;
        ss2 << refAlignmentGraph;
        SEQAN_ASSERT_EQ(ss1.str(), ss2.str());
    }

    // Compute with the variant that does not have an alignment tag.
    {
        TAlign align(strings);
        res1 = globalAlignment(align, scoringScheme, alignConfig);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(align == refAlign);

        TGaps gapsH(strH), gapsV(strV);
        res1 = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(gapsH == refGapsH);
        SEQAN_ASSERT(gapsV == refGapsV);

        TFragmentString fragments;
        res1 = globalAlignment(fragments, strings, scoringScheme, alignConfig);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(fragments == refFragments);

        TAlignmentGraph alignmentGraph(strings);
        res1 = globalAlignment(alignmentGraph, scoringScheme, alignConfig);
        SEQAN_ASSERT_EQ(res1, refRes);
        std::stringstream ss1, ss2;
        ss1 << alignmentGraph;
        ss2 << refAlignmentGraph;
        SEQAN_ASSERT_EQ(ss1.str(), ss2.str());
    }

    // Compute with the variant that has no AlignConfig<> and algorithm tag.
    {
        TAlign align(strings);
        res1 = globalAlignment(align, scoringScheme);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(align == refAlign);

        TGaps gapsH(strH), gapsV(strV);
        res1 = globalAlignment(gapsH, gapsV, scoringScheme);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(gapsH == refGapsH);
        SEQAN_ASSERT(gapsV == refGapsV);

        TFragmentString fragments;
        res1 = globalAlignment(fragments, strings, scoringScheme);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(fragments == refFragments);

        TAlignmentGraph alignmentGraph(strings);
        res1 = globalAlignment(alignmentGraph, scoringScheme);
        SEQAN_ASSERT_EQ(res1, refRes);
        std::stringstream ss1, ss2;
        ss1 << alignmentGraph;
        ss2 << refAlignmentGraph;
        SEQAN_ASSERT_EQ(ss1.str(), ss2.str());
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_shorter_interfaces_affine)
{
    using namespace seqan;

    AlignConfig<> alignConfig;

    Dna5String strH = "ATGT";
    Dna5String strV = "ATAGAT";

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Gaps<Dna5String, ArrayGaps> TGaps;
    typedef Align<Dna5String> TAlign;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    TStringSet strings;
    appendValue(strings, strH);
    appendValue(strings, strV);

    TAlign refAlign(strings);

    TGaps refGapsH(strH);
    TGaps refGapsV(strV);

    Score<int, Simple> scoringScheme(2, -1, -1, -3);

    // Compute reference results from long interfaces.

    int refRes = globalAlignment(refAlign, scoringScheme, alignConfig, Gotoh());

    int res1 = globalAlignment(refGapsH, refGapsV, scoringScheme, alignConfig, Gotoh());
    SEQAN_ASSERT_EQ(res1, refRes);

    TFragmentString refFragments;
    int res2 = globalAlignment(refFragments, strings, scoringScheme, alignConfig, Gotoh());
    SEQAN_ASSERT_EQ(res2, refRes);

    TAlignmentGraph refAlignmentGraph(strings);
    int res3 = globalAlignment(refAlignmentGraph, scoringScheme, alignConfig, Gotoh());
    SEQAN_ASSERT_EQ(res3, refRes);

    // Compute with the variant that does not have an AlignConfig parameter.
    {
        TAlign align(strings);
        res1 = globalAlignment(align, scoringScheme, Gotoh());
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(align == refAlign);

        TGaps gapsH(strH), gapsV(strV);
        res1 = globalAlignment(gapsH, gapsV, scoringScheme, Gotoh());
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(gapsH == refGapsH);
        SEQAN_ASSERT(gapsV == refGapsV);

        TFragmentString fragments;
        res1 = globalAlignment(fragments, strings, scoringScheme, Gotoh());
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(fragments == refFragments);

        TAlignmentGraph alignmentGraph(strings);
        res1 = globalAlignment(alignmentGraph, scoringScheme, Gotoh());
        SEQAN_ASSERT_EQ(res1, refRes);
        std::stringstream ss1, ss2;
        ss1 << alignmentGraph;
        ss2 << refAlignmentGraph;
        SEQAN_ASSERT_EQ(ss1.str(), ss2.str());
    }

    // Compute with the variant that does not have an alignment tag.
    {
        TAlign align(strings);
        res1 = globalAlignment(align, scoringScheme, alignConfig);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(align == refAlign);

        TGaps gapsH(strH), gapsV(strV);
        res1 = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(gapsH == refGapsH);
        SEQAN_ASSERT(gapsV == refGapsV);

        TFragmentString fragments;
        res1 = globalAlignment(fragments, strings, scoringScheme, alignConfig);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(fragments == refFragments);

        TAlignmentGraph alignmentGraph(strings);
        res1 = globalAlignment(alignmentGraph, scoringScheme, alignConfig);
        SEQAN_ASSERT_EQ(res1, refRes);
        std::stringstream ss1, ss2;
        ss1 << alignmentGraph;
        ss2 << refAlignmentGraph;
        SEQAN_ASSERT_EQ(ss1.str(), ss2.str());
    }

    // Compute with the variant that has no AlignConfig<> and algorithm tag.
    {
        TAlign align(strings);
        res1 = globalAlignment(align, scoringScheme);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(align == refAlign);

        TGaps gapsH(strH), gapsV(strV);
        res1 = globalAlignment(gapsH, gapsV, scoringScheme);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(gapsH == refGapsH);
        SEQAN_ASSERT(gapsV == refGapsV);

        TFragmentString fragments;
        res1 = globalAlignment(fragments, strings, scoringScheme);
        SEQAN_ASSERT_EQ(res1, refRes);
        SEQAN_ASSERT(fragments == refFragments);

        TAlignmentGraph alignmentGraph(strings);
        res1 = globalAlignment(alignmentGraph, scoringScheme);
        SEQAN_ASSERT_EQ(res1, refRes);
        std::stringstream ss1, ss2;
        ss1 << alignmentGraph;
        ss2 << refAlignmentGraph;
        SEQAN_ASSERT_EQ(ss1.str(), ss2.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_global_different_container)
{
    using namespace seqan;

    // Global alignment with Segment and DnaString
    {
        Dna5String strH = "ATGTAT";
        Dna5String strV = "ATAGAT";

        Segment<Dna5String, PrefixSegment> prefixSegment(strH, 4);
        SEQAN_ASSERT_EQ(prefixSegment, "ATGT");

        Gaps<Segment<Dna5String, PrefixSegment>, ArrayGaps> gapsH(prefixSegment);
        Gaps<Dna5String> gapsV(strV);


        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(gapsH, gapsV, scoringScheme, AlignConfig<>(), NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Global alignment with two Segments
    {
        Dna5String strH = "ATGTAT";
        Dna5String strV = "ATATAGAT";

        Segment<Dna5String, PrefixSegment> prefixSegment(strH, 4);
        SEQAN_ASSERT_EQ(prefixSegment, "ATGT");

        Segment<Dna5String, SuffixSegment> suffixSegment(strV, 2);
        SEQAN_ASSERT_EQ(suffixSegment, "ATAGAT");

        Gaps<Segment<Dna5String, PrefixSegment>, ArrayGaps> gapsH(prefixSegment);
        Gaps<Segment<Dna5String, SuffixSegment>, ArrayGaps> gapsV(suffixSegment);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(gapsH, gapsV, scoringScheme, AlignConfig<>(), NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AT-G-T");
        SEQAN_ASSERT_EQ(ssV.str(), "ATAGAT");
    }

    // Global alignment score with Segment and DnaString
    {
        Dna5String strH = "ATGTAT";
        Dna5String strV = "ATAGAT";

        Segment<Dna5String, PrefixSegment> prefixSegment(strH, 4);
        SEQAN_ASSERT_EQ(prefixSegment, "ATGT");

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignmentScore(prefixSegment, strV, scoringScheme, AlignConfig<>(), NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);
    }

    // Global alignment with two Segments
    {
        Dna5String strH = "ATGTAT";
        Dna5String strV = "ATATAGAT";

        Segment<Dna5String, PrefixSegment> prefixSegment(strH, 4);
        SEQAN_ASSERT_EQ(prefixSegment, "ATGT");

        Segment<Dna5String, SuffixSegment> suffixSegment(strV, 2);
        SEQAN_ASSERT_EQ(suffixSegment, "ATAGAT");

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignmentScore(prefixSegment, suffixSegment, scoringScheme, AlignConfig<>(), NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);
    }
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_GLOBAL_H_
