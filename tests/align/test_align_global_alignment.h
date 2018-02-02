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
// Tests for the various global alignment algorithms.
//
// We test all result types with the "full" global alignment interface (using
// the maximal number of parameters) with comprehensive tests, i.e. with
// multiple inputs
//
// Then, we test the shorter interfaces with one input value only.
// ==========================================================================

#ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_GLOBAL_ALIGNMENT_H_
#define SEQAN_TESTS_ALIGN_TEST_ALIGN_GLOBAL_ALIGNMENT_H_

#include <sstream>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/score.h>

// ==========================================================================
// Long Interfaces
// ==========================================================================

SEQAN_DEFINE_TEST(test_align_global_alignment_align_gaps_free_top_left_right_bottom_nw)
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
        // std::cerr << align << "\n";
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "GGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, 7);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "----GGGGG----");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTTTT";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, 16);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTTTT--------");
        SEQAN_ASSERT_EQ(ssV.str(), "------TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_align_gaps_free_notop_left_noright_bottom_nw)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, 14);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTTTT----G");
        SEQAN_ASSERT_EQ(ssV.str(), "--AATTTTTTTTTTGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_gaps_gaps_free_top_left_right_bottom_nw)
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
        Dna5String strV = "GGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 7);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAGGGGTTTT");
        SEQAN_ASSERT_EQ(ssV.str(), "----GGGGG----");
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        DnaString strH = "AAAAAATTTTTTTT";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 16);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTTTT--------");
        SEQAN_ASSERT_EQ(ssV.str(), "------TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_gaps_gaps_free_notop_left_noright_bottom_nw)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;

    // More or less simple alignment.
    {
        DnaString strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 14);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTTTT----G");
        SEQAN_ASSERT_EQ(ssV.str(), "--AATTTTTTTTTTGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_graph_gaps_free_top_left_right_bottom_nw)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        // TODO(holtgrew): We'd rather have a function to extract the alignment in a more compact form, e.g. get an Align<> object?
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
        Dna5String strV = "GGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 7);

        // TODO(holtgrew): We'd rather have a function to extract the alignment in a more compact form, e.g. get an Align<> object?
        std::stringstream ss;
        ss << alignmentGraph;
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    \n"
                   << "        AAAAAGGGGTTTT\n"
                   << "             ||||    \n"
                   << "        ----GGGGG----\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTTTT";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 16);

        // TODO(holtgrew): We'd rather have a function to extract the alignment in a more compact form, e.g. get an Align<> object?
        std::stringstream ss;
        ss << alignmentGraph;
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :   \n"
                   << "        AAAAAATTTTTTTT--------\n"
                   << "              ||||||||        \n"
                   << "        ------TTTTTTTTGGGGGGGG\n\n\n";
        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_graph_gaps_free_notop_left_noright_bottom_nw)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 14);

        // TODO(holtgrew): We'd rather have a function to extract the alignment in a more compact form, e.g. get an Align<> object?
        std::stringstream ss;
        ss << alignmentGraph;
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .     \n"
                   << "        AAAAAATTTTTTTT----G\n"
                   << "          ||  ||||||||    |\n"
                   << "        --AATTTTTTTTTTGGGGG\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_fragments_gaps_free_top_left_right_bottom_nw)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        TFragmentString fragments;

        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 6);

        SEQAN_ASSERT_EQ(length(fragments), 3u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 5, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 3, 1));
        SEQAN_ASSERT(fragments[2] == TFragment(0, 0, 1, 0, 2));
    }

    // Alignment with both leading and trailing gaps in one row.
    {
        Dna5String strH = "AAAAAGGGGTTTT";
        Dna5String strV = "GGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        TFragmentString fragments;

        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 7);

        SEQAN_ASSERT_EQ(length(fragments), 1u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 4, 1, 0, 5));
    }

    // Alignment with both leading and trailing gaps in different rows.
    {
        Dna5String strH = "AAAAAATTTTTTTT";
        Dna5String strV = "TTTTTTTTGGGGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        TFragmentString fragments;

        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 16);

        SEQAN_ASSERT_EQ(length(fragments), 1u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 6, 1, 0, 8));
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_fragments_gaps_free_notop_left_noright_bottom_nw)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        TFragmentString fragments;

        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, NeedlemanWunsch());

        SEQAN_ASSERT_EQ(score, 14);

        // TODO(holtgrew): Implement operator<< for Fragment<>?
        SEQAN_ASSERT_EQ(length(fragments), 2u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 14, 1, 16, 1));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 0, 12));
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_align_gaps_free_top_left_right_bottom_gotoh)
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
        // std::cerr << align << "\n";

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
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, 13);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTGGG-----");
        SEQAN_ASSERT_EQ(ssV.str(), "---TTTTTTTTGGGGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_align_gaps_free_notop_left_noright_bottom_gotoh)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -4);

        int score = globalAlignment(align, scoringScheme, alignConfig, Gotoh());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, 8);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTTTTG----");
        SEQAN_ASSERT_EQ(ssV.str(), "--AATTTTTTTTTTGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_gaps_gaps_free_top_left_right_bottom_gotoh)
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

SEQAN_DEFINE_TEST(test_align_global_alignment_gaps_gaps_free_notop_left_noright_bottom_gotoh)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;

    // More or less simple alignment.
    {
        DnaString strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        Gaps<DnaString, ArrayGaps> gapsH(strH);
        Gaps<Dna5String, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -4);

        int score = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 8);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTTTTG----");
        SEQAN_ASSERT_EQ(ssV.str(), "--AATTTTTTTTTTGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_graph_gaps_free_top_left_right_bottom_gotoh)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);

        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 4);

        // TODO(holtgrew): We'd rather have a function to extract the alignment in a more compact form, e.g. get an Align<> object?
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

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);

        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 7);

        // TODO(holtgrew): We'd rather have a function to extract the alignment in a more compact form, e.g. get an Align<> object?
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

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);

        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 13);

        // TODO(holtgrew): We'd rather have a function to extract the alignment in a more compact form, e.g. get an Align<> object?
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

SEQAN_DEFINE_TEST(test_align_global_alignment_graph_gaps_free_notop_left_noright_bottom_gotoh)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        Score<int, Simple> scoringScheme(2, -1, -1, -4);

        int score = globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 8);

        // TODO(holtgrew): We'd rather have a function to extract the alignment in a more compact form, e.g. get an Align<> object?
        std::stringstream ss;
        ss << alignmentGraph;
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .     \n"
                   << "        AAAAAATTTTTTTTG----\n"
                   << "          ||  |||||||||    \n"
                   << "        --AATTTTTTTTTTGGGGG\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_fragments_gaps_free_top_left_right_bottom_gotoh)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    // Simple alignment without any leading or trailing gaps.
    {
        Dna5String strH = "ATGT";
        Dna5String strV = "ATAGAT";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -3);

        TFragmentString fragments;

        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 4);

        SEQAN_ASSERT_EQ(length(fragments), 1u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 0, 1, 4, 2));
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

        SEQAN_ASSERT_EQ(score, 7);

        SEQAN_ASSERT_EQ(length(fragments), 2u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 8, 1, 3, 3));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 2, 1, 0, 3));
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

        SEQAN_ASSERT_EQ(score, 13);

        SEQAN_ASSERT_EQ(length(fragments), 1u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 3, 1, 0, 11));
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_fragments_gaps_free_notop_left_noright_bottom_gotoh)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;

    typedef StringSet<Dna5String> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        Score<int, Simple> scoringScheme(2, -1, -1, -4);

        TFragmentString fragments;

        int score = globalAlignment(fragments, strings, scoringScheme, alignConfig, Gotoh());

        SEQAN_ASSERT_EQ(score, 8);

        SEQAN_ASSERT_EQ(length(fragments), 1u);
        SEQAN_ASSERT(fragments[0] == TFragment(0, 2, 1, 0, 13));
    }
}

// ==========================================================================
// Shorter Interfaces
// ==========================================================================

// We only test the resulting scores and compare with the result from the
// longer interfaces (those are already tested above).

SEQAN_DEFINE_TEST(test_align_global_alignment_shorter_interfaces_nw)
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

SEQAN_DEFINE_TEST(test_align_global_alignment_shorter_interfaces_gotoh)
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

#endif  // #ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_GLOBAL_ALIGNMENT_H_
