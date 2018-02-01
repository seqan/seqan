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
// Tests for the more specialized global alignment algorithms Hirschberg
// and Myers-Hirschberg.
// ==========================================================================

#ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_GLOBAL_ALIGNMENT_SPECIALIZED_H_
#define SEQAN_TESTS_ALIGN_TEST_ALIGN_GLOBAL_ALIGNMENT_SPECIALIZED_H_

SEQAN_DEFINE_TEST(test_align_global_alignment_hirschberg_single_character)
{
    using namespace seqan;

    // Horizontal sequence has length 1.
    {
        Dna5String strH = "T";
        Dna5String strV = "AAT";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(3, -1, -1);

        int score = globalAlignment(align, scoringScheme, Hirschberg());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, 1);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "--T");
        SEQAN_ASSERT_EQ(ssV.str(), "AAT");
    }

    // Vertical sequence has length 1.
    {

        Dna5String strH = "AAT";
        Dna5String strV = "A";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(3, -1, -1);

        int score = globalAlignment(align, scoringScheme, Hirschberg());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, 1);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAT");
        SEQAN_ASSERT_EQ(ssV.str(), "A--");
    }

    // Both sequences have length 1.
    {

        Dna5String strH = "T";
        Dna5String strV = "A";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(3, -1, -1);

        int score = globalAlignment(align, scoringScheme, Hirschberg());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, -1);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "T");
        SEQAN_ASSERT_EQ(ssV.str(), "A");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_hirschberg_align)
{
    using namespace seqan;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(align, scoringScheme, Hirschberg());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, 14);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTTTT----G");
        SEQAN_ASSERT_EQ(ssV.str(), "A--ATTTTTTTTTTGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_hirschberg_gaps)
{
    using namespace seqan;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        DnaString strV = "AATTTTTTTTTTGGGGG";

        Gaps<Dna5String, ArrayGaps> gapsH(strH);
        Gaps<DnaString, ArrayGaps> gapsV(strV);

        Score<int, Simple> scoringScheme(2, -1, -1);

        int score = globalAlignment(gapsH, gapsV, scoringScheme, Hirschberg());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, 14);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTTTT----G");
        SEQAN_ASSERT_EQ(ssV.str(), "A--ATTTTTTTTTTGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_hirschberg_fragments)
{
    // TODO(holtgrew): Test when implemented!
}

SEQAN_DEFINE_TEST(test_align_global_alignment_hirschberg_graph)
{
    // TODO(holtgrew): Test when implemented!
}

SEQAN_DEFINE_TEST(test_align_global_alignment_myers_hirschberg_align)
{
    using namespace seqan;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        Dna5String strV = "AATTTTTTTTTTGGGGG";

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        int score = globalAlignment(align, MyersHirschberg());
        // std::cerr << align << "\n";

        SEQAN_ASSERT_EQ(score, -8);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTT--TTG");
        SEQAN_ASSERT_EQ(ssV.str(), "AATTTTTTTTTTGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_myers_hirschberg_gaps)
{
    using namespace seqan;

    // More or less simple alignment.
    {
        Dna5String strH = "AAAAAATTTTTTTTG";
        DnaString strV = "AATTTTTTTTTTGGGGG";

        Gaps<Dna5String, ArrayGaps> gapsH(strH);
        Gaps<DnaString, ArrayGaps> gapsV(strV);

        int score = globalAlignment(gapsH, gapsV, MyersHirschberg());
        // std::cerr << gapsH << "\n" << gapsV << "\n";

        SEQAN_ASSERT_EQ(score, -8);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "AAAAAATTTTTT--TTG");
        SEQAN_ASSERT_EQ(ssV.str(), "AATTTTTTTTTTGGGGG");
    }
}

SEQAN_DEFINE_TEST(test_align_global_alignment_myers_hirschberg_fragments)
{
    // TODO(holtgrew): Test when implemented!
}

SEQAN_DEFINE_TEST(test_align_global_alignment_myers_hirschberg_graph)
{
    // TODO(holtgrew): Test when implemented!
}

SEQAN_DEFINE_TEST(test_align_global_alignment_score_hirschberg)
{
    using namespace seqan;

    Dna5String strH = "ATGT";
    DnaString strV = "ATAGAT";

    StringSet<Dna5String> strings;
    appendValue(strings, strH);
    appendValue(strings, strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    int res = 0;

    res = globalAlignmentScore(strH, strV, scoringScheme, Hirschberg());
    SEQAN_ASSERT_EQ(res, 6);


    res = globalAlignmentScore(strings, scoringScheme, Hirschberg());
    SEQAN_ASSERT_EQ(res, 6);
}

SEQAN_DEFINE_TEST(test_align_global_alignment_score_myers)
{
    using namespace seqan;

    Dna5String strH = "AAAAAATTTTTTTTG";
    DnaString strV = "AATTTTTTTTTTGGGGG";

    StringSet<Dna5String> strings;
    appendValue(strings, strH);
    appendValue(strings, strV);

    int res = 0;

    res = globalAlignmentScore(strH, strV, MyersBitVector());
    SEQAN_ASSERT_EQ(res, -8);

    res = globalAlignmentScore(strings, MyersBitVector());
    SEQAN_ASSERT_EQ(res, -8);
}

SEQAN_DEFINE_TEST(test_align_global_alignment_score_myers_hirschberg)
{
    using namespace seqan;

    Dna5String strH = "AAAAAATTTTTTTTG";
    DnaString strV = "AATTTTTTTTTTGGGGG";

    StringSet<Dna5String> strings;
    appendValue(strings, strH);
    appendValue(strings, strV);

    int res = 0;

    res = globalAlignmentScore(strH, strV, MyersHirschberg());
    SEQAN_ASSERT_EQ(res, -8);

    res = globalAlignmentScore(strings, MyersHirschberg());
    SEQAN_ASSERT_EQ(res, -8);
}

#endif  // #ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_GLOBAL_ALIGNMENT_SPECIALIZED_H_
