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
// Tests for the various variants of globalAlignmentScore().
// ==========================================================================

#ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_GLOBAL_ALIGNMENT_SCORE_H_
#define SEQAN_TESTS_ALIGN_TEST_ALIGN_GLOBAL_ALIGNMENT_SCORE_H_

#include <sstream>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/score.h>


SEQAN_DEFINE_TEST(test_align_global_alignment_score_nw)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    Dna5String strH = "ATGT";
    DnaString strV = "ATAGAT";

    StringSet<Dna5String> strings;
    appendValue(strings, strH);
    appendValue(strings, strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    int res = 0;

    res = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strH, strV, scoringScheme, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strH, strV, scoringScheme, alignConfig);
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strH, strV, scoringScheme);
    SEQAN_ASSERT_EQ(res, 6);


    res = globalAlignmentScore(strings, scoringScheme, alignConfig, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strings, scoringScheme, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strings, scoringScheme, alignConfig);
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strings, scoringScheme);
    SEQAN_ASSERT_EQ(res, 6);
}

SEQAN_DEFINE_TEST(test_align_global_alignment_score_banded_nw)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    Dna5String strH = "ATGT";
    DnaString strV = "ATAGAT";

    StringSet<Dna5String> strings;
    appendValue(strings, strH);
    appendValue(strings, strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    int res = 0;

    res = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, -2, 2, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strH, strV, scoringScheme, -2, 2, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, -2, 2);
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strH, strV, scoringScheme, -2, 2);
    SEQAN_ASSERT_EQ(res, 6);


    res = globalAlignmentScore(strings, scoringScheme, alignConfig, -2, 2, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strings, scoringScheme, -2, 2, NeedlemanWunsch());
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strings, scoringScheme, alignConfig, -2, 2);
    SEQAN_ASSERT_EQ(res, 6);

    res = globalAlignmentScore(strings, scoringScheme, -2, 2);
    SEQAN_ASSERT_EQ(res, 6);
}

SEQAN_DEFINE_TEST(test_align_global_alignment_score_gotoh)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    Dna5String strH = "ATGT";
    DnaString strV = "ATAGAT";

    StringSet<Dna5String> strings;
    appendValue(strings, strH);
    appendValue(strings, strV);

    Score<int, Simple> scoringScheme(2, -1, -1, -3);

    int res = 0;

    res = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, Gotoh());
    SEQAN_ASSERT_EQ(res, 4);

    res = globalAlignmentScore(strH, strV, scoringScheme, Gotoh());
    SEQAN_ASSERT_EQ(res, 2);

    res = globalAlignmentScore(strH, strV, scoringScheme, alignConfig);
    SEQAN_ASSERT_EQ(res, 4);

    res = globalAlignmentScore(strH, strV, scoringScheme);
    SEQAN_ASSERT_EQ(res, 2);


    res = globalAlignmentScore(strings, scoringScheme, alignConfig, Gotoh());
    SEQAN_ASSERT_EQ(res, 4);

    res = globalAlignmentScore(strings, scoringScheme, Gotoh());
    SEQAN_ASSERT_EQ(res, 2);

    res = globalAlignmentScore(strings, scoringScheme, alignConfig);
    SEQAN_ASSERT_EQ(res, 4);

    res = globalAlignmentScore(strings, scoringScheme);
    SEQAN_ASSERT_EQ(res, 2);
}

SEQAN_DEFINE_TEST(test_align_global_alignment_score_banded_gotoh)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    Dna5String strH = "ATGT";
    DnaString strV = "ATAGAT";

    StringSet<Dna5String> strings;
    appendValue(strings, strH);
    appendValue(strings, strV);

    Score<int, Simple> scoringScheme(2, -1, -1, -3);

    int res = 0;

    res = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, -2, 2, Gotoh());
    SEQAN_ASSERT_EQ(res, 2);

    res = globalAlignmentScore(strH, strV, scoringScheme, -2, 2, Gotoh());
    SEQAN_ASSERT_EQ(res, 2);

    res = globalAlignmentScore(strH, strV, scoringScheme, alignConfig, -2, 2);
    SEQAN_ASSERT_EQ(res, 2);

    res = globalAlignmentScore(strH, strV, scoringScheme, -2, 2);
    SEQAN_ASSERT_EQ(res, 2);


    res = globalAlignmentScore(strings, scoringScheme, alignConfig, -2, 2, Gotoh());
    SEQAN_ASSERT_EQ(res, 2);

    res = globalAlignmentScore(strings, scoringScheme, -2, 2, Gotoh());
    SEQAN_ASSERT_EQ(res, 2);

    res = globalAlignmentScore(strings, scoringScheme, alignConfig, -2, 2);
    SEQAN_ASSERT_EQ(res, 2);

    res = globalAlignmentScore(strings, scoringScheme, -2, 2);
    SEQAN_ASSERT_EQ(res, 2);
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

#endif  // #ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_GLOBAL_ALIGNMENT_SCORE_H_
