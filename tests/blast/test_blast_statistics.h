// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Tests for the blast module
// ==========================================================================

#ifndef SEQAN_TESTS_TEST_BLAST_STATISTICS_H_
#define SEQAN_TESTS_TEST_BLAST_STATISTICS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/blast.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_blast_scoring_scheme_conversion)
{
    //TODO replace
//     Blosum62 scheme;
//     setScoreGapOpen(scheme, -11);
//     setScoreGapExtend(scheme, -1);
//
//     blastScoringScheme2seqanScoringScheme(scheme);
//
//     SEQAN_ASSERT_EQ(scoreGapOpen(scheme), -12);
//     SEQAN_ASSERT_EQ(scoreGapExtend(scheme), -1);
//
//     seqanScoringScheme2blastScoringScheme(scheme);
//
//     SEQAN_ASSERT_EQ(scoreGapOpen(scheme), -11);
//     SEQAN_ASSERT_EQ(scoreGapExtend(scheme), -1);
}

SEQAN_DEFINE_TEST(test_blast_scoring_adapter)
{
    // Blosum62 and some general stuff
    {
        typedef Blosum62 TScheme;
        BlastScoringScheme<TScheme> scheme;
        setScoreGapOpenBlast(scheme, -11);
        setScoreGapExtendBlast(scheme, -1);

        // TEST isValid()
        SEQAN_ASSERT(isValid(scheme));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(scheme), 0.267);
        SEQAN_ASSERT_EQ(getKappa(scheme),  0.041);
        SEQAN_ASSERT_EQ(getH(scheme),      0.140);
        SEQAN_ASSERT_EQ(getAlpha(scheme),  1.900);
        SEQAN_ASSERT_EQ(getBeta(scheme), -30.000);

        setScoreGapOpen(scheme, -1);
        setScoreGapExtend(scheme, -11);
        SEQAN_ASSERT(!isValid(scheme));// invalid parameters, because not the setScoreGap*Blast version

        // reset to valid
        setScoreGapOpenBlast(scheme, -11);
        setScoreGapExtendBlast(scheme, -1);

        // now check values again
        SEQAN_ASSERT(isValid(scheme));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(scheme), 0.267);
        SEQAN_ASSERT_EQ(getKappa(scheme),  0.041);
        SEQAN_ASSERT_EQ(getH(scheme),      0.140);
        SEQAN_ASSERT_EQ(getAlpha(scheme),  1.900);
        SEQAN_ASSERT_EQ(getBeta(scheme), -30.000);
    }

    // Blosum45
    {
        typedef Blosum45 TScheme;
        BlastScoringScheme<TScheme> scheme;
        setScoreGapOpenBlast(scheme, -11);
        setScoreGapExtendBlast(scheme, -3);

        // TEST isValid()
        SEQAN_ASSERT(isValid(scheme));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(scheme), 0.190);
        SEQAN_ASSERT_EQ(getKappa(scheme),  0.031);
        SEQAN_ASSERT_EQ(getH(scheme),      0.095);
        SEQAN_ASSERT_EQ(getAlpha(scheme),  2.000);
        SEQAN_ASSERT_EQ(getBeta(scheme), -38.000);
    }

    // Blosum80
    {
        typedef Blosum80 TScheme;
        BlastScoringScheme<TScheme> scheme;
        setScoreGapOpenBlast(scheme, -11);
        setScoreGapExtendBlast(scheme, -1);

        // TEST isValid()
        SEQAN_ASSERT(isValid(scheme));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(scheme), 0.314);
        SEQAN_ASSERT_EQ(getKappa(scheme),  0.095);
        SEQAN_ASSERT_EQ(getH(scheme),      0.350);
        SEQAN_ASSERT_EQ(getAlpha(scheme),  0.900);
        SEQAN_ASSERT_EQ(getBeta(scheme),  -9.000);
    }

    // Pam250
    {
        typedef Pam250 TScheme;
        BlastScoringScheme<TScheme> scheme;
        setScoreGapOpenBlast(scheme, -11);
        setScoreGapExtendBlast(scheme, -3);

        // TEST isValid()
        SEQAN_ASSERT(isValid(scheme));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(scheme), 0.174);
        SEQAN_ASSERT_EQ(getKappa(scheme),  0.020);
        SEQAN_ASSERT_EQ(getH(scheme),      0.070);
        SEQAN_ASSERT_EQ(getAlpha(scheme),  2.500);
        SEQAN_ASSERT_EQ(getBeta(scheme), -48.000);
    }

    // SimpleScore
    {
        typedef Score<int, Simple> TScheme;
        BlastScoringScheme<TScheme> scheme;
        setScoreMatch(scheme, 2);
        setScoreMismatch(scheme, -3);
        setScoreGapOpenBlast(scheme, -4);
        setScoreGapExtendBlast(scheme, -2);

        // TEST isValid()
        SEQAN_ASSERT(isValid(scheme));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(scheme), 0.610);
        SEQAN_ASSERT_EQ(getKappa(scheme),  0.350);
        SEQAN_ASSERT_EQ(getH(scheme),      0.680);
        SEQAN_ASSERT_EQ(getAlpha(scheme),  0.900);
        SEQAN_ASSERT_EQ(getBeta(scheme),  -3.000);
    }
}

SEQAN_DEFINE_TEST(test_blast_blastmatch_stats_and_score)
{
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<CharString, CharString, uint32_t, TAlign> TBlastMatch;


    TBlastMatch m;

    String<AminoAcid> src0 = "ARNDAYVBRNDCQFGCYVBQARNDCQEGEG";
    String<AminoAcid> src1 = "ARNAYVBRNDCCYCYVBQARNQEGEG";

    resize(rows(m.align), 2);
    assignSource(row(m.align, 0), src0);
    assignSource(row(m.align, 1), src1);

    typedef Blosum62 TScheme;
    BlastScoringScheme<TScheme> scheme;
    setScoreGapOpenBlast(scheme, -11);
    setScoreGapExtendBlast(scheme, -1);
    SEQAN_ASSERT(isValid(scheme));

    int score = globalAlignment(m.align, static_cast<TScheme>(scheme));
//     ARNDAYVBRNDCQFGCYVBQARNDCQEGEG
//     ||| ||||||||   ||||||||  |||||
//     ARN-AYVBRNDCCY-CYVBQARN--QEGEG

    SEQAN_ASSERT_EQ(score, 94);

    computeAlignmentStats(m.alignStats, m.align, static_cast<TScheme>(scheme));

    SEQAN_ASSERT_EQ(m.alignStats.alignmentScore, score);
//     SEQAN_ASSERT_EQ(m.alignLength, 30u);
    SEQAN_ASSERT_EQ(m.alignStats.numMatches, 24u);
    SEQAN_ASSERT_EQ(m.alignStats.numPositiveScores, 25u);
    SEQAN_ASSERT_EQ(m.alignStats.numMismatches, 2u);
    SEQAN_ASSERT_EQ(m.alignStats.numGapExtensions, 1u);
    SEQAN_ASSERT_EQ(m.alignStats.numGapOpens, 3u);
}

SEQAN_DEFINE_TEST(test_blast_blastmatch_bit_score_e_value)
{
    typedef Align<String<AminoAcid>, ArrayGaps> TAlign;
    typedef BlastMatch<CharString, CharString, uint32_t, TAlign> TBlastMatch;
    typedef Blosum62 TScheme;

    TBlastMatch m;

    String<AminoAcid> src0 =
    "SSITEEKHIPHKEQDKDAEFLSKEALKTHMTENVLQMDRRAVQDPSTSFLQLLKAKGLLG"
    "LPDYEVNLADVNSPGFRKVAYAQTKPRRLCFPNGGTRRGSFIMDTAVVVMVSLRYVNIGK"
    "VIFPGATDVSEGEDEFWAGLPQAYGCLATEFLCIHIAIYSWIHVQSSRYDDMNASVIRAK"
    "LNLPGGLTLIQARGNEKETI";
    String<AminoAcid> src1 = "VAYAQPRKLCYP";

    resize(rows(m.align), 2);
    assignSource(row(m.align, 0), src0);
    assignSource(row(m.align, 1), src1);

    BlastScoringScheme<TScheme> scheme;
    setScoreGapOpenBlast(scheme, -11);
    setScoreGapExtendBlast(scheme, -1);
    SEQAN_ASSERT(isValid(scheme));

    int score = localAlignment(m.align, static_cast<TScheme>(scheme));
//         VAYAQTKPRRLCFP
//         |||||  || || |
//         VAYAQ--PRKLCYP

    SEQAN_ASSERT_EQ(score, 48);

    computeAlignmentStats(m.alignStats, m.align, static_cast<TScheme>(scheme));

    SEQAN_ASSERT_EQ(m.alignStats.alignmentScore, score);
    SEQAN_ASSERT_EQ(m.alignStats.alignmentLength, 14u);
    SEQAN_ASSERT_EQ(m.alignStats.numMatches, 10u);
    SEQAN_ASSERT_EQ(m.alignStats.numPositiveScores, 12u);
    SEQAN_ASSERT_EQ(m.alignStats.numMismatches, 2u);
    SEQAN_ASSERT_EQ(m.alignStats.numGapExtensions, 1u);
    SEQAN_ASSERT_EQ(m.alignStats.numGapOpens, 1u);

    // same as previous test until here (just local)

    m.bitScore = computeBitScore(m.alignStats.alignmentScore, scheme);
    double epsilon = 1e-4;
    SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 23.0978), epsilon);

    m.eValue = computeEValue(m.alignStats.alignmentScore, length(src1), length(src0), scheme);
    epsilon = 1e-8; // evalues are smaller
    SEQAN_ASSERT_LEQ(std::abs(m.eValue - 0.000267348), epsilon);
}

#endif  // SEQAN_TESTS_TEST_BLAST_STATISTICS_H_
