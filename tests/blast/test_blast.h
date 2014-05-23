// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2014, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Tests for basic/translation.h
// ==========================================================================

#ifndef SEQAN_EXTRAS_TESTS_BASIC_TEST_BLAST_H_
#define SEQAN_EXTRAS_TESTS_BASIC_TEST_BLAST_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/blast.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_blast_scoring_scheme_conversion)
{
    Blosum62 scheme;
    setScoreGapOpen(scheme, -11);
    setScoreGapExtend(scheme, -1);

    blastScoringScheme2seqanScoringScheme(scheme);

    SEQAN_ASSERT_EQ(scoreGapOpen(scheme), -12);
    SEQAN_ASSERT_EQ(scoreGapExtend(scheme), -1);

    seqanScoringScheme2blastScoringScheme(scheme);

    SEQAN_ASSERT_EQ(scoreGapOpen(scheme), -11);
    SEQAN_ASSERT_EQ(scoreGapExtend(scheme), -1);
}

SEQAN_DEFINE_TEST(test_blast_scoring_adapter)
{
    // Blosum62 and some general stuff
    {
        typedef Blosum62 TScheme;
        TScheme scheme;
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -1);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.267);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.041);
        SEQAN_ASSERT_EQ(getH(adapter),      0.140);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  1.900);
        SEQAN_ASSERT_EQ(getBeta(adapter), -30.000);

        setScoreGapOpen(scheme, -1);
        setScoreGapExtend(scheme, -11); // invalid parameters
        // TEST assignScoreScheme
        bool res = assignScoreScheme(adapter, scheme);
        SEQAN_ASSERT(!res);

        // TEST getScoreScheme
        Blosum62 const & schemeInt = getScoreScheme(adapter);
        SEQAN_ASSERT_EQ(scoreGapOpen(schemeInt), -1);
        SEQAN_ASSERT_EQ(scoreGapExtend(schemeInt), -11);

        // reset to valid
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -1);
        blastScoringScheme2seqanScoringScheme(scheme);
        assignScoreScheme(adapter, scheme);
        // now check values again
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.267);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.041);
        SEQAN_ASSERT_EQ(getH(adapter),      0.140);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  1.900);
        SEQAN_ASSERT_EQ(getBeta(adapter), -30.000);
    }

    // Blosum45
    {
        typedef Blosum45 TScheme;
        TScheme scheme;
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -3);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.190);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.031);
        SEQAN_ASSERT_EQ(getH(adapter),      0.095);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  2.000);
        SEQAN_ASSERT_EQ(getBeta(adapter), -38.000);
    }

    // Blosum80
    {
        typedef Blosum80 TScheme;
        TScheme scheme;
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -1);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.314);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.095);
        SEQAN_ASSERT_EQ(getH(adapter),      0.350);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  0.900);
        SEQAN_ASSERT_EQ(getBeta(adapter),  -9.000);
    }

    // Pam250
    {
        typedef Pam250 TScheme;
        TScheme scheme;
        setScoreGapOpen(scheme, -11);
        setScoreGapExtend(scheme, -3);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.174);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.020);
        SEQAN_ASSERT_EQ(getH(adapter),      0.070);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  2.500);
        SEQAN_ASSERT_EQ(getBeta(adapter), -48.000);
    }

    // SimpleScore
    {
        typedef Score<int, Simple> TScheme;
        TScheme scheme;
        setScoreMatch(scheme, 2);
        setScoreMismatch(scheme, -3);
        setScoreGapOpen(scheme, -4);
        setScoreGapExtend(scheme, -2);
        blastScoringScheme2seqanScoringScheme(scheme);

        // TEST CONSTRUCTOR
        BlastScoringAdapter<TScheme> adapter(scheme);
        // TEST isValid()
        SEQAN_ASSERT(isValid(adapter));
        // TEST get*
        SEQAN_ASSERT_EQ(getLambda(adapter), 0.610);
        SEQAN_ASSERT_EQ(getKappa(adapter),  0.350);
        SEQAN_ASSERT_EQ(getH(adapter),      0.680);
        SEQAN_ASSERT_EQ(getAlpha(adapter),  0.900);
        SEQAN_ASSERT_EQ(getBeta(adapter),  -3.000);
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
    TScheme scheme;
    setScoreGapOpen(scheme, -11);
    setScoreGapExtend(scheme, -1);
    blastScoringScheme2seqanScoringScheme(scheme);

    int score = globalAlignment(m.align, scheme);
//     ARNDAYVBRNDCQFGCYVBQARNDCQEGEG
//     ||| ||||||||   ||||||||  |||||
//     ARN-AYVBRNDCCY-CYVBQARN--QEGEG

    SEQAN_ASSERT_EQ(score, 94);

    calcStatsAndScore(m, scheme);

    SEQAN_ASSERT_EQ(m.score, score);
    SEQAN_ASSERT_EQ(m.aliLength, 30u);
    SEQAN_ASSERT_EQ(m.identities, 24u);
    SEQAN_ASSERT_EQ(m.positives, 25u);
    SEQAN_ASSERT_EQ(m.mismatches, 2u);
    SEQAN_ASSERT_EQ(m.gaps, 4u);
    SEQAN_ASSERT_EQ(m.gapOpenings, 3u);
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

    TScheme scheme;
    setScoreGapOpen(scheme, -11);
    setScoreGapExtend(scheme, -1);
    blastScoringScheme2seqanScoringScheme(scheme);

    int score = localAlignment(m.align, scheme);
//         VAYAQTKPRRLCFP
//         |||||  || || |
//         VAYAQ--PRKLCYP

    SEQAN_ASSERT_EQ(score, 48);

    calcStatsAndScore(m, scheme);

    SEQAN_ASSERT_EQ(m.score, score);
    SEQAN_ASSERT_EQ(m.aliLength, 14u);
    SEQAN_ASSERT_EQ(m.identities, 10u);
    SEQAN_ASSERT_EQ(m.positives, 12u);
    SEQAN_ASSERT_EQ(m.mismatches, 2u);
    SEQAN_ASSERT_EQ(m.gaps, 2u);
    SEQAN_ASSERT_EQ(m.gapOpenings, 1u);

    // same as previous test until here (just local)

    BlastScoringAdapter<TScheme> adapter(scheme);
    SEQAN_ASSERT(isValid(adapter));

    calcBitScoreAndEValue(m.bitScore, m.eVal, m.score, length(src0),
                          length(src1), adapter);

    double epsilon = 1e-4;
    SEQAN_ASSERT_LEQ(std::abs(m.bitScore - 23.0978), epsilon);
    epsilon = 1e-8; // more important on evalue
    SEQAN_ASSERT_LEQ(std::abs(m.eVal - 0.000267348), epsilon);
}




#endif  // SEQAN_EXTRAS_TESTS_BASIC_TEST_BLAST_H_
