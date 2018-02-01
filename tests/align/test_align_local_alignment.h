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
// Tests for the local alignment algorithms, both unbanded and banded.
// ==========================================================================

#ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_LOCAL_ALIGNMENT_H_
#define SEQAN_TESTS_ALIGN_TEST_ALIGN_LOCAL_ALIGNMENT_H_

#include <sstream>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/score.h>

// TODO(holtgrew): Tests for local alignment using Alignment Graphs.

// ==========================================================================
// Local Alignment Computation
// ==========================================================================

SEQAN_DEFINE_TEST(test_align_local_alignment_align)
{
    using namespace seqan;

    {
        Dna5String strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -2, -2);

        int score = localAlignment(align, scoringScheme);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTTA-GCT");
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_gaps)
{
    using namespace seqan;

    {
        DnaString strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Gaps<DnaString> gapsH(strH);
        Gaps<Dna5String> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2, -2);

        int score = localAlignment(gapsH, gapsV, scoringScheme);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTTA-GCT");
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_graph)
{
    // TODO(holtgrew): Test after this is written.
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    {
        Dna5String strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        SimpleScore scoringScheme(2, -1, -2, -2);

        int score = localAlignment(alignmentGraph, scoringScheme);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ss;
        ss << alignmentGraph;

        // Note that the non-overlapping of the segments not part of the
        // sequence is intended: They simply do not take part in the local
        // alignment.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :    .   \n"
                   << "        ----GGGGCTTAAGCT------TGGGG\n"
                   << "                ||| ||||           \n"
                   << "        AAAA----CTT-AGCTCTAAAA-----\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_fragment)
{
    using namespace seqan;

    // TODO(holtgrew): Test after this is written.
}

SEQAN_DEFINE_TEST(test_align_local_alignment_banded_align)
{
    using namespace seqan;

    {
        Dna5String strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -2, -2);

        int score = localAlignment(align, scoringScheme, -10, 10);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTT-AGCT");
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_banded_gaps)
{
    using namespace seqan;

    {
        DnaString strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Gaps<DnaString> gapsH(strH);
        Gaps<Dna5String> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2, -2);

        int score = localAlignment(gapsH, gapsV, scoringScheme, -6, 6);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTT-AGCT");
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_banded_graph)
{
    using namespace seqan;

    // TODO(holtgrew): Test after this is written.
}

SEQAN_DEFINE_TEST(test_align_local_alignment_banded_fragment)
{
    using namespace seqan;

    // TODO(holtgrew): Test after this is written.
}

// ==========================================================================
// Local Alignment Enumeration
// ==========================================================================

SEQAN_DEFINE_TEST(test_align_local_alignment_enumeration_align)
{
    using namespace seqan;

    // Example from original testLocalAlign2 test.
    {
        std::stringstream ssH, ssV;

        Dna5String strH("ATAAGCGTCTCG");
        Dna5String strV("TCATAGAGTTGC");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -2, -2);
        int cutoff = 5;

        LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 9);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "ATAAGCGT");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "ATA-GAGT");

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 5);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "TC-TCG");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "TCATAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(align, enumerator));
    }

    // First example from original testBandedLocalAlign test, presumably from Birte.
    {
        std::stringstream ssH, ssV;

        Dna5String strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -2, -2);
        int cutoff = 5;

        LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 12);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT(ssV.str() == "CTTA-GCT" || ssV.str() == "CTT-AGCT");

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 10);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "GCT-TAA");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "GCTCTAA");

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 10);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "AAGCTTGG");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "AAACTTAG");

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 7);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "CTTAA");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "CTAAA");

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 5);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "AAGCTTG");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "AACTTAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(align, enumerator));
    }

    // Second example from original testBandedLocalAlign test, presumably from Birte.
    {

        std::stringstream ssH, ssV;

        Dna5String strH("GCAGAATTAAGGAGGATTACAAGTGGGAATTTGAAGAGCTTTTGAAATCC");
        Dna5String strV("CGGTTGAGCAGAACTTGGGCTACGAGACTCCCCCCGAGGAATTTGAAGGCTTTCTTCAAATCCAAAAGCA");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(1, -9, -9, -9);
        int cutoff = 7;

        LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 11);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "GGAATTTGAAG");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "GGAATTTGAAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(align, enumerator));
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_enumeration_gaps)
{
    using namespace seqan;

    // Example from original testLocalAlign2 test.
    {
        std::stringstream ssH, ssV;

        Dna5String strH("ATAAGCGTCTCG");
        DnaString strV("TCATAGAGTTGC");

        Gaps<Dna5String, ArrayGaps> gapsH(strH);
        Gaps<DnaString, ArrayGaps> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2, -2);
        int cutoff = 5;

        LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 9);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "ATAAGCGT");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "ATA-GAGT");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 5);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "TC-TCG");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "TCATAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(gapsH, gapsV, enumerator));
    }

    // First example from original testBandedLocalAlign test, presumably from Birte.
    {
        std::stringstream ssH, ssV;

        Dna5String strH("GGGGCTTAAGCTTGGGG");
        DnaString strV("AAAACTTAGCTCTAAAA");

        Gaps<Dna5String, ArrayGaps> gapsH(strH);
        Gaps<DnaString, ArrayGaps> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2, -2);
        int cutoff = 5;

        LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 12);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT(ssV.str() == "CTTA-GCT" || ssV.str() == "CTT-AGCT");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 10);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "GCT-TAA");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "GCTCTAA");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 10);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "AAGCTTGG");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "AAACTTAG");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 7);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "CTTAA");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "CTAAA");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 5);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "AAGCTTG");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "AACTTAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(gapsH, gapsV, enumerator));
    }

    // Second example from original testBandedLocalAlign test, presumably from Birte.
    {

        std::stringstream ssH, ssV;

        Dna5String strH("GCAGAATTAAGGAGGATTACAAGTGGGAATTTGAAGAGCTTTTGAAATCC");
        DnaString strV("CGGTTGAGCAGAACTTGGGCTACGAGACTCCCCCCGAGGAATTTGAAGGCTTTCTTCAAATCCAAAAGCA");

        Gaps<Dna5String, ArrayGaps> gapsH(strH);
        Gaps<DnaString, ArrayGaps> gapsV(strV);

        SimpleScore scoringScheme(1, -9, -9, -9);
        int cutoff = 7;

        LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 11);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "GGAATTTGAAG");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "GGAATTTGAAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(gapsH, gapsV, enumerator));
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_enumeration_fragment)
{
    using namespace seqan;

    // TODO(holtgrew): Test after this is written.
}

SEQAN_DEFINE_TEST(test_align_local_alignment_enumeration_banded_align)
{
    using namespace seqan;

    // Example from original testLocalAlign2 test.
    {
        std::stringstream ssH, ssV;

        Dna5String strH("ATAAGCGTCTCG");
        Dna5String strV("TCATAGAGTTGC");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -2, -2);
        int cutoff = 5;

        LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, -4, 4, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 9);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "ATAAGCGT");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT(ssV.str() == "ATA-GAGT" || ssV.str() == "AT-AGAGT");

        SEQAN_ASSERT_NOT(nextLocalAlignment(align, enumerator));
    }

    // First example from original testBandedLocalAlign test, presumably from Birte.
    {
        std::stringstream ssH, ssV;

        Dna5String strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -2, -2);
        int cutoff = 5;

        LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, -10, 10, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 12);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT(ssV.str() == "CTTA-GCT" || ssV.str() == "CTT-AGCT");

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 10);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "AAGCTTGG");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "AAACTTAG");

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 10);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "GCT-TAA");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "GCTCTAA");

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 7);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "CTTAA");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "CTAAA");

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 5);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "AAGCTTG");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "AACTTAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(align, enumerator));
    }

    // Second example from original testBandedLocalAlign test, presumably from Birte.
    {

        std::stringstream ssH, ssV;

        Dna5String strH("GCAGAATTAAGGAGGATTACAAGTGGGAATTTGAAGAGCTTTTGAAATCC");
        Dna5String strV("CGGTTGAGCAGAACTTGGGCTACGAGACTCCCCCCGAGGAATTTGAAGGCTTTCTTCAAATCCAAAAGCA");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(1, -9, -9, -9);
        int cutoff = 7;

        LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, -20, 0, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(align, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 11);
        ssH.clear(); ssH.str(""); ssH << row(align, 0);
        SEQAN_ASSERT_EQ(ssH.str(), "GGAATTTGAAG");
        ssV.clear(); ssV.str(""); ssV << row(align, 1);
        SEQAN_ASSERT_EQ(ssV.str(), "GGAATTTGAAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(align, enumerator));
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_enumeration_banded_gaps)
{
    using namespace seqan;

    // Example from original testLocalAlign2 test.
    {
        std::stringstream ssH, ssV;

        Dna5String strH("ATAAGCGTCTCG");
        DnaString strV("TCATAGAGTTGC");

        Gaps<Dna5String, ArrayGaps> gapsH(strH);
        Gaps<DnaString, ArrayGaps> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2, -2);
        int cutoff = 5;

        LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, -4, 4, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 9);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "ATAAGCGT");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "AT-AGAGT");

        SEQAN_ASSERT_NOT(nextLocalAlignment(gapsH, gapsV, enumerator));
    }

    // First example from original testBandedLocalAlign test, presumably from Birte.
    {
        std::stringstream ssH, ssV;

        Dna5String strH("GGGGCTTAAGCTTGGGG");
        DnaString strV("AAAACTTAGCTCTAAAA");

        Gaps<Dna5String, ArrayGaps> gapsH(strH);
        Gaps<DnaString, ArrayGaps> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2, -2);
        int cutoff = 5;

        LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, -10, 10, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 12);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT(ssV.str() == "CTTA-GCT" || ssV.str() == "CTT-AGCT");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 10);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "AAGCTTGG");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "AAACTTAG");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 10);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "GCT-TAA");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "GCTCTAA");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 7);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "CTTAA");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "CTAAA");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 5);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "AAGCTTG");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "AACTTAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(gapsH, gapsV, enumerator));
    }

    // Second example from original testBandedLocalAlign test, presumably from Birte.
    {

        std::stringstream ssH, ssV;

        Dna5String strH("GCAGAATTAAGGAGGATTACAAGTGGGAATTTGAAGAGCTTTTGAAATCC");
        DnaString strV("CGGTTGAGCAGAACTTGGGCTACGAGACTCCCCCCGAGGAATTTGAAGGCTTTCTTCAAATCCAAAAGCA");

        Gaps<Dna5String, ArrayGaps> gapsH(strH);
        Gaps<DnaString, ArrayGaps> gapsV(strV);

        SimpleScore scoringScheme(1, -9, -9, -9);
        int cutoff = 7;

        LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, -20, 0, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 11);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "GGAATTTGAAG");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "GGAATTTGAAG");

        SEQAN_ASSERT_NOT(nextLocalAlignment(gapsH, gapsV, enumerator));
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_enumeration_banded_fragment)
{
    using namespace seqan;

    // TODO(holtgrew): Test after this is written.
}

#endif  // #ifndef SEQAN_TESTS_ALIGN_TEST_ALIGN_LOCAL_ALIGNMENT_H_
