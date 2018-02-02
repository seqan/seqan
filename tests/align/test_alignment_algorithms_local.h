// ==========================================================================
//                     test_alignment_algorithms_local.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_LOCAL_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_LOCAL_H_

#include <sstream>
#include <seqan/basic.h>

#include <seqan/score.h>
#include <seqan/align.h>


SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_linear)
{
    using namespace seqan;

    {
        Dna5String strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -2);
        int score = localAlignment(align, scoringScheme);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTT-AGCT");
    }

    {
        Dna5String strH("GGGGGGGGG");
        Dna5String strV("CCCCCCCCC");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -2);
        int score = localAlignment(align, scoringScheme);

        SEQAN_ASSERT_EQ(score, 0);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_local_linear)
{
    using namespace seqan;

    {
        DnaString strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Gaps<DnaString> gapsH(strH);
        Gaps<Dna5String> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2);

        String<TraceSegment_<unsigned, unsigned> > traces;
        int score = localAlignment(gapsH, gapsV, scoringScheme);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTT-AGCT");
    }

    {
        DnaString strH("GGGGGGGGG");
        Dna5String strV("CCCCCCCCC");

        Gaps<DnaString> gapsH(strH);
        Gaps<Dna5String> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2);
        int score = localAlignment(gapsH, gapsV, scoringScheme);

        SEQAN_ASSERT_EQ(score, 0);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_local_linear)
{
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

        SimpleScore scoringScheme(2, -1, -2);
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
                   << "        GGGG----CTTAAGCTTGGGG------\n"
                   << "                ||| ||||           \n"
                   << "        ----AAAACTT-AGCT-----CTAAAA\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }

    {
        Dna5String strH("GGGGGGGG");
        Dna5String strV("CCCCCCCC");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        SimpleScore scoringScheme(2, -1, -2);
        int score = localAlignment(alignmentGraph, scoringScheme);

        SEQAN_ASSERT_EQ(score, 0);

        std::stringstream ss;
        ss << alignmentGraph;

        // Note that for an empty alignment the output of the graph is simply
        // adjacency list and the empty ege list.

        std::stringstream expectedSS;
        expectedSS << "Adjacency list:\n"
                   << "Edge list:\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_local_linear)
{
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        Dna5String strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        SimpleScore scoringScheme(2, -1, -2);
        TFragmentString fragments;
        int score = localAlignment(fragments, strings, scoringScheme);

        SEQAN_ASSERT_EQ(score, 12);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 8, 1, 7, 4));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 4, 1, 4, 3));
    }

    {
        Dna5String strH("GGGGGGGGGG");
        Dna5String strV("CCCCCCCCCC");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        SimpleScore scoringScheme(2, -1, -2);
        TFragmentString fragments;
        int score = localAlignment(fragments, strings, scoringScheme);

        SEQAN_ASSERT_EQ(score, 0);

        SEQAN_ASSERT_EQ(empty(fragments), true);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_affine)
{
    using namespace seqan;

    {
        Dna5String strH("CACACTTAACTTCACAA");
        Dna5String strV("GGGGCTTGAGAGCTTGGGG");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -1, -3);
        int score = localAlignment(align, scoringScheme);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(score, 8);

        SEQAN_ASSERT_EQ(ssH.str(), "CTT---AACTT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTTGAGAGCTT");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_local_affine)
{
    using namespace seqan;

    {
        DnaString strH("CACACTTAACTTCACAA");
        Dna5String strV("GGGGCTTGAGAGCTTGGGG");

        Gaps<DnaString> gapsH(strH);
        Gaps<Dna5String> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -1, -3);
        int score = localAlignment(gapsH, gapsV, scoringScheme);

        SEQAN_ASSERT_EQ(score, 8);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "CTT---AACTT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTTGAGAGCTT");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_local_affine)
{
    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    {
        Dna5String strH("CACACTTAACTTCACAA");
        Dna5String strV("GGGGCTTGAGAGCTTGGGG");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);
        TAlignmentGraph alignmentGraph(strings);

        SimpleScore scoringScheme(2, -1, -1, -3);
        int score = localAlignment(alignmentGraph, scoringScheme);

        SEQAN_ASSERT_EQ(score, 8);

        std::stringstream ss;
        ss << alignmentGraph;

        // Note that the non-overlapping of the segments not part of the
        // sequence is intended: They simply do not take part in the local
        // alignment.
        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :    .    \n"
                   << "        CACA----CTT---AACTTCACAA----\n"
                   << "                |||   | |||         \n"
                   << "        ----GGGGCTTGAGAGCTT-----GGGG\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_local_affine)
{
    using namespace seqan;

    using namespace seqan;

    typedef StringSet<Dna5String, Dependent<> > TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        Dna5String strH("CACACTTAACTTCACAA");
        Dna5String strV("GGGGCTTGAGAGCTTGGGG");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        SimpleScore scoringScheme(2, -1, -1, -3);
        TFragmentString fragments;
        int score = localAlignment(fragments, strings, scoringScheme);

        SEQAN_ASSERT_EQ(score, 8);

        SEQAN_ASSERT(fragments[0] == TFragment(0, 7, 1, 10, 5));
        SEQAN_ASSERT(fragments[1] == TFragment(0, 4, 1, 4, 3));
    }
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

    // Test scoring matrix
    {

        std::stringstream ssH, ssV;

        String<AminoAcid> strH("IGYELAPIPHTRTMDDFGNWWWKKWIHDDELNYFGTQLLIWHLQEKEGEQ");
        String<AminoAcid> strV("KHSDQGQIALLIHNTLQDWRPKVECDSPRTMIRRDFDDPQLAPPPHTNHRGNM");

        Gaps<String<AminoAcid>, ArrayGaps> gapsH(strH);
        Gaps<String<AminoAcid>, ArrayGaps> gapsV(strV);

        Blosum62 scoringScheme;
        int cutoff =  40;

        LocalAlignmentEnumerator<Blosum62, Unbanded> enumerator(scoringScheme, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 69);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "GYELAP--IPHTRTMDDFGNWWWK-KWIH-DD-E-L---NYFGT-QLLIW---HLQEKEG");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "G-QIA-LLI-HN-TLQD---W--RPK-VECDSPRTMIRRD-FDDPQLA--PPPHTNHR-G");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 57);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "KWIHDDELNYFGTQ--LLIWH--LQE---K-E");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT(ssV.str() == "K--HSDQ----G-QIALLI-HNTLQDWRPKVE");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 51);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "IGYE-LA---P-I----PHTRTMD---DFGNWWWKKWIHDD-EL------NYF-GTQL");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "I-HNTLQDWRPKVECDSP--RTM-IRRDF----------DDPQLAPPPHTNH-RGN-M");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 46);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "I--GYE---LAPIP-HT--RTMDDFGN");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "IRRDFDDPQLAP-PPHTNHR-----GN");

        SEQAN_ASSERT_NOT(nextLocalAlignment(gapsH, gapsV, enumerator));
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_enumeration_fragment)
{
    using namespace seqan;

    // TODO(holtgrew): Test after this is written.
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_LOCAL_H_
