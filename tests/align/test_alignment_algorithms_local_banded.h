// ==========================================================================
//                  test_alignment_algorithms_local_banded.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_LOCAL_BANDED_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_LOCAL_BANDED_H_

#include <seqan/basic.h>
#include <seqan/score.h>
#include <seqan/align.h>

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_linear_banded)
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

        int score = localAlignment(align, scoringScheme, -2, 2);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTT-AGCT");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_local_linear_banded)
{
    using namespace seqan;

    {
        DnaString strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Gaps<DnaString> gapsH(strH);
        Gaps<Dna5String> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2, -2);

        int score = localAlignment(gapsH, gapsV, scoringScheme, -2, 2);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTT-AGCT");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_local_linear_banded)
{
    using namespace seqan;

    typedef Graph<Alignment<StringSet<DnaString, Dependent<> > > > TGraph;
    typedef StringSet<DnaString> TStringSet;

    {
        DnaString strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        TGraph alignmentGraph(strings);

        SimpleScore scoringScheme(2, -1, -2, -2);

        int score = localAlignment(alignmentGraph, scoringScheme, -2, 2);

        SEQAN_ASSERT_EQ(score, 12);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :    .   \n"
                   << "        GGGG----CTTAAGCTTGGGG------\n"
                   << "                ||| ||||           \n"
                   << "        ----AAAACTT-AGCT-----CTAAAA\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_local_linear_banded)
{
    using namespace seqan;

    typedef StringSet<DnaString> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        DnaString strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        SimpleScore scoringScheme(2, -1, -2, -2);
        TFragmentString fragments;
        int score = localAlignment(fragments, strings, scoringScheme, -2, 2);

        SEQAN_ASSERT_EQ(score, 12);

        SEQAN_ASSERT_EQ(length(fragments), 2u);
        SEQAN_ASSERT_EQ(fragments[0] == TFragment(0, 8, 1, 7, 4), true);
        SEQAN_ASSERT_EQ(fragments[1] == TFragment(0, 4, 1, 4, 3), true);
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_affine_banded)
{
    using namespace seqan;

    {
        Dna5String strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), strH);
        assignSource(row(align, 1), strV);

        SimpleScore scoringScheme(2, -1, -2, -4);

        int score = localAlignment(align, scoringScheme, -2, 2);
        SEQAN_ASSERT_EQ(score, 10);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTT-AGCT");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_gaps_local_affine_banded)
{
    using namespace seqan;

    {
        DnaString strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        Gaps<DnaString> gapsH(strH);
        Gaps<Dna5String> gapsV(strV);

        SimpleScore scoringScheme(2, -1, -2, -4);

        int score = localAlignment(gapsH, gapsV, scoringScheme, -2, 2);

        SEQAN_ASSERT_EQ(score, 10);

        std::stringstream ssH, ssV;
        ssH << gapsH;
        ssV << gapsV;

        SEQAN_ASSERT_EQ(ssH.str(), "CTTAAGCT");
        SEQAN_ASSERT_EQ(ssV.str(), "CTT-AGCT");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_graph_local_affine_banded)
{
    using namespace seqan;

    typedef Graph<Alignment<StringSet<DnaString, Dependent<> > > > TGraph;
    typedef StringSet<DnaString> TStringSet;

    {
        DnaString strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        TGraph alignmentGraph(strings);

        SimpleScore scoringScheme(2, -1, -2, -4);

        int score = localAlignment(alignmentGraph, scoringScheme, -2, 2);

        SEQAN_ASSERT_EQ(score, 10);

        std::stringstream ss;
        ss << alignmentGraph;

        std::stringstream expectedSS;
        expectedSS << "Alignment matrix:\n"
                   << "      0     .    :    .    :    .   \n"
                   << "        GGGG----CTTAAGCTTGGGG------\n"
                   << "                ||| ||||           \n"
                   << "        ----AAAACTT-AGCT-----CTAAAA\n\n\n";

        SEQAN_ASSERT_EQ(ss.str(), expectedSS.str());
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_fragments_local_affine_banded)
{
    using namespace seqan;

    typedef StringSet<DnaString> TStringSet;
    typedef Fragment<unsigned> TFragment;
    typedef String<TFragment>  TFragmentString;

    {
        DnaString strH("GGGGCTTAAGCTTGGGG");
        Dna5String strV("AAAACTTAGCTCTAAAA");

        TStringSet strings;
        appendValue(strings, strH);
        appendValue(strings, strV);

        SimpleScore scoringScheme(2, -1, -2, -4);
        TFragmentString fragments;
        int score = localAlignment(fragments, strings, scoringScheme, -2, 2);

        SEQAN_ASSERT_EQ(score, 10);

        SEQAN_ASSERT_EQ(length(fragments), 2u);
        SEQAN_ASSERT_EQ(fragments[0] == TFragment(0, 8, 1, 7, 4), true);
        SEQAN_ASSERT_EQ(fragments[1] == TFragment(0, 4, 1, 4, 3), true);
    }
}

// ==========================================================================
// Local Alignment Enumeration
// ==========================================================================

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

    // Test scoring matrix
    {
        std::stringstream ssH, ssV;

        String<AminoAcid> strH("IGYELAPIPHTRTMDDFGNWWWKKWIHDDELNYFGTQLLIWHLQEKEGEQ");
        String<AminoAcid> strV("KHSDQGQIALLIHNTLQDWRPKVECDSPRTMIRRDFDDPQLAPPPHTNHRGNM");

        Gaps<String<AminoAcid>, ArrayGaps> gapsH(strH);
        Gaps<String<AminoAcid>, ArrayGaps> gapsV(strV);

        Blosum62 scoringScheme;
        int cutoff =  40;

        LocalAlignmentEnumerator<Blosum62, Banded> enumerator(scoringScheme, -20, 20, cutoff);

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 69);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "GYELAP--IPHTRTMDDFGNWWWK-KWIH-DD-E--L--NYF-GTQLLIW---HLQEKEG");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT_EQ(ssV.str(), "G-QIA-LLI-H-NTLQD-----WRPK-VECDSPRTMIRRD-FDDPQ-LA-PPPHTNHR-G");

        SEQAN_ASSERT(nextLocalAlignment(gapsH, gapsV, enumerator));
        SEQAN_ASSERT_EQ(getScore(enumerator), 51);
        ssH.clear(); ssH.str(""); ssH << gapsH;
        SEQAN_ASSERT_EQ(ssH.str(), "IGYE-L---AP-I----PHTRTM--DDFGNWWWKKWIHDD-EL------NYF-GTQL");
        ssV.clear(); ssV.str(""); ssV << gapsV;
        SEQAN_ASSERT(ssV.str() == "I-HNTLQDWRPKVECDSP--RTMIRRDF----------DDPQLAPPPHTNH-RG-NM");

        SEQAN_ASSERT_NOT(nextLocalAlignment(gapsH, gapsV, enumerator));
    }
}

SEQAN_DEFINE_TEST(test_align_local_alignment_enumeration_banded_fragment)
{
    using namespace seqan;

    // TODO(holtgrew): Test after this is written.
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_LOCAL_BANDED_H_
