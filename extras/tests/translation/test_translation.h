// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SEQAN_EXTRAS_TESTS_BASIC_TEST_TRANSLATION_H_
#define SEQAN_EXTRAS_TESTS_BASIC_TEST_TRANSLATION_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/translation.h>

using namespace seqan;

template <typename TString>
inline void
test_translation_onestring_singleframe_impl(TString const & str)
{
    {
        String<AminoAcid> outStr;
        translate(outStr, str);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<VertMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<YeastMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<InvertMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<Ciliate>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<FlatwormMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<Euplotid>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<Prokaryote>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<AltYeast>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<AscidianMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<Blepherisma>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<ChlorophyceanMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<TrematodeMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<ScenedesmusMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<ThraustochytriumMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<PterobranchiaMitochondrial>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<Gracilibacteria>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF");
    }
}


SEQAN_DEFINE_TEST(test_translation_onestring_singleframe_allcodes)
{
    const char * dna = "aaaaacaagaatacaaccacgactagaagcaggagtataatcatgattcaacacc"
    "agcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcggg"
    "ggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgttt";
    const char * rna = "aaaaacaagaauacaaccacgacuagaagcaggaguauaaucaugauucaacacc"
    "agcauccacccccgccucgacgccggcgucuacuccugcuugaagacgaggaugcagccgcggcuggaggcggg"
    "gguguagucgugguuuaauacuaguauucauccucgucuugaugcugguguuuauucuuguuu";


    // dna
    {
        DnaString str(dna);
        test_translation_onestring_singleframe_impl(str);
    }

    // dna5
    {
        Dna5String str(dna);
        test_translation_onestring_singleframe_impl(str);
    }

    // rna
    {
        RnaString str(rna);
        test_translation_onestring_singleframe_impl(str);
    }

    // rna5
    {
        Rna5String str(rna);
        test_translation_onestring_singleframe_impl(str);
    }

    // CharString
    {
        CharString str(dna);
        test_translation_onestring_singleframe_impl(str);
    }

    // 0-terminated
    {
        test_translation_onestring_singleframe_impl(dna);
    }

    // std::string
    {
        std::string str(dna);
        test_translation_onestring_singleframe_impl(str);
    }
}

template <typename TTargetSet>
inline void
test_translation_onestring_multiframe_impl()
{
    Dna5String str("acgtnncgtaaaccgttaaaccgnntaagtnnaccccggtaccgataan");
    //              __|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
    //               __|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
    //                __|__|__|__|__|__|__|__|__|__|__|__|__|__|__|

    StringSet<String<AminoAcid> > trans;
    appendValue(trans, "TXRKPLNXXSXPRYR*");
    appendValue(trans, "RXVNR*TX*XXPGTDX");
    appendValue(trans, "XX*TVKPXKXTPVPI" );
    appendValue(trans, "XIGTGXXLXGLTVYXT");
    appendValue(trans, "LSVPGXTXXV*RFTXR");
    appendValue(trans, "YRYRGXLXRFNGLXX" );

    TTargetSet res;

    {
        translate(res, str, SINGLE_FRAME);

        SEQAN_ASSERT_EQ(length(res), 1);
        SEQAN_ASSERT_EQ(res[0], trans[0]);
    }

    {
        translate(res, str, WITH_REV_COMP);

        SEQAN_ASSERT_EQ(length(res), 2);
        SEQAN_ASSERT_EQ(res[0], trans[0]);
        SEQAN_ASSERT_EQ(res[1], trans[3]);
    }

    {
        translate(res, str, WITH_FRAME_SHIFT);

        SEQAN_ASSERT_EQ(length(res), 3);
        SEQAN_ASSERT_EQ(res[0], trans[0]);
        SEQAN_ASSERT_EQ(res[1], trans[1]);
        SEQAN_ASSERT_EQ(res[2], trans[2]);
    }

    {
        translate(res, str, SIX_FRAME);

        SEQAN_ASSERT_EQ(length(res), 6);
        SEQAN_ASSERT_EQ(res[0], trans[0]);
        SEQAN_ASSERT_EQ(res[1], trans[1]);
        SEQAN_ASSERT_EQ(res[2], trans[2]);
        SEQAN_ASSERT_EQ(res[3], trans[3]);
        SEQAN_ASSERT_EQ(res[4], trans[4]);
        SEQAN_ASSERT_EQ(res[5], trans[5]);
    }
}


SEQAN_DEFINE_TEST(test_translation_onestring_multiframe)
{
    typedef StringSet<String<AminoAcid> > TRegSet;
    test_translation_onestring_multiframe_impl<TRegSet>();
}

SEQAN_DEFINE_TEST(test_translation_onestring_multiframe_concatdirect)
{
    typedef StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > TConcatSet;
    test_translation_onestring_multiframe_impl<TConcatSet>();
}

template <typename TSourceSet, typename TResultSet>
inline void
test_translation_stringset_multiframe_impl(TResultSet const & comp,
                                                      TSourceSet const & source)
{
    TResultSet res;
    unsigned l = length(source);

    {
        translate(res, source, SINGLE_FRAME);
        unsigned r = length(res) / l;

        SEQAN_ASSERT_EQ(r, 1);
        for (unsigned i = 0; i < l; ++i)
        {
            SEQAN_ASSERT_EQ(res[0+i*r], comp[0+i*6]);
        }
    }

    {
        translate(res, source, WITH_REV_COMP);
        unsigned r = length(res) / l;

        SEQAN_ASSERT_EQ(r, 2);
        for (unsigned i = 0; i < l; ++i)
        {
            SEQAN_ASSERT_EQ(res[0+i*r], comp[0+i*6]);
            SEQAN_ASSERT_EQ(res[1+i*r], comp[3+i*6]);
        }
    }

    {
        translate(res, source, WITH_FRAME_SHIFT);
        unsigned r = length(res) / l;

        SEQAN_ASSERT_EQ(r, 3);
        for (unsigned i = 0; i < l; ++i)
        {
            SEQAN_ASSERT_EQ(length(res), 3*length(source));
            SEQAN_ASSERT_EQ(res[0+i*r], comp[0+i*6]);
            SEQAN_ASSERT_EQ(res[1+i*r], comp[1+i*6]);
            SEQAN_ASSERT_EQ(res[2+i*r], comp[2+i*6]);
        }
    }

    {
        translate(res, source, SIX_FRAME);
        unsigned r = length(res) / l;

        SEQAN_ASSERT_EQ(r, 6);
        for (unsigned i = 0; i < l; ++i)
        {

            SEQAN_ASSERT_EQ(res[0+i*r], comp[0+i*6]);
            SEQAN_ASSERT_EQ(res[1+i*r], comp[1+i*6]);
            SEQAN_ASSERT_EQ(res[2+i*r], comp[2+i*6]);
            SEQAN_ASSERT_EQ(res[3+i*r], comp[3+i*6]);
            SEQAN_ASSERT_EQ(res[4+i*r], comp[4+i*6]);
            SEQAN_ASSERT_EQ(res[5+i*r], comp[5+i*6]);
        }
    }
}

template <typename TSetSpec>
inline void
test_translation_stringset_multiframe_impl0()
{
    StringSet<String<AminoAcid>, TSetSpec> comp;
    appendValue(comp, "TXRKPLNXXSXPRYR*");
    appendValue(comp, "RXVNR*TX*XXPGTDX");
    appendValue(comp, "XX*TVKPXKXTPVPI" );
    appendValue(comp, "XIGTGXXLXGLTVYXT");
    appendValue(comp, "LSVPGXTXXV*RFTXR");
    appendValue(comp, "YRYRGXLXRFNGLXX" );

    appendValue(comp, "GYVXYRLVLGASRX");
    appendValue(comp, "VTYXTG*YLGRVXX" );
    appendValue(comp, "LRXXPVSTWGE*XV" );
    appendValue(comp, "NXLLAPSTNRXXT*" );
    appendValue(comp, "XXYSPQVLTGXIRN" );
    appendValue(comp, "XSTRPKY*PVXYVT" );

    //DNA
    {
        StringSet<Dna5String, TSetSpec> source;
        Dna5String str("acgtnncgtaaaccgttaaaccgnntaagtnnaccccggtaccgataan");
        //              __|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
        //               __|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
        //                __|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
        appendValue(source, str);
        str = "ggttacgtatnntaccggttagtacttggggcgagtaganngtt";
        //     __|__|__|__|__|__|__|__|__|__|__|__|__|__|__
        //      __|__|__|__|__|__|__|__|__|__|__|__|__|__|
        //       __|__|__|__|__|__|__|__|__|__|__|__|__|__|
        appendValue(source, str);


        test_translation_stringset_multiframe_impl(comp, source);
    }

    //RNA
    {
        StringSet<Rna5String, TSetSpec> source;
        Rna5String str("acgunncguaaaccguuaaaccgnnuaagunnaccccgguaccgauaan");
        //              __|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
        //               __|__|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
        //                __|__|__|__|__|__|__|__|__|__|__|__|__|__|__|
        appendValue(source, str);
        str = "gguuacguaunnuaccgguuaguacuuggggcgaguagannguu";
        //     __|__|__|__|__|__|__|__|__|__|__|__|__|__|__
        //      __|__|__|__|__|__|__|__|__|__|__|__|__|__|
        //       __|__|__|__|__|__|__|__|__|__|__|__|__|__|
        appendValue(source, str);


        test_translation_stringset_multiframe_impl(comp, source);
    }
}

SEQAN_DEFINE_TEST(test_translation_stringset_multiframe)
{

    test_translation_stringset_multiframe_impl0<Owner<> >();
}

SEQAN_DEFINE_TEST(test_translation_stringset_multiframe_concatdirect)
{

    test_translation_stringset_multiframe_impl0<Owner<ConcatDirect<> > >();
}

#endif  // SEQAN_EXTRAS_TESTS_BASIC_TEST_TRANSLATION_H_
