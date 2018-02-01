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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Tests for basic/translation.h
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_TRANSLATION_H_
#define SEQAN_TESTS_BASIC_TEST_TRANSLATION_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
// #include <seqan/parallel.h>

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
        translate(outStr, str, SINGLE_FRAME, GeneticCode<VERT_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<YEAST_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<INVERT_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<CILIATE>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<FLATWORM_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<EUPLOTID>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<PROKARYOTE>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<ALT_YEAST>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<ASCIDIAN_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<BLEPHARISMA>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<CHLOROPHYCEAN_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<TREMATODE_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<SCENEDESMUS_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<THRAUSTOCHYTRIUM_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<PTEROBRANCHIA_MITOCHONDRIAL>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GeneticCode<GRACILIBACTERIA>());
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF");
    }
}

template <typename TString>
inline void
test_translation_onestring_singleframe_impl_runtime(TString const & str)
{
    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, CANONICAL);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, VERT_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, YEAST_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, INVERT_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, CILIATE);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, FLATWORM_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, EUPLOTID);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, PROKARYOTE);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, ALT_YEAST);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, ASCIDIAN_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, BLEPHARISMA);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, CHLOROPHYCEAN_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, TREMATODE_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, SCENEDESMUS_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, THRAUSTOCHYTRIUM_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, PTEROBRANCHIA_MITOCHONDRIAL);
        SEQAN_ASSERT_EQ(outStr,
        "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF");
    }

    {
        String<AminoAcid> outStr;
        translate(outStr, str, SINGLE_FRAME, GRACILIBACTERIA);
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

SEQAN_DEFINE_TEST(test_translation_onestring_singleframe_allcodes_runtime)
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
        test_translation_onestring_singleframe_impl_runtime(str);
    }

    // dna5
    {
        Dna5String str(dna);
        test_translation_onestring_singleframe_impl_runtime(str);
    }

    // rna
    {
        RnaString str(rna);
        test_translation_onestring_singleframe_impl_runtime(str);
    }

    // rna5
    {
        Rna5String str(rna);
        test_translation_onestring_singleframe_impl_runtime(str);
    }

    // CharString
    {
        CharString str(dna);
        test_translation_onestring_singleframe_impl_runtime(str);
    }

    // 0-terminated
    {
        test_translation_onestring_singleframe_impl_runtime(dna);
    }

    // std::string
    {
        std::string str(dna);
        test_translation_onestring_singleframe_impl_runtime(str);
    }
}


template <typename TTargetSet, typename TParallelism>
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
        translate(res, str, SINGLE_FRAME, TParallelism());

        SEQAN_ASSERT_EQ(length(res), 1u);
        SEQAN_ASSERT_EQ(res[0], trans[0]);
    }

    {
        translate(res, str, WITH_REVERSE_COMPLEMENT, TParallelism());

        SEQAN_ASSERT_EQ(length(res), 2u);
        SEQAN_ASSERT_EQ(res[0], trans[0]);
        SEQAN_ASSERT_EQ(res[1], trans[3]);
    }

    {
        translate(res, str, WITH_FRAME_SHIFTS, TParallelism());

        SEQAN_ASSERT_EQ(length(res), 3u);
        SEQAN_ASSERT_EQ(res[0], trans[0]);
        SEQAN_ASSERT_EQ(res[1], trans[1]);
        SEQAN_ASSERT_EQ(res[2], trans[2]);
    }

    {
        translate(res, str, SIX_FRAME, TParallelism());

        SEQAN_ASSERT_EQ(length(res), 6u);
        SEQAN_ASSERT_EQ(res[0], trans[0]);
        SEQAN_ASSERT_EQ(res[1], trans[1]);
        SEQAN_ASSERT_EQ(res[2], trans[2]);
        SEQAN_ASSERT_EQ(res[3], trans[3]);
        SEQAN_ASSERT_EQ(res[4], trans[4]);
        SEQAN_ASSERT_EQ(res[5], trans[5]);
    }
}


SEQAN_DEFINE_TEST(test_translation_onestring_multiframe_serial)
{
    typedef StringSet<String<AminoAcid> > TRegSet;
    test_translation_onestring_multiframe_impl<TRegSet, Serial>();
}

SEQAN_DEFINE_TEST(test_translation_onestring_multiframe_concatdirect_serial)
{
    typedef StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > TConcatSet;
    test_translation_onestring_multiframe_impl<TConcatSet, Serial>();
}

SEQAN_DEFINE_TEST(test_translation_onestring_multiframe_parallel)
{
    typedef StringSet<String<AminoAcid> > TRegSet;
    test_translation_onestring_multiframe_impl<TRegSet, Parallel>();
}

SEQAN_DEFINE_TEST(test_translation_onestring_multiframe_concatdirect_parallel)
{
#ifndef __alpha__ // NOTE(h-2): fails on alpha for unknown reasons
    typedef StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > TConcatSet;
    test_translation_onestring_multiframe_impl<TConcatSet, Parallel>();
#endif
}

template <typename TSourceSet, typename TResultSet, typename TParallelism>
inline void
test_translation_stringset_multiframe_impl(TResultSet const & comp,
                                           TSourceSet const & source,
                                           TParallelism const & /**/)
{
    TResultSet res;
    unsigned l = length(source);

    {
        translate(res, source, SINGLE_FRAME, TParallelism());
        unsigned r = length(res) / l;

        SEQAN_ASSERT_EQ(r, 1u);
        for (unsigned i = 0; i < l; ++i)
        {
            SEQAN_ASSERT_EQ(res[0+i*r], comp[0+i*6]);
        }
    }

    {
        translate(res, source, WITH_REVERSE_COMPLEMENT,
                  TParallelism());
        unsigned r = length(res) / l;

        SEQAN_ASSERT_EQ(r, 2u);
        for (unsigned i = 0; i < l; ++i)
        {
            SEQAN_ASSERT_EQ(res[0+i*r], comp[0+i*6]);
            SEQAN_ASSERT_EQ(res[1+i*r], comp[3+i*6]);
        }
    }

    {
        translate(res, source, WITH_FRAME_SHIFTS,
                  TParallelism());
        unsigned r = length(res) / l;

        SEQAN_ASSERT_EQ(r, 3u);
        for (unsigned i = 0; i < l; ++i)
        {
            SEQAN_ASSERT_EQ(length(res), 3*length(source));
            SEQAN_ASSERT_EQ(res[0+i*r], comp[0+i*6]);
            SEQAN_ASSERT_EQ(res[1+i*r], comp[1+i*6]);
            SEQAN_ASSERT_EQ(res[2+i*r], comp[2+i*6]);
        }
    }

    {
        translate(res, source, SIX_FRAME, TParallelism());
        unsigned r = length(res) / l;

        SEQAN_ASSERT_EQ(r, 6u);
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

template <typename TSetSpec, typename TParallelism>
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


        test_translation_stringset_multiframe_impl(comp, source, TParallelism());
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


        test_translation_stringset_multiframe_impl(comp, source, TParallelism());
    }

    // handle empty input sequences correctly, see #2228
    {
        // one empty input results in 6 empty output frames for consistency
        insertValue(comp,  0, "");
        insertValue(comp,  1, "");
        insertValue(comp,  2, "");
        insertValue(comp,  3, "");
        insertValue(comp,  4, "");
        insertValue(comp,  5, "");

        // too short for anything also results in six empty
        insertValue(comp,  6, "");
        insertValue(comp,  7, "");
        insertValue(comp,  8, "");
        insertValue(comp,  9, "");
        insertValue(comp, 10, "");
        insertValue(comp, 11, "");

        // too short for shift results in some empty frames
        insertValue(comp, 12, "T");
        insertValue(comp, 13, "");
        insertValue(comp, 14, "");
        insertValue(comp, 15, "R");
        insertValue(comp, 16, "");
        insertValue(comp, 17, "");

        StringSet<Dna5String, TSetSpec> source;
        appendValue(source, ""); // empty
        appendValue(source, "a"); // empty
        appendValue(source, "acg"); // some empty frames
        appendValue(source, "acgtnncgtaaaccgttaaaccgnntaagtnnaccccggtaccgataan");
        appendValue(source, "ggttacgtatnntaccggttagtacttggggcgagtaganngtt");

        test_translation_stringset_multiframe_impl(comp, source, TParallelism());
    }
}

SEQAN_DEFINE_TEST(test_translation_stringset_multiframe_serial)
{

    test_translation_stringset_multiframe_impl0<Owner<>, Serial >();
}

SEQAN_DEFINE_TEST(test_translation_stringset_multiframe_concatdirect_serial)
{

    test_translation_stringset_multiframe_impl0<Owner<ConcatDirect<> >,Serial>();
}

SEQAN_DEFINE_TEST(test_translation_stringset_multiframe_parallel)
{

    test_translation_stringset_multiframe_impl0<Owner<>, Parallel>();
}

SEQAN_DEFINE_TEST(test_translation_stringset_multiframe_concatdirect_parallel)
{
#ifndef __alpha__ // NOTE(h-2): fails on alpha for unknown reasons
    test_translation_stringset_multiframe_impl0<Owner<ConcatDirect<> >,
                                                Parallel>();
#endif
}

#endif  // SEQAN_TESTS_BASIC_TEST_TRANSLATION_H_
