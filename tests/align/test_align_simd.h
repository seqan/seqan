// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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

#ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_H_
#define TESTS_ALIGN_TEST_ALIGN_SIMD_H_

#if SEQAN_SIMD_ENABLED

#include <seqan/basic.h>
#include <seqan/align.h>

template <typename TAlphabet>
struct TestSequences_;

template <>
struct TestSequences_<seqan::Dna>
{
    using TSeq = seqan::String<seqan::Dna>

    static std::tuple<seqan::StringSet<TSeq>, seqan::StringSet<TSeq> >
    getSequencesEqualLength()
    {
        seqan::StringSet<TSeq> set;
        appendValue(set, "AGCGACTGCAAACATCAGATCAGAG");
        appendValue(set, "TAATACTAGCATGCGATAAGTCCCT");
        appendValue(set, "GGCACGTGGATGGTTTAGAGGAATC");
        appendValue(set, "AGATTCAAGTCTGGTTAACCATCAA");
        appendValue(set, "ACAGGTCTTGAGTCTAAAATTGTCG");
        appendValue(set, "TCTCCTGCGTACGAGATGGAAATAC");
        appendValue(set, "TAGGTAACTACAGGGACTCCGACGT");
        appendValue(set, "TATGTACGTTGCTCCGTCAGAGGCG");

        appendValue(set, "CCATTCAGGATCACGTTACCGCGAA");
        appendValue(set, "AAAAAGGGACCAGGAGCTCTTCTCC");
        appendValue(set, "CCTGCGGTCACGTCTATAGAAATTA");
        appendValue(set, "CACCATTAACCCTCCTGAGAACCGG");
        appendValue(set, "GAGGCGGGAATCCGTCACGTATGAG");
        appendValue(set, "AAGGTATTTGCCCGATAATCAATAC");
        appendValue(set, "CCCAGGCTTCTAACTTTTTCCACTC");
        appendValue(set, "GCTTGAGCCGGCTAGGCCTTTCTGC");

        appendValue(set, "ATCTCGGGTCCTGCCCAACCGGTCT");
        appendValue(set, "AAAAAGGGACCAGGAGCTCTTCTCC");
        appendValue(set, "ACACGCTAATATAGCGAATCACCGA");
        appendValue(set, "GAACCCGGCGCCACGCAATGGAACG");
        appendValue(set, "TCCTTAACTCCGGCAGGCAATTAAA");
        appendValue(set, "ACAGAAAAATAGGCGAATGAATCTT");
        appendValue(set, "GGGAACGTATGTATAACGCAAAAAA");
        appendValue(set, "TTCTCTGTGTATCGAAGAATGGCCT");

        appendValue(set, "CCGAAGTTTCGATGGACTGGTGCCA");
        appendValue(set, "ACGCGCAGGCATAGTTTTAGGAGAA");
        appendValue(set, "TTATTCGGGGGCAGTGACAACCAAC");

        decltype(set) set2(set);
        std::next_permutation(seqan::begin(set, seqan::Standard()), seqan::end(set, seqan::Standard()),
                              [](auto& strA, auto& strB){ return seqan::isLess(strA, strB); });
        return std::make_tuple(set, set2);
    }

    static std::tuple<seqan::StringSet<TSeq>, seqan::StringSet<TSeq> >
    getSequencesVariableLength()
    {
        seqan::StringSet<TSeq> set;
        appendValue(set, "AGCGACTGCAAACATCAGATCAGAGGTAGAG");
        appendValue(set, "TAATACTAGCATGCGATAAGTCCCT");
        appendValue(set, "GGCACGTGTGGTTTAGAGGAATC");
        appendValue(set, "AGATTCAAGTCTGGTTAACCATCAA");
        appendValue(set, "ACAGGTCTTGAGTCTAAAATTGTCGAA");
        appendValue(set, "TCTCCTGCGTACGAGATGGAAATAC");
        appendValue(set, "TAGGTAACTACAGGGACACGT");
        appendValue(set, "TATGTACGTCTCCGTCAGAGGCG");

        appendValue(set, "CCATTCAGGATCACGTTACCGCGAAGTACCC");
        appendValue(set, "AAGGGACCAGGAGCTCTTCTCC");
        appendValue(set, "CCTGCGGTCACGTCTATAGAAATT");
        appendValue(set, "CACCATTAACCCTCCTGAGAACCGAGTAGG");
        appendValue(set, "GAGGCGGGAATCCGTCACGTATGAG");
        appendValue(set, "AAGGTATTTGCCCGATAATCAATACGATGAGATAGAGAGATAGAATAGAGAAGGGACCGCGCATGACTACGATCGACTGACTACGA");
        appendValue(set, "CGAGTATATCGAGAGAGGTCACG");
        appendValue(set, "GCTTGAGCCGGCTAGGCTCTGC");

        appendValue(set, "ATCTCGGGTCCTGCCAACCGGTCT");
        appendValue(set, "AAAAAGGGACCAGGAGCTCTTCTCC");
        appendValue(set, "ACACGCTAATATAGCGAATCACCGA");
        appendValue(set, "AATGGAACG");
        appendValue(set, "TCCTTAACTCCGGCAGGCAATTATACCGGACTGACACTTAAA");
        appendValue(set, "ACAGAAAAATAGGCGAATGAAACACTCTT");
        appendValue(set, "GGGAACGTATGTATAACGCAAAAA");
        appendValue(set, "TTCTCTGTGTATCGAAGAATGCT");

        appendValue(set, "CCGAAGTTTCGATGGATGGATTCCACACACCTGGTGCCA");
        appendValue(set, "ACGCGCAGGCATAGTTGGAGAA");
        appendValue(set, "TTATTCGGGGGCAGTGACAACACTTAGCGACTAC");

        auto set2(set);
        std::next_permutation(seqan::begin(set, seqan::Standard()), seqan::end(set, seqan::Standard()),
                              [](auto& strA, auto& strB){ return seqan::isLess(strA, strB); });
        return std::make_tuple(set, set2);
    }
};

template <>
struct TestSequences_<seqan::AminoAcid>
{
    using TSeq = seqan::String<seqan::AminoAcid>

    static std::tuple<seqan::StringSet<TSeq>, seqan::StringSet<TSeq> >
    getSequencesEqualLength()
    {
        seqan::StringSet<TSeq> set;
        appendValue(set, "FNQSAEYPDISLHCGVLKWRATLGT");
        appendValue(set, "EIKSDVLLHRPGNIGMQVAESYFAT");
        appendValue(set, "PIIMWSMKNRTIERLPTGVLMISHT");
        appendValue(set, "FMATNEKVHCACGADYQMIIDCNEA");
        appendValue(set, "MFHQTSNANWMFVSNKFHIKFGTLD");
        appendValue(set, "SNEMGQCFPHEPACFFDKDFRLFIN");
        appendValue(set, "FPWAHYVVHTLREHRKDANHRSTSY");
        appendValue(set, "QYRNTESMGCEMRCFTETIMIAGVA");

        appendValue(set, "VVRMDGKEVLKQHVPTYADKHPTGQ");
        appendValue(set, "TMLKWCEWCFAEFPPFASEPKFPPN");
        appendValue(set, "GTWGWVDGVHHTMGEQCGPGRACWG");
        appendValue(set, "ECDFQTWYFYCVNQEIFELFICCMG");
        appendValue(set, "KRRELNGQERGGWWTVDGPGVSMGT");
        appendValue(set, "CWAAHYVCWRTKQKQLVAFQRLNCI");
        appendValue(set, "NRLVGFQIHCFLIRCVEPGQTHTID");
        appendValue(set, "AYYVRGFMMGQMYGRPVILMTFTKP");

        appendValue(set, "SFTQPVELHIPHYWWHLAYFMIMFY");
        appendValue(set, "PMNKMFDFNNHQDLLTFTKRFPTPW");
        appendValue(set, "VIPMIYHDWSIISALMMQKDIYYIA");
        appendValue(set, "TPGMWGMATLTGNFNSIFVSKYVKN");
        appendValue(set, "GKELWGMVIARAGMAVQNMYSRDTF");
        appendValue(set, "VHASDLYAKCYSNCVYQENIDIAEV");
        appendValue(set, "KQSGTLSGPQYWENVHRVLEDYPKE");
        appendValue(set, "DPHGYCFYEGTFAWDVEVHEFNNKD");

        appendValue(set, "NMQDVIGGKSLAQHSSVTYKAQQEH");
        appendValue(set, "CQTPRWECSLNFDEKEAADLMIDVS");
        appendValue(set, "PMMDLDHCMLIECLRPHNRDNCARH");

        decltype(set) set2(set);
        std::next_permutation(seqan::begin(set, seqan::Standard()), seqan::end(set, seqan::Standard()),
                              [](auto& strA, auto& strB){ return seqan::isLess(strA, strB); });
        return std::make_tuple(set, set2);
    }

    static std::tuple<seqan::StringSet<TSeq>, seqan::StringSet<TSeq> >
    getSequencesVariableLength()
    {
        seqan::StringSet<TSeq> set;
        appendValue(set, "FNQSAEYPDISHCGVMQLKWRATLGT");
        appendValue(set, "EIKSDVLLHRWSMKNPGNILMIDVGMQVAESYFAT");
        appendValue(set, "PIIMWSMKNRTIEPTGLMISHT");
        appendValue(set, "FMATNEKVHCACGWSMKNADLMIDVYQMIIDCNEA");
        appendValue(set, "MWSMKNFHQTSNANWMFVSNMQKFHIKFGTLD");
        appendValue(set, "SNEGQCFPHEPACFWSMKFDKDFRLFIN");
        appendValue(set, "FPWAHYVVHTLREHMQRKDANHRSTSY");
        appendValue(set, "QYRNTWSMKNESMGCEMRFLMIVTETIMIAGVA");

        appendValue(set, "VVRMDGKEVLWSMKNKQHVPTYADKHPTGQ");
        appendValue(set, "TMLKWCEWCFALMIDVEFPPFASEPKFPPN");
        appendValue(set, "GTWGVDGVHHWSMWSMKNLMIDVTMGEQCGPGRACWG");
        appendValue(set, "ECDFQTWYFYCVNQMQEIFELFICCMG");
        appendValue(set, "KRREWSMKNLNGQERGGWWTVDGPGVSMGT");
        appendValue(set, "CWAAHYCWRWWSMKNSMKNTKLMIDVQMQKQLVAFQRLNCI");
        appendValue(set, "NRLVGFQIHCFIRCVEPGQTHTID");
        appendValue(set, "AYYVRGFMMGQMMQYGRPVILMTFTKP");

        appendValue(set, "SFTQPVELHIPHYWLMIDVWHLAYFMIMFY");
        appendValue(set, "PMNKMFDFNHQMQDLLTFTKPTPW");
        appendValue(set, "VIPMIYHDWSIISALMMLMIDVQKDIYYIA");
        appendValue(set, "TPGMWGMATLTGMQNFNSFVSKYVKN");
        appendValue(set, "GKELWGMVIARAGMAVQNLMIDVMYSRDTF");
        appendValue(set, "VHASDLWSNYAKCYSNCVYQEIDIAEV");
        appendValue(set, "KQSGTLSMQGPYWENVHRVLLMIDVEDYPKE");
        appendValue(set, "DPHGYCFMQYEGTFAWDVEVHEFNNKD");

        appendValue(set, "NMQDVIGGKSLAQHSSVTYAQQEH");
        appendValue(set, "CQTPRWECMQSLNFDEKEAADLMIDVS");
        appendValue(set, "PMMDLDWSMKNMLIECLRPHNRMQDNLMIDVCARH");
        
        auto set2(set);
        std::next_permutation(seqan::begin(set, seqan::Standard()), seqan::end(set, seqan::Standard()),
                              [](auto& strA, auto& strB){ return seqan::isLess(strA, strB); });
        return std::make_tuple(set, set2);
    }
};

// These are interface tests.
struct LocalAlignTester_
{
    template <typename TAlign, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static seqan::String<TScoreValue>
    run(TAlign & align, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const &,
        int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return localAlignment(align, score);
        else
            return localAlignment(align, score, lDiag, uDiag);
    }

    template <typename TAlign, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static TScoreValue
    gold(TAlign & align, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const &,
         int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return localAlignment(align, score);
        else
            return localAlignment(align, score, lDiag, uDiag);

    }
};

struct LocalAlignScoreTester_
{
    template <typename TStringsH, typename TStringsV, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static seqan::String<TScoreValue>
    run(TStringsH const & strH, TStringsV const & strV, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const &)
    {
        return localAlignmentScore(strH, strV, score);
    }

    template <typename TSeqH, typename TSeqV, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static TScoreValue
    gold(TSeqH const & seqH, TSeqV const & seqV, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const &)
    {
        return localAlignmentScore(seqH, seqV, score);
    }
};

struct GlobalAlignTester_
{
    template <typename TAlign, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static seqan::String<TScoreValue>
    run(TAlign & align, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const & config,
        int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return globalAlignment(align, score, config);
        else
            return globalAlignment(align, score, config, lDiag, uDiag);
    }

    template <typename TAlign, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static TScoreValue
    gold(TAlign & align, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const & config,
         int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return globalAlignment(align, score, config);
        else
            return globalAlignment(align, score, config, lDiag, uDiag);
    }
};

struct GlobalAlignScoreTester_
{
    template <typename TStringsH, typename TStringsV, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static seqan::String<TScoreValue>
    run(TStringsH const & strH, TStringsV const & strV, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const & config,
        int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return globalAlignmentScore(strH, strV, score, config);
        else
            return globalAlignmentScore(strH, strV, score, config, lDiag, uDiag);
    }

    template <typename TSeqH, typename TSeqV, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static TScoreValue
    gold(TSeqH const & seqH, TSeqV const & seqV, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const & config,
         int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return globalAlignmentScore(seqH, seqV, score, config);
        else
            return globalAlignmentScore(seqH, seqV, score, config, lDiag, uDiag);
    }
};

template <typename TAlphabet, typename TScoreValue, typename TScoreSpec, typename TAlignConfig, typename TFunctor>
void testAlignSimd(TFunctor const &,
                   seqan::Score<TScoreValue, TScoreSpec> const & score,
                   TAlignConfig const & config,
                   int const lDiag = seqan::MinValue<int>::VALUE,
                   int const uDiag = seqan::MaxValue<int>::VALUE)
{
    using namespace seqan;

    for (unsigned i = 0; i < 2; ++i)
    {
        typedef typename
        StringSet<Align<String<TAlphabet> > > alignments;
        for(unsigned i = 0; i < 34; ++i)
        {
            Align<DnaString> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), seqH);
            assignSource(row(align, 1), seqV);
            appendValue(alignments, align);
        }

        String<TScoreValue> scores = TFunctor::run(alignments, score, config, lDiag, uDiag);

        SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));

        Align<DnaString> goldAlign;
        resize(rows(goldAlign), 2);
        assignSource(row(goldAlign, 0), seqH);
        assignSource(row(goldAlign, 1), seqV);

        TScoreValue goldScore = TFunctor::gold(goldAlign, score, config, lDiag, uDiag);

        for(size_t i = 0; i < 34; ++i)
        {
            SEQAN_ASSERT_EQ(scores[i], goldScore);
            SEQAN_ASSERT(row(alignments[i], 0) == row(goldAlign, 0));
            SEQAN_ASSERT(row(alignments[i], 1) == row(goldAlign, 1));
        }
    }
}

template <typename TSeq, typename TScoreValue, typename TScoreSpec, typename TAlignConfig, typename TTester>
void testAlignSimdScore(TTester const &,
                        TSeq const & seqH,
                        TSeq const & seqV,
                        seqan::Score<TScoreValue, TScoreSpec> const & score,
                        TAlignConfig const & config,
                        int const lDiag = seqan::MinValue<int>::VALUE,
                        int const uDiag = seqan::MaxValue<int>::VALUE)
{
    using namespace seqan;

    StringSet<DnaString> stringsH, stringsV;

    for(unsigned i = 0; i < 34; ++i)
    {
        appendValue(stringsH, seqH);
        appendValue(stringsV, seqV);
    }

    String<TScoreValue> scores = TTester::run(stringsH, stringsV, score, config, lDiag, uDiag);

    SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));


    TScoreValue goldScore = TTester::gold(seqH, seqV, score, config, lDiag, uDiag);

    for(size_t i = 0; i < 34; ++i)
        SEQAN_ASSERT_EQ(scores[i], goldScore);
}

// Problem is clearly that the result is different.
SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_global_linear)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    testAlignSimd<Dna>(GlobalAlignTester_(), Score<int, Simple>(2, -1, -1), alignConfig);
    testAlignSimd<AminoAcid>(GlobalAlignTester_(), Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_global_linear)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    auto sets = getDnaSequencesEqualLength();
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), std::get<0>(sets), std::get<1>(sets), Score<int, Simple>(2, -1, -1), alignConfig);
    sets = getDnaSequencesVariableLength();
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), std::get<0>(sets), std::get<1>(sets), Score<int, Simple>(2, -1, -1), alignConfig);

    auto setsAA = getAninoAcidSequencesEqualLength();
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), std::get<0>(setsAA), std::get<1>(setsAA), Score<int, Simple>(2, -1, -1), alignConfig);
    setsAA = getAninoAcidSequencesVariableLength();
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), std::get<0>(setsAA), std::get<1>(setsAA), Score<int, Simple>(2, -1, -1), alignConfig);

}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_global_affine)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_global_affine)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_overlap_linear)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with no overlap - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "GGGGGGGGG", "CCCCCCCCC", Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_overlap_linear)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with no overlap - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "GGGGGGGGG", "CCCCCCCCC", Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_overlap_affine)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_overlap_affine)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_semi_global_linear)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_semi_global_linear)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_semi_global_affine)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_semi_global_affine)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore<DnaString>(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

// ----------------------------------------------------------------------------
// Global Alignments Banded.
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_global_linear_banded)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    Score<int, Simple> score(2, -1, -1);
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ATGT", "ATAGAT", score, alignConfig, -3, 2);
    // Alignment with both leading and trailing gaps in one row - Simd Version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAGGGGTTTT", "GGGGG", score, alignConfig, -2, 8);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTTTT", "TTTTTTTTGGGGGGGG", score, alignConfig, -4, 4);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_global_affine_banded)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1, -3);
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ATGT", "ATAGAT", scoringScheme, alignConfig, -3, 2);
    // Alignment with both leading and trailing gaps in one row - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAGGGGTTTT", "GGGGG", scoringScheme, alignConfig, -2, 8);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTTTT", "TTTTTTTTGGGGGGGG", scoringScheme, alignConfig, -4, 4);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_overlap_linear_banded)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1);
    // Simple alignment without any leading or trailing gaps.
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ATGT", "ATAGAT", scoringScheme, alignConfig, -2, 2);
    // Alignment with both leading and trailing gaps in one row - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAGGGGTTTT", "GGGGG", scoringScheme, alignConfig, -2, 2);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTTTT", "TTTTTTTTGGGGGGGG", scoringScheme, alignConfig, -2, 2);
    // Alignment that starts at first position where the band crosses the bottom of the matrix - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AACGCATTTTT", "TTTACGCA", scoringScheme, alignConfig, -2, 2);
    // Alignment that starts at first position where the band crosses the bottom of the matrix - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AACGCA", "TTTACGCA", scoringScheme, alignConfig, -2, 4);
    // Alignment that starts at first position where the band crosses the bottom of the matrix - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ACGAGTGTTTGCC", "TTTTTACGA", scoringScheme, alignConfig, -5, 7);
    // Alignment that starts at first position where the band crosses the bottom of the matrix - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ACGA", "TTTTTACGA", scoringScheme, alignConfig, -5, 4);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_overlap_affine_banded)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1, -3);

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "ATGT", "ATAGAT", scoringScheme, alignConfig, -2, 2);
    // Alignment with both leading and trailing gaps in one row - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", scoringScheme, alignConfig, -2, 2);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", scoringScheme, alignConfig, -2, 2);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_semi_global_linear_banded)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1);
    // More or less simple alignment - Simd version
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTTTTG", "AATTTTTTTTTTGGGGG", scoringScheme, alignConfig, -2, 2);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_semi_global_affine_banded)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1, -4);
    testAlignSimd<DnaString>(GlobalAlignTester_(), "AAAAAATTTTTTTTG", "AATTTTTTTTTTGGGGG", scoringScheme, alignConfig, -2, 2);
}

// ----------------------------------------------------------------------------
// Local Alignments
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_linear)
{
    testAlignSimd<seqan::DnaString>(LocalAlignTester_(), "GGGGCTTAAGCTTGGGG", "AAAACTTAGCTCTAAAA", seqan::SimpleScore(2, -1, -2), seqan::Nothing());
    testAlignSimd<seqan::DnaString>(LocalAlignTester_(), "GGGGGGGGG", "CCCCCCCCC", seqan::SimpleScore(2, -1, -2), seqan::Nothing());
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_affine)
{
    testAlignSimd<seqan::DnaString>(LocalAlignTester_(), "CACACTTAACTTCACAA", "GGGGCTTGAGAGCTTGGGG", seqan::SimpleScore(2, -1, -1, -3), seqan::Nothing());
}

// ----------------------------------------------------------------------------
// Local Alignments Banded
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_linear_banded)
{
    using namespace seqan;

    testAlignSimd<seqan::DnaString>(LocalAlignTester_(), "GGGGCTTAAGCTTGGGG", "AAAACTTAGCTCTAAAA", SimpleScore(2, -1, -2, -2), Nothing(), -2, 2);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_affine_banded)
{
    using namespace seqan;

    testAlignSimd<seqan::DnaString>(LocalAlignTester_(), "GGGGCTTAAGCTTGGGG", "AAAACTTAGCTCTAAAA", SimpleScore(2, -1, -2, -4), Nothing(), -2, 2);
}
#endif  // SEQAN_SIMD_ENABLED

#endif  // #ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_H_
