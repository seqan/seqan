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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_ALIGN_TEST_MOCK_H_
#define TESTS_ALIGN_TEST_MOCK_H_

namespace impl
{
namespace test_align_mock
{

struct TestAlignSimdVariableLength_;
using VariableLengthSimd = seqan::Tag<TestAlignSimdVariableLength_>;

struct TestAlignSimdEqualLength_;
using EqualLengthSimd = seqan::Tag<TestAlignSimdEqualLength_>;

template <typename TAlphabet, typename TSimdLength>
struct TestSequences_;

template <>
struct TestSequences_<seqan::Dna, EqualLengthSimd>
{
    using TSeq = seqan::String<seqan::Dna>;

    static auto
    getSequences()
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
        appendValue(set, "AACAAGGGACCAGGAGCTCTTCTCC");
        appendValue(set, "ACACGCTAATATAGCGAATCACCGA");
        appendValue(set, "GAACCCGGCGCCACGCAATGGAACG");
        appendValue(set, "TCCTTAACTCCGGCAGGCAATTAAA");
        appendValue(set, "ACAGAAAAATAGGCGAATGAATCTT");
        appendValue(set, "GGGAACGTATGTATAACGCAAAAAA");
        appendValue(set, "TTCTCTGTGTATCGAAGAATGGCCT");

        appendValue(set, "CCGAAGTTTCGATGGACTGGTGCCA");
        appendValue(set, "ACGCGCAGGCATAGTTTTAGGAGAA");
        appendValue(set, "TTATTCGGGGGCAGTGACAACCAAC");

        seqan::StringSet<TSeq>  set2(set);
        std::sort(seqan::begin(set2, seqan::Standard()), seqan::end(set2, seqan::Standard()),
                  [](auto& strA, auto& strB){ return seqan::isLess(strA, strB); });
        return std::make_tuple(set, set2);
    }
};

template <>
struct TestSequences_<seqan::Dna, VariableLengthSimd>
{
    using TSeq = seqan::String<seqan::Dna>;

    static auto
    getSequences()
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
        std::sort(seqan::begin(set2, seqan::Standard()), seqan::end(set2, seqan::Standard()),
                  [](auto& strA, auto& strB){ return seqan::isLess(strA, strB); });
        return std::make_tuple(set, set2);
    }
};

template <>
struct TestSequences_<seqan::AminoAcid, EqualLengthSimd>
{
    using TSeq = seqan::String<seqan::AminoAcid>;

    static auto
    getSequences()
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
        std::sort(seqan::begin(set2, seqan::Standard()), seqan::end(set2, seqan::Standard()),
                  [](auto& strA, auto& strB){ return seqan::isLess(strA, strB); });
        return std::make_tuple(set, set2);
    }
};

template <>
struct TestSequences_<seqan::AminoAcid, VariableLengthSimd>
{
    using TSeq = seqan::String<seqan::AminoAcid>;

    static auto
    getSequences()    {
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
        std::sort(seqan::begin(set2, seqan::Standard()), seqan::end(set2, seqan::Standard()),
                  [](auto& strA, auto& strB){ return seqan::isLess(strA, strB); });
        return std::make_tuple(set, set2);
    }
};
} // namespace test_align_mock
} // namespace impl
#endif // TEST_MOCK_H_
