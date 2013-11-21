// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_WAVELT_TREE_BETA_H_
#define TESTS_WAVELT_TREE_BETA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

template <typename TAlphabetSpecPair>
class RankDictionaryTest : public seqan::Test
{
public:
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 0>::Type TRankDictionarySpec;
};

// TODO(singer): MMap causes lots of seg faults
typedef seqan::TagList<
            seqan::TagList<seqan::WaveletTree<seqan::Dna> >, seqan::TagList<
            seqan::TagList<seqan::WaveletTree<seqan::Dna5> >, seqan::TagList<
            seqan::TagList<seqan::WaveletTree<seqan::AminoAcid> >, seqan::TagList<
            seqan::TagList<seqan::WaveletTree<char> >, seqan::TagList<
            seqan::TagList<seqan::WaveletTree<unsigned char> >, seqan::TagList<
            seqan::TagList<seqan::WaveletTree<signed char> >
            > > > > > >
        RankDictionaryTestTypes;


template <typename T>
class RankDictionaryTestCommon : public RankDictionaryTest<T>
{};

SEQAN_TYPED_TEST_CASE(RankDictionaryTestCommon, RankDictionaryTestTypes);

using namespace seqan;

template <typename TRankDictionary>
void dictionaryClear(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
	TRankDictionary dictionary(text);

    clear(dictionary);

	SEQAN_ASSERT_EQ(empty(dictionary), true);
}

SEQAN_TYPED_TEST(RankDictionaryTestCommon, Clear)
{
    using namespace seqan;

    RankDictionary<typename TestFixture::TRankDictionarySpec> seg;
    dictionaryClear(seg);
}

template <typename TRankDictionary>
void rankDictionaryConstructor(TRankDictionary & /*tag*/)
{
	{
		TRankDictionary rankDictionary;
	}
	{
		String<typename Value<TRankDictionary>::Type> text;        
        generateText(text);
		TRankDictionary rankDictionary(text);
		TRankDictionary rankDictionary2(rankDictionary);
        
//        SEQAN_ASSERT(rankDictionary == rankDictionary2);

		for (unsigned i = 0; i < length(text); ++i)
        {
		    SEQAN_ASSERT_EQ(getValue(rankDictionary, i), text[i]);
		    SEQAN_ASSERT_EQ(getValue(rankDictionary2, i), text[i]);
        }
    }
}

SEQAN_TYPED_TEST(RankDictionaryTestCommon, Constuctor)
{
    using namespace seqan;

    RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
    rankDictionaryConstructor(dictionay);
}

// template <typename TRankDictionary>
// void rankDictionaryDollarPosition(TRankDictionary & /*tag*/)
// {
// 	typedef typename Fibre<TRankDictionary, FibreTreeStructure>::Type TRankDictionaryStructure;
// 	typedef typename Fibre<TRankDictionaryStructure, FibreTreeStructureEncoding>::Type TRankDictionaryVertices;
// 	typedef typename Value<TRankDictionaryVertices>::Type TRankDictionaryVertex;
// 	typedef typename Value<TRankDictionary>::Type TChar;
// 	typedef typename Value<TRankDictionaryVertex, 2>::Type TPos;
// 	typedef typename Fibre<TRankDictionary, FibreDollarPosition>::Type TDollarPos;
// 
// 	TRankDictionary rankDictionary;
// 
// 	TDollarPos dollarPos = 0u;
// 
// 	setDollarPosition(rankDictionary, dollarPos);
// 	SEQAN_ASSERT_EQ(getDollarPosition(rankDictionary), dollarPos);
// }
// 
// template <typename TRankDictionary>
// void rankDictionaryDollarSubstitute(TRankDictionary & /*tag*/)
// {
// 	typedef typename Fibre<TRankDictionary, FibreTreeStructure>::Type TRankDictionaryStructure;
// 	typedef typename Fibre<TRankDictionaryStructure, FibreTreeStructureEncoding>::Type TRankDictionaryVertices;
// 	typedef typename Value<TRankDictionaryVertices>::Type TRankDictionaryVertex;
// 	typedef typename Value<TRankDictionary>::Type TChar;
// 	typedef typename Value<TRankDictionaryVertex, 2>::Type TPos;
// 
// 	TRankDictionary rankDictionary;
// 
// 	setDollarSubstitute(rankDictionary, 'C');
// 	SEQAN_ASSERT_EQ(getDollarSubstitute(rankDictionary), 'C');
// }

template <typename TRankDictionary>
void rankDictionaryEmpty(TRankDictionary & /*tag*/)
{
    {
        TRankDictionary rankDictionary;
        
        SEQAN_ASSERT_EQ(empty(rankDictionary), true);
    }
    {
        String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
        TRankDictionary rankDictionary(text);
        
        SEQAN_ASSERT_EQ(empty(rankDictionary), false);
    }
}

SEQAN_TYPED_TEST(RankDictionaryTestCommon, Empty)
{
    using namespace seqan;

    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
        rankDictionaryEmpty(dictionay);
    }
    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> const dictionay;
        rankDictionaryEmpty(dictionay);
    }
}

template <typename TRankDictionary>
void rankDictionaryGetValue(TRankDictionary & /*tag*/)
{
    String<typename Value<TRankDictionary>::Type> text;
    generateText(text);

    TRankDictionary rankDictionary(text); 

    for (unsigned i = 0; i < length(text); ++i)
    {
        SEQAN_ASSERT_EQ(getValue(rankDictionary, i), text[i]);
    }
}

SEQAN_TYPED_TEST(RankDictionaryTestCommon, GetValue)
{
    using namespace seqan;

    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
        rankDictionaryGetValue(dictionay);
    }
    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> const dictionay;
        rankDictionaryGetValue(dictionay);
    }
}

// NOTE(esiragusa): GetFibre test makes no sense.
//template <typename TRankDictionary>
//void rankDictionaryGetFibre(TRankDictionary & /*tag*/)
//{
//	String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
//	TRankDictionary rankDictionary(text);
//
//    typename Fibre<TRankDictionary, FibreBitStrings>::Type & tempBitStrings = getFibre(rankDictionary, FibreBitStrings());
//    typename Fibre<TRankDictionary, FibreTreeStructure>::Type & tempWaveletTreeStructure = getFibre(rankDictionary, FibreTreeStructure());
//
//    resize(tempBitStrings, 110);
//    _resize(tempWaveletTreeStructure, 100, Exact());
//
//	SEQAN_ASSERT_EQ(length(getFibre(rankDictionary, FibreBitStrings())), 110u);
//	SEQAN_ASSERT_EQ(_length(getFibre(rankDictionary, FibreTreeStructure())), 100u);
//}
//
//SEQAN_TYPED_TEST(RankDictionaryTestCommon, GetFibre)
//{
//    using namespace seqan;
//
//    {
//        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
//        rankDictionaryGetFibre(dictionay);
//    }
//}

template <typename TRankDictionary>
void rankDictionaryGetRank(TRankDictionary & /*tag*/)
{
    typedef typename Value<TRankDictionary>::Type TChar;
    String<TChar> text;
    generateText(text);
    resize(text, 1000);

    TRankDictionary rankDictionary(text);

    for (int i = MinValue<TChar>::VALUE; i <= MaxValue<TChar>::VALUE; ++i)
    {
        unsigned counter = 0;
        for (unsigned j = 0; j < length(text); ++j)
        {
            if (text[j] == (TChar)i)
                ++counter;
            SEQAN_ASSERT_EQ(getRank(rankDictionary, j, (TChar)i), counter);
        }
    }
}

SEQAN_TYPED_TEST(RankDictionaryTestCommon, GetRank)
{
    using namespace seqan;
    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
        rankDictionaryGetRank(dictionay);
    }
    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> const dictionay;
        rankDictionaryGetRank(dictionay);
    }
}

// NOTE(esiragusa): Fill test makes no sense.
//template <typename TRankDictionary>
//void _rankDictionaryFill(TRankDictionary & /*tag*/)
//{
//    String<typename Value<TRankDictionary>::Type> text;
//    generateText(text);
//    resize(text, 1000);
//
//    TRankDictionary rankDictionary(text);
//
//    clear(getFibre(rankDictionary, FibreBitStrings()));
//
//    _fillWaveletTree(rankDictionary, text);
//
//    for (unsigned i = 0; i < length(text); ++i)
//    {
//        SEQAN_ASSERT_EQ(getValue(rankDictionary, i), text[i]);
//    }
//}
//
//SEQAN_TYPED_TEST(RankDictionaryTestCommon, Fill)
//{
//    using namespace seqan;
//    {
//        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
//        _rankDictionaryFill(dictionay);
//    }
//}

// NOTE(esiragusa): OpenSave test shouldn't use operator==().
//template <typename TRankDictionary>
//void rankDictionaryOpenSave(TRankDictionary & /*tag*/)
//{
//	String<typename Value<TRankDictionary>::Type> text;
//    generateText(text);
//    resize(text, 1000);
//
//    CharString tempFilename = SEQAN_TEMP_FILENAME();
//
//    TRankDictionary rankDictionary(text);
//    save(rankDictionary, toCString(tempFilename));
//
//    TRankDictionary rankDictionaryOpen;
//    open(rankDictionaryOpen, toCString(tempFilename));
//    SEQAN_ASSERT(rankDictionary == rankDictionaryOpen);
//}
//
//SEQAN_TYPED_TEST(RankDictionaryTestCommon, OpenSave)
//{
//    using namespace seqan;
//    {
//        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
//        rankDictionaryOpenSave(dictionay);
//    }
//}


#endif  // TESTS_WAVELT_TREE_STRUCTURE_BETA_H_
