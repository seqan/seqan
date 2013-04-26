// ==========================================================================
//                 seqan - the library for sequence analysis
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
            seqan::TagList<seqan::WaveletTree<signed char> >, seqan::TagList<
            seqan::TagList<seqan::SequenceBitMask<seqan::Dna> >, seqan::TagList<
            seqan::TagList<seqan::SequenceBitMask<seqan::Dna5> >, seqan::TagList<
            seqan::TagList<seqan::SequenceBitMask<seqan::AminoAcid> >
            > > > > > >
            > > >
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
        
        SEQAN_ASSERT(rankDictionary == rankDictionary2);

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
// 
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
// 
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

template <typename TRankDictionary>
void rankDictionaryGetFibre(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
	TRankDictionary rankDictionary(text);

    typename Fibre<TRankDictionary, FibreBitStrings>::Type & tempBitStrings = getFibre(rankDictionary, FibreBitStrings());
    typename Fibre<TRankDictionary, FibreTreeStructure>::Type & tempWaveletTreeStructure = getFibre(rankDictionary, FibreTreeStructure());

    resize(tempBitStrings, 110);
    _resize(tempWaveletTreeStructure, 100, Exact());

	SEQAN_ASSERT_EQ(length(getFibre(rankDictionary, FibreBitStrings())), 110u);
	SEQAN_ASSERT_EQ(_length(getFibre(rankDictionary, FibreTreeStructure())), 100u);
}

template <typename TValue>
void rankDictionaryGetFibre(RankDictionary<SequenceBitMask<TValue> > & /*tag*/)
{
    String<typename Value<RankDictionary<SequenceBitMask<TValue> > >::Type> text = "ACGTNACGTNACGTN";
	RankDictionary<SequenceBitMask<TValue> > rankDictionary(text);

    typename Fibre<RankDictionary<SequenceBitMask<TValue> >, FibreBitStrings>::Type & tempBitStrings = getFibre(rankDictionary, FibreBitStrings());

    resize(tempBitStrings, 110);

	SEQAN_ASSERT_EQ(length(getFibre(rankDictionary, FibreBitStrings())), 110u);
}


SEQAN_TYPED_TEST(RankDictionaryTestCommon, GetFibre)
{
    using namespace seqan;

    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
        rankDictionaryGetFibre(dictionay);
    }
}

template <typename TRankDictionary>
void rankDictionaryCountOcc(TRankDictionary & /*tag*/)
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
            if(text[j] == (TChar)i)
                ++counter;
            SEQAN_ASSERT_EQ(countOccurrences(rankDictionary, (TChar)i, j), counter);
        }
    }
}

SEQAN_TYPED_TEST(RankDictionaryTestCommon, CountOcc)
{
    using namespace seqan;
    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
        rankDictionaryCountOcc(dictionay);
    }
    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> const dictionay;
        rankDictionaryCountOcc(dictionay);
    }
}


// // template <typename TString, typename TSpec>
// // void rankDictionaryCountOcc(WaveletTree<TString, FmiDollarSubstituted<TSpec> > & /*tag*/)
// // {
// //     typedef WaveletTree<TString, FmiDollarSubstituted<TSpec> > TRankDictionary;
// // 	typedef typename Fibre<TRankDictionary, FibreTreeStructure>::Type TRankDictionaryStructure;
// // 	typedef typename Fibre<TRankDictionaryStructure, FibreTreeStructureEncoding>::Type TRankDictionaryVertices;
// // 	typedef typename Value<TRankDictionaryVertices>::Type TRankDictionaryVertex;
// // 	typedef typename Value<TRankDictionary>::Type TChar;
// // 	typedef typename Value<TRankDictionaryVertex, 2>::Type TPos;
// // 
// // 	{
// // 		String<TChar> text;
// //  		generateText(text);
// //  		resize(text, 1000);
// // 
// // 		TRankDictionary rankDictionary(text);
// // 		setDollarSubstitute(rankDictionary, getCharacter(rankDictionary, 0u));
// // 		setDollarPosition(rankDictionary, 0u);
// // 
// //         for (int i = MinValue<TChar>::VALUE; i <= MaxValue<TChar>::VALUE; ++i)
// //         {
// //             unsigned counter = 0;
// //             for (unsigned j = 0; j < length(text); ++j)
// //             {
// //                 if(text[j] == (TChar)i)
// //                 {
// //                     ++counter;
// //                     if(text[j] == getDollarSubstitute(rankDictionary) && j == getDollarPosition(rankDictionary))
// //                         --counter;
// //                 }
// // 		        SEQAN_ASSERT_EQ(countOccurrences(rankDictionary, (TChar)i, j), counter);
// //             }
// //         }
// //     }
// // }
// 
template <typename TRankDictionary>
void _rankDictionaryFill(TRankDictionary & /*tag*/)
{
    String<typename Value<TRankDictionary>::Type> text;
    generateText(text);
    resize(text, 1000);

    TRankDictionary rankDictionary(text);

    clear(getFibre(rankDictionary, FibreBitStrings()));

    _fillWaveletTree(rankDictionary, text);

    for (unsigned i = 0; i < length(text); ++i)
    {
        SEQAN_ASSERT_EQ(getValue(rankDictionary, i), text[i]);
    }
}

template <typename TValue>
void _rankDictionaryFill(RankDictionary<SequenceBitMask<TValue> > & /*tag*/) {}

SEQAN_TYPED_TEST(RankDictionaryTestCommon, Fill)
{
    using namespace seqan;
    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
        _rankDictionaryFill(dictionay);
    }
}

template <typename TRankDictionary>
void rankDictionaryOpenSave(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text;
    generateText(text);
    resize(text, 1000);

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    TRankDictionary rankDictionary(text);
    save(rankDictionary, toCString(tempFilename));

    TRankDictionary rankDictionaryOpen;
    open(rankDictionaryOpen, toCString(tempFilename));
    SEQAN_ASSERT(rankDictionary == rankDictionaryOpen);
}

template <typename TValue>
void rankDictionaryOpenSave(RankDictionary<SequenceBitMask<TValue> > & /*tag*/) {}

SEQAN_TYPED_TEST(RankDictionaryTestCommon, OpenSave)
{
    using namespace seqan;
    {
        RankDictionary<typename TestFixture::TRankDictionarySpec> dictionay;
        rankDictionaryOpenSave(dictionay);
    }
}

// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_clear)
// {
//     using namespace seqan;
// 
//     RankDictionary<WaveletTree<String<Dna> > > dnaTag;
//     WaveletTree<String<Dna5>, void> dna5Tag;
//     WaveletTree<String<AminoAcid>, void> asTag;
//     WaveletTree<String<char>, void> charTag;
//     WaveletTree<String<unsigned char>, void> uCharTag;
//     rankDictionaryClear(dnaTag);
//     rankDictionaryClear(dna5Tag);
//     rankDictionaryClear(asTag);
//     rankDictionaryClear(charTag);
//     rankDictionaryClear(uCharTag);
// }
// 
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_constructor)
// {
//     using namespace seqan;
// 
//     WaveletTree<String<Dna>, void> dnaTag;
//     WaveletTree<String<Dna5>, void> dna5Tag;
//     WaveletTree<String<AminoAcid>, void> asTag;
//     WaveletTree<String<char>, void> charTag;
//     WaveletTree<String<unsigned char>, void> uCharTag;
//     WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharDollarTag;
//     rankDictionaryConstructor(dnaTag);
//     rankDictionaryConstructor(dna5Tag);
//     rankDictionaryConstructor(asTag);
//     rankDictionaryConstructor(charTag);
//     rankDictionaryConstructor(uCharTag);
//     rankDictionaryConstructor(uCharDollarTag);
// }
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_dollar_position)
// {
//     using namespace seqan;
// 
//     WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
//     WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
//     WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
//     WaveletTree<String<char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
//     WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;
//     rankDictionaryDollarPosition(dnaTag);
//     rankDictionaryDollarPosition(dna5Tag);
//     rankDictionaryDollarPosition(asTag);
//     rankDictionaryDollarPosition(charTag);
//     rankDictionaryDollarPosition(uCharTag);
// }
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_dollar_substitute)
// {
//     using namespace seqan;
// 
//     WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
//     WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
//     WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
//     WaveletTree<String<char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
//     WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;
//     rankDictionaryDollarSubstitute(dnaTag);
//     rankDictionaryDollarSubstitute(dna5Tag);
//     rankDictionaryDollarSubstitute(asTag);
//     rankDictionaryDollarSubstitute(charTag);
//     rankDictionaryDollarSubstitute(uCharTag);
// }
// 
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_empty)
// {
//     using namespace seqan;
// 
//     WaveletTree<String<Dna>, void> dnaTag;
//     WaveletTree<String<Dna5>, void> dna5Tag;
//     WaveletTree<String<AminoAcid>, void> asTag;
//     WaveletTree<String<char>, void> charTag;
//     WaveletTree<String<unsigned char>, void> uCharTag;
//     rankDictionaryEmpty(dnaTag);
//     rankDictionaryEmpty(dna5Tag);
//     rankDictionaryEmpty(asTag);
//     rankDictionaryEmpty(charTag);
//     rankDictionaryEmpty(uCharTag);
// }
// 
// 
// 
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_get_character)
// {
//     using namespace seqan;
//     {   
//         WaveletTree<String<Dna>, void> dnaTag;
//         WaveletTree<String<Dna5>, void> dna5Tag;
//         WaveletTree<String<AminoAcid>, void> asTag;
//         WaveletTree<String<signed char>, void> charTag;
//         WaveletTree<String<unsigned char>, void> uCharTag;
// 
//         rankDictionaryGetCharacter(dnaTag);
//         rankDictionaryGetCharacter(dna5Tag);
//         rankDictionaryGetCharacter(asTag);
//         rankDictionaryGetCharacter(charTag);
//         rankDictionaryGetCharacter(uCharTag);
//     }
//     {
//         WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
//         WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
//         WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
//         WaveletTree<String<signed char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
//         WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;
// 
//         rankDictionaryGetCharacter(dnaTag);
//         rankDictionaryGetCharacter(dna5Tag);
//         rankDictionaryGetCharacter(asTag);
//         rankDictionaryGetCharacter(charTag);
//         rankDictionaryGetCharacter(uCharTag);
//     }
// }
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_get_fibre)
// {
//     using namespace seqan;
// 
//     WaveletTree<String<Dna>, void> dnaTag;
//     WaveletTree<String<Dna5>, void> dna5Tag;
//     WaveletTree<String<AminoAcid>, void> asTag;
//     WaveletTree<String<char>, void> charTag;
//     WaveletTree<String<unsigned char>, void> uCharTag;
//     rankDictionaryGetFibre(dnaTag);
//     rankDictionaryGetFibre(dna5Tag);
//     rankDictionaryGetFibre(asTag);
//     rankDictionaryGetFibre(charTag);
//     rankDictionaryGetFibre(uCharTag);
// }
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_get_occ)
// {
//     using namespace seqan;
// 
//     {
//         WaveletTree<String<Dna>, void> dnaTag;
//         WaveletTree<String<Dna5>, void> dna5Tag;
//         WaveletTree<String<AminoAcid>, void> asTag;
//         WaveletTree<String<signed char>, void> charTag;
//         WaveletTree<String<unsigned char>, void> uCharTag;
// 
//         rankDictionaryCountOcc(dnaTag);
//         rankDictionaryCountOcc(dna5Tag);
//         rankDictionaryCountOcc(asTag);
//         rankDictionaryCountOcc(charTag);
//         rankDictionaryCountOcc(uCharTag);
//     }
// //     {
// //         WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
// //         WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
// //         WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
// //         WaveletTree<String<signed char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
// //         WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;
// // 
// //         rankDictionaryCountOcc(dnaTag);
// //         rankDictionaryCountOcc(dna5Tag);
// //         rankDictionaryCountOcc(asTag);
// //         rankDictionaryCountOcc(charTag);
// //         rankDictionaryCountOcc(uCharTag);
// //     }
// }
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_fill_wavelet_tree_)
// {
//     using namespace seqan;
// 
//     WaveletTree<String<Dna>, void> dnaTag;
//     WaveletTree<String<Dna5>, void> dna5Tag;
//     WaveletTree<String<AminoAcid>, void> asTag;
//     WaveletTree<String<signed char>, void> charTag;
//     WaveletTree<String<unsigned char>, void> uCharTag;
// 
//     _rankDictionaryFillWaveletTree(dnaTag);
//     _rankDictionaryFillWaveletTree(dna5Tag);
//     _rankDictionaryFillWaveletTree(asTag);
//     _rankDictionaryFillWaveletTree(charTag);
//     _rankDictionaryFillWaveletTree(uCharTag);
// }
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_open_save)
// {
//     using namespace seqan;
// 
//     WaveletTree<String<Dna5>, void> dna5Tag;
// 
//     rankDictionaryOpenSave(dna5Tag);
// }

#endif  // TESTS_WAVELT_TREE_STRUCTURE_BETA_H_
