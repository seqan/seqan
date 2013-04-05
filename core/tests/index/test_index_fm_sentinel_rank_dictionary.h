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

#ifndef TEST_INDEX_FM_SENTINEL_RANK_DICTIONARY
#define TEST_INDEX_FM_SENTINEL_RANK_DICTIONARY

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

template <typename TAlphabetSpecPair>
class SentinelRankDictionaryTest : public seqan::Test
{
public:
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 0>::Type TSentinelRankDictionarySpec;
};

typedef seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::Dna> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::Dna5> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::AminoAcid> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<char> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<unsigned char> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<signed char> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::Dna> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::Dna5> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::AminoAcid> >, Sentinel> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::Dna> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::Dna5> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<seqan::AminoAcid> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<char> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<unsigned char> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::WaveletTree<signed char> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::Dna> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::Dna5> >, Sentinels> >, seqan::TagList<
            seqan::TagList<SentinelRankDictionary<RankDictionary<seqan::SequenceBitMask<seqan::AminoAcid> >, Sentinels> >
            > > > > > >
            > > >
            > > > > > >
            > > >
        SentinelRankDictionaryTestTypes;


template <typename T>
class SentinelRankDictionaryTestCommon : public SentinelRankDictionaryTest<T>
{};

SEQAN_TYPED_TEST_CASE(SentinelRankDictionaryTestCommon, SentinelRankDictionaryTestTypes);

using namespace seqan;


template <typename TRankDictionary>
void sentinelRankDictionaryConstructor(TRankDictionary & /*tag*/)
{
	{
		TRankDictionary sentinelRankDictionary;
	}
	{
		String<typename Value<TRankDictionary>::Type> text;        
        generateText(text);
		TRankDictionary sentinelRankDictionary(text);
		TRankDictionary sentinelRankDictionary2(sentinelRankDictionary);
        
        SEQAN_ASSERT(sentinelRankDictionary == sentinelRankDictionary2);

		for (unsigned i = 0; i < length(text); ++i)
        {
		    SEQAN_ASSERT_EQ(getValue(sentinelRankDictionary, i), text[i]);
		    SEQAN_ASSERT_EQ(getValue(sentinelRankDictionary2, i), text[i]);
        }
    }
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, Constuctor)
{
    using namespace seqan;

    typename TestFixture::TSentinelRankDictionarySpec dictionay;
    sentinelRankDictionaryConstructor(dictionay);
}

template <typename TRankDictionary>
void sentinelDictionaryClear(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
	TRankDictionary sentinelDictionary(text);

    clear(sentinelDictionary);

	SEQAN_ASSERT_EQ(empty(sentinelDictionary), true);
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, Clear)
{
    using namespace seqan;

    typename TestFixture::TSentinelRankDictionarySpec dictionay;
    sentinelDictionaryClear(dictionay);
}

template <typename TRankDictionary>
void sentinelDictionarySentinelPosition(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
	TRankDictionary sentinelDictionary(text);

    SEQAN_ASSERT_EQ(_getSentinelPosition(sentinelDictionary), length(text));

    setSentinelPosition(sentinelDictionary, 2u);
    SEQAN_ASSERT_EQ(_getSentinelPosition(sentinelDictionary), 2u);
}
    
template <typename TRankDictionary>
void sentinelDictionarySentinelPosition(SentinelRankDictionary<TRankDictionary, Sentinels> & /*tag*/)
{
	String<typename Value<SentinelRankDictionary<TRankDictionary, Sentinels> >::Type> text = "ACGTNACGTNACGTN";
	SentinelRankDictionary<TRankDictionary, Sentinels> sentinelDictionary(text);

    String<unsigned> sentinelPos;

    SEQAN_ASSERT_EQ(_getSentinelPosition(sentinelDictionary), sentinelPos);

    RankSupportBitString<> bitString;
    resize(bitString, length(text), 0);
    setBitTo(bitString, 2u, 1);
    setSentinelPosition(sentinelDictionary, bitString);
    String<unsigned> temp;
    appendValue(temp, 2);
    SEQAN_ASSERT_EQ(_getSentinelPosition(sentinelDictionary), temp);
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, SentinelPosition)
{
    using namespace seqan;

    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelDictionarySentinelPosition(dictionay);
    }
}

template <typename TRankDictionary>
void sentinelRankDictionaryEmpty(TRankDictionary & /*tag*/)
{
    {
        TRankDictionary sentinelRankDictionary;
        
        SEQAN_ASSERT_EQ(empty(sentinelRankDictionary), true);
    }
    {
        String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
        TRankDictionary sentinelRankDictionary(text);
        
        SEQAN_ASSERT_EQ(empty(sentinelRankDictionary), false);
    }
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, Empty)
{
    using namespace seqan;

    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryEmpty(dictionay);
    }
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryEmpty(dictionay);
    }
}

template <typename TRankDictionary>
void sentinelRankDictionaryGetValue(TRankDictionary & /*tag*/)
{
    String<typename Value<TRankDictionary>::Type> text;
    generateText(text);

    TRankDictionary sentinelRankDictionary(text); 

    for (unsigned i = 0; i < length(text); ++i)
    {
        SEQAN_ASSERT_EQ(getValue(sentinelRankDictionary, i), text[i]);
    }
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, GetValue)
{
    using namespace seqan;

    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetValue(dictionay);
    }
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetValue(dictionay);
    }
}

template <typename TRankDictionary>
void sentinelRankDictionaryGetFibre(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text = "ACGTNACGTNACGTN";
	TRankDictionary sentinelRankDictionary(text);

    typename Fibre<TRankDictionary, FibreRankDictionary>::Type & dictionary = getFibre(sentinelRankDictionary, FibreRankDictionary());

	SEQAN_ASSERT_EQ(empty(getFibre(sentinelRankDictionary, FibreRankDictionary())), false);

    clear(dictionary);

	SEQAN_ASSERT_EQ(empty(getFibre(sentinelRankDictionary, FibreRankDictionary())), true);
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, GetFibre)
{
    using namespace seqan;

    using namespace seqan;
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetFibre(dictionay);
    }
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryGetFibre(dictionay);
    }
}

template <typename TRankDictionary>
void sentinelRankDictionaryCountOcc(TRankDictionary & /*tag*/)
{
    typedef typename Value<TRankDictionary>::Type TChar;
    String<TChar> text;
    generateText(text);
    text[0] = 'A';

    resize(text, 1000);

    TRankDictionary sentinelRankDictionary(text);
    setSentinelSubstitute(sentinelRankDictionary, 'A');
    setSentinelPosition(sentinelRankDictionary, 0);

    for (int i = MinValue<TChar>::VALUE; i <= MaxValue<TChar>::VALUE; ++i)
    {
        unsigned counter = 0;
        for (unsigned j = 0; j < length(text); ++j)
        {
            if (text[j] == (TChar)i)
                ++counter;
            if ((TChar)i == 'A' && j == 0)
                --counter;
            SEQAN_ASSERT_EQ(countOccurrences(sentinelRankDictionary, (TChar)i, j), counter);
        }
    }
}
template <typename TRankDictionary>
void sentinelRankDictionaryCountOcc(SentinelRankDictionary<TRankDictionary, Sentinels> & /*tag*/)
{
    typedef typename Value<TRankDictionary>::Type TChar;
    String<TChar> text;
    generateText(text);
    text[0] = 'A';
    text[99] = 'A';
    text[999] = 'A';

    resize(text, 1000);
    
    SentinelRankDictionary<TRankDictionary, Sentinels> sentinelRankDictionary(text);
    setSentinelSubstitute(sentinelRankDictionary, 'A');
    setBitTo(sentinelRankDictionary.sentinelPosition, 0, 1);
    setBitTo(sentinelRankDictionary.sentinelPosition, 99, 1);
    setBitTo(sentinelRankDictionary.sentinelPosition, 999, 1);
    _updateRanks(sentinelRankDictionary.sentinelPosition);


    for (int i = MinValue<TChar>::VALUE; i <= MaxValue<TChar>::VALUE; ++i)
    {
        unsigned counter = 0;
        for (unsigned j = 0; j < length(text); ++j)
        {
            if (text[j] == (TChar)i)
                ++counter;
            if ((TChar)i == 'A' && (j == 0 || j == 99 || j == 999))
                --counter;
            SEQAN_ASSERT_EQ(countOccurrences(sentinelRankDictionary, (TChar)i, j), counter);
        }
    }
}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, CountOcc)
{
    using namespace seqan;
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryCountOcc(dictionay);
    }
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryCountOcc(dictionay);
    }
}


// // template <typename TString, typename TSpec>
// // void sentinelRankDictionaryCountOcc(WaveletTree<TString, FmiDollarSubstituted<TSpec> > & /*tag*/)
// // {
// //     typedef WaveletTree<TString, FmiDollarSubstituted<TSpec> > TRankDictionary;
// // 	typedef typename Fibre<TRankDictionary, FibreTreeStructure>::Type TRankDictionaryStructure;
// // 	typedef typename Fibre<TRankDictionaryStructure, FibreTreeStructureEncoding>::Type TRankDictionaryVertieces;
// // 	typedef typename Value<TRankDictionaryVertieces>::Type TRankDictionaryVertex;
// // 	typedef typename Value<TRankDictionary>::Type TChar;
// // 	typedef typename Value<TRankDictionaryVertex, 2>::Type TPos;
// // 
// // 	{
// // 		String<TChar> text;
// //  		generateText(text);
// //  		resize(text, 1000);
// // 
// // 		TRankDictionary sentinelRankDictionary(text);
// // 		setDollarSubstitute(sentinelRankDictionary, getCharacter(sentinelRankDictionary, 0u));
// // 		setDollarPosition(sentinelRankDictionary, 0u);
// // 
// //         for (int i = MinValue<TChar>::VALUE; i <= MaxValue<TChar>::VALUE; ++i)
// //         {
// //             unsigned counter = 0;
// //             for (unsigned j = 0; j < length(text); ++j)
// //             {
// //                 if(text[j] == (TChar)i)
// //                 {
// //                     ++counter;
// //                     if(text[j] == getDollarSubstitute(sentinelRankDictionary) && j == getDollarPosition(sentinelRankDictionary))
// //                         --counter;
// //                 }
// // 		        SEQAN_ASSERT_EQ(countOccurrences(sentinelRankDictionary, (TChar)i, j), counter);
// //             }
// //         }
// //     }
// // }
// 

template <typename TRankDictionary>
void sentinelRankDictionaryOpenSave(TRankDictionary & /*tag*/)
{
	String<typename Value<TRankDictionary>::Type> text;
    generateText(text);
    resize(text, 1000);

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    TRankDictionary sentinelRankDictionary(text);
    sentinelRankDictionary.sentinelSubstitute = 'A';
    save(sentinelRankDictionary, toCString(tempFilename));

    TRankDictionary sentinelRankDictionaryOpen;
    open(sentinelRankDictionaryOpen, toCString(tempFilename));
    SEQAN_ASSERT(sentinelRankDictionary == sentinelRankDictionaryOpen);
}

template <typename TValue>
void sentinelRankDictionaryOpenSave(RankDictionary<SequenceBitMask<TValue> > & /*tag*/) {}

SEQAN_TYPED_TEST(SentinelRankDictionaryTestCommon, OpenSave)
{
    using namespace seqan;
    {
        typename TestFixture::TSentinelRankDictionarySpec dictionay;
        sentinelRankDictionaryOpenSave(dictionay);
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
//     sentinelRankDictionaryClear(dnaTag);
//     sentinelRankDictionaryClear(dna5Tag);
//     sentinelRankDictionaryClear(asTag);
//     sentinelRankDictionaryClear(charTag);
//     sentinelRankDictionaryClear(uCharTag);
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
//     sentinelRankDictionaryConstructor(dnaTag);
//     sentinelRankDictionaryConstructor(dna5Tag);
//     sentinelRankDictionaryConstructor(asTag);
//     sentinelRankDictionaryConstructor(charTag);
//     sentinelRankDictionaryConstructor(uCharTag);
//     sentinelRankDictionaryConstructor(uCharDollarTag);
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
//     sentinelRankDictionaryDollarPosition(dnaTag);
//     sentinelRankDictionaryDollarPosition(dna5Tag);
//     sentinelRankDictionaryDollarPosition(asTag);
//     sentinelRankDictionaryDollarPosition(charTag);
//     sentinelRankDictionaryDollarPosition(uCharTag);
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
//     sentinelRankDictionaryDollarSubstitute(dnaTag);
//     sentinelRankDictionaryDollarSubstitute(dna5Tag);
//     sentinelRankDictionaryDollarSubstitute(asTag);
//     sentinelRankDictionaryDollarSubstitute(charTag);
//     sentinelRankDictionaryDollarSubstitute(uCharTag);
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
//     sentinelRankDictionaryEmpty(dnaTag);
//     sentinelRankDictionaryEmpty(dna5Tag);
//     sentinelRankDictionaryEmpty(asTag);
//     sentinelRankDictionaryEmpty(charTag);
//     sentinelRankDictionaryEmpty(uCharTag);
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
//         sentinelRankDictionaryGetCharacter(dnaTag);
//         sentinelRankDictionaryGetCharacter(dna5Tag);
//         sentinelRankDictionaryGetCharacter(asTag);
//         sentinelRankDictionaryGetCharacter(charTag);
//         sentinelRankDictionaryGetCharacter(uCharTag);
//     }
//     {
//         WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
//         WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
//         WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
//         WaveletTree<String<signed char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
//         WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;
// 
//         sentinelRankDictionaryGetCharacter(dnaTag);
//         sentinelRankDictionaryGetCharacter(dna5Tag);
//         sentinelRankDictionaryGetCharacter(asTag);
//         sentinelRankDictionaryGetCharacter(charTag);
//         sentinelRankDictionaryGetCharacter(uCharTag);
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
//     sentinelRankDictionaryGetFibre(dnaTag);
//     sentinelRankDictionaryGetFibre(dna5Tag);
//     sentinelRankDictionaryGetFibre(asTag);
//     sentinelRankDictionaryGetFibre(charTag);
//     sentinelRankDictionaryGetFibre(uCharTag);
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
//         sentinelRankDictionaryCountOcc(dnaTag);
//         sentinelRankDictionaryCountOcc(dna5Tag);
//         sentinelRankDictionaryCountOcc(asTag);
//         sentinelRankDictionaryCountOcc(charTag);
//         sentinelRankDictionaryCountOcc(uCharTag);
//     }
// //     {
// //         WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
// //         WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
// //         WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
// //         WaveletTree<String<signed char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
// //         WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;
// // 
// //         sentinelRankDictionaryCountOcc(dnaTag);
// //         sentinelRankDictionaryCountOcc(dna5Tag);
// //         sentinelRankDictionaryCountOcc(asTag);
// //         sentinelRankDictionaryCountOcc(charTag);
// //         sentinelRankDictionaryCountOcc(uCharTag);
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
//     _sentinelRankDictionaryFillWaveletTree(dnaTag);
//     _sentinelRankDictionaryFillWaveletTree(dna5Tag);
//     _sentinelRankDictionaryFillWaveletTree(asTag);
//     _sentinelRankDictionaryFillWaveletTree(charTag);
//     _sentinelRankDictionaryFillWaveletTree(uCharTag);
// }
// 
// SEQAN_DEFINE_TEST(test_wavelet_tree_open_save)
// {
//     using namespace seqan;
// 
//     WaveletTree<String<Dna5>, void> dna5Tag;
// 
//     sentinelRankDictionaryOpenSave(dna5Tag);
// }

#endif  // TEST_INDEX_FM_SENTINEL_RANK_DICTIONARY

