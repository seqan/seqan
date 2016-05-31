// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/reduced_aminoacid.h>
#include <seqan/index.h>

#include "test_index_helpers.h"

using namespace seqan;

// ==========================================================================
// Metafunctions
// ==========================================================================

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

//namespace seqan {
//template <typename TValue>
//struct Size<RankDictionary<TValue, Levels<unsigned> > >
//{
//    typedef unsigned    Type;
//};
//}

// ==========================================================================
// Types
// ==========================================================================

// --------------------------------------------------------------------------
// RankDictionary Types
// --------------------------------------------------------------------------

typedef
        //TagList<RankDictionary<bool,             Levels<void, LevelsPrefixRDConfig<> > >,
        TagList<RankDictionary<Dna,             Levels<void, LevelsPrefixRDConfig<> > >,
        TagList<RankDictionary<Rna,             Levels<void, LevelsPrefixRDConfig<> > >,
        TagList<RankDictionary<Dna5,             Levels<void, LevelsPrefixRDConfig<> > >,
        TagList<RankDictionary<Rna5,             Levels<void, LevelsPrefixRDConfig<> > >,
        TagList<RankDictionary<SimpleType<unsigned char, ReducedAminoAcid_<Murphy10> >,             Levels<void, LevelsPrefixRDConfig<> > >,
        TagList<RankDictionary<AminoAcid,             Levels<void, LevelsPrefixRDConfig<> > >,
        TagList<RankDictionary<char,             Levels<void, LevelsPrefixRDConfig<> > >
    > > > > > > >// >
    RankDictionaryPrefixTypes;

typedef
/*TagList<RankDictionary<bool,            Naive<> >,,
TagList<RankDictionary<Rna5,            Levels<> >,
TagList<RankDictionary<Dna,             Levels<> >,
TagList<RankDictionary<Rna,             Levels<> >,
TagList<RankDictionary<Dna5,            Levels<> >,
TagList<RankDictionary<Rna5,            Levels<> >,
TagList<RankDictionary<SimpleType<unsigned char, ReducedAminoAcid_<Murphy10> >,            Levels<> >,
TagList<RankDictionary<AminoAcid,       Levels<> >,
TagList<RankDictionary<char,            Levels<> >,
TagList<RankDictionary<Dna,             WaveletTree<> >,
TagList<RankDictionary<Dna5,            WaveletTree<> >,
TagList<RankDictionary<DnaQ,            WaveletTree<> >,
TagList<RankDictionary<Dna5Q,           WaveletTree<> >,
TagList<RankDictionary<AminoAcid,       WaveletTree<> >,
TagList<RankDictionary<char,            WaveletTree<> >,
TagList<RankDictionary<unsigned char,   WaveletTree<> >,*/
    RankDictionaryPrefixTypes
    //> > > > > > > > > > > > > > > >
    RankDictionaryTypes;

// ==========================================================================
// Test Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class RankDictionaryTest
// --------------------------------------------------------------------------


template <typename TRankDictionary>
class RankDictionaryTest : public Test
{
public:
    typedef TRankDictionary                             TRankDict;
    typedef typename Value<TRankDict>::Type             TValue;
    typedef typename Size<TValue>::Type                 TValueSize;
    typedef String<TValue>                              TText;
    typedef typename Iterator<TText, Standard>::Type    TTextIterator;

    TValueSize      alphabetSize;
    TText           text;
    TTextIterator   textBegin;
    TTextIterator   textEnd;

    RankDictionaryTest() :
        alphabetSize(ValueSize<TValue>::VALUE), text(), textBegin(), textEnd()
    {}

    void setUp()
    {
        //createText(text, TValue());
        generateText(text, 50000);
        //text = "ATAGACTCTCGCAATTGGAAGCCTAATAAT"; // AUAGACUCUCGCAAUUGGAAUG
        // AUAGACUCUCGCAAUUGGAAGCCUAAUAAUG CCGG CACUCCGUGGAUCAA
        //std::cout << text << std::endl;
        textBegin = begin(text, Standard());
        textEnd = end(text, Standard());
    }
};

/*template <typename TRankDictionary>
class RankDictionaryPrefixTest : public RankDictionaryTest<TRankDictionary> {};*/

SEQAN_TYPED_TEST_CASE(RankDictionaryTest, RankDictionaryTypes);
//SEQAN_TYPED_TEST_CASE(RankDictionaryPrefixTest, RankDictionaryPrefixTypes);

// ==========================================================================
// Tests
// ==========================================================================

// ----------------------------------------------------------------------------
// Test RankDictionary()
// ----------------------------------------------------------------------------

/*SEQAN_TYPED_TEST(RankDictionaryTest, Constructor)
{
    typename TestFixture::TRankDict dict(this->text);
}*/

// ----------------------------------------------------------------------------
// Test createRankDictionary()
// ----------------------------------------------------------------------------

/*SEQAN_TYPED_TEST(RankDictionaryTest, CreateRankDictionary)
{
    typename TestFixture::TRankDict dict;
    createRankDictionary(dict, this->text);
}*/

// ----------------------------------------------------------------------------
// Test clear() and empty()
// ----------------------------------------------------------------------------

/*SEQAN_TYPED_TEST(RankDictionaryTest, ClearEmpty)
{
    typename TestFixture::TRankDict dict;

    SEQAN_ASSERT(empty(dict));
    createRankDictionary(dict, this->text);
    SEQAN_ASSERT_NOT(empty(dict));
    clear(dict);
    SEQAN_ASSERT(empty(dict));
}*/

// ----------------------------------------------------------------------------
// Test getValue()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(RankDictionaryTest, GetValue)
{
    typedef typename TestFixture::TTextIterator         TTextIterator;

    typename TestFixture::TRankDict dict(this->text);
    /*std::cout << "TEXT:    " << this->text << std::endl;

    std::cout << std::bitset<64>(dict.ranks[0].values[0].i) << std::endl;
    std::cout << std::bitset<64>(dict.ranks[0].values[1].i) << std::endl;
    std::cout << std::bitset<64>(dict.ranks[0].values[2].i) << std::endl;
    std::cout << std::bitset<64>(dict.ranks[0].values[3].i) << std::endl;*/

    for (TTextIterator textIt = this->textBegin; textIt != this->textEnd; ++textIt)
        SEQAN_ASSERT_EQ(getValue(dict, (unsigned long)(textIt - this->textBegin)), value(textIt));
}

// ----------------------------------------------------------------------------
// Test getRank()
// ----------------------------------------------------------------------------

/*SEQAN_TYPED_TEST(RankDictionaryTest, GetRank)
{
    typedef typename TestFixture::TValueSize            TValueSize;
    typedef typename TestFixture::TText                 TText;
    typedef typename TestFixture::TTextIterator         TTextIterator;
    typedef typename Size<TText>::Type                  TTextSize;
    typedef String<TTextSize>                           TPrefixSum;

    typename TestFixture::TRankDict dict(this->text);

    // The prefix sum is built while scanning the text.
    TPrefixSum prefixSum;

    resize(prefixSum, this->alphabetSize, 0);

    // Scan the text.
    for (TTextIterator textIt = this->textBegin; textIt != this->textEnd; ++textIt)
    {
        // Update the prefix sum.
        prefixSum[ordValue(value(textIt))]++;

        // Check the rank for all alphabet symbols.
        for (TValueSize c = 0; c < this->alphabetSize; ++c)
        {
            unsigned long pos = textIt - this->textBegin;
            SEQAN_ASSERT_EQ(getRank(dict, pos, c), prefixSum[c]);
        }
    }
}*/

/*SEQAN_TYPED_TEST(RankDictionaryPrefixTest, GetCumulativeRank)
{
    typedef typename TestFixture::TValueSize            TValueSize;
    typedef typename TestFixture::TText                 TText;
    typedef typename TestFixture::TTextIterator         TTextIterator;
    typedef typename Size<TText>::Type                  TTextSize;
    typedef String<TTextSize>                           TPrefixSum;

    typename TestFixture::TRankDict dict(this->text);

    // The prefix sum is built while scanning the text.
    TPrefixSum prefixSum;

    typedef typename Value<TText>::Type TValue;

    resize(prefixSum, this->alphabetSize, 0);

    // Scan the text.
    for (TTextIterator textIt = this->textBegin; textIt != this->textEnd; ++textIt)
    {
        // Update the prefix sum.
        prefixSum[ordValue(value(textIt))]++;

        // Check the rank for all alphabet symbols.
        unsigned long smallerNaive = 0;
        for (TValueSize c = 0; c < this->alphabetSize; ++c)
        {
            unsigned long smaller;
            unsigned long rank = getRank(dict, (unsigned long)(textIt - this->textBegin), TValue((uint16_t) c), smaller); // TODO: getRank!!! uint16_t cast???
            //SEQAN_ASSERT_EQ(smaller, smallerNaive);
            SEQAN_ASSERT_EQ(rank, prefixSum[c]);
            smallerNaive += prefixSum[c];
        }
    }
}*/

/*SEQAN_TYPED_TEST(RankDictionaryTest, GetRankWithPrefix)
{
    typedef typename TestFixture::TValueSize            TValueSize;
    typedef typename TestFixture::TText                 TText;
    typedef typename TestFixture::TTextIterator         TTextIterator;
    typedef typename Size<TText>::Type                  TTextSize;
    typedef String<TTextSize>                           TPrefixSum;

    typename TestFixture::TRankDict dict(this->text);

    // The prefix sum is built while scanning the text.
    TPrefixSum prefixSum;

    resize(prefixSum, this->alphabetSize, 0);

    // Scan the text.
    for (TTextIterator textIt = this->textBegin; textIt != this->textEnd; ++textIt)
    {
        // Update the prefix sum.
        prefixSum[ordValue(value(textIt))]++;

        // Check the rank for all alphabet symbols.
        unsigned long smallerNaive = 0;
        for (TValueSize c = 0; c < this->alphabetSize; ++c)
        {
            //std::cout << "prefixSum: " << prefixSum[c] << std::endl;
            //std::cout << "unsigned long: " << (unsigned long)(textIt - this->textBegin) << std::endl;
            unsigned long pos = textIt - this->textBegin;
            //if (pos == 128)
            //    std::cout << "xxx" << std::endl;
            unsigned long smaller = 0;
            //std::cout << prefixSum[c] << std::endl;
            //unsigned long rank = getRank(dict, pos, c);

            unsigned long rank = getRank(dict, pos, c, smaller);
            SEQAN_ASSERT_EQ(rank, prefixSum[c]);
            SEQAN_ASSERT_EQ(smaller, smallerNaive);
            smallerNaive += prefixSum[c];
        }
        //std::cout << prefixSum[0] << " " << prefixSum[1] << " " << prefixSum[2] << " " << prefixSum[3] << std::endl;
    }
}*/

// ----------------------------------------------------------------------------
// Test setValue()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): rename it as assignValue()

// ----------------------------------------------------------------------------
// Test length()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Test resize()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Test reserve()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Test open() and save()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Test Size<>
// ----------------------------------------------------------------------------

//SEQAN_DEFINE_TEST(test_rss_sizeof)
//{
//    typedef Dna                                         TAlphabet;
//    typedef Alloc<unsigned>                             TTextSpec;
//    typedef String<TAlphabet, TTextSpec>                TText;
//
//    typedef Levels<TAlphabet, unsigned>              TRankDictionarySpec;
//    typedef RankDictionary<TRankDictionarySpec>         TRankDictionary;
//
//    TRankSupport rs;
//
//    std::cout << "sizeof(Block): " << sizeof(rs.block) << std::endl;
//
//    std::cout << "bits(Block): " << BitsPerValue<TRankSupport::TBlock>::VALUE << std::endl;
//    std::cout << "length(Block): " << length(rs.block) << std::endl;
//    std::cout << "capacity(Block): " << capacity(rs.block) << std::endl;
//    std::cout << std::endl;
//
//    std::cout << "sizeof(SuperBlock): " << sizeof(rs.sblock) << std::endl;
//    std::cout << "bits(SuperBlock): " << BitsPerValue<TRankSupport::TSuperBlock>::VALUE << std::endl;
//    std::cout << "length(SuperBlock): " << length(rs.sblock) << std::endl;
//    std::cout << "capacity(SuperBlock): " << capacity(rs.sblock) << std::endl;
//    std::cout << std::endl;
//
//    std::cout << "sizeof(RankSupport): " << sizeof(rs) << std::endl;
//    std::cout << "bits(RankSupport): " << BitsPerValue<TRankSupport>::VALUE << std::endl;
//    std::cout << std::endl;
//}

//SEQAN_DEFINE_TEST(test_rss_resize)
//{
//    typedef Dna                                         TAlphabet;
//    typedef Alloc<unsigned>                             TTextSpec;
//    typedef String<TAlphabet, TTextSpec>                TText;
//
//    typedef Levels<TAlphabet, unsigned>              TRankDictionarySpec;
//    typedef RankDictionary<TRankDictionarySpec>         TRankDictionary;
//
////    TText text = "ACGTNACGTNACGTNACGTNA";
//    TText text = "ACGTACGTACGTACGTACGTACGTACGTACGT";
////    TText text = "ACGTACGTACGTACGTACGTACGTACGTACGTCCCCCCCCCCCCCCC";
//
//    TRankDictionary dict(text);
////    createRankDictionary(dict, text);
//
//    std::cout << "Text: " << text << std::endl;
////    std::cout << "Block: " << rs.block << std::endl;
//
//    for (unsigned i = 0; i < 10; i++)
//        for (unsigned char c = 0; c < 4; c++)
//            std::cout << "getRank(" << Dna(c) << ", " << i << "): " << getRank(dict, i, Dna(c)) << std::endl;
//
//    std::cout << std::endl;
//}

// ========================================================================== 
// Functions
// ========================================================================== 

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
