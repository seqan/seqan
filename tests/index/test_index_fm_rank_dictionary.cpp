// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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

typedef SimpleType<unsigned char, ReducedAminoAcid_<Murphy10> >         ReducedMurphy10;

typedef Levels<void, LevelsPrefixRDConfig<uint32_t, Alloc<>, 1, 1> > Prefix1Level;
typedef Levels<void, LevelsPrefixRDConfig<uint32_t, Alloc<>, 2, 2> > Prefix2Level;
typedef Levels<void, LevelsPrefixRDConfig<uint32_t, Alloc<>, 3, 3> > Prefix3Level;

typedef Levels<void, LevelsRDConfig<uint32_t, Alloc<>, 1, 1> >       Default1Level;
typedef Levels<void, LevelsRDConfig<uint32_t, Alloc<>, 2, 2> >       Default2Level;
typedef Levels<void, LevelsRDConfig<uint32_t, Alloc<>, 3, 3> >       Default3Level;

typedef
    TagList<RankDictionary<bool,            Prefix1Level>,
    TagList<RankDictionary<Dna,             Prefix1Level>,
    TagList<RankDictionary<Dna5,            Prefix1Level>,
    TagList<RankDictionary<Dna5Q,           Prefix1Level>,
    TagList<RankDictionary<ReducedMurphy10, Prefix1Level>,
    TagList<RankDictionary<AminoAcid,       Prefix1Level>,
    TagList<RankDictionary<bool,            Prefix2Level>,
    TagList<RankDictionary<Dna,             Prefix2Level>,
    TagList<RankDictionary<Dna5,            Prefix2Level>,
    TagList<RankDictionary<Dna5Q,           Prefix2Level>,
    TagList<RankDictionary<ReducedMurphy10, Prefix2Level>,
    TagList<RankDictionary<AminoAcid,       Prefix2Level>,
    TagList<RankDictionary<bool,            Prefix3Level>,
    TagList<RankDictionary<Dna,             Prefix3Level>,
    TagList<RankDictionary<Dna5,            Prefix3Level>,
    TagList<RankDictionary<Dna5Q,           Prefix3Level>,
    TagList<RankDictionary<ReducedMurphy10, Prefix3Level>,
    TagList<RankDictionary<AminoAcid,       Prefix3Level>,
    TagList<RankDictionary<bool,            WaveletTree<> >,
    TagList<RankDictionary<Dna,             WaveletTree<> >,
    TagList<RankDictionary<Dna5,            WaveletTree<> >,
    TagList<RankDictionary<Dna5Q,           WaveletTree<> >,
    TagList<RankDictionary<ReducedMurphy10, WaveletTree<> >,
    TagList<RankDictionary<AminoAcid,       WaveletTree<> >
    > > > > > > > > > > > > > > > > > > > > > > > >
    RankDictionaryPrefixSumTypes;

typedef
    TagList<RankDictionary<bool,            Naive<> >,
    TagList<RankDictionary<bool,            Default1Level>,
    TagList<RankDictionary<Dna,             Default1Level>,
    TagList<RankDictionary<Dna5,            Default1Level>,
    TagList<RankDictionary<Dna5Q,           Default1Level>,
    TagList<RankDictionary<ReducedMurphy10, Default1Level>,
    TagList<RankDictionary<AminoAcid,       Default1Level>,
    TagList<RankDictionary<bool,            Default2Level>,
    TagList<RankDictionary<Dna,             Default2Level>,
    TagList<RankDictionary<Dna5,            Default2Level>,
    TagList<RankDictionary<Dna5,            Default2Level>,
    TagList<RankDictionary<ReducedMurphy10, Default2Level>,
    TagList<RankDictionary<AminoAcid,       Default2Level>,
    TagList<RankDictionary<bool,            Default3Level>,
    TagList<RankDictionary<Dna,             Default3Level>,
    TagList<RankDictionary<Dna5,            Default3Level>,
    TagList<RankDictionary<Dna5Q,           Default3Level>,
    TagList<RankDictionary<ReducedMurphy10, Default3Level>,
#ifndef __alpha__ // NOTE(h-2): fails on alpha for unknown reasons
    TagList<RankDictionary<AminoAcid,       Default3Level>,
#endif
    RankDictionaryPrefixSumTypes
    > > > > > > > > > > > > > > > > > >
#ifndef __alpha__
    >
#endif
    RankDictionaryAllTypes;

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
        generateText(text, 3947);
        textBegin = begin(text, Standard());
        textEnd = end(text, Standard());
    }
};

template <typename TRankDictionary>
class RankDictionaryPrefixTest : public RankDictionaryTest<TRankDictionary> {};

SEQAN_TYPED_TEST_CASE(RankDictionaryTest, RankDictionaryAllTypes);
SEQAN_TYPED_TEST_CASE(RankDictionaryPrefixTest, RankDictionaryPrefixSumTypes);

// ==========================================================================
// Tests
// ==========================================================================

// ----------------------------------------------------------------------------
// Test RankDictionary()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(RankDictionaryTest, Constructor)
{
    typename TestFixture::TRankDict dict(this->text);
}

// ----------------------------------------------------------------------------
// Test createRankDictionary()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(RankDictionaryTest, CreateRankDictionary)
{
    typename TestFixture::TRankDict dict;
    createRankDictionary(dict, this->text);
}

// ----------------------------------------------------------------------------
// Test clear() and empty()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(RankDictionaryTest, ClearEmpty)
{
    typename TestFixture::TRankDict dict;

    SEQAN_ASSERT(empty(dict));
    createRankDictionary(dict, this->text);
    SEQAN_ASSERT_NOT(empty(dict));
    clear(dict);
    SEQAN_ASSERT(empty(dict));
}

// ----------------------------------------------------------------------------
// Test getValue()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(RankDictionaryTest, GetValue)
{
    typedef typename TestFixture::TTextIterator         TTextIterator;

    typename TestFixture::TRankDict dict(this->text);

    for (TTextIterator textIt = this->textBegin; textIt != this->textEnd; ++textIt)
        SEQAN_ASSERT_EQ(getValue(dict, (unsigned long)(textIt - this->textBegin)), value(textIt));
}

// ----------------------------------------------------------------------------
// Test getRank()
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(RankDictionaryTest, GetRank)
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
            SEQAN_ASSERT_EQ(getRank(dict, (unsigned long)(textIt - this->textBegin), c), prefixSum[c]);
    }
}

SEQAN_TYPED_TEST(RankDictionaryPrefixTest, GetPrefixRank)
{
    typedef typename TestFixture::TValueSize TValueSize;
    typedef typename TestFixture::TText TText;
    typedef typename TestFixture::TTextIterator TTextIterator;
    typedef typename Size<TText>::Type TTextSize;
    typedef String<TTextSize> TPrefixSum;

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
            unsigned long smaller;
            SEQAN_ASSERT_EQ(getRank(dict, (unsigned long)(textIt - this->textBegin), c, smaller), prefixSum[c]);
            SEQAN_ASSERT_EQ(smaller, smallerNaive);
            smallerNaive += prefixSum[c];
        }
    }
}

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
