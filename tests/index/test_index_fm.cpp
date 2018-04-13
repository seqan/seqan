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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/index.h>

#include "test_index_helpers.h"

using namespace seqan;

// ==========================================================================
// Types
// ========================================================================== 

// --------------------------------------------------------------------------
// FMIndex Configs
// --------------------------------------------------------------------------

template <typename TSpec = void, typename TLengthSum = size_t>
struct WTFMIndexConfig : FMIndexConfig<TSpec, TLengthSum> {};

template <typename TSpec = void, typename TLengthSum = size_t>
struct SmallWTFMIndexConfig : FMIndexConfig<TSpec, TLengthSum>
{
    typedef TLengthSum                                           LengthSum;
    typedef Naive<TSpec, RDConfig<TLengthSum, Alloc<>, 1, 0> >   Sentinels;
};

template <typename TSpec = void, typename TLengthSum = size_t>
struct SmallLVFMIndexConfig : FMIndexConfig<TSpec, TLengthSum>
{
    typedef TLengthSum                                                  LengthSum;
    typedef Levels<TSpec, LevelsRDConfig<TLengthSum, Alloc<>, 1, 0> >   Bwt;
    typedef Naive<TSpec, RDConfig<TLengthSum, Alloc<>, 1, 0> >          Sentinels;
};

template <typename TSpec = void, typename TLengthSum = size_t>
struct PrefixLVFMIndexConfig : FMIndexConfig<TSpec, TLengthSum>
{
    typedef TLengthSum                                                       LengthSum;
    typedef Levels<TSpec, LevelsPrefixRDConfig<TLengthSum, Alloc<>, 1, 0> >  Bwt;
    typedef Naive<TSpec, RDConfig<TLengthSum, Alloc<>, 1, 0> >               Sentinels;
};

// --------------------------------------------------------------------------
// FMIndex Specs
// --------------------------------------------------------------------------

typedef FMIndex<void, WTFMIndexConfig<> >       WTFMIndex;
typedef FMIndex<void, SmallWTFMIndexConfig<> >  SmallWTFMIndex;
typedef FMIndex<void, SmallLVFMIndexConfig<> >  SmallLVFMIndex;
typedef FMIndex<void, PrefixLVFMIndexConfig<> > PrefixLVFMIndex;

// --------------------------------------------------------------------------
// FMIndex Types
// --------------------------------------------------------------------------

typedef
    TagList<Index<DnaString, WTFMIndex>,
    TagList<Index<CharString, WTFMIndex>,
    TagList<Index<StringSet<CharString>, WTFMIndex>,
    TagList<Index<StringSet<CharString>, SmallWTFMIndex>,
    TagList<Index<StringSet<DnaString>, SmallLVFMIndex>,
    TagList<Index<String<bool>, PrefixLVFMIndex>,
    TagList<Index<DnaString, PrefixLVFMIndex>,
    TagList<Index<CharString, PrefixLVFMIndex>,
    TagList<Index<StringSet<DnaString>, PrefixLVFMIndex>
    > > > > > > > > >
    FMIndexTypes2;

// ========================================================================== 
// Test Classes
// ========================================================================== 

// --------------------------------------------------------------------------
// Class LFTest
// --------------------------------------------------------------------------

template <typename TFMIndex>
class LFTest : public FibreTest<TFMIndex, FibreLF> {};

SEQAN_TYPED_TEST_CASE(LFTest, FMIndexTypes2);

// --------------------------------------------------------------------------
// Class CSATest
// --------------------------------------------------------------------------

template <typename TFMIndex>
class CSATest : public FibreTest<TFMIndex, FibreSA> {};

SEQAN_TYPED_TEST_CASE(CSATest, FMIndexTypes2);

// ==========================================================================
// LFTable Tests
// ========================================================================== 

// --------------------------------------------------------------------------
// Test LF(pos)
// --------------------------------------------------------------------------

//SEQAN_TYPED_TEST(LFTest, LFPos)
//{
//    typedef typename Iterator<typename TestFixture::TText>::Type    TIter;
//
//    unsigned pos = 0;
//    TIter textBegin = begin(this->text, Standard());
//    for (TIter textIt = end(this->text, Standard()) - 1; textIt >= textBegin; goPrevious(textIt))
//    {
//        pos = this->fibre(pos);
//    }
//}

// --------------------------------------------------------------------------
// Test LF(pos, val)
// --------------------------------------------------------------------------

//SEQAN_TYPED_TEST(LFTest, LFPosVal)
//}

// --------------------------------------------------------------------------
// Test isSentinel()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(LFTest, IsSentinel)
{
    typedef typename Size<typename TestFixture::TFibre>::Type   TSize;

    TSize sentinels = 0;
    for (TSize pos = 0; pos < bwtLength(this->text); ++pos)
        sentinels += isSentinel(this->fibre, pos);

    // Assert that there is exactly one sentinel per text in the collection.
    SEQAN_ASSERT_EQ(sentinels, countSequences(this->text));
}

// ==========================================================================
// CompressedSA Tests
// ==========================================================================

// --------------------------------------------------------------------------
// Test indexSA()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(CSATest, IndexSA)
{
// NOTE(esiragusa): Actual behavior is:
    SEQAN_ASSERT_EQ(length(this->fibre), bwtLength(this->text));
// NOTE(esiragusa): Correct behavior should be:
//    SEQAN_ASSERT_EQ(length(this->fibre), lengthSum(this->text));

// NOTE(esiragusa): this must work.
//    SEQAN_ASSERT(isSuffixArray(this->fibre, this->text));
}

// --------------------------------------------------------------------------
// Test getValue()
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Test begin() and end()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(CSATest, BeginEnd)
{
    typedef typename TestFixture::TFibre            TSA;
    typedef typename Iterator<TSA, Standard>::Type  TIter;
    typedef typename Position<TIter>::Type          TPos;
    typedef typename Difference<TIter>::Type        TDiff;

    TIter itBeg = begin(this->fibre, Standard());
    TIter itEnd = end(this->fibre, Standard());

    SEQAN_ASSERT_EQ(itEnd - itBeg, static_cast<TDiff>(length(this->fibre)));
    SEQAN_ASSERT_EQ(position(itBeg), static_cast<TPos>(0));
    SEQAN_ASSERT_EQ(position(itEnd), static_cast<TPos>(length(this->fibre)));
}

// ========================================================================== 
// Functions
// ========================================================================== 

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
