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
// Author: Christopher Pockrandt <christopher.pockrandt@fu-berlin.de>
// ==========================================================================

#include <seqan/index.h>
#include <ctime>

#include "test_index_helpers.h"

using namespace seqan;

template <typename TSpec = void, typename TLengthSum = size_t>
struct FMIndexConfigLevelsPrefix
{
    typedef TLengthSum                                                      LengthSum;
    typedef Levels<TSpec, LevelsPrefixRDConfig<TLengthSum, Alloc<>, 1, 0> > Bwt;
    typedef Levels<TSpec, LevelsRDConfig<TLengthSum, Alloc<>, 1, 0> >       Sentinels;

    static const unsigned SAMPLING =                                        10;
};

template <typename TSpec = void, typename TLengthSum = size_t>
struct FMIndexWTConfig
{
    typedef TLengthSum                                                      LengthSum;
    typedef WaveletTree<TSpec, WTRDConfig<TLengthSum, Alloc<>, 1, 0> >      Bwt;
    typedef Levels<TSpec, LevelsRDConfig<TLengthSum, Alloc<>, 1, 0> >       Sentinels;

    static const unsigned SAMPLING =                                        10;
};

typedef
    TagList<Index<String<bool>,   BidirectionalIndex<FMIndex<void, FMIndexConfigLevelsPrefix<> > > >,
    TagList<Index<DnaString,      BidirectionalIndex<FMIndex<void, FMIndexConfigLevelsPrefix<> > > >,
    TagList<Index<Dna5String,     BidirectionalIndex<FMIndex<void, FMIndexConfigLevelsPrefix<> > > >,
    TagList<Index<String<bool>,   BidirectionalIndex<FMIndex<void, FMIndexWTConfig<> > > >,
    TagList<Index<DnaString,      BidirectionalIndex<FMIndex<void, FMIndexWTConfig<> > > >,
    TagList<Index<Dna5String,     BidirectionalIndex<FMIndex<void, FMIndexWTConfig<> > > >
    > > > > > >
    FMIndices;

// ==========================================================================
// Test Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class RankDictionaryTest
// --------------------------------------------------------------------------

template <typename TIndex_>
class BidirectionalFMIndexTest : public Test
{
public:
    typedef TIndex_ TIndex;
};

SEQAN_TYPED_TEST_CASE(BidirectionalFMIndexTest, FMIndices);

// testing the bidirectional FM index by comparing ranges and hits against two stand-alone
// FM indices of the original and the reversed text
template <typename TBiFMIndex, typename TText, typename TPattern>
inline bool
testBidirectionalIndex(TBiFMIndex & bifmIndex, TText & text, TText & revText, TPattern & pattern)
{
    typedef Index<TText, FMIndex<> >                           TFMIndex;
    typedef typename Iterator<TFMIndex, TopDown<> >::Type      TFMIter;
    typedef typename Iterator<TBiFMIndex, TopDown<> >::Type    TBiFMIter;

    ModifiedString<TPattern, ModReverse> revPattern(pattern);

    TFMIndex indexFwd(text);
    TFMIndex indexRev(revText);
    TFMIter itFwd(indexFwd);
    TFMIter itRev(indexRev);

    TBiFMIter bifm(bifmIndex);

    bool res1 = goDown(itFwd, revPattern);
    bool res2 = goDown(itRev, pattern);

    std::mt19937 rng(time(nullptr));
    unsigned left = rng() % length(pattern);
    unsigned right = left;

    bool res3 = goDown(bifm, pattern[left]);
    while (res3 && (0 < left || right < length(pattern) - 1))
    {
        if (rng() % 2 && 0 < left)
        {
            --left;
            res3 = goDown(bifm, pattern[left], Fwd());
        }
        else if (right < length(pattern) - 1)
        {
            ++right;
            res3 = goDown(bifm, pattern[right], Rev());
        }
    }

    SEQAN_ASSERT_EQ(res1, res2);
    SEQAN_ASSERT_EQ(res1, res3);

    if (res1) // if pattern was found in index
    {
        SEQAN_ASSERT(getOccurrences(itFwd) == getOccurrences(bifm, Fwd()));
        SEQAN_ASSERT(getOccurrences(itRev) == getOccurrences(bifm, Rev()));
    }

    return 0;
}

SEQAN_TYPED_TEST(BidirectionalFMIndexTest, SearchInString)
{
    typedef typename TestFixture::TIndex                        TIndex;
    typedef typename Host<TIndex>::Type                         TText;

    std::mt19937 rng(time(nullptr));

    TText text;
    generateText(rng, text, 3947);
    TText revText(text);
    reverse(revText);

    TIndex index(text);

    for (unsigned patternLength = 1; patternLength <= 10; ++patternLength)
    {
        TText pattern;
        generateText(rng, pattern, patternLength);

        testBidirectionalIndex(index, text, revText, pattern);
    }
}

#ifndef __alpha__ // NOTE(h-2): fails on alpha for unknown reasons
SEQAN_TYPED_TEST(BidirectionalFMIndexTest, SearchInStringSet)
{
    typedef typename TestFixture::TIndex                        TIndex;
    typedef typename Host<TIndex>::Type                         TText;
    typedef typename Spec<TIndex>::Type                         TIndexSpec;
    typedef StringSet<TText, Owner<ConcatDirect<void> > >       TStringSet;
    typedef Index<TStringSet, TIndexSpec>                       TStringSetIndex;

    std::mt19937 rng(time(nullptr));

    unsigned textLength = 3947;

    TStringSet stringSet;
    TStringSet revStringSet;
    for (unsigned stringSetSize = 1; stringSetSize <= 3; ++stringSetSize)
    {
        TText text;
        generateText(rng, text, textLength);
        ModifiedString<TText, ModReverse> revText(text);

        appendValue(stringSet, text);
        appendValue(revStringSet, revText);

        TStringSetIndex index(stringSet);
        for (unsigned patternLength = 1; patternLength <= 20; ++patternLength)
        {
            TText pattern;
            if (rng() % 2) // guaranteed hit
                pattern = infixWithLength(text, rng() % (textLength - patternLength), patternLength);
            else // likely to have no hits (for longer pattern and short texts)
                generateText(rng, pattern, patternLength);

            testBidirectionalIndex(index, stringSet, revStringSet, pattern);
        }
    }
}
#endif // __alpha__

// ==========================================================================
// Functions
// ==========================================================================

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
