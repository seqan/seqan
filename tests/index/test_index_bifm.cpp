// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

#include <seqan/basic.h>
#include <seqan/index.h>

#include "test_index_helpers.h"

using namespace seqan;

// generating a random text and reversing it
template <typename TString>
inline void
generateText(StringSet<TString, Owner<ConcatDirect<void> > > &text,
             StringSet<TString, Owner<ConcatDirect<void> > > &revText,
             const unsigned int textLength,
             const unsigned int seqNumb)
{
    for (unsigned int i = 0; i < seqNumb; ++i) {
        TString _text, _revText;
        generateText<TString>(_text, textLength);
        appendValue(text, _text);

        resize(_revText, textLength);
        for (unsigned int j = 0; j < textLength; ++j)
        {
            _revText[j] = _text[textLength - j - 1];
        }
        appendValue(revText, _revText);
    }
}

// generating a random text and reversing it
template <typename TText>
inline void
generateText(TText &text, TText &revText, const unsigned int textLength, const unsigned int /*seqNumb*/)
{
    generateText<TText>(text, textLength);
    revText = text;
    reverse(revText);
}

// comparing two lists of occurrences
template <typename TIter1, typename TIter2>
inline void
compareOccurrences(TIter1 &iter1, TIter2 &iter2)
{
    typedef typename Container<TIter1>::Type  TIndex1;
    typedef typename Container<TIter2>::Type  TIndex2;
    typedef typename Infix<typename Fibre<TIndex1, FibreSA>::Type const>::Type  TOccs1;
    typedef typename Infix<typename Fibre<TIndex2, FibreSA>::Type const>::Type  TOccs2;
    TOccs1 occs1 = getOccurrences(iter1);
    TOccs2 occs2 = getOccurrences(iter2);

    SEQAN_ASSERT_EQ(length(occs1), length(occs2));
    for (unsigned int i = 0; i < length(occs1); ++i)
    {
        SEQAN_ASSERT_EQ(occs1[i], occs2[i]);
    }
}

// testing the bidirectional FM index by comparing ranges and hits against two stand-alone
// FM indices of the original and the reversed text
template <typename TBiFMIndex, typename TPattern>
inline bool
testBidirectionalIndex(const unsigned int textLength, const unsigned int patternLength, const unsigned int seqNumb = 1)
{
    typedef typename Fibre<TBiFMIndex, FibreText>::Type        TText;

    typedef Index<TText, FMIndex<> >                           TFMIndex;
    typedef typename Iterator<TFMIndex, TopDown<> >::Type      TFMIter;
    typedef typename Iterator<TBiFMIndex, TopDown<> >::Type    TBiFMIter;

    TText text, revText;
    generateText(text, revText, textLength, seqNumb);

    TPattern pattern;
    generateText<TPattern>(pattern, patternLength);
    ModifiedString<TPattern, ModReverse> revPattern(pattern);

    TFMIndex indexFwd(text);
    TFMIndex indexRev(revText);
    TFMIter itFwd(indexFwd);
    TFMIter itRev(indexRev);

    TBiFMIndex bifmIndex(text);
    TBiFMIter bifm1(bifmIndex);
    TBiFMIter bifm2(bifmIndex);

    bool res1 = goDown(itFwd, revPattern);
    bool res2 = goDown(itRev, pattern);
    bool res3 = goDown(bifm1.bwdIter, pattern);
    bool res4 = goDown(bifm2.fwdIter, revPattern);

    SEQAN_ASSERT_EQ(res1, res2);
    SEQAN_ASSERT_EQ(res1, res3);
    SEQAN_ASSERT_EQ(res1, res4);

    if (res1) // if all are true
    {
        SEQAN_ASSERT_EQ(countOccurrences(itFwd), countOccurrences(bifm1.fwdIter));
        SEQAN_ASSERT_EQ(countOccurrences(itFwd), countOccurrences(bifm2.fwdIter));
        SEQAN_ASSERT_EQ(countOccurrences(itRev), countOccurrences(bifm1.bwdIter));
        SEQAN_ASSERT_EQ(countOccurrences(itRev), countOccurrences(bifm2.bwdIter));

        compareOccurrences(itFwd, bifm1.fwdIter);
        compareOccurrences(itFwd, bifm2.fwdIter);
        compareOccurrences(itRev, bifm1.bwdIter);
        compareOccurrences(itRev, bifm2.bwdIter);
    }

    return 0;
}

SEQAN_DEFINE_TEST(bifm_index_iterator_range_check)
{
    using namespace seqan;

    typedef DnaString TText;
    typedef Index<TText, BidirectionalFMIndex<> > TIndex;
    typedef Index<StringSet<TText, Owner<ConcatDirect<void> > >, BidirectionalFMIndex<> >  TStringSetIndex;

    for (int textLength = 1; textLength < 15; ++textLength)
    {
        for (int patternLength = 0; patternLength <= textLength; ++patternLength)
        {
            testBidirectionalIndex<TIndex, TText>(textLength, patternLength);
            for (int stringSetSize = 1; stringSetSize <= 3; ++stringSetSize)
            {
                testBidirectionalIndex<TStringSetIndex, TText>(textLength, patternLength, stringSetSize);
            }
        }
    }
}

// ========================================================================== 
// Functions
// ========================================================================== 

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    SEQAN_CALL_TEST(bifm_index_iterator_range_check);
    return TestSystem::runAll();
}
