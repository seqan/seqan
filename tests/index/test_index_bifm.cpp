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

// testing the bidirectional FM index by comparing ranges and hits against two stand-alone
// FM indices of the original and the reversed text
template <typename TBiFMIndex, typename TText, typename TPattern>
inline bool
testBidirectionalIndex(TText & text, TPattern & pattern)
{
    typedef Index<TText, FMIndex<> >                           TFMIndex;
    typedef typename Iterator<TFMIndex, TopDown<> >::Type      TFMIter;
    typedef typename Iterator<TBiFMIndex, TopDown<> >::Type    TBiFMIter;

    TText revText(text);
    reverse(revText);

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
    bool res3 = goDown(bifm1, pattern, Rev());
    bool res4 = goDown(bifm2, revPattern, Fwd());

    SEQAN_ASSERT_EQ(res1, res2);
    SEQAN_ASSERT_EQ(res1, res3);
    SEQAN_ASSERT_EQ(res1, res4);

    if (res1) // if pattern was found in string
    {
        getOccurrences(bifm2, Rev());
        SEQAN_ASSERT(getOccurrences(itFwd) == getOccurrences(bifm1, Fwd()));
        SEQAN_ASSERT(getOccurrences(itFwd) == getOccurrences(bifm2, Fwd()));
        SEQAN_ASSERT(getOccurrences(itRev) == getOccurrences(_iter(bifm1, Rev())));
        SEQAN_ASSERT(getOccurrences(itRev) == getOccurrences(bifm2, Rev()));
    }

    return 0;
}

SEQAN_DEFINE_TEST(bifm_index_iterator_range_check)
{
    using namespace seqan;

    typedef DnaString                                           TText;
    typedef StringSet<TText, Owner<ConcatDirect<void> > >       TStringSet;

    typedef Index<TText, BidirectionalIndex<FMIndex<> > >       TIndex;
    typedef Index<TStringSet, BidirectionalIndex<FMIndex<> > >  TStringSetIndex;

    const unsigned int textLength = 10;

    for (unsigned int patternLength = 0; patternLength <= textLength; ++patternLength)
    {

        TText text, pattern;
        generateText(text, textLength);
        generateText(pattern, patternLength);

        testBidirectionalIndex<TIndex>(text, pattern);

        TStringSet stringSet;
        for (unsigned int stringSetSize = 1; stringSetSize <= 3; ++stringSetSize)
        {
            generateText(text, textLength);
            appendValue(stringSet, text);

            testBidirectionalIndex<TStringSetIndex>(stringSet, pattern);
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
