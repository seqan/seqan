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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/index.h>

#include "test_index_helpers.h"

using namespace seqan;

template <typename TBiFMIndex>
inline bool
testBidirectionalIndex(unsigned int textLength, unsigned int patternLength)
{
	typedef DnaString											TText;

	typedef Index<TText, FMIndex<> >							TFMIndex;
	typedef typename Iterator<TFMIndex, TopDown<> >::Type		TFMIter;
	typedef typename Iterator<TBiFMIndex, TopDown<> >::Type		TBiFMIter;

	TText text, pattern;
	generateText<TText>(text, textLength);
	generateText<TText>(pattern, patternLength);

	TText revText(text);
	reverse(revText);

	TFMIndex fmIndex1(text);
	TFMIndex fmIndex2(revText);
	TFMIter fm1(fmIndex1);
	TFMIter fm2(fmIndex2);

	TBiFMIndex bifmIndex(text);

	TBiFMIter bifm1(bifmIndex);
	TBiFMIter bifm2(bifmIndex);

	bool res1 = true, res2 = true, res3 = true, res4 = true;

	unsigned int i;
	for (i = 0; i < patternLength; ++i)
	{
		res1 &= goDown(fm1, pattern[patternLength - i - 1]);
		res2 &= goDown(fm2, pattern[i]);
		res3 &= rightExtend(bifm1, pattern[i]);
		res4 &= leftExtend(bifm2, pattern[patternLength - i - 1]);
	}

	if (i == patternLength)
	{
	    SEQAN_ASSERT_EQ(res1, res2);
	    SEQAN_ASSERT_EQ(res1, res3);
	    SEQAN_ASSERT_EQ(res1, res4);

	    if (res1) // if all are true
	    {
	    	SEQAN_ASSERT_EQ(fm1.vDesc.range, bifm1.fwdIter.vDesc.range);
	    	SEQAN_ASSERT_EQ(fm1.vDesc.range, bifm2.fwdIter.vDesc.range);
	    	SEQAN_ASSERT_EQ(fm2.vDesc.range, bifm1.bwdIter.vDesc.range);
	    	SEQAN_ASSERT_EQ(fm2.vDesc.range, bifm2.bwdIter.vDesc.range);
	    }
	}

	return 0;
}

template <typename TBiFMIndex>
inline bool
testBidirectionalIndex(unsigned int textLength, unsigned int patternLength, unsigned int seqNumb)
{
	typedef DnaString											TString;
	typedef StringSet<TString, Owner<ConcatDirect<void> > >		TText;

	typedef Index<TText, FMIndex<> >							TFMIndex;
	typedef typename Iterator<TFMIndex, TopDown<> >::Type		TFMIter;
	typedef typename Iterator<TBiFMIndex, TopDown<> >::Type		TBiFMIter;

	TText text, revText;
	TString pattern;

	unsigned int i;
	for (i = 0; i < seqNumb; ++i) {
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

	generateText<TString>(pattern, patternLength);

	TFMIndex fmIndex1(text);
	TFMIndex fmIndex2(revText);
	TFMIter fm1(fmIndex1);
	TFMIter fm2(fmIndex2);

	TBiFMIndex bifmIndex(text);

	TBiFMIter bifm1(bifmIndex);
	TBiFMIter bifm2(bifmIndex);

	bool res1 = true, res2 = true, res3 = true, res4 = true;

	for (i = 0; i < patternLength; ++i)
	{
		res1 &= goDown(fm1, pattern[patternLength - i - 1]);
		res2 &= goDown(fm2, pattern[i]);
		res3 &= rightExtend(bifm1, pattern[i]);
		res4 &= leftExtend(bifm2, pattern[patternLength - i - 1]);
	}

	if (i == patternLength)
	{
	    SEQAN_ASSERT_EQ(res1, res2);
	    SEQAN_ASSERT_EQ(res1, res3);
	    SEQAN_ASSERT_EQ(res1, res4);

	    if (res1) // if all are true
	    {
	    	SEQAN_ASSERT_EQ(fm1.vDesc.range, bifm1.fwdIter.vDesc.range);
	    	SEQAN_ASSERT_EQ(fm1.vDesc.range, bifm2.fwdIter.vDesc.range);
	    	SEQAN_ASSERT_EQ(fm2.vDesc.range, bifm1.bwdIter.vDesc.range);
	    	SEQAN_ASSERT_EQ(fm2.vDesc.range, bifm2.bwdIter.vDesc.range);
	    }
	}

	return 0;
}

SEQAN_DEFINE_TEST(bifm_index_iterator_range_check)
{
    using namespace seqan;
    {
    	typedef DnaString TText;
    	typedef Index<TText, BidirectionalFMIndex<> > TIndex;
    	typedef Index<StringSet<TText, Owner<ConcatDirect<void> > >, BidirectionalFMIndex<> >  TStringSetIndex;
		for (int textLength = 0; textLength < 15; ++textLength)
		{
			for (int patternLength = 0; patternLength <= textLength; ++patternLength)
			{
				testBidirectionalIndex<TIndex>(textLength, patternLength);
				for (int stringSetSize = 1; stringSetSize <= 15; ++stringSetSize)
				{
					testBidirectionalIndex<TStringSetIndex>(textLength, patternLength, stringSetSize);
				}
			}
		}
    }
    {
    	testBidirectionalIndex<Index<DnaString, BidirectionalFMIndex<> > >(10000, 100);
    	testBidirectionalIndex<Index<StringSet<DnaString, Owner<ConcatDirect<void> > >, BidirectionalFMIndex<> > >(10000, 100, 10);
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
