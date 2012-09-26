// ==========================================================================
//                                  FMIndex
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <your.email@example.net>
// ==========================================================================

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_PREFIX_SUM_TABLE_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_PREFIX_SUM_TABLE_H_

#include <iostream>
#include <algorithm>

namespace seqan {



template <typename TChar>
struct PST
{
	typedef unsigned TPos;

	String<Pair<TChar, TPos> > entries;
	PST() :
		entries()
	{}

	PST(String<TChar> & input) :
		entries()
	{
		Pair<TChar, unsigned> temp(input[0], 1);
		appendValue(entries, temp);
		for(unsigned i = 1; i < length(input); ++i)
		{
			//Pair<TChar, unsigned> * entrie_;
			unsigned pos;
			if(binarySearch_(*this, input[i], pos))
				++(entries[pos].i2);
			else
			{
				Pair<TChar, unsigned> temp(input[i], 1);
				appendValue(entries, temp);
				std::sort(begin(entries), end(entries));
			}
		}
	}
};

template <typename TChar>
struct Iterator<PST<TChar> const, Standard>
{
    typedef Iter<String<Pair<TChar, unsigned> > const, PositionIterator> Type;
};

template <typename TChar>
struct Iterator<PST<TChar>, Standard>
{
	typedef Iter<String<Pair<TChar, unsigned> >, PositionIterator> Type;
};

template <typename TChar>
struct Iterator<PST<TChar> const, Rooted>:
	Iterator<PST<TChar> const, Standard>{};

template <typename TChar>
struct Iterator<PST<TChar>, Rooted>:
	Iterator<PST<TChar>, Standard>{};

template <typename TChar>
bool binarySearch_(PST<TChar> & pst, TChar & character, unsigned & pos)
{
	Pair<TChar, unsigned> character_(character, 0);
	typename Iterator<String<Pair<TChar, unsigned> >, Rooted>::Type low, up;
	low = std::lower_bound(begin(pst.entries, Rooted()), end(pst.entries, Rooted()), character_);
	up = std::upper_bound(begin(pst.entries, Rooted()), end(pst.entries, Rooted()), character_);
	if(getValue(low).i1 != character)
	{
		return false;
		pos = position(up);
	}
	pos = position(low);
	return true;

//	Pair<TChar, unsigned> low(pst.entries[0].i1,0);
//	Pair<TChar, unsigned> up(pst.entries[length(pst.entries) -1].i1, length(pst.entries) -1);
//	Pair<TChar, unsigned> middle(pst.entries[low.i2 + (up.i2 - low.i2) / 2].i1, low.i2 + (up.i2 - low.i2) / 2);
//
//	do
//	{
//		std::cerr << "low: " << low << " middle: " << middle << " up: " << up << std::endl;
//		if(middle.i1 <= character)
//			low = middle;
//		else
//			up = middle;
//		middle = Pair<TChar, unsigned>(pst.entries[low.i2 + (up.i2 - low.i2) / 2].i1, low.i2 + (up.i2 - low.i2) / 2);
//	} while(low.i2 < up.i2);
//
//	if(low.i1 != character)
//		return false;
//	pos = low.i2;

}

}


#endif // SANDBOX_MY_SANDBOX_APPS_FMINDEX_PREFIX_SUM_TABLE_H_
