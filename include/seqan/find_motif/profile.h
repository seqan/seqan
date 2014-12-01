// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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

#ifndef SEQAN_HEADER_PROFILE_H
#define SEQAN_HEADER_PROFILE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template<typename TProfile, typename TIterStr>
void 
convertPatternToProfile(TProfile & profile,
						TIterStr str_begin,
						TIterStr str_end)
{
//IOREV _notio_
	typedef typename Position<TProfile>::Type TPos;
	unsigned int str_size = str_end-str_begin;
	resize(profile, str_size);
	TPos pos = 0;
	while(str_begin!=str_end)
	{
		convertResidueToFrequencyDist(profile[pos], *str_begin);
		++str_begin;
		++pos;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TProfile, typename TStrings, typename TPseudocountMode>
void
convertSetOfPatternsToProfile(
        TProfile & profile,
        TStrings const & l_mers, 
        TPseudocountMode & pseudocount)
{
//IOREV _notio_
	typedef typename Value<TStrings const>::Type TString;
	//typedef typename Value<TProfile>::Type TFrequencyDistribution;

	typename Size<TString>::Type l = length(l_mers[0]);
	resize(profile, l);

	typename Iterator<TStrings const>::Type l_mers_iter, l_mers_end;
	typename Iterator<TString const>::Type l_mer_iter, l_mer_end;
	l_mers_iter = begin(l_mers);
	l_mers_end = end(l_mers);
	while(l_mers_iter!=l_mers_end)
	{
		l_mer_iter = begin(*l_mers_iter);
		l_mer_end = end(*l_mers_iter);
		while(l_mer_iter!=l_mer_end)
		{
			++profile[(int)(l_mer_iter-begin(*l_mers_iter))][(int)*l_mer_iter];
			++l_mer_iter;
		}
		++l_mers_iter;
	}
	normalize(profile, pseudocount);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TProfile>
void 
normalize(TProfile & profile)
{
	typename Iterator<TProfile>::Type iter = begin(profile);
	typename Iterator<TProfile>::Type iter_end = end(profile);
	while(iter!=iter_end)
	{
		normalize(*iter);
		++iter;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TProfile>
void 
completeProfile(TProfile & profile,
				typename Value<TProfile>::Type & background_distribution)
{
//IOREV _notio_
	TProfile copy(profile);
	resize(profile, length(copy)+1);
	profile[0] = background_distribution;

	typename Iterator<TProfile>::Type iter = begin(copy);
	int counter = 1;
	for(; !atEnd(iter, copy); goNext(iter))
	{
		profile[counter] = *iter;
		++counter;
	}
}

//////////////////////////////////////////////////////////////////////////////


template<typename TStrings>
void 
display(TStrings & strings)
{
	if(length(strings)!=0)
	{
		typename Iterator<TStrings>::Type iter = begin(strings);
		int counter = 0;
		for(; !atEnd(iter, strings); goNext(iter))
		{
			std::cout << "[" << counter << "]: " << *iter << "\n";
			++counter;
		}
		std::cout << "\n";
	}
	else
	{
		std::cout << "EMPTY STRINGS !!!\n";
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
