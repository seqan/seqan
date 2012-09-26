// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_CHAIN_GENERIC_H
#define SEQAN_HEADER_CHAIN_GENERIC_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TWeight>
struct ChainGenericEntry_
{
	TPos me;		//position of fragment here (within Source)
	TPos pre;		//position of precursor or -1 for top (within Frags)
	TWeight weight;	//weight for best chain until here
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSource>
struct ChainGenericSortFragsPredFunctional_
{
	TSource & source;
	ChainGenericSortFragsPredFunctional_(TSource & src)
		: source(src)
	{
	}
	ChainGenericSortFragsPredFunctional_(ChainGenericSortFragsPredFunctional_ const & other)
		: source(other.source)
	{
	}
	inline ChainGenericSortFragsPredFunctional_ &
	operator = (ChainGenericSortFragsPredFunctional_ const & other)
	{
		source = other.source;
		return *this;
	}
	~ChainGenericSortFragsPredFunctional_()
	{
	}

	template <typename TFrags>
	inline bool 
	operator() (TFrags const & left, TFrags const & right) const
	{
		return leftPosition(source[left.me], 0) < leftPosition(source[right.me], 0);
	}
};

//____________________________________________________________________________


template <typename TSource, typename TFrags, typename TScoring>
inline void
_chainGenericInitFrags(TSource & source,
						 TFrags & frags,
						 TScoring scoring)
{
	SEQAN_ASSERT_GT(length(source), 0u);

	typedef typename Position<TSource>::Type TPos;
	typedef typename Value<TSource>::Type TFragment;
	typedef typename Value<TFrags>::Type TFrag;

	//create top fragment
	TFragment top(dimension(source[0]));
	makeBeginFragment(top);

	//create entry in frags for each item in source
	resize(frags, length(source));
	for (TPos i = beginPosition(source); i < endPosition(source); ++i)
	{
		TFrag & frag = frags[i];
		frag.me = i;
		frag.pre = ~0UL; //link it with the top fragment
		frag.weight = scoreChainGap(scoring, top, source[i]) + weight(source[i]);
	}

	std::sort(begin(frags, Standard()), end(frags, Standard()), ChainGenericSortFragsPredFunctional_<TSource>(source));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFrag>
inline bool
_chainGenericChainable(TFrag & f1,
						 TFrag & f2)
{
	SEQAN_ASSERT_EQ(dimension(f1), dimension(f2));

	unsigned int dim = dimension(f1); 
	while (dim > 0)
	{
		--dim;
		if (rightPosition(f1, dim) > leftPosition(f2, dim))
		{
			return false;
		}
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TFrags, typename TIterator, typename TScoring>
inline void
_chainGenericFindBest(TSource & source,
						TFrags & frags,
						TIterator & it_act,
						TScoring scoring)
{
	typedef typename Iterator<TFrags, Standard>::Type TFragsIterator;
	typedef typename Value<TScoring>::Type TWeight;
	typedef typename Reference<TIterator>::Type TFragRef;

	TFragRef act = *it_act;
	TWeight act_weight = weight(source[act.me]);
	TFragsIterator it_begin = begin(frags, Standard());
	for (TFragsIterator it = it_begin; it < it_act; ++it)
	{
		TFragRef frag = *it;
		if (_chainGenericChainable(source[frag.me], source[act.me]))
		{
			TWeight score = frag.weight + scoreChainGap(scoring, source[frag.me], source[act.me]) + act_weight;
			if (score > act.weight)
			{//better predecessor found
				act.pre = it - it_begin;
				act.weight = score;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TFrags, typename TIterator, typename TDest, typename TFragment>
inline void
_chainGenericTraceback(TSource & source,
						 TFrags & frags,
						 TIterator & it_best,
						 TDest & dest,
						 TFragment & bottom)
{
	typedef typename Position<TSource>::Type TPos;

	clear(dest);
	//build chain in reverse order
	appendValue(dest, bottom); //chain will end with bottom fragment

	appendValue(dest, source[(*it_best).me]);
	for (TPos pos = (*it_best).pre; pos != ~0UL; pos = frags[pos].pre)
	{
		appendValue(dest, source[frags[pos].me]);
	}

	//chain will start with top fragment
	TFragment top(dimension(bottom));
	makeBeginFragment(top);
	appendValue(dest, top);

	//reverse chain
	std::reverse(begin(dest, Standard()), end(dest, Standard()));
}


//////////////////////////////////////////////////////////////////////////////

//spec for GenericChaining
template< typename TSource, typename TDest, typename TScoring>
inline typename Value<TScoring>::Type
globalChaining(TSource & source, 
	  TDest & dest, 
	  TScoring const & scoring, 
	  GenericChaining)
{
	typedef typename Value<TSource>::Type TFragment;
	typedef typename Weight<TFragment>::Type TWeight;
	typedef typename Position<TSource>::Type TSourcePosition;
	typedef ChainGenericEntry_<TSourcePosition, TWeight> TFrag;
	typedef String<TFrag> TFrags;
	typedef typename Iterator<TFrags, Standard>::Type TFragsIterator;

	//initialize fragments
	TFrags frags;
	_chainGenericInitFrags(source, frags, scoring);

	TFragsIterator it_begin = begin(frags, Standard());
	TFragsIterator it_end = end(frags, Standard());
	TFragsIterator it_best = it_begin;
	TWeight weight_best = MinValue<TWeight>::VALUE;


	//create bottom fragment
	unsigned int dim = dimension(source[0]);
	TFragment bottom(dim);
	makeEndFragment(bottom, source);

	//iterate all fragments
	for (TFragsIterator it = it_begin; it < it_end; ++it)
	{
		//find best predecessor for *it
		_chainGenericFindBest(source, frags, it, scoring);

		//determine heaviest fragment
		TWeight weight_it = (*it).weight + scoreChainGap(scoring, source[(*it).me], bottom);

		if (weight_it > weight_best)
		{
			it_best = it;
			weight_best = weight_it;
		}
	}

	//follow best fragment back to the beginning of the chain
	_chainGenericTraceback(source, frags, it_best, dest, bottom);

	return weight_best;
}



//____________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
