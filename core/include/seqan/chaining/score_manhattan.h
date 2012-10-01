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


#ifndef SEQAN_HEADER_SCORE_MANHATTAN_H
#define SEQAN_HEADER_SCORE_MANHATTAN_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Score spec for manhattan distance

/**.Spec.Score Manhattan
..summary:Scoring scheme for chaining that computes gap scores using manhattan distance.
..cat:Chaining
..general:Class.Score
..signature:Score<TValue, Manhattan>
..param.TValue:Type of the score values.
..remarks:The manhattan distance between two n-dimensional points is defined is the sum of the (absolute) differences of their coordinates. 
*/
template <typename TValue>
class Score<TValue, Manhattan>
{
public:
	TValue data_match;
	TValue data_mismatch;
	TValue data_gap;

public:
	Score( TValue _match = 0, TValue _misalign = 1 ):
		data_match( _match ),
		data_mismatch( _misalign ),
		data_gap( _misalign )
	{
	}

	Score( TValue score ):
		data_match( 0 ),
		data_mismatch( score ),
		data_gap( score )
	{
	}

	Score(Score const & other):
		data_match(other.data_match),
		data_mismatch(other.data_mismatch),
		data_gap(other.data_gap)
	{
	}

	~Score()
	{
	}

	Score & operator = (Score const & other)
	{
		data_match = other.data_match;
		data_mismatch = other.data_mismatch;
		data_gap = other.data_gap;
		return *this;
	}



//____________________________________________________________________________
};


//____________________________________________________________________________

template <typename TValue>
inline TValue 
scoreMatch(Score<TValue, Manhattan> & me)
{
	return me.data_match;
}
template <typename TValue>
inline TValue
scoreMatch(Score<TValue, Manhattan> const & me)
{
	return me.data_match;
}

template <typename TValue>
inline TValue 
scoreMismatch(Score<TValue, Manhattan> & me)
{
	return me.data_mismatch;
}
template <typename TValue>
inline TValue
scoreMismatch(Score<TValue, Manhattan> const & me)
{
	return me.data_mismatch;
}

template <typename TValue>
inline TValue 
scoreGap(Score<TValue, Manhattan> & me)
{
	return me.data_gap;
}
template <typename TValue>
inline TValue const &
scoreGap(Score<TValue, Manhattan> const & me)
{
	return me.data_gap;
}
//////////////////////////////////////////////////////////////////////////////

//Shortcut:


//template< typename TValue >
//typedef typename Score< TValue, Manhattan > ManhattanScore;

//////////////////////////////////////////////////////////////////////////////

// compute score of aligning two characters

//template <typename TValue, typename T>
//inline TValue
//score(Score<TValue, Manhattan> const & me,
//	  T const & left,
//	  T const & right)
//{
//	if (left == right) return scoreMatch(me);
//	else return scoreMismatch(me);
//}

//////////////////////////////////////////////////////////////////////////////
//compute score for chaining two fragments 
//return value this is only valid for f1 < f2, 
//that is f2 can be appended to f1

template <typename TValue, typename TFragment>
inline TValue
scoreChainGap(Score<TValue, Manhattan> const & me,
			  TFragment & f1,
			  TFragment & f2)
{
	SEQAN_ASSERT_EQ(dimension(f1), dimension(f2));

	unsigned int dim = dimension(f1);
	TValue score = 0;
	TValue score_gap = scoreGap(me);
	for (unsigned int i = 0; i < dim; ++i)
	{
		score -= score_gap * (leftPosition(f2, i) - rightPosition(f1, i));
	}
	return score;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
