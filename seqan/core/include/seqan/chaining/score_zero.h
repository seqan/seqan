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


#ifndef SEQAN_HEADER_SCORE_ZERO_H
#define SEQAN_HEADER_SCORE_ZERO_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
/**.Spec.Score Zero
..summary:Scoring scheme for chaining that set gap scores to 0
..cat:Chaining
..general:Class.Score
..signature:Score<TValue, Zero>
..param.TValue:Type of the score values.
*/

template <typename TValue>
class Score<TValue, Zero>
{
private:

public:
	Score()
	{
	}

	Score(Score const & )
	{
	}

	~Score()
	{
	}

	Score & operator = (Score const & other)
	{
		if( this == & other )
			return *this;
		return *this;
	}

//____________________________________________________________________________
};

template <typename TValue>
inline TValue 
scoreMatch(Score<TValue, Zero> &)
{
	return 0;
}
template <typename TValue>
inline TValue const 
scoreMatch(Score<TValue, Zero> const &)
{
	return 0;
}

template <typename TValue>
inline TValue 
scoreMismatch(Score<TValue, Zero> & /*me*/)
{
	return 0;
}
template <typename TValue>
inline TValue const 
scoreMismatch(Score<TValue, Zero> const &)
{
	return 0;
}

template <typename TValue>
inline TValue 
scoreGap(Score<TValue, Zero> &)
{
	return 0;
}
template <typename TValue>
inline TValue const 
scoreGap(Score<TValue, Zero> const &)
{
	return 0;
}
	
//////////////////////////////////////////////////////////////////////////////

//Shortcut:

//template< typename TValue >
//typedef typename Score<TValue, Zero> ZeroScore;

//////////////////////////////////////////////////////////////////////////////

// compute score of aligning two characters

template <typename TValue, typename T>
inline TValue
score(Score<TValue, Zero> const & /*me*/,
	  T const & /*left*/,
	  T const & /*right*/)
{
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
//compute score for chaining two fragments

template <typename TValue, typename TFragment>
inline TValue
scoreChainGap(Score<TValue, Zero> const &,
			  TFragment &,
			  TFragment &)
{
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
