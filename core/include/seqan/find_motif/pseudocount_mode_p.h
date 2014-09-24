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

#ifndef SEQAN_HEADER_PSEUDOCOUNT_MODE_P_H
#define SEQAN_HEADER_PSEUDOCOUNT_MODE_P_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// PMode
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.PMode:
..summary: Represents the P computation scheme for handling "zero" probabilities.
..general:Class.Pseudocount
..cat:Motif Search
..signature:Pseudocount<TValue, PMode>
..param.TValue:The type of sequence which is considered.
...type:Spec.Dna
...type:Spec.AminoAcid
..remarks:The P mode computation scheme distributes the pseudocounts among the various residue
          according to their background probabilities.
..include:seqan/find_motif.h
*/

///.Class.Pseudocount.param.TSpec.type:Spec.PMode

struct PMode_;
typedef Tag<PMode_> PMode;


template<typename TValue>
class Pseudocount<TValue, PMode>
{
	enum { SIZE = ValueSize<TValue>::VALUE };

//_________________________________________________________________________________________________

public:
	double pseudocounts[SIZE];
	double epsilon;

//_________________________________________________________________________________________________

	template<typename TFrequencyDistribution>
	Pseudocount(double epsilon_, TFrequencyDistribution & background_):
		epsilon(epsilon_)
	{
		_computePseudocount(background_); 
	}
	Pseudocount(Pseudocount const & other_):
		pseudocounts(other_.pseudocounts),
		epsilon(other_.epsilon)
	{
	}
	~Pseudocount()
	{
	}

	Pseudocount const &
	operator = (Pseudocount const & other_)
	{
		this->pseudocount = other_.pseudocounts;
		this->epsilon = other_.epsilon;

		return *this;
	}

//_________________________________________________________________________________________________

private:

	template<typename TFrequencyDistribution>
	void
	_computePseudocount(TFrequencyDistribution & background_frequency) 
	{
		typedef typename Position< TFrequencyDistribution >::Type TPos;
		
		// calculating pseudocounts for each residue i.
		// pseudocount = epsilon*fi
		for(TPos i=0; i<length(background_frequency); ++i)
		{
			pseudocounts[i] = (double)(epsilon*background_frequency[i]);
		}
	}
//_________________________________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

// Function.normalize (s. profile.h)

template<typename TProfile, typename TValue>
void 
normalize(TProfile & profile, Pseudocount<TValue, PMode> & mode)
{
	typedef typename Value<TProfile>::Type TFreqDist;
	typedef typename Spec<TFreqDist>::Type TFrequencyType;

	typename Size<TProfile>::Type profile_size = length(profile);
	for(typename Position<TProfile>::Type i=0; 
		i<profile_size; 
		++i)
	{
		if(std::find(begin(profile[i]), end(profile[i]), 0)!=end(profile[i]))
		{
			// N:=row sum
			TFrequencyType N = sum(profile[i]);

			// add pseudocounts
			for(typename Position<TFreqDist>::Type j=0; 
				j<length(profile[i]); 
				++j)
			{
				profile[i][j] = 
					((TFrequencyType)(profile[i][j]+mode.pseudocounts[j]))/
					((TFrequencyType)(N+mode.epsilon));
			}
		}
		else
		{
			normalize(profile[i]);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
