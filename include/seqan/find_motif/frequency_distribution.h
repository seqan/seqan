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

#ifndef SEQAN_HEADER_FREQUENCY_DISTRIBUTION_H
#define SEQAN_HEADER_FREQUENCY_DISTRIBUTION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec = double>
class FrequencyDistribution
{
//____________________________________________________________________________________________

	enum { SIZE = ValueSize<TValue>::VALUE }; 

//____________________________________________________________________________________________

public:
	String<TSpec> frequency_list;

	// constructor & destructor
	FrequencyDistribution()
	{
        resize(frequency_list, static_cast<unsigned>(SIZE), 0);
	}
//	FrequencyDistribution(TValue const & letter_)
//	{
//		resize(frequency_list, (unsigned int) SIZE);
//		convertResidueToFrequencyDist(*this, letter_);
//	}
//	FrequencyDistribution(FrequencyDistribution const & other_)
//	{
//		frequency_list = other_.frequency_list; 
//	}
//	~FrequencyDistribution()
//	{
//	}

//	// overloading operators
//	FrequencyDistribution & 
//	operator = (FrequencyDistribution const & other_)
//	{
//		if(this!=&other_)
//		{
//			clear(frequency_list);
//			frequency_list = other_.frequency_list; 
//		}
//		return *this;
//	}

	FrequencyDistribution &
	operator += (FrequencyDistribution const & other_)
	{
		for(unsigned int i=0; i<SIZE; ++i)
		{
			frequency_list[i]+=other_.frequency_list[i];
		}

		return *this;
	}



	FrequencyDistribution &
	operator -= (FrequencyDistribution const & other_)
	{
		for(unsigned int i=0; i<SIZE; ++i)
		{
			frequency_list[i]-=other_.frequency_list[i];
		}

		return *this;
	}

	template<typename TType>
	FrequencyDistribution &
	operator *= (TType value_)
	{
		for(unsigned int i=0; i<SIZE; ++i)
		{
			frequency_list[i]*= (TSpec)value_;
		}

		return *this;
	}

	template<typename TPos>
	inline TSpec &
	operator [] (TPos index_)
	{
		return frequency_list[index_];
	}

	template<typename TPosition>
	inline TSpec const & 
	operator [] (TPosition const index_) const
	{
		return frequency_list[index_];
	}


//____________________________________________________________________________________________

};

template <typename TValue, typename TSpec>
FrequencyDistribution<TValue, TSpec> 
operator + (FrequencyDistribution<TValue, TSpec> const & lhs_, FrequencyDistribution<TValue, TSpec> const & rhs_)
{
	FrequencyDistribution<TValue, TSpec> ret(lhs_);
	ret+=rhs_;

	return ret;
}


template <typename TValue, typename TSpec>
FrequencyDistribution<TValue, TSpec> 
operator - (FrequencyDistribution<TValue, TSpec> const & lhs_, FrequencyDistribution<TValue, TSpec> const & rhs_)
{
	FrequencyDistribution<TValue, TSpec> ret(lhs_);
	ret-=rhs_;

	return ret;
}

template <typename TValue, typename TSpec, typename TType>
FrequencyDistribution<TValue, TSpec> 
operator * (FrequencyDistribution<TValue, TSpec> const & fd_, TType value_)
{
	FrequencyDistribution<TValue, TSpec> ret(fd_);
	ret*=value_;

	return ret;
}

template <typename TTarget, typename TValue, typename TSpec>
inline void
write(TTarget &target, FrequencyDistribution<TValue, TSpec> & fd_)
{
    write(target, fd_.frequency_list);
}

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator<<(TStream & target,
           FrequencyDistribution<TValue, TSpec> const & source)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, source);
    return target;
}

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.FrequencyDistribution

template <typename TValue, typename TSpec>
struct Spec< FrequencyDistribution<TValue, TSpec> >
{
	typedef TSpec Type;
};

template <typename TValue, typename TSpec>
struct Spec< FrequencyDistribution<TValue, TSpec> const>
{
	typedef TSpec Type;
};

//the following is a workaround for an error in GCC 4.1.2
template <typename TValue>
struct Spec< FrequencyDistribution<TValue, double> >
{
	typedef double Type;
};
template <typename TValue>
struct Spec< FrequencyDistribution<TValue, double> const>
{
	typedef double Type;
};


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.FrequencyDistribution

template <typename TValue, typename TSpec, typename TIteratorSpec>
struct Iterator< FrequencyDistribution<TValue, TSpec>, TIteratorSpec >
{
	typedef String<TSpec> TString_;
	typedef typename Iterator<TString_, TIteratorSpec>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Position.param.T.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
struct Position< FrequencyDistribution<TValue, TSpec> >
{
	typedef String<TSpec> TString_;
	typedef typename Position<TString_>::Type Type;
};
template<typename TValue, typename TSpec>
struct Position< FrequencyDistribution<TValue, TSpec> const>
{
	typedef String<TSpec> TString_;
	typedef typename Position<TString_ const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
struct Size< FrequencyDistribution<TValue, TSpec> >
{
	typedef String<TSpec> TString_;
	typedef typename Size<TString_>::Type Type;
};
template<typename TValue, typename TSpec>
struct Size<FrequencyDistribution<TValue, TSpec> const>
{
	typedef String<TSpec> TString_;
	typedef typename Size<TString_ const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.FrequencyDistribution
/*
.Metafunction.Value:
..summary:Returns the sequence type of a @Class.FrequencyDistribution@ type
	     (TValue for FrequencyDistribution<TValue, TSpec>).
..include:seqan/find_motif.h
*/

template<typename TValue, typename TSpec>
struct Value< FrequencyDistribution<TValue, TSpec> >
{
	typedef TValue Type;
};
template<typename TValue, typename TSpec>
struct Value< FrequencyDistribution<TValue, TSpec> const>
{
	typedef TValue const Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSeqIter> 
void 
absFreqOfLettersInSeq(FrequencyDistribution<TValue, TSpec> & fd,
					  TSeqIter seq_start,
					  TSeqIter seq_end) 
{	
	while(seq_start!=seq_end)
	{
		++fd[(int)*seq_start];
		++seq_start;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TIter>
void
absFreqOfLettersInSetOfSeqs(FrequencyDistribution<TValue, TSpec> & fd,
							TIter seq_start,
							TIter seq_end)
{
	while(seq_start!=seq_end)
	{
		absFreqOfLettersInSeq(fd, begin(*seq_start), end(*seq_start));
		++seq_start;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TType>
void 
addValue(FrequencyDistribution<TValue, TSpec> & fd, TType const & val)
{
	typedef typename Position< FrequencyDistribution<TValue, TSpec> >::Type TPos;
	for(TPos i=0; i<length(fd); ++i)
	{
		fd[i]+= (TSpec)val;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec,typename TDatasetIter> 
void 
backgroundFrequency(FrequencyDistribution<TValue, TSpec> & fd,
					TDatasetIter dataset_start,
					TDatasetIter dataset_end)
{
	absFreqOfLettersInSetOfSeqs(fd, dataset_start, dataset_end);

	// check for zero entries
	if(std::find(begin(fd), end(fd), (TSpec)0)!= end(fd))
	{
		// add pseudocounts
		double epsilon = 0.1;
		seqan::Pseudocount<TValue, CMode> p(epsilon);
		addValue(fd, p.pseudocount);
	}
	normalize(fd);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.begin.param.object.type:Class.FrequencyDistribution
///.Function.begin.class:Class.FrequencyDistribution

template <typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> >::Type
begin(FrequencyDistribution<TValue, TSpec> & me)
{
	return begin(me.frequency_list);
}
template <typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> >::Type
begin(FrequencyDistribution<TValue, TSpec> const & me)
{
	return begin(me.frequency_list);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.clear.param.object.type:Class.FrequencyDistribution
///.Function.clear.class:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
void 
clear(FrequencyDistribution<TValue, TSpec> & fd) 
{
	fd = FrequencyDistribution<TValue, TSpec>();
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec>
void 
convertResidueToFrequencyDist(FrequencyDistribution<TValue, TSpec> & fd, TValue const & residue)
{
	typedef typename Position< FrequencyDistribution<TValue, TSpec> >::Type TPos;
	TSpec probability = 
		(TSpec)(0.5/(ValueSize<TValue>::VALUE-1));

	for(TPos i=0; i<length(fd); ++i)
	{
		if(i==residue)
		{
			fd[i] = 0.5;
		}
		else
		{
			fd[i] = probability;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

///.Function.end.param.object.type:Class.FrequencyDistribution
///.Function.end.class:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> >::Type
end(FrequencyDistribution<TValue, TSpec> & me)
{
	return begin(me)+length(me);
}
template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> const >::Type
end(FrequencyDistribution<TValue, TSpec> const & me)
{
	return begin(me)+length(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.FrequencyDistribution
///.Function.length.class:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
inline typename Size< FrequencyDistribution<TValue, TSpec> >::Type
length(FrequencyDistribution<TValue, TSpec> & me)
{
	return length(me.frequency_list);
}

template<typename TValue, typename TSpec>
inline typename Size< FrequencyDistribution<TValue, TSpec> >::Type
length(FrequencyDistribution<TValue, TSpec> const & me)
{
	return length(me.frequency_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec>
void 
logarithmize(FrequencyDistribution<TValue, TSpec> & fd) 
{
	typedef FrequencyDistribution<TValue, TSpec> TFrequencyDistribution;
	typedef typename Position<TFrequencyDistribution>::Type TPos;
	
	for(TPos i=0; i<length(fd); ++i)
	{
		fd[i] = (TSpec)log(fd[i]);
	}
}

//////////////////////////////////////////////////////////////////////////////

/* s. normalize() (profile.h)
.Function.normalize:
..summary:Determines the normalized frequencies.
..cat:Motif Search
..signature:normalize(frequencies)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
..include:seqan/find_motif.h
*/

template<typename TValue, typename TSpec>
void 
normalize(FrequencyDistribution<TValue, TSpec> & fd)
{
	typedef FrequencyDistribution<TValue, TSpec> TFrequencyDistribution;
	typedef typename Position<TFrequencyDistribution>::Type TPos;
	
	TSpec amount = sum(fd);
	for(TPos i=0; i<length(fd); ++i)
	{
		fd[i] = (TSpec)(fd[i]/amount);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec>
typename Position< FrequencyDistribution<TValue, TSpec> >::Type
posOfMax(FrequencyDistribution<TValue, TSpec> & me)
{
	typedef FrequencyDistribution<TValue, TSpec> TFrequencyDistribution; 
	typedef typename Position<TFrequencyDistribution>::Type TPos;

	TPos position = 0;
	TSpec max_value = (TSpec)0;
	for(TPos i=0; i<length(me); ++i)
	{
		if(me[i]>max_value)
		{
			position = i;
			max_value = me[i];
		}
	}

	return position;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec>
TSpec 
sum(FrequencyDistribution<TValue, TSpec> & me)
{
	TSpec amount = 
		std::accumulate(begin(me), end(me), (TSpec)0);

	return amount;
}
								
//////////////////////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
