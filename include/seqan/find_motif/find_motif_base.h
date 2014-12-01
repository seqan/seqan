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

#ifndef SEQAN_HEADER_FIND_MOTIF_BASE_H
#define SEQAN_HEADER_FIND_MOTIF_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

///////////////////////////////////////////////////////////////////////////////////////////////////

struct MotifFinderClass_;
typedef Tag<MotifFinderClass_> MotifFinderClass;
    
template <typename TValue, typename TSpec, typename TRng = typename GetDefaultRng<MotifFinderClass>::Type>
class MotifFinder;

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TRng>
struct Value< MotifFinder<TValue, TSpec, TRng> >
{
	typedef TValue Type;
};
template<typename TValue, typename TSpec, typename TRng>
struct Value< MotifFinder<TValue, TSpec, TRng> const>
{
	typedef TValue const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////


template<typename TType>
TType factorial(TType n)
{
    SEQAN_CHECKPOINT;

	TType result = 0;

	if(n==0)
	{
		result = 1;
	}
	else
	{
		result = n*factorial(n-1);
	}
   
	return result;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TType>
TType binomialCoefficient(TType n, TType k)
{
    SEQAN_CHECKPOINT;

	//SEQAN_ASSERT(!(n<0) & !(k<0));
	TType result = 1;
	for(TType i=(n-k+1); i<=n; ++i)
	{
		result*=i;
	}
	result = result/factorial(k);
	
	return result;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TType, typename TStringIterator>
TType hammingDistance(TStringIterator start1, TStringIterator end1, TStringIterator start2)
{
    SEQAN_CHECKPOINT;

	TType num_of_mismatches = 0;
	while(start1!=end1)
	{
		if(*start1!=*start2)
		{
			++num_of_mismatches;
		}
		++start1;
		++start2;
	}

	return num_of_mismatches;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TType>
String<TValue>
inverseHash(TType const & hash_value, 
			typename Size<TValue>::Type const & alp_size, 
			typename Size< String<TValue> >::Type const & seq_size)
{
    SEQAN_CHECKPOINT;

	typedef String<TValue> TString;
	TString seq;
	resize(seq, seq_size);

	TType hash_val = hash_value;
	typedef typename Position<TString>::Type TPos;
	for(TPos i=0; i<seq_size; ++i)
	{
		int letter = hash_val%alp_size;
		seq[i] = (TValue)letter;
		hash_val = (hash_val-letter)/alp_size;
	}

	std::reverse(begin(seq), end(seq));
	return seq;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TAlgorithm, typename TRng>
void
displayResult(MotifFinder<TValue, TAlgorithm, TRng> & finder)
{
    SEQAN_CHECKPOINT;

	typedef String<TValue> TString;
	typedef String<TString> TStrings;

	if(length(finder.set_of_motifs)!=0)
	{
		unsigned int counter = 0;
		typename Iterator<TStrings>::Type iter = begin(finder.set_of_motifs);
		for(; !atEnd(iter, finder.set_of_motifs); goNext(iter))
		{
			std::cout << "[" << counter << "]: " << *iter << "\n";
			++counter;
		}
		std::cout << "\n";
	}
	else
	{
		std::cout << "NO MOTIF HAS BEEN FOUND!!!\n";
	}
}

/////////////////////////////////////////////////////////////////////////
template <typename T>
struct Motif;

template <typename TValue, typename TSpec, typename TRng>
struct Motif< MotifFinder<TValue, TSpec, TRng> >
{
	typedef String<TValue> Type;
};

/////////////////////////////////////////////////////////////////////////
template <typename TValue, typename TSpec, typename TPosition, typename TRng>
inline typename Motif<MotifFinder<TValue, TSpec, TRng> >::Type &
getMotif(MotifFinder<TValue, TSpec, TRng> & me,
		 TPosition pos)
{
    SEQAN_CHECKPOINT;
	return me.set_of_motifs[pos];
}

template <typename TValue, typename TSpec, typename TRng>
inline typename Motif<MotifFinder<TValue, TSpec, TRng> >::Type &
getMotif(MotifFinder<TValue, TSpec, TRng> & me)
{
    SEQAN_CHECKPOINT;
	return me.set_of_motifs[0];
}

/////////////////////////////////////////////////////////////////////////
template <typename TValue, typename TSpec, typename TRng>
inline size_t
motifCount(MotifFinder<TValue, TSpec, TRng> const & me)
{
    SEQAN_CHECKPOINT;
	return length(me.set_of_motifs);
}


/////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
