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

/**
.Class.MotifFinder:
..cat:Motif Search
..summary:Holds the algorithm parameter values and the motif instance(s) found by the appropriate
          motif discovery algorithm.
..cat:Motif Search
..signature:MotifFinder<TValue, TSpec, TRng>
..param.TValue:The type of sequences to be analyzed.
...metafunction:Metafunction.Value
...type:Spec.Dna
...type:Spec.AminoAcid
..param.TSpec:The motif finding algorithm to search with.
...type:Spec.Projection
...type:Spec.EPatternBranching
...type:Spec.Pms1
...type:Spec.Pmsp
..param.TRng:The @Class.Rng@ specialization to use for random number generation.
...default:$GetDefaultRng<MotifFinderClass>::Type$
..include:seqan/find_motif.h
*/

struct MotifFinderClass_;
typedef Tag<MotifFinderClass_> MotifFinderClass;
    
template <typename TValue, typename TSpec, typename TRng = typename GetDefaultRng<MotifFinderClass>::Type>
class MotifFinder;

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.MotifFinder
///.Metafunction.Value.class:Class.MotifFinder

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

/**
.Function.findMotif
..class:Class.MotifFinder
..summary:Represents the main function which is used to start the search for noticeable motif patterns.
..cat:Motif Search
..signature:findMotif(finder,dataset,seq_model)
..param.finder:The @Class.MotifFinder@ object.
...type:Class.MotifFinder
..param.dataset:The dataset object representing the input sequences.
...type:Class.StringSet
..param.seq_model:The seq_model object.
...type:Tag.Oops
...type:Tag.Omops
...type:Tag.Zoops
...type:Tag.Tcm
...remarks:The sequence models rely on different assumptions about the distribution of motif occurrences
           across the sample sequences. 
..remarks:The PROJECTION algorithm is able to run in @Tag.Oops@, @Tag.Zoops@ and @Tag.Tcm@ mode.
..remarks:The ePatternBranching algorithm is able to run in @Tag.Oops@ and @Tag.Omops@ mode.
..remarks:The Pms1 and Pmsp algorithm is able to run in  @Tag.Oops@,  @Tag.Omops@,  @Tag.Zoops@ and  
          @Tag.Tcm@ mode.
..include:seqan/find_motif.h
*/

/**
.Function.factorial:
..summary:Calculates the factorial value of any integer number.
..cat:Motif Search
..signature:factorial(value)
..param.value:The value object.
...remarks:$value$ must be a positive integer.
..remarks:The factorial of a non-negative integer $value$ is 
          the product of all positive integers less than or equal to $value$.  
..include:seqan/find_motif.h
*/

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

/**
.Function.binomialCoefficient:
..summary:Calculates the binomial coefficient C(n,k).
..cat:Motif Search
..signature:binomialCoefficient(n,k)
..param.n:The first parameter object.
...remarks:$n$ must be a positive integer.
..param.k:The second parameter object.
...remarks:$k$ must be a positive integer.
..remarks:The binomial coefficient of $n$ and $k$ is equal to zero 
          if $k$ is greater than $n$.   
..include:seqan/find_motif.h
*/

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

/**
.Function.hammingDistance:
..summary:Determines the Hamming distance between two sequences.
..cat:Motif Search
..signature:hammingDistance<TType>(begin1,end1,begin2)
..param.TType:Distance type.
..param.begin1:An iterator pointing to the beginning of the first sequence which is either
              a @Shortcut.DnaString@ or a @Shortcut.Peptide@. 
...type:Concept.RandomAccessIteratorConcept
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..param.end1:An iterator pointing to the end of the first sequence which is either
            a @Shortcut.DnaString@ or a @Shortcut.Peptide@. 
...type:Concept.RandomAccessIteratorConcept
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..param.begin2:An iterator pointing to the beginning of the second sequence which is either
              a @Shortcut.DnaString@ or a @Shortcut.Peptide@. 
...type:Concept.RandomAccessIteratorConcept
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..include:seqan/find_motif.h
*/

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

/**
.Function.inverseHash:
..summary:Determines the corresponding sequence pattern given the hash value.
..cat:Motif Search
..signature:inverseHash<TValue>(hash_value,alphabet_size,seq_size)
..param.hash_value:The hash_value object.
..param.alphabet_size:The alphabet_size object.
...remarks:$alphabet_size$ is four for nucleotide sequences and twenty for amino acid sequences.
..param.seq_size:The seq_size object representing the size of the corresponding sequence.
..include:seqan/find_motif.h
*/

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

/**
.Function.displayResult:
..class:Class.MotifFinder
..summary:Displays all found motif candidates. In the case of the Projection Motif Finder
          the function displays the consensus pattern of the found motif candidate.
..cat:Motif Search
..signature:displayResult(motif_finder)
..param.motif_finder:The @Class.MotifFinder@ object.
...type:Class.MotifFinder
..include:seqan/find_motif.h
*/

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
/**
.Metafunction.Motif:
..cat:Motif Search
..summary:The string type of the finder.
..signature:Motif<T>::Type
..param.T:Finder for which the string type is determined.
...type:Class.String
..returns.param.Type:Underlying sequence type of finder $T$.
..include:seqan/find_motif.h
 */
template <typename T>
struct Motif;

template <typename TValue, typename TSpec, typename TRng>
struct Motif< MotifFinder<TValue, TSpec, TRng> >
{
	typedef String<TValue> Type;
};

/////////////////////////////////////////////////////////////////////////
/**
.Function.getMotif:
..class:Class.MotifFinder
..summary:Gets the motif out of a @Class.MotifFinder@.  If pos is given, the pos-th motif is returned, otherwise the first motif is returned.
..cat:Motif Search
..signature:getMotif(motifFinder, pos)
..param.motifFinder:
...type:Class.MotifFinder
..param.pos:Position 
..include:seqan/find_motif.h
*/

template <typename TValue, typename TSpec, typename TPosition, typename TRng>
inline typename Motif<MotifFinder<TValue, TSpec, TRng> >::Type &
getMotif(MotifFinder<TValue, TSpec, TRng> & me,
		 TPosition pos)
{
    SEQAN_CHECKPOINT;
	return me.set_of_motifs[pos];
}

///.Function.getMotif.signature:getMotif(motifFinder)
template <typename TValue, typename TSpec, typename TRng>
inline typename Motif<MotifFinder<TValue, TSpec, TRng> >::Type &
getMotif(MotifFinder<TValue, TSpec, TRng> & me)
{
    SEQAN_CHECKPOINT;
	return me.set_of_motifs[0];
}

/////////////////////////////////////////////////////////////////////////
/**
.Function.motifCount:
..class:Class.MotifFinder
..summary:Gets number of motifs in the @Class.MotifFinder@.
..cat:Motif Search
..signature:motifCount(motifFinder)
..param.motifFinder:
...type:Class.MotifFinder
..include:seqan/find_motif.h
*/

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
