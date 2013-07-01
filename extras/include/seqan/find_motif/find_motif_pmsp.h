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

#ifndef SEQAN_HEADER_FIND_PMSP_H
#define SEQAN_HEADER_FIND_PMSM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Pmsp
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Pmsp:
..summary: Represents the Pmsp algorithm of Davila et al.
..description.note:There are various known problems with the motif finding in SeqAn. We plan to fix this in an upcoming release.
..general:Class.MotifFinder
..cat:Motif Search
..signature:MotifFinder<TValue, Pmsp, TRng>
..param.TValue:The type of sequences to be analyzed.
...type:Spec.Dna
...type:Spec.AminoAcid
..remarks:The @Spec.Pmsp|Pmsp algorithm@ is an improvement of the @Spec.Pms1|Pms1 algorithm@.
          It examines each possible l-mer of the first input sequence, explores its neighborhood
		  and finally checks whether an l-mer in the neighborhood is a motif instance.
..param.TRng:The @Class.Rng@ specialization to use for random number generation.
...default:$GetDefaultRng<MotifFinderClass>::Type$
..include:seqan/find_motif.h
*/

///.Class.MotifFinder.param.TSpec.type:Spec.Pmsp

struct Pmsp_;
typedef Tag<Pmsp_> const Pmsp;

//////////////////////////////////////////////////////////////////////////////
// MotifFinder - Pmsp Spec
//
// t:=dataset size (number of sequences)
// n:=average sequence size
// l:=motif size
// d:=number of substitutions
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TRng>
class MotifFinder<TValue, Pmsp, TRng>
{
//_________________________________________________________________________________

public:
	typedef unsigned int TSize;
	typedef String<TValue> TString;
	typedef String<TString> TStrings;
	
	TSize motif_size;
	TSize num_of_substitutions;
	bool has_exact_substitutions;
	TStrings set_of_motifs; // result set

//_________________________________________________________________________________

	MotifFinder()
	{
    SEQAN_CHECKPOINT;
	}
	MotifFinder(TSize const & l_, TSize const & d_, bool const & is_exact_):
		motif_size(l_),
		num_of_substitutions(d_),
		has_exact_substitutions(is_exact_)
	{
    SEQAN_CHECKPOINT;
	}
	MotifFinder(MotifFinder const & other_):
		motif_size(other_.motif_size),
		num_of_substitutions(other_.num_of_substitutions),
		has_exact_substitutions(other_.has_exact_substitutions)
	{
    SEQAN_CHECKPOINT;
	}
	~MotifFinder()
	{
    SEQAN_CHECKPOINT;
	}

	MotifFinder const &
	operator = (MotifFinder const & other_)
	{
    SEQAN_CHECKPOINT;
		if(this!=&other_)
		{
			motif_size = other_.motif_size;
			num_of_substitutions = other_.num_of_substitutions;
			has_exact_substitutions = other_.has_exact_substitutions;
		}

		return *this;
	}

//_________________________________________________________________________________

}; // class MotifFinder<TValue, Pmsp>

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template<typename TSeqType, typename TStrings, typename TModel, typename TRng>
inline void
findMotif(MotifFinder<TSeqType, Pmsp, TRng> & finder, 
		  TStrings & dataset, 
		  TModel seq_model)
{
    SEQAN_CHECKPOINT;
	pmsp(finder.set_of_motifs,
		 dataset, 
		 finder.motif_size, 
		 finder.num_of_substitutions, 
		 finder.has_exact_substitutions,
		 seq_model);
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function.pmsp:
..summary:Represents the Pmsp algorithm.
          For each l-mer x in the first sequence (s0) Pmsp generates the set of neighbors of x (V)
          and tries to guess if an l-mer y in that neighborhood is a motif by checking
          whether there are l-mers in the rest sequences (s1,s2,...,st-1) that are at/at most 
          distance d from it.
..cat:Motif Search
..signature:pms1(result_set,dataset,l,d,is_exact,h,model_type)
..param.result_set:The result_set object.
..param.dataset:The dataset object representing the input sequences.
...type:Class.String
...signature:String<TString>
...param.TString:A @Class.String@ type
....type:Class.String
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.is_exact:The size of Hamming distance.
...type:$bool$
..param.model_type:The model_type object.
...type:Tag.Oops
...type:Tag.Omops
...type:Tag.Zoops
...type:Tag.Tcm
..remarks:The Pmsp algorithm is able to run in Oops, Omops, Zoops and Tcm mode.
..remarks:The resulted motif candidates found by the Pmsp algorithm will be stored in the result_set object.
..include:seqan/find_motif.h
*/

//////////////////////////////////////////////////////////////////////////////
//	Oops model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TType>
void 
pmsp(TStrings & result,
	  TStrings & dataset,
	  TType const & l,
	  TType const & d, 
	  bool const & is_exact,
	  Oops const & /*model_type*/)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TString>::Type TPos;

	typename Size<TStrings>::Type t = length(dataset);
	typename Iterator<TStrings>::Type ds_iter1 = begin(dataset);

	typename Size<TString>::Type seq_len1 = length(*ds_iter1); 
	typename Iterator<TString>::Type seq_iter1 = begin(*ds_iter1);
	typename Iterator<TString>::Type seq_end1 = begin(*ds_iter1)+(seq_len1-l+1);
	Shape<TValue> shape(l);
	std::set<TString> result_set;
	// ----------------------------------------------------------------------------
	// STEP 1:
	// create set of d-variants (set V) for every l-mer x in the first sequence.
	// create set L[i] := set of l-mers in sequence s(i+1) 
	//                    having a hamming distance of <=2d to x, i=0,...,(t-2).
	// guess if an l-mer y in of neighborhood V is a motif by checking
	// whether there are l-mers in the rest sequences (s1,s2,...,st-1) that are 
	// at distance d from it.
	// ----------------------------------------------------------------------------
	std::vector<int> V;
	String< String<TPos> > L;
	//unsigned int count = 0;
	while(seq_iter1!=seq_end1)
	{
		//construct set V
		createDVariants(V, seq_iter1, l, d, is_exact, shape);
		//construct set L
		resize(L, t-1);
		for(typename Position< String< String<TPos> > >::Type i=0; i<length(L); ++i)
		{
			String<TPos> pos_ar;
			std::vector<TPos> pos_vect;

			typename Size<TString>::Type seq_len2 = length(dataset[i+1]);
			typename Iterator<TString>::Type seq_iter2 = begin(dataset[i+1]);
			typename Iterator<TString>::Type seq_end2 = begin(dataset[i+1])+seq_len2-l+1;
			int seq_pos = 0;
			while(seq_iter2!=seq_end2)
			{
				if( hammingDistance<TType>(seq_iter2, seq_iter2+l, seq_iter1)<=2*d )
				{
					pos_vect.push_back(seq_pos);
				}
				++seq_iter2;
				++seq_pos;
			}
			resize(pos_ar, pos_vect.size());
			std::copy(pos_vect.begin(), pos_vect.end(), begin(pos_ar));
			L[i] = pos_ar;
			pos_vect.clear();
		}
		std::vector<int>::iterator v_iter = V.begin();
		std::vector<int>::iterator v_end = V.end();
		while(v_iter!=v_end)
		{
			TString l_mer = inverseHash<TValue>(*v_iter,ValueSize<TValue>::VALUE,l); 
			if(hasExactOneOccurrence(begin(l_mer),begin(dataset[0]),end(dataset[0]),l,d,is_exact))
			{
				unsigned int number = 0;
				bool isMotif = true;
				for(typename Position< String< String<TPos> > >::Type i=0; i<length(L); ++i)
				{
					typename Iterator< String<TPos> >::Type l_iter = begin(L[i]);
					typename Iterator< String<TPos> >::Type l_end = end(L[i]);
					while((l_iter!=l_end) && (number<2))
					{
						TPos pos = *l_iter;
						if(is_exact && hammingDistance<TType>(begin(dataset[i+1])+pos, begin(dataset[i+1])+pos+l, 
							                          begin(l_mer))==d)
						{
							++number;
						}
						else if(!(is_exact) && 
							      hammingDistance<TType>(begin(dataset[i+1])+pos, begin(dataset[i+1])+pos+l,
								                        begin(l_mer))<=d)                  
						{
							++number;
						}
						++l_iter;
					}
					if(number!=1)
					{
						isMotif = false;
						i = length(L);
					}
					number = 0;
				}
				if(isMotif)
				{
					result_set.insert(l_mer);
				}
			}
			++v_iter;
		}
		V.clear();
		++seq_iter1;
		clear(L);
	}
	// ----------------------------------------------------------------------------
	// STEP 2:
	// collect motif candidates
	// ----------------------------------------------------------------------------
	resize(result, result_set.size());
	std::copy(result_set.begin(), result_set.end(), begin(result));
}

//////////////////////////////////////////////////////////////////////////////
//	Omops model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TType>
void 
pmsp(TStrings & result,
	  TStrings & dataset,
	  TType  const & l,
	  TType const & d, 
	  bool const & is_exact,
	  Omops const & /*model_type*/)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TString>::Type TPos;

	typename Size<TStrings>::Type t = length(dataset);
	typename Iterator<TStrings>::Type ds_iter1 = begin(dataset);

	typename Size<TString>::Type seq_len1 = length(*ds_iter1); 
	typename Iterator<TString>::Type seq_iter1 = begin(*ds_iter1);
	Shape<TValue> shape(l);
	std::set<TString> result_set;
	// ----------------------------------------------------------------------------
	// STEP 1:
	// create set of d-variants (set V) for every l-mer x in the first sequence.
	// create set L[i] := set of l-mers in sequence s(i+1) 
	//                    having a hamming distance of <=2d to x, i=0,...,(t-2).
	// guess if an l-mer y in of neighborhood V is a motif by checking
	// whether there are l-mers in the rest sequences (s1,s2,...,st-1) that are 
	// at distance d from it.
	// ----------------------------------------------------------------------------
	std::vector<int> V;
	String< String<TPos> > L;
	while(seq_iter1!=(begin(*ds_iter1)+seq_len1-l+1))
	{		
		//construct set V
		createDVariants(V, seq_iter1, l, d, is_exact, shape);	
		//construct set L
		resize(L, t-1);
		for(typename Position< String< String<TPos> > >::Type i=0; i<length(L); ++i)
		{
			String<TPos> pos_ar;
			std::vector<TPos> pos_vect;

			typename Size<TString>::Type seq_len2 = length(dataset[i+1]);
			typename Iterator<TString>::Type seq_iter2 = begin(dataset[i+1]);
			typename Iterator<TString>::Type seq_end2 = begin(dataset[i+1])+seq_len2-l+1;
			int seq_pos = 0;
			while(seq_iter2!=seq_end2)
			{
				if( hammingDistance<TType>(seq_iter2, seq_iter2+l, seq_iter1)<=2*d )
				{
					pos_vect.push_back(seq_pos);
				}
				++seq_iter2;
				++seq_pos;
			}
			resize(pos_ar, pos_vect.size());
			std::copy(pos_vect.begin(), pos_vect.end(), begin(pos_ar));
			L[i] = pos_ar;
			pos_vect.clear();
		}
		std::vector<int>::iterator v_iter = V.begin();
		std::vector<int>::iterator v_end = V.end();
		while(v_iter!=v_end)
		{
			TString l_mer = inverseHash<TValue>(*v_iter,ValueSize<TValue>::VALUE,l);
			unsigned int number = 0;
			bool isMotif = true;
			for(typename Position< String< String<TPos> > >::Type i=0; i<length(L); ++i)
			{
				typename Iterator< String<TPos> >::Type l_iter = begin(L[i]);
				typename Iterator< String<TPos> >::Type l_end = end(L[i]);
				while((l_iter!=l_end) && (number<1))
				{
					TPos pos = *l_iter;
					if(is_exact && hammingDistance<TType>(begin(dataset[i+1])+pos, begin(dataset[i+1])+pos+l, 
												  begin(l_mer))==d)
					{
						++number;
					}
					else if(!(is_exact) && 
							  hammingDistance<TType>(begin(dataset[i+1])+pos, begin(dataset[i+1])+pos+l,
													begin(l_mer))<=d)                              
					{
						++number;
					}
					++l_iter;
				}
				if(number!=1)
				{
					isMotif = false;
					i=length(L);
				}
				number = 0;
			}
			if(isMotif)
			{
				result_set.insert(l_mer);
			}

			++v_iter;
		}
		V.clear();
		++seq_iter1;
		clear(L);
	}
	// ----------------------------------------------------------------------------
	// STEP 2:
	// collect motif candidates
	// ----------------------------------------------------------------------------
	resize(result, result_set.size());
	std::copy(result_set.begin(), result_set.end(), begin(result));
}

//////////////////////////////////////////////////////////////////////////////
//	Zoops model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TType>
void 
pmsp(TStrings & result,
	 TStrings & dataset,
	 TType const & l,
	 TType const & d, 
	 bool const & is_exact,
	 Zoops const & model_type)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TStrings>::Type TPos1;
	typedef typename Position<TString>::Type TPos2;

	typename Size<TStrings>::Type t = length(dataset);
	typename Iterator<TStrings>::Type ds_iter1 = begin(dataset);

	std::set<TString> result_set;
	TPos1 seq_nr = 0;
	for(; !atEnd(ds_iter1, dataset); goNext(ds_iter1))
	{
		std::vector<TPos1> relevant_pos_vect;
		for(TPos1 i=0; i<t; ++i)
		{
			if(i!=seq_nr)
			{
				relevant_pos_vect.push_back(i);
			}
		}
		typename Size<TString>::Type seq_len1 = length(*ds_iter1);
		typename Iterator<TString>::Type seq_iter1 = begin(*ds_iter1);
		typename Iterator<TString>::Type seq_end1 = begin(*ds_iter1)+(seq_len1-l+1);
		Shape<TValue> shape(l);
		while(seq_iter1!=seq_end1)
		{
			//construct set V
			std::vector<int> V;
			createDVariants(V, seq_iter1, l, d, is_exact, shape);
			//construct set L
			String< String<TPos2> > L;
			resize(L, t-1);
			for(typename Position< String< String<TPos2> > >::Type i=0; i<length(L); ++i)
			{
				String<TPos2> pos_ar;
				std::vector<TPos2> pos_vect;
				typename Size<TString>::Type seq_len2 = length(dataset[relevant_pos_vect[i]]);
				typename Iterator<TString>::Type seq_iter2 = begin(dataset[relevant_pos_vect[i]]);
				typename Iterator<TString>::Type seq_end2 = begin(dataset[relevant_pos_vect[i]])+(seq_len2-l+1);
				int seq_pos = 0;
				while(seq_iter2!=seq_end2)
				{
					if( hammingDistance<TType>(seq_iter2, seq_iter2+l, seq_iter1)<=2*d )
					{
						pos_vect.push_back(seq_pos);
					}
					++seq_iter2;
					++seq_pos;
				}
				resize(pos_ar, pos_vect.size());
				std::copy(pos_vect.begin(), pos_vect.end(), begin(pos_ar));
				L[i] = pos_ar;
				pos_vect.clear();
				clear(pos_ar);
			}

			std::vector<int>::iterator V_iter = V.begin();
			std::vector<int>::iterator V_end = V.end();
			while(V_iter!=V_end)
			{
				TString l_mer = inverseHash<TValue>(*V_iter, ValueSize<TValue>::VALUE, l);
				typename std::set<TString>::iterator iter = result_set.find(l_mer);
				if( (iter==result_set.end()) &&
					hasExactOneOccurrence(begin(l_mer), 
					                     begin(dataset[seq_nr]), 
										 end(dataset[seq_nr]), l, d, is_exact) )
				{
					int lower_limit = (int) floor(t*(model_type.threshold)+0.5)-1;
					bool isMotif = true;
					unsigned int number = 0;
					for(typename Position< String< String<TPos2> > >::Type i=0; i<length(L); ++i)
					{
						typename Iterator< String<TPos2> >::Type L_iter = begin(L[i]);
						typename Iterator< String<TPos2> >::Type L_end = end(L[i]);
						while((L_iter!=L_end) && (number<2))
						{
							TPos2 pos = *L_iter;
							if(is_exact && 
								hammingDistance<TType>(begin(dataset[relevant_pos_vect[i]])+pos, 
												begin(dataset[relevant_pos_vect[i]])+pos+l, 
												begin(l_mer))==d)
							{
								++number;
							}
							else
							{
								if((!is_exact) && 
									hammingDistance<TType>(begin(dataset[relevant_pos_vect[i]])+pos, 
													begin(dataset[relevant_pos_vect[i]])+pos+l, 
													begin(l_mer))<=d)
								{
									++number;
								}}

							++L_iter;
						}
						if(number>1)
						{
							isMotif = false;
							i = length(L);
						}
						else if(number==1)
						{
							--lower_limit;
						}
						if(lower_limit<=0)
						{
							i = length(L);
						}
						number = 0;
					}
					if(isMotif && (lower_limit<=0))
					{
						result_set.insert(l_mer);
					}
				}

				++V_iter;
			}
			V.clear();
			clear(L);
			++seq_iter1;
		}
		relevant_pos_vect.clear();
		++seq_nr;
	}
	// ----------------------------------------------------------------------------
	// STEP 2:
	// collect motif candidates
	// ----------------------------------------------------------------------------
	resize(result, result_set.size());
	std::copy(result_set.begin(), result_set.end(), begin(result));
}

//////////////////////////////////////////////////////////////////////////////
//	Tcm model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TType>
void 
pmsp(TStrings & result,
	 TStrings & dataset,
	 TType const & l,
	 TType const & d, 
	 bool const & is_exact,
	 Tcm const & model_type)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TStrings>::Type TPos1;
	typedef typename Position<TString>::Type TPos2;

	typename Size<TStrings>::Type t = length(dataset);
	typename Iterator<TStrings>::Type ds_iter1 = begin(dataset);

	std::set<TString> result_set;
	TPos1 seq_nr = 0;
	for(; !atEnd(ds_iter1, dataset); goNext(ds_iter1))
	{
		std::vector<TPos1> relevant_pos_vect;
		for(TPos1 i=0; i<t; ++i)
		{
			if(i!=seq_nr)
			{
				relevant_pos_vect.push_back(i);
			}
		}
		typename Size<TString>::Type seq_len1 = length(*ds_iter1);
		typename Iterator<TString>::Type seq_iter1 = begin(*ds_iter1);
		typename Iterator<TString>::Type seq_end1 = begin(*ds_iter1)+(seq_len1-l+1);
		Shape<TValue> shape(l);
		while(seq_iter1!=seq_end1)
		{
			//construct set V
			std::vector<int> V;
			createDVariants(V, seq_iter1, l, d, is_exact, shape);
			//construct set L
			String< String<TPos2> > L;
			resize(L, t-1);
			for(typename Position< String< String<TPos2> > >::Type i=0; i<length(L); ++i)
			{
				String<TPos2> pos_ar;
				std::vector<TPos2> pos_vect;
				typename Size<TString>::Type seq_len2 = length(dataset[relevant_pos_vect[i]]);
				typename Iterator<TString>::Type seq_iter2 = begin(dataset[relevant_pos_vect[i]]);
				typename Iterator<TString>::Type seq_end2 = begin(dataset[relevant_pos_vect[i]])+(seq_len2-l+1);
				int seq_pos = 0;
				while(seq_iter2!=seq_end2)
				{
					if( hammingDistance<TType>(seq_iter2, seq_iter2+l, seq_iter1)<=2*d )
					{
						pos_vect.push_back(seq_pos);
					}
					++seq_iter2;
					++seq_pos;
				}
				resize(pos_ar, pos_vect.size());
				std::copy(pos_vect.begin(), pos_vect.end(), begin(pos_ar));
				L[i] = pos_ar;
				pos_vect.clear();
				clear(pos_ar);
			}

			std::vector<int>::iterator V_iter = V.begin();
			std::vector<int>::iterator V_end = V.end();
			while(V_iter!=V_end)
			{
				TString l_mer = inverseHash<TValue>(*V_iter, ValueSize<TValue>::VALUE, l);
				typename std::set<TString>::iterator iter = result_set.find(l_mer);
				if( iter==result_set.end() )
				{
					int lower_limit = (int) floor(t*(model_type.threshold)+0.5)-1;
					bool isMotif = false;
					unsigned int number = 0;
					for(typename Position< String< String<TPos2> > >::Type i=0; i<length(L); ++i)
					{
						typename Iterator< String<TPos2> >::Type L_iter = begin(L[i]);
						typename Iterator< String<TPos2> >::Type L_end = end(L[i]);
						while((L_iter!=L_end) && (number<1))
						{
							TPos2 pos = *L_iter;
							if((is_exact) && 
								(hammingDistance<TType>(begin(dataset[relevant_pos_vect[i]])+pos, 
												begin(dataset[relevant_pos_vect[i]])+pos+l, 
												begin(l_mer))==d))
							{
								++number;
							}
							else if((!is_exact) && 
								(hammingDistance<TType>(begin(dataset[relevant_pos_vect[i]])+pos, 
												begin(dataset[relevant_pos_vect[i]])+pos+l, 
												begin(l_mer))<=d))
							{
								++number;
							}

							++L_iter;
						}
						if(number==1)
						{
							--lower_limit;
						}

						if(lower_limit==0)
						{
							isMotif = true;
							i = length(L);
						}
						number = 0;
					}
					if(isMotif)
					{
						result_set.insert(l_mer);
					}
				}

				++V_iter;
			}
			V.clear();
			clear(L);
			++seq_iter1;
		}
		relevant_pos_vect.clear();
		++seq_nr;
	}
	// ----------------------------------------------------------------------------
	// STEP 2:
	// collect motif candidates
	// ----------------------------------------------------------------------------
	resize(result, result_set.size());
	std::copy(result_set.begin(), result_set.end(), begin(result));
}

//////////////////////////////////////////////////////////////////////////////
//Subfunctions
//////////////////////////////////////////////////////////////////////////////

/*
.Function.hasExactOneOccurrence:
..summary:Checks if a given l-mer occurs exactly once in a given sequence
..cat:Motif Search
..signature:hasExactOneOccurrence(l_mer_begin,seq_begin,seq_end,l,d,is_exact)
..param.l_mer_begin:An iterator pointing to the beginning of a given l-mer pattern.
...type:Concept.RandomAccessIteratorConcept
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.seq_begin:An iterator pointing to the beginning of a given sequence.
...type:Concept.RandomAccessIteratorConcept
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.seq_end:An iterator pointing to the end of a given sequence.
...type:Concept.RandomAccessIteratorConcept
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.is_exact:The size of Hamming distance
...type:$bool$
..include:seqan/find_motif.h
*/

template<typename TStringIter, typename TType>
bool
hasExactOneOccurrence(TStringIter l_mer_begin,
					 TStringIter seq_begin,
					 TStringIter seq_end,
					 TType const & l,
					 TType const & d,
					 bool const & is_exact)
{
    SEQAN_CHECKPOINT;

	bool result = false;
	TType counter = 0;
	while( (seq_begin!=(seq_end-l+1)) )
	{
		if(is_exact)
		{
			counter += (hammingDistance<TType>(seq_begin, seq_begin+l, l_mer_begin)==d) ? 1 : 0;
		}
		else
		{
			counter += (hammingDistance<TType>(seq_begin, seq_begin+l, l_mer_begin)<=d) ? 1 : 0;
		}		
		++seq_begin;
	}

	if(counter==1)
	{
		result = true;
	}
	return result;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
