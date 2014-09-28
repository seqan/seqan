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

#ifndef SEQAN_HEADER_FIND_MOTIF_PMS1_H
#define SEQAN_HEADER_FIND_MOTIF_PMS1_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Pms1
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Pms1:
..summary: Represents the Pms1 algorithm developed by Rajasekaran et al.
..description.note:There are various known problems with the motif finding in SeqAn. We plan to fix this in an upcoming release.
..general:Class.MotifFinder
..cat:Motif Search
..signature:MotifFinder<TValue, Pms1, TRng>
..param.TValue:The type of sequences to be analyzed.
...type:Spec.Dna
...type:Spec.AminoAcid
..remarks:The exact @Spec.Pms1|Pms1 algorithm@ (Planted Motif Search 1) searches in the space
          of possible motifs such as @Spec.EPatternBranching@. The procedure of the algorithm
		  is quite simple. For every l-mer in each input sequence the algorithm generates
		  all possible length-l patterns in the Hamming distance $d$-neighborhood of $x$.
		  The neighbor sets for each sequence are then intersected so that at the end of the process
		  we get a group of l-mers or a single l-mer that occur(s) in each input sequence with $d$
		  substitutions.
..param.TRng:The @Class.Rng@ specialization to use for random number generation.
...default:$GetDefaultRng<MotifFinderClass>::Type$
..include:seqan/find_motif.h
*/

///.Class.MotifFinder.param.TSpec.type:Spec.Pms1

struct Pms1_;
typedef Tag<Pms1_> const Pms1;

//////////////////////////////////////////////////////////////////////////////
// MotifFinder - Pms1 Spec
//
// t:=dataset size (number of sequences)
// n:=average sequence size
// l:=motif size
// d:=number of substitutions
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TRng>
class MotifFinder<TValue, Pms1, TRng>
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

}; // class MotifFinder<TValue, Pms1>

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template<typename TSeqType, typename TStrings, typename TModel, typename TRng>
inline void
findMotif(MotifFinder<TSeqType, Pms1, TRng> & finder, 
		  TStrings & dataset, 
		  TModel seq_model)
{
    SEQAN_CHECKPOINT;
	pms1(finder.set_of_motifs,
		 dataset, 
		 finder.motif_size, 
		 finder.num_of_substitutions, 
		 finder.has_exact_substitutions,
		 seq_model);
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function.pms1:
..summary:Represents the Pms1 algorithm.
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
..remarks:The Pms1 algorithm is able to run in Oops, Omops, Zoops and Tcm mode.
..remarks:The resulted motif candidates found by the Pms1 algorithm will be stored in the result_set object.
..include:seqan/find_motif.h
*/

//////////////////////////////////////////////////////////////////////////////
//	Oops model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TType>
void 
pms1(TStrings & result_set,
	 TStrings & dataset,
	 TType const & l,
	 TType const & d, 
	 bool const & is_exact,
	 Oops const & /*model_type*/)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	Shape<TValue> shape(l);
	// ----------------------------------------------------------------------------
	// STEP 1:
	// processing the first sequence.
	// ----------------------------------------------------------------------------
	typename Iterator<TStrings>::Type ds_iter = begin(dataset);
	TString l_mer;
	std::vector<int> var_of_l_mer, var_all;

	typename Size<TString>::Type seq_len = length(*ds_iter);
	typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
	typename Iterator<TString>::Type seq_end = begin(*ds_iter)+(seq_len-l+1);
	while(seq_iter!=seq_end)
	{
		createDVariants(var_of_l_mer, seq_iter, l, d, is_exact, shape);
		std::copy(var_of_l_mer.begin(), var_of_l_mer.end(), std::back_inserter(var_all));
		var_of_l_mer.clear();
		++seq_iter;
	}

	// sort var_all and filter d-variants that occur exactly once in the sequence
	std::sort(var_all.begin(), var_all.end());
	std::vector<int>::iterator vect_iter = var_all.begin();
	std::vector<int>::iterator vect_end = var_all.end();
	std::vector<int> var_unique;
	while(vect_iter!=vect_end)
	{
		std::vector<int>::iterator upper_iter =
			std::upper_bound(var_all.begin(),var_all.end(),*vect_iter);
		int count = 
			std::count(vect_iter, upper_iter, *vect_iter);
		if(count==1)
		{
			var_unique.push_back(*vect_iter);
		}
		vect_iter+=count;
	}
	var_all.clear();
	// ----------------------------------------------------------------------------
	// STEP 2:
	// processing the remaining input sequences
	// ----------------------------------------------------------------------------
	std::vector<int> var_unique2, result;
	++ds_iter;
	for(; !atEnd(ds_iter, dataset); goNext(ds_iter))
	{
		seq_len = length(*ds_iter);
		seq_iter = begin(*ds_iter);
		seq_end = begin(*ds_iter)+(seq_len-l+1);
		while(seq_iter!=seq_end)
		{
			createDVariants(var_of_l_mer, seq_iter, l, d, is_exact, shape);
			std::copy(var_of_l_mer.begin(), var_of_l_mer.end(), std::back_inserter(var_all));
			var_of_l_mer.clear();
			++seq_iter;
		}
		// sort var_all and filter d-variants that occur exactly once in the sequence
		std::sort(var_all.begin(), var_all.end());
		vect_iter = var_all.begin();
		vect_end = var_all.end();
		while(vect_iter!=vect_end)
		{
			std::vector<int>::iterator upper_iter =
				std::upper_bound(var_all.begin(),var_all.end(),*vect_iter);
			int count = 
				std::count_if(vect_iter, upper_iter, std::bind2nd(std::equal_to<int>(), *vect_iter));
			if(count==1)
			{
				var_unique2.push_back(*vect_iter);
			}
			vect_iter+=count;
		}
		var_all.clear();

		//create an insert_iterator for int-vector result
		std::insert_iterator< std::vector<int> > res_ins(result, result.begin());
		std::set_intersection(var_unique.begin(), var_unique.end(),
							  var_unique2.begin(), var_unique2.end(), res_ins);
		var_unique2.clear();
		var_unique.clear();
		var_unique.resize(result.size());
		std::copy(result.begin(),result.end(),var_unique.begin());
		result.clear();
	}
	// ----------------------------------------------------------------------------
	// STEP 3:
	// insert the relevant d-variants into result set
	// ----------------------------------------------------------------------------
	resize(result_set, var_unique.size());
	vect_iter = var_unique.begin();
	vect_end = var_unique.end();
	unsigned int variant_nr =0;
	while(vect_iter!=vect_end)
	{
		l_mer = inverseHash<TValue>(*vect_iter, ValueSize<TValue>::VALUE, l);
		result_set[variant_nr] = l_mer;
		++variant_nr;
		++vect_iter;
	}
}

//////////////////////////////////////////////////////////////////////////////
//	Omops model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TType>
void 
pms1(TStrings & result_set,
	 TStrings & dataset,
	 TType const & l,
	 TType const & d, 
	 bool const & is_exact,
	 Omops const & /*model_type*/)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	Shape<TValue> shape(l);
	// ----------------------------------------------------------------------------
	// STEP 1:
	// processing the first sequence.
	// ----------------------------------------------------------------------------
	typename Iterator<TStrings>::Type ds_iter = begin(dataset);
	TString l_mer;
	std::vector<int> var_of_l_mer, var_all;

	typename Size<TString>::Type seq_len = length(*ds_iter);
	typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
	typename Iterator<TString>::Type seq_end = begin(*ds_iter)+(seq_len-l+1);
	while(seq_iter!=seq_end)
	{
		createDVariants(var_of_l_mer, seq_iter, l, d, is_exact, shape);
		std::copy(var_of_l_mer.begin(), var_of_l_mer.end(), std::back_inserter(var_all));
		var_of_l_mer.clear();
		++seq_iter;
	}
	// sort var_all and filter d-variants that occur exactly once in the sequence
	std::sort(var_all.begin(), var_all.end());
	std::vector<int> var_unique;
	std::unique_copy(var_all.begin(), var_all.end(), std::back_inserter(var_unique));
	var_all.clear();
	// ----------------------------------------------------------------------------
	// STEP 2:
	// processing remaining input sequences
	// ----------------------------------------------------------------------------
	std::vector<int> var_unique2, result;
	++ds_iter;
	for(; !atEnd(ds_iter, dataset); goNext(ds_iter))
	{
		seq_len = length(*ds_iter);
		seq_iter = begin(*ds_iter);
		while(seq_iter!=begin(*ds_iter)+(seq_len-l+1))
		{
			createDVariants(var_of_l_mer, seq_iter, l, d, is_exact, shape);
			std::copy(var_of_l_mer.begin(), var_of_l_mer.end(), std::back_inserter(var_all));
			var_of_l_mer.clear();
			++seq_iter;
		}
		// sort var_all and eliminate duplicates
		std::sort(var_all.begin(), var_all.end());
		std::unique_copy(var_all.begin(), var_all.end(), std::back_inserter(var_unique2));
		var_all.clear();

		//create an insert_iterator for int-vector result
		std::insert_iterator< std::vector<int> > res_ins(result, result.begin());
		std::set_intersection(var_unique.begin(), var_unique.end(),
							  var_unique2.begin(), var_unique2.end(), res_ins);
		var_unique2.clear();
		var_unique.clear();
		var_unique.resize(result.size());
		std::copy(result.begin(),result.end(),var_unique.begin());
		result.clear();
	}
	// ----------------------------------------------------------------------------
	// STEP 3:
	// insert relevant d-variants into result set
	// ----------------------------------------------------------------------------
	resize(result_set, var_unique.size());
	typename std::vector<int>::iterator vect_iter = var_unique.begin();
	typename std::vector<int>::iterator vect_end = var_unique.end();
	unsigned int variant_nr =0;
	while(vect_iter!=vect_end)
	{
		l_mer = inverseHash<TValue>(*vect_iter, ValueSize<TValue>::VALUE, l);
		result_set[variant_nr] = l_mer;
		++variant_nr;
		++vect_iter;
	}
}

//////////////////////////////////////////////////////////////////////////////
//	Zoops model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TType>
void 
pms1(TStrings & result_set,
	 TStrings & dataset,
	 TType const & l,
	 TType const & d, 
	 bool const & is_exact,
	 Zoops const & model_type)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	Shape<TValue> shape(l);

	// ----------------------------------------------------------------------------
	// STEP 1:
	// building d-variants for all l-mers from the input sequences and storing their
	// corresponding hash-values (values of variants which occur exactly once 
	// in an input sequence) in an int-string 'var_ar'.
	// index_mat contains the variant's hash values for each sequence.
	// ----------------------------------------------------------------------------
	typename Size<TStrings>::Type t = length(dataset);
	typename Iterator<TStrings>::Type ds_iter = begin(dataset);
	TString l_mer;
	std::vector<int> var_of_l_mer, var_all, result_vect;
	for(; !atEnd(ds_iter, dataset); goNext(ds_iter))
	{
		typename Size<TString>::Type seq_len = length(*ds_iter);
		typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
		typename Iterator<TString>::Type seq_end = begin(*ds_iter)+(seq_len-l+1);
		while(seq_iter!=seq_end)
		{
			createDVariants(var_of_l_mer, seq_iter, l, d, is_exact, shape);
			std::copy(var_of_l_mer.begin(), var_of_l_mer.end(), std::back_inserter(var_all));
			var_of_l_mer.clear();
			++seq_iter;
		}
		// sort var_all and filter d-variants that occur exactly once in the sequence
		std::sort(var_all.begin(), var_all.end());
		std::vector<int>::iterator vect_iter = var_all.begin();
		std::vector<int>::iterator vect_end = var_all.end();
		std::vector<int> var_unique;
		while(vect_iter!=vect_end)
		{
			std::vector<int>::iterator upper_iter =
				std::upper_bound(var_all.begin(),var_all.end(),*vect_iter);
			int count = 
				std::count(vect_iter, upper_iter, *vect_iter);
			if(count==1)
			{
				var_unique.push_back(*vect_iter);
			}
			vect_iter+=count;
		}
		var_all.clear();

		//create an insert_iterator for int-vector var_all
		std::insert_iterator< std::vector<int> > res_ins(var_all, var_all.begin());
		std::merge(result_vect.begin(), result_vect.end(),
			       var_unique.begin(), var_unique.end(), res_ins);
		var_unique.clear();
		result_vect.clear();
		result_vect.resize(var_all.size());
		std::copy(var_all.begin(),var_all.end(),result_vect.begin());
		var_all.clear();
	}
	// ----------------------------------------------------------------------------
	// STEP 2:
	// count votes of l-mers (d-variants) and insert relevant d-variants into result set
	// ----------------------------------------------------------------------------
	int lower_limit = (int) floor(t*(model_type.threshold)+0.5);
	// count votes of relevant l-mers (d-variants)
	std::vector<int>::iterator iter = result_vect.begin();
	std::vector<int>::iterator end_iter = result_vect.end();
	unsigned int num_of_results = 0;
	while(iter!=end_iter)
	{
		std::vector<int>::iterator upper_iter =
			std::upper_bound(result_vect.begin(),result_vect.end(),*iter);
		int count = 
			std::count(iter, upper_iter, *iter);
		if( count>=lower_limit )
		{
			var_all.push_back(*iter);
			++num_of_results;
		}
		iter+=count;
	}
	result_vect.clear();

	resize(result_set, num_of_results);
	num_of_results = 0;
	iter = var_all.begin();
	end_iter = var_all.end();
	while(iter!=end_iter)
	{
		l_mer = inverseHash<TValue>(*iter, ValueSize<TValue>::VALUE, l);
		result_set[num_of_results] = l_mer;
		++num_of_results;
		++iter;
	}
	var_all.clear();
}

//////////////////////////////////////////////////////////////////////////////
//	Tcm model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TType>
void 
pms1(TStrings & result_set,
	 TStrings & dataset,
	 TType const & l,
	 TType const & d, 
	 bool const & is_exact,
	 Tcm const & model_type)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	Shape<TValue> shape(l);
	// ----------------------------------------------------------------------------
	// STEP 1:
	// building d-variants for all l-mers from the input sequences and storing their
	// corresponding hash-values (values of variants which occur exactly once 
	// in an input sequence) in an int-string 'var_ar'.
	// index_mat contains the variant's hash values for each sequence.
	// ----------------------------------------------------------------------------
	typename Size<TStrings>::Type t = length(dataset);
	typename Iterator<TStrings>::Type ds_iter = begin(dataset);
	TString l_mer;
	std::vector<int> var_of_l_mer, var_all, result_vect, unique_vect;
	for(; !atEnd(ds_iter, dataset); goNext(ds_iter))
	{
		typename Size<TString>::Type seq_len = length(*ds_iter);
		typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
		while(seq_iter!=begin(*ds_iter)+(seq_len-l+1))
		{
			createDVariants(var_of_l_mer, seq_iter, l, d, is_exact, shape);
			std::copy(var_of_l_mer.begin(), var_of_l_mer.end(), std::back_inserter(var_all));
			var_of_l_mer.clear();
			++seq_iter;
		}
	
		// sort var_all and filter d-variants that occur exactly once in the sequence
		std::sort(var_all.begin(), var_all.end());
		std::insert_iterator< std::vector<int> > res_ins1(unique_vect, unique_vect.begin());
		std::unique_copy(var_all.begin(), var_all.end(), res_ins1);
		var_all.clear();

		//create an insert_iterator for int-vector var_all
		std::insert_iterator< std::vector<int> > res_ins2(var_all, var_all.begin());
		std::merge(result_vect.begin(), result_vect.end(),
			       unique_vect.begin(), unique_vect.end(), res_ins2);
		unique_vect.clear();
		result_vect.clear();
		result_vect.resize(var_all.size());
		std::copy(var_all.begin(),var_all.end(),result_vect.begin());
		var_all.clear();
	}
	// ----------------------------------------------------------------------------
	// STEP 2:
	// count votes of l-mers (d-variants) and insert relevant d-variants into result set
	// ----------------------------------------------------------------------------
	int lower_limit = (int) floor(t*(model_type.threshold)+0.5);

	// count votes of relevant l-mers (d-variants)
	std::vector<int>::iterator iter = result_vect.begin();
	std::vector<int>::iterator end_iter = result_vect.end();
	unsigned int num_of_results = 0;
	while(iter!=end_iter)
	{
		std::vector<int>::iterator upper_iter =
			std::upper_bound(result_vect.begin(),result_vect.end(),*iter);
		int count = 
			std::count(iter, upper_iter, *iter);
		if( count>=lower_limit )
		{
			var_all.push_back(*iter);
			++num_of_results;
		}
		iter+=count;
	}
	result_vect.clear();

	resize(result_set, num_of_results);
	num_of_results = 0;
	iter = var_all.begin();
	end_iter = var_all.end();
	while(iter!=end_iter)
	{
		l_mer = inverseHash<TValue>(*iter, ValueSize<TValue>::VALUE, l);
		result_set[num_of_results] = l_mer;
		++num_of_results;
		++iter;
	}
	var_all.clear();
}

//////////////////////////////////////////////////////////////////////////////
//Subfunctions
//////////////////////////////////////////////////////////////////////////////

/*
.Function.createDVariants:
..summary:Creates the d-variants of a given l-mer and computes their hash-values which will be
          finally stored in array 'variants'.
..cat:Motif Search
..signature:createDVariants(variants,l_mer_begin,l,d,is_exact,shape)
..param.variants:The set of d-variants of a given l-mer.
...type:String<int>
...remarks:The string holds the hash-values of each pattern being in the d-vicinity of a given l-mer.
..param.l_mer_begin:An iterator pointing to the beginning of a given l-mer pattern.
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
..param.shape:The @Class.Shape@ object.
...type:Class.Shape
...signature:Shape<TValue, SimpleShape>
...remarks:Is used for computing the hash-value of each variant.
..remarks: #d-variants = 
                  - sum(i,0,d,binomialCoefficient(l,i)*(alp_size-1)^i), if is_exact=false
                  - binomialCoefficient(l.d)*(alp_size-1)^d),           else
..include:seqan/find_motif.h
*/

template<typename TIntVect, typename TIterString, typename TType, typename TShape>
void 
createDVariants(TIntVect & variants,
				TIterString l_mer_begin,
				TType const & l,
				TType const & d,
				bool is_exact,
				TShape & shape)
{
    SEQAN_CHECKPOINT;

	typedef String<int> TIntArray;
	String<TIntArray> bitsets;
	if(is_exact)
	{
		// get all variants of bitsets consisting of d zeros(pos of mismatches)
		// and (l-d) ones
		_getVariantsOfBitset(bitsets, l, d);

		typename Iterator< String<TIntArray> >::Type bitsets_iter = begin(bitsets);
		typename Iterator< String<TIntArray> >::Type bitsets_end = end(bitsets);
		while(bitsets_iter!=bitsets_end)
		{
			TIntVect hash_val_of_variants;
			_buildVariants(hash_val_of_variants, l_mer_begin, d, *bitsets_iter, shape);
			std::copy(hash_val_of_variants.begin(), hash_val_of_variants.end(), std::back_inserter(variants));
			++bitsets_iter;
		}
		clear(bitsets);
	}
	else
	{
		for(unsigned int i=0; i<=d; ++i)
		{
			_getVariantsOfBitset(bitsets, l, i);
			typename Iterator< String<TIntArray> >::Type bitsets_iter = begin(bitsets);
			typename Iterator< String<TIntArray> >::Type bitsets_end = end(bitsets);
			while(bitsets_iter!=bitsets_end)
			{
				TIntVect hash_val_of_variants;
				_buildVariants(hash_val_of_variants, l_mer_begin, i, *bitsets_iter, shape);
				std::copy(hash_val_of_variants.begin(), hash_val_of_variants.end(), std::back_inserter(variants));
				++bitsets_iter;
			}
			clear(bitsets);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._getVariantsOfBitset:
..summary:Builds all variants of bitsets consisting of d zeros & (l-d) ones which
          will be stored in an seqan::String.
..cat:Motif Search
..signature:_getVariantsOfBitset(bitsets,l,d)
..param.bitsets:The collection of bitsets.
...type:Class.String
...signature:String< String<int> >
...remarks:Bitsets contains different length-l arrays consisting of d zero- (position of mismatches) 
           and (l-d) one-values
..param.l:The size of the motif.
..param.d:The number of substitutions.
..remarks: #bitsets = binomialCoefficient(l,d).
..include:seqan/find_motif.h
*/

template<typename TStrings, typename TType>
void 
_getVariantsOfBitset(TStrings & bitsets, 
					 TType const & l, 
					 TType const & d)
{
    SEQAN_CHECKPOINT;

	unsigned int num_of_bitsets = binomialCoefficient(l, d);
	resize(bitsets, num_of_bitsets);

	typename Value<TStrings>::Type bitset;
	resize(bitset, l);
	std::fill(begin(bitset), end(bitset), 1);
	std::fill_n(begin(bitset), d, 0);

	for(typename Position<TStrings>::Type i=0;i<num_of_bitsets; ++i)
	{
		bitsets[i] = bitset;
		std::next_permutation(begin(bitset), end(bitset));
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._buildVariants:
..summary:Builds all d-variants of a given l-mer given the position(s) of mismatch(es) and 
          then computes the hash-values for each variant which will be finally stored 
          in an int-vector variants.
..cat:Motif Search
..signature:_buildVariants(variants,l_mer_begin,d,bitset,shape)
..param.variants:The int-vector of d-variants of a given l-mer 
...remarks:The vector holds the hash-values of each d-variant. 
...type:vector<int>
..param.l_mer_begin:An iterator pointing to the beginning of a given l-mer pattern.
...type:Concept.RandomAccessIteratorConcept
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.d:The number of substitutions.
..param.bitset:The bitset object.
...type:Class.String
...signature:String<int>
...remarks:The length-l int-array consisting of d zeros (for mismatch positions)
           & (l-d) ones)
..param.shape:The @Class.Shape@ object.
...type:Class.Shape
...signature:Shape<TValue, SimpleShape>
...remarks:Is used for computing the hash-value of each variant.
..remarks: #variants = (alp_size-1)^d for a given bitset, d-value and l-mer.
..include:seqan/find_motif.h
*/

template<typename TIntVect, typename TStringIter, typename TType, typename TBitset, typename TValue, typename TSpec>
void
_buildVariants(TIntVect & variants,
			   TStringIter l_mer_begin,
			   TType const & d,
			   TBitset const & bitset,
			   Shape<TValue, TSpec> & shape)
{
    SEQAN_CHECKPOINT;

	typedef String<TValue> TString;
	typedef typename Position<TString>::Type TPos;

	//build l-mer
	TString l_mer;
	resize(l_mer, length(bitset));
	for(TPos i=0; i<length(l_mer); ++i)
	{
		l_mer[i] = *l_mer_begin;
		++l_mer_begin;
	}

	// d-mers
	unsigned int num_of_d_mers = 
		(unsigned int) pow((double)ValueSize<TValue>::VALUE, (int)d);
	for(unsigned int i=0; i<num_of_d_mers; ++i)
	{
		TString d_mer = inverseHash<TValue>(i, ValueSize<TValue>::VALUE, d);
		TString variant = l_mer;

		typename Position<TString>::Type pos_of_d_mer = 0;
		typename Position<TString>::Type pos_of_variant = 0;
		typename Iterator<TString>::Type iter = begin(variant);
		typename Iterator<TString>::Type iter_end = end(variant);
		while(iter!=iter_end)
		{
			if( bitset[pos_of_variant]==0)
			{
				if(l_mer[pos_of_variant]!=d_mer[pos_of_d_mer])
				{
					variant[pos_of_variant] = d_mer[pos_of_d_mer];
					++pos_of_d_mer;
					++pos_of_variant;
					++iter;
				}
				else
				{
					iter = end(variant);
				}
			}
			else
			{
				++iter;
				++pos_of_variant;
			}
		}
		if(pos_of_d_mer==d)
		{
			variants.push_back(hash(shape, begin(variant)));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

