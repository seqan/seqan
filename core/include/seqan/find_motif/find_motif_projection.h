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

#ifndef SEQAN_HEADER_FIND_MOTIF_PROJECTION_H
#define SEQAN_HEADER_FIND_MOTIF_PROJECTION_H

namespace SEQAN_NAMESPACE_MAIN
{

// Manual Forwards
template <typename TType> TType _computeProjectionSize(TType const & alp_size, TType const & l, TType const & d,
                                                       TType const & m_total);
template <typename TType> TType _computeBucketThreshold(TType const & alp_size, TType const & l, TType const & d,
                                                        TType const & m_total, TType const & k);
template <typename TBucketAr, typename TArray, typename TType, typename TStrings, typename TPositions>
void _filteringStep(TBucketAr & buckets, TArray & count_ar, TType & num_of_relevant_buckets, TStrings & dataset,
                    TPositions & positions, TType const & l, TType const & s);
template<typename TType> 
TType _computeNumOfTrials(TType const & t, TType const & l, TType const & d, TType const & k, TType const & s,
                          double const & prob_q);

//////////////////////////////////////////////////////////////////////////////
// Projection
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Projection:
..summary: Represents the PROJECTION algorithm of Buhler and Tompa.
..description.note:There are various known problems with the motif finding in SeqAn. We plan to fix this in an upcoming release.
..general:Class.MotifFinder
..cat:Motif Search
..signature:MotifFinder<TValue, Projection, TRng>
..param.TValue:The type of sequences to be analyzed.
...type:Spec.Dna
...type:Spec.AminoAcid
..remarks:The @Spec.Projection|Projection algorithm@ is a heuristic algorithm that does not guarantee
          that the unknown motif will be found every time. We can increase the chance of success
		  by performing a large number of independent trials to generate multiple guesses.
		  In each trial, @Spec.Projection@ makes a preselection of sets of length-l patterns called l-mers
		  which are likely to be a collection of motif instances (filtering step) and 
		  refines them by some local searching techniques, e.g. @Function.em|EM algorithm@, Gibbs Sampling etc (refinement step).
..param.TRng:The @Class.Rng@ specialization to use for random number generation.
...default:$GetDefaultRng<MotifFinderClass>::Type$
..include:seqan/find_motif.h
*/

struct Projection_;
typedef Tag<Projection_> const Projection;

//////////////////////////////////////////////////////////////////////////////
// MotifFinder - Projection Spec
//
// t:=dataset size (number of sequences)
// l:=motif size
// m:=number of possible l-mers
// d:=number of substitutions
// k:=projection size
// s:=bucket size
// tr:=number of independent trials
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TRng>
class MotifFinder<TValue, Projection, TRng>
{
	enum { ALPHABET_SIZE = ValueSize<TValue>::VALUE };
//____________________________________________________________________________________________
	/*
	enum { ALPHABET_SIZE = ValueSize<TValue>::VALUE };
	typedef String<TValue> TString;
	typedef String<TString> TStrings;
	typedef typename Size<TStrings>::Type TSize1;
	typedef typename Size<TString>::Type TSize2;*/

//____________________________________________________________________________________________

public:
	typedef unsigned int TSize;
	typedef String<TValue> TString;
	typedef String<TString> TStrings;

    Holder<TRng> _rng;
	TSize dataset_size;
	TSize motif_size;
	TSize total_num_of_l_mers;
	TSize num_of_substitutions;
	bool has_exact_substitutions;
	TSize projection_size;
	TSize bucket_threshold;
	TSize num_of_trials;
	TStrings set_of_motifs;
	int score;

//____________________________________________________________________________________________

	MotifFinder():
        _rng(defaultRng<TRng>()),
		dataset_size(0),
		motif_size(0),
		total_num_of_l_mers(0),
		num_of_substitutions(0),
		has_exact_substitutions(false),
		projection_size(0),
		bucket_threshold(0),
		num_of_trials(0)
	{
    SEQAN_CHECKPOINT;
	}
	MotifFinder(TSize const & t_, 
				TSize const & l_, 
				TSize const & m_total_,
				TSize const & d_,
				bool const & is_exact_,
				TSize const & k_,
				TSize const & s_,
				TSize const & tr_):
        _rng(defaultRng<TRng>()),
		dataset_size(t_),
		motif_size(l_),
        total_num_of_l_mers(m_total_),
		num_of_substitutions(d_),
		has_exact_substitutions(is_exact_),
		projection_size(k_),
		bucket_threshold(s_),
		num_of_trials(tr_)
	{
    SEQAN_CHECKPOINT;
	}
	MotifFinder(TSize const & t_, 
				TSize const & l_, 
				TSize const & m_total_,
				TSize const & d_,
				bool const & is_exact_):
        _rng(defaultRng<TRng>()),
		dataset_size(t_),
		motif_size(l_),
		total_num_of_l_mers(m_total_),
		num_of_substitutions(d_),
		has_exact_substitutions(is_exact_),
		projection_size(0),
		bucket_threshold(0),
		num_of_trials(0)
	{
    SEQAN_CHECKPOINT;

		projection_size = 
			_computeProjectionSize<TSize>(ALPHABET_SIZE,
			                                     motif_size,
								                 num_of_substitutions,
								                 total_num_of_l_mers);
		
		bucket_threshold = 
			_computeBucketThreshold<TSize>(ALPHABET_SIZE,
			                                      motif_size,
								                  num_of_substitutions,
								                  total_num_of_l_mers,
												  projection_size);

		double prob_q = static_cast<double>(0.95);
		num_of_trials = _computeNumOfTrials(dataset_size,
								            motif_size,
								            num_of_substitutions,
								            projection_size,
								            bucket_threshold,
								            prob_q);
	}
	MotifFinder(MotifFinder const & other_):
        _rng(other_.rng),
		dataset_size(other_.dataset_size),
		motif_size(other_.motif_size),
		total_num_of_l_mers(other_.total_num_of_l_mers),
		num_of_substitutions(other_.num_of_substitutions),
		has_exact_substitutions(other_.has_exact_substitutions),
		projection_size(other_.projection_size),
		bucket_threshold(other_.bucket_threshold),
		num_of_trials(other_.num_of_trials)
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
			this->dataset_size = other_.dataset_size;
			this->motif_size = other_.motif_size;
			this->total_num_of_l_mers = other_.total_num_of_l_mers;
			this->num_of_substitutions = other_.number_of_substitutions;
			this->has_exact_substitutions = other_.has_exact_substitutions;
			this->projection_size = other_.projection_size;
			this->bucket_threshold = other_.bucket_threshold;
			this->num_of_trials = other_.num_of_trials;
		}

		return *this;
	}

//____________________________________________________________________________________________

}; // class MotifFinder<TValue, Projection>


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeProjectionSize:
..summary:Computes the projection size (k).
..cat:Motif Search
..signature:_computeProjectionSize(alp_size,l,d,m)
..param.alp_size:The size of the sequence alphabet.
...remarks:The alp_size object is four for nucleotide sequences and twenty for amino acid sequences.
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.m:The total number of possible l-mers of a given dataset.
..include:seqan/find_motif.h
*/

template<typename TType> 
TType
_computeProjectionSize(TType const & alp_size,
					   TType const & l, 
					   TType const & d,
					   TType const & m_total)
{
    SEQAN_CHECKPOINT;

	TType result = static_cast<TType>(0);
	double numerator = log(static_cast<double>(m_total));
	double denominator = log(static_cast<double>(alp_size));
	double fraction = (numerator/denominator);

	if( fraction<static_cast<double>(l-d-1) )
	{
		result = static_cast<TType>(floor(fraction)+1);
	}
	else
	{
		result = l-d-1;
	}

	return result;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeBucketThreshold:
..summary:Computes the bucket threshold size (s).
..cat:Motif Search
..signature:_computeBucketThreshold(alp_size,l,d,m,k)
..param.alp_size:The size of the sequence alphabet.
...remarks:The alp_size object is four for nucleotide sequences and twenty for amino acid sequences.
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.m:The total number of possible l-mers of a given dataset.
..param.k:The projection size.
..include:seqan/find_motif.h
*/

template<typename TType> 
TType
_computeBucketThreshold(TType const & alp_size,
					    TType const & l, 
					    TType const & d,
					    TType const & m_total, 
						TType const & k)
{
    SEQAN_CHECKPOINT;

	TType result = 3; // or 4
	double numerator = log(static_cast<double>(m_total));
	double denominator = log(static_cast<double>(alp_size));
	double fraction = (numerator/denominator);

	if( fraction>=static_cast<double>(l-d-1) )
	{
		result =
			static_cast<TType>(floor(static_cast<double>(m_total/
			                          pow(static_cast<double>(alp_size), static_cast<double>(k))*2)));
		if(result<1)
		{
			result = 1;
		}
	}

	return result;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeNumOfTrials:
..summary:Computes the number of independent trials (tr).
..cat:Motif Search
..signature:_computeNumOfTrials(t,l,d,k,s,prob_q)
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.k:The projection size.
..param.s:The bucket threshold size.
..param.prob_q:
...remarks:The prob_q object represents the probability that the planted bucket contains "s" or 
           more planted motif instances in at least one of the "tr" trials. Normally, we use 
		   prob_q=0.95.	
...type:$double$
..remarks:tr>= log(1-q)/log(B), where p is the probability that each motif occurrence hashes 
          to the planted bucket and B is the probability that fewer than s planted occurrences hash
          to the planted buckes in a given trial
..include:seqan/find_motif.h
*/

template<typename TType> 
TType
_computeNumOfTrials(TType const & t,
					TType const & l,
					TType const & d, 
					TType const & k, 
					TType const & s,
					double const & prob_q)
{
    SEQAN_CHECKPOINT;

	double prob_p =
		static_cast<double>(binomialCoefficient( (l-d), k ))
	   /static_cast<double>(binomialCoefficient(l,k));
	
	double prob_B = static_cast<double>(0);
	for(unsigned int i=0; i<s; ++i)
	{
		prob_B+=
			static_cast<double>(binomialCoefficient(t,i))
		   *pow(prob_p, static_cast<double>(i))
		   *pow(static_cast<double>(1)-prob_p, static_cast<double>(t-i));
	}

	double numerator = log(static_cast<double>(1)-prob_q);
	double denominator = log(static_cast<double>(prob_B));
	TType result = 
		static_cast<TType>(ceil(static_cast<double>(numerator/denominator)-static_cast<double>(0.5)));

	if(result<1)
	{
		result = 1;
	}
	
	return result;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function.findMotif (for PROJECTION)
..summary:Represents the PROJECTION algorithm.
..cat:Motif Search
..signature:findMotif(finder,dataset,seq_model)
..param.finder:The @Class.MotifFinder@ object.
...type:Class.MotifFinder
..param.dataset:The dataset object representing the input sequences.
...type:Class.StringSet
..param.seq_model:The seq_model object.
...type:Tag.Oops
...type:Tag.Zoops
...type:Tag.Tcm
...remarks:The sequence models rely on different assumptions about the distribution of motif occurrences
           across the sample sequences. 
..remarks:The PROJECTION algorithm which consists of two steps, the filtering and the refinement step,
          is able to run in Oops, Zoops and Tcm mode.
..remarks:The algorithm uses the EM procedure during the refinement phase which was introduced by Bailey 
          and Elkan.
..include:seqan/find_motif.h
*/

template<typename TSeqType, typename TStrings, typename TModel, typename TRng>
void
findMotif(MotifFinder<TSeqType, Projection, TRng> & finder, 
		  TStrings & dataset,
		  TModel seq_model)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Position<TString>::Type TPos;
	typedef String<int> TArray; 
	typedef std::vector<int> TBucket;
	typedef typename Size<TArray>::Type TArSize;
	typedef String<TBucket> TBuckets;

	// dataset information
	typedef typename Size<TStrings>::Type TStringsSize;
	TStringsSize t = length(dataset);

	// count_ar:=array of votes for each h(k-mer)
	TArSize ar_size =
		static_cast<TArSize>(pow(static_cast<double>(ValueSize<TSeqType>::VALUE), static_cast<int>(finder.projection_size)));
	finder.score = -1;
	//std::set< String<int> > occurred_positions;
	for(unsigned int trial=0; trial<finder.num_of_trials; ++trial)
	{
		//
		//std::cout << " . ";
		//

		TArray count_ar;
		resize(count_ar, ar_size);

		// ----------------------------------------------------------------------------
		// STEP 1:
		// filtering phase (:=key random projection phase)
		// ----------------------------------------------------------------------------
		// choose randomly k different positions
		std::set<int> positions;
		choosePositions(positions,finder.motif_size,finder.projection_size, value(finder._rng));

		// array of collection of l-mers
		TBuckets bucket_ar;
		unsigned int num_of_relevant_buckets = 0;

		_filteringStep(bucket_ar,
					   count_ar,
					   num_of_relevant_buckets,
					   dataset,
					   positions,
					   finder.motif_size,
					   finder.bucket_threshold);

		// ----------------------------------------------------------------------------
		// STEP 2:
		// checking phase (:= local search-based refinement procedure)
		// ----------------------------------------------------------------------------
		TPos i = 0;
		TPos j = 0;
		while( (j<num_of_relevant_buckets) & (i<ar_size) )
		{
			unsigned int bucket_size = (bucket_ar[i]).size();
			if(bucket_size>=finder.bucket_threshold)
			{
				++j;

				TStrings l_mers;
				resize(l_mers, bucket_size);
				TBucket::iterator bucket_iter, bucket_end;
				bucket_iter = (bucket_ar[i]).begin();
				bucket_end = (bucket_ar[i]).end();
				int bucket_element = 0;
				while(bucket_iter!=bucket_end)
				{
					TString l_mer = 
						inverseHash<TSeqType>(*bucket_iter, valueSize<TSeqType>(), finder.motif_size);
					l_mers[bucket_element] = l_mer;
					++bucket_element;
					++bucket_iter;
				}

				TString consensus_pat;
				TStrings motif_instances;
				int score = 
					_refinementStep(consensus_pat,
									l_mers,dataset,
									finder.motif_size,
									finder.num_of_substitutions,
									finder.has_exact_substitutions,
									seq_model);

				if(score>finder.score)
				{
					finder.score = score;
					if (!length(finder.set_of_motifs)) appendValue(finder.set_of_motifs, consensus_pat);
					else finder.set_of_motifs[0] = consensus_pat;

					if((TStringsSize) finder.score==t)
					{
						trial = finder.num_of_trials;
						j = num_of_relevant_buckets;
					}
				}
			}
			++i;
		}// end while( (j<num_of_relevant_buckets) & (i<ar_size) )
	}// end for(unsigned int trial=0; trial<finder.num_of_trials; ++trial)
	
	//
	//std::cout << "\n";
	//
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function.projectLMer:
..summary:Based on set "positions" the function uses the letters of a given l-mer at
          chosen k positions to compute an appropriate hash value of the new k-mer.
..cat:Motif Search
..signature:projectLMer(positions,l,k)
..param.positions:The set of k chosen positions.
...remarks:$positions$ is of type $set<int>$
..param.k:An iterator pointing to the first positions of a given sequence.
..include:seqan/find_motif.h
*/

template<typename TValue, typename TIter>
inline std::set<int>::value_type
projectLMer(std::set<int> & positions, TIter it)
{
    SEQAN_CHECKPOINT;

	typedef std::set<int>::value_type THValue;
	THValue prev_pos; //, cur_pos;

	std::set<int>::iterator positions_iter, positions_end;
	positions_iter = positions.begin();
	positions_end = positions.end();
	prev_pos = *positions_iter;
	it += prev_pos;
	THValue hValue = ordValue(*it);
	++positions_iter;
	while(positions_iter!=positions_end)
	{
		THValue cur_pos = *positions_iter;
		goFurther(it, (cur_pos-prev_pos));
		hValue = hValue * ValueSize<TValue>::VALUE + ordValue(*it);
		prev_pos = cur_pos;
		++positions_iter;
	}

	return hValue;
}
//////////////////////////////////////////////////////////////////////////////

/*
.Function._filteringStep:
..summary:Given a position set with k different positions we compute a projection value
          for each l-mer in the input sequences and store the specific l-mer in the appropriate
		  bucket which is labeled with the specific projection value.
..cat:Motif Search
..signature:_filteringStep(buckets,count_ar,num_of_relevant_buckets,dataset,shape,l,s)
..include:seqan/find_motif.h
*/

template<typename TBucketAr, typename TArray, typename TType, typename TStrings, typename TPositions>
void 
_filteringStep(TBucketAr & buckets, 
			   TArray & count_ar,
			   TType & num_of_relevant_buckets,
			   TStrings & dataset,
			   TPositions & positions,
			   TType const & l,
			   TType const & s)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	//typedef typename Value<TBucketAr>::Type TBucket;
	typename Iterator<TStrings>::Type ds_iter = begin(dataset);
	typename Size<TArray>::Type ar_size = length(count_ar);
	Shape<TValue> shape(l); //to compute hash value of l-mer x
 
	// initialize pointer by setting it to null 
	// (=std::fill(begin(count_ar),end(count_ar),0))
	typename Iterator<TArray>::Type count_ar_iter = begin(count_ar) ;
	typename Iterator<TArray>::Type count_ar_end = end(count_ar);
	while(count_ar_iter!=count_ar_end)
	{
		*count_ar_iter = 0;
		++count_ar_iter;
	}

	// go over input sequences & increment corresponding counter in count_ar
	// fill l_mer_index with entries
	resize(buckets, ar_size);
	int y = 0; //hash-value of created k-mer
	int x = 0; //hash-value of l-mer
	for(; !atEnd(ds_iter, dataset); goNext(ds_iter))
	{
		typename Size<TString>::Type seq_length = length(*ds_iter);
		typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
		typename Iterator<TString>::Type seq_end = begin(*ds_iter)+(seq_length-l+1);
		while( seq_iter!=seq_end )
		{
			x = hash(shape, seq_iter);
			y = projectLMer<TValue>(positions, seq_iter);
	    	++count_ar[y];
			(buckets[y]).push_back(x);
			++seq_iter;
		}
	}

	num_of_relevant_buckets = 
		std::count_if(begin(count_ar),
		              end(count_ar),
					  bind2nd(std::greater_equal<int>(),static_cast<int>(s)));
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._refinementStep:
..summary:Refines the collection of l-mers in each relevant bucket which contains at least s l-mers.
..cat:Motif Search
..signature:_refinementStep(consensus_seq,positions,l_mers,dataset,t,l,d,is_exact,model_type)
..include:seqan/find_motif.h
*/

//////////////////////////////////////////////////////////////////////////////
//	Oops model
//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TType>
int
_refinementStep(TString & consensus_seq,
			    String<TString> const & l_mers,
			    String<TString> & dataset,
				TType const & l,
				TType const & d,
				bool const & is_exact,
				Oops const & /*oops*/)
{
    SEQAN_CHECKPOINT;

	typedef String<TString> TStrings;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TString>::Type TPos;
	typedef FrequencyDistribution<TValue> TFrequencyDistribution;
	TType t = length(dataset);
	int score = 0;

	// compute background frequency
	TFrequencyDistribution background;
	backgroundFrequency(background, begin(dataset), end(dataset));

	// step1: initial guess (profile) from bucket
	double epsilon = 0.1;
	Pseudocount<TValue, CMode> pseudocount_mode_c(epsilon);
	String<TFrequencyDistribution> profile;
	convertSetOfPatternsToProfile(profile, l_mers, pseudocount_mode_c);
	completeProfile(profile, background);

	// step2: refinement of initial profile with em: 5 trials
	// double likelihood_score = 0;  // TODO(holtgrew): Why is this read nowhere?
	int iterations = 3; //5

	while(iterations>0)
	{
		// likelihood_score = em(profile, begin(dataset), t, l, oops);
		--iterations;
	}

	determineConsensusSeq(consensus_seq, profile, l);
	typename Iterator<TStrings>::Type ds_iter, ds_end;
	ds_iter = begin(dataset);
	ds_end = end(dataset);
	typename Position<TStrings>::Type seq_nr;
	do
	{
		seq_nr = t-(ds_end-ds_iter);
		TPos m = (TPos)(length(*ds_iter)-l+1);
		int * hd_ar = new int[m];
		typename Iterator<TString>::Type seq_iter, seq_end, consensus_begin;
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			consensus_begin = begin(consensus_seq);
			hd_ar[m-(seq_end-seq_iter)] = hammingDistance<int>(seq_iter, seq_iter+l, consensus_begin);
			++seq_iter;
		}

		if( (!is_exact) & 
			(count_if(hd_ar,hd_ar+m,bind2nd(std::less_equal<TType>(), d))==1) )
		{
			++score;
		}
		else if( is_exact & 
			     (count_if(hd_ar,hd_ar+m,bind2nd(std::equal_to<TType>(), d))==1) )
		{
			++score;
		}

		delete[] hd_ar;
		++ds_iter;
	}
	while( (ds_iter!=ds_end) & (score==(int) seq_nr+1) );

	return score;
}

//////////////////////////////////////////////////////////////////////////////
//	Omops model
//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TType>
int
_refinementStep(TString & consensus_seq,
			    String<TString> const & l_mers,
			    String<TString> & dataset,
				TType const & l,
				TType const & d,
				bool const & is_exact,
				Omops const & /*omops*/)
{
	typedef String<TString> TStrings;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TString>::Type TPos;
	typedef FrequencyDistribution<TValue> TFrequencyDistribution;
	TType t = length(dataset);
	int score = 0;

	// compute background frequency
	TFrequencyDistribution background;
	backgroundFrequency(background, begin(dataset), end(dataset));

	// step1: initial guess (profile) from bucket
	double epsilon = 0.1;
	Pseudocount<TValue, CMode> pseudocount_mode_c(epsilon);
	String<TFrequencyDistribution> profile;
	convertSetOfPatternsToProfile(profile, l_mers, pseudocount_mode_c);
	completeProfile(profile, background);

	// step2: refinement of initial profile with em: 5 trials
	// double likelihood_score = 0;  // TODO(holtgrew): Why is this read nowhere?
	int iterations = 3; //5

	while(iterations>0)
	{
		// likelihood_score = em(profile, begin(dataset), t, l, Oops());
		--iterations;
	}

	determineConsensusSeq(consensus_seq, profile, l);
	typename Iterator<TStrings>::Type ds_iter, ds_end;
	ds_iter = begin(dataset);
	ds_end = end(dataset);
	typename Position<TStrings>::Type seq_nr;
	do
	{
		seq_nr = t-(ds_end-ds_iter);
		TPos m = (TPos)(length(*ds_iter)-l+1);
		typename Iterator<TString>::Type seq_iter, seq_end, consensus_begin;
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			consensus_begin = begin(consensus_seq);
			TType hd = hammingDistance<TType>(seq_iter, seq_iter+l, consensus_begin);

			if( (!is_exact) & (hd<=d) )
			{
				++score;
				seq_iter = seq_end-1;
			}
			else if( is_exact & (hd==d) )
			{
				++score;
				seq_iter = seq_end-1;
			}
			++seq_iter;
		}
		++ds_iter;
	}
	while( (ds_iter!=ds_end) & (score== (int) seq_nr+1) );

	return score;
}

//////////////////////////////////////////////////////////////////////////////
//	Zoops model
//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TType>
int
_refinementStep(TString & consensus_seq,
			    String<TString> const & l_mers,
			    String<TString> & dataset,
				TType const & l,
				TType const & d,
				bool const & is_exact,
				Zoops const & zoops)
{
    SEQAN_CHECKPOINT;

	typedef String<TString> TStrings;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TString>::Type TPos;
	typedef FrequencyDistribution<TValue> TFrequencyDistribution;
	TType t = length(dataset);
	int score = 0;
	int lower_limit = (int) floor(t*(zoops.threshold)+0.5);

	// compute background frequency
	TFrequencyDistribution background;
	backgroundFrequency(background, begin(dataset), end(dataset));

	// step1: initial guess (profile) from bucket
	double epsilon = 0.1;
	Pseudocount<TValue, CMode> pseudocount_mode_c(epsilon);
	String<TFrequencyDistribution> profile;
	convertSetOfPatternsToProfile(profile, l_mers, pseudocount_mode_c);
	completeProfile(profile, background);

	// step2: refinement of initial profile with em: 5 trials
	// double likelihood_score = 0;  // TODO(holtgrew): Why is this read nowhere?
	// double gamma = static_cast<double>(1)/sqrt(static_cast<double>(t));
	int iterations = 3; //5
	while(iterations>0)
	{
		// likelihood_score = em(profile, begin(dataset), t, l, gamma, zoops);
		--iterations;
	}

	determineConsensusSeq(consensus_seq, profile, l);
	typename Iterator<TStrings>::Type ds_iter, ds_end;
	ds_iter = begin(dataset);
	ds_end = end(dataset);
	// typename Position<TStrings>::Type seq_nr;  // TODO(holtgrew): Why is this read nowhere?
	do
	{
		// seq_nr = t-(ds_end-ds_iter);
		TPos m = (TPos)(length(*ds_iter)-l+1);
		int * hd_ar = new int[m];
		typename Iterator<TString>::Type seq_iter, seq_end, consensus_begin;
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			consensus_begin = begin(consensus_seq);
			hd_ar[m-(seq_end-seq_iter)] = hammingDistance<int>(seq_iter, seq_iter+l, consensus_begin);
			++seq_iter;
		}

		if(!is_exact)
		{
			TType num = count_if(hd_ar,hd_ar+m,bind2nd(std::less_equal<TType>(), d));
			if(num>1)
			{
				ds_iter = ds_end-1;
			}
			else if(num==1)
			{
				++score;
			}
		}
		else
		{
			TType num = count_if(hd_ar,hd_ar+m,bind2nd(std::equal_to<TType>(), d));
			if(num>1)
			{
				ds_iter = ds_end-1;
			}
			else if(num==1)
			{
				++score;
			}
		}
		delete[] hd_ar;
		++ds_iter;
	}
	while(ds_iter!=ds_end);

	if(score<lower_limit)
	{
		score = 0;
	}

	return score;
}

//////////////////////////////////////////////////////////////////////////////
//	Tcm model
//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TType>
int 
_refinementStep(TString & consensus_seq,
			    String<TString> const & l_mers,
			    String<TString> & dataset,
				TType const & l,
				TType const & d,
				bool const & is_exact,
				Tcm const & tcm)
{
    SEQAN_CHECKPOINT;

	typedef String<TString> TStrings;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TString>::Type TPos;
	typedef FrequencyDistribution<TValue> TFrequencyDistribution;
	TType t = length(dataset);
	int score = 0;
	int lower_limit = (int) floor(t*(tcm.threshold)+0.5);

	// compute background frequency
	TFrequencyDistribution background;
	backgroundFrequency(background, begin(dataset), end(dataset));

	// step1: initial guess (profile) from bucket
	double epsilon = 0.1;
	Pseudocount<TValue, CMode> pseudocount_mode_c(epsilon);
	String<TFrequencyDistribution> profile;
	convertSetOfPatternsToProfile(profile, l_mers, pseudocount_mode_c);
	completeProfile(profile, background);

	// step2: refinement of initial profile with em: 5 trials
	double gamma = static_cast<double>(1)/sqrt(static_cast<double>(t));
	double lambda = 0;
	typename Iterator<TStrings>::Type ds_iter, ds_end;
	ds_iter = begin(dataset);
	ds_end = end(dataset);
	while(ds_iter!=ds_end)
	{
		TPos m = (TPos)(length(*ds_iter)-l+1);
		lambda += (gamma/((double)m));
		++ds_iter;
	}
	lambda = lambda/((double)t);

	double likelihood_score = 0;
	int iterations = 3; //3
	while(iterations>0)
	{
		likelihood_score = em(profile, begin(dataset), t, l, lambda, tcm);
		--iterations;
	}
    (void)likelihood_score;  // TODO(holtgrew): Why is this never read?

	determineConsensusSeq(consensus_seq, profile, l);
	ds_iter = begin(dataset);
	// typename Position<TStrings>::Type seq_nr;
	do
	{
		// seq_nr = t-(ds_end-ds_iter);
		TPos m = (TPos)(length(*ds_iter)-l+1);
		typename Iterator<TString>::Type seq_iter, seq_end, consensus_begin;
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			consensus_begin = begin(consensus_seq);
			TType hd = hammingDistance<TType>(seq_iter, seq_iter+l, consensus_begin);

			if( (!is_exact) & (hd<=d) )
			{
				++score;
				seq_iter = seq_end-1;
			}
			else if( is_exact & (hd==d) )
			{
				++score;
				seq_iter = seq_end-1;
			}
			++seq_iter;
		}
		++ds_iter;
	}
	while(ds_iter!=ds_end);

	if(score<lower_limit)
	{
		score = 0;
	}

	return score;
}

//////////////////////////////////////////////////////////////////////////////
//Subfunctions
//////////////////////////////////////////////////////////////////////////////

/*
.Function.choosePositions:
..summary:Chooses randomly k different positions from {0,1,...,(l-1)}
..cat:Motif Search
..signature:choosePositions(positions,l,k,rng)
..param.positions:The set of k chosen positions.
...type:$set<int>$
..param.l:The size of the motif.
..param.k:The projection size.
..param.rng:The @Class.Rng@ object to get random numbers from.
..include:seqan/find_motif.h
*/

template<typename TAssociativeContainer, typename TType, typename TRng>
void
choosePositions(TAssociativeContainer & positions, TType const & l, TType const & k, TRng & rng)
{
    SEQAN_CHECKPOINT;
    
    SEQAN_ASSERT_GT(l, static_cast<TType>(0));
    Pdf<Uniform<int> > positionPdf(0, l - 1);  // PDF for picking positions.
    
	while (positions.size() < k)
	{
		int position = pickRandomNumber(rng, positionPdf);
		positions.insert(position);
	}
}


//////////////////////////////////////////////////////////////////////////////

/*
.Function._getLMersWithTheLargestLikelihoodRatio:
..summary:Forms a guess for the planted motif by selecting from each input sequence
		  the l-mer x with the largest likelihood ratio.
..signature:_getLMersWithTheLargestLikelihoodRatio(l_mers,positions,dataset_start,dataset_end,profile,l)
..param.l_mers:The collection of t l-mers.
..param.dataset_start:An iterator pointing to the first input sequence of a given dataset.
..param.dataset_end:An iterator pointing to the last input sequence of a given dataset.
..param.t:The number of input sequences.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.l:The size of the motif.
..include:seqan/find_motif.h
*/

template<typename TStrings, typename TIter, typename TType, typename TProfile>
void
_getLMersWithTheLargestLikelihoodRatio(TStrings & l_mers,
									   TIter dataset_start,
									   TIter dataset_end,
									   TProfile const & profile,
								       TType const & l)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TStrings>::Type TString;
	//typedef typename Position<TStrings>::Type TPos1;
	typedef typename Position<TString>::Type TPos2;
	typename Size<TStrings>::Type t = (dataset_end-dataset_start);
	resize(l_mers, t);
	while(dataset_start!=dataset_end)
	{
		typename Position<TStrings>::Type seq_nr = t-(dataset_end-dataset_start);
		TPos2 m = (TPos2)(length(*dataset_start)-l+1);
		double * likelihood_ratio_mat = new double[m];
		typename Iterator<TString>::Type seq_iter, seq_end;
		seq_iter = begin(*dataset_start);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			likelihood_ratio_mat[m-(seq_end-seq_iter)] = 
				_computeLikelihoodRatioOfLMer(seq_iter, seq_iter+l, profile);
			++seq_iter;
		}

		TPos2 max_pos = 
			(std::max_element(likelihood_ratio_mat, likelihood_ratio_mat+m)-likelihood_ratio_mat);
		l_mers[seq_nr] = infix(*dataset_start, max_pos, max_pos+l);

		delete[] likelihood_ratio_mat;
		++dataset_start;
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeLikelihoodRatioOfLMer:
..summary:Computes the likelihood ratio of a given l-mer.
..signature:_computeLikelihoodRatioOfLMer(l_mer_begin,l_mer_end,profile)
..param.l_mer_begin:An iterator pointing to the beginning of a given l-mer pattern.
...type:Concept.RandomAccessIteratorConcept
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.l_mer_end:An iterator pointing to the end of a given l-mer pattern.
...type:Concept.RandomAccessIteratorConcept
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..remarks:Computes the sum of log probabilites instead of the product of probabilites
..include:seqan/find_motif.h
*/

template<typename TStrIter, typename TProfile>
double
_computeLikelihoodRatioOfLMer(TStrIter l_mer_begin, 
							  TStrIter l_mer_end,
							  TProfile const & profile)
{
    SEQAN_CHECKPOINT;

	double result = 0;
	unsigned int l = (unsigned int)(l_mer_end-l_mer_begin);
	typedef typename Position<TProfile>::Type TPos;
	TProfile log_profile = profile;
	for(TPos i=0; i<length(log_profile); ++i)
	{
		logarithmize(log_profile[i]);
	}

	double motif_component = 0;
	double backgr_component = 0;
	while(l_mer_begin!=l_mer_end)
	{
		motif_component += log_profile[l-(l_mer_end-l_mer_begin)+1][(int)*l_mer_begin];
		backgr_component += log_profile[0][(int)*l_mer_begin];
		++l_mer_begin;
	}
	result = motif_component-backgr_component;

	return exp(result);
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeLikelihoodRatioOfLMers:
..summary:Computes the likelihood ratio of a given set of l-mers.
..signature:_computeLikelihoodRatioOfLMers(l_mers,profile)
..param.l_mers:The collection of l-mers.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..include:seqan/find_motif.h
*/

template<typename TStrings, typename TProfile>
double
_computeLikelihoodRatioOfLMers(TStrings const & l_mers, 
							   TProfile const & profile)
{
    SEQAN_CHECKPOINT;

	typedef typename Position<TStrings>::Type TPos;
	double score = 1;
	typename Size<TStrings>::Type num_of_l_mers = length(l_mers);
	for(TPos i=0; i<num_of_l_mers; ++i)
	{
		score *= _computeLikelihoodRatioOfLMer(begin(l_mers[i]), end(l_mers[i]), profile);
	}

	return score;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.determineConsensusSeq:
..summary:Determines the consensus pattern of a given profile.
..cat:Motif Search
..signature:determineConsensusSeq(consensus_seq,profile,l)
..param.consensus_seq:The consensus pattern.
...type:Class.String
...type:Shortcut.DnaString
...type:Shortcut.Peptide
..param.profile:A StringSet of @Class.FrequencyDistribution|frequency distributions@.
...type:Class.StringSet
..param.l:The size of the motif.
..include:seqan/find_motif.h
*/

template<typename TString, typename TProfile>
void
determineConsensusSeq(TString & consensus_seq,
					  TProfile & profile,
					  typename Size<TString>::Type const & l)
{
    SEQAN_CHECKPOINT;

	typedef typename Value<TString>::Type TValue;
	typename Position<TString>::Type i;

	resize(consensus_seq, l);
	if(length(profile)==l) //profile only consists of the motif component
	{
		for(i=0; i<l; ++i)
		{
			consensus_seq[i] = 
				static_cast<TValue>(posOfMax(profile[i]));
		}
	}
	else
	{
		for(i=1; i<=l; ++i)
		{
			consensus_seq[i-1] = 
				static_cast<TValue>(posOfMax(profile[i]));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function.displayResult:
..summary:Displays the consensus pattern of the found motif candidate.
..cat:Motif Search
..signature:displayResult(motif_finder)
..param.motif_finder:The @Class.MotifFinder@ object.
...type:Class.MotifFinder
..param.dataset:The dataset object representing the input sequences.
...type:Class.String
...signature:String<TString>
...param.TString:A @Class.String@ type
....type:Class.String
..include:seqan/find_motif.h
*/

template<typename TValue, typename TRng>
void
displayResult(MotifFinder<TValue, Projection, TRng> & projection)
{
    SEQAN_CHECKPOINT;

	if(length(projection.set_of_motifs)!=0)
	{
		std::cout << projection.set_of_motifs[0] << "\n";
	}
	else
	{
		std::cout << "NO MOTIF HAS BEEN FOUND!!!\n";
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Access Functions
//////////////////////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TRng>
inline int
getScore(MotifFinder<TValue, Projection, TRng> & me)
{
    SEQAN_CHECKPOINT;
	return me.score;
}


//////////////////////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
