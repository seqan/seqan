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

#ifndef SEQAN_HEADER_EM_ALGORITHM_H
#define SEQAN_HEADER_EM_ALGORITHM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.em:
..summary:Represents the EM algorithm as used by MEME.
..cat:Motif Search
..signature:em(profile,dataset_start,t,l,oops_model)
..signature:em(profile,dataset_start,t,l,gamma,zoops_model)
..signature:em(profile,dataset_start,t,l,lambda,tcm_model)
..param.profile:A StringSet of @Class.FrequencyDistribution|frequency distributions@.
...type:Class.StringSet
..param.dataset_start:An iterator pointing to the first input sequence of a given dataset.
...type:Concept.RandomAccessIteratorConcept
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.oops_model:The oops_model object.
...type:Tag.Oops
..param.zoops_model:The zoops_model object.
...type:Tag.Zoops
..param.tcm_model:The tcm_model object.
...type:Tag.Tcm
..param.gamma:The probability of sequence having a motif occurrence.
..param.lambda:The probability of starting a motif occurrence 
...remarks:$lambda$ is calculated by dividing $gamma$ by the length of the corresponding sequence.
..remarks:This version of EM is used in the MEME program of Bailey and Elkan. It is a Bayesian
          variant of the basic EM which allows multiple occurrences of a motif in any sequence and can 
		  therefore be performed on sequences of one of the model types @Tag.Oops@, @Tag.Zoops@ and 
		  @Tag.Tcm@. We use the EM algorithm of MEME for the refinement step of PROJECTION.
..include:seqan/find_motif.h
*/

//////////////////////////////////////////////////////////////////////////////
//	Oops model
//////////////////////////////////////////////////////////////////////////////

template<typename TProfile, typename TIter, typename TType>
double
em(TProfile & profile, 
   TIter dataset_start, 
   TType const & t,
   TType const & l,
   Oops const & oops)
{
	// matrix w - allocate space (memory)
	TType row_size = t;
	double ** matrix_w = new double*[row_size];
	for(TType pos=0; pos<row_size; ++pos)
	{
		TType col_size = length(*(dataset_start+pos))-l+1;
		matrix_w[pos] = new double[col_size];
	}

	// E-step: compute matrix z and joint log likelihood
	double log_likelihood = 0;
	_computeEStep(matrix_w,log_likelihood,profile,dataset_start,t,l,oops);

	// M-step: refine profile
	_computeMStep(profile,dataset_start,matrix_w,t,l,oops);

	return log_likelihood;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeEStep:
..summary:Represents the E-step of the EM algorithm for Oops models.
..cat:Motif Search
..signature:_computeEStep(matrix_w,joint_log_likelihood,profile,dataset_start,t,l,oops_model)
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.joint_log_likelihood:The joint log-likelihood.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.oops_model:The oops_model object.
...type:Tag.Oops
..remarks:The joint log likelihood is not computed in the M-step as MEME does, but rather in the E-step
          of the algorithm.
..include:seqan/find_motif.h
*/

template<typename TMatrix, typename TProfile, typename TIter, typename TType>
void
_computeEStep(TMatrix & matrix_w, 
	double & joint_log_likelihood,
	TProfile & profile, 
	TIter dataset_start, //start iterator of dataset 
	TType const & t,
	TType const & l,
	Oops const & /*oops*/)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Value<TFrequencyDist>::Type TValue;
	typedef String<TValue> TString;
	typedef typename Position<TProfile>::Type TPos;
	for(TPos pos=0; pos<length(profile); ++pos)
	{
		logarithmize(profile[pos]); 
	}

	//compute matrix w
	TIter dataset_end = dataset_start+t;
	while(dataset_start!=dataset_end)
	{
		TType seq_nr = t-(dataset_end-dataset_start);
		TType m = length(*dataset_start)-l+1;
		TFrequencyDist absLetterFrequencies;
		absFreqOfLettersInSeq(absLetterFrequencies,begin(*dataset_start),end(*dataset_start));
		typename Iterator<TString>::Type seq_iter, seq_end;
		seq_iter = begin(*dataset_start);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			TFrequencyDist fd = absLetterFrequencies;
			double sum_of_log_probs = 0;
			for(TType h=0; h<l; ++h)
			{
				sum_of_log_probs += profile[h+1][(int)*(seq_iter+h)];
				--fd[(int)*(seq_iter+h)];
			}
			sum_of_log_probs +=
				(double)std::inner_product(begin(fd),end(fd),begin(profile[0]),(double)0);
			matrix_w[seq_nr][m-(seq_end-seq_iter)] = sum_of_log_probs;
			++seq_iter;
		}
		++dataset_start;
	}

	//
	dataset_start -= t;
	while(dataset_start!=dataset_end)
	{
		TType seq_nr = t-(dataset_end-dataset_start);
		TType m = length(*dataset_start)-l+1;
		joint_log_likelihood += (log((double)1)-log((double)m));

		double log_of_sums = matrix_w[seq_nr][0];
		TType j;
		for(j=1; j<m; ++j)
		{
			if( (matrix_w[seq_nr][j]-log_of_sums)>DBL_MIN_EXP )
			{
				log_of_sums += log(1+exp(matrix_w[seq_nr][j]-log_of_sums));
			}
		}

		for(j=0; j<m; ++j)
		{
			double prob = matrix_w[seq_nr][j];
			matrix_w[seq_nr][j] -= log_of_sums;
			matrix_w[seq_nr][j] = exp(matrix_w[seq_nr][j]);
			joint_log_likelihood += matrix_w[seq_nr][j]*prob;
		}

		++dataset_start;
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeStep_M:
..summary:Represents the M-step of the EM algorithm for Oops models.
          Refines the background and motif component.
..cat:Motif Search
..signature:_computeStep_M(profile,dataset_start,matrix_w,t,l,oops_model)
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.oops_model:The oops_model object.
...type:Tag.Oops
..remarks:The joint log likelihood is computed in the E-step of the algorithm.
..include:seqan/find_motif.h
*/

template<typename TProfile, typename TIter, typename TMatrix, typename TType>
void
_computeMStep(TProfile & profile, 
	TIter dataset_start,
	TMatrix const & matrix_w,
	TType const & t,
	TType const & l,
	Oops const & /*oops*/)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Value<TFrequencyDist>::Type TValue;
	typedef typename Position<TProfile>::Type TPos;
	TPos k = 0;
	typedef String<TValue> TString;
	
	// total_counts_of_letters is used for the computation of the background frequency
	TFrequencyDist total_counts_of_letters; //<-c
	absFreqOfLettersInSetOfSeqs(total_counts_of_letters,dataset_start, dataset_start+t);

	//
	for(k=0; k<length(profile); ++k)
	{
		profile[k] = TFrequencyDist();
	}

	// refine motif component
	TIter dataset_end = dataset_start+t;
	while(dataset_start!=dataset_end)
	{
		TType seq_nr = t-(dataset_end-dataset_start);
		TType m = length(*dataset_start)-l+1;
		typename Iterator<TString>::Type seq_iter, seq_end, ptr;
		seq_iter = begin(*dataset_start);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			ptr = seq_iter;
			for(TType i=0; i<l; ++i)
			{
				profile[i+1][(int)*(ptr+i)] += matrix_w[seq_nr][m-(seq_end-seq_iter)];
			}
			++seq_iter;
		}
		++dataset_start;
	}

	for(k=1; k<length(profile); ++k)
	{
		total_counts_of_letters -= profile[k];
	}

	// refine background component
	profile[0] = total_counts_of_letters;
	
	// addPseudocount (if necessary) & normalize
	double epsilon = 0.1;
	normalize(profile, Pseudocount<TValue, CMode>(epsilon));

	// matrix w - deallocate space (memory)
	for(TType pos=0; pos<t; ++pos)
	{
		delete[] matrix_w[pos];
	}
	delete[] matrix_w;
}

//////////////////////////////////////////////////////////////////////////////
//	Zoops model
//////////////////////////////////////////////////////////////////////////////

// gamma=(1/t)*sum(i,1,t,Qi), Qi=sum(j,1,m,zij) (i=1,...,t)

template<typename TProfile, typename TIter, typename TType>
double
em(TProfile & profile, 
   TIter dataset_start, 
   TType const & t,
   TType const & l,
   double & gamma,
   Zoops const & zoops)
{
	// matrix w - allocate space (memory)
	TType row_size = t;
	double ** matrix_w = new double*[row_size];
	for(TType pos=0; pos<row_size; ++pos)
	{
		TType col_size = length(*(dataset_start+pos))-l+1;
		matrix_w[pos] = new double[col_size];
	}

	// compute matrix w and joint log likelihood, refine gamma
	double log_likelihood = 0;
	_computeEStep(matrix_w,log_likelihood,profile,dataset_start,gamma,t,l,zoops);

	// refine profile
	_computeMStep(profile,dataset_start,matrix_w,t,l,zoops);

	return log_likelihood;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeEStep:
..summary:Represents the E-step of the EM algorithm for Zoops models.
..cat:Motif Search
..signature:_computeEStep(matrix_w,joint_log_likelihood,profile,dataset_start,gamma,t,l,zoops_model)
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.joint_log_likelihood:The joint log-likelihood.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.gamma:The probability of sequence having a motif occurrence.
...type:$double$
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.zoops_model:The oops_model object.
...type:Tag.Zoops
..remarks:The joint log likelihood is not computed in the M-step as MEME does, but rather in the E-step
          of the algorithm.
..remarks:The parameter gamma is reestimated during the E-step to save calculation time.
..include:seqan/find_motif.h
*/

template<typename TMatrix, typename TProfile, typename TIter, typename TType>
void
_computeEStep(TMatrix & matrix_w,
			  double & joint_log_likelihood,
			  TProfile & profile,
			  TIter dataset_start,
			  double & gamma,
			  TType const & t,
			  TType const & l,
			  Zoops const & /*zoops*/)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Value<TFrequencyDist>::Type TValue;
	typedef String<TValue> TString;
	typedef typename Position<TProfile>::Type TPos;
	for(TPos pos=0; pos<length(profile); ++pos)
	{
		logarithmize(profile[pos]); 
	}

	double * log_probs_of_motifless_sequences = new double[t];

	//compute matrix w
	TIter dataset_end = dataset_start+t;
	while(dataset_start!=dataset_end)
	{
		TType seq_nr = t-(dataset_end-dataset_start);
		TType m = length(*dataset_start)-l+1;
		TFrequencyDist absLetterFrequencies;
		absFreqOfLettersInSeq(absLetterFrequencies,begin(*dataset_start),end(*dataset_start));
		typename Iterator<TString>::Type seq_iter, seq_end;
		seq_iter = begin(*dataset_start);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			TFrequencyDist fd = absLetterFrequencies;
			double sum_of_log_probs = 0;
			for(TType h=0; h<l; ++h)
			{
				sum_of_log_probs += profile[h+1][(int)*(seq_iter+h)];
				--fd[(int)*(seq_iter+h)];
			}
			sum_of_log_probs +=
				(double)std::inner_product(begin(fd),end(fd),begin(profile[0]),(double)0);
			matrix_w[seq_nr][m-(seq_end-seq_iter)] = sum_of_log_probs;
			++seq_iter;
		}
		log_probs_of_motifless_sequences[seq_nr] =
			(double)std::inner_product(begin(absLetterFrequencies),
									   end(absLetterFrequencies),
									   begin(profile[0]),(double)0);
		++dataset_start;
	}

	//
	double lambda_i = 0; 
	double sum_of_Q_i = 0;
	dataset_start -= t;
	while(dataset_start!=dataset_end)
	{
		TType seq_nr = t-(dataset_end-dataset_start);
		TType m = length(*dataset_start)-l+1;
		lambda_i = gamma/((double)m);
		joint_log_likelihood += (log((double)1)-log((double)m));

		double log_of_sums = matrix_w[seq_nr][0];
		TType j;
		for(j=1; j<m; ++j)
		{
			if( (matrix_w[seq_nr][j]-log_of_sums)>DBL_MIN_EXP )
			{
				log_of_sums += log(1+exp(matrix_w[seq_nr][j]-log_of_sums));
			}
		}
		for(j=0; j<m; ++j)
		{
			double prob = matrix_w[seq_nr][j];
			matrix_w[seq_nr][j] += log(lambda_i);
			double exponent =
				log_of_sums+log(lambda_i)-log_probs_of_motifless_sequences[seq_nr]-log(1-gamma);
			if( exponent>DBL_MIN_EXP )
			{
				matrix_w[seq_nr][j] -=
					(log_probs_of_motifless_sequences[seq_nr]+log(1-gamma)+log(1+exp(exponent)));
			}
			else
			{
				matrix_w[seq_nr][j] -= 
					(log_probs_of_motifless_sequences[seq_nr]+log(1-gamma));
			}
			matrix_w[seq_nr][j] = exp(matrix_w[seq_nr][j]);
			joint_log_likelihood += matrix_w[seq_nr][j]*prob;
		}
		double Q_i = 
			std::accumulate(matrix_w[seq_nr],matrix_w[seq_nr]+m,(double)0);
		sum_of_Q_i += Q_i;
		joint_log_likelihood +=
			((1-Q_i)*log_probs_of_motifless_sequences[seq_nr])
		   +(Q_i*log(lambda_i))
		   +((1-Q_i)*log(1-gamma));

		++dataset_start;
	}
	delete[] log_probs_of_motifless_sequences;

	// refine value of gamma
	gamma =  sum_of_Q_i*((double)1/(double)t);
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeStep_M:
..summary:Represents the M-step of the EM algorithm for Zoops models.
          Refines the background and motif component.
..cat:Motif Search
..signature:_computeStep_M(profile,dataset_start,matrix_w,t,l,zoops_model)
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.zoops_model:The zoops_model object.
...type:Tag.Zoops
..remarks:The joint log likelihood is computed in the E-step of the algorithm.
..include:seqan/find_motif.h
*/

template<typename TProfile, typename TIter, typename TMatrix, typename TType>
void
_computeMStep(TProfile & profile,
			  TIter dataset_start,
			  TMatrix const & matrix_w,
			  TType const & t,
	          TType const & l,
			  Zoops const & /*zoops*/)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Value<TFrequencyDist>::Type TValue;
	typedef typename Position<TProfile>::Type TPos;
	TPos k = 0;
	typedef String<TValue> TString;
	
	// total_counts_of_letters is used for the computtation of the background frequency
	TFrequencyDist total_counts_of_letters; //<-c
	absFreqOfLettersInSetOfSeqs(total_counts_of_letters,dataset_start, dataset_start+t);

	//
	for(k=0; k<length(profile); ++k)
	{
		profile[k] = TFrequencyDist();
	}

	// refine motif component
	TIter dataset_end = dataset_start+t;
	while(dataset_start!=dataset_end)
	{
		TType seq_nr = t-(dataset_end-dataset_start);
		TType m = length(*dataset_start)-l+1;
		typename Iterator<TString>::Type seq_iter, seq_end, ptr;
		seq_iter = begin(*dataset_start);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			ptr = seq_iter;
			for(TType i=0; i<l; ++i)
			{
				profile[i+1][(int)*(ptr+i)] += matrix_w[seq_nr][m-(seq_end-seq_iter)];
			}
			++seq_iter;
		}
		++dataset_start;
	}

	for(k=1; k<length(profile); ++k)
	{
		total_counts_of_letters -= profile[k];
	}

	// refine background component
	profile[0] = total_counts_of_letters;
	
	// addPseudocount (if necessary) & normalize
	double epsilon = 0.1;
	normalize(profile, Pseudocount<TValue, CMode>(epsilon));

	// matrix w - deallocate space (memory)
	for(TType pos=0; pos<t; ++pos)
	{
		delete[] matrix_w[pos];
	}
	delete[] matrix_w;
}

//////////////////////////////////////////////////////////////////////////////
//	Tcm model
//////////////////////////////////////////////////////////////////////////////

/*
..remarks:Dataset X is converted into a new pseudo-dataset consisting of all the width-l
          overlapping subsequences by running a window of width-l along each sequence Xi 
          and writing down the string contained in the window. Then the new dataset is 
          modeled as though it were generated by a two-component mixture model where each 
          component generates a string of width-l.
..remarks:Xij:=a width-l pseudo-sequence (=[Xij,Xij+1,...,Xij+l-1])
..remarks:lambda=gamma/m
*/

template<typename TProfile, typename TIter, typename TType>
double
em(TProfile & profile,
   TIter dataset_start,
   TType const & t,
   TType const & l,
   double & lambda,
   Tcm const & tcm)
{
	// matrix z - allocate space (memory)
	TType row_size = t;
	double ** matrix_w = new double*[row_size];
	for(TType pos=0; pos<row_size; ++pos)
	{
		TType col_size = length(*(dataset_start+pos))-l+1;
		matrix_w[pos] = new double[col_size];
	}

	// E-step: compute matrix w and the joint log likelihood
	double log_likelihood = 0;
	_computeEStep(matrix_w,log_likelihood,profile,dataset_start,lambda,t,l,tcm);

	// M-step: refine profile
	_computeMStep(profile,dataset_start,matrix_w,t,l,tcm);

	return log_likelihood;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeEStep:
..summary:Represents the E-step of the EM algorithm for Tcm models.
..cat:Motif Search
..signature:_computeEStep(matrix_w,joint_log_likelihood,profile,dataset_start,lambda,t,l,tcm_model)
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.joint_log_likelihood:The joint log-likelihood.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.lambda:The probability of starting a motif occurrence 
...type:$double$
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.tcm_model:The tcm_model object.
...type:Tag.Tcm
..remarks:The joint log likelihood is not computed in the M-step as MEME does, but rather in the E-step
          of the algorithm.
..remarks:The parameter lambda is reestimated during the E-step to save calculation time.
..include:seqan/find_motif.h
*/

template<typename TMatrix, typename TProfile, typename TIter, typename TType>
void
_computeEStep(TMatrix & matrix_w,
			  double & joint_log_likelihood,
			  TProfile & profile, 
			  TIter dataset_start,
			  double & lambda,
			  TType const & t,
			  TType const & l,
			  Tcm const & /*tcm*/)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Value<TFrequencyDist>::Type TValue;
	typedef String<TValue> TString;
	typedef typename Position<TProfile>::Type TPos;
	for(TPos pos=0; pos<length(profile); ++pos)
	{
		logarithmize(profile[pos]); 
	}

	//compute matrix w
	double prev_lambda = lambda;
	TIter dataset_end = dataset_start+t;
	while(dataset_start!=dataset_end)
	{
		TType seq_nr = t-(dataset_end-dataset_start);
		TType m = length(*dataset_start)-l+1;
		double Q_i = 0;
		typename Iterator<TString>::Type seq_iter, seq_end, ptr;
		seq_iter = begin(*dataset_start);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			TType seq_pos = m-(seq_end-seq_iter);
			ptr = seq_iter;
			double log_prob_given_theta0 = 0;
			double log_prob_given_theta1 = 0;
			for(TType i=0; i<l; ++i)
			{
				log_prob_given_theta0 += profile[0][(int)*(ptr+i)];
				log_prob_given_theta1 += profile[i+1][(int)*(ptr+i)];
			}
			double minuend = log_prob_given_theta1+log(prev_lambda);
			double subtrahend = 0;
			double exponent = minuend-log_prob_given_theta0-log((double)(1-prev_lambda));
			if(exponent>DBL_MIN_EXP)
			{
				subtrahend =
					log_prob_given_theta0+log((double)(1-prev_lambda))
				   +log((double)(1+exp(exponent)));
			}
			else
			{
				subtrahend =
					log_prob_given_theta0+log((double)(1-prev_lambda));
			}

			matrix_w[seq_nr][seq_pos] = exp(minuend-subtrahend);
			Q_i += matrix_w[seq_nr][seq_pos];
			joint_log_likelihood += ((double)(1-matrix_w[seq_nr][seq_pos]))*log_prob_given_theta0
									+matrix_w[seq_nr][seq_pos]*log_prob_given_theta1
									+((double)(1-matrix_w[seq_nr][seq_pos]))*log((double)(1-prev_lambda))
									+matrix_w[seq_nr][seq_pos]*log(prev_lambda);
			++seq_iter;
		}
		lambda += (Q_i/((double)m));
		++dataset_start;
	}
	lambda = lambda/((double)t);

	// apply a smoothing step to reduce the degree to which any two overlapping 
	// subsequences can both be assigned to the motif component
	dataset_start -= t;
	_smoothingStep(matrix_w, dataset_start, t, l);
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeStep_M:
..summary:Represents the M-step of the EM algorithm for Tcm models.
          Refines the background and motif component.
..cat:Motif Search
..signature:_computeStep_M(profile,dataset_start,matrix_w,t,l,tcm_model)
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.dataset_start:n iterator pointing to the first input sequence of a given dataset.
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.tcm_model:The tcm_model object.
...type:Tag.Tcm
..remarks:The joint log likelihood is computed in the E-step of the algorithm.
..include:seqan/find_motif.h
*/

template<typename TProfile, typename TIter, typename TMatrix, typename TType>
void
_computeMStep(TProfile & profile,
			  TIter dataset_start,
			  TMatrix const & matrix_w,
			  TType const & t,
	          TType const & l,
			  Tcm const & /*tcm*/)
{
	typedef typename Value<TProfile>::Type TFrequencyDist;
	typedef typename Value<TFrequencyDist>::Type TValue;
	typedef typename Position<TProfile>::Type TPos;
	TPos k = 0;
	typedef String<TValue> TString;
	
	//
	for(k=0; k<length(profile); ++k)
	{
		profile[k] = TFrequencyDist();
	}

	// refine motif component
	TIter dataset_end = dataset_start+t;
	while(dataset_start!=dataset_end)
	{
		TType seq_nr = t-(dataset_end-dataset_start);
		TType m = length(*dataset_start)-l+1;
		typename Iterator<TString>::Type seq_iter, seq_end, ptr;
		seq_iter = begin(*dataset_start);
		seq_end = seq_iter+m;
		while(seq_iter!=seq_end)
		{
			ptr = seq_iter;
			for(TType i=0; i<l; ++i)
			{
				profile[i+1][(int)*(ptr+i)] += matrix_w[seq_nr][m-(seq_end-seq_iter)];
				profile[0][(int)*(ptr+i)] += (double)(1-matrix_w[seq_nr][m-(seq_end-seq_iter)]);
			}
			++seq_iter;
		}
		++dataset_start;
	}

	// addPseudocount (if necessary) & normalize
	double epsilon = 0.1;
	normalize(profile, Pseudocount<TValue, CMode>(epsilon));

	// matrix w - deallocate space (memory)
	for(TType pos=0; pos<t; ++pos)
	{
		delete[] matrix_w[pos];
	}
	delete[] matrix_w;
}

//////////////////////////////////////////////////////////////////////////////
//Subfunctions
//////////////////////////////////////////////////////////////////////////////

/*
.Function._smoothingStep:
..summary:Applies a smoothing step to reduce the degree to which any two overlapping 
          subsequences can both be assigned to the motif component.
         (We do not want the model to predict that two overlapping substrings are both 
          motif occurrences.)
..cat:Motif Search
..signature:_smoothingStep(matrix_w,dataset_start,t,l)
..param.matrix_w:The matrix_w object.
...remarks:w(i,j) is the probability that the motif pattern 
            starts at position j in the i-th input sequence.
..param.dataset_start:An iterator pointing to the first input sequence of a given dataset.
..param.t:The number of input sequences.
..param.l:The size of the motif.
..remarks:Function is used by the EM algorithm for Tcm models.
..include:seqan/find_motif.h
*/

template<typename TMatrix, typename TIter, typename TType>
void
_smoothingStep(TMatrix & matrix_w,
			   TIter dataset_start,
			   TType const & t,
			   TType const & l)
{
	TType h = 0;
	for(TType i=0; i<t; ++i)
	{
		TType m = length(*(dataset_start+i))-l+1;
		for(TType s=0; s<l; ++s)
		{
			TType k = s;
			double w_sum;
			double w_max;
			double factor;
			while( ((m-k)/l)>=1 )
			{
				w_sum = 0;
				w_max = 0;
				factor = 0;

				// w_sum:=sum of Zij in the current window of size l (l:=l)
				w_sum = (double)(std::accumulate(matrix_w[i]+k, 
										matrix_w[i]+k+l, (double)0));
				// w_max:=the largest Zij in the current window
				w_max =
					(double)(*std::max_element(matrix_w[i]+k,matrix_w[i]+k+l));
				if(w_sum>1.0)
				{
					factor = ((double)1-w_max)/(w_sum-w_max);
					for(h=k; h<k+l; ++h)
					{
						if(matrix_w[i][h]!=w_max)
						{
							matrix_w[i][h] *= factor;
						}
					}
				}
				k += l;
			}

			if( (m-k)>1 )
			{
				w_sum = 0;
				w_max = 0;
				factor = 0;

				// w_sum:=sum of Zij in the current window of size l (l:=l)
				w_sum = (double)(std::accumulate(matrix_w[i]+k, 
										matrix_w[i]+m, (double)0));
				// w_max:=the largest Zij in the current window
				w_max =
					(double)(*std::max_element(matrix_w[i]+k,matrix_w[i]+m));
				if(w_sum>1.0)
				{
					factor = ((double)1-w_max)/(w_sum-w_max);
					for(h=k; h<m; ++h)
					{
						if(matrix_w[i][h]!=w_max)
						{
							matrix_w[i][h] *= factor;
						}
					}
				}
			}
		}
	}
}

} // SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
