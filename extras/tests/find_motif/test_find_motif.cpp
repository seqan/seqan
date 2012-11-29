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

#include <iostream>

#include <seqan/find_motif.h>

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template<typename TIter, typename TString, typename TType>
bool
isOOPSMotif(TIter ds_iter, TIter ds_end, TString const & motif, TType const & d, bool & is_exact)
{
	bool isRelevantMotif = false;
	int t = (int) (ds_end-ds_iter);
	int counter = 0;
	int hd = 0;
	typename Size<TString>::Type l = length(motif);
	typename Iterator<TString>::Type seq_iter, seq_end, l_mer_begin;
	while( (ds_iter!=ds_end) & (((int)(ds_end-ds_iter))==t) )
	{
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+length(*ds_iter)-l+1;
		counter = 0;
		while( (seq_iter!=seq_end) & (counter<2) )
		{
			l_mer_begin = begin(motif);
			hd = hammingDistance<int>(seq_iter, seq_iter+l, l_mer_begin);
			if( is_exact & (hd == (int) d) )
			{
				++counter;
			}
			else if( !is_exact & (hd<=(int) d))
			{
				++counter;
			}
			++seq_iter;
		}
		if(counter==1)
		{
			--t;
		}
		++ds_iter;
	}
	if(t==0)
	{
		isRelevantMotif = true;
	}
	return isRelevantMotif;
}

template<typename TIter, typename TString, typename TType>
bool
isOMOPSMotif(TIter ds_iter, TIter ds_end, TString const & motif, TType const & d, bool & is_exact)
{
	bool isRelevantMotif = false;
	int t = (int) (ds_end-ds_iter);
	int counter = 0;
	int hd = 0;
	typename Size<TString>::Type l = length(motif);
	typename Iterator<TString>::Type seq_iter, seq_end, l_mer_begin;
	while( (ds_iter!=ds_end) & (((int)(ds_end-ds_iter))==t) )
	{
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+length(*ds_iter)-l+1;
		counter = 0;
		while( (seq_iter!=seq_end) & (counter<1) )
		{
			l_mer_begin = begin(motif);
			hd = hammingDistance<int>(seq_iter, seq_iter+l, l_mer_begin);
			if( is_exact & (hd== (int) d) )
			{
				++counter;
			}
			else if( !is_exact & (hd<=(int) d))
			{
				++counter;
			}
			++seq_iter;
		}
		if(counter==1)
		{
			--t;
		}
		++ds_iter;
	}
	if(t==0)
	{
		isRelevantMotif = true;
	}
	return isRelevantMotif;
}

template<typename TIter, typename TString, typename TType>
bool
isZOOPSMotif(TIter ds_iter, TIter ds_end, TString const & motif, TType const & d, bool & is_exact)
{
	bool isRelevantMotif = false;
	int t = (int) (ds_end-ds_iter);
	int upper_limit = t-((int) floor(t*(Zoops().threshold)+0.5)-1);
	int counter = 0;
	int hd = 0;
	typename Size<TString>::Type l = length(motif);
	typename Iterator<TString>::Type seq_iter, seq_end, l_mer_begin;
	while( (ds_iter!=ds_end) )
	{
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+length(*ds_iter)-l+1;
		counter = 0;
		while( (seq_iter!=seq_end) & (counter<2) )
		{
			l_mer_begin = begin(motif);
			hd = hammingDistance<int>(seq_iter, seq_iter+l, l_mer_begin);
			if( is_exact & (hd== (int) d) )
			{
				++counter;
			}
			else if( !is_exact & (hd<= (int) d))
			{
				++counter;
			}
			++seq_iter;
		}
		if(counter==1)
		{
			--t;
		}
		++ds_iter;
	}
	if(t<=upper_limit)
	{
		isRelevantMotif = true;
	}
	return isRelevantMotif;
}

template<typename TIter, typename TString, typename TType>
bool
isTCMMotif(TIter ds_iter, TIter ds_end, TString const & motif, TType const & d, bool & is_exact)
{
	bool isRelevantMotif = false;
	int t = (int) (ds_end-ds_iter);
	int upper_limit = t-((int) floor(t*(Tcm().threshold)+0.5)-1);
	int counter = 0;
	int hd = 0;
	typename Size<TString>::Type l = length(motif);
	typename Iterator<TString>::Type seq_iter, seq_end, l_mer_begin;
	while( (ds_iter!=ds_end) )
	{
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+length(*ds_iter)-l+1;
		counter = 0;
		while( (seq_iter!=seq_end) & (counter<1) )
		{
			l_mer_begin = begin(motif);
			hd = hammingDistance<int>(seq_iter, seq_iter+l, l_mer_begin);
			if( is_exact & (hd==(int) d) )
			{
				++counter;
			}
			else if( !is_exact & (hd<=(int) d))
			{
				++counter;
			}
			++seq_iter;
		}
		if(counter==1)
		{
			--t;
		}
		++ds_iter;
	}
	if(t<=upper_limit)
	{
		isRelevantMotif = true;
	}
	return isRelevantMotif;
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_find_motif_exact_algorithms) {
	//Testing Pms1 & Pmsp algorithm

	unsigned int t = 0;      //number of sequences
	unsigned int l = 0;		//length of motif
	unsigned int d = 0;		//number of substitutions
	bool is_exact = false;	//size of Hamming distance
	unsigned int i = 0;

    // Initialize random number generator.
    typedef typename GetDefaultRng<MotifFinderClass>::Type TRng;
    TRng & rng = defaultRng(MotifFinderClass());
    reSeed(rng, 0);

//____________________________________________________________________________
// Test1 - Search for Oops motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	t = 3;
	l = 4;		
	d = 1;		
	is_exact = true;	

	String<DnaString> dataset1;
	resize(dataset1, t);
	dataset1[0] = "ACAGCA";
	dataset1[1] = "AGGCAG";
	dataset1[2] = "TCAGTC";

	//Application of Pms1-Oops
	MotifFinder<Dna, Pms1> motif_finder1(l,d,is_exact);
	findMotif(motif_finder1,dataset1,Oops());

	//Application of Pmsp-Oops
	MotifFinder<Dna, Pmsp> motif_finder2(l,d,is_exact);
	findMotif(motif_finder2,dataset1,Oops());

	SEQAN_ASSERT(length(motif_finder1.set_of_motifs)==length(motif_finder2.set_of_motifs));
	for(i=0; i<length(motif_finder1.set_of_motifs); ++i)
	{
		SEQAN_ASSERT(motif_finder1.set_of_motifs[i]==motif_finder2.set_of_motifs[i]);
	}

//____________________________________________________________________________
// Test2 - Search for Omops motifs on a small set of nucleotide sequences
//         given the inexact Hamming distance (<=d)

	l = 6;		
	d = 2;		
	is_exact = false;	

	String<DnaString> dataset2;
	appendValue(dataset2,DnaString("GCTGGACGTG"));
	appendValue(dataset2,DnaString("TCTAGACATA"));
	appendValue(dataset2,DnaString("AGTGGGGGAC"));
	appendValue(dataset2,DnaString("CTAGTCAAGA"));
	appendValue(dataset2,DnaString("CTCGAGGGGT"));

	//Application of Pms1-Omops
	MotifFinder<Dna, Pms1> motif_finder3(l,d,is_exact);
	findMotif(motif_finder3,dataset2,Omops());

	//Application of Pmsp-Omops
	MotifFinder<Dna, Pmsp> motif_finder4(l,d,is_exact);
	findMotif(motif_finder4,dataset2,Omops());

	SEQAN_ASSERT(length(motif_finder3.set_of_motifs)==length(motif_finder4.set_of_motifs));
	for(i=0; i<length(motif_finder3.set_of_motifs); ++i)
	{
		SEQAN_ASSERT(motif_finder3.set_of_motifs[i]==motif_finder4.set_of_motifs[i]);
	}

//____________________________________________________________________________
// Test3 - Search for Zoops motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	l = 6;		
	d = 1;		
	is_exact = true;	

	String<DnaString> dataset3;
	appendValue(dataset3,DnaString("AGCCGTCTGA"));
	appendValue(dataset3,DnaString("TCCAGGCAAG"));
	appendValue(dataset3,DnaString("GAACGTCCAA"));
	appendValue(dataset3,DnaString("GCTTTCTAAC"));
	appendValue(dataset3,DnaString("AGTAGCTCGC"));

	//Application of Pms1-Zoops
	MotifFinder<Dna, Pms1> motif_finder5(l,d,is_exact);
	findMotif(motif_finder5,dataset3,Zoops());

	//Application of Pmsp-Zoops
	MotifFinder<Dna, Pmsp> motif_finder6(l,d,is_exact);
	findMotif(motif_finder6,dataset3,Zoops());

	SEQAN_ASSERT(length(motif_finder5.set_of_motifs)==length(motif_finder6.set_of_motifs));
	for(i=0; i<length(motif_finder5.set_of_motifs); ++i)
	{
		SEQAN_ASSERT(motif_finder5.set_of_motifs[i]==motif_finder6.set_of_motifs[i]);
	}

//____________________________________________________________________________
// Test4 - Search for Tcm motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	l = 5;		
	d = 1;		

	String<DnaString> dataset4;
	appendValue(dataset4,DnaString("CAATTAACTC"));
	appendValue(dataset4,DnaString("ATAAACAGTG"));
	appendValue(dataset4,DnaString("GAATGCATTG"));

	//Application of Pms1-Tcm
	MotifFinder<Dna, Pms1> motif_finder7(l,d,is_exact);
	findMotif(motif_finder7,dataset4,Tcm());

	//Application of Pmsp-Tcm
	MotifFinder<Dna, Pmsp> motif_finder8(l,d,is_exact);
	findMotif(motif_finder8,dataset4,Tcm());
	
	SEQAN_ASSERT(length(motif_finder7.set_of_motifs)==length(motif_finder8.set_of_motifs));
	for(i=0; i<length(motif_finder7.set_of_motifs); ++i)
	{
		SEQAN_ASSERT(motif_finder7.set_of_motifs[i]==motif_finder8.set_of_motifs[i]);
	}

}


SEQAN_DEFINE_TEST(test_find_motif_approximation_algorithms) {
	//Testing Projection & ePatternBranching algorithm
	
	unsigned int t = 0;      //number of sequences
	unsigned int n = 0;		//length of sequence
	unsigned int l = 0;		//length of motif
	unsigned int d = 0;		//number of substitutions
	bool is_exact = false;	//size of Hamming distance
	unsigned int m =0;		//total number of possible l-mers
	unsigned int h = 0;		//size of the neighborhood considering at first
	unsigned int i = 0;

    // Initialize random number generator.
    typedef typename GetDefaultRng<MotifFinderClass>::Type TRng;
    TRng & rng = defaultRng(MotifFinderClass());
    reSeed(rng, 0);

//____________________________________________________________________________
// Test1 - Search for Oops motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	t = 3;  
	n = 6;		
	l = 4;		
	d = 1;		
	is_exact = true;	
	m = t*(n-l+1);

	String<DnaString> dataset1;
	appendValue(dataset1,DnaString("ACAGCA"));
	appendValue(dataset1,DnaString("AGGCAG"));
	appendValue(dataset1,DnaString("TCAGTC"));

	//Application of PROJECTION-Oops
    MotifFinder<Dna, Projection> motif_finder1(t,l,m,d,is_exact);
	findMotif(motif_finder1, dataset1, Oops());
	//check whether found motif is really an Oops motif
	SEQAN_ASSERT(isOOPSMotif(begin(dataset1),
				  end(dataset1),
				  getMotif(motif_finder1,0),
				  d,
				  is_exact)==true);

//____________________________________________________________________________
//
	//Application of ePatternBranching-Oops
    MotifFinder<Dna, EPatternBranching> motif_finder2(t,l,d,is_exact,h);
	findMotif(motif_finder2, dataset1, Oops());
	//check whether found motif is really an Oops motif
	for(i=0; i<length(motif_finder2.set_of_motifs); ++i)
	{
		SEQAN_ASSERT(isOOPSMotif(begin(dataset1),
					  end(dataset1),
				      motif_finder2.set_of_motifs[i],
				      d,
				      is_exact)==true);
	}

//____________________________________________________________________________
// Test2 - Search for Omops motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (<=d)

	is_exact = false;

	//Application of PROJECTION-Omops
    MotifFinder<Dna, Projection> motif_finder3(t,l,m,d,is_exact);
	findMotif(motif_finder3, dataset1, Omops());
	//check whether found motif is really an Omops motif
	SEQAN_ASSERT(isOMOPSMotif(begin(dataset1),
				  end(dataset1),
				  getMotif(motif_finder3, 0),
				  d,
				  
				  is_exact)==true);

//____________________________________________________________________________
//
	//Application of ePatternBranching-Omops
	MotifFinder<Dna, EPatternBranching> motif_finder4(t,l,d,is_exact,h);
	findMotif(motif_finder4, dataset1, Omops());
	//check whether found motif is really an Omops motif
	for(i=0; i<length(motif_finder4.set_of_motifs); ++i)
	{
		SEQAN_ASSERT(isOMOPSMotif(begin(dataset1),
					  end(dataset1),
				      motif_finder4.set_of_motifs[i],
				      d,
				      is_exact)==true);
	}

//____________________________________________________________________________
// Test3 - Search for Zoops motifs on a set of small nucleotide sequences
//         given the inexact Hamming distance (<=d)

	//Application of PROJECTION-Zoops
    MotifFinder<Dna, Projection> motif_finder5(t,l,m,d,is_exact);
	findMotif(motif_finder5, dataset1, Zoops());
	//check whether found motif is really a Zoops motif
	SEQAN_ASSERT(isZOOPSMotif(begin(dataset1),
				  end(dataset1),
				  getMotif(motif_finder5, 0),
				  d,
				  is_exact)==true);

//____________________________________________________________________________
// Test4 - Search for Tcm motifs on a set of small nucleotide sequences
//         given the exact Hamming distance (=d)

	is_exact = true;

	//Application of PROJECTION-Tcm
    MotifFinder<Dna, Projection> motif_finder6(t,l,m,d,is_exact);
	findMotif(motif_finder6, dataset1, Tcm());
	//check whether found motif is really a Tcm motif
	SEQAN_ASSERT(isTCMMotif(begin(dataset1),
				  end(dataset1),
				  getMotif(motif_finder6),
				  d,
				  is_exact)==true);
}

SEQAN_BEGIN_TESTSUITE(test_find_motif)
{
    SEQAN_CALL_TEST(test_find_motif_exact_algorithms);
    SEQAN_CALL_TEST(test_find_motif_approximation_algorithms);
}
SEQAN_END_TESTSUITE
