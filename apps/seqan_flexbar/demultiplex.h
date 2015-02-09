// ==========================================================================
//                               demultiplex.h
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// ==========================================================================
// This file provides the demultiplexing functionality of seqan-flexbar
// which is based in the implementation of the original flexbar program
// in [1].
// [1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, C.  FLEXBAR—Flexible
// Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
// Biology 2012, 1, 895-905.
// ==========================================================================


#ifndef SANDBOX_GROUP3_APPS_SEQDPT_DEMULTIPLEX_H_
#define SANDBOX_GROUP3_APPS_SEQDPT_DEMULTIPLEX_H_

#ifdef _OPENMP
#include <omp.h>
#endif
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include "general_processing.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef String<Dna5Q> TAlphabet;

struct DemultiplexStats
{
	String<unsigned> groups;
};


// ============================================================================
// Functions
// ============================================================================

template <typename TSeqs, typename TIds, typename TBarcodes> //This version works with paired-end data
bool check(TSeqs& seqs,  TIds& ids, TSeqs& seqsRev, TIds& idsRev, TBarcodes& barcodes, GeneralStats& stats)
{
	unsigned len = length(barcodes[0]);
	for (unsigned i = 1; i < length(barcodes); ++i)
	{
		if (len != length(barcodes[i]))
		{
			std::cerr << "ERROR: Barcodes differ in length. All barcodes must be of equal length.\n";
			return false;
		}
	}   //Iterating backward to avoid error after deletion of a sequence
    unsigned ex = 0;
	for (int i = length(seqs) - 1; i >= 0; --i)  
	{   //integer necessary because unsigned would cause error after last iteration
		if (length(seqs[i]) <= len)	
		{
            //Divinding the comming jobs between the threads
            //sorting all  sequences which shall be deleted to the end of the container
            swap(seqs[i], seqs[len - ex - 1]);
            swap(seqsRev[i], seqsRev[len - ex - 1]);
            swap(ids[i], ids[len - ex - 1]);
            swap(idsRev[i], idsRev[len - ex - 1]);
            ++ex;
		}
	}
    if (ex != 0)
    {
        resize(seqs, len - ex);
        resize(seqsRev, len - ex);
        resize(ids, len - ex);
        resize(idsRev, len - ex);
        stats.removedSeqsShort += (2 * ex);
    }
	return true;
}
//Overload for single-end data
template <typename TSeqs, typename TIds, typename TBarcodes>
bool check(TSeqs& seqs, TIds& ids,TBarcodes& barcodes, GeneralStats& stats) 
{
	unsigned len = length(barcodes[0]);
	for (unsigned i = 1; i < length(barcodes); ++i)
	{
		if (len != length(barcodes[i]))
		{
			std::cerr << "ERROR: Barcodes differ in length. All barcodes must be of equal length.\n";
			return false;
		}
	} //Iterating backward to avoid error after deletion of a sequence
    unsigned ex = 0;
    unsigned lenSeq = length(seqs);
    for (int i = lenSeq - 1; i >= 0; --i)
	{
		if (length(seqs[i]) <= len)
		{
            //sorting all  sequences which shall be deleted to the end of the container
            swap(seqs[i], seqs[lenSeq - ex - 1]);
            swap(ids[i], ids[lenSeq - ex - 1]);
            ++ex;
		}
    }
    if (ex != 0)
    {
        resize(seqs, lenSeq - ex);
        resize(ids, lenSeq - ex);
        stats.removedSeqsShort += ex;
    }
	return true;
}

template <typename TSeqs>
void getPrefix(TSeqs& prefices, TSeqs& seqs, unsigned len)
{
	int limit = length(seqs);
	resize(prefices, limit);
	SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static))
	for (int i = 0; i < limit; ++i) //Remark: OMP requires an integer as loop-variable
	{
		prefices[i] = prefix(seqs[i], len);
	}
}

template <typename TBarcode>
void buildVariations(StringSet<Dna5String>& variations, const TBarcode& barcode)
{
	int limit = (length(barcode))*5;		    //possible number of variations with one error (A,T,G,C,N)
	resize(variations, limit);				    //resizes according to calculated number of variations
	SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static))
	for (int i = 0; i < limit; ++i)
	{
        assign(variations[i], barcode);	        //fills resultset with original barcode
	}
	limit = limit/5;
	SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static))
	for (int i = 0; i < limit; ++i)
	{										    //each run through the loop modifies 1 position of the barcode 4 times
		move(variations[5*i][i], 'A');	    	//multiplication with 5 and addition of 0...4 prevents index conflicts
		move(variations[5*i+1][i], 'C');
		move(variations[5*i+2][i], 'G');
		move(variations[5*i+3][i], 'T');
		move(variations[5*i+4][i], 'N');
	}
}

template <typename TBarcodes>
void buildAllVariations(TBarcodes& barcodes)
{
	StringSet<Dna5String> newbarcodes;			//stores the new barcodes
	for (unsigned i = 0; i < length(barcodes); ++i)
	{
		StringSet<Dna5String> tempbarcodes;
        buildVariations(tempbarcodes, barcodes[i]);
		for (unsigned j = 0; j < length(tempbarcodes); ++j)
        {
			appendValue(newbarcodes, tempbarcodes[j]);
        }
	}
	clear(barcodes);
	barcodes = newbarcodes;
}

template <typename TPrefix, typename TFinder>
int findExactIndex(const TPrefix& prefix, TFinder& finder)
{
	clear(finder);								//resets finder
	if (find(finder, prefix))
    {
		return getSeqNo(position(finder));		//returns index of barcode. ONLY THE FIRST HIT!
    }
	else return -1;								//return -1 if no hit occured
}

template <typename TPrefices, typename TFinder, typename TStats>
void findAllExactIndex(String<int>& matches, const TPrefices& prefices, const TFinder& finder, TStats& stats)
{
    resize(matches, length(prefices));

	int tnum = 1;
#ifdef _OPENMP
    tnum = omp_get_max_threads();
#endif
    StringSet<TFinder> finderSet;
    resize(finderSet, tnum);            //creating a set of finders to prevent their continous re-initialisation
    for (int i = 0; i < tnum; ++i)
    {
        finderSet[i] = finder;
    }
	int limit = length(prefices);
	SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static))
	for (int i = 0; i < limit; ++i)
	{
		int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        int hit = findExactIndex(prefices[i], finderSet[tid]);
		matches[i] = hit;
	}
    for (unsigned i = 0; i < length(matches); ++i) //outside of parallel loop to avoid use of "atomic"
    {
        ++stats.groups[matches[i]+1];
    }
}

template <typename TSeqs>
void clipBarcodes(TSeqs& seqs, const String<int>& matches, unsigned len)
{
	int limit = length(matches);
	SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static))
	for (int i = 0; i < limit; ++i)
		if (matches[i] != -1)					//only erases barcode from sequence if it could be matched
        {
			erase(seqs[i], 0 , len);
        }
}

//Overload for deleting the barcodes in any case.
template <typename TSeqs>
void clipBarcodes(TSeqs& seqs, int len)
{
	int limit = length(seqs);
	SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static))
	for (int i = 0; i < limit; ++i)
    {
		erase(seqs[i], 0 , len);
    }
}

template <typename TMatches, typename TBarcodes>
void group(StringSet<String<int> >& sortedSequences, const TMatches& matches, const TBarcodes& barcodes, bool exclude)
{
	resize(sortedSequences, length(barcodes)+1);
	for (unsigned i = 0; i < length(matches); ++i)
    {                                                       //adds index of sequence to respective group.
		if ((!exclude) || (matches[i] != -1))                //Check if unidentified seqs have to be excluded
        {                                                   //offset by 1 is necessary since group 0 is...
            appendValue(sortedSequences[matches[i]+1], i);  //...reserved for unidentified sequences)
        }
    }                                                  
}

//Overload if approximate search has been used.
template <typename TMatches, typename TBarcodes, typename TApprox>
void group(StringSet<String<int> >& sortedSequences, const TMatches& matches,
    const TBarcodes& barcodes, TApprox const &, bool exclude)
{
	resize(sortedSequences, length(barcodes)/5+1);
	float dividend = float(length(barcodes[0])*5.0);		//value by which the index will be corrected.
	for (unsigned i = 0; i < length(matches); ++i)			//adds index of sequence to respective group.
    {
		if ((!exclude) || (matches[i] != -1))                //Check if unidentified seqs have to be excluded
        {
            appendValue(sortedSequences[int(floor(float(matches[i])/dividend))+1], i);
        }
    }
}

//Using exact search and multiplex barcodes.
template<typename TBarcodes, typename TMultiplex ,typename TFinder>
void doAll(StringSet<String<int> >& sortedSequences, TMultiplex& multiplex, TBarcodes& barcodes,
    TFinder& esaFinder, DemultiplexStats& stats, bool exclude)
{
	String<int> matches;
	findAllExactIndex(matches, multiplex, esaFinder, stats);
    group(sortedSequences, matches, barcodes, exclude);
}
// Using approximate search and multiplex barcodes.
template<typename TBarcodes, typename TMultiplex ,typename TFinder, typename TApprox>
void doAll(StringSet<String<int> >& sortedSequences, TMultiplex& multiplex, TBarcodes& barcodes, TFinder& esaFinder,
    DemultiplexStats& stats, const TApprox approximate, bool exclude)
{
	String<int> matches;
	findAllExactIndex(matches, multiplex, esaFinder, stats);
    group(sortedSequences, matches, barcodes, approximate, exclude);
}
//Using exact search and inline barcodes.
template<typename TSeqs, typename TBarcodes, typename TFinder>
void doAll(StringSet<String<int> >& sortedSequences, TSeqs& seqs, TBarcodes& barcodes, TFinder& esaFinder, bool hardClip,
    DemultiplexStats& stats, bool exclude)
{
	TSeqs prefices;
    getPrefix(prefices, seqs, length(barcodes[0]));
	String<int> matches;
	findAllExactIndex(matches, prefices, esaFinder, stats);
	if (hardClip)		//clip barcodes according to selected method
	{
		clipBarcodes(seqs, length(barcodes[0]));
	}
	else
	{
        clipBarcodes(seqs, matches, length(barcodes[0]));
	}
    group(sortedSequences, matches, barcodes, exclude);
}
// Using approximate search and inline barcodes.
template<typename TSeqs, typename TBarcodes, typename TFinder, typename TApprox>
void doAll(StringSet<String<int> >& sortedSequences, TSeqs& seqs, TBarcodes& barcodes, TFinder& esaFinder,
    bool hardClip, DemultiplexStats& stats, const TApprox approximate, bool exclude)
{
	TSeqs prefices;
    getPrefix(prefices, seqs, length(barcodes[0]));
	String<int> matches;
	findAllExactIndex(matches, prefices, esaFinder, stats);
	if (hardClip)		//clip barcodes according to selected method
	{
	    clipBarcodes(seqs, length(barcodes[0]));
	}
	else
	{
        clipBarcodes(seqs, matches, length(barcodes[0]));
	}
    group(sortedSequences, matches, barcodes, approximate, exclude);
}

//Version for paired-end data
template<typename TSeqs, typename TIds>
void buildSets(TSeqs& seqs, TSeqs& seqsRev, TIds& ids, TIds& idsRev, const StringSet<String<int> >& groups,
		String<TSeqs>& gSeqs, String<TSeqs>& gSeqsRev, String<TIds>& gIds, String<TIds>& gIdsRev)
{
	unsigned len = length(groups);
    resize(gSeqs, len);
	resize(gSeqsRev, len);
	resize(gIds, len);
	resize(gIdsRev, len);
    unsigned k = 0;
	for (unsigned i = 0; i < length(groups); ++i)
	{
		for (unsigned j = 0; j < length(groups[i]); ++j)
		{
			appendValue(gSeqs[k], seqs[groups[i][j]]);
			appendValue(gSeqsRev[k], seqsRev[groups[i][j]]);
			appendValue(gIds[k], ids[groups[i][j]]);
			appendValue(gIdsRev[k], idsRev[groups[i][j]]);
		}
		if (length(groups[i]) != 0)
		{
			++k;
		}
	}
	resize(gSeqs, k);
	resize(gSeqsRev, k);
	resize(gIds, k);
	resize(gIdsRev, k);
	clear(seqs);
	clear(seqsRev);
	clear(ids);
	clear(idsRev);
}
//Overload for single-end data.
template<typename TSeqs, typename TIds>
void buildSets(TSeqs& seqs, TIds& ids, const StringSet<String<int> >& groups, String<TSeqs>& gSeqs, String<TIds>& gIds)
{
	resize(gSeqs, length(groups));
	resize(gIds, length(groups));
	unsigned k = 0;
	for (unsigned i = 0; i < length(groups); ++i)
	{
		for (unsigned j = 0; j < length(groups[i]); ++j)
		{
            appendValue(gSeqs[k], seqs[groups[i][j]]);
		    appendValue(gIds[k],ids[groups[i][j]]);
		}
		if (length(groups[i]) != 0)
		{
			++k;
		}
	}
    resize(gSeqs, k);
	resize(gIds, k);
	clear(seqs);
 	clear(ids);
}
#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_DEMULTIPLEX_H_
