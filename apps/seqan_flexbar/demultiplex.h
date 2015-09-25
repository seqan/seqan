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

#include "helper_functions.h"


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef seqan::String<seqan::Dna5Q> TAlphabet;

struct DemultiplexStats
{
	seqan::String<unsigned> groups;
};

struct DemultiplexingParams
{
	seqan::String<char> barcodeFile;
	seqan::StringSet<seqan::String<seqan::Dna5> > barcodes;
	seqan::StringSet<seqan::String<char> > barcodeIds;
	seqan::String<char> multiplexFile;
	bool approximate;
	bool hardClip;
	bool run;
	bool runx;
	bool exclude;
	DemultiplexStats stats;

	DemultiplexingParams() :
		approximate(false),
		hardClip(false),
		run(false),
		runx(false),
		exclude(false)
	{};
};

// ============================================================================
// Functions
// ============================================================================

template <typename TReads, typename TBarcodes, typename TStats>
bool check(TReads& reads, TBarcodes& barcodes, TStats& stats) noexcept
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
    auto it = std::remove_if(reads.begin(), reads.end(), [len](auto& read) {return length(read.seq) <= len;});
    stats.removedSeqsShort += std::distance(it, reads.end());
    reads.erase(it, reads.end());
    return true;
}

template <template <typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>,Read<TSeq>>::value || std::is_same<TRead<TSeq>,ReadPairedEnd<TSeq>>::value>>
void getPrefix(std::vector<TSeq>& prefices, std::vector<TRead<TSeq>>& reads, unsigned len)
{
    int limit = reads.size();
    assert(prefices.size() == reads.size());
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static))
        for (int i = 0; i < limit; ++i) //Remark: OMP requires an integer as loop-variable
        {
            prefices[i] = prefix(reads[i].seq, len);
        }
}

template <template <typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value>>
void getPrefix(std::vector<TSeq>& prefices, std::vector<TRead<TSeq>>& reads, unsigned len, bool = false)
{
    (void)len;
    int limit = reads.size();
    prefices.resize(limit);
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static))
        for (int i = 0; i < limit; ++i) //Remark: OMP requires an integer as loop-variable
        {
            prefices[i] = reads[i].demultiplex;
        }
}

template <typename TBarcode>
void buildVariations(seqan::StringSet<seqan::Dna5String>& variations, const TBarcode& barcode)
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
	seqan::StringSet<seqan::Dna5String> newbarcodes;			//stores the new barcodes
	for (unsigned i = 0; i < length(barcodes); ++i)
	{
		seqan::StringSet<seqan::Dna5String> tempbarcodes;
        buildVariations(tempbarcodes, barcodes[i]);
		for (unsigned j = 0; j < length(tempbarcodes); ++j)
        {
			seqan::appendValue(newbarcodes, tempbarcodes[j]);
        }
	}
	clear(barcodes);
	barcodes = newbarcodes;
}

template <typename TPrefix, typename TFinder>
int findExactIndex(const TPrefix& prefix, TFinder& finder) noexcept
{
	clear(finder);								//resets finder
	if (find(finder, prefix))
    {
		return getSeqNo(position(finder));		//returns index of barcode. ONLY THE FIRST HIT!
    }
	else return -1;								//return -1 if no hit occured
}

template <typename TMatches, typename TPrefices, typename TFinder, typename TStats>
void findAllExactIndex(TMatches& matches, const TPrefices& prefices, const TFinder& finder, TStats& stats) noexcept
{
    assert(length(matches) == length(prefices));

	int tnum = 1;
#ifdef _OPENMP
    tnum = omp_get_max_threads();
#endif
    seqan::StringSet<TFinder> finderSet;
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

template <typename TRead>
void clipBarcodes(std::vector<TRead>& reads, const seqan::String<int>& matches, unsigned len) noexcept
{
    int limit = reads.size();
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static))
        for (int i = 0; i < limit; ++i)
            if (matches[i] != -1)					//only erases barcode from sequence if it could be matched
            {
                erase(reads[i].seq, 0, len);
            }
}

//Overload for deleting the barcodes in any case.
template<typename TRead>
void clipBarcodes(std::vector<TRead>& reads, int len) noexcept
{
    int limit = reads.size();
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static))
        for (int i = 0; i < limit; ++i)
        {
            erase(reads[i].seq, 0, len);
        }
}

struct ApproximateBarcodeMatching {};
struct ExactBarcodeMatching {};

template <typename TRead, typename TMatches>
void group(std::vector<TRead>& reads, const TMatches& matches, 
    const ExactBarcodeMatching&, bool exclude) noexcept
{
	for (unsigned i = 0; i < length(matches); ++i)
    {                                                       //adds index of sequence to respective group.
		if ((!exclude) || (matches[i] != -1))                //Check if unidentified seqs have to be excluded
        {                                                   //offset by 1 is necessary since group 0 is...
            reads[i].demuxResult = matches[i] + 1;               //...reserved for unidentified sequences)
        }
    }                                                  
}

template <typename TRead, typename TMatches, typename TBarcodes>
void group(std::vector<TRead>& reads, const TMatches& matches,
    const TBarcodes& barcodes, const ExactBarcodeMatching& dummy, bool exclude) noexcept
{
    (void)barcodes;
    group(reads, matches, dummy, exclude);
}

//Overload if approximate search has been used.
template <typename TRead, typename TMatches, typename TBarcodes>
void group(std::vector<TRead>& reads, const TMatches& matches,
    const TBarcodes& barcodes, const ApproximateBarcodeMatching&, bool exclude) noexcept
{
    float dividend = float(length(barcodes[0])*5.0);		//value by which the index will be corrected.
    for (unsigned i = 0; i < length(matches); ++i)			//adds index of sequence to respective group.
    {
        if ((!exclude) || (matches[i] != -1))                //Check if unidentified seqs have to be excluded
        {
            reads[i].demuxResult = int(floor(float(matches[i]) / dividend)) + 1;
        }
    }
}

template<template <typename> class TRead, typename TSeq, typename TBarcodes, typename TFinder, typename TApprox>
void doAll(std::vector<TRead<TSeq>>& reads, TBarcodes& barcodes, TFinder& esaFinder,
    bool hardClip, DemultiplexStats& stats, const TApprox approximate, bool exclude)
{
    std::vector<TSeq> prefices(length(reads));
    getPrefix(prefices, reads, length(barcodes[0]));
    std::vector<int> matches(length(prefices));
    findAllExactIndex(matches, prefices, esaFinder, stats);
    if (hardClip)		//clip barcodes according to selected method
    {
        clipBarcodes(reads, length(barcodes[0]));
    }
    else
    {
        clipBarcodes(reads, matches, length(barcodes[0]));
    }
    group(reads, matches, barcodes, approximate, exclude);
}

//Overload for single-end data.
template<typename TSeqs, typename TIds>
void buildSets(TSeqs& seqs, TIds& ids, const seqan::StringSet<seqan::String<int> >& groups, seqan::String<TSeqs>& gSeqs, seqan::String<TIds>& gIds)
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
