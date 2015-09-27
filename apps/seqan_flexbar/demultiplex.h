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
#include "general_stats.h"


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct DemultiplexingParams
{
	std::string barcodeFile;
	seqan::StringSet<seqan::String<seqan::Dna5> > barcodes;
	seqan::StringSet<seqan::CharString> barcodeIds;
	std::string multiplexFile;
	bool approximate;
	bool hardClip;
	bool run;
	bool runx;
	bool exclude;

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
    for (const auto& barcode : barcodes)
    {
        if (len != length(barcode))
        {
            std::cerr << "ERROR: Barcodes differ in length. All barcodes must be of equal length.\n";
            return false;
        }
    }
    auto it = std::remove_if(reads.begin(), reads.end(), [len](auto& read) {return length(read.seq) <= len;});
    stats.removedShort += std::distance(it, reads.end());
    reads.erase(it, reads.end());
    return true;
}

// always use the forward read for barcode detection
template <template <typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>,Read<TSeq>>::value || std::is_same<TRead<TSeq>,ReadPairedEnd<TSeq>>::value>>
void getPrefix(std::vector<TSeq>& prefices, std::vector<TRead<TSeq>>& reads, unsigned len)
{
    assert(prefices.size() == reads.size());

    std::transform(reads.begin(), reads.end(), prefices.begin(), [len](const auto& read) {
        return prefix(read.seq, len);});
}

template <template <typename> class TRead, typename TSeq, typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value>>
void getPrefix(std::vector<TSeq>& prefices, std::vector<TRead<TSeq>>& reads, unsigned len, bool = false)
{
    (void)len;
    std::transform(reads.begin(), reads.end(), prefices.begin(), [](const auto& read) {
        return read.demultiplex;
    });
}

template <typename TBarcode>
void buildVariations(seqan::StringSet<seqan::Dna5String>& variations, const TBarcode& barcode)
{
	int limit = (length(barcode))*5;		    //possible number of variations with one error (A,T,G,C,N)
	resize(variations, limit);				    //resizes according to calculated number of variations
	for (int i = 0; i < limit; ++i)
	{
        assign(variations[i], barcode);	        //fills resultset with original barcode
	}
	limit = limit/5;
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
inline int findExactIndex(const TPrefix& prefix, TFinder& finder) noexcept
{
	clear(finder);								//resets finder
	if (find(finder, prefix))
    {
		return getSeqNo(position(finder));		//returns index of barcode. ONLY THE FIRST HIT!
    }
	else return -1;								//return -1 if no hit occured
}

template <typename TMatches, typename TPrefices, typename TFinder>
void findAllExactIndex(TMatches& matches, const TPrefices& prefices, TFinder finder) noexcept
{
    assert(length(matches) == length(prefices));

    std::transform(seqan::begin(prefices), seqan::end(prefices), seqan::begin(matches), [&finder](const auto& prefix)->auto {
        return findExactIndex(prefix, finder);});
}

template <typename TRead>
void clipBarcodes(std::vector<TRead>& reads, const seqan::String<int>& matches, const unsigned len) noexcept
{
    int limit = reads.size();
    for (int i = 0; i < limit; ++i)
        if (matches[i] != -1)					//only erases barcode from sequence if it could be matched
        {
            erase(reads[i].seq, 0, len);
        }
}

//Overload for deleting the barcodes in any case.
template<typename TRead>
void clipBarcodes(std::vector<TRead>& reads, const int len) noexcept
{
    std::transform(reads.begin(), reads.end(), reads.begin(), [len](auto& read)->auto {
        erase(read.seq, 0, len);
        return read;});
}

struct ApproximateBarcodeMatching {};
struct ExactBarcodeMatching {};

template <typename TRead, typename TMatches, typename TStats>
void group(std::vector<TRead>& reads, const TMatches& matches, TStats& stats,
    const ExactBarcodeMatching&, const bool exclude) noexcept
{
    unsigned i = 0;
	for (const auto matchResult : matches)
    {                                                       //adds index of sequence to respective group.
		if ((!exclude) || (matchResult != -1))                //Check if unidentified seqs have to be excluded
        {                                                   //offset by 1 is necessary since group 0 is...
            reads[i++].demuxResult = matchResult + 1;               //...reserved for unidentified sequences)
            ++stats.matchedBarcodeReads[matchResult + 1];
        }
    }                                                  
}

template <typename TRead, typename TMatches, typename TBarcodes, typename TStats>
void group(std::vector<TRead>& reads, const TMatches& matches,
    const TBarcodes& barcodes, TStats& stats, const ExactBarcodeMatching& dummy, bool exclude) noexcept
{
    (void)barcodes;
    group(reads, matches, stats, dummy, exclude);
}

//Overload if approximate search has been used.
template <typename TRead, typename TMatches, typename TBarcodes, typename TStats>
void group(std::vector<TRead>& reads, const TMatches& matches,
    const TBarcodes& barcodes, TStats& stats, const ApproximateBarcodeMatching&, bool exclude) noexcept
{
    unsigned i = 0;
    float dividend = float(length(barcodes[0])*5.0);		//value by which the index will be corrected.
    for (const int matchResult : matches)			        //adds index of sequence to respective group.
    {
        if ((!exclude) || (matchResult != -1))                //Check if unidentified reads have to be excluded
        {
            reads[i++].demuxResult = int(floor(float(matchResult) / dividend)) + 1;
            ++stats.matchedBarcodeReads[static_cast<int>(floor(float(matchResult) / dividend)) + 1];
        }
    }
}

template<template <typename> class TRead, typename TSeq, typename TBarcodes, typename TFinder, typename TStats, typename TApprox>
void doAll(std::vector<TRead<TSeq>>& reads, const TBarcodes& barcodes, const TFinder& esaFinder,
    const bool hardClip, TStats& stats, const TApprox& approximate, bool exclude)
{
    std::vector<TSeq> prefices(length(reads));
    getPrefix(prefices, reads, length(barcodes[0]));
    std::vector<int> matches(length(prefices));
    findAllExactIndex(matches, prefices, esaFinder);
    if (std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value)   // clipping is not done for multiplex barcodes, only for inline barcodes
    {
        if (hardClip)		//clip barcodes according to selected method
        {
            clipBarcodes(reads, length(barcodes[0]));
        }
        else
        {
            clipBarcodes(reads, matches, length(barcodes[0]));
        }
    }
    group(reads, matches, barcodes, stats, approximate, exclude);
}
#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_DEMULTIPLEX_H_
