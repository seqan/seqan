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
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================


#ifndef DEMULTIPLEX_H
#define DEMULTIPLEX_H

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
	std::vector<std::string> barcodes;
	std::vector<std::string> barcodeIds;
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
// Objects
// ============================================================================

struct BarcodeMatcher
{
    template <typename TContainer>
    BarcodeMatcher(const TContainer& patterns)
        : _patterns(patterns)
    {
        if (patterns.empty())
            return;
        _barcodeLength = _patterns[0].size();;
        for (const auto& pattern : patterns)
        {
            (void)pattern;
            assert(pattern.size() == _barcodeLength);
        }
    }
    template <template <typename> class TRead, typename TSeq>
    int getMatchIndex(const TRead<TSeq>& read) const noexcept
    {
        unsigned int index = 0;
        const std::string prefix = getPrefix(read, getBarcodeLength());
        for (const auto& pattern : _patterns)
        {
            if (pattern == prefix)   // likely
                return index;
            ++index;
        }
        return -1;
    }
    inline unsigned int getBarcodeLength() const noexcept
    {
        return _barcodeLength;
    }
private:
    const std::vector<std::string> _patterns;
    unsigned int _barcodeLength;
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

template <typename TBarcodes, typename TBarcode>
void buildVariations(TBarcodes& variations, const TBarcode& barcode)
{
	int limit = (length(barcode))*5;		    //possible number of variations with one error (A,T,G,C,N)
	resize(variations, limit);				    //resizes according to calculated number of variations
	for (auto& variation : variations)
        variation = barcode;	        //fills resultset with original barcode
	limit = limit/5;
	for (int i = 0; i < limit; ++i)
	{										    
		variations[5*i][i] = 'A';	    	
		variations[5*i+1][i] = 'C';
		variations[5*i+2][i] = 'G';
		variations[5*i+3][i] = 'T';
		variations[5*i+4][i] = 'N';
	}
}

template <typename TBarcodes>
void buildAllVariations(TBarcodes& barcodes)
{
    TBarcodes newbarcodes;
	for(const auto& barcode : barcodes)
	{
        TBarcodes tempBarcodes;
        buildVariations(tempBarcodes, barcode);
        for (const auto& tempBarcode : tempBarcodes)
			seqan::appendValue(newbarcodes, tempBarcode);
	}
	barcodes = newbarcodes;
}

// only used for testing
template <template <typename> class TRead, typename TSeq, typename TBarcodeFinder>
void MatchBarcodes(std::vector<TRead<TSeq>>& reads, const TBarcodeFinder& finder) noexcept
{
    std::for_each(reads.begin(), reads.end(), [&finder](auto& read){
        read.demuxResult = finder.getMatchIndex(read);});
}

struct ClipHard {};
struct ClipSoft {};

template <typename TRead>
void clipBarcodes(std::vector<TRead>& reads, const unsigned len, const ClipHard&) noexcept
{
    for (auto& read : reads)
        erase(read.seq, 0, len);
}

//Overload for deleting only matched barcodes 
template<typename TRead>
void clipBarcodes(std::vector<TRead>& reads, const int len, const ClipSoft&) noexcept
{
    // unmatched reads are partitioned to the end
    std::for_each(reads.begin(), std::find_if(reads.begin(), reads.end(), [](const auto& read)->auto {return read.demuxResult == 0;}), [len](auto& read) {
        erase(read.seq, 0, len);
    });
}

struct ApproximateBarcodeMatching {};
struct ExactBarcodeMatching {};

template <typename TRead, typename TFinder, typename TStats>
void MatchBarcodes(std::vector<TRead>& reads, const TFinder& finder, TStats& stats, const ExactBarcodeMatching&)
{
    for (auto& read : reads)
    {
        read.demuxResult = finder.getMatchIndex(read) + 1;
        if (stats.matchedBarcodeReads.size() < static_cast<unsigned int>(read.demuxResult) + 1)
        {
            std::cout << "error: matchedBarcodeReads too small!" << std::endl;
            throw(std::runtime_error("error: matchedBarcodeReads too small!"));
        }
        ++stats.matchedBarcodeReads[read.demuxResult];
    }
}

//Overload if approximate search has been used.
template <typename TRead, typename TFinder, typename TStats>
void MatchBarcodes(std::vector<TRead>& reads, const TFinder& finder, TStats& stats, const ApproximateBarcodeMatching&)
{
    const float dividend = float(finder.getBarcodeLength()*5.0);		//value by which the index will be corrected.
    for (auto& read: reads)			             
    {
        read.demuxResult = finder.getMatchIndex(read) ;
        if (read.demuxResult != -1)
            read.demuxResult = int(floor(float(read.demuxResult) / dividend));
        ++read.demuxResult;
        if (stats.matchedBarcodeReads.size() < static_cast<unsigned int>(read.demuxResult) + 1)
        {
            std::cout << "error: matchedBarcodeReads too small!" << std::endl;
            throw(std::runtime_error("error: matchedBarcodeReads too small!"));
        }
        ++stats.matchedBarcodeReads[read.demuxResult];
    }
}

template<template <typename> class TRead, typename TSeq, typename TFinder, typename TStats, typename TApprox>
void demultiplex(std::vector<TRead<TSeq>>& reads, const TFinder& finder,
    const bool hardClip, TStats& stats, const TApprox& approximate, const bool exclude)
{
    MatchBarcodes(reads, finder, stats, approximate);
    if (exclude)
        reads.erase(std::remove_if(reads.begin(), reads.end(), [](const auto& read)->auto {return read.demuxResult == 0;}), reads.end());
    else
        std::partition(reads.begin(), reads.end(), [](const auto& read)->auto {return read.demuxResult != 0;});

    if (std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value)   // clipping is not done for multiplex barcodes, only for inline barcodes
    {
        if (hardClip)
            clipBarcodes(reads, finder.getBarcodeLength(), ClipHard());
        else
            clipBarcodes(reads, finder.getBarcodeLength(), ClipSoft());
    }
}
#endif 
