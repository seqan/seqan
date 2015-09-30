// ==========================================================================
//                               readTrimming.h
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
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// ==========================================================================

#ifndef READTRIMMING_H
#define READTRIMMING_H

#include "helper_functions.h"
#include "general_stats.h"

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Tags for choosing the appropriate trimming method.
//brief The tagging structure for the Tail trimming algorithm.
struct Tail {};
//brief The tagging structure for the BWA trimming algorithm.
struct BWA {};
//brief The tagging structure for the window trimming algorithm.
struct Mean {
	unsigned window;	// The window size used for quality calculations.
	Mean(unsigned w) : window(w) {}
};

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

inline unsigned getQuality(const seqan::String<seqan::Dna5Q>& seq, unsigned i)
{
	return seqan::getQualityValue(seq[i]);
}

// Trimming methods
// ----------------------------------------------------------------------------
template <typename TSeq>
unsigned _trimRead(const TSeq& seq, unsigned const cutoff, Tail const &) noexcept
{
	for (int i = length(seq) - 1; i >= 0; --i)
    {
		if (getQuality(seq, i) >= cutoff)
        {
			return i + 1;
        }
    }
	return 0;
}

//Trimming mechanism using BWA. Trim to argmax_x sum_{i=x+1}^l {cutoff - q_i}
template <typename TSeq>
unsigned _trimRead(const TSeq& seq, unsigned const cutoff, BWA const &) noexcept
{
	int max_arg = length(seq) - 1, sum = 0, max = 0;
	for (int i = length(seq) - 1; i >= 0; --i)
	{
		sum += cutoff - getQuality(seq, i);
		if (sum < 0)
        {
			break;
        }
		if (sum > max)
		{
			max = sum;
			max_arg = i;
		}
	}
	return max_arg + 1;
}

// Trim by shifting a window over the sequence and cut where avg. qual. in window turns bad the first time.
template <typename TSeq>
unsigned _trimRead(const TSeq& seq, unsigned const _cutoff, Mean const & spec) noexcept
{
	unsigned window = spec.window;
	unsigned avg = 0, i = 0;
	// Work with absolute cutoff in window to avoid divisions.
	unsigned cutoff = _cutoff*window;
	// Calculate average quality of initial window.
	for (i = 0; i < window; ++i)
    {
		avg += getQuality(seq, i);
    }
	// Shift window over read and keep mean quality, update in constant time.
	for (i = 0; i < length(seq) && avg >= cutoff; ++i)
	{   // Take care only not to go over the end of the sequence. Shorten window near the end.
		avg -= getQuality(seq, i);
		avg += i + window < length(seq) ? getQuality(seq, i + window) : 0;
	}
	return i;   // i now holds the start of the first window that turned bad.
}

template <typename TSeq, typename TSpec>
unsigned trimRead(TSeq& seq, unsigned const cutoff, TSpec const & spec) noexcept
{
	unsigned ret, cut_pos;
	cut_pos = _trimRead(seqan::Dna5QString(seq), cutoff, spec);
	ret = length(seq) - cut_pos;
	erase(seq, cut_pos, length(seq));
    return ret;
}

template<bool tag>
struct TagTrimming
{
    static const bool value = tag;
};

template <typename TRead, typename TSpec, typename TTagTrimming>
unsigned _trimReads(std::vector<TRead>& reads, unsigned const cutoff, const TSpec& spec, TTagTrimming) noexcept(!TTagTrimming::value)
{
    int trimmedReads = 0;
    std::transform(reads.begin(),reads.end(),reads.begin(),[&trimmedReads, cutoff, &spec](auto& read)->auto
    {
        if (trimRead(read.seq, cutoff, spec))
        {
            ++trimmedReads;
            if (TTagTrimming::value)
                append(read.id, "[Trimmed]");
        }
        return read;
    });
    return trimmedReads;
}

template <typename TRead, typename TSpec>
unsigned trimBatch(std::vector<TRead>& reads, unsigned const cutoff, TSpec const& spec, bool tagOpt)
{
    unsigned trimmedReads;
    if(tagOpt)
        trimmedReads = _trimReads(reads, cutoff, spec, TagTrimming<true>());
    else
        trimmedReads = _trimReads(reads, cutoff, spec, TagTrimming<false>());
    return trimmedReads;
}
#endif
