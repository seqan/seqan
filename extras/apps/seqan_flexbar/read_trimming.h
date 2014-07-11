// ==========================================================================
//                               readTrimming.h
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Benjamin Strauch <b.strauch@fu-berlin.de>
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// ==========================================================================
// This file provides the quality control functionality of seqan-flexbar
// which is based in the implementation of the original flexbar program
// in [1].
// [1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, C.  FLEXBAR—Flexible
// Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
// Biology 2012, 1, 895-905.
// ==========================================================================


#ifndef SANDBOX_GROUP3_APPS_SEQDPT_READTRIMMING_H_
#define SANDBOX_GROUP3_APPS_SEQDPT_READTRIMMING_H_

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

/**
.Class.QualityTrimmingStats:
..summary:Struct holding information about certain quality trimming statistics.
.Memvar.QualityTrimmingStats#dropped_1
..class:Class.QualityTrimmingStats
..summary:Unsigned int holding the number of dropped forward reads
..type:nolink:unsigned
.Memvar.QualityTrimmingStats#dropped_2
..class:Class.QualityTrimmingStats
..summary:Unsigned int holding the number of dropped backward reads
..type:nolink:unsigned

.Memfunc.QualityTrimmingStats#clear:
..class:Class.QualityTrimmingStats
..summary:Resets the object.
..signature:clear(void)
..remarks:Both member variabels are set to 0.
*/
struct QualityTrimmingStats
{
	unsigned dropped_1, dropped_2;
	QualityTrimmingStats() : dropped_1(0), dropped_2(0) {};
	void clear()
	{
		dropped_1 = 0;
		dropped_2 = 0;
	}
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

/**
.Function.getQuality:
..summary:Determines the quality at position i of the Dna5QString.
..signature:loadSeqs(seq, i)
..param.seq:The Dna5QString containing bases and qualities.
...type:Class.String
..param.i:The index of the base whose quality shall be returned.
...type:nolink:unsigned
..returns:Phred quality of the base at position i.

 */
inline unsigned getQuality(const seqan::String<seqan::Dna5Q>& seq, unsigned i)
{
	return seqan::getQualityValue(seq[i]);
}

// Trimming methods
// ----------------------------------------------------------------------------
/**
.Function._trimRead:
..summary:Funtion for finding the trimming position for a given read.
..signature:_trimRead(seq, cutoff, spec)
..param.seq:The sequence that shall be trimmed..
...type:Class.String
..param.cutoff:The minimum quality required..
...type:nolink:unsigned
..param.spec:Specialisation Tag
...type:nolink:Taik
...type:nolink:BWA
...type:nolink:Mean
..returns:The trimming position, i.e. the first base to be removed.
...type:nolink:unsigned
..remarks:This function has three specialisations: The "Tail" method simply cuts off as many low quality bases from
 the end as possible before finding a good base with good quality. The "BWA" method uses the formular 
 argmax_x sum_{i=x+1}^l {cutoff - q_i} to trim to argmax. The "Mean" method uses a sliding window to calculate
 mean qualities.
 */
template <typename TSeq>
unsigned _trimRead(const TSeq& seq, unsigned const cutoff, Tail const &)
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
unsigned _trimRead(const TSeq& seq, unsigned const cutoff, BWA const &)
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
unsigned _trimRead(const TSeq& seq, unsigned const _cutoff, Mean const & spec)
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

/**
.Function.trimRead:
..summary: Interface that trims a sequence.
..signature:trimRead(seq, qual, cutoff, spec)
..param.seq:The sequence which shall be trimmed.
...type:Class.String
..param.qual:The quality scores of the sequence.
...type:nolink:String
..param.cutoff:The minimum quality required from a base.
...type:nolink:unsigned
..param.spec:The trimming algorithm used for trimming.
...type:nolink:Tail
...type:nolink:BWA
...type:nolink:Mean
..returns:The number of bases trimmed from the sequence.
...type:nolink:unsigned
 */
template <typename TSeq, typename TSpec>
unsigned trimRead(TSeq& seq, unsigned const cutoff, TSpec const & spec)
{
	unsigned ret, cut_pos;
	cut_pos = _trimRead(seq, cutoff, spec);
	ret = length(seq) - cut_pos;
	erase(seq, cut_pos, length(seq));
    return ret;
}
/**
.Function._trimReads:
..summary:Trims a set of reads.
..signature:_trimReads(seqSet, cutoff, spec)
..param.seqSet: A collection of sequences that shall be trimmed.
...type:Class.StringSet
..param.cutoff:The minimum quality required from a base.
...type:nolink:unsigned int
..param.spec:The trimming algorithm used for trimming.
...type:nolink:Tail
...type:nolink:BWA
...type:nolink:Mean
..returns:The number of reads which had bases removed.
...type:nolink:unsigned
 */
template <typename TSet, typename TIdSet, typename TSpec>
unsigned _trimReads(TSet & seqSet, TIdSet& idSet, unsigned const cutoff, TSpec const & spec, bool tagOpt)
{
	typedef typename seqan::Value<TSet>::Type TSeq;
	int trimmedReads = 0;
	int len = length(seqSet);
	SEQAN_OMP_PRAGMA(parallel for schedule(static) reduction(+:trimmedReads))
	for (int i=0; i < len; ++i)
	{
		TSeq& read = value(seqSet, i);
		unsigned trimmed = trimRead(read, cutoff, spec);
		if (trimmed > 0)
        {
            ++trimmedReads;
            if (tagOpt)
            {
                append(idSet[i], "[Trimmed]");
            }
        }
	}
	return trimmedReads;
}

/**
.Function.dropReads:
..summary:Drops reads which are too short. This is done in a way
 the pair structure is conserved.
..signature:dropReads(idSet1, seqSet1, min_length, stats)
..signature:dropReads(idSet1, seqSet1, idSet2, seqSet2, min_length, stats)
..param.idSet1:Set of FastA-IDs for the set of forward reads.
...type:Class.StringSet
..param.seqSet1:Set containing the forward reads.
...type:Class.StringSet
..param.idSet2:Set of FastA-IDs for the set of backward reads.
...type:Class.StringSet
..param.seqSet2:Set containing the backward reads.
...type:Class.StringSet
..param.min_length:The minimum length required after trimming. Shorter sequences will be deleted.
...type:nolink:unsigned
..param.stats:QualityTrimmingStats object storing the number of dropped reads.
...type:Class.QualityTrimmingStats
..returns:void
..remarks:In paired end case: If one read is marked for removal, while its sibling is not, the read will
 be replaced by a single 'N'. If both reads are marked for removal, they are removed.
 This way the relation between paired reads stays intact.
 */
template <typename TId, typename TSeq>
unsigned dropReads(seqan::StringSet<TId> & idSet, seqan::StringSet<TSeq> & seqSet,
		unsigned const min_length, QualityTrimmingStats& stats)
{
	int len = length(seqSet);
    seqan::StringSet<bool> rem;
    resize(rem, len);
    SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static))
    for (int i = 0; i < len; ++i)
    {
        if (length(seqSet[i]) < min_length)
        {
            rem[i] = true;
        }
        else
        {
            rem[i] = false;
        }
    }

    unsigned ex = 0;
	for (int i = len - 1; i >= 0; --i)
	{
		if (rem[i])
        {
            seqan::swap(seqSet[i], seqSet[len - ex - 1]);
            seqan::swap(idSet[i], idSet[len - ex - 1]);
            ++ex;
        }
	}
    if (ex != 0)
    {
        seqan::resize(seqSet, len - ex);
        seqan::resize(idSet, len - ex);
        stats.dropped_1 += ex;
    }
	return 0;
}

//overload for paired end data
template <typename TId, typename TSeq>
unsigned dropReads(seqan::StringSet<TId> & idSet1, seqan::StringSet<TSeq> & seqSet1,
		  seqan::StringSet<TId> & idSet2, seqan::StringSet<TSeq> & seqSet2, unsigned const min_length, QualityTrimmingStats& stats)
{
	// Iterate over the set and drop filter out those reads that are
	// too short. If only one read of the pair is too short, mark it
	// with a single N. If both reads are too short, remove them.
	// Possible feature for the future generation: Write out orphaned reads to extra file.
    int len = length(seqSet1);
    seqan::StringSet<bool> rem;
    resize(rem, len);
    unsigned dropped1 = 0;
    unsigned dropped2 = 0;
    SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static) reduction(+:dropped1, dropped2))
    for (int i = 0; i < len; ++i)
    {
        bool drop1 = length(seqSet1[i]) < min_length;
		bool drop2 = length(seqSet2[i]) < min_length;
        if (drop1 && drop2)
        {
            rem[i] = true;
            ++dropped1;
            ++dropped2;
        }
        else if (drop1)
        {
            rem[i] = false;
            ++dropped1;
			seqan::moveValue(seqSet1, i, TSeq("N"));
        }
        else if (drop2)
        {
            rem[i] = false;
            ++dropped2;
			seqan::moveValue(seqSet2, i, TSeq("N"));
        } 
        else
        {
            rem[i] = false;
        }
    }
    stats.dropped_1 += dropped1;
    stats.dropped_2 += dropped2;
    unsigned ex = 0;
    for (int i = len - 1; i >= 0; --i)
    {
        if (rem[i])
        {
            seqan::swap(seqSet1[i], seqSet1[len - ex - 1]);
            seqan::swap(idSet1[i], idSet1[len - ex - 1]);
            seqan::swap(seqSet2[i], seqSet2[len - ex - 1]);
            seqan::swap(idSet2[i], idSet2[len - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqSet1, len - ex);
        resize(idSet1, len - ex);
        resize(seqSet2, len - ex);
        resize(idSet2, len - ex);
    }
	return 0;
}

/**
.Function.trimBatch:
..summary:Trims bad quality bases from a set of sequences.
..signature:trimBatch(seqSet, cutoff, spec)
..param.seqSet:StringSet containing the reads.
...type:Class.StringSet
..param.cutoff:The minimum quality required of a base.
...type:nolink:unsigned int
..param.spec:The trimming algorithm used for trimming.
...type:nolink:Tail
...type:nolink:BWA
...type:nolink:Mean
..returns:The number of sequences which had bases removed from.
...type:nolink:unsigned
 */
template <typename TSeq, typename TId, typename TSpec>
unsigned trimBatch(seqan::StringSet<TSeq>& seqSet, seqan::StringSet<TId>& idSet, unsigned const cutoff,
    TSpec const& spec, bool tagOpt)
{
	unsigned trimmedReads = _trimReads(seqSet, idSet, cutoff, spec, tagOpt);
	return trimmedReads;
}

/**
.Function.trimPairBatch:
..summary:Trims bad quality bases from two sets of sequences.
..signature:trimPairBatch(seqSet1, seqSet2, cutoff, spec)
..param.seqSet1:StringSet containing the forward reads of the paired sequences.
...type:Class.StringSet
..param.seqSet2:StringSet containing the backward reads of the paired sequences.
...type:Class.StringSet
..param.cutoff:The minimum quality required of a base.
...type:nolink:unsigned int
..param.spec:The trimming algorithm used for trimming.
...type:nolink:Tail
...type:nolink:BWA
...type:nolink:Mean
..returns:A pair of unsigned ints containing the number of reads which had bases removed from.
...type:nolink:pair
 */
template <typename TSeq, typename TId, typename TSpec>
seqan::Pair<unsigned, unsigned> trimPairBatch(seqan::StringSet<TSeq>& seqSet1, seqan::StringSet<TId>& idSet1,
    seqan::StringSet<TSeq> & seqSet2, seqan::StringSet<TId>& idSet2, unsigned const cutoff,
    TSpec const & spec, bool tagOpt)
{
	unsigned trimmedReads1 = _trimReads(seqSet1, idSet1, cutoff, spec, tagOpt);
	unsigned trimmedReads2 = _trimReads(seqSet2, idSet2, cutoff, spec, tagOpt);
	return seqan::Pair<unsigned, unsigned>(trimmedReads1, trimmedReads2);
}
#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_READTRIMMING_H_
