// ==========================================================================
//                             adapterTrimming.h
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
// Author: Benjamin Strauch <b.strauch@fu-berlin.de>
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================
// This file provides the adapter trimming functionality of seqan-flexbar
// which is based in the implementation of the original flexbar program in 
// [1].
// [1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, C.  FLEXBAR—Flexible
// Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
// Biology 2012, 1, 895-905.
// ==========================================================================

#ifndef SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_
#define SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_

#include <seqan/align.h>
#include <seqan/find.h>
#include "helper_functions"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Custom scoring matrix for ignoring alignments with N.
namespace seqan{
	//Struct used to define a new custom scoring matrix.
	struct AdapterScoringMatrix {};
     //Struct containing data for the seqan::AdapterScoringMatrix custom scoring matrix.
	 //Matches score 1, mismatches -1 and matches against N with 0.
	template <>
	struct ScoringMatrixData_<int, Dna5, AdapterScoringMatrix> {
		enum {
			VALUE_SIZE = ValueSize<Dna5>::VALUE,
			TAB_SIZE = VALUE_SIZE * VALUE_SIZE
		};
		static inline int const * getData() {
			static int const _data[TAB_SIZE] = {
			   1, -1, -1, -1, 0,
			  -1,  1, -1, -1, 0,
			  -1, -1,  1, -1, 0,
			  -1, -1, -1,  1, 0,
			   0,  0,  0,  0, 0
			};
			return _data;
		}
	};
}

struct AdapterItem;
typedef seqan::Dna5QString TAdapter;
typedef std::vector< AdapterItem > TAdapterSet;
struct AdapterItem
{
    AdapterItem() : anchored(false) {};
    AdapterItem(const TAdapter &adapter) : seq(adapter), anchored(false) {};

    TAdapter seq;
    // the anchored matching mode should have the same functionality as cutadapt's anchored mode
    // see http://cutadapt.readthedocs.org/en/latest/guide.html
    // In anchored mode, the adapter has to start (3" adapter) or has to end (5" adapter) with the sequence
    // The anchored mode is rarely used, at least for 3" adapters
    // Todo: implement rooted mode
    bool anchored;
};

// Define scoring function type.
typedef Score<int, ScoreMatrix<Dna5, AdapterScoringMatrix> > TScore;

//Tagging struct representing the the match algorithm working
//with values supplied by the user. Saves those  values as members.
struct AdapterMatchSettings
{
    AdapterMatchSettings(int m, int e, double er) : min_length(m), errors(e), errorRate(er), erMode(false), modeAuto(false), singleMatch(false)
    {
        erMode = ((e == 0) && (er != 0));
    }
    AdapterMatchSettings() : min_length(0), errors(0), errorRate(0), erMode(false), modeAuto(true), singleMatch(false) {};
   
    int min_length; //The minimum length of the overlap.
	int errors;     //The maximum number of errors we allow.
	double errorRate;  //The maximum number of errors allowed per overlap
	bool erMode;
    bool modeAuto;
    bool singleMatch;
};

struct AdapterTrimmingStats
{
	AdapterTrimmingStats() : a1count(0), a2count(0), overlapSum(0),
			minOverlap(std::numeric_limits<unsigned>::max()), maxOverlap(0) {};

    AdapterTrimmingStats& operator+= (AdapterTrimmingStats const& rhs)
    {
        a1count += rhs.a1count;
        a2count += rhs.a2count;
        overlapSum += rhs.overlapSum;
        minOverlap = minOverlap < rhs.minOverlap ? minOverlap : rhs.minOverlap;
        maxOverlap = maxOverlap < rhs.maxOverlap ? rhs.maxOverlap : maxOverlap;
        return *this;
    }
	void clear()
	{
		a1count = 0;
		a2count = 0;
		overlapSum = 0;
		minOverlap = std::numeric_limits<unsigned>::max();
		maxOverlap = 0;
	}

    unsigned a1count, a2count;
    unsigned overlapSum;
    unsigned minOverlap, maxOverlap;
};

// ============================================================================
// Metafunctions
// ============================================================================


//brief A metafunction which constructs a ModifiedString type for an Alphabet.
template <class TValue>
struct STRING_REVERSE_COMPLEMENT
{
	typedef ModifiedString<
			ModifiedString<	String<TValue>, ModView<FunctorComplement<TValue> > >,
			ModReverse>	Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TSeq1, typename TSeq2, bool TTop, bool TLeft, bool TRight, bool TBottom>
void alignPair(seqan::Pair<unsigned, seqan::Align<TSeq1> >& ret, TSeq1& seq1, TSeq2& seq2,
		const seqan::AlignConfig<TTop, TLeft, TRight, TBottom>& config, bool band = false)
{
    seqan::resize(rows(ret.i2), 2);
    seqan::assignSource(row(ret.i2, 0), seq1);
    seqan::assignSource(row(ret.i2, 1), seq2);
	// Overlap alignment. We again don't allow gaps, since they are unlikely in Illumina data.
	// We mostly use AlignConfigs true,false,true,false and true,true,true,true here.
	TScore adapterScore(-100);
	// If we can do a banded alignment, restrict ourselves to the upper half of the matrix.
    if (band)
    {
        ret.i1 = globalAlignment(ret.i2, adapterScore, config, -2, length(seq1));
    }
    else
    {
        ret.i1 = globalAlignment(ret.i2, adapterScore, config);
    }
}

template <typename TRow>
unsigned countTotalGaps(TRow& row)
{
	return length(row) - length(source(row));
}

template <typename TAlign>
unsigned getOverlap(TAlign& align)
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	TRow &row1 = seqan::row(align,0);
	TRow &row2 = seqan::row(align,1);
	return length(source(row1)) - countTotalGaps(row2);
}

template <typename TAlign>
unsigned getInsertSize(TAlign& align)
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	TRow &row1 = seqan::row(align,0);
	TRow &row2 = seqan::row(align,1);
	unsigned seq1_length = length(source(row1));
	unsigned seq2_length = length(source(row2));
	// Calculate overlap and overhangs.
	unsigned overlap = seq1_length - countTotalGaps(row2);
	unsigned seq2l = seqan::countGaps(seqan::begin(row1)); // Overhang of sequence 2 = Gaps at start of row 1.
	unsigned seq1r = countTotalGaps(row2) - seqan::countGaps(seqan::begin(row2)); // Overhang of sequence 1 = Gaps at end of row 2.
	// Insert size: Add sequence lengths, subtract common region (overlap)
	// and subtract overhangs left and right of insert.
	return seq1_length + seq2_length - overlap - (seq2l + seq1r);
}

template <typename TSeq>
unsigned stripPair(TSeq& seq1, TSeq& seq2)
{
	// When aligning the two sequences, the complementary sequence is reversed and
	// complemented, so we have an overlap alignment with complementary bases being the same.
	typedef typename Value<TSeq>::Type TAlphabet;
	typedef typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type TReverseComplement;
	TReverseComplement mod(seq2);
	typedef seqan::Align<TSeq> TAlign;
	seqan::Pair<unsigned, TAlign> ret;
    alignPair(ret, seq1, mod, seqan::AlignConfig<true, true, true, true>());
	unsigned score = ret.i1;
	TAlign align = ret.i2;
	// Use the overlap of the two sequences to determine the end position.
	unsigned overlap = getOverlap(align);
	unsigned mismatches = (overlap-score) / 2;
	// We require a certain correct overlap to exclude spurious hits.
	// (Especially reverse 3'->5' alignments not caused by adapters.)
	if (overlap <= 5 || mismatches > overlap * 0.15)
    {
		return 0;
    }
	// Get actual size of the insert (possible to determine from overlap etc.).
	unsigned insert = getInsertSize(align);
	// Now cut both sequences to insert size (no cuts happen if they are smaller)
	if (length(seq1) > insert)
    {
		seqan::erase(seq1, insert, length(seq1));
    }
	if (length(seq2) > insert)
    {
		seqan::erase(seq2, insert, length(seq2));
    }
	return insert;
}
 //Overload using adapter information.
 //Currently not in use, because not as good as the version without adapters.
template <typename TSeq>
unsigned stripPair(TSeq& seq1, TSeq const& adapter1, TSeq& seq2, TSeq const& adapter2)
{
	typedef typename Value<TSeq>::Type TAlphabet;
	typedef typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type TReverseComplement;
	// Add adapters to the sequences.
	seqan::insert(seq1, 0, adapter2);
	seqan::insert(seq2, 0, adapter1);
	TReverseComplement mod(seq2);
	typedef seqan::Align<TSeq> TAlign;
	seqan::Pair<unsigned, TAlign> ret;
    alignPair(ret, seq1, mod, seqan::AlignConfig<true, false, true, false>(), true);
	TAlign align = ret.i2;
	unsigned score = ret.i1;
	unsigned overlap = getOverlap(align);
	unsigned mismatches = (overlap-score) / 2;
	// We can still get the proper insert size from the overlap, since we only
	// get properly overlapping alignments in this alignment mode. We just have
	// to subtract the adapter additions. One exception is the case where they
	// don't overlap at all, but in this case getInsertSize returns 0.
	unsigned insert = getInsertSize(align);
	if (insert > 0)
    {
		insert = insert - length(adapter1) - length(adapter2);
    }
	// Remove adapter sequences again.
	seqan::erase(seq1, 0, length(adapter2));
	seqan::erase(seq2, 0, length(adapter1));
	// Check the quality of the overlap.
	if (mismatches >= overlap * 0.15)
    {
		return 0;
    }
	// Cut down to insert size, if necessary.
	if (length(seq1) > insert)
    {
		seqan::erase(seq1, insert, length(seq1));
    }
	if (length(seq2) > insert)
    {
		seqan::erase(seq2, insert, length(seq2));
    }
	return insert;
}

template <typename TSeq, typename TId>
unsigned stripPairBatch(seqan::StringSet<TSeq>& set1, seqan::StringSet<TId>& idSet1,
    seqan::StringSet<TSeq>& set2, seqan::StringSet<TId>& idSet2, AdapterTrimmingStats& stats, bool tagOpt)
{
    int t_num = 1;
#ifdef _OPENMP
    t_num = omp_get_max_threads();
#endif
	// Create local counting variables to avoid concurrency problems.
	unsigned a1count = 0, a2count = 0, overlapSum = 0;
	seqan::String<unsigned> minOverlap;
	seqan::String<unsigned> maxOverlap;
	seqan::resize(minOverlap, t_num, std::numeric_limits<unsigned>::max());
	seqan::resize(maxOverlap, t_num, 0);
	int len = length(set1);
	SEQAN_OMP_PRAGMA(parallel for schedule(static) reduction(+:a1count, a2count, overlapSum))
	for (int i = 0; i < len; ++i)
	{
		// The reads processed in this iteration. The lengths will change after stripPair().
		TSeq& read1 = value(set1, i);
		TSeq& read2 = value(set2, i);
		unsigned len1 = length(read1);
		unsigned len2 = length(read2);
		stripPair(read1, read2);
		// Determine how much was cut off. (Length before and after.)
		unsigned over1 = len1 - length(read1);
		unsigned over2 = len2 - length(read2);
        // Apply tags
        if (over1 != 0 && tagOpt)
        {
            append(idSet1[i], "[AdapterRemoved]");
        }
        if (over2 != 0 && tagOpt)
        {
            append(idSet2[i], "[AdapterRemoved]");
        }
		// Count cuts.
		a1count += (over1 != 0);
		a2count += (over2 != 0);
		overlapSum += over1 + over2; // Count for each mate.
		// Thread saves local min/max seen in the batch it processed.
		int t_id = omp_get_thread_num();
		if (over1 > 0 && std::min(over1,over2) < minOverlap[t_id])
        {
            minOverlap[t_id] = std::min(over1,over2);
        }
		if (over2 > 0 && std::max(over1,over2) > maxOverlap[t_id])
        {
            maxOverlap[t_id] = std::max(over1,over2);
        }
	}
	// Update statistics using information gathered by the threads.
	stats.a1count += a1count;
	stats.a2count += a2count;
	stats.overlapSum += overlapSum;
	unsigned batch_min = *std::min_element(begin(minOverlap), end(minOverlap));
	unsigned batch_max = *std::max_element(begin(maxOverlap), end(maxOverlap));
	if (batch_min < stats.minOverlap) 
    {
        stats.minOverlap = batch_min;
    }
	if (batch_max > stats.maxOverlap)
    {
        stats.maxOverlap = batch_max;
    }
	return a1count + a2count;
}

template <typename TSeq, typename TAdapterItem>
void alignAdapter(seqan::Pair<unsigned, seqan::Align<TSeq> >& ret, TSeq& seq, TAdapterItem const& adapterItem)
{
	// Gaps are not allowed by setting the gap penalty high.
    // Changed the AlignConfig to behave like cutadapt's non anchored mode
    // examples using ADAPTER as adapter, required overlap 4, allowed error rate 0.2
    //
    // NNADAPTERMYSEQUENCE -> MYSEQUENCE
    // MYSEQUENCEADAPTERXX -> MYSEQUENCE
    // DAPTERMYSEQUENCE -> MYSEQUENCE
    // MYSEQUENCEDAPTER -> MYSEQUENC
    //
    // anchored mode would look like
    // NNADAPTERMYSEQUENCE -> NNADAPTERMYSEQUENCE
    // MYSEQUENCEADAPTERXX -> NNADAPTERMYSEQUENCE
    // DAPTERMYSEQUENCE -> NNADAPTERMYSEQUENCE
    // MYSEQUENCEDAPTER -> MYSEQUENC
    alignPair(ret, seq, adapterItem.seq, seqan::AlignConfig<true, true, true, true>(), true);
}

//Version for automatic matching options
inline bool isMatch(const int overlap, const int mismatches, const AdapterMatchSettings &adatperMatchSettings)
{
    if (overlap == 0)
        return false;
    if (adatperMatchSettings.modeAuto)
    {
        int errors = 0; // No errors for overlap up to 5 bases.
        if (overlap > 5)
        {
            errors = 1; // One error for overlap up to 10 bases.
        }
        if (overlap > 10)
        {
            errors = int(0.33 * overlap); //33% of overlap length otherwise.
        }
        return mismatches <= errors;
    }
    else 
    {
        if (adatperMatchSettings.erMode)
            return overlap >= adatperMatchSettings.min_length && (static_cast<double>(mismatches) / static_cast<double>(overlap)) <= adatperMatchSettings.errorRate;
        else
            return overlap >= adatperMatchSettings.min_length && mismatches <= adatperMatchSettings.errors;
    }
}

template <typename TSeq, typename TAdapterSet, typename TSpec>
unsigned stripAdapter(TSeq& seq, AdapterTrimmingStats& stats, TAdapterSet const& adapterSet, TSpec const& spec, const bool reverse)
{
	typedef typename seqan::Align<TSeq> TAlign;
    typedef typename seqan::Row<TAlign>::Type TRow;
	seqan::Pair<unsigned, TAlign> ret;

    unsigned removed{ 0 };
    for (AdapterItem const& adapterItem : adapterSet)
    {
        //std::cout << "seq: " << seq << std::endl;
        //std::cout << "adapter: " << adapterItem.seq << std::endl;
        alignAdapter(ret, seq, adapterItem);    // align crashes if seq is an empty string!
        const unsigned int overlap = getOverlap(ret.i2);
        const int score = ret.i1;
        const int mismatches = (overlap - score) / 2;
        if (isMatch(overlap, mismatches, spec))
        {
            //std::cout << "score: " << ret.i1 << " overlap: " << overlap << " mismatches: " << mismatches << std::endl;
            //std::cout << ret.i2 << std::endl;
            TRow row2 = row(ret.i2, 1);
            //std::cout << "adapter start position: " << toViewPosition(row2, 0) << std::endl;
            removed += length(seq) - toViewPosition(row2, 0);
            seqan::erase(seq, toViewPosition(row2, 0), length(seq));
            //std::cout << "stripped seq: " << seq << std::endl;
            stats.overlapSum += overlap;
            if (reverse)
                ++stats.a1count;
            else
                ++stats.a2count;

            stats.maxOverlap = stats.maxOverlap > overlap ? stats.maxOverlap : overlap;
            stats.minOverlap = stats.minOverlap < overlap ? stats.minOverlap : overlap;
            if(spec.singleMatch || empty(seq))
                return removed;
        }
    }
    return removed;
}

template <typename TSeq, typename TAdapterSet, typename TSpec>
unsigned stripAdapter(TSeq& seq, AdapterTrimmingStats& stats, TAdapterSet const& adapterSet, TSpec const& spec)
{
    return stripAdapter(seq, stats, adapterSet, spec, false);
}

template <typename TSeq, typename TId, typename TAdapterSet, typename TSpec>
unsigned stripAdapterBatch(seqan::StringSet<TSeq>& set, seqan::StringSet<TId>& idSet, TAdapterSet const& adapterSet, TSpec const& spec,
		AdapterTrimmingStats& stats, bool reverse = false, bool tagOpt = false)
{
	int t_num = omp_get_max_threads();
	// Create local counting variables to avoid concurrency problems.
	seqan::String<unsigned> minOverlap;
	seqan::String<unsigned> maxOverlap;
    std::vector<AdapterTrimmingStats> adapterTrimmingStatsVector;
	seqan::resize(minOverlap, t_num, std::numeric_limits<unsigned>::max());
	seqan::resize(maxOverlap, t_num, 0);
    seqan::resize(adapterTrimmingStatsVector, t_num);
	int len = length(set);
	SEQAN_OMP_PRAGMA(parallel for schedule(static))
    for (int i = 0; i < len; ++i)
    {
        if (empty(value(set, i)))
            continue;
        const int t_id = omp_get_thread_num();
        // Every thread has its own adapterTrimmingStatsVector
        const unsigned over = stripAdapter(value(set, i), adapterTrimmingStatsVector[t_id], adapterSet, spec, reverse);
        if (tagOpt && over != 0)
            insertAfterFirstToken(idSet[i], TId(":AdapterRemoved"));
	}
    std::for_each(adapterTrimmingStatsVector.begin(), adapterTrimmingStatsVector.end(), 
        [&stats](AdapterTrimmingStats const& _stats) {stats += _stats;});
	return stats.a1count + stats.a2count;
}

template <typename TSeq, typename TId, typename TSpec>
unsigned stripReverseAdapterBatch(seqan::StringSet<TSeq>& set, seqan::StringSet<TId>& IdSet, TSeq const& adapter, TSpec const& spec,
		AdapterTrimmingStats& stats, bool tagOpt)
{
	typedef typename Value<TSeq>::Type TAlphabet;
	typedef typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type TReverseComplement;
	TReverseComplement mod(adapter);
	return stripAdapterBatch(set, IdSet, mod, spec, stats, true, tagOpt);
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_
