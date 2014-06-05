// ==========================================================================
//                             adapterTrimming.h
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


#ifndef SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_
#define SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_

#include <seqan/align.h>
#include <seqan/find.h>

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

// Define scoring function type.
typedef Score<int, ScoreMatrix<Dna5, AdapterScoringMatrix> > TScore;


//brief A struct encapsulating information about the match algorithm.
struct Mode {};
//brief Tagging struct for the automatic match algorithm.
struct Auto : Mode {};

//Tagging struct representing the the match algorithm working
//with values supplied by the user. Saves those values as members.
struct User : Mode
{
	int min_length; //The minimum length of the overlap.
	int errors;     //The maximum number of errors we allow.
	User (int m, int e): min_length(m), errors(e) {}
};

/**
.Class.AdapterTrimmingStats:
..summary:Struct holding information about certain adapter trimming statistics.
..include:seqan/adapterTrimming.h
.Memvar.AdapterTrimmingStats#a1count
..class:Class.AdapterTrimmingStats
..summary:Unsigned int holding the number of adapters present in forward reads.
..type:nolink:unsigned
.Memvar.AdapterTrimmingStats#a2count
..class:Class.AdapterTrimmingStats
..summary:Unsigned int holding the number of adapters present in backwardreads.
..type:nolink:unsigned
.Memvar.AdapterTrimmingStats#overlapSum
..class:Class.AdapterTrimmingStats
..summary:Unsigned int holding the sum of the removed adapter bases. Later used to calculate mean adapter size.
..type:nolink:unsigned
.Memvar.AdapterTrimmingStats#minOverlap
..class:Class.AdapterTrimmingStats
..summary:Unsigned int holding the minimum overlap between adapter and reads.
..type:nolink:unsigned
.Memvar.AdapterTrimmingStats#maxOverlap
..class:Class.AdapterTrimmingStats
..summary:Unsigned int holding the maximum overlap between adapter and reads.
..type:nolink:unsigned
 */
struct AdapterTrimmingStats
{
	unsigned a1count, a2count;
	unsigned overlapSum;
	unsigned minOverlap, maxOverlap;
	AdapterTrimmingStats() : a1count(0), a2count(0), overlapSum(0),
			minOverlap(std::numeric_limits<unsigned>::max()), maxOverlap(0) {};
	void clear()
	{
		a1count = 0;
		a2count = 0;
		overlapSum = 0;
		minOverlap = std::numeric_limits<unsigned>::max();
		maxOverlap = 0;
	}
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

/**
.Function.alignPair:
..summary:Aligns two sequences, returning the score and the align object of the alignment.
..signature:alignPair(ret, seq1, seq2, config, band)
..param.ret:A pair which will be used to store the score of the alignment and the align object.
...type:nolink:Pair<unsigned, Align<TSeq1> >
..param.seq1:The first sequence of the alignment.
...type:nolink:TSeq1
..param.seq2:The second sequence of the alignment.
...type:nolink:TSeq2
..param.config:The alignment configuration used by the Seqan's globalAlignment method.
...type:Class.AlignConfig
..param.band:Bool inidicating if the alignment's lower diagonal should be banded at -2. (Two diagonals
 below the main diagonal.) Useful to speed up computation for 5'-3' overlap alignments.
...type:nolink:bool
..remarks:The main intended uses are AlignConfig<true, false, true, false> with band = true
 and AlignConfig<true, true, true, true> with band = false
..remarks:We use a custom scoring scheme which scores 1 for match, -1 for mismatch, 0 for match with N.
..returns:void
*/
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

/**
.Function.countTotalGaps:
..summary:Returns the total number of gaps in a row object.
..signature:countTotalGaps(row)
..param.row:The gapped row object used by seqan alignments.
...type:Class.Align
..returns:The number of gaps.
...type:nolink:unsigned
*/
template <typename TRow>
unsigned countTotalGaps(TRow& row)
{
	return length(row) - length(source(row));
}

/**
.Function.getOverlap:
..summary:Determines the overlap between two (free-shift) aligned gapless sequences.
..signature:getOverlap(align)
..param.align:An alignment object with two aligned overlapping sequences.
...type:Class.Align
..remarks:The aligned sequences must not contain any internal gaps.
..returns:The number of overlapping positions.
...type:nolink:unsigned
*/
template <typename TAlign>
unsigned getOverlap(TAlign& align)
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	TRow &row1 = seqan::row(align,0);
	TRow &row2 = seqan::row(align,1);
	return length(source(row1)) - countTotalGaps(row2);
}

/**
.Function.getInsertSize:
..summary:Given an alignment of two overlapping (forward and reverse) sequences,
 this function determines the actual size of the insert they cover.
..signature:getInsertSize(align)
..param.align:An alignment object with two aligned overlapping sequences.
...type:Class.Align
..remarks:The aligned sequences must not contain any internal gaps.
..remarks:This method can reconstruct the insert perfectly, provided the alignment
 object represents the real alignment. There can be small discrepancies if indels occurred.
..returns:The insert size that was determined from the overlap alignment.
...type:nolink:unsigned
 */
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

/**
.Function.stripPair:
..summary:Removes adapter contamination from paired-end reads.
..signature:stripPair(seq1, seq2)
..signature:stripPair(seq1, adapter1, seq2, adapter2)
..param.seq1:The forward read.
...type:nolink:String
..param.adapter1:The adapter that contaminates the forward read.
...type:nolink:String
..param.seq2:The backward read.
...type:nolink:String
..param.adapter2:The adapter whose reverse complement contaminates the backward read.
...type:nolink:String
..remarks:This method is very accurate and does not need any knowledge of the specific
 adapter sequences contaminating the reads.
..remarks:The number of trimmed bases is len(seq) - insert, if the sequence was greater than the insert size.
..remarks: The method using adapters can be less accurate than the overload which doesn't use any adapter
 sequence information, since it is more constrained in how it tries to overlap the sequences. Therefore the overload
 using adapter information is currently not in use.
..returns:The determined actual insert size or 0 if the sequences don't overlap.
...type:nolink:unsigned
*/
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
unsigned stripPair(TSeq& seq1, TSeq& adapter1, TSeq& seq2, TSeq& adapter2)
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

/**
.Function.stripPairBatch:
..summary:Removes adapter contamination from paired-end reads and tags their IDs.
..signature:stripPairBatch(set1, idSet1, set2, idSet2, stats, tagOpt)
..param.set1:The set of forward reads.
...type:Class.StringSet
..param.idSet1:The set of IDs associated with the forward reads.
...type:Class.StringSet
..param.set2:The set of backward reads.
...type:Class.StringSet
..param.idSet2:The set of IDs associated with the backward reads.
...type:Class.StringSet
..param.stats:The AdapterTrimmingStats object.
...type:Class.AdapterTrimmingStats
..param.tagOpt:Bool incidating that reads with removed adapter shall be tagged (i.e their IDs).
...type:nolink:bool
..remarks:This method is very accurate and does not need any knowledge of
 the specific adapter sequences contaminating the reads.
..returns:The number of sequences that had some bases removed.
...type:nolink:unsigned
 */
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

/**
.Function.alignAdapter:
..summary:Aligns a sequence to an adapter.
..signature:alignAdapter(ret, alignPair, seq, adapter)
..param.ret:A pair which will be used to store the score of the alignment and the align object.
...type:nolink:Pair<unsigned, Align<TSeq1> >
..param.seq:The sequence whose adapter contamination shall be removed.
...type:Class.String
..param.adapter:The adapter sequence that might contaminate the sequence.
...type:Class.String
..remarks:The aligned sequences must not contain any internal gaps.
..returns:void
 */
template <typename TSeq, typename TAdapter>
void alignAdapter(seqan::Pair<unsigned, seqan::Align<TSeq> >& ret, TSeq& seq, TAdapter& adapter)
{
	// Global free-end alignment. The Alignment configuration specifies that gaps on the
	// start of the read and at the end of the adapter are not penalized -> 5'-3' overlap.
	// We also don't allow gaps by setting the gap penalty high.
	alignPair(ret, seq, adapter, seqan::AlignConfig<true, false, true, false>(), true);
}

/**
.Function.isMatch:
..summary:Checks if a overlap of an alignment is accepted, based on mismatches and the length of the overlap.
..signature:isMatch(overlap, mismatches)
..signature:isMatch(overlap, mismatches, userOptions)
..param.overlap:The number of overlapping positions in the overlap alignment.
...type:nolink:int
..param.mismatches:The number of allowed mismatches in the overlapping region.
...type:nolink:int
..param.userOptions:Parameters specifying requirements for the overlap.
...type:nolink:int
..remarks:The method using two arguments automatically uses a very simple heuristic to determine matches.
..remarks:The overload using "userOptions" lets the user specify the match requirements.
..returns:A bool indicating if the alignment is significant.
...type:nolink:bool
 */
//Version for automatic matching options
inline bool isMatch(int overlap, int mismatches, const Auto &)
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
//Overload for user-definied matching options
inline bool isMatch(int overlap, int mismatches, const User& userOptions)
{
	return overlap >= userOptions.min_length && mismatches <= userOptions.errors;
}

/**
.Function.stripAdapter:
..summary:Removes adapter sequence from a sequence.
..signature:stripAdapter(seq, adapter, spec)
..param.seq:The sequence whose adapter contamination should be removed.
...type:Class.String
..param.adapter:The adapter sequence that might contaminate the sequence.
...type:Class.String
..param.spec:Tag determining the match algorithm that decides whether a match was significant.
...type:nolink:TSpec
..returns:The overlap of the sequence with the adapter.
...type:nolink:unsigned
 */
template <typename TSeq, typename TAdapter, typename TSpec>
unsigned stripAdapter(TSeq& seq, TAdapter& adapter, TSpec const & spec)
{
	typedef seqan::Align<TSeq> TAlign;
	seqan::Pair<unsigned, TAlign> ret;
    alignAdapter(ret, seq, adapter);
	int overlap = getOverlap(ret.i2);
	int score = ret.i1;
	int mismatches = (overlap-score) / 2;
	if (isMatch(overlap, mismatches, spec))
	{
		seqan::erase(seq, length(seq) - overlap, length(seq));
		return overlap;
	}
    else
    {
	    return 0;
    }
}
/**
.Function.stripAdapterBatch:
..summary:Remove adapter sequence from a set of sequences and Tag their IDs.
..signature:stripAdapterBatch(set, idSet, adapter, spec, stats, reverse, tagOpt)
..param.set: A StringSet of sequences whose adapter contaminations shall be removed.
...type:Class.StringSet
..param.idSet: StringSet of IDs associated with the reads.
...type:Class.StringSet
..param.adapter:The adapter sequence that might contaminate the sequences.
...type:Class.String
..param.spec:Tag determining the match algorithm that decides whether a match was significant.
...type:nolink:TSpec
..param.stats: AdapterTrimmingStats object holding statistic information on the Adaptertrimming process.
...type:Class.AdapterTrimmingStats
..param.reverse:Bool indicating that the reverse complement of the sequence shall be used.
...type:nolink:bool
..param.tagOpt:Bool incidating that reads with removed adapter shall be tagged (i.e their IDs).
...type:nolink:bool
..remarks:This method simply applies stripAdapter to all sequences in the set and adds a Tag to the IDs.
..returns:The number of sequences which had bases removed.
...type:nolink:unsigned
..see:Function.stripAdapter
 */
template <typename TSeq, typename TId, typename TAdapter, typename TSpec>
unsigned stripAdapterBatch(seqan::StringSet<TSeq>& set, seqan::StringSet<TId>& idSet, TAdapter& adapter, TSpec const & spec,
		AdapterTrimmingStats& stats, bool reverse = false, bool tagOpt = false)
{
	int t_num = omp_get_max_threads();
	// Create local counting variables to avoid concurrency problems.
	unsigned a_count = 0, overlapSum = 0;
	seqan::String<unsigned> minOverlap;
	seqan::String<unsigned> maxOverlap;
	seqan::resize(minOverlap, t_num, std::numeric_limits<unsigned>::max());
	seqan::resize(maxOverlap, t_num, 0);
	int len = length(set);
	SEQAN_OMP_PRAGMA(parallel for schedule(static) reduction(+:a_count, overlapSum))
	for (int i=0; i < len; ++i)
	{
        unsigned over = stripAdapter(value(set, i), adapter, spec);
		overlapSum += over;
		a_count += (over != 0);
		// Thread saves local min/max seen in the batch it processed.
		int t_id = omp_get_thread_num();
		if (over > 0 && over < minOverlap[t_id])
        {
            minOverlap[t_id] = over;
        }
		if (over > 0 && over > maxOverlap[t_id])
        {
            maxOverlap[t_id] = over;
        }
        if (tagOpt && over != 0)
        {
            append(idSet[i], "[AdapterRemoved]");
        }
	}
	if (!reverse)
    {
		stats.a1count += a_count;
    }
	else
    {
		stats.a2count += a_count;
    }
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
	return a_count;
}

/**
.Function.stripReverseAdapterBatch:
..summary:A simple interface to align the reverse complement of an adapter to a set of sequences
 and removing significant matches.
..signature:stripReverseAdapterBatch(set, idSet, adapter, spec, stats, tagOpt)
..param.set: A StringSet of sequences whose adapter contaminations shall be removed.
...type:Class.StringSet
..param.idSet: StringSet of IDs associated with the reads.
...type:Class.StringSet
..param.adapter:The adapter sequence whose reverse complement might contaminate the sequences.
...type:Class.String
..param.spec:Tag determining the match algorithm that decides whether a match was significant.
...type:nolink:TSpec
..param.stats: AdapterTrimmingStats object holding statistic information on the Adaptertrimming process.
...type:Class.AdapterTrimmingStats
..param.tagOpt:Bool incidating that reads with removed adapter shall be tagged (i.e their IDs).
...type:nolink:bool
..remarks:This interface can be used to detect contamination in reverse reads of paired-end
 reads. Those reads might read into the reverse complement of an adapter.
..returns:The number of sequences which had bases removed.
...type:nolink:unsigned
..see:Function.stripAdapter
..see:Function.stripAdapterBatch
 */
template <typename TSeq, typename TId, typename TSpec>
unsigned stripReverseAdapterBatch(seqan::StringSet<TSeq>& set, seqan::StringSet<TId>& IdSet, TSeq& adapter, TSpec const & spec,
		AdapterTrimmingStats& stats, bool tagOpt)
{
	typedef typename Value<TSeq>::Type TAlphabet;
	typedef typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type TReverseComplement;
	TReverseComplement mod(adapter);
	return stripAdapterBatch(set, IdSet, mod, spec, stats, true, tagOpt);
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_
