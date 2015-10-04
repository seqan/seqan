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
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// ==========================================================================W

#ifndef ADAPTERTRIMMING_H
#define ADAPTERTRIMMING_H

#include <seqan/align.h>
#include "helper_functions.h"
#include "general_stats.h"

// ============================================================================
// Metafunctions
// ============================================================================


//brief A metafunction which constructs a ModifiedString type for an Alphabet.
template <class TValue>
struct STRING_REVERSE_COMPLEMENT
{
    typedef seqan::ModifiedString<
        seqan::ModifiedString<	seqan::String<TValue>, seqan::ModView<seqan::FunctorComplement<TValue> > >,
        seqan::ModReverse>	Type;
};


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

using TAdapterAlphabet = seqan::Dna5Q;
using TAdapterSequence = seqan::String<TAdapterAlphabet>;

struct AdapterItem;
typedef std::vector< AdapterItem > TAdapterSet;
using TReverseComplement = STRING_REVERSE_COMPLEMENT<TAdapterAlphabet>::Type;

struct AdapterItem
{
    AdapterItem() : anchored(false) {};
    AdapterItem(const TAdapterSequence &adapter) : seq(adapter), anchored(false) {};
    AdapterItem(const TAdapterSequence &adapter, const bool anchored) : seq(adapter), anchored(anchored) {};

    AdapterItem getReverseComplement() const noexcept
    {
        auto seqCopy = seq;
        return AdapterItem(TReverseComplement(seqCopy), anchored);
    }

    TAdapterSequence seq;
    // the anchored matching mode should have the same functionality as cutadapt's anchored mode
    // see http://cutadapt.readthedocs.org/en/latest/guide.html
    // In anchored mode, the adapter has to start (3" adapter) or has to end (5" adapter) with the sequence
    // The anchored mode is rarely used, at least for 3" adapters
    // Todo: implement rooted mode
    bool anchored;
};

// Define scoring function type.
typedef seqan::Score<int, seqan::ScoreMatrix<seqan::Dna5, seqan::AdapterScoringMatrix> > TScore;

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


// ============================================================================
// Functions
// ============================================================================


template <typename TRow>
unsigned countTotalGaps(TRow& row) noexcept
{
	return length(row) - length(source(row));
}

template <typename TAlign>
unsigned getOverlap(TAlign& align) noexcept
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	TRow &row1 = seqan::row(align,0);
	TRow &row2 = seqan::row(align,1);
	return length(source(row1)) - countTotalGaps(row2);
}

template <typename TAlign>
unsigned getInsertSize(TAlign& align) noexcept
{
	typedef typename seqan::Row<TAlign>::Type TRow;
	const TRow &row1 = seqan::row(align,0);
	const TRow &row2 = seqan::row(align,1);
	const unsigned seq1_length = length(source(row1));
	const unsigned seq2_length = length(source(row2));
	// Calculate overlap and overhangs.
	const unsigned overlap = seq1_length - countTotalGaps(row2);
	const unsigned seq2l = seqan::countGaps(seqan::begin(row1)); // Overhang of sequence 2 = Gaps at start of row 1.
	const unsigned seq1r = countTotalGaps(row2) - seqan::countGaps(seqan::begin(row2)); // Overhang of sequence 1 = Gaps at end of row 2.
	// Insert size: Add sequence lengths, subtract common region (overlap)
	// and subtract overhangs left and right of insert.
	return seq1_length + seq2_length - overlap - (seq2l + seq1r);
}


template <typename TSeq1, typename TSeq2, bool TTop, bool TLeft, bool TRight, bool TBottom>
void alignPair(std::pair<unsigned, seqan::Align<TSeq1> >& ret, const TSeq1& seq1, const TSeq2& seq2,
		const seqan::AlignConfig<TTop, TLeft, TRight, TBottom>& config, bool band = false) noexcept
{
    seqan::resize(rows(ret.second), 2);
    seqan::assignSource(row(ret.second, 0), seq1);
    seqan::assignSource(row(ret.second, 1), seq2);
	// Overlap alignment. We again don't allow gaps, since they are unlikely in Illumina data.
	// We mostly use AlignConfigs true,false,true,false and true,true,true,true here.
	const TScore adapterScore(-100);
	// If we can do a banded alignment, restrict ourselves to the upper half of the matrix.
    if (band)
    {
        ret.first = globalAlignment(ret.second, adapterScore, config, -2, length(seq1), seqan::LinearGaps());
    }
    else
    {
        ret.first = globalAlignment(ret.second, adapterScore, config, seqan::LinearGaps());
    }
}

template <typename TSeq>
unsigned stripPair(TSeq& seq1, TSeq& seq2) noexcept
{
    // When aligning the two sequences, the complementary sequence is reversed and
    // complemented, so we have an overlap alignment with complementary bases being the same.
    typedef typename seqan::Value<TSeq>::Type TAlphabet;
    typedef typename STRING_REVERSE_COMPLEMENT<TAlphabet>::Type TReverseComplement;
    TReverseComplement mod(seq2);
    typedef seqan::Align<TSeq> TAlign;
    std::pair<unsigned, TAlign> ret;
    alignPair(ret, seq1, mod, seqan::AlignConfig<true, true, true, true>());
    const auto& score = ret.first;
    const TAlign& align = ret.second;
    // Use the overlap of the two sequences to determine the end position.
    const unsigned overlap = getOverlap(align);
    const unsigned mismatches = (overlap - score) / 2;
    // We require a certain correct overlap to exclude spurious hits.
    // (Especially reverse 3'->5' alignments not caused by adapters.)
    if (overlap <= 5 || mismatches > overlap * 0.15)
    {
        return 0;
    }
    // Get actual size of the insert (possible to determine from overlap etc.).
    const unsigned insert = getInsertSize(align);
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

template <typename TSeq, typename TAdapterItem>
void alignAdapter(std::pair<unsigned, seqan::Align<TSeq> >& ret, const TSeq& seq, TAdapterItem const& adapterItem) noexcept
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
inline bool isMatch(const int overlap, const int mismatches, const AdapterMatchSettings &adatperMatchSettings) noexcept
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

enum adapterDirection : bool
{
    reverse, 
    forward
};

template<bool _val>
struct TagAdapter
{
    static const bool value = _val;
};

template<bool _direction>
struct StripAdapterDirection
{
    static const bool value = _direction;
};

template <typename TSeq, typename TAdapters, typename TSpec, typename TStripAdapterDirection>
unsigned stripAdapter(TSeq& seq, AdapterTrimmingStats& stats, TAdapters const& adapters, TSpec const& spec,
    const TStripAdapterDirection&) noexcept
{
    typedef typename seqan::Align<TSeq> TAlign;
    typedef typename seqan::Row<TAlign>::Type TRow;
    std::pair<unsigned, TAlign> ret;

    unsigned removed{ 0 };

    for (auto const& adapterItem : adapters)
    {
        //std::cout << "seq: " << seq << std::endl;
        //std::cout << "adapter: " << adapterItem.seq << std::endl;
        if (length(adapterItem.seq) == 0)
            continue;

        if (TStripAdapterDirection::value == adapterDirection::reverse)
        {
            alignAdapter(ret, seq, adapterItem.getReverseComplement());    // align crashes if seq is an empty string!
        }
        else
            alignAdapter(ret, seq, adapterItem);    // align crashes if seq is an empty string!

        const unsigned int overlap = getOverlap(ret.second);
        const int score = ret.first;
        const int mismatches = (overlap - score) / 2;
        if (isMatch(overlap, mismatches, spec))
        {
            //std::cout << "score: " << ret.i1 << " overlap: " << overlap << " mismatches: " << mismatches << std::endl;
            //std::cout << ret.i2 << std::endl;
            const TRow row2 = row(ret.second, 1);
            //std::cout << "adapter start position: " << toViewPosition(row2, 0) << std::endl;
            const auto seqLen = length(seq);
            removed += seqLen - toViewPosition(row2, 0);
            seqan::erase(seq, toViewPosition(row2, 0), seqLen);
            //std::cout << "stripped seq: " << seq << std::endl;
            stats.overlapSum += overlap;
            if (TStripAdapterDirection::value == adapterDirection::forward)
                ++stats.a2count;
            else
                ++stats.a1count;

            stats.maxOverlap = stats.maxOverlap > overlap ? stats.maxOverlap : overlap;
            stats.minOverlap = stats.minOverlap < overlap ? stats.minOverlap : overlap;
            if (spec.singleMatch || empty(seq))
                return removed;
        }
    }
    return removed;
}

template < template <typename> class TRead, typename TSeq, typename TAdaptersArray, typename TSpec, typename TTagAdapter,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value> >
unsigned stripAdapterBatch(std::vector<TRead<TSeq>>& reads, TAdaptersArray const& adapters, TSpec const& spec, const bool pairedNoAdapterFile,
    AdapterTrimmingStats& stats, TTagAdapter, bool = false) noexcept(!TTagAdapter::value)
{
    (void)pairedNoAdapterFile;
    for(auto& read: reads)
    {
        if (seqan::empty(read.seq))
            continue;
        const unsigned over = stripAdapter(read.seq, stats, adapters[1], spec, StripAdapterDirection<adapterDirection::forward>());
        if (TTagAdapter::value && over != 0)
            insertAfterFirstToken(read.id, ":AdapterRemoved");
    }
    return stats.a1count + stats.a2count;
}

// pairedEnd adapters will be trimmed in single mode, each seperately
template < template <typename> class TRead, typename TSeq, typename TAdaptersArray, typename TSpec, typename TTagAdapter,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq>>::value> >
    unsigned stripAdapterBatch(std::vector<TRead<TSeq>>& reads, TAdaptersArray const& adapters, TSpec const& spec, const bool pairedNoAdapterFile,
        AdapterTrimmingStats& stats, TTagAdapter) noexcept(!TTagAdapter::value)
{
    for (auto& read : reads)
    {
        if (seqan::empty(read.seq))
            continue;
        unsigned over = 0;
        if (pairedNoAdapterFile)
        {
            stripPair(read.seq, read.seqRev);
        }
        else
        {
            over = stripAdapter(read.seq, stats, adapters[1], spec, StripAdapterDirection<adapterDirection::forward>());
            if (!seqan::empty(read.seqRev))
                over += stripAdapter(read.seqRev, stats, adapters[0], spec, StripAdapterDirection<adapterDirection::reverse>());
        }
        if (TTagAdapter::value && over != 0)
            insertAfterFirstToken(read.id, ":AdapterRemoved");
    }
    return stats.a1count + stats.a2count;
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_ADAPTERTRIMMING_H_
