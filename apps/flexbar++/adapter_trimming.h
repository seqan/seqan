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

struct AdapterMatchSettings
{
    AdapterMatchSettings(const int m, const int e, const double er, const unsigned int lo, const unsigned int times) : min_length(m), errors(e), leftOverhang(lo), errorRate(er), erMode(false), modeAuto(false), times(times)
    {
        erMode = ((e == 0) && (er != 0));
    }
    AdapterMatchSettings() : min_length(0), errors(0), leftOverhang(0), errorRate(0), erMode(false), modeAuto(true), times(1) {};
   
    int min_length; //The minimum length of the overlap.
	int errors;     //The maximum number of errors we allow.
    unsigned int leftOverhang;
	double errorRate;  //The maximum number of errors allowed per overlap
	bool erMode;
    bool modeAuto;
    unsigned int times;
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


template <typename TSeq, typename TAdapter, bool TTop, bool TLeft, bool TRight, bool TBottom>
void alignPair(std::pair<unsigned, seqan::Align<TSeq> >& ret, const TSeq& seq1, const TAdapter& seq2,
    const seqan::AlignConfig<TTop, TLeft, TRight, TBottom>& config, const unsigned int leftOverhang, const unsigned int min_overlap) noexcept
{
    seqan::resize(rows(ret.second), 2);
    seqan::assignSource(row(ret.second, 0), seq1);
    seqan::assignSource(row(ret.second, 1), seq2);

    // dont allow gaps by setting the gap penalty high
    const TScore adapterScore(-100);

    // banded alignment
    ret.first = globalAlignment(ret.second, adapterScore, config, -static_cast<int>(leftOverhang), length(seq1)-min_overlap, seqan::LinearGaps());
}

// used only for testing and paired end data
template <typename TSeq, typename TAdapter, bool TTop, bool TLeft, bool TRight, bool TBottom>
void alignPair(std::pair<unsigned, seqan::Align<TSeq> >& ret, const TSeq& seq1, const TAdapter& seq2,
		const seqan::AlignConfig<TTop, TLeft, TRight, TBottom>& config) noexcept
{
    seqan::resize(rows(ret.second), 2);
    seqan::assignSource(row(ret.second, 0), seq1);
    seqan::assignSource(row(ret.second, 1), seq2);

    // dont allow gaps by setting the gap penalty high
    const TScore adapterScore(-100);

	// global alignment, not banded
    ret.first = globalAlignment(ret.second, adapterScore, config, seqan::LinearGaps());
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

//Version for automatic matching options
inline bool isMatch(const int overlap, const int mismatches, const AdapterMatchSettings &adatperMatchSettings) noexcept
{
    if (overlap == 0)
        return false;
    if (adatperMatchSettings.erMode)
        return overlap >= adatperMatchSettings.min_length && (static_cast<double>(mismatches) / static_cast<double>(overlap)) <= adatperMatchSettings.errorRate;
    else
        return overlap >= adatperMatchSettings.min_length && mismatches <= adatperMatchSettings.errors;
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

template <typename TSeq, typename TAdapters, typename TStripAdapterDirection>
unsigned stripAdapter(TSeq& seq, AdapterTrimmingStats& stats, TAdapters const& adapters, AdapterMatchSettings const& spec,
    const TStripAdapterDirection&) noexcept
{
    typedef typename seqan::Align<TSeq> TAlign;
    typedef typename seqan::Row<TAlign>::Type TRow;
    std::pair<unsigned, TAlign> ret;

    unsigned removed{ 0 };
    std::vector<std::tuple<unsigned int, TRow, unsigned int>> matches;
    matches.reserve(spec.times);
    for (unsigned int n = 0;n < spec.times; ++n)
    {
        for (auto const& adapterItem : adapters)
        {
            //std::cout << "seq: " << seq << std::endl;
            //std::cout << "adapter: " << adapterItem.seq << std::endl;
            if (length(adapterItem.seq) == 0)
                continue;

            // always use banded alignment
            if (TStripAdapterDirection::value == adapterDirection::reverse)
            {
                alignPair(ret, seq, adapterItem.getReverseComplement().seq, seqan::AlignConfig<true, true, true, true>(), spec.leftOverhang, spec.min_length);    // align crashes if seq is an empty string!
            }
            else
                alignPair(ret, seq, adapterItem.seq, seqan::AlignConfig<true, true, true, true>(), spec.leftOverhang, spec.min_length);    // align crashes if seq is an empty string!

            const int score = ret.first;
            if (score < 0)
                continue;
            const unsigned int overlap = getOverlap(ret.second);
            const int mismatches = (overlap - score) / 2;
            //std::cout << "score: " << ret.first << " overlap: " << overlap << " mismatches: " << mismatches << std::endl;
            //std::cout << ret.second << std::endl;
            //std::cout << "adapter start position: " << toViewPosition(row(ret.second, 1), 0) << std::endl;
            if (isMatch(overlap, mismatches, spec))
            {
                const TRow row2 = row(ret.second, 1);
                matches.push_back(std::make_tuple(score,row2, overlap));
            }
        }
        if (matches.empty())
            return removed;

        // find best matching adapter
        auto maxIt = matches.cbegin();
        for (auto it = matches.cbegin(); it != matches.cend();++it)
            if (std::get<0>(*it)>std::get<0>(*maxIt))
                maxIt = it;

        // erase best matching adapter from sequence
        const auto seqLen = length(seq);
        const auto erasePosition = toViewPosition(std::get<1>(*maxIt), 0);
        const auto overlap = std::get<2>(*maxIt);
        removed += seqLen - erasePosition;
        //std::cout << "unstripped seq: " << seq << std::endl;
        seqan::erase(seq, erasePosition, seqLen);
        //std::cout << "stripped seq  : " << seq << std::endl;
        stats.overlapSum += overlap;
        if (TStripAdapterDirection::value == adapterDirection::forward)
            ++stats.a2count;
        else
            ++stats.a1count;

        stats.maxOverlap = std::max(stats.maxOverlap, overlap);
        stats.minOverlap = std::min(stats.minOverlap, overlap);
        if (empty(seq))
            return removed;
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
