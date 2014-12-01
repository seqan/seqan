// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "overlapper.h"

#include <seqan/align.h>
#include <seqan/score.h>

namespace {  // anonymous namespace

// --------------------------------------------------------------------------
// Function _fixBandSize()
// --------------------------------------------------------------------------

// Fix band size for the given alignment.

template <typename TSequenceH, typename TSequenceV, typename TAlignConfig, typename TAlgoTag>
void _fixBandSize(int & lDiag,
                  int & uDiag,
                  TSequenceH const & seqH,
                  TSequenceV const & seqV,
                  TAlignConfig const & /*alignConfig*/,
                  TAlgoTag const & /*algoTag*/)
{
    using namespace seqan;

    // typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef typename If<typename IsSameType<TAlgoTag, Gotoh>::Type, AffineGaps, LinearGaps>::Type TGapsType;
    typedef typename SetupAlignmentProfile_<DPGlobal, TAlignConfig, TGapsType, TracebackConfig_<SingleTrace, GapsLeft> >::Type TDPProfile;

    if (uDiag < -(int)length(seqV))
        uDiag = -(int)length(seqV);
    if (lDiag > (int)length(seqH))
        lDiag = length(seqV);

    if (uDiag < 0 && !IsFreeEndGap_<TDPProfile, DPFirstColumn>::VALUE)
        uDiag = 0;

    if (lDiag > 0 && !IsFreeEndGap_<TDPProfile, DPFirstRow>::VALUE)
        lDiag = 0;

    if (uDiag + (int)length(seqV) < (int)length(seqH) && !IsFreeEndGap_<TDPProfile, DPLastRow>::VALUE)
        uDiag = (int)length(seqH) - (int)length(seqV);

    if (lDiag + (int)length(seqV) > (int)length(seqH) && !IsFreeEndGap_<TDPProfile, DPLastColumn>::VALUE)
        lDiag = (int)length(seqH) - (int)length(seqV);
}

}  // anonymous namespace

namespace assembler {

// ----------------------------------------------------------------------------
// Class OverlapperImpl
// ----------------------------------------------------------------------------

class OverlapperImpl
{
public:
    OverlapperImpl(OverlapperOptions options = OverlapperOptions()) : options(options)
    {}

    bool computeOverlap(Overlap & overlap,
                        TFragments & alignment,
                        seqan::Dna5String const & seqH,
                        seqan::Dna5String const & seqV,
                        OverlapCandidate const & candidate) const;

    bool computeOverlap(OverlapStore & store,
                        seqan::Dna5String const & seqH,
                        seqan::Dna5String const & seqV,
                        OverlapCandidate const & candidate) const;

    Overlap computeOverlap(seqan::Dna5String const & seqH,
                           seqan::Dna5String const & seqV,
                           OverlapCandidate const & candidate) const;
    
private:
    // Generate Overlap record from the alignment stored in fragments.  Length information is taken from seqs.
    Overlap overlapFromAlignment(
            seqan::String<seqan::Fragment<>> const & fragments,
            seqan::StringSet<seqan::Dna5String, seqan::Dependent<>> const & strings) const;

    // Configuration to use.
    OverlapperOptions options;
};

Overlap OverlapperImpl::overlapFromAlignment(seqan::String<seqan::Fragment<> > const & fragments,
                                             seqan::StringSet<seqan::Dna5String, seqan::Dependent<>> const & strings) const
{
    typedef seqan::StringSet<seqan::Dna5String, seqan::Dependent<> > TStringSet;
    TStringSet & stringsNC = const_cast<seqan::StringSet<seqan::Dna5String, seqan::Dependent<> > &>(strings);

    if (options.logging)
    {
        std::cerr << "FRAGMENTS\n";
        for (unsigned i = 0; i < length(fragments); ++i)
            std::cerr << "  Fragment(" << fragments[i].seqId1 << ", " << fragments[i].begin1
                      << ", " << fragments[i].seqId2 << ", " << fragments[i].begin2 << ", "
                      << fragments[i].len << ")\n";
    }

    // TODO(holtgrew): overlap length should actually be the length of the alignment

    auto frags = fragments;
    std::sort(begin(frags, seqan::Standard()), end(frags, seqan::Standard()));

    auto frag0 = front(frags);  // first
    auto id0 = sequenceId(frag0, 0);
    auto id1 = sequenceId(frag0, 1);
    auto len0 = length(strings[0]);
    auto len1 = length(strings[1]);

    typedef decltype(+fragmentBegin(frag0, 0)) TPos;
    std::pair<TPos, TPos> range0(seqan::maxValue<TPos>(), seqan::minValue<TPos>());
    auto range1 = range0;
    int errors = 0;
    for (auto itF = begin(frags, seqan::Standard()); itF != end(frags, seqan::Standard()); ++itF)
    {
        // Get some shortcuts.
        auto fLen = fragmentLength(*itF);
        auto begin0 = fragmentBegin(*itF, id0);
        auto begin1 = fragmentBegin(*itF, id1);

        // Count indels.
        if (itF != begin(frags, seqan::Standard()))
        {
            SEQAN_ASSERT_LEQ(range0.second, begin0);
            SEQAN_ASSERT_LEQ(range1.second, begin1);
            SEQAN_ASSERT_NEQ((begin0 != range0.second), (begin1 != range1.second));
            errors += (begin0 - range0.second);
            errors += (begin1 - range1.second);
        }

        // Update begin/end position in either read.
        range0.first = std::min(range0.first, begin0);
        range0.second = std::max(range0.second, begin0 + fLen);
        range1.first = std::min(range1.first, begin1);
        range1.second = std::max(range1.second, begin1 + fLen);

        // Count matches/mismatches.
        auto label0 = label(*itF, stringsNC, id0);
        auto label1 = label(*itF, stringsNC, id1);
        for (auto it0 = begin(label0, seqan::Standard()), it1 = begin(label1, seqan::Standard());
             it0 != end(label0, seqan::Standard()); ++it0, ++it1)
            errors += ((seqan::Dna5)*it0 == 'N' ||
                       (seqan::Dna5)*it1 == 'N' ||
                       (seqan::Dna5)*it0 != (seqan::Dna5)*it1);
    }

    // In case that the alignment to the right aligns to a gap, flush left.
    int delta = std::min(range0.first, range1.first);
    range0.first -= delta;
    range1.first -= delta;

    SEQAN_ASSERT_MSG(range0.first == 0 || range1.first == 0, "One must start at beginning");
    // NB: Do not activate the following, does not have to be true, can end in alignment to gap.
    // SEQAN_ASSERT_MSG(range0.second == len0 || range1.second == len1, "One must end at last");

    auto begin1 = range0.first, begin0 = range1.first;
    // int overlapLen = std::max(range0.second - range0.first, range1.second - range1.first);

    if (options.logging)
    {
        std::cerr << "range0 = (" << range0.first << ", " << range0.second << ")\n"
                  << "range1 = (" << range1.first << ", " << range1.second << ")\n";
    }

    return Overlap(id0, id1, len0, len1, begin0, begin1, errors);
}

bool OverlapperImpl::computeOverlap(Overlap & overlap,
                                    TFragments & frags,
                                    seqan::Dna5String const & seqH,
                                    seqan::Dna5String const & seqV,
                                    OverlapCandidate const & candidate) const
{
    clear(frags);

    if (options.logging)
        std::cerr << "Computing overlap\n"
                  << "  seqH: " << seqH << "\n"
                  << "  seqV: " << seqV << "\n"
                  << "  cand: " << candidate << "\n";

    seqan::StringSet<seqan::Dna5String, seqan::Dependent<>> pairSet;
    appendValue(pairSet, const_cast<seqan::Dna5String &>(seqH));  // id 0
    appendValue(pairSet, const_cast<seqan::Dna5String &>(seqV));  // id 1

    seqan::Score<int, seqan::Simple> scoringScheme(1000, -1000, -1001);

    seqan::AlignConfig<true, true, true, true> alignConfig;
    int uDiag = candidate.uDiag;
    int lDiag = candidate.lDiag;

    _fixBandSize(lDiag, uDiag, pairSet[0], pairSet[1], alignConfig, seqan::Gotoh());

    // DEBUG
    if (options.logging)
    {
        std::cerr << "\n\n(alignment of " << candidate << ")\n"
                  << "0:\t" << pairSet[0] << "\n"
                  << "1:\t" << pairSet[1] << "\n";
        typedef seqan::Value<seqan::StringSet<seqan::Dna5String>>::Type TStringSeq;
        seqan::Align<TStringSeq> align;
        resize(rows(align), 2);
        setSource(row(align, 0), pairSet[0]);
        setSource(row(align, 1), pairSet[1]);
        globalAlignment(align, scoringScheme, alignConfig, lDiag, uDiag, seqan::Gotoh());
        std::cerr << "\n" << align << "\n";
    }
    // /DEBUG

    int overlapScore = globalAlignment(frags, pairSet, scoringScheme, alignConfig, lDiag, uDiag,
                                       seqan::NeedlemanWunsch());
    (void)overlapScore;
    if (empty(frags))
        return false;

    overlap = overlapFromAlignment(frags, pairSet);

    // Replace ids in fragments and overlap.
    overlap.seq0 = candidate.seq0;
    overlap.seq1 = candidate.seq1;
    std::for_each(begin(frags, seqan::Standard()), end(frags, seqan::Standard()),
                   [candidate](seqan::Fragment<> & frag) {
                      sequenceId(frag, 0) = candidate.seq0;
                      sequenceId(frag, 1) = candidate.seq1;
                  });

    int ovlLen = overlap.length();
    bool ok = ((options.overlapMinLength < 0 || ovlLen >= options.overlapMinLength) &&
               ((options.overlapErrorRate < 0 || ovlLen == 0)
                || (1.0 * overlap.errors / ovlLen) <= options.overlapErrorRate + 0.00001));
    if (options.logging)
        std::cerr << "Resulting overlap:\t" << overlap << " (passes quality? " << ok << ")\n";
    return ok;
}

bool OverlapperImpl::computeOverlap(OverlapStore & store,
                                    seqan::Dna5String const & seqH,
                                    seqan::Dna5String const & seqV,
                                    OverlapCandidate const & candidate) const
{
    Overlap ovl;
    TFragments ali;
    bool result = computeOverlap(ovl, ali, seqH, seqV, candidate);
    if (result)
    {
        if (options.logging)
            std::cerr << "Inserting into store: " << ovl.normalize() << "\n";
        store.insert(ovl, ali);
    }
    return result;
}

Overlap OverlapperImpl::computeOverlap(seqan::Dna5String const & seqH,
                                       seqan::Dna5String const & seqV,
                                       OverlapCandidate const & candidate) const
{
    Overlap ovl;
    TFragments ali;
    computeOverlap(ovl, ali, seqH, seqV, candidate);

    return ovl;
}

// ----------------------------------------------------------------------------
// Class Overlapper
// ----------------------------------------------------------------------------

Overlapper::Overlapper(OverlapperOptions options) : impl(new OverlapperImpl(options))
{}

Overlapper::~Overlapper()
{}

bool Overlapper::computeOverlap(Overlap & overlap,
                                TFragments & alignment,
                                seqan::Dna5String const & seqH,
                                seqan::Dna5String const & seqV,
                                OverlapCandidate const & candidate) const
{
    return impl->computeOverlap(overlap, alignment, seqH, seqV, candidate);
}

bool Overlapper::computeOverlap(OverlapStore & store,
                                seqan::Dna5String const & seqH,
                                seqan::Dna5String const & seqV,
                                OverlapCandidate const & candidate) const
{
    return impl->computeOverlap(store, seqH, seqV, candidate);
}

Overlap Overlapper::computeOverlap(seqan::Dna5String const & seqH,
                                   seqan::Dna5String const & seqV,
                                   OverlapCandidate const & candidate) const
{
    return impl->computeOverlap(seqH, seqV, candidate);
}

}  // namespace assembler
