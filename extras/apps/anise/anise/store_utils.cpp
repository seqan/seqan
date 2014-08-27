// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#include "store_utils.h"

#include <seqan/realign.h>

#include "realignment.h"
#include "sliding_window.h"

// --------------------------------------------------------------------------
// Class DifferingReadsRemover
// --------------------------------------------------------------------------

// Removes reads from fragment store that have a too high error rate compared to the contig.
//
// These reads are moved to their own contig (with contigStore and contigNameStore resized appropriately.  Note that we
// we also remove gaps and shift alignments appropriately.

class DifferingReadsRemoverImpl
{
public:

    DifferingReadsRemoverImpl(double maxErrorRate) : maxErrorRate(maxErrorRate)
    {}

    // The IDs of the removed reads are written to readIDs.
    //
    // Will resort the read alignments.
    void run(seqan::FragmentStore<void, MyStoreConfig> & store,
             std::vector<unsigned> & readIDs) const
    {
        // Compute coverage on each contig, so we can comfortably shift reads and gaps later.
        std::vector<std::vector<unsigned> > coverages;
        computeCoverages(coverages, store);

        // Move reads alignments with an error rate of more than maxErrorRate to their own contig.
        removeDifferingReads(readIDs, coverages, store);

        // Remove paddings, leading, and trailing gaps.
        //
        // Afterwards, reads will be sorted by (contig id, begin pos).
        updateAtZeroCoverage(store, coverages);

        // Build profiles form read alignments and recreated contigs.
        fixContigs(store);
    }

private:

    // Write coverage at position pos on contig rID into coverages[rID][pos].
    void computeCoverages(std::vector<std::vector<unsigned> > & coverages,
                          TFragmentStore /*const*/ & store) const;

    // Update coverage for given contig with given read element.
    void updateCoverage(std::vector<unsigned> & coverage,
                        TAlignedReadStoreElement /*const*/ & el,
                        TFragmentStore /*const*/ store,
                        int delta) const;

    // Remove differing reads, updating the coverage in coverages.
    void removeDifferingReads(std::vector<unsigned> & readIDs,
                              std::vector<std::vector<unsigned> > & coverages,
                              TFragmentStore /*const*/ & store) const;

    // Update the store at zero coverage positions by shifting alignments and removing padding gaps.
    void updateAtZeroCoverage(TFragmentStore /*const*/ & store,
                              std::vector<std::vector<unsigned> > & coverages) const;

    // true if removed gap that is not leading/trailing gap
    template <typename TPos, typename TGapAnchor, typename TSpec, typename TGapPos>
    bool removeGap2(seqan::AlignedReadStoreElement<TPos, TGapAnchor, TSpec> & alignedRead,
                    TGapPos const gapPos) const;

    // Returns number of removed gaps.
    template <typename TAlignedReads, typename TGapPos>
    int removeGap(TAlignedReads & alignedReadStore, TGapPos const gapPos) const;

    // For each contig: Build contig profile and update sequence and contig gaps for this.
    void fixContigs(TFragmentStore /*const*/ & store) const;

    // Maximal read error rate.
    double maxErrorRate;
};

void DifferingReadsRemoverImpl::computeCoverages(
        std::vector<std::vector<unsigned> > & coverages,
        TFragmentStore /*const*/ & store) const
{
    // Initialize coverages.
    coverages.resize(length(store.contigStore));
    for (unsigned i = 0; i < length(store.contigStore); ++i)
    {
        typedef decltype(store.contigStore[0].seq) TContigSeq;
        typedef typename std::remove_reference<decltype(store.contigStore[0].gaps)>::type TContigGapAnchors;
        typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<TContigGapAnchors> > TContigGaps;
        TContigGaps contigGaps(store.contigStore[i].seq, store.contigStore[i].gaps);

        coverages[i].resize(length(contigGaps), 0);
    }

    for (auto /*const*/ & el : store.alignedReadStore)
        updateCoverage(coverages[el.contigId], el, store, +1);
}

void DifferingReadsRemoverImpl::updateCoverage(
        std::vector<unsigned> & coverage,
        TAlignedReadStoreElement /*const*/ & el,
        TFragmentStore /*const*/ store,
        int delta) const
{
    typedef typename std::remove_reference<decltype(store.readSeqStore[0])>::type TReadSeq;
    typedef typename std::remove_reference<decltype(el.gaps)>::type TReadGapAnchors;
    typedef seqan::Gaps<TReadSeq, seqan::AnchorGaps<TReadGapAnchors> > TReadGaps;

    // Create read gaps for coverage computation.
    TReadGaps readGaps(store.readSeqStore[el.readId], el.gaps);
    int pos = std::min(el.beginPos, el.endPos);
    for (auto it = begin(readGaps, seqan::Standard()), itEnd = end(readGaps, seqan::Standard());
         it != itEnd; ++it, ++pos)
        coverage.at(pos) += delta * !isGap(it);
}

void DifferingReadsRemoverImpl::removeDifferingReads(
        std::vector<unsigned> & readIDs,
        std::vector<std::vector<unsigned> > & coverages,
        TFragmentStore /*const*/ & store) const
{
    unsigned nextContigID = length(store.contigStore);

    // Count alignments on the contigs.
    seqan::String<unsigned> numAlis;
    resize(numAlis, length(store.contigStore), 0);
    for (auto it = begin(store.alignedReadStore, seqan::Rooted()); !atEnd(it); ++it)
        numAlis[it->contigId] += 1;

    // Distribute non-singletons to new contigs.
    for (auto it = begin(store.alignedReadStore, seqan::Rooted()); !atEnd(it); ++it)
    {
        if (numAlis[it->contigId] == 1u)
            continue;  // Skip if singleton.
        if (countErrors(store, *it) > (int)floor(maxErrorRate * length(store.readSeqStore[it->readId])))
        {
            updateCoverage(coverages[it->contigId], *it, store, -1);

            readIDs.push_back(it->readId);
            it->contigId = nextContigID++;
            it->beginPos = 0;
            it->endPos = length(store.readSeqStore[it->readId]);
            clear(it->gaps);  // remove gaps

            resize(store.contigStore, nextContigID);
            back(store.contigStore).seq = store.readSeqStore[it->readId];  // contig seq == read seq
            resize(store.contigNameStore, nextContigID);

            // Creage coverage = 1 for new contig.
            coverages.resize(coverages.size() + 1);
            coverages.back().resize(it->endPos, 1);
        }
    }
}

// Update the store at zero coverage positions by shifting alignments and removing padding gaps.
void DifferingReadsRemoverImpl::updateAtZeroCoverage(
        TFragmentStore /*const*/ & store,
        std::vector<std::vector<unsigned> > & coverages) const
{
    sortAlignedReads(store.alignedReadStore, seqan::SortBeginPos());
    sortAlignedReads(store.alignedReadStore, seqan::SortContigId());

    for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
    {
        auto itBegin = lowerBoundAlignedReads(store.alignedReadStore, contigID, seqan::SortContigId());
        auto itEnd = upperBoundAlignedReads(store.alignedReadStore, contigID, seqan::SortContigId());
        auto alignedReads = infix(store.alignedReadStore,
                                  itBegin - begin(store.alignedReadStore, seqan::Standard()),
                                  itEnd - begin(store.alignedReadStore, seqan::Standard()));
        int pos = coverages[contigID].size() - 1;
        for (auto it = coverages[contigID].rbegin(); it != coverages[contigID].rend(); ++it, --pos)
            if (*it == 0u)
                removeGap(alignedReads, pos);
    }
}

template <typename TPos, typename TGapAnchor, typename TSpec, typename TGapPos>
bool DifferingReadsRemoverImpl::removeGap2(
        seqan::AlignedReadStoreElement<TPos, TGapAnchor, TSpec> & alignedRead,
        TGapPos const gapPos) const
{
    if (gapPos < (TGapPos) alignedRead.beginPos)
    {
        --alignedRead.beginPos; --alignedRead.endPos;
    }
    else if (gapPos < (TGapPos) alignedRead.endPos)
    {
        --alignedRead.endPos;
        auto gapIt = upperBoundGapAnchor(alignedRead.gaps, gapPos - alignedRead.beginPos, seqan::SortGapPos());
        auto gapItEnd = end(alignedRead.gaps, seqan::Standard());
        // Note: We might create empty gaps here
        for (; gapIt != gapItEnd; ++gapIt)
            --(gapIt->gapPos);
        return true;
    }
    return false;
}

template <typename TAlignedReads, typename TGapPos>
int DifferingReadsRemoverImpl::removeGap(TAlignedReads & alignedReadStore, TGapPos const gapPos) const
{
    int result = 0;

    auto alignIt = begin(alignedReadStore, seqan::Standard());
    auto alignItEnd = end(alignedReadStore, seqan::Standard());
    for (; alignIt != alignItEnd; ++alignIt)
        result += removeGap2(*alignIt, gapPos);
    return result;
}

void DifferingReadsRemoverImpl::fixContigs(TFragmentStore /*const*/ & store) const
{
    seqan::String<TProfileChar> profile;  // the current contig's profile

    for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
    {
        auto itBegin = lowerBoundAlignedReads(store.alignedReadStore, contigID, seqan::SortContigId());
        auto itEnd = upperBoundAlignedReads(store.alignedReadStore, contigID, seqan::SortContigId());
        if (itBegin == itEnd)
            continue;  // no such contig

        // Get profile length.
        typedef decltype(*itBegin) TAlignedRead;
        auto itMax = std::max_element(itBegin, itEnd, [](TAlignedRead lhs, TAlignedRead rhs) {
                return std::max(lhs.beginPos, lhs.endPos) < std::max(rhs.beginPos, rhs.endPos);
            });
        unsigned profileLen = std::max(itMax->beginPos, itMax->endPos);
        clear(profile);

        // Build the profile.
        resize(profile, profileLen, TProfileChar());
        for (auto it = itBegin; it != itEnd; ++it)
        {
            typedef std::remove_reference<decltype(store.readSeqStore[0])>::type TReadSeq;
            typedef typename std::remove_reference<decltype(it->gaps)>::type TReadGapAnchors;
            typedef seqan::Gaps<TReadSeq, seqan::AnchorGaps<TReadGapAnchors> > TReadGaps;

            TReadGaps readGaps(store.readSeqStore[it->readId], it->gaps);
            auto pos = std::min(it->beginPos, it->endPos);
            for (auto it = begin(readGaps, seqan::Standard()); it != end(readGaps, seqan::Standard()); ++it, ++pos)
                if (isGap(it))
                    profile[pos].count[5] += 1;
                else
                    profile[pos].count[ordValue(*it)] += 1;
        }

        // Build the sequence if consensus says it's not a gap character.
        clear(store.contigStore[contigID].seq);
        unsigned pos = 0;
        std::vector<unsigned> gapPositions;
        for (auto const & c : profile)
        {
            if (_getMaxIndex(c) != 5)
            {
                appendValue(store.contigStore[contigID].seq, seqan::Dna5(_getMaxIndex(c)));
                ++pos;
            }
            else
            {
                gapPositions.push_back(pos);
            }
        }

        // Insert the gaps.
        clear(store.contigStore[contigID].gaps);
        typedef decltype(store.contigStore[0].seq) TContigSeq;
        typedef typename std::remove_reference<decltype(store.contigStore[0].gaps)>::type TContigGapAnchors;
        typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<TContigGapAnchors> > TContigGaps;
        TContigGaps contigGaps(store.contigStore[contigID].seq, store.contigStore[contigID].gaps);
        for (auto it = gapPositions.rbegin(); it != gapPositions.rend(); ++it)
            insertGap(contigGaps, *it);
    }
}

// ----------------------------------------------------------------------------
// Class DifferingReadsRemover
// ----------------------------------------------------------------------------

DifferingReadsRemover::DifferingReadsRemover(double maxErrorRate) :
        impl(new DifferingReadsRemoverImpl(maxErrorRate))
{}

DifferingReadsRemover::~DifferingReadsRemover()
{}

void DifferingReadsRemover::run(TFragmentStore & store, std::vector<unsigned> & readIDs) const
{
    impl->run(store, readIDs);
}

// ----------------------------------------------------------------------------
// Class ContigSplitterImpl
// ----------------------------------------------------------------------------

// Split contigs of fragment store at zero-coverage positions and at positions where a split would improve the pairwise
// alignment score (given a scheme with match=1, mismatch=-2, gaps=-2).

class ContigSplitterImpl
{
public:
    ContigSplitterImpl(TFragmentStore & store, bool debug = false) : debug(debug), store(store)
    {}

    // Perform the separation
    void run();

private:

    // Move alignments to new contig.
    void moveToNewContig(TAlignedReadIter itBegin, TAlignedReadIter itEnd, unsigned contigID);

    // Print a fragment store to a std::ostream.
    void printStore(std::ostream & out, TFragmentStore & store) const;

    // Prints a profile.
    template <typename TProfile>
    void printProfile(std::ostream & out, char const * label, TProfile profile);

    // Splite the given contig at the given positions.
    std::vector<TAlignedReadIter> collectSplitPoints(TAlignedReadIter aliItBegin, TAlignedReadIter aliItEnd);

    // Return SOP score of profile.
    template <typename TProfileSeq>
    int profileSopScore(TProfileSeq const & profile,
                        seqan::Score<int, seqan::Simple> const & scoringScheme) const;

    // Split the profile.
    void splitProfile(seqan::String<TProfileChar> & leftProf,
                      seqan::String<TProfileChar> & rightProf,
                      typename seqan::Infix<seqan::String<TProfileChar>>::Type source,
                      TAlignedReadIter itBegin,
                      TAlignedReadIter itEnd);

    // Whether or not to enable debugging.
    bool debug;
    // The FragmentStore to work on.
    TFragmentStore & store;
};

void ContigSplitterImpl::run()
{
    if (debug)
    {
        std::cerr << "BEFORE SPLITTING CONTIGS\n";
        printStore(std::cerr, store);
    }
    // Sort aligned reads by coordinate.
    sortAlignedReads(store.alignedReadStore, seqan::SortBeginPos());
    sortAlignedReads(store.alignedReadStore, seqan::SortContigId());

    // Collect split points for each contig.
    std::vector<std::vector<TAlignedReadIter>> splitPoints;
    for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
    {
        auto itBegin = lowerBoundAlignedReads(store.alignedReadStore, contigID, seqan::SortContigId());
        auto itEnd = upperBoundAlignedReads(store.alignedReadStore, contigID, seqan::SortContigId());
        splitPoints.emplace_back(collectSplitPoints(itBegin, itEnd));
        auto it = std::unique(splitPoints.back().begin(), splitPoints.back().end());
        splitPoints.back().resize(it - splitPoints.back().begin());
    }

    // Iterate over the split points for each contig and split contig.
    std::vector<unsigned> contigIDs;
    unsigned nextContigID = length(store.contigStore);  // next contig id to use
    for (unsigned contigID = 0; contigID < length(splitPoints); ++contigID)
    {
        if (splitPoints[contigID].size() == 2u)
            continue;  // no work, ignore
        contigIDs.push_back(contigID);  // realign current
        for (unsigned i = 0; i + 1 < length(splitPoints[contigID]); ++i)
            moveToNewContig(splitPoints[contigID][i], splitPoints[contigID][i + 1], nextContigID++);
    }
    resize(store.contigStore, nextContigID);  // resize store for realigner
    for (unsigned contigID = contigIDs.size(); contigID < nextContigID; ++contigID)
    {
        std::stringstream ss;
        ss << "contig_" << contigID;
        appendValue(store.contigNameStore, ss.str());
        contigIDs.push_back(contigID);
    }

    if (debug)
    {
        std::cerr << "BEFORE REALIGNING CONTIGS\n";
        printStore(std::cerr, store);
    }

    realign(store, contigIDs, 40, 30, false);  // TODO(holtgrew): Hand through ANISE options for this.

    if (debug)
    {
        std::cerr << "AFTER SPLITTING CONTIGS\n";
        printStore(std::cerr, store);
    }
}

// Move alignments to new contig.
void ContigSplitterImpl::moveToNewContig(TAlignedReadIter itBegin, TAlignedReadIter itEnd, unsigned contigID)
{
    if (itBegin == itEnd)
        return;

    int shift = std::min(itBegin->beginPos, itBegin->endPos);

    for (auto it = itBegin; it != itEnd; ++it)
    {
        if (debug)
            std::cerr << "MOVING " << it->contigId << ", " << it->beginPos << ", " << it->endPos << ", " << store.readSeqStore[it->readId] << " TO " << contigID << " WITH SHIFT " << shift << "\n";
        it->contigId = contigID;
        it->beginPos -= shift;
        it->endPos -= shift;
    }
}

// Print a fragment store to a std::ostream.
void ContigSplitterImpl::printStore(std::ostream & out, TFragmentStore & store) const
{
    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    for (unsigned i = 0; i < length(layout.contigRows); ++i)
    {
        // Get coordinates to plot.
        __int64 l = 0;
        __int64 r = l;
        for (unsigned j = 0; j < length(layout.contigRows[i]); ++j)
        {
            unsigned id = back(layout.contigRows[i][j]);
            if (r < store.alignedReadStore[id].beginPos)
                r = store.alignedReadStore[id].beginPos;
            if (r < store.alignedReadStore[id].endPos)
                r = store.alignedReadStore[id].endPos;
        }

        out << ">multi-read-alignment: contig_" << i << "\n";
        printAlignment(out, seqan::Raw(), layout, store, i, l, r, 0, 1000);
    }
}

// Prints a profile.
template <typename TProfile>
void ContigSplitterImpl::printProfile(std::ostream & out, char const * label, TProfile profile)
{
    out << ">" << label << "\n";

    out << "    ";
    for (unsigned i = 0; i < length(profile); ++i)
        if (i % 10 == 0)
            out << ": ";
        else if (i % 5 == 0)
            out << ". ";
        else
            out << "  ";
    out << "\n";

    char labels[6] = {'A', 'C', 'G', 'T', 'N', '-'};
    for (unsigned i = 0; i < 6; ++i)
    {
        out << labels[i] << " ";
        for (unsigned j = 0; j < length(profile); ++j)
        {
            int count = profile[j].count[i];
            if (count < 10)
                out << "  " << count;
            else
                out << " " << count;
        }
        out << "\n";
    }
}

// Return SOP score of profile.
template <typename TProfileSeq>
int ContigSplitterImpl::profileSopScore(TProfileSeq const & profile,
                                        seqan::Score<int, seqan::Simple> const & scoringScheme) const
{
    int result = 0;
    for (TProfileChar c : profile)
    {
        // Count matches.
        for (unsigned i = 0; i < 5; ++i)
            if (c.count[i] > 1)
                result += scoreMatch(scoringScheme) * (c.count[i] - 1) * (c.count[i] - 1);
        // Count mismatches.
        for (unsigned i = 0; i < 5; ++i)
            for (unsigned j = i + 1; j < 5; ++j)
                result += scoreMismatch(scoringScheme) * c.count[i] * c.count[j];
        // Count gaps.
        for (unsigned i = 0; i < 5; ++i)
            result += scoreGapOpen(scoringScheme) * c.count[i] * c.count[5/* => gap*/];
    }
    return result;
}

// Split the profile.
void ContigSplitterImpl::splitProfile(seqan::String<TProfileChar> & leftProf,
                                      seqan::String<TProfileChar> & rightProf,
                                      typename seqan::Infix<seqan::String<TProfileChar>>::Type source,
                                      TAlignedReadIter itBegin,
                                      TAlignedReadIter itEnd)
{
    // Get begin and end position of MSA/profile to extract.
    int beginPos = beginPosition(source);
    int endPos = endPosition(source);

    // Build right profile.
    resize(rightProf, length(source), TProfileChar());

    for (auto it = itBegin; it != itEnd; ++it)
    {
        typedef std::remove_reference<decltype(store.readSeqStore[0])>::type TReadSeq;
        typedef typename std::remove_reference<decltype(it->gaps)>::type TReadGapAnchors;
        typedef seqan::Gaps<TReadSeq, seqan::AnchorGaps<TReadGapAnchors> > TReadGaps;

        TReadGaps readGaps(store.readSeqStore[it->readId], it->gaps);
        auto pos = std::min(it->beginPos, it->endPos);
        SEQAN_ASSERT_GEQ(pos, beginPos);
        for (auto itR = begin(readGaps, seqan::Standard()); itR != end(readGaps, seqan::Standard()) && pos < endPos; ++itR, ++pos)
            if (isGap(itR))
                rightProf[pos - beginPos].count[5] += 1;
            else
                rightProf[pos - beginPos].count[ordValue(*itR)] += 1;
    }

    // leftProf is source minus counts from rightProf.
    leftProf = source;
    for (unsigned i = 0; i < length(source); ++i)
        for (unsigned j = 0; j < 6; ++j)
            leftProf[i].count[j] = source[i].count[j] - rightProf[i].count[j];
}

// TODO(holtgrew): This function needs some love.
std::vector<TAlignedReadIter> ContigSplitterImpl::collectSplitPoints(
    TAlignedReadIter aliItBegin, TAlignedReadIter aliItEnd)
{
    if (aliItBegin == aliItEnd)
        return { aliItBegin, aliItEnd };  // no work
    std::vector<TAlignedReadIter> result;
    result.push_back(aliItBegin);

    // Alignment iterator, indexed by read id.
    std::map<unsigned, TAlignedReadIter> aliMap;
    for (auto it = aliItBegin; it != aliItEnd; ++it)
        aliMap[it->readId] = it;

    // Build profile for these alignments.
    //
    // Get alignment length and initialize profile.
    typedef decltype(*aliItBegin) TAlignedRead;
    auto getEndPos = [](TAlignedRead const & el) { return std::max(el.beginPos, el.endPos); };
    auto itMax = std::max_element(aliItBegin, aliItEnd, [getEndPos](TAlignedRead lhs, TAlignedRead rhs) {
                                      return getEndPos(lhs) < getEndPos(rhs);
                                  });
    int maxPos = getEndPos(*itMax);
    seqan::String<TProfileChar> profile;
    resize(profile, maxPos, TProfileChar());
    // Create the actual profile.
    for (auto it = aliItBegin; it != aliItEnd; ++it)
    {
        typedef std::remove_reference<decltype(store.readSeqStore[0])>::type TReadSeq;
        typedef typename std::remove_reference<decltype(it->gaps)>::type TReadGapAnchors;
        typedef seqan::Gaps<TReadSeq, seqan::AnchorGaps<TReadGapAnchors> > TReadGaps;

        TReadGaps readGaps(store.readSeqStore[it->readId], it->gaps);
        auto pos = std::min(it->beginPos, it->endPos);
        for (auto it = begin(readGaps, seqan::Standard()); it != end(readGaps, seqan::Standard()); ++it, ++pos)
            if (isGap(it))
                profile[pos].count[5] += 1;
            else
                profile[pos].count[ordValue(*it)] += 1;
    }
    // Print complete profile.
    if (debug)
        printProfile(std::cerr, "complete profile", profile);

    // TODO(holtgrew): Add sliding window over aligned pairs to check that we do not cut through too many links.

    // We keep an iterator into the aligned reads and keep it at >= of the current sliding position.
    TAlignedReadIter it = aliItBegin;

    // Create intervals from alignments.
    std::vector<SlidingWindowInterval> intervals;
    std::transform(aliItBegin, aliItEnd, std::back_inserter(intervals),
                   [](TAlignedRead el) {
                        return SlidingWindowInterval(el.contigId, std::min(el.beginPos, el.endPos),
                                                     std::max(el.beginPos, el.endPos), el.readId);
                   });
    // Create sliding windows algorithm.
    SlidingWindowAlgo slidingWindows(intervals);

    for (; !slidingWindows.atEnd(); slidingWindows.advance())
    {
        if (std::get<2>(slidingWindows.location()))
            continue;  // is before end, we only look after end/on begin
        // Get current position.
        int rID = std::get<0>(slidingWindows.location());
        int pos = std::get<1>(slidingWindows.location());

        // Advance it such that it is not left of the current position (should be on current).
        auto getPos = [](TAlignedReadIter it) {
            return std::pair<int, int>(it->contigId, std::min(it->beginPos, it->endPos));
        };
        while (it != aliItEnd && getPos(it) < std::make_pair(rID, pos))
            ++it;
        if (it == aliItEnd)
            continue;  // behind end position in sliding windows, no alignment here

        // Check if coverage is 0 and add location (need to subtract alignments starting here).
        unsigned x = 0;
        for (auto it = slidingWindows.activeBegin(); it != slidingWindows.activeEnd(); ++it)
        {
            std::cerr << "..aliMap[" << slidingWindows.getItemID(*it) << "]->beginPos == "
                      << aliMap[slidingWindows.getItemID(*it)]->beginPos << ", pos == " << pos << "\n";
            x += (aliMap[slidingWindows.getItemID(*it)]->endPos == pos);
        }
        if (debug)
            std::cerr << "rID=" << rID << ", pos=" << pos << ", slidingWindows.activeCount() == "
                      << slidingWindows.activeCount() << ", x == " << x << "\n";
        if (slidingWindows.activeCount() == x)
        {
            if (debug)
                std::cerr << "A SPLIT POINT IS " << rID << ", " << pos << "\n";
            result.push_back(it);
            continue;
        }

        // Check SOP change if we split at (rID, pos).
        //
        // First, find out the largest extent into the profile by currently active reads.
        int endPos = pos;
        for (auto it = slidingWindows.activeBegin(); it != slidingWindows.activeEnd(); ++it)
            endPos = std::max(endPos, (int)std::max(aliMap[*it]->beginPos, aliMap[*it]->endPos));
        seqan::Score<int, seqan::Simple> scoringScheme(1, -2, -2);
        // Compute SOP score of the profile infix [pos, endPos).
        typename seqan::Infix<seqan::String<TProfileChar> >::Type inf(profile, pos, endPos);
        int scoreBefore = profileSopScore(inf, scoringScheme);

        // Compute end of overlap, i.e. iterator right of it that does not overlap with endPos.
        TAlignedReadIter itOvlEnd = it;
        while (itOvlEnd != aliItEnd && getPos(itOvlEnd) < std::make_pair(rID, endPos))
            ++itOvlEnd;

        // Split the profile infix; alignments in [it, itOvlEnd) will be subtracted to leftProf and added to rightProf.
        seqan::String<TProfileChar> leftProf, rightProf;
        splitProfile(leftProf, rightProf, inf, it, itOvlEnd);
        int scoreAfter = profileSopScore(leftProf, scoringScheme) + profileSopScore(rightProf, scoringScheme);

        if (debug)
        {
            std::cerr << "rID == " << rID << "\n";
            std::cerr << "pos == " << pos << "\n";
            std::cerr << "original score = " << profileSopScore(inf, scoringScheme) << "\n";
            std::cerr << "left score     = " << profileSopScore(leftProf, scoringScheme) << "\n";
            std::cerr << "right score    = " << profileSopScore(rightProf, scoringScheme) << "\n";
            printProfile(std::cerr, "original profile", inf);
            printProfile(std::cerr, "left profile", leftProf);
            printProfile(std::cerr, "right profile", rightProf);
        }

        // If the SOP score improved then we add (rID, pos) as a split position.
        if (scoreAfter > scoreBefore)
        {
            result.push_back(it);
            if (debug)
            {
                std::cerr << "A SPLIT POINT IS " << rID << ", " << pos << "\n"
                          << "  score before: " << scoreBefore << ", score after: " << scoreAfter << "\n";
                profileSopScore(inf, scoringScheme);
                printProfile(std::cerr, "original profile", inf);
                printProfile(std::cerr, "left profile", leftProf);
                printProfile(std::cerr, "right profile", rightProf);
            }

            // Advance behind the last/rightmost current overlap.
            slidingWindows.advanceTo(rID, endPos, true);
        }
    }

    result.push_back(aliItEnd);
    return result;
}

// ----------------------------------------------------------------------------
// Class ContigSplitter
// ----------------------------------------------------------------------------

ContigSplitter::ContigSplitter(TFragmentStore & store, bool debug) :
        impl(new ContigSplitterImpl(store, debug))
{}

ContigSplitter::~ContigSplitter()
{}

void ContigSplitter::run()
{
    impl->run();
}

// ----------------------------------------------------------------------------
// Function countErrors()
// ----------------------------------------------------------------------------

// Number of errors from read alignment vs. contig.

int countErrors(TFragmentStore /*const*/ & store, TAlignedReadStoreElement /*const*/ & el)
{
    typedef typename TFragmentStore::TContigSeq TContigSeq;
    typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
    typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<seqan::String<TContigGapAnchor> > > TContigGaps;
    TContigGaps contigGaps(store.contigStore[el.contigId].seq, store.contigStore[el.contigId].gaps);

    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename TFragmentStore::TReadGapAnchor TReadGapAnchor;
    typedef seqan::Gaps<TReadSeq, seqan::AnchorGaps<seqan::String<TReadGapAnchor> > > TReadGaps;
    SEQAN_CHECK(el.beginPos <= el.endPos, "Everything should be forward here.");
    TReadGaps readGaps(store.readSeqStore[el.readId], el.gaps);

    SEQAN_CHECK((int)length(readGaps) == el.endPos - el.beginPos, "Must match.");

    auto itC = iter(contigGaps, el.beginPos, seqan::Standard());
    auto itCEnd = iter(contigGaps, el.endPos, seqan::Standard());
    auto itR = begin(readGaps, seqan::Standard());
    auto itREnd = end(readGaps, seqan::Standard());

    int result = 0;
    for (; itC != itCEnd && itR != itREnd; ++itC, ++itR)
    {
        if (isGap(itC) != isGap(itR))
            result += 1;
        else if (!isGap(itC) && !isGap(itR) && (*itC != *itR))
            result += 1;
    }

    return result;
}

// ----------------------------------------------------------------------------
// Function printStore()
// ----------------------------------------------------------------------------

void printStore(std::ostream & out, TFragmentStore const & storeC)
{
    TFragmentStore store(storeC);

    seqan::AlignedReadLayout layout;
    layoutAlignment(layout, store);
    for (unsigned i = 0; i < length(layout.contigRows); ++i)
    {
        // Get coordinates to plot.
        __int64 l = 0;
        __int64 r = l;
        for (unsigned j = 0; j < length(layout.contigRows[i]); ++j)
        {
            unsigned id = back(layout.contigRows[i][j]);
            if (r < store.alignedReadStore[id].beginPos)
                r = store.alignedReadStore[id].beginPos;
            if (r < store.alignedReadStore[id].endPos)
                r = store.alignedReadStore[id].endPos;
        }

        out << ">multi-read-alignment: contig_" << i << "\n";
        printAlignment(out, seqan::Raw(), layout, store, i, l, r, 0, 1000);
    }
}

