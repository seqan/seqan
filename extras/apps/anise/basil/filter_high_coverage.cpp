// ==========================================================================
//                                 BASIL
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

#include "filter_high_coverage.h"

#include <queue>

#include <seqan/seq_io.h>

namespace  // anonymous
{

// ----------------------------------------------------------------------------
// Class CoverageAlgorithm
// ----------------------------------------------------------------------------

class CoverageAlgorithm
{
public:
    CoverageAlgorithm(seqan::String<seqan::GenomicRegion> & result, int maxCoverage)
            : result(result), maxCoverage(maxCoverage), rID(-1), coverage(0), beginPos(-1), pos(-1),
              doneUpTo(-1)
    {}

    // Register the next BamAlignmentRecord with algorithm.
    void push(seqan::BamAlignmentRecord const & record);

    // Finish running the algorithm.
    void finish();

private:
    // The resulting regions.
    seqan::String<seqan::GenomicRegion> & result;

    // Maximum coverage to allow.
    int maxCoverage;

    // The current reference id.
    int rID;
    // The current coverage.
    int coverage;
    // The first position where maxCoverage was reached.
    int beginPos;
    // The current position.
    int pos;
    // The position up towards we guarantee there either is an entry in results or there is no too high coverage.
    int doneUpTo;

    // A list with (pos, isBegin) events at each begin and end position of a read alignment.  The interpretation of
    // integers gives an ordering where end events come before begin events.
    std::deque<std::pair<int, bool> > events;

    // Go through events until we reach behind the last end event at pos.
    void advanceTo(int toPos);

    // Push the current region in (rID, beginPos, pos) to result.
    void pushInterval();

    // TODO(holtgrew): Clean up HighCoverageFilterImpl so this is not necessary any more.
    friend class ::HighCoverageFilterImpl;
};

// ----------------------------------------------------------------------------
// Class CoverageAlgorithm
// ----------------------------------------------------------------------------

void CoverageAlgorithm::push(seqan::BamAlignmentRecord const & record)
{
    // Handle transition between contigs.
    if (record.rID != rID)
    {
        finish();  // finish contig
        rID = record.rID;
        SEQAN_ASSERT_EQ(coverage, 0);
        beginPos = -1;
        pos = -1;
    }

    // Break if record is an orphan.
    if (record.rID == seqan::BamAlignmentRecord::INVALID_REFID)
        return;

    // Otherwise, insert the begin and end events and advance to the current begin position.
    std::pair<int, bool> beginEvent(record.beginPos, true);
    std::pair<int, bool> endEvent(record.beginPos + getAlignmentLengthInRef(record), false);
    events.insert(std::upper_bound(events.begin(), events.end(), beginEvent), beginEvent);
    events.insert(std::upper_bound(events.begin(), events.end(), endEvent), endEvent);
    advanceTo(record.beginPos);
}

void CoverageAlgorithm::finish()
{
    if (!events.empty())
        advanceTo(events.back().first);
}

void CoverageAlgorithm::advanceTo(int toPos)
{
    std::pair<int, bool> toPair(toPos, true);
    while (!events.empty() && toPair > events.front())
    {
        std::pair<int, bool> const & event = events.front();
        pos = event.first;
        if (event.second)  // OPEN
        {
            coverage += 1;
            if (coverage == maxCoverage)
                beginPos = pos;
        }
        else  // CLOSE
        {
            if (coverage == maxCoverage)
                pushInterval();
            coverage -= 1;
        }
        events.pop_front();

        // If we are below the coverage for the current position then we can store it as doneUpTo.
        if (coverage < maxCoverage)
            doneUpTo = pos;
    }
}

void CoverageAlgorithm::pushInterval()
{
    seqan::GenomicRegion region;
    region.seqId = rID;
    region.beginPos = beginPos;
    region.endPos = pos;
    appendValue(result, region);
    // std::cerr << "PUSHING REGION beginPos = " << beginPos << ", endPos = " << pos << "\n";
    doneUpTo = pos;  // Can process all records up to doneUpTo.
}

}  // namespace anonymous

// ----------------------------------------------------------------------------
// Class HighCoverageFilterImpl
// ----------------------------------------------------------------------------

class HighCoverageFilterImpl
{
public:
    HighCoverageFilterImpl(int maxCoverage) :
            coverageAlgorithm(highCoverageRegions, maxCoverage),
            currentCoverageRegion(0), maxBufferSize(0)
    {}

    void filter(std::vector<seqan::BamAlignmentRecord *> & out,
                std::vector<seqan::BamAlignmentRecord *> const & in);

    void finish(std::vector<seqan::BamAlignmentRecord *> & out);

private:
    // The algorithm object to use for computing high-coverage region and a result storage.
    CoverageAlgorithm coverageAlgorithm;
    seqan::String<seqan::GenomicRegion> highCoverageRegions;
    // Index of high coverage region that has not completely passed yet.
    int currentCoverageRegion;

    // Buffer for the alignment records up to the current position.
    std::queue<seqan::BamAlignmentRecord *> buffer;

    void advance(std::vector<seqan::BamAlignmentRecord *> & out);

    unsigned maxBufferSize;
};

void HighCoverageFilterImpl::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                                    std::vector<seqan::BamAlignmentRecord *> const & in)
{
    // Push everything from the input buffer into the internal buffer and the coverage algorithm.
    for (auto ptr : in)
    {
        buffer.push(ptr);
        coverageAlgorithm.push(*ptr);
    }

    maxBufferSize = std::max(maxBufferSize, (unsigned)buffer.size());

    // Advance in buffers.
    advance(out);
}

void HighCoverageFilterImpl::finish(std::vector<seqan::BamAlignmentRecord *> & out)
{
    // Set into state that will make advance() write out everything.
    coverageAlgorithm.finish();
    if (!buffer.empty())
        coverageAlgorithm.rID = buffer.front()->rID;
    coverageAlgorithm.doneUpTo = seqan::maxValue<int>();

    // Advance in buffers.
    advance(out);

    SEQAN_CHECK(buffer.empty(), "Buffer must be empty at end. Size was: %d.",
                (int)buffer.size());
    //fprintf(stderr, "HIGH COVERAGE MAX BUFFER SIZE\t%d\n", maxBufferSize);
}

void HighCoverageFilterImpl::advance(std::vector<seqan::BamAlignmentRecord *> & out)
{
    // We can now write out records from the buffer.
    while (!buffer.empty())
    {
        // Get shortcut to first record in buffer and its positions.
        seqan::BamAlignmentRecord const * record = buffer.front();
        int beginPos = record->beginPos;
        int endPos = beginPos + getAlignmentLengthInRef(*record);

        // Check whether we (1) drop the record, (2) write it to out, (3) advance the current region, (4) break because
        // we are behind coverageAlgorithm.doneUpTo.

        // Break if we would go further than the coverage algorithms' doneUpTo.
        if (std::make_pair(record->rID, endPos) >= std::make_pair(coverageAlgorithm.rID, coverageAlgorithm.doneUpTo))
            break;

        // Advance current high coverage region until record overlaps with the current one or the record is in front of
        // the coverage region.
        while (currentCoverageRegion < (int)length(highCoverageRegions) &&
               std::make_pair(record->rID, record->beginPos) >= std::make_pair(highCoverageRegions[currentCoverageRegion].seqId,
                                                                               highCoverageRegions[currentCoverageRegion].endPos))
            ++currentCoverageRegion;

        // Discard record if it overlaps with current region (if any).
        if (currentCoverageRegion < (int)length(highCoverageRegions))
        {
            seqan::GenomicRegion const & region = highCoverageRegions[currentCoverageRegion];
            if (beginPos < region.endPos && region.beginPos < endPos)
            {
                delete buffer.front();
                buffer.pop();
                continue;
            }
        }

        // Otherwise, move record from buffer to out buffer.
        out.push_back(buffer.front());
        buffer.pop();
    }
}

// ----------------------------------------------------------------------------
// Class HighCoverageFilter
// ----------------------------------------------------------------------------

HighCoverageFilter::HighCoverageFilter(int maxCoverage) : impl(new HighCoverageFilterImpl(maxCoverage))
{}

HighCoverageFilter::~HighCoverageFilter()
{}

void HighCoverageFilter::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                                std::vector<seqan::BamAlignmentRecord *> const & in)
{
    impl->filter(out, in);
}

void HighCoverageFilter::finish(std::vector<seqan::BamAlignmentRecord *> & out)
{
    impl->finish(out);
}
