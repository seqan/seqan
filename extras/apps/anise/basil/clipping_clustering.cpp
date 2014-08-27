// ==========================================================================
//                                 BASIL
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

#include "clipping_clustering.h"

#include "utils.h"

namespace  // anonymous
{

// ----------------------------------------------------------------------------
// Enum ClippingSide
// ----------------------------------------------------------------------------

// Identifies the side of the clipping.
enum class ClippingSide { INVALID, LEFT, RIGHT };

// ----------------------------------------------------------------------------
// Function clippingSide
// ----------------------------------------------------------------------------

ClippingSide clippingSide(seqan::BamAlignmentRecord const & record)
{
    return startsWithClipping(record.cigar) ? ClippingSide::RIGHT : ClippingSide::LEFT;
}

// ----------------------------------------------------------------------------
// Class FilterMinCoverage
// ----------------------------------------------------------------------------

// Input is a stream of strictly sorted positions.  Each position/item is assigned an id.  The output is a stream of
// pairs (item id, coverage) of items that fall into regions with a min coverage (in windows of a given size).  This can
// be used to remove spurious hits and determine the coverage by points in a streaming fashion.

template <typename TCargo>
class FilterMinCoverage
{
public:

    typedef std::pair<int, int> TPosition;

    // Constructor, default constructed filter just passes through the ids.
    FilterMinCoverage(int windowSpan = 100, int minCoverage = 1, bool debug = false) :
            _windowSpan(windowSpan), _minCoverage(minCoverage),
            _lastYieldedItem(std::make_pair(-1, 0), TCargo()),
            debug(debug),
            numIn(0),
            numOut(0)
    {}

    // Push new position, returning new id for this position.
    //
    // Write out the ids of all positions whose coverage cannot be influenced by positions pushed in the future.
    template <typename TOutIt>
    void push(TPosition const & pos, TCargo const & cargo, TOutIt it)
    {
        numIn += 1;
        TItem item(pos, cargo);

        // Fast forward through items, stop if current position empty (will be set to pos below) or it is within
        // _windowSpan of pos.
        while (!_queueCurrent.empty() && (unsigned)_distance(pos, _queueCurrent.front().first) > (unsigned)_windowSpan)
            _advanceByOnePosition(it);

        // If pos was a sentinel value then stop here.  We will have written out everything above.
        if (pos.first == -1)
            return;

        // Insert item with new pos into appropriate queue.
        if (_queueCurrent.empty() || (_queueRight.empty() && _queueCurrent.front().first == pos))
        {
            if (debug)
                printf("pushing to current: (%u, %u)\n", (unsigned)item.first.first, (unsigned)_queueRight.front().first.second);
            _queueCurrent.push_back(item);
        }
        else
        {
            if (debug)
                printf("pushing to right: (%u, %u)\n", (unsigned)item.first.first, (unsigned)_queueRight.front().first.second);
            _queueRight.push_back(item);
        }
    }

    // Trigger writing out the remaining ids.
    template <typename TOutIt>
    void finish(TOutIt it)
    {
        push(std::make_pair(-1, -1), TCargo(), it);
    }

private:

    // Advance queue to next position from right queue, writing out all positions from items that go from the current
    // position to the left queue.
    template <typename TOutIt>
    void _advanceByOnePosition(TOutIt outIt)
    {
        // Write out coverage.
        int coverage = _queueLeft.size() + _queueCurrent.size() + _queueRight.size();
        for (typename std::deque<TItem>::iterator it = _queueCurrent.begin(); it != _queueCurrent.end(); ++it)
            if (coverage >= _minCoverage)
            {
                numOut += 1;
                if (debug)
                    std::cout << "Filter min coverage yields " << it->first.first << ", " << it->first.second << ", " << coverage << "\n";
                *outIt++ = std::make_pair(*it, coverage);
            }
        if (!_queueCurrent.empty())
            _lastYieldedItem = _queueCurrent.back();
        // Move items from current queue to left queue.
        while (!_queueCurrent.empty())
        {
            if (debug)
                printf("from current to left: (%u, %u)\n", (unsigned)_queueCurrent.front().first.first, (unsigned)_queueCurrent.front().first.second);
            _queueLeft.push_back(_queueCurrent.front());
            _queueCurrent.pop_front();
        }
        // Move items from right queue into current queue.
        while ((_queueCurrent.empty() && !_queueRight.empty()) ||
               (!_queueRight.empty() && _queueCurrent.back().first == _queueRight.front().first))
        {
            if (debug)
                printf("from right to current: (%u, %u)\n", (unsigned)_queueRight.front().first.first, (unsigned)_queueRight.front().first.second);
            _queueCurrent.push_back(_queueRight.front());
            _queueRight.pop_front();
        }
        // Handle case of queue running empty.
        if (_queueCurrent.empty())
        {
            SEQAN_ASSERT(_queueRight.empty());
            if (debug)
                printf("clearing left queue\n");
            _queueLeft.clear();
            return;
        }
        // Let items fall out of left queue that are too far away from current position.
        //
        // We use conversion to unsigned to get maxValue in case of distance -1.
        while (!_queueLeft.empty() &&
               (unsigned)_distance(_queueCurrent.front().first, _queueLeft.front().first) > (unsigned)_windowSpan)
        {
            if (debug)
                printf("popping from left queue: (%u, %u)\n", (unsigned)_queueLeft.front().first.first, (unsigned)_queueLeft.front().first.second);
            _queueLeft.pop_front();
        }
    }

    // Return distance of both positions.  Two positions that lie on different contigs have distance -1, indicating
    // infinity.
    int _distance(TPosition const & lhs, TPosition const & rhs)
    {
        if (lhs.first != rhs.first)
            return -1;
        return lhs.second - rhs.second;
    }

    typedef std::pair<TPosition, TCargo> TItem;

    // The length of the window to use towards each side to use.
    int _windowSpan;

    // The minimal coverage to require.
    int _minCoverage;

    // The last yielded item.  No item with a smaller id can be yielded.
    TItem _lastYieldedItem;

    // We manage three queues, one for the left side of the window, one for the current position, one for the right
    // side.  Items are yielded once they go from the current queue to the left queue.
    std::deque<TItem> _queueLeft;
    std::deque<TItem> _queueCurrent;
    std::deque<TItem> _queueRight;

    // Whether or not to enable debugging.
    bool debug;

    // Number of items that went into/out of the filter.
    unsigned numIn;
    unsigned numOut;
};

// ----------------------------------------------------------------------------
// Class FilterIntervalMerger
// ----------------------------------------------------------------------------

// Reads through a stream of intervals with increasing begin position, yields a stream of merged intervals.

class FilterIntervalMerger
{
public:
    typedef std::pair<int, int> TPosition;
    typedef std::pair<TPosition, TPosition> TInterval;

    FilterIntervalMerger(bool debug = false) : debug(debug)
    {}

    template <typename TOutIt>
    void push(TPosition const & beginPos, TPosition const & endPos, ClippingSide side, TOutIt outIt)
    {
        SEQAN_ASSERT_EQ(beginPos.first, endPos.first);  // must be on same contig

        // Write out any intervals that cannot be merged with any future intervals any more.
        TIntervalWithCounts newInterval(std::make_pair(beginPos, endPos),
                                        std::make_pair((side == ClippingSide::LEFT), (side == ClippingSide::RIGHT)));
        if (debug)
            std::cerr << "Pushed into merger: [(" << newInterval.first.first.first << ", "
                      << newInterval.first.first.second << "), ("
                      << newInterval.first.second.first << ", "
                      << newInterval.first.second.second << ")]\n";
        while (!_intervalQueue.empty() && _beginsAfter(newInterval.first, _intervalQueue.front().first))
        {
            if (debug)
                std::cerr << "Merger yields [(" << _intervalQueue.front().first.first.first << ", "
                          << _intervalQueue.front().first.first.second << "), ("
                          << _intervalQueue.front().first.second.first << ", "
                          << _intervalQueue.front().first.second.second << ")]\n";
            *outIt++ = _intervalQueue.front();
            _intervalQueue.pop_front();
        }

        // We are done if this was a sentinel value.
        if (beginPos.first == -1)
            return;

        // Merge the new interval into the sorted queue of intervals.
        //
        // We search for a range of intervals in _intervalQueue that overlap with newInterval.
        typedef std::deque<TIntervalWithCounts>::iterator TIter;
        bool foundAny = false;
        std::pair<TIter, TIter> overlapRange(_intervalQueue.end(), _intervalQueue.end());
        for (auto it = _intervalQueue.begin(); it != _intervalQueue.end(); ++it)
            if (_canJoin(it->first, newInterval.first))
            {
                if (!foundAny)
                    overlapRange.first = it;
                foundAny = true;
                overlapRange.second = it;
                ++overlapRange.second;
            }

        // Add interval if no overlap was found.  Can always be appended since the intervals are ordered by start
        // position.
        if (overlapRange.first == overlapRange.second)
        {
            if (debug)
                std::cerr << "Merger gets new interval [(" << newInterval.first.first.first << ", "
                          << newInterval.first.first.second << "), (" << newInterval.first.second.first << ", "
                          << newInterval.first.second.second << "]\n";
            _intervalQueue.push_back(newInterval);
            return;
        }

        // Merge with overlapping intervals.
        TIter it = --overlapRange.second;
        it->second.first += newInterval.second.first;
        it->second.second += newInterval.second.second;
        for (TIter it2 = overlapRange.first; it2 != it; ++it2)
        {
            it->second.first += it2->second.first;
            it->second.second += it2->second.second;
        }
        it->first.first = std::min(beginPos, overlapRange.first->first.first);
        it->first.second = std::max(endPos, it->first.second);
        _intervalQueue.erase(overlapRange.first, overlapRange.second);
        if (debug)
            std::cerr << "Merger merged into [(" << it->first.first.first << ", "
                      << it->first.first.second << "), (" << it->first.second.first << ", "
                      << it->first.second.second << " {" << it->second.first << ", " << it->second.second << "}]\n";
    }

    // Trigger writing out the remaining ids.
    template <typename TOutIt>
    void finish(TOutIt it)
    {
        push(std::make_pair(-1, -1), std::make_pair(-1, -1), ClippingSide::INVALID, it);
    }

private:

    typedef std::pair<TInterval, std::pair<unsigned, unsigned> > TIntervalWithCounts;  // [(interval, (countL, countR))]

    std::deque<TIntervalWithCounts> _intervalQueue;

    // Whether or not to enable debugging.
    bool debug;

    // Whether or not lhs begins after rhs.
    bool _beginsAfter(TInterval const & lhs, TInterval const & rhs) const
    {
        if (lhs.first.first == -1)
            return true;  // Stop marker.
        if (lhs.first.first > rhs.second.first)
            return true;  // Is on next contig.
        return (lhs.first.second > rhs.second.second);
    }

    // Whether or not the intervals can be joined (overlap or right one is directly next to the next).
    bool _canJoin(TInterval const & lhs, TInterval const & rhs) const
    {
        // All intervals should be on the same contig.
        SEQAN_ASSERT_EQ(lhs.first.first, lhs.second.first);
        SEQAN_ASSERT_EQ(rhs.first.first, rhs.second.first);
        SEQAN_ASSERT_EQ(lhs.first.first, rhs.second.first);

        // lhs = [a, b), rhs = [x, y), lhs and rhs overlap iff (b >= x) && (a <= y)
        return (lhs.second.second >= rhs.first.second) && (lhs.first.second <= rhs.second.second);
    }
};

// ----------------------------------------------------------------------------
// Class FilterReorder
// ----------------------------------------------------------------------------

// Input is a stream of positions (contig id, position on contig) that is roughly sorted (the distance between two
// positions on the same contig is limited and once a new contig is started, no position on an old contig may be added).
// When pushing, each position/item is assigned an id.  The output is a reordered stream of the position ids.

template <typename TCargo>
class FilterReorder
{
public:
    typedef std::pair<int, int> TPosition;
    typedef std::pair<TPosition, TCargo> TItem;

    // Constructor, default constructed filter just passes through the ids.
    FilterReorder(int maxSpan = 0, bool debug = false) : _maxSpan(maxSpan), debug(debug), numIn(0), numOut(0)
    {}

    // Push a position into the reordering filter, returning its id.  A position on contig -1 is interpreted as a
    // sentinel value for "at end" and the filter writes out all remaining ids.
    template <typename TOutIt>
    void push(TPosition const & pos, TCargo const & cargo, TOutIt it)
    {
        if (!_queue.empty() && pos.first != -1)
            SEQAN_ASSERT_MSG(pos >= _lastYieldedItem.first, "Cannot push a position left of the last yielded position.");
        numIn += 1;

        // Enqueue position as next item if it is not a sentinel value.
        if (pos.first != -1)
        {
            // Build item that is to be inserted and isert it into the queue such that the queue stays sorted.
            TItem item(pos, cargo);
            _queue.insert(std::upper_bound(_queue.begin(), _queue.end(), item), item);
        }

        // Write out ids of all enqueued item that are sufficiently far away from the queue tail's position.
        //
        // We use the trick of converting from int to unsigned to get maxValue<unsigned>() for distance -1;
        while (!_queue.empty() && ((pos.first == -1) || ((unsigned)_distance(_queue.back().first, _queue.front().first) > (unsigned)_maxSpan)))
        {
            _lastYieldedItem = _queue.front();
            if (debug)
                std::cout << "Yielding (" << _lastYieldedItem.first.first << ", " << _lastYieldedItem.first.second << ")\n";
            _queue.pop_front();
            numOut += 1;
            *it++ = _lastYieldedItem;
        }
    }

    // Trigger writing out the remaining ids.
    template <typename TOutIt>
    void finish(TOutIt it)
    {
        push(std::make_pair(-1, -1), TCargo(), it);
    }

private:

    // The largest deviation that two positions can have.  If a position of the next contig has been pushed then no more
    // position from the previous contig can be added.  The contig of id -1 is the largest one and marks the end.
    int _maxSpan;

    // The last yielded item.  No item with a smaller id can be yielded.
    TItem _lastYieldedItem;

    // The queue of unprocessed items.
    std::deque<TItem> _queue;

    // Whether or not to enable debugging.
    bool debug;

    // Number of items that went into the filter and out again.
    unsigned numIn;
    unsigned numOut;

    // Return distance of both positions.  Two positions that lie on different contigs have distance -1, indicating
    // infinity.
    int _distance(TPosition const & lhs, TPosition const & rhs)
    {
        if (lhs.first != rhs.first)
            return -1;
        return lhs.second - rhs.second;
    }
};

// ----------------------------------------------------------------------------
// Class ClippingScanner
// ----------------------------------------------------------------------------

// Scan through a stream of BamStreamRecord objects and generate a list of positions with sufficiently clustered
// clipping positions, intervals around them, and a list of the records that have alignments on them.

class ClippingScanner
{
public:
    ClippingScanner(ClippingClusterOptions const & options) :
            options(options), reorderFilter(options.maxAlignmentLength),
            minCoverageFilter(options.minCoverageWindowLength, 0), minCoverage(options.minCoverage),
            minCoverageWindowLength(options.minCoverageWindowLength), minClippingLength(options.minClippingLength)
    {}
    
    void push(std::vector<seqan::BamAlignmentRecord *> const & records);
    
    void finish();

private:

    typedef FilterReorder<ClippingSide>::TPosition TPosition;
    typedef FilterIntervalMerger::TInterval TInterval;

    ClippingClusterOptions options;

    // Reorder the records by the clipping position.
    FilterReorder<ClippingSide> reorderFilter;
    // A buffer of ids, sorted by position.
    std::vector<std::pair<TPosition, ClippingSide> > orderedClippingSites;

    // We do not filter out records in this filter but only use it to determine the coverage of each position.
    FilterMinCoverage<ClippingSide> minCoverageFilter;
    // Items yielded by minCoverageFilter with a smaller coverage are filtered out.
    int minCoverage;
    // A buffer of ((position, id), coverage), sorted by position and having a given minimal coverage within a window around the
    // position.
    typedef std::pair<TPosition, ClippingSide> TPosWithCargo;
    std::vector<std::pair<TPosWithCargo, unsigned> > coverageResults;

    // Generate intervals, given a stream of intervals with increasing begin positions.
    FilterIntervalMerger intervalMerger;
public:  // TODO(holtgrew): Create wrapper?
    // The collected intervals being generated around the sufficiently covered clipping positions.
    std::vector<std::pair<TInterval, std::pair<unsigned, unsigned> > > breakpointWindows;
private:
    // The window length to count clippings in.
    int minCoverageWindowLength;
    // The minimal length of considered clippings.
    int minClippingLength;
};

void ClippingScanner::push(std::vector<seqan::BamAlignmentRecord *> const & in)
{
    for (auto ptr : in)
    {
        seqan::BamAlignmentRecord const & record = *ptr;
        int const INVALID_REFID = seqan::BamAlignmentRecord::INVALID_REFID;
        if (record.rID != INVALID_REFID)  // We allow INVALID_REFID as a sentinel value.
        {
            if (hasFlagUnmapped(record))
                continue;  // Ignore unaligned record.
            // TODO(holtgrew): Ignore 1-2 clipped bases.
            if (startsWithClipping(record.cigar) == endsWithClipping(record.cigar))
                continue;  // Ignore record with no clipping or local alignment only.
            if (clippingLength(record) < minClippingLength)
                continue;  // Ignore too short clippings.
        }

        // Add clipping position to reordering filter.
        if (record.rID == INVALID_REFID)
            reorderFilter.finish(std::back_inserter(orderedClippingSites));
        else
            reorderFilter.push(std::make_pair(record.rID, clippingPosition(record)),
                               clippingSide(record),
                               std::back_inserter(orderedClippingSites));

        // Move all elements falling out of the reordering filter into the min coverage filter.
        for (auto const & site : orderedClippingSites)
            minCoverageFilter.push(site.first, site.second, std::back_inserter(coverageResults));

        // Stop minCoverageFilter in case of sentinel value.
        if (record.rID == INVALID_REFID)
            minCoverageFilter.finish(std::back_inserter(coverageResults));
        orderedClippingSites.clear();

        // Process coverage computation results:
        //
        // - Push intervals around sufficiently covered positions into interval merger.
        // - Remove records that are on insufficiently covered positions.
        for (auto const & siteWithCoverage : coverageResults)
        {
            std::pair<TPosition, ClippingSide> const & site = siteWithCoverage.first;
            unsigned coverage = siteWithCoverage.second;
            if (coverage >= (unsigned)minCoverage)
            {
                // Push intervals into interval merger.
                int rID = site.first.first;
                int pos = site.first.second;
                TPosition beginPos(rID, std::max(0, pos - minCoverageWindowLength));
                TPosition endPos(rID, pos + 1 + minCoverageWindowLength);
                intervalMerger.push(beginPos, endPos, site.second, std::back_inserter(breakpointWindows));
            }
        }
        coverageResults.clear();

        // Finish interval merging in case of sentinel.
        if (record.rID == INVALID_REFID)
            intervalMerger.finish(std::back_inserter(breakpointWindows));
    }
}

void ClippingScanner::finish()
{
    seqan::BamAlignmentRecord record;
    record.rID = seqan::BamAlignmentRecord::INVALID_REFID;
    std::vector<seqan::BamAlignmentRecord *> buffer;
    buffer.push_back(&record);
    push(buffer);
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class ClippingClusterAlgoImpl
// ----------------------------------------------------------------------------

class ClippingClusterAlgoImpl
{
public:
    ClippingClusterAlgoImpl(ClippingClusterOptions const & options) : clippingScanner(options)
    {}

    void push(std::vector<seqan::BamAlignmentRecord *> const & records);

    void finish();

    std::vector<ClippingClusterRecord> result();

private:
    ClippingScanner clippingScanner;

    std::vector<ClippingClusterRecord> _result;
};


void ClippingClusterAlgoImpl::push(std::vector<seqan::BamAlignmentRecord *> const & records)
{
    clippingScanner.push(records);
}

void ClippingClusterAlgoImpl::finish()
{
    clippingScanner.finish();

    // Convert breakpoint windows to clipping records.
    for (auto & window : clippingScanner.breakpointWindows)
    {
        int rID = window.first.first.first;
        int beginPos = window.first.first.second;
        int endPos = window.first.second.second;
        int leftWeight = window.second.first;
        int rightWeight = window.second.second;
        _result.push_back(ClippingClusterRecord(rID, beginPos, endPos, beginPos + (endPos - beginPos) / 2,
                                                leftWeight, rightWeight, leftWeight + rightWeight));
    }

    std::sort(_result.begin(), _result.end());
}

std::vector<ClippingClusterRecord> ClippingClusterAlgoImpl::result()
{
    return _result;
}

// ----------------------------------------------------------------------------
// Class ClippingClusterAlgo
// ----------------------------------------------------------------------------

ClippingClusterAlgo::ClippingClusterAlgo(ClippingClusterOptions const & options) : impl(new ClippingClusterAlgoImpl(options))
{}

ClippingClusterAlgo::~ClippingClusterAlgo()
{}

void ClippingClusterAlgo::push(std::vector<seqan::BamAlignmentRecord *> const & records)
{
    impl->push(records);
}

void ClippingClusterAlgo::finish()
{
    impl->finish();
}

std::vector<ClippingClusterRecord> ClippingClusterAlgo::result()
{
    return impl->result();
}
