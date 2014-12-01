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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_ASM_BEST_OVERLAP_GRAPH_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_ASM_BEST_OVERLAP_GRAPH_H_

#include <cstddef>
#include <functional>
#include <iosfwd>
#include <memory>
#include <map>
#include <vector>

namespace assembler {

// ============================================================================
// Forwards
// ============================================================================

class ContainmentFilterImpl;
class Overlap;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// UnitiggerOptions
// ----------------------------------------------------------------------------

struct UnitiggerOptions
{
    explicit UnitiggerOptions(unsigned longBogPathMinLength = 200,
                              unsigned longBogPathMinReadCount = 3) :
            longBogPathMinLength(longBogPathMinLength),
            longBogPathMinReadCount(longBogPathMinReadCount)
    {}

    // Minimal length for a path to be called "long".
    unsigned longBogPathMinLength { 200 };
    // Minimal number of reads to be in one path to consider it "long".
    unsigned longBogPathMinReadCount { 3 };

    // Recursion depth to explore in MIRA-like approach.
    unsigned miraRecursionDepth { 4 };
    // Number of best edges to traverse.
    unsigned miraNumEdges { 4 };
};

// ----------------------------------------------------------------------------
// Class ContainmentFilter
// ----------------------------------------------------------------------------

class ContainmentFilter
{
public:
    ContainmentFilter(bool logging);
    ~ContainmentFilter();  // for pimpl

    // Run the filter and filter out reads that are contained in others.  In case of "stacking" overlaps, complete
    // overlaps, pick the one with the smallest ID.
    std::vector<Overlap> run(size_t numReads, std::vector<Overlap> const & in);

    // Returns a vector with the contained sorted IDs.
    std::vector<unsigned> const & containedIDs() const;
    // Returns true iff the read is contained.
    bool isContained(unsigned readID) const;
    // Returns the id of the container.  This will not be a contained read.
    unsigned container(unsigned readID) const;
    // Returns the overlap with the container.
    Overlap overlap(unsigned readID) const;

private:
    std::unique_ptr<ContainmentFilterImpl> impl;
};

// ----------------------------------------------------------------------------
// Class BestOverlapGraph
// ----------------------------------------------------------------------------

// Stores the best overlap to the left and to the right.

struct BestOverlapGraph
{
    // Store information about a (tentative) best overlap.  Allows comparison by length, ties are broken by score and
    // readID.
    struct OverlapInfo
    {
        // Constant for invalid reads.
        static const unsigned INVALID_READ_ID = -1;

        // The read id for the overlap.
        unsigned readID { INVALID_READ_ID };
        // The length of the overlap.
        int length { 0 };
        // The score of the overlap (-errors).
        int score { 0 };
        // The begin position of the right sequence (seq1) in the overlap.
        int begin1 { 0 };

        OverlapInfo() = default;

        OverlapInfo(unsigned readID, int length, int score, int begin1 = 0) :
                readID(readID), length(length), score(score), begin1(begin1)
        {}

        // Returns whether the overlap info refers to a valid or an invalid read.
        bool valid() const
        {
            return (readID != INVALID_READ_ID);
        }

        // Reset to invalid state.
        void reset()
        {
            readID = INVALID_READ_ID;
            length = 0;
            score = 0;
            begin1 = 0;
        }

        bool operator<(OverlapInfo const & other) const
        {
            return make_tuple() < other.make_tuple();
        }

        bool operator==(OverlapInfo const & other) const
        {
            return (make_tuple() == other.make_tuple()) && (begin1 == other.begin1);
        }

    private:
        std::tuple<int, int, unsigned> make_tuple() const
        {
            // "-begin1" since smalles first.
            return std::make_tuple(-begin1, score, readID);
        }
    };

    // Construct BestOverlapGraph with a given read/sequence count.
    explicit BestOverlapGraph(size_t numReads = 0) : bestLeft(numReads), bestRight(numReads)
    {}

    // Returns number of reads this BestOverlapGraph is for.
    size_t numReads() const
    {
        return bestLeft.size();
    }

    // "Pointer" to the best overlap to the left or right.
    std::vector<OverlapInfo> bestLeft, bestRight;
};

void print(std::ostream & out, BestOverlapGraph const & g);

// ----------------------------------------------------------------------------
// Class BogPathLengths
// ----------------------------------------------------------------------------

// Best path lengths in BestOverlapGraph.
//
// The counts do not include the current node.

struct BogPathLengths
{
    static const unsigned UNSET_PATH_LENGTH;

    explicit BogPathLengths(size_t numReads = 0) :
            left(numReads, UNSET_PATH_LENGTH), right(numReads, UNSET_PATH_LENGTH)
    {}

    void reinit(size_t numReads);

    unsigned score(size_t readID) const
    {
        return 1 + left[readID] + right[readID];
    }

    std::vector<unsigned> left, right;
};

// ----------------------------------------------------------------------------
// Class BogPath
// ----------------------------------------------------------------------------

// A path in the BestOverlapGraph.

struct BogPath
{
    struct OverlapInfo
    {
        OverlapInfo() = default;

        OverlapInfo(unsigned id0, unsigned id1, int begin1) :
                id0(id0), id1(id1), begin1(begin1)
        {}

        // Ids of first/left and second/right sequence.
        unsigned id0 { 0 };
        unsigned id1 { 0 };
        // Begin position of the right sequence on the path.
        int begin1 { -1 };

        bool operator==(OverlapInfo const & other) const
        {
            return makeTuple() == other.makeTuple();
        }

    private:
        std::pair<std::pair<unsigned, unsigned>, int> makeTuple() const
        {
            return std::make_pair(std::make_pair(id0, id1), begin1);
        }
    };

    // Update length from the overlap data and readLengths.
    void updateLength();

    // Return number of reads.
    size_t size() const
    {
        return reads.size();
    }

    // Returns whether path contains read in O(n).
    bool contains(unsigned readID)
    {
        for (auto id : reads)
            if (readID == id)
                return true;
        return false;
    }

    // The length of the path (in approximate bases).
    unsigned length { 0 };

    // The lengths of the reads.
    std::vector<unsigned> readLengths;
    // The ids of the reads along the path.
    std::vector<unsigned> reads;
    // Information about the overlaps along the path.
    std::vector<OverlapInfo> overlaps;
};

std::ostream & operator<<(std::ostream & out, BogPath::OverlapInfo const & info);
std::ostream & operator<<(std::ostream & out, BogPath const & bogPath);

// ----------------------------------------------------------------------------
// Class BogPathSet
// ----------------------------------------------------------------------------

// A set of paths in the BestOverlapGraph, also allow for splitting of paths.

struct BogPathSet
{
    unsigned registerPath(BogPath && path);

    size_t size() const
    {
        return paths.size();
    }

    std::vector<BogPath>::iterator begin() { return paths.begin(); }
    std::vector<BogPath>::const_iterator begin() const { return paths.begin(); }
    std::vector<BogPath>::iterator end() { return paths.end(); }
    std::vector<BogPath>::const_iterator end() const { return paths.end(); }

    // The set of paths in the ste.
    std::vector<BogPath> paths;
    // The id of the path for the given read.
    //
    // Note that we store only the id of the first path that this read was assigned to.  Thus, it is elementary to
    // perform path splitting later that ensures each read is only on one path.
    std::map<unsigned, unsigned> pathForRead;
};

// ----------------------------------------------------------------------------
// Class RemovedOverlapStore
// ----------------------------------------------------------------------------

// Store information about overlaps removed from BestOverlapGraph.  Required for adding them back later to path
// splitting is less biased.

struct RemovedOverlapStore
{
    enum Direction { INVALID, LEFT, RIGHT };

    // Stores information about one removed overlap info entry.
    struct Entry {
        unsigned readID { (unsigned)-1 };
        Direction direction { INVALID };
        BestOverlapGraph::OverlapInfo info;

        Entry() = default;
        Entry(unsigned readID, Direction direction, BestOverlapGraph::OverlapInfo const & info) :
                readID(readID), direction(direction), info(info)
        {}
    };

    void addEntry(unsigned readID, Direction direction, BestOverlapGraph::OverlapInfo const & info)
    {
        entries.push_back(Entry(readID, direction, info));
    }

    std::vector<Entry>::const_iterator begin() const { return entries.begin(); }
    std::vector<Entry>::const_iterator end() const { return entries.end(); }

    std::vector<Entry> entries;

    void print(std::ostream & out) const;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function operator<<()                        [BestOverlapGraph::OverlapInfo]
// ----------------------------------------------------------------------------

std::ostream & operator<<(std::ostream & out, BestOverlapGraph::OverlapInfo info);

// ----------------------------------------------------------------------------
// Function buildBestOverlapGraph()
// ----------------------------------------------------------------------------

// Build a best overlap graph from overlaps.

BestOverlapGraph buildBestOverlapGraph(size_t numReads,
                                       std::vector<Overlap> const & overlaps,
                                       bool logging = false);

// ----------------------------------------------------------------------------
// Function printDot()
// ----------------------------------------------------------------------------

void printDot(std::ostream & out, BestOverlapGraph const & bog);

// ----------------------------------------------------------------------------
// Function removeCycles()
// ----------------------------------------------------------------------------

// Remove cycles from BestOverlapGraph.
//
// The resulting missing overlaps are stored in ovls and can be added back later before path splitting.

void removeCycles(RemovedOverlapStore & ovls,
                  BestOverlapGraph & og);

// ----------------------------------------------------------------------------
// Function completePaths()
// ----------------------------------------------------------------------------

// Consider the removed overlaps in ovl.  In the case that they belong to a path ending at the read for which the
// overlap was removed, add back the other read together with the overlap so they can be used for path splitting or
// joining the contigs downstream through this overlap.

void completePaths(BogPathSet & paths,
                   RemovedOverlapStore const & ovls,
                   std::vector<unsigned> const & readLengths);

// ----------------------------------------------------------------------------
// Function computePathLengths()
// ----------------------------------------------------------------------------

BogPathLengths computePathLengths(BestOverlapGraph const & bog);

// ----------------------------------------------------------------------------
// Function buildInitialPaths()
// ----------------------------------------------------------------------------

// Build initial/"promiscuous" paths through BOG graph.

BogPathSet buildInitialPaths(std::vector<unsigned> const & readLengths,
                             BestOverlapGraph const & bog,
                             BogPathLengths const & pathLengths,
                             bool logging = false);

// ----------------------------------------------------------------------------
// Function performPathSplitting()
// ----------------------------------------------------------------------------

// Perform path splitting as described in CABOG paper.
//
// Also ensures that each read is only on one read.

void performPathSplitting(BogPathSet & pathSet,
                          ContainmentFilter const & filter,
                          UnitiggerOptions const & options,
                          bool logging = false);

}  // namespace assembler

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_ASM_BEST_OVERLAP_GRAPH_H_
