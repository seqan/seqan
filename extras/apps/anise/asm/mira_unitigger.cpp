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

#include "mira_unitigger.h"

#include <iostream>
#include <vector>

#include <seqan/file.h>

#include "asm/best_overlap_graph.h"

namespace assembler {

namespace {  // anonymous namespace

// --------------------------------------------------------------------------
// Class AnnotatedOverlap
// --------------------------------------------------------------------------

// Wrapper for normalized overlap that adds an id to the overlap (for "used" flag).

struct AnnotatedOverlap
{
    static const unsigned INVALID = (unsigned)-1;

    AnnotatedOverlap() = default;
    AnnotatedOverlap(Overlap const & ovl, unsigned id) : ovl(ovl.normalize()), id(id) {}

    // Returns approximate score from overlap length and error count (approximately +1 for match, -1 for any error).
    int score() const
    {
        // Overlap is normalized, no need to call ovl.length() which would get normalized copy.
        return length();// - ovl.errors;
    }

    // Returns overlap length.
    int length() const
    {
        return ovl.len0 - ovl.begin1;
    }

    // Returns whether is a containment or stack overlap.
    bool isContainment() const
    {
        return (ovl.begin1 + ovl.len1 <= ovl.len0);  // INVARIANT: (ovl.begin1 >= 0)
    }

    Overlap ovl;
    unsigned id { INVALID };
};

// Write to out stream.
std::ostream & operator<<(std::ostream & out, AnnotatedOverlap const & ann)
{
    return out << "AnnotatedOverlap(ovl=" << ann.ovl << ", id=" << ann.id << ", score()=" << ann.score() << ")";
}

// Compares two annotated overlaps (seq0, -score) such that they can be used in an adjacency array representation
// for "forward overlap" edges.
bool ltSeq0(AnnotatedOverlap const & lhs, AnnotatedOverlap const & rhs)
{
    return (std::make_tuple(lhs.ovl.seq0, -lhs.score()) <
            std::make_tuple(rhs.ovl.seq0, -rhs.score()));
}

// Compares two annotated overlaps (seq1, -score) such that they can be used in an adjacency array representation
// for "backward overlap" edges.
bool ltSeq1(AnnotatedOverlap const & lhs, AnnotatedOverlap const & rhs)
{
    return (std::make_tuple(lhs.ovl.seq1, -lhs.score()) <
            std::make_tuple(rhs.ovl.seq1, -rhs.score()));
}

// Compares two annotated overlap by (score, key) such that the highest-scoring overlaps come first.
bool ltScore(AnnotatedOverlap const & lhs, AnnotatedOverlap const & rhs)
{
    return (std::make_pair(-lhs.score(), lhs.ovl.key()) < std::make_pair(-rhs.score(), rhs.ovl.key()));
}

// --------------------------------------------------------------------------
// Class MiraGraph
// --------------------------------------------------------------------------

// Compact representation for the MIRA-like unitigging.

class MiraGraph
{
public:
    MiraGraph(std::vector<Overlap> const & overlaps, unsigned numSeqs) : _numSeqs(numSeqs)
    {
        buildLeftBuckets(overlaps);
        buildRightBuckets(overlaps);
    }

    // Print MIRA Graph.
    void print(std::ostream & out) const;

    unsigned numSeqs() const { return _numSeqs; }

    // Return range of overlaps to the left/right for the given read.
    typedef std::vector<AnnotatedOverlap>::const_iterator TIter;
    std::pair<TIter, TIter> leftBucket(unsigned readID) const
    {
        return std::make_pair(leftOverlaps.begin() + leftBuckets.at(readID),
                              leftOverlaps.begin() + leftBuckets.at(readID + 1));
    }
    std::pair<TIter, TIter> rightBucket(unsigned readID) const
    {
        return std::make_pair(rightOverlaps.begin() + rightBuckets.at(readID),
                              rightOverlaps.begin() + rightBuckets.at(readID + 1));
    }

private:

    // Build leftBuckets and leftOverlaps.
    void buildLeftBuckets(std::vector<Overlap> const & overlaps);
    // Build rightBuckets and rightOverlaps.
    void buildRightBuckets(std::vector<Overlap> const & overlaps);

    // Number of sequences.
    unsigned _numSeqs;

    // Begin/end borders of buckets in leftOverlaps.
    std::vector<unsigned> leftBuckets;
    // Overlaps to the left.
    std::vector<AnnotatedOverlap> leftOverlaps;
    // Begin/end borders of buckets in rightOverlaps.
    std::vector<unsigned> rightBuckets;
    // Overlaps to the right.
    std::vector<AnnotatedOverlap> rightOverlaps;
};

void MiraGraph::buildRightBuckets(std::vector<Overlap> const & overlaps)
{
    // Copy into AnnotatedOverlaps.
    for (unsigned ovlID = 0; ovlID < overlaps.size(); ++ovlID)
        rightOverlaps.push_back(AnnotatedOverlap(overlaps[ovlID], ovlID));
    // Sort into correct order.
    std::sort(rightOverlaps.begin(), rightOverlaps.end(), ltSeq0);

    auto it = rightOverlaps.begin();
    for (unsigned seqID = 0; seqID < _numSeqs; ++seqID)
    {
        for (; it != rightOverlaps.end() && it->ovl.seq0 < seqID; ++it)
            continue;  // Skip over.
        rightBuckets.push_back(it - rightOverlaps.begin());
    }
    while (rightBuckets.size() <= _numSeqs)
        rightBuckets.push_back(rightOverlaps.size());
}

void MiraGraph::buildLeftBuckets(std::vector<Overlap> const & overlaps)
{
    // Copy into AnnotatedOverlaps.
    for (unsigned ovlID = 0; ovlID < overlaps.size(); ++ovlID)
        leftOverlaps.push_back(AnnotatedOverlap(overlaps[ovlID], ovlID));
    // Sort into correct order.
    std::sort(leftOverlaps.begin(), leftOverlaps.end(), ltSeq1);

    auto it = leftOverlaps.begin();
    for (unsigned seqID = 0; seqID < _numSeqs; ++seqID)
    {
        for (; it != leftOverlaps.end() && it->ovl.seq1 < seqID; ++it)
            continue;  // Skip over.
        leftBuckets.push_back(it - leftOverlaps.begin());
    }
    while (leftBuckets.size() <= _numSeqs)
        leftBuckets.push_back(leftOverlaps.size());
}

void MiraGraph::print(std::ostream & out) const
{
    out << "MIRA GRAPH\n"
        << "\n"
        << "LEFT\n";
    for (unsigned i = 0; i < _numSeqs; ++i)
    {
        out << i << ":\n";
        for (unsigned j = leftBuckets[i]; j < leftBuckets[i + 1]; ++j)
            out << "\t" << leftOverlaps[j] << "\n";
    }
    out << "RIGHT\n";
    for (unsigned i = 0; i < _numSeqs; ++i)
    {
        out << i << ":\n";
        for (unsigned j = rightBuckets[i]; j < rightBuckets[i + 1]; ++j)
            out << "\t" << rightOverlaps[j] << "\n";
    }
}

// --------------------------------------------------------------------------
// Class VectorSet
// --------------------------------------------------------------------------

// Small helper of unordered vector set.

struct VectorSet
{
    VectorSet() = default;

    void push(unsigned v) { items.push_back(v); }
    void pop() { items.pop_back(); }

    size_t size() const { return items.size(); }
    size_t empty() const { return items.empty(); }
    bool contains(unsigned v) const { return (find(v) != end()); }
    std::vector<unsigned>::const_iterator find(unsigned v) const { return std::find(begin(), end(), v); }

    std::vector<unsigned>::iterator begin() { return items.begin(); }
    std::vector<unsigned>::iterator end()   { return items.end(); }
    std::vector<unsigned>::const_iterator begin() const { return items.begin(); }
    std::vector<unsigned>::const_iterator end() const   { return items.end(); }

    std::vector<unsigned> items;
};

// --------------------------------------------------------------------------
// Class MiraBogBuilder
// --------------------------------------------------------------------------

// Build BOG from MiraGraph.

class MiraBogBuilder
{
public:
    MiraBogBuilder(MiraGraph const & miraGraph,
                   std::vector<Overlap> const & overlaps,
                   UnitiggerOptions const & options,
                   bool logging) :
            miraGraph(miraGraph), options(options), logging(logging), result(miraGraph.numSeqs())
    {
        init(overlaps);
    }

    BestOverlapGraph run();

private:

    // Initialize sortedOverlaps and reached.
    void init(std::vector<Overlap> const & overlaps);

    // Explore graph starting at the given overlap.
    void exploreGraphStartingAt(AnnotatedOverlap const & ann);

    // Recurse left / right from the given read.
    void recurseLeftFrom(unsigned readID);
    void recurseRightFrom(unsigned readID);

    // Recurse left / right to compute best score.  tmpReached is used to keep a patch to the "reached" map, all nodes
    // marked as reached and that are in tmpReached are not touched.
    int leftScoreRecursion(unsigned readID, unsigned counter);
    int rightScoreRecursion(unsigned readID, unsigned counter);

    // Input

    // MiraGraph to explore.
    MiraGraph const & miraGraph;
    // Configuration.
    UnitiggerOptions const & options;
    // Whether to enable/disable logging.
    bool logging;

    // Result
    
    BestOverlapGraph result;

    // State

    // Overlaps, sorted by score.
    std::vector<AnnotatedOverlap> sortedOverlaps;
    // Map for each read whether it has been reached yet.
    std::vector<bool> reached;

    // Temporary patch to reached map.
    VectorSet tmpReached;
};

BestOverlapGraph MiraBogBuilder::run()
{
    for (auto const & ann : sortedOverlaps)
        exploreGraphStartingAt(ann);

    return std::move(result);
}

int MiraBogBuilder::leftScoreRecursion(unsigned readID, unsigned counter)
{
    if (!counter)
        return 0;  // Break recursion for depth.

    // Compute best score when proceeding left.
    auto range = miraGraph.leftBucket(readID);
    unsigned tries = 0;  // number of tried edges
    int bestScore = 0;
    for (auto it = range.first; it != range.second && tries < options.miraNumEdges; ++it)
    {
        SEQAN_ASSERT_EQ(it->ovl.seq1, readID);
        if (tmpReached.contains(it->ovl.seq0))
            continue;  // Skip, in current recursion
        if (!it->isContainment())
            ++tries;  // considering edge, do not increase counter in case of containment

        // Compute best score along this edge.
        tmpReached.push(it->ovl.seq0);
        int score = it->score() + leftScoreRecursion(it->ovl.seq0, counter - 1);
        tmpReached.pop();

        bestScore = std::max(score, bestScore);
    }

    return bestScore;
}

int MiraBogBuilder::rightScoreRecursion(unsigned readID, unsigned counter)
{
    if (!counter)
        return 0;  // Break recursion for depth.

    // Compute best score when proceeding right.
    auto range = miraGraph.rightBucket(readID);
    unsigned tries = 0;  // number of tried edges
    int bestScore = 0;
    for (auto it = range.first; it != range.second && tries < options.miraNumEdges; ++it)
    {
        SEQAN_ASSERT_EQ(it->ovl.seq0, readID);
        if (tmpReached.contains(it->ovl.seq1))
            continue;  // Skip, in current recursion
        if (!it->isContainment())
            ++tries;  // considering edge, do not increase counter in case of countainment

        // Compute best score along this edge.
        tmpReached.push(it->ovl.seq1);
        int score = it->score() + rightScoreRecursion(it->ovl.seq1, counter - 1);
        tmpReached.pop();

        bestScore = std::max(score, bestScore);
    }

    return bestScore;
}

void MiraBogBuilder::exploreGraphStartingAt(AnnotatedOverlap const & ann)
{
    if (logging)
        std::cerr << "Starting exploration from " << ann << "\n";

    if (reached.at(ann.ovl.seq0) && reached.at(ann.ovl.seq1))
        return;  // Nothing to explore if connecting two already reached reads.

    if (reached.at(ann.ovl.seq0))
    {
        if (logging)
            std::cerr << "Only recursing right\n";
        // Left sequence in ovl has already been reached.  Update best left of seq1 to seq0 and start recursion left.
        result.bestLeft[ann.ovl.seq1] = BestOverlapGraph::OverlapInfo(
                ann.ovl.seq0, ann.length(), ann.score(), ann.ovl.begin1);
        reached.at(ann.ovl.seq1) = true;
        recurseRightFrom(ann.ovl.seq1);
    }
    else if (reached.at(ann.ovl.seq1))
    {
        if (logging)
            std::cerr << "Only recursing left\n";
        // Right sequence has already been reached.  Update best right of seq0 to seq1 and start recursion right.
        result.bestRight[ann.ovl.seq0] = BestOverlapGraph::OverlapInfo(
                ann.ovl.seq1, ann.length(), ann.score(), ann.ovl.begin1);
        reached.at(ann.ovl.seq0) = true;
        recurseLeftFrom(ann.ovl.seq0);
    }
    else
    {
        if (logging)
            std::cerr << "Recursion in both directions\n";
        // Both left and right are unreached.  Add bidirectional link and recurse both to the left and to the right.
        result.bestLeft[ann.ovl.seq1] = BestOverlapGraph::OverlapInfo(
                ann.ovl.seq0, ann.length(), ann.score(), ann.ovl.begin1);
        result.bestRight[ann.ovl.seq0] = BestOverlapGraph::OverlapInfo(
                ann.ovl.seq1, ann.length(), ann.score(), ann.ovl.begin1);
        reached.at(ann.ovl.seq0) = true;
        reached.at(ann.ovl.seq1) = true;
        recurseLeftFrom(ann.ovl.seq0);
        recurseRightFrom(ann.ovl.seq1);
    }
}

void MiraBogBuilder::init(std::vector<Overlap> const & overlaps)
{
    reached.resize(miraGraph.numSeqs(), false);

    // Created list of overlaps sorted by score.
    for (unsigned i = 0; i < overlaps.size(); ++i)
        sortedOverlaps.push_back(AnnotatedOverlap(overlaps[i], i));
    std::sort(sortedOverlaps.begin(), sortedOverlaps.end(), ltScore);
}

void MiraBogBuilder::recurseLeftFrom(unsigned readID)
{
    SEQAN_ASSERT(tmpReached.empty());
    if (logging)
        std::cerr << "Recursing left from " << readID << "\n";

    // Pick best overlap to the left.
    auto range = miraGraph.leftBucket(readID);
    unsigned tries = 0;  // number of tried edges
    auto best = range.second;  // no matching found
    int bestScore = 0;
    for (auto it = range.first; it != range.second && tries < options.miraNumEdges; ++it)
    {
        SEQAN_ASSERT_EQ(it->ovl.seq1, readID);
        if (!it->isContainment())
            ++tries;

        // Compute best score along this edge.
        SEQAN_ASSERT(tmpReached.empty());
        tmpReached.push(it->ovl.seq0);
        int score = it->score() + leftScoreRecursion(it->ovl.seq0, options.miraRecursionDepth);
        tmpReached.pop();

        // Pick current as best if non picked yet or better than previously best.
        if (best == range.second || score > bestScore)
        {
            bestScore = score;
            best = it;
        }
    }

    if (best == range.second)  // found no destination, break recursion
        return;

    // Recurse along best overlap to the left.
    if (reached.at(best->ovl.seq0))
    {
        if (logging)
            std::cerr << "  => stopping, already reached " << best->ovl.seq0 << "\n";
        // Already reached, simply update the left-to-left link for readID and break recursion.
        result.bestLeft[readID] = BestOverlapGraph::OverlapInfo(
                best->ovl.seq0, best->length(), best->score(), best->ovl.begin1);
    }
    else
    {
        if (logging)
            std::cerr << "  => proceeding left (with link)\n";
        // Otherwise, update left-to-left link of this and left-to-left link of best->ovl.seq0, thus fixing this
        // overlap.  Then, mark best->ovl.seq0 as reached and recurse along edge.
        result.bestRight[best->ovl.seq0] = BestOverlapGraph::OverlapInfo(
                readID, best->length(), best->score(), best->ovl.begin1);
        result.bestLeft[readID] = BestOverlapGraph::OverlapInfo(
                best->ovl.seq0, best->length(), best->score(), best->ovl.begin1);
        reached.at(best->ovl.seq0) = true;
        recurseLeftFrom(best->ovl.seq0);
    }
}

void MiraBogBuilder::recurseRightFrom(unsigned readID)
{
    SEQAN_ASSERT(tmpReached.empty());
    if (logging)
        std::cerr << "Recursing right from " << readID << "\n";

    // Pick best overlap to the right.
    auto range = miraGraph.rightBucket(readID);
    unsigned tries = 0;  // number of tried edges
    auto best = range.second;  // no matching found
    int bestScore = 0;
    for (auto it = range.first; it != range.second && tries < options.miraNumEdges; ++it, ++tries)
    {
        SEQAN_ASSERT_EQ(it->ovl.seq0, readID);
        if (!it->isContainment())
            ++tries;

        // Compute best score along this edge.
        SEQAN_ASSERT(tmpReached.empty());
        tmpReached.push(it->ovl.seq1);
        int score = it->score() + rightScoreRecursion(it->ovl.seq1, options.miraRecursionDepth);
        tmpReached.pop();

        // Pick current as best if non picked yet or better than previously best.
        if (best == range.second || score > bestScore)
        {
            bestScore = score;
            best = it;
        }
    }

    if (best == range.second)  // found no destination, break recursion
        return;

    // Recurse along best overlap to the right.
    if (reached.at(best->ovl.seq1))
    {
        if (logging)
            std::cerr << "  => stopping, already reached " << best->ovl.seq0 << "\n";
        // Already reached, simply update the left-to-right link for readID and break recursion.
        result.bestRight[readID] = BestOverlapGraph::OverlapInfo(
                best->ovl.seq1, best->length(), best->score(), best->ovl.begin1);
    }
    else
    {
        if (logging)
            std::cerr << "  => proceeding right (with link)\n";
        // Otherwise, update left-to-right link of this and right-to-left link of best->ovl.seq1, thus fixing this
        // overlap.  Then, mark best->ovl.seq1 as reached and recurse along edge.
        result.bestRight[readID] = BestOverlapGraph::OverlapInfo(
                best->ovl.seq1, best->length(), best->score(), best->ovl.begin1);
        result.bestLeft[best->ovl.seq1] = BestOverlapGraph::OverlapInfo(
                readID, best->length(), best->score(), best->ovl.begin1);
        reached.at(best->ovl.seq1) = true;
        recurseRightFrom(best->ovl.seq1);
    }
}

}  // anonymous namespace

// --------------------------------------------------------------------------
// Class MiraUnitigger
// --------------------------------------------------------------------------

std::unique_ptr<ContigGraph> MiraUnitigger::run(seqan::StringSet<seqan::Dna5String> const & seqs,
                                                std::vector<Overlap> const & overlaps) const
{
    if (logging)
        std::cerr << "Filtering out containment overlaps...\n";
    ContainmentFilter filter(logging);
    std::vector<Overlap> filteredOverlaps = filter.run(length(seqs), overlaps);

    if (logging)
        std::cerr << "MIRA-like Unitigging...\n";

    // Build MIRA graph.
    if (logging)
        std::cerr << "MIRA graph\n";
    MiraGraph mg(filteredOverlaps, length(seqs));
    if (logging)
        mg.print(std::cerr);

    // Perform MIRA path exploration, yielding a BestOverlapGraph with paths.
    UnitiggerOptions uo;
    MiraBogBuilder builder(mg, filteredOverlaps, uo, logging);
    BestOverlapGraph bog = builder.run();

    // Now, perform BOG path splitting.
    if (logging)
    {
        std::cerr << "Before cycle removal\n";
        print(std::cerr, bog);
    }
    RemovedOverlapStore removedOvls;
    removeCycles(removedOvls, bog);
    if (logging)
    {
        std::cerr << "After cycle removal\n";
        print(std::cerr, bog);
        std::cerr << "REMOVED OVERLAP STORE\n";
        removedOvls.print(std::cerr);
    }
    BogPathLengths pathLengths = computePathLengths(bog);
    auto readLengths = computeReadLengths(seqs);
    BogPathSet bps = buildInitialPaths(readLengths, bog, pathLengths, logging);
    if (logging)
    {
        std::cerr << "Paths before completion (adding formerly removed overlaps)\n";
        for (unsigned i = 0; i < bps.size(); ++i)
        {
            std::cerr << "path #" << i << "\n";
            std::cerr << bps.paths[i] << "\n";
        }
    }
    completePaths(bps, removedOvls, readLengths);
    // Perform path splitting, we can build our results from that.
    if (logging)
        std::cerr << "Performing Path Splitting...\n";
    if (logging)
    {
        std::cerr << "Paths before splitting\n";
        for (unsigned i = 0; i < bps.size(); ++i)
        {
            std::cerr << "path #" << i << "\n";
            std::cerr << bps.paths[i] << "\n";
        }
    }
    UnitiggerOptions options;  // TODO(holtgrew): Make configurable from the outside.
    performPathSplitting(bps, filter, options, logging);
    if (logging)
    {
        std::cerr << "Paths after splitting\n";
        for (unsigned i = 0; i < bps.size(); ++i)
        {
            std::cerr << "path #" << i << "\n";
            std::cerr << bps.paths[i] << "\n";
        }
    }

    // Compute contig graph without edges.
    auto result = augmentContigGraph(computeContigGraph(bps, filter, filteredOverlaps, seqs), filter, seqs);
    auto & out = *result;

    // Add arcs/contig overlaps to ContigGraph.
    //
    // First collect first/last read of for each contig.
    std::map<unsigned, unsigned> firstReadOf, lastReadOf;
    for (unsigned i = 0; i < out.node.size(); ++i)
    {
        // std::cerr << "firstReadOf[" << out.contig[out.node[i]].posReads.front().readID << "] = " << i << "\n"
        //           << "lastReadOf[" << out.contig[out.node[i]].posReads.back().readID << "] = " << i << "\n";
        firstReadOf[out.contig[out.node[i]].posReads.front().readID] = i;
        lastReadOf[out.contig[out.node[i]].posReads.back().readID] = i;
    }
    // Consider all overlaps.  Collect overlaps between first/last reads on a FCFS basis.
    std::map<std::pair<unsigned, unsigned>, Overlap> firstLastOverlap;
    auto handleOverlap = [&](Overlap const & ovl)  // ovl left-to-right or stack
    {
        auto itL = lastReadOf.find(ovl.seq0);
        auto itF = firstReadOf.find(ovl.seq1);
        if (itL == lastReadOf.end() || itF == firstReadOf.end() || itL->second == itF->second)
            return;  // not first and last or of the same contig
        if (firstLastOverlap.count(std::make_pair(itL->second, itF->second)))
            return;  // we already have an overlap for this pair
        firstLastOverlap[std::make_pair(itL->second, itF->second)] = ovl;
    };
    for (auto const & ovl : overlaps)
    {
        if (ovl.begin0 == 0)
            handleOverlap(ovl);
        if (ovl.begin1 == 0)
            handleOverlap(ovl.flip());
    }
    for (auto const & contigOverlaps : firstLastOverlap)
    {
        auto idU = contigOverlaps.first.first;
        auto idV = contigOverlaps.first.second;
        SEQAN_CHECK(idU != idV, "idU = %u, idV = %u", idU, idV);
        auto u = out.node.at(idU);
        auto v = out.node.at(idV);
        SEQAN_CHECK(u != v, "idU = %u, idV = %u", idU, idV);
        auto arc = out.graph.addArc(u, v);
        out.overlapLength[arc] = contigOverlaps.second.len0 - contigOverlaps.second.begin1;
    }

    return result;
}

}  // namespace assembler
