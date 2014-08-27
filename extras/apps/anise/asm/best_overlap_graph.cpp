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

#include "best_overlap_graph.h"

#include <iterator>
#include <map>

#include <seqan/file.h>  // for printing seqan::String

#include <seqan/graph_types.h>  // ugly dependency of misc_union_find.h
#include <seqan/misc/misc_union_find.h>

#include <lemon/smart_graph.h>
#include <lemon/dfs.h>

#include "asm/overlapper.h"

namespace assembler {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Function eraseIf()
// ----------------------------------------------------------------------------

template <typename TValue, typename TPred>
void eraseIf(std::vector<TValue> & v, TPred pred)
{
    v.erase(std::remove_if(v.begin(), v.end(), pred), v.end());
}

// ----------------------------------------------------------------------------
// Function isContainment()
// ----------------------------------------------------------------------------

// Returns true if seq1 is contained in seq0 in the overlap.

inline bool isContainment(Overlap const & ovl)
{
    SEQAN_ASSERT_EQ(ovl.begin0, 0u);
    return (ovl.begin1 + ovl.len1 <= ovl.begin0 + ovl.len0);
}

// --------------------------------------------------------------------------
// Function infoFromOverlap()
// --------------------------------------------------------------------------

BestOverlapGraph::OverlapInfo infoFromOverlap(Overlap const & ovl, bool left)
{
    return BestOverlapGraph::OverlapInfo(left ? ovl.seq0 : ovl.seq1, ovl.length(), -ovl.errors, ovl.begin1);
}

// --------------------------------------------------------------------------
// Class BogCycleRemover
// --------------------------------------------------------------------------

// Cycle removal using DFS.  DFS edges first taken are kept (and their parallel back edge from the "opposite site" if
// any exists.

class BogCycleRemover
{
public:
    // An edge in the Best Overlap Graph.  Sequence SEQ is represented by the nodes 2*SEQ and 2*SEQ+1 which are
    // connected by an edge.  Links to the right go from 2*FROM+1 to 2*TO while links to the left go from 2*FROM to
    // 2*TO+1.
    //
    // We store the link such that node0 < node1.
    struct Link
    {
        unsigned node0  { 0 };
        unsigned node1  { 0 };
        unsigned weight { 0 };

        Link() = default;

        Link(unsigned node0, unsigned node1, unsigned weight = 1) :
                node0(std::min(node0, node1)), node1(std::max(node0, node1)), weight(weight)
        {}

        std::pair<unsigned, unsigned> key() const
        {
            return std::make_pair(node0, node1);
        }

        bool operator<(Link const & other) const
        {
            return makeTuple() > other.makeTuple();
        }

    private:

        std::tuple<unsigned, unsigned, unsigned> makeTuple() const
        {
            return std::make_tuple(weight, node0, node1);
        }
    };

    BogCycleRemover(RemovedOverlapStore & removedOvls, BestOverlapGraph & graph) :
            removedOvls(removedOvls), graph(graph)
    {}

    // Perform the removal.
    void run();

private:

    // Compute links sorted descendingly by weight.
    std::vector<Link> sortedLinks() const;

    // Overlaps that were previously removed.
    RemovedOverlapStore & removedOvls;
    // Reference to the BestOverlapGraph to work on.
    BestOverlapGraph & graph;
};

// std::ostream & operator<<(std::ostream & out, BogCycleRemover::Link const & link)
// {
//     return out << "BogCycleRemover::Link(node0=" << link.node0 << ", node1=" << link.node1
//                << ", weight=" << link.weight << ")";
// }

std::vector<BogCycleRemover::Link> BogCycleRemover::sortedLinks() const
{
    std::vector<Link> result;

    // Collect links.
    for (unsigned i = 0; i < graph.bestLeft.size(); ++i)
    {
        if (graph.bestLeft[i].valid())
            result.push_back(Link(2 * i, 2 * graph.bestLeft[i].readID + 1));
        if (graph.bestRight[i].valid())
            result.push_back(Link(2 * i + 1, 2 * graph.bestRight[i].readID));
    }

    // Sort links.
    if (result.size() <= 1u)
        return result;
    std::stable_sort(result.begin(), result.end());

    // Merge adjacent links.
    auto it = result.begin(), itNext = std::next(it);
    for (; itNext != result.end(); ++itNext)
    {
        if (it->key() == itNext->key())
            ++(it->weight);
        else
            *(++it) = *itNext;
    };
    result.resize(std::next(it) - result.begin());

    std::stable_sort(result.begin(), result.end());

    return result;
}

void BogCycleRemover::run()
{
    // Obtain links sorted descendingly by weight.
    auto links = sortedLinks();

    // Find links closing a cycle and remove them from the graph.
    seqan::UnionFind<unsigned> uf;
    resize(uf, 2 * graph.bestRight.size());

    // Pre-join all sequence edges.
    for (unsigned i = 0; i < graph.bestRight.size(); ++i)
        joinSets(uf, findSet(uf, 2 * i), findSet(uf, 2 * i + 1));

    // Perform joining with the links.
    for (Link link : links)
    {
        unsigned repr0 = findSet(uf, link.node0);
        unsigned repr1 = findSet(uf, link.node1);
        if (repr0 != repr1)
        {
            // No cycles closed, join sets.
            joinSets(uf, repr0, repr1);
        }
        else
        {
            SEQAN_ASSERT_NEQ(link.node0 % 2u, link.node1 % 2u);
            unsigned read0 = link.node0 / 2;
            unsigned read1 = link.node1 / 2;
            if (link.node0 % 2u)
                std::swap(read0, read1);

            // Cycle closed, remove link.
            if (graph.bestLeft[read0].valid() && (graph.bestLeft[read0].readID == read1))
            {
                removedOvls.addEntry(read0, RemovedOverlapStore::LEFT, graph.bestLeft[read0]);
                graph.bestLeft[read0].reset();
            }
            if (graph.bestRight[read1].valid() && (graph.bestRight[read1].readID == read0))
            {
                removedOvls.addEntry(read1, RemovedOverlapStore::RIGHT, graph.bestRight[read1]);
                graph.bestRight[read1].reset();
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Class InitialBogPathBuilder
// ----------------------------------------------------------------------------

class InitialBogPathBuilder
{
public:
    InitialBogPathBuilder(std::vector<unsigned> const & readLengths,
                          BestOverlapGraph const & bog,
                          BogPathLengths const & pathLengths,
                          bool logging) :
            readLengths(readLengths), bog(bog), pathLengths(pathLengths), logging(logging),
            nodeToPath(bog.numReads(), INVALID_PATH)
    {}

    BogPathSet run();

private:
    // Compute sorted read ids by their score.
    std::vector<unsigned> rankedReadIDs() const;

    // Collect path starting at readID.
    BogPath collectPathFrom(unsigned readID) const;

    // Input

    // The length of each read.
    std::vector<unsigned> const & readLengths;
    // Best overlap graph.
    BestOverlapGraph const & bog;
    // Path lengths info.
    BogPathLengths const & pathLengths;
    // Enable/disable logging.
    bool logging;

    // State

    // Mapping of the nodes to the path if any.
    static const unsigned INVALID_PATH;
    std::vector<unsigned> nodeToPath;
};

const unsigned InitialBogPathBuilder::INVALID_PATH = -1;

std::vector<unsigned> InitialBogPathBuilder::rankedReadIDs() const
{
    std::vector<unsigned> result(bog.numReads());
    std::iota(result.begin(), result.end(), 0);

    auto cmp = [&](unsigned lhs, unsigned rhs) { return pathLengths.score(lhs) > pathLengths.score(rhs); };
    std::stable_sort(result.begin(), result.end(), cmp);

    return result;
}

BogPath InitialBogPathBuilder::collectPathFrom(unsigned readID) const
{
    if (logging)
        std::cerr <<  "collectPathFrom(" << readID << ")\n";
    BogPath path;
    SEQAN_ASSERT_EQ(nodeToPath[readID], INVALID_PATH);

    // Get read ids and overlaps to the left of readID.
    for (unsigned current = readID, left = bog.bestLeft.at(readID).readID;
         left != (unsigned)-1;
         current = left, left = bog.bestLeft.at(left).readID)
    {
        path.reads.push_back(left);
        path.readLengths.push_back(readLengths[left]);
        path.overlaps.push_back(BogPath::OverlapInfo(left, current, bog.bestLeft.at(current).begin1));
        if (nodeToPath[left] != INVALID_PATH)
            break;  // reached node that is part of another path
    }
    std::reverse(path.reads.begin(), path.reads.end());
    std::reverse(path.readLengths.begin(), path.readLengths.end());
    std::reverse(path.overlaps.begin(), path.overlaps.end());

    // Add current read.
    path.readLengths.push_back(readLengths[readID]);
    path.reads.push_back(readID);

    // Get read ids and overlaps to the right of readID.
    for (unsigned current = readID, right = bog.bestRight.at(readID).readID;
         right != (unsigned)-1;
         current = right, right = bog.bestRight.at(right).readID)
    {
        path.reads.push_back(right);
        path.readLengths.push_back(readLengths[right]);
        path.overlaps.push_back(BogPath::OverlapInfo(current, right, bog.bestRight.at(current).begin1));
        if (nodeToPath[right] != INVALID_PATH)
            break;  // reached node that is part of another path
    }

    path.updateLength();

    return path;
}

BogPathSet InitialBogPathBuilder::run()
{
    BogPathSet result;

    auto sortedReadIDs = rankedReadIDs();
    for (auto readID : sortedReadIDs)
        if (nodeToPath[readID] == INVALID_PATH)
        {
            unsigned pathID = result.registerPath(collectPathFrom(readID));
            for (auto readID : result.paths[pathID].reads)
                nodeToPath[readID] = pathID;
        }

    return result;
}

// ----------------------------------------------------------------------------
// Class BogPathSplitter
// ----------------------------------------------------------------------------

class BogPathSplitter
{
public:
    BogPathSplitter(BogPathSet & pathSet, ContainmentFilter const & filter,
                    UnitiggerOptions const & options, bool logging) :
            pathSet(pathSet), filter(filter), options(options), logging(logging)
    {}

    void run();

private:

    // Remove first read and overlap if any.
    void removeFirstRead(BogPath & path)
    {
        if (path.reads.empty())
            return;
        path.reads.erase(path.reads.begin());
        path.readLengths.erase(path.readLengths.begin());
        if (!path.overlaps.empty())
            path.overlaps.erase(path.overlaps.begin());
        path.updateLength();
    }

    // Remove last read and overlap if any.
    void removeLastRead(BogPath & path)
    {
        if (path.reads.empty())
            return;
        path.reads.pop_back();
        path.readLengths.pop_back();
        if (!path.overlaps.empty())
            path.overlaps.pop_back();
        path.updateLength();
    }

    // Split the path with pathID at the read with readID when the converging path comes from the right.
    void splitFromRight(unsigned pathID, unsigned readID);

    // Split the path with pathID at the read with readID when the converging path comes from the left.
    void splitFromLeft(unsigned pathID, unsigned readID);

    // Split the path at the given iterator, returns new path ID.
    unsigned splitPathAt(BogPath & path, std::vector<unsigned>::const_iterator it);

    // Update pathForRead for the just-created path with pathID.
    void updatePathForRead(unsigned pathID);

    // Returns true if the path is long.
    inline bool isLong(BogPath const & path);

    // The BogPathSet to modify.
    BogPathSet & pathSet;
    // The containment filter to use for resolving containment overlaps.
    ContainmentFilter const & filter;
    // Configuration for the Unitigging.
    UnitiggerOptions const & options;
    // Enable/disable logging.
    bool logging;
};

void BogPathSplitter::splitFromRight(unsigned pathID, unsigned readID)
{
    auto & path = pathSet.paths.at(pathID);
    if (path.reads.back() == readID)
        return;  // no splitting necessary, at end

    SEQAN_ASSERT(std::find(path.reads.begin(), path.reads.end(), readID) != path.reads.end());
    std::vector<unsigned>::const_iterator it = std::find(path.reads.begin(), path.reads.end(), readID);
    updatePathForRead(splitPathAt(path, std::next(it)));
}

void BogPathSplitter::splitFromLeft(unsigned pathID, unsigned readID)
{
    auto & path = pathSet.paths.at(pathID);
    if (path.reads.front() == readID)
        return;  // no splitting necessary, at beginning

    SEQAN_ASSERT(std::find(path.reads.begin(), path.reads.end(), readID) != path.reads.end());
    std::vector<unsigned>::const_iterator it = std::find(path.reads.begin(), path.reads.end(), readID);
    SEQAN_ASSERT(it != path.reads.begin());
    updatePathForRead(splitPathAt(path, it));
}

unsigned BogPathSplitter::splitPathAt(BogPath & path, std::vector<unsigned>::const_iterator it)
{
    unsigned offset = (it - path.reads.begin());  // number of reads left of it
    SEQAN_ASSERT_GT(offset, 0u);

    if (logging)
        std::cerr << "Splitting BOG path at position " << offset << "\n"
                  << "PATH\n" << path << "\n";

    // Create new path.
    //
    // NB: We have to be careful about the order of operations here since pushing to pathSet.paths invalidates path and
    // it since the path can change places in memory.
    BogPath newPath;
    std::copy(path.reads.begin() + offset, path.reads.end(), std::back_inserter(newPath.reads));
    std::copy(path.readLengths.begin() + offset, path.readLengths.end(), std::back_inserter(newPath.readLengths));
    std::copy(path.overlaps.begin() + offset, path.overlaps.end(), std::back_inserter(newPath.overlaps));
    newPath.updateLength();

    // Update old path.
    path.reads.resize(offset);
    path.readLengths.resize(offset);
    path.overlaps.resize(offset - 1);
    path.updateLength();

    if (logging)
        std::cerr << "UPDATED PATH\n"
                  << path
                  << "NEW PATH (id: " << pathSet.paths.size() << ")\n"
                  << newPath
                  << "\n";

    // Push path to pathSet.
    pathSet.paths.emplace_back(newPath);
    return pathSet.paths.size() - 1;
}

void BogPathSplitter::updatePathForRead(unsigned pathID)
{
    auto const & path = pathSet.paths.at(pathID);
    for (auto readID : path.reads)
        pathSet.pathForRead[readID] = pathID;
}

bool BogPathSplitter::isLong(BogPath const & path)
{
    return (path.length >= options.longBogPathMinLength) && (path.size() >= options.longBogPathMinReadCount);
}

void BogPathSplitter::run()
{
    // We look for paths whose first/last node are node mapped to the same path by pathForRead.  This is caused by the
    // path stopping at a read that is already incorporated in a previous, possibly longer path.  If it turns out that
    // the current path is not long enough for splitting then we remove the leading/trailing node again.  If we can
    // split the other path then shared begin node will also remain with the one part that is created during splitting.
    for (unsigned pathID = 0, numPaths = pathSet.size(); pathID < numPaths; ++pathID)
    {
        auto & path = pathSet.paths.at(pathID);
        SEQAN_ASSERT_NOT(path.reads.empty());
        unsigned firstReadID = path.reads.front();
        unsigned lastReadID = path.reads.back();

        // Get IDs of other path.
        SEQAN_ASSERT(pathSet.pathForRead.count(firstReadID));
        SEQAN_ASSERT(pathSet.pathForRead.count(lastReadID));
        auto otherID0 = pathSet.pathForRead[firstReadID];
        SEQAN_ASSERT(pathSet.paths[otherID0].contains(firstReadID));

        // Path length is decreased by one read after splitting but take length at beginning.
        bool isPathLong = isLong(pathSet.paths[pathID]);

        // If the path is long then perform the splitting of the other path (we have to differentiate between splitting
        // for the beginning and end of path) and remove first/last node.
        //
        // NB: The splitFromRight() and splitFromLeft() calls invalidate pathSet.paths, so path is invalid from here!
        if (otherID0 != pathID)
        {
            SEQAN_ASSERT_MSG(pathSet.paths[otherID0].contains(firstReadID),
                             "otherID0=%u, firstReadID=%u", otherID0, firstReadID);
            if (isPathLong)
            {
                if (logging)
                    std::cerr << "Splitting path " << otherID0 << " from right at id " << firstReadID << "\n"
                              << "USING THE FOLLOWING PATH (#" << pathID << ") FOR SPLITTING\n"
                              << path << "\n";
                splitFromRight(otherID0, firstReadID);
            }
            removeFirstRead(pathSet.paths[pathID]);
        }

        // Get id of second path here since it may be a newly created one and getting it from above would get old,
        // invalid state.
        auto otherID1 = pathSet.pathForRead[lastReadID];
        SEQAN_ASSERT(pathSet.paths[otherID1].contains(lastReadID));

        if (otherID1 != pathID && firstReadID != lastReadID)  // ignore in case one-read paths
        {
            SEQAN_ASSERT_MSG(pathSet.paths[otherID1].contains(lastReadID),
                             "otherID1=%u, lastReadID=%u", otherID1, lastReadID);
            if (isPathLong)
            {
                if (logging)
                    std::cerr << "Splitting path " << otherID1 << " from left at id " << lastReadID << "\n"
                              << "USING THE FOLLOWING PATH (#" << pathID << ") FOR SPLITTING\n"
                              << path << "\n";
                splitFromLeft(otherID1, lastReadID);
            }
            removeLastRead(pathSet.paths[pathID]);
        }
    }

    // TODO(holtgrew): Remove empty paths!
    // TODO(holtgrew): Register overlap between reads that started a split.
}

}  // anonymous namespace


// ----------------------------------------------------------------------------
// Class RemovedOverlapStore
// ----------------------------------------------------------------------------

void RemovedOverlapStore::print(std::ostream & out) const
{
    char const * label[3] = { "INVALID", "LEFT", "RIGHT" };

    out << "RemovedOverlapStore\n";
    for (auto const & entry : entries)
        out << "\tRemovedOverlapStore::Entry(readID=" << entry.readID << ", direction=" << label[entry.direction]
            << ", info=" << entry.info << ")\n";
}

// ----------------------------------------------------------------------------
// Function operator<<()
// ----------------------------------------------------------------------------

void print(std::ostream & out, BestOverlapGraph const & g)
{
    out << "Best Overlap Graph\n\n";
    for (unsigned i = 0; i < g.bestLeft.size(); ++i)
        std::cerr << i << "\t" << g.bestLeft[i] << "\t" << g.bestRight[i] << "\n";
}

// ----------------------------------------------------------------------------
// Class ContainmentFilterImpl
// ----------------------------------------------------------------------------

class ContainmentFilterImpl
{
public:
    ContainmentFilterImpl(bool logging) :
            logging(logging), onContig([](unsigned /*readID*/) { return false; })
    {}

    ContainmentFilterImpl(bool logging, std::function<bool(unsigned readID)> onContig) :
            logging(logging), onContig(onContig)
    {}

    std::vector<Overlap> run(size_t numReads, std::vector<Overlap> const & in);
    std::vector<unsigned> const & containedIDs() const;
    unsigned container(unsigned readID) const;

    Overlap overlap(unsigned readID) const;

    bool isContained(unsigned readID) const
    {
        return containers_.count(readID);
    }

private:
    bool logging;
    // Callback that allows to ignore overlaps with reads that are on contigs.
    std::function<bool(unsigned readID)> onContig;

    // Returns true if contained and stacking with other read.
    bool isStackContained(unsigned readID) const
    {
        if (!isContained(readID))
            return false;
        auto ovl = overlap(readID);
        return (ovl.begin0 == 0 && ovl.begin1 == 0 && ovl.len0 == ovl.len1);
    }

    // Collect start nodes for DFS in containment computation.  The result will be the nodes in a topological sorting
    // (ignoring cyclic links obtained in between).
    std::vector<lemon::SmartDigraph::Node> collectStartNodes(lemon::SmartDigraph const & graph) const
    {
        std::vector<lemon::SmartDigraph::Node> result;
        result.reserve(countNodes(graph));

        enum Mark { UNMARKED, TEMPORARY, PERMANENT };
        lemon::SmartDigraph::NodeMap<Mark> mark(graph, UNMARKED);

        std::function<void(lemon::SmartDigraph::Node)> visit;
        visit = [&](lemon::SmartDigraph::Node u) {
            if (mark[u] == TEMPORARY)
                return;  // traversed non-DAG edge
            if (mark[u] == PERMANENT)
                return;  // traversed into already processed part of graph
            mark[u] = TEMPORARY;
            for (lemon::SmartDigraph::OutArcIt arc(graph, u); arc != lemon::INVALID; ++arc)
            {
                SEQAN_ASSERT(graph.target(arc) != u);
                visit(graph.target(arc));
            }
            mark[u] = PERMANENT;
            result.push_back(u);
        };

        for (lemon::SmartDigraph::NodeIt u(graph); u != lemon::INVALID; ++u)
            if (mark[u] == UNMARKED)
                visit(u);

        std::reverse(result.begin(), result.end());
        return result;
    }

    // The read contained readIDs.
    std::vector<unsigned> containedIDs_;
    // Overlaps, indexed by contained readID (contained => (container, overlap)).
    std::map<unsigned, std::pair<unsigned, Overlap>> containers_;
};

std::vector<Overlap> ContainmentFilterImpl::run(size_t numReads, std::vector<Overlap> const & in)
{
    std::vector<Overlap> result(in);
    (void)numReads;  // TODO(holtgrew): Remove?

    // Assign reads to container in FCFS mode until no more assignments can be performed.
    std::map<unsigned, std::pair<unsigned, Overlap>> mapping;  // mapping done in this round
    std::set<unsigned> newContainers;  // new containers in this round
    do
    {
        mapping.clear();
        newContainers.clear();

        // Collect mapping.
        //
        // Note that overlaps are assumed to be normalized, seq1 is the contained.
        for (auto const & ovl : result)
            if (isContainment(ovl) && !mapping.count(ovl.seq0) && !mapping.count(ovl.seq1))
            {
                if (newContainers.count(ovl.seq1))
                    continue;  // Marked as contained itself in this round, skip.
                newContainers.insert(ovl.seq0);
                if (logging)
                    std::cerr << "mapping[" << ovl.seq1 << "] = (" << ovl.seq0 << ", " << ovl << ")\n";
                mapping[ovl.seq1] = std::make_pair(ovl.seq0, ovl);
            }

        // If we are mapping a container then update this to a mapping and overlap to the new container.
        if (!mapping.empty())
            for (auto & containerPair : containers_)
            {
                auto containedID = containerPair.first;
                (void)containedID;  // only used in debug mode
                unsigned containerID = containerPair.second.first;
                auto const & ovl = containerPair.second.second;
                (void)ovl;  // only used for assertion
                SEQAN_ASSERT_EQ(ovl.seq0, containerID);
                SEQAN_ASSERT_EQ(ovl.seq1, containedID);
                if (mapping.count(containerID))
                {
                    // if (logging)
                    //     std::cerr << "Updating from " << ovl << "\n";
                    unsigned newContainerID = mapping[containerID].first;
                    // if (logging)
                    //     std::cerr << "\tUpdating with " << mapping[containerID].second << "\n";
                    containerPair.second.first = newContainerID;
                    containerPair.second.second.seq0 = newContainerID;
                    containerPair.second.second.begin1 += mapping[containerID].second.begin1;
                    // if (logging)
                    //     std::cerr << "\tUpdating to " << ovl << "\n";
                }
            }

        // Enter new mappings into containedIDs and containers.
        for (auto const & pair : mapping)
        {
            containedIDs_.push_back(pair.first);
            if (logging)
                std::cerr << "container_.insert(" << pair.first << ", " << pair.second.first << ", "
                          << pair.second.second << ")\n";
            SEQAN_ASSERT_NOT(containers_.count(pair.first));
            containers_.insert(pair);
        }

        // Re-sort containedIDs_.
        std::sort(containedIDs_.begin(), containedIDs_.end());
        containedIDs_.resize(std::unique(containedIDs_.begin(), containedIDs_.end()) - containedIDs_.begin());

        // Translate overlaps replacing ids of contained reads with id of their container.
        std::for_each(result.begin(), result.end(),
                       [&](Overlap & ovl) {
                           bool contained0 = isContained(ovl.seq0);
                           bool contained1 = isContained(ovl.seq1);
                           // if (logging && (contained0 || contained1))
                           //     std::cerr << "Updating from " << ovl << "\n";
                           if (contained0)
                           {
                               auto it = containers_.find(ovl.seq0);
                               ovl.seq0 = it->second.first;
                               ovl.begin1 += it->second.second.begin1;
                           }
                           if (contained1)
                           {
                               auto it = containers_.find(ovl.seq1);
                               ovl.seq1 = it->second.first;
                               ovl.begin1 += it->second.second.begin1;
                           }
                           // if (logging && (contained0 || contained1))
                           //     std::cerr << "Updating to " << ovl << "\n";
                           SEQAN_ASSERT_NOT_MSG(isContained(ovl.seq0), "ovl.seq0 = %u", ovl.seq0);
                           SEQAN_ASSERT_NOT_MSG(isContained(ovl.seq1), "ovl.seq1 = %u", ovl.seq1);
                       });

        // Remove self-links.
        eraseIf(result, [](Overlap const & ovl) { return (ovl.seq0 == ovl.seq1); });

        // Sort result by (seq0, seq1, begin1, errors) such that std::unique() gives us the best overlap.
        auto lt = [](Overlap const & lhs, Overlap const & rhs) {
            return (std::make_tuple(lhs.seq0, lhs.seq1, lhs.begin1, lhs.errors) <
                    std::make_tuple(rhs.seq0, rhs.seq1, rhs.begin1, rhs.errors));
        };
        std::sort(result.begin(), result.end(), lt);
        auto eq = [](Overlap const & lhs, Overlap const & rhs) {
            return (std::make_pair(lhs.seq0, lhs.seq1) == std::make_pair(lhs.seq0, rhs.seq1));
        };
        auto itEnd = std::unique(result.begin(), result.end(), eq);
        result.resize(itEnd - result.begin());
    }
    while (!mapping.empty());

    return result;
}

std::vector<unsigned> const & ContainmentFilterImpl::containedIDs() const
{
    return containedIDs_;
}

unsigned ContainmentFilterImpl::container(unsigned readID) const
{
    SEQAN_ASSERT(containers_.count(readID));
    return containers_.find(readID)->second.first;
}

Overlap ContainmentFilterImpl::overlap(unsigned readID) const
{
    SEQAN_ASSERT(containers_.count(readID));
    return containers_.find(readID)->second.second;
}

// --------------------------------------------------------------------------
// Class ContainmentFilter
// --------------------------------------------------------------------------

ContainmentFilter::ContainmentFilter(bool logging) :
        impl(new ContainmentFilterImpl(logging))
{}

ContainmentFilter::~ContainmentFilter()
{}

std::vector<Overlap> ContainmentFilter::run(size_t numReads, std::vector<Overlap> const & in)
{
    return impl->run(numReads, in);
}

std::vector<unsigned> const & ContainmentFilter::containedIDs() const
{
    return impl->containedIDs();
}

unsigned ContainmentFilter::container(unsigned readID) const
{
    return impl->container(readID);
}

bool ContainmentFilter::isContained(unsigned readID) const
{
    return impl->isContained(readID);
}

Overlap ContainmentFilter::overlap(unsigned readID) const
{
    return impl->overlap(readID);
}

// ----------------------------------------------------------------------------
// Class BogPathLengths
// ----------------------------------------------------------------------------

unsigned const BogPathLengths::UNSET_PATH_LENGTH = -1;

void BogPathLengths::reinit(size_t numReads)
{
    left.clear();
    left.resize(numReads, UNSET_PATH_LENGTH);
    right.clear();
    right.resize(numReads, UNSET_PATH_LENGTH);
}

// ----------------------------------------------------------------------------
// Function operator<<()                        [BestOverlapGraph::OverlapInfo]
// ----------------------------------------------------------------------------

std::ostream & operator<<(std::ostream & out, BestOverlapGraph::OverlapInfo info)
{
    return out << "BestOverlapGraph::OverlapInfo(" << info.readID << ", " << info.length
               << ", " << info.score << ", " << info.begin1 << ")";
}

// --------------------------------------------------------------------------
// Function buildBestOverlapGraph()
// --------------------------------------------------------------------------

BestOverlapGraph buildBestOverlapGraph(size_t numReads,
                                       std::vector<Overlap> const & overlaps,
                                       bool logging)
{
    BestOverlapGraph result(numReads);

    for (auto const & ovl : overlaps)
    {
        SEQAN_ASSERT_EQ(ovl.begin0, 0u);

        auto leftInfo = infoFromOverlap(ovl, /*left=*/true);
        if (!result.bestLeft.at(ovl.seq1).valid() || result.bestLeft[ovl.seq1] < leftInfo)
            result.bestLeft[ovl.seq1] = leftInfo;
        auto rightInfo = infoFromOverlap(ovl, /*left=*/false);
        if (!result.bestRight.at(ovl.seq0).valid() || result.bestRight[ovl.seq0] < rightInfo)
            result.bestRight[ovl.seq0] = rightInfo;
    }

    if (logging)
        print(std::cerr, result);

    return result;
}

// ----------------------------------------------------------------------------
// Function printDot()
// ----------------------------------------------------------------------------

void printDot(std::ostream & out, BestOverlapGraph const & bog)
{
    out << "digraph BOG {\n"
        << "/* nodes */\n";

    for (unsigned i = 0; i < bog.bestLeft.size(); ++i)
        out << i << ";\n";

    out << "\n/* arcs */\n";
    for (unsigned i = 0; i < bog.bestLeft.size(); ++i)
        if (bog.bestLeft[i].valid())
            out << i << " -> " << bog.bestLeft[i].readID << " [color=blue,label=L];\n";
    for (unsigned i = 0; i < bog.bestRight.size(); ++i)
        if (bog.bestRight[i].valid())
            out << i << " -> " << bog.bestRight[i].readID << " [color=red,label=R];\n";

    out << "}\n";
}

// ----------------------------------------------------------------------------
// Function removeCycles()
// ----------------------------------------------------------------------------

void removeCycles(RemovedOverlapStore & removedOvls, BestOverlapGraph & og)
{
    BogCycleRemover remover(removedOvls, og);
    remover.run();
}

// ----------------------------------------------------------------------------
// Function computePathLengths()
// ----------------------------------------------------------------------------

BogPathLengths computePathLengths(BestOverlapGraph const & bog)
{
    BogPathLengths pathLengths;
    pathLengths.reinit(bog.numReads());

    // TODO(holtgrew): Shuffle start vertices?

    for (size_t readID = 0; readID < bog.numReads(); ++readID)
    {
        size_t countRight = 0;
        for (unsigned right = bog.bestRight.at(readID).readID;
             right != (unsigned)-1;
             right = bog.bestRight.at(right).readID)
            if (pathLengths.right[right] != BogPathLengths::UNSET_PATH_LENGTH)
            {
                // std::cerr << "right=" << right << "\tcountRight=" << countRight << "\n";
                countRight += pathLengths.right[right] + 1;
                break;
            }
            else
            {
                // std::cerr << "right=" << right << "\tcountRight=" << countRight << "\n";
                ++countRight;
            }
        pathLengths.right[readID] = countRight;

        size_t countLeft = 0;
        for (unsigned left = bog.bestLeft.at(readID).readID;
             left != (unsigned)-1;
             left = bog.bestLeft.at(left).readID)
            if (pathLengths.left[left] != BogPathLengths::UNSET_PATH_LENGTH)
            {
                countLeft += pathLengths.left[left] + 1;
                break;
            }
            else
            {
                ++countLeft;
            }
        pathLengths.left[readID] = countLeft;
    }

    return pathLengths;
}

// ----------------------------------------------------------------------------
// Function buildInitialPaths()
// ----------------------------------------------------------------------------

BogPathSet buildInitialPaths(std::vector<unsigned> const & readLengths,
                             BestOverlapGraph const & bog,
                             BogPathLengths const & pathLengths,
                             bool logging)
{
    InitialBogPathBuilder builder(readLengths, bog, pathLengths, logging);
    return builder.run();
}

// ----------------------------------------------------------------------------
// Class BogPath
// ----------------------------------------------------------------------------

void BogPath::updateLength()
{
    length = 0;
    if (reads.empty())
        return;
    for (auto const & ovl : overlaps)
        length += ovl.begin1;
    length += readLengths.back();
}

std::ostream & operator<<(std::ostream & out, BogPath::OverlapInfo const & info)
{
    return out << "BogPath::OverlapInfo(id0=" << info.id0 << ", id1=" << info.id1
               << ", begin1=" << info.begin1 << ")";
}

std::ostream & operator<<(std::ostream & out, BogPath const & bogPath)
{
    out << "BOG Path\n"
        << "reads:";
    for (auto readID : bogPath.reads)
        out << " " << readID;
    out << "\nlengths:";
    for (auto len : bogPath.readLengths)
        out << " " << len;
    out << "\n";
    out << "overlaps: ";
    std::copy(bogPath.overlaps.begin(), bogPath.overlaps.end(),
              std::ostream_iterator<BogPath::OverlapInfo>(std::cerr, " "));
    out << "\n";

    return out;
}

// ----------------------------------------------------------------------------
// Class BogPathSet
// ----------------------------------------------------------------------------

unsigned BogPathSet::registerPath(BogPath && path)
{
    for (auto readID : path.reads)
        if (!pathForRead.count(readID))
            pathForRead[readID] = paths.size();
    paths.emplace_back(path);
    return paths.size() - 1;
}

// ----------------------------------------------------------------------------
// Function performPathSplitting()
// ----------------------------------------------------------------------------

void performPathSplitting(BogPathSet & pathSet,
                          ContainmentFilter const & filter,
                          UnitiggerOptions const & options,
                          bool logging)
{
    BogPathSplitter splitter(pathSet, filter, options, logging);
    splitter.run();
}

// ----------------------------------------------------------------------------
// Function completePaths()
// ----------------------------------------------------------------------------

void completePaths(BogPathSet & paths,
                   RemovedOverlapStore const & ovls,
                   std::vector<unsigned> const & readLengths)
{
    for (auto const & entry : ovls)
    {
        SEQAN_CHECK(paths.pathForRead.count(entry.readID), "Must exist in set.");
        unsigned pathID = paths.pathForRead.find(entry.readID)->second;
        auto & path = paths.paths[pathID];
        unsigned readLen = readLengths.at(entry.info.readID);
        if (std::find(path.reads.begin(), path.reads.end(), entry.readID) != path.reads.end())
            continue;  // Skip, already on path
        if (entry.direction == RemovedOverlapStore::RIGHT && path.reads.back() == entry.readID)
        {
            path.readLengths.push_back(readLen);
            path.reads.push_back(entry.info.readID);
            BogPath::OverlapInfo ovl(entry.readID, entry.info.readID, entry.info.begin1);
            path.overlaps.push_back(ovl);
        }
        else if (entry.direction == RemovedOverlapStore::LEFT && path.reads.front() == entry.readID)
        {
            path.readLengths.insert(path.readLengths.begin(), readLen);
            path.reads.insert(path.reads.begin(), entry.info.readID);
            BogPath::OverlapInfo ovl(entry.info.readID, entry.readID, entry.info.begin1);
            path.overlaps.insert(path.overlaps.begin(), ovl);
        }
    }
}

}  // namespace assembler
