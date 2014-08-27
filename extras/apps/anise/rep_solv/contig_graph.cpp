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

#include "contig_graph.h"

#include <iterator>

#include <lemon/dfs.h>

#include <seqan/sequence.h>

#include "asm/overlapper.h"
#include "asm/contig_graph.h"

#include "scaffolder/gpm.h"
#include "scaffolder/gpm_options.h"
#include "scaffolder/mate_link.h"

namespace rep_solv {

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
// Class ContigGraphBuilder
// ----------------------------------------------------------------------------

class ContigGraphBuilder
{
public:
    ContigGraphBuilder(ContigGraph & out,
                       assembler::ContigGraph const & in,
                       TFragmentStore const & fragStore) :
            out(out), in(in), fragStore(fragStore)
    {}

    void run()
    {
        copyNodes();
        copyOverlapLinks();
    }

private:
    // Copy over the nodes (adds nodes to graph, contig labels to contig, and nodes to node).
    void copyNodes();
    // Copy over the links from assembler ContigGraph (represents overlaps) to repeat resolution contig graph.
    void copyOverlapLinks();

    // The resulting contig graph is written here.
    ContigGraph & out;
    // The input contig graph is available here.
    assembler::ContigGraph const & in;
    // FragmentStore.
    TFragmentStore const & fragStore;
};

void ContigGraphBuilder::copyNodes()
{
    out.node.resize(in.node.size());

    // Compute mate indices and compute number of pairs for each contig.
    seqan::String<unsigned> mateIdx;
    calculateMateIndices(mateIdx, fragStore);
    std::vector<unsigned> numPairs(length(fragStore.contigStore), 0u);
    for (auto const & el : fragStore.alignedReadStore)
        if (mateIdx[el.readId] != TAlignedReadStoreElement::INVALID_ID)
        {
            auto const & otherEl = fragStore.alignedReadStore[mateIdx[el.readId]];
            if (el.contigId == otherEl.contigId && el.readId < otherEl.readId)
                numPairs[el.contigId] += 1;
        }

    // Copy over nodes and assign labels.
    for (auto u : in.node)
    {
        auto n = out.graph.addNode();
        auto contigID = in.contig[u].id;
        out.node[contigID] = n;
        unsigned contigLen = length(fragStore.contigStore[contigID].seq);
        out.contig[n] = ContigLabel(contigID, contigLen, numPairs.at(contigID));
    }

    // Copy over read to contig mapping.
    for (auto const & pair : in.readToContig)
    {
        SEQAN_ASSERT_NOT(out.readToContigs.count(pair.first));
        out.readToContigs[pair.first] = { in.contig[pair.second].id };
    }
}

void ContigGraphBuilder::copyOverlapLinks()
{
    for (lemon::SmartDigraph::ArcIt arc(in.graph); arc != lemon::INVALID; ++arc)
    {
        auto left = in.graph.source(arc);
        auto right = in.graph.target(arc);
        auto leftID = in.contig[left].id;
        auto rightID = in.contig[right].id;

        auto e = out.graph.addEdge(out.node[leftID], out.node[rightID]);
        auto & label = out.link[e];
        label = EdgeLabel(leftID, rightID);
        label.overlapDistance.mean = -in.overlapLength[arc];
    }
}

// ----------------------------------------------------------------------------
// Class Dagifier
// ----------------------------------------------------------------------------

// Helper class that explores the ContigGraph from s and t (towards the left and right, respectively) and marks non-DAG
// edges with false flag.
//
// Also flags unreachable edges as non-DAG edges.

class Dagifier
{
private:
    enum Marker { UNREACHED, TEMPORARY, PERMANENT };

public:
    Dagifier(lemon::SmartGraph::EdgeMap<bool> & dagEdge,
             ContigGraph const & in,
             lemon::SmartGraph::NodeMap<bool> const & doRemove) :
            dagEdge(dagEdge), graph(in), doRemove(doRemove), usedEdge(in.graph, false)
    {}

    void run();

private:

    // Start DFS from graph.s to the right.
    void dfsRight(lemon::SmartGraph::NodeMap<Marker> & reached);
    // Recursion helper for DFS from graph.s to the right.
    void dfsRightRec(lemon::SmartGraph::NodeMap<Marker> & reached, lemon::SmartGraph::Node u);
    // Start DFS from graph.t to the left.
    void dfsLeft(lemon::SmartGraph::NodeMap<Marker> & reached);
    // Recursion helper for DFS from graph.t to the left.
    void dfsLeftRec(lemon::SmartGraph::NodeMap<Marker> & reached, lemon::SmartGraph::Node u);

    // Input / Output

    // Edge map for flagging edges as not being part of the DAGified version.
    lemon::SmartGraph::EdgeMap<bool> & dagEdge;
    // Input graph.
    ContigGraph const & graph;
    // Flags for whether or not to a node is flagged for removal.
    lemon::SmartGraph::NodeMap<bool> const & doRemove;

    // Whether or not an edge is traversed in the DFS.
    lemon::SmartGraph::EdgeMap<bool> usedEdge;
};

void Dagifier::run()
{
    lemon::SmartGraph::NodeMap<Marker> markR(graph.graph, UNREACHED), markL(graph.graph, UNREACHED);
    dfsRight(markR);
    dfsLeft(markL);

    // Mark unreached edges as non-DAG edge.
    for (lemon::SmartGraph::EdgeIt edge(graph.graph); edge != lemon::INVALID; ++edge)
        if (!usedEdge[edge])
            dagEdge[edge] = false;
}

void Dagifier::dfsRight(lemon::SmartGraph::NodeMap<Dagifier::Marker> & mark)
{
    mark[graph.s] = TEMPORARY;

    for (lemon::SmartGraph::OutArcIt arc(graph.graph, graph.s); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(graph.s));
        if (graph.link[arc].leftID != graph.contig[graph.s].id)
            continue;  // skip, link has wrong direction
        if (mark[graph.graph.target(arc)] == TEMPORARY)
        {
            // std::cerr << "==== NOT DAG EDGE\t" << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
            dagEdge[arc] = false;  // can read node marked as temporary over this edge, mark for removal
        }
        else
        {
            // std::cerr << "==== USING\t" << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
            usedEdge[arc] = true;
            dfsRightRec(mark, graph.graph.target(arc));
        }
    }

    mark[graph.s] = PERMANENT;
}

void Dagifier::dfsRightRec(lemon::SmartGraph::NodeMap<Dagifier::Marker> & mark, lemon::SmartGraph::Node u)
{
    if (doRemove[u])
        return;  // break recursion if to be removed, link does not exist.
    if (u == graph.t || mark[u] != UNREACHED)
        return;  // break recursion at s and marked nodes
    mark[u] = TEMPORARY;

    // Iterate over left-to-right edges.
    for (lemon::SmartGraph::OutArcIt arc(graph.graph, u); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(u));
        if (graph.link[arc].leftID != graph.contig[u].id)
            continue;  // skip, link has wrong direction
        // std::cerr << "traversing " << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
        if (mark[graph.graph.target(arc)] == TEMPORARY)
        {
            // std::cerr << "==== NOT DAG EDGE\t" << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
            dagEdge[arc] = false;  // can read node marked as temporary over this edge, mark for removal
        }
        else
        {
            // std::cerr << "==== USING\t" << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
            usedEdge[arc] = true;
            dfsRightRec(mark, graph.graph.target(arc));
        }
    }

    mark[u] = PERMANENT;
}

void Dagifier::dfsLeft(lemon::SmartGraph::NodeMap<Dagifier::Marker> & mark)
{
    mark[graph.t] = TEMPORARY;

    for (lemon::SmartGraph::OutArcIt arc(graph.graph, graph.t); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(graph.t));
        if (graph.link[arc].rightID != graph.contig[graph.t].id)
            continue;  // skip, link has wrong direction
        if (mark[graph.graph.target(arc)] == TEMPORARY)
        {
            // std::cerr << "==== NOT DAG EDGE\t" << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
            dagEdge[arc] = false;  // can read node marked as temporary over this edge, mark for removal
        }
        else
        {
            // std::cerr << "==== USING\t" << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
            usedEdge[arc] = true;
            dfsLeftRec(mark, graph.graph.target(arc));
        }
    }

    mark[graph.t] = PERMANENT;
}

void Dagifier::dfsLeftRec(lemon::SmartGraph::NodeMap<Dagifier::Marker> & mark, lemon::SmartGraph::Node u)
{
    if (doRemove[u])
        return;  // break recursion if to be removed, link does not exist.
    if (u == graph.t || mark[u] != UNREACHED)
        return;  // break recursion at t and mark nodes
    mark[u] = TEMPORARY;

    // Iterate over right-to-left edges.
    for (lemon::SmartGraph::OutArcIt arc(graph.graph, u); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(u));
        if (graph.link[arc].rightID != graph.contig[u].id)
            continue;  // skip, link has wrong direction
        if (mark[graph.graph.target(arc)] == TEMPORARY)
        {
            // std::cerr << "==== NOT DAG EDGE\t" << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
            dagEdge[arc] = false;  // can read node marked as temporary over this edge, mark for removal
        }
        else
        {
            // std::cerr << "==== USING\t" << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
            usedEdge[arc] = true;
            dfsLeftRec(mark, graph.graph.target(arc));
        }
    }

    mark[u] = PERMANENT;
}

// ----------------------------------------------------------------------------
// Class ContigRemover
// ----------------------------------------------------------------------------

// Helper for removeContigs().

class ContigRemover
{
public:
    ContigRemover(ContigGraph & out, ContigGraph const & in, lemon::SmartGraph::NodeMap<bool> const & doRemove) :
            out(out), in(in), doRemove(doRemove), dagEdge(in.graph, true), newToOld(buildNewToOld()),
            oldToNew(buildOldToNew())
    {}

    void run()
    {
        // std::cerr << "ContigRemover Input\n";
        // in.print(std::cerr);

        // Compute DAG edges.
        //
        // This is required since adding paired end links might lead to DAG edges.
        Dagifier helper(dagEdge, in, doRemove);
        helper.run();

        copyContigs();
        copyLinks();
        copyReadInfo();
    }

private:

    // Marker for invalid values.
    static const unsigned INVALID;

    // Build mapping from new contig ID to old contig ID.
    std::vector<unsigned> buildNewToOld() const;
    // Build mapping from old contig ID to new contig ID, INVALID if removed.
    std::vector<unsigned> buildOldToNew() const;

    // Copy out contig.
    void copyContigs();
    // Copy out links.
    void copyLinks();
    // Copy out read info (member readToContigs).
    void copyReadInfo();

    // Input / Output

    ContigGraph & out;
    ContigGraph const & in;
    lemon::SmartGraph::NodeMap<bool> const & doRemove;

    // State
    
    // Flag for whether the an edge in in remains part of the output graph (used in Dagifier).
    lemon::SmartGraph::EdgeMap<bool> dagEdge;

    // Mapping new contigID to old contigID and vice versa.
    std::vector<unsigned> newToOld, oldToNew;
};

const unsigned ContigRemover::INVALID = (unsigned)-1;

std::vector<unsigned> ContigRemover::buildNewToOld() const
{
    std::vector<unsigned> result;

    for (unsigned oldID = 0; oldID < in.node.size(); ++oldID)  // old contig ID
    {
        auto u = in.node[oldID];
        if (!doRemove[u])
            result.push_back(in.contig[u].id);
    }

    return result;
}

std::vector<unsigned> ContigRemover::buildOldToNew() const
{
    std::vector<unsigned> result(in.node.size(), INVALID);

    unsigned newID = 0;
    for (unsigned oldID = 0; oldID < in.node.size(); ++oldID)  // old contig ID
    {
        auto u = in.node[oldID];
        if (!doRemove[u])
        {
            // std::cerr << "ContigRemover\toldToNew[" << in.contig[u].id << "] = " << newID << "\n";
            result.at(in.contig[u].id) = newID++;
        }
    }

    return result;
}

void ContigRemover::copyContigs()
{
    SEQAN_ASSERT_NOT(oldToNew.empty());
    SEQAN_ASSERT_NOT(newToOld.empty());

    for (unsigned contigID = 0; contigID < newToOld.size(); ++contigID)  // new contig ID
    {
        SEQAN_ASSERT_EQ(out.node.size(), contigID);
        auto u = out.graph.addNode();
        out.node.push_back(u);
        out.contig[u] = in.contig[in.node[newToOld.at(contigID)]];
        SEQAN_ASSERT_EQ(newToOld.at(contigID), out.contig[u].id);
        out.contig[u].id = contigID;
    }
}

void ContigRemover::copyLinks()
{
    for (lemon::SmartGraph::EdgeIt edge(in.graph); edge != lemon::INVALID; ++edge)
    {
        if (!dagEdge[edge])
            continue;  // Skip if not part of DAG edge.

        // Get shortcut to label from input.
        auto const & label = in.link[edge];
        SEQAN_ASSERT_NEQ(label.rightID, +EdgeLabel::SOURCE);
        SEQAN_ASSERT_NEQ(label.leftID, +EdgeLabel::TARGET);

        // std::cerr << "\tContigRemover::copyLinks()\tConsidering\t" << label.leftID
        //           << " -- " << label.rightID << "\n";

        if (label.leftID == EdgeLabel::SOURCE)
        {
            if (oldToNew.at(label.rightID) == INVALID)
                continue;  // skip, other vertex removed
            auto e = out.graph.addEdge(out.s, out.node.at(oldToNew.at(label.rightID)));
            out.link[e] = label;
            out.link[e].rightID = oldToNew[label.rightID];
            continue;
        }
        else if (label.rightID == EdgeLabel::TARGET)
        {
            if (oldToNew.at(label.leftID) == INVALID)
                continue;  // skip, other vertex removed
            auto e = out.graph.addEdge(out.node.at(oldToNew.at(label.leftID)), out.t);
            out.link[e] = label;
            out.link[e].leftID = oldToNew[label.leftID];
            continue;
        }

        // std::cerr << "\tContigRemover::copyLinks()\tno special case\n";

        if (oldToNew.at(label.leftID) == INVALID || oldToNew.at(label.rightID) == INVALID)
            continue;  // skip, one or both vertices are removed

        // Add edge and copy link label.
        auto e = out.graph.addEdge(out.node.at(oldToNew[label.leftID]), out.node.at(oldToNew[label.rightID]));
        out.link[e] = label;
        // std::cerr << "\tContigRemover::copyLinks() => ADDING LINK\t" << out.link[e].leftID
        //           << " -- " << out.link[e].rightID << "\n";
        // Modify the link label.
        if (out.link[e].leftID != EdgeLabel::SOURCE && out.link[e].leftID != EdgeLabel::TARGET)
            out.link[e].leftID = oldToNew.at(out.link[e].leftID);
        if (out.link[e].rightID != EdgeLabel::SOURCE && out.link[e].rightID != EdgeLabel::TARGET)
            out.link[e].rightID = oldToNew.at(out.link[e].rightID);
    }
}

void ContigRemover::copyReadInfo()
{
    for (auto const & pair : in.readToContigs)
        for (auto const contigID : pair.second)  // old contig ID
            if (oldToNew.at(contigID) != INVALID)  // valid mapping
                out.readToContigs[/*readID=*/pair.first].push_back(oldToNew[contigID]);
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class ContigGraph
// ----------------------------------------------------------------------------

ContigGraph::ContigGraph(ContigGraph const & other) :
        contig(graph), link(graph), readToContigs(other.readToContigs)
{
    // D'oh, copy constructors are deleted for LEMON graphs, so we have to copy and reconstruct on our own.

    // Copy over the graph manually.
    lemon::SmartGraph::NodeMap<lemon::SmartGraph::Node> nr(other.graph);
    lemon::SmartGraph::EdgeMap<lemon::SmartGraph::Edge> er(other.graph);
    graphCopy(other.graph, graph).nodeRef(nr).edgeRef(er).run();

    // Copy s, t.
    s = nr[other.s];
    t = nr[other.t];
    contig[s].id = EdgeLabel::SOURCE;
    contig[t].id = EdgeLabel::TARGET;

    // Copy contig and link labels.
    for (lemon::SmartGraph::NodeIt u(other.graph); u != lemon::INVALID; ++u)
        contig[nr[u]] = other.contig[u];
    for (lemon::SmartGraph::EdgeIt e(other.graph); e != lemon::INVALID; ++e)
        link[er[e]] = other.link[e];

    // Copy over node mapping.
    node.resize(other.node.size());
    for (unsigned contigID = 0; contigID < other.node.size(); ++contigID)
        node[contigID] = nr[other.node[contigID]];
}

void ContigGraph::print(std::ostream & out) const
{
    out << "ContigGraph\n"
        << "NODES\n";
    for (lemon::SmartGraph::NodeIt u(graph); u != lemon::INVALID; ++u)
        out << "\t" << contig[u].id << " (id=" << contig[u].id << ", len=" << contig[u].length << ")\n";
    out << "\nEDGES\n";
    for (lemon::SmartGraph::EdgeIt edge(graph); edge != lemon::INVALID; ++edge)
        out << "\t" << link[edge].leftID << " -- " << link[edge].rightID << "\n";
    // out << "READ TO CONTIGS\n";
    // for (auto const & pair : readToContigs)
    // {
    //     out << pair.first << "\t->\t";
    //     std::copy(pair.second.begin(), pair.second.end(),
    //               std::ostream_iterator<unsigned>(out, " "));
    //     out << "\n";
    // }
}

// ----------------------------------------------------------------------------
// Class MateInfo
// ----------------------------------------------------------------------------

std::ostream & operator<<(std::ostream & out, MateInfo const & info)
{
    return out << "MateInfo(readL=" << info.readL << ", readR=" << info.readR << ", leftID=" << info.leftID
               << ", rightID=" << info.rightID << ", beginPos=" << info.beginPos
               << ", endPos=" << info.endPos << ", leftLen=" << info.leftLen
               << ", rightLen=" << info.rightLen << ", libraryID=" << info.libraryID
               << ", numCopies=" << info.numCopies << ")";
}

bool lt(MateInfo const & lhs, MateInfo const & rhs)
{
    return (std::make_tuple(lhs.readL, lhs.readR) < std::make_tuple(rhs.readL, rhs.readR));
}

// ----------------------------------------------------------------------------
// Class MateInfos
// ----------------------------------------------------------------------------

void MateInfos::print(std::ostream & out)
{
    out << "MateInfos\n";
    std::copy(records.begin(), records.end(), std::ostream_iterator<MateInfo>(out, "\n"));
}

void MateInfos::insert(MateInfo const & info)
{
    readToRecords[info.readL].push_back(records.size());
    readToRecords[info.readR].push_back(records.size());
    records.push_back(info);
}

// ----------------------------------------------------------------------------
// Function buildContigGraph()
// ----------------------------------------------------------------------------

void buildContigGraph(ContigGraph & out,
                      assembler::ContigGraph const & in,
                      TFragmentStore const & fragStore)
{
    ContigGraphBuilder builder(out, in, fragStore);
    builder.run();
}

// ----------------------------------------------------------------------------
// Function removeContigs()
// ----------------------------------------------------------------------------

// Remove the contigs that are flagged so in doRemove.
//
// Note that this can yield unmapped reads (marked by out.readToContigs[readID].empty()).

void removeContigs(ContigGraph & out,
                   ContigGraph const & in,
                   lemon::SmartGraph::NodeMap<bool> const & doRemove)
{
    ContigRemover helper(out, in, doRemove);
    helper.run();
}

// ----------------------------------------------------------------------------
// Function stPathExists()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Limit traversed edges to left-to-right?

bool stPathExists(ContigGraph const & in)
{
    lemon::Dfs<lemon::SmartGraph> dfs(in.graph);
    dfs.run(in.s, in.t);
    return dfs.reached(in.t);
}

// ----------------------------------------------------------------------------
// Function addPairedEndLinks()
// ----------------------------------------------------------------------------

namespace {  // anonymous namespace

class PairedEndLinkAdder
{
public:
    PairedEndLinkAdder(ContigGraph & cg,
                       MateInfos const & mateInfos,
                       seqan::StringSet<seqan::Dna5String> const * contigSeqs = nullptr,
                       bool logging = false) :
            cg(cg), mateInfos(mateInfos), contigSeqs(contigSeqs), logging(logging)
    {}

    void run()
    {
        if (logging)
        {
            std::cerr << "PAIRED-LINK ADDING FOR\n";
            cg.print(std::cerr);
        }
        collectContigInfos();
        collectMateLinks();
        buildReducedLinks();
        augmentGraph();
    }

private:

    void collectContigInfos();
    void collectMateLinks();
    void buildReducedLinks();
    void removeContainmentLinks();
    void augmentGraph();

    // Input.
    ContigGraph & cg;
    MateInfos const & mateInfos;
    // Contig sequences to use for containment overlap removal.
    seqan::StringSet<seqan::Dna5String> const * contigSeqs;
    bool logging;

    // Links and contig infos collected from input.
    std::vector<scaffolder::MateLink> mateLinks;
    std::vector<scaffolder::ContigEdgeLabel> contigInfos;

    // Reduced links computed from mateLinks using link bundling and transitive reduction.
    std::vector<scaffolder::MateLink> reducedLinks;
};

void PairedEndLinkAdder::collectContigInfos()
{
    for (auto node : cg.node)
        contigInfos.push_back(scaffolder::ContigEdgeLabel(cg.contig[node].length));
}

void PairedEndLinkAdder::collectMateLinks()
{
    for (auto const & info : mateInfos.records)
    {
        // std::cerr << "mateInfos.libraries[" << info.libraryID << "].mean == "
        //           << mateInfos.libraries[info.libraryID].mean
        //           << "\tcontigInfos[" << info.leftID << "].length == "
        //           << contigInfos[info.leftID].length
        //           << "\tinfo.beginPos == " << info.beginPos
        //           << "\tinfo.endPos == " << info.endPos << "\n";
        double mean = mateInfos.libraries[info.libraryID].mean;
        mean -= contigInfos[info.leftID].length - info.beginPos;
        mean -= info.endPos;
        double sd = mateInfos.libraries[info.libraryID].sd;
        mateLinks.push_back(scaffolder::MateLink(
                info.leftID, info.rightID, scaffolder::MateEdgeLabel(mean, sd, 1, 1.0 / info.numCopies)));
        if (logging)
            std::cerr << "COLLECTING\t" << mateLinks.back() << "\n";
    }
}

void PairedEndLinkAdder::buildReducedLinks()
{
    using scaffolder::MateLink;

    // Bundle links and perform transitive reduction.
    scaffolder::PathMergingOptions options;
    std::vector<MateLink> bundledLinks = bundleLinks(mateLinks, options);
    if (logging)
    {
        std::cerr << "BUNDLED LINKS\n";
        std::copy(bundledLinks.begin(), bundledLinks.end(),
                  std::ostream_iterator<MateLink>(std::cerr, "\n"));
    }
    transitiveReduction(reducedLinks, bundledLinks, contigInfos, options);

    // Remove links with too low weight.
    eraseIf(reducedLinks, [&options](MateLink const & mateLink) {
            return (mateLink.label.weight < options.minWeight);
        });

    // Pick label with highest count of all pairs.
    auto lt = [](MateLink const & lhs, MateLink const & rhs) {
        return (std::make_tuple(lhs.source, lhs.target, -(int)lhs.label.weight) <
                std::make_tuple(rhs.source, rhs.target, -(int)rhs.label.weight));
    };
    std::sort(reducedLinks.begin(), reducedLinks.end(), lt);

    auto eq = [](MateLink const & lhs, MateLink const & rhs) {
        return (std::make_pair(lhs.source, lhs.target) == std::make_pair(rhs.source, rhs.target));
    };
    reducedLinks.resize(std::distance(
            reducedLinks.begin(),
            std::unique(reducedLinks.begin(), reducedLinks.end(), eq)));

    if (logging)
    {
        std::cerr << "REDUCED LINKS\n";
        std::copy(reducedLinks.begin(), reducedLinks.end(),
                  std::ostream_iterator<MateLink>(std::cerr, "\n"));
    }

    if (contigSeqs)
    {
        removeContainmentLinks();  // also removes links not verified by overlap!
        if (logging)
        {
            std::cerr << "REDUCED LINKS WITHOUT CONTAINMENTS\n";
            std::copy(reducedLinks.begin(), reducedLinks.end(),
                      std::ostream_iterator<MateLink>(std::cerr, "\n"));
        }
    }
}

void PairedEndLinkAdder::removeContainmentLinks()
{
    SEQAN_CHECK(contigSeqs, "Must not be null.");
    auto const & contigs = *contigSeqs;
    assembler::Overlapper overlapper;  // default options are OK for now

    eraseIf(reducedLinks, [&](scaffolder::MateLink const & mateLink) {
            int const MIN_BAND = 20;
            int const k = 3;
            int diagonal = (int)length(contigs[mateLink.source]) - mateLink.label.lengthMean;
            int band = std::max((int)(k * mateLink.label.lengthStdDev), MIN_BAND);
            auto ovl = overlapper.computeOverlap(contigs[mateLink.source], contigs[mateLink.target],
                                                 mateLink.source, mateLink.target,
                                                 diagonal, band);
            if (ovl.errors == assembler::Overlap::INVALID)
                return true;  // could not verify overlap, remove
            ovl = ovl.normalize();  // make sure it is normalized
            bool result = (ovl.begin1 + ovl.len1 <= ovl.len0);
            // std::cerr << "CONSIDERING CONTAINMENT OVERLAP\t" << ovl << " => " << result << "\n";
            return result;
        });
}

void PairedEndLinkAdder::augmentGraph()
{
    // At this point, reducedLinks contains non-redundant links.  We can simply add edge to the graph or add a label to
    // existing ones.

    // TODO(holtgrew): Differentiate unique and redundant links!

    // Build map of (left, right) => edge.
    std::map<std::pair<unsigned, unsigned>, lemon::SmartGraph::Edge> edges;
    for (lemon::SmartGraph::EdgeIt edge(cg.graph); edge != lemon::INVALID; ++edge)
    {
        auto const & label = cg.link[edge];
        std::pair<unsigned, unsigned> key(label.leftID, label.rightID);
        if (key.first > key.second)
            std::swap(key.first, key.second);
        edges[key] = edge;
    }

    // Now, add the links indicated by reducedLinks to cg.graph with a new label or update the given edge label.
    for (auto const & mateLink : reducedLinks)
    {
        unsigned leftID = mateLink.source;
        unsigned rightID = mateLink.target;
        std::pair<unsigned, unsigned> key(leftID, rightID);
        if (key.first > key.second)
            std::swap(key.first, key.second);
        auto it = edges.find(key);
        if (key.first != leftID)
            continue;  // ignore, link for reverse order already in graph
        lemon::SmartGraph::Edge e;
        if (it != edges.end())
        {
            e = it->second;
            cg.link[it->second].uniqueLinks = mateLink.label.weight;
            cg.link[it->second].duplicateLinks = mateLink.label.count;
            cg.link[it->second].uniqueLinkDistance.mean = mateLink.label.lengthMean;
            cg.link[it->second].uniqueLinkDistance.sd = mateLink.label.lengthStdDev;
        }
        else
        {
            EdgeLabel label(leftID, rightID, mateLink.label.weight, mateLink.label.count);
            label.uniqueLinkDistance.mean = mateLink.label.lengthMean;
            label.uniqueLinkDistance.sd = mateLink.label.lengthStdDev;
            cg.addEdge(cg.node[leftID], cg.node[rightID], label);
        }
    }
}

}  // anonymous namespace

void addPairedEndLinks(ContigGraph & cg,
                       MateInfos const & mateInfos,
                       seqan::StringSet<seqan::Dna5String> const * contigSeqs)
{
    PairedEndLinkAdder helper(cg, mateInfos, contigSeqs, /*logging=*/false);
    helper.run();
}

}  // namespace rep_solv
