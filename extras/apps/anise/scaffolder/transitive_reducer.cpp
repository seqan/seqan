// ==========================================================================
//                                  ANISE
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

#include "gpm.h"

#include <iostream>
#include <vector>

#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>

#include "scaffolder/internal_shared.h"
#include "scaffolder/gpm_options.h"
#include "scaffolder/mate_link.h"
#include "scaffolder/utils.h"

namespace scaffolder {

namespace {  // anonymous namespace

// --------------------------------------------------------------------------
// Class ReductionGraph
// --------------------------------------------------------------------------

struct ReductionGraph
{
    // The graph with the underlying structure.
    lemon::SmartGraph graph;
    // The edge map for labeling the edges.
    lemon::SmartGraph::EdgeMap<ReductionGraphLabel> labels;
    // A mapping from idx to the given node.  Contig idx has two nodes (2 * idx) and (2 * idx + 1).
    std::vector<lemon::SmartGraph::Node> contigNodes;
    // Assigns nodes to the biconnected components and number of biconnected components.
    unsigned numBiComponents;
    lemon::SmartGraph::NodeMap<unsigned> biComponents;

    ReductionGraph() : labels(graph), numBiComponents(0), biComponents(graph)
    {}
};

// --------------------------------------------------------------------------
// Class PathEnumerationState
// --------------------------------------------------------------------------

template <typename TSubGraph>
struct PathEnumerationState
{
    enum  // Previous edge type.
    {
        MATE,
        CONTIG
    } edgeType;

    // The ReductionGraph.
    ReductionGraph const & graph;
    // The current path.
    lemon::SimplePath<TSubGraph> path;
    // Ids of the nodes on the current path.
    std::set<int> nodesOnPath;
    // These edges are removed in the following steps.

    PathEnumerationState(ReductionGraph const & graph, TSubGraph const & /*tag*/) :
            edgeType(CONTIG), graph(graph)
    {}

    void toggleEdgeType()
    {
        edgeType = (edgeType == CONTIG) ? MATE : CONTIG;
    }

    void printReachedPath(std::ostream & out, const char * label = nullptr) const
    {
        out << ",-- REACHED PATH (" << label << ")\n";
        for (int i = 0; i < path.length(); ++i)
        {
            out << "|  " << i << "\t::\t"
                << "" << graph.graph.id(graph.graph.source(path.nth(i))) << ""
                << " -> " << graph.graph.id(graph.graph.target(path.nth(i))) << "\t"
                << "IS CONTIG: " << graph.labels[path.nth(i)].isContig() << "\t"
                << "MEAN: " << graph.labels[path.nth(i)].lengthMean << "\t"
                << "DEV: " << graph.labels[path.nth(i)].lengthStdDev << "\t"
                << "IDX: " << graph.labels[path.nth(i)].idx << "\n";
        }
        out << "`---\n";
    }

    template <typename TArcIt>
    void pushEdge(TArcIt arc)
    {
        path.addBack(arc);
        toggleEdgeType();
        if (nodesOnPath.size() == 0u)
            nodesOnPath.insert(graph.graph.id(graph.graph.source(arc)));
        nodesOnPath.insert(graph.graph.id(graph.graph.target(arc)));
    }

    template <typename TArcIt>
    void popEdge(TArcIt arc)
    {
        nodesOnPath.erase(graph.graph.id(graph.graph.target(arc)));
        if (nodesOnPath.size() == 1u)
            nodesOnPath.erase(graph.graph.id(graph.graph.source(arc)));
        toggleEdgeType();
        path.eraseBack();
    }

    int firstNodeID() const
    {
        return graph.graph.id(graph.graph.source(path.front()));
    }

    int lastNodeID() const
    {
        return graph.graph.id(graph.graph.target(path.back()));
    }

    // Helper function for easier lookup in TransitiveReducer.
    std::pair<int, int> mapKey() const
    {
        std::pair<int, int> result(firstNodeID(), lastNodeID());
        if (result.first > result.second)
            std::swap(result.first, result.second);
        return result;
    }
};

template <typename TSubGraph>
PathEnumerationState<TSubGraph> pathEnumerationState(ReductionGraph const & graph, TSubGraph const & tag)
{
    return PathEnumerationState<TSubGraph>(graph, tag);
}

// --------------------------------------------------------------------------
// Class TransitiveReducer
// --------------------------------------------------------------------------

class TransitiveReducer
{
    // Input.
    //
    // Information about the links.
    std::vector<MateLink> const & links;
    // Information about the contig.
    std::vector<ContigEdgeLabel> const & contigInfos;
    // Configuration for path merging.
    PathMergingOptions const & options;

    // Result.
    //
    // The updated links.  We will first modify weights during updates and mark weight with 0 if removed.  In the end,
    // we will purge zero-weight links.
    std::vector<MateLink> result;

    typedef lemon::SmartGraph const TConstSmartGraph;

public:

    TransitiveReducer(std::vector<MateLink> const & links,
                      std::vector<ContigEdgeLabel> const & contigInfos,
                      PathMergingOptions const & options) :
            links(links), contigInfos(contigInfos), options(options)
    {}

    void run(std::vector<MateLink> & reduced)
    {
        // Start out with original links.
        result = links;

        // Build graph with from info in links and contigInfos.
        ReductionGraph graph;
        buildReductionGraph(graph);

        // Compute biconnected components.  We only have to enumerate paths in each component below.
        graph.numBiComponents = biEdgeConnectedComponents(graph.graph, graph.biComponents);
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            for (lemon::SmartGraph::NodeIt n(graph.graph); n != lemon::INVALID; ++n)
                std::cerr << graph.graph.id(n) << " -- " << graph.biComponents[n] << "\n";

        // The reduce edges are flagged with true.
        TConstSmartGraph::EdgeMap<bool> reducedEdges(graph.graph, false);

        // Enumerate paths in the biconnected components.  Paths are enumerated into incPaths, the to be deleted links
        // are collected in removeLinks.
        for (unsigned i = 0; i < graph.numBiComponents; ++i)
            enumeratePaths(reducedEdges, graph, i);

        // Purge links with empty weight and write out as result.
        auto itEnd = std::copy_if(result.begin(), result.end(), result.begin(),
                                  [](MateLink link) { return link.label.count > 0; });
        result.resize(itEnd - result.begin());
        swap(reduced, result);
    }

private:

    // Enumerate all paths in graph with the given component index.
    void enumeratePaths(TConstSmartGraph::EdgeMap<bool> & reducedEdges, ReductionGraph const & graph, unsigned idx)
    {
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "ENUMERATING PATHS FOR BI-COMPONENT " << idx << "\n"
                      << "========================================\n";

        lemon::SmartGraph::NodeMap<bool> enabled(graph.graph, false);
        for (typename lemon::SmartGraph::NodeIt n(graph.graph); n != lemon::INVALID; ++n)
            if (graph.biComponents[n] == idx)
                enabled[n] = true;
        lemon::FilterNodes<lemon::SmartGraph const> subgraph(graph.graph, enabled);

        auto state = pathEnumerationState(graph, subgraph);
        for (decltype(subgraph)::NodeIt n(subgraph); n != lemon::INVALID; ++n)
            enumeratePathRec(reducedEdges, state, subgraph, n, options.reductionPathLength);
    }

    // Enumerate paths of length up to num more vertices.
    template <typename TSubGraph, typename TNodeIt>
    void enumeratePathRec(TConstSmartGraph::EdgeMap<bool> & reducedEdges,
                          PathEnumerationState<TSubGraph> & state, TSubGraph const & subGraph, TNodeIt nIt,
                          unsigned num)
    {
        if (num == 0u)
            return;  // Add no more vertices to path.
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "ENUMERATE PATH REC FROM " << subGraph.id(nIt) << " NUM = " << num << "\n";
        for (typename TSubGraph::OutArcIt aIt(subGraph, nIt); aIt != lemon::INVALID; ++aIt)
        {
            // Skip ignored edges.
            if (reducedEdges[aIt])
                continue;  // Edge is already reduced.
            // Ensure the edge is not already on the path.
            if (state.nodesOnPath.count(subGraph.id(subGraph.target(aIt))))
                continue;  // Already on path.
            // Ensure we traverse the path alternatingly.
            bool isContigEdge = state.graph.labels[aIt].isContig();
            bool needContigEdge = (state.edgeType == PathEnumerationState<TSubGraph>::MATE);
            if (isContigEdge != needContigEdge)
                continue;  // No alternating traversal.

            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "TRAVERSING EDGE\t"
                          << state.graph.graph.id(state.graph.graph.source(aIt)) << "\t"
                          << state.graph.graph.id(state.graph.graph.target(aIt)) << "\n";

            // If we reach here then we select this edge.
            state.pushEdge(aIt);
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                state.printReachedPath(std::cerr, "after pushing");

            // Check if we should look at this edge.
            if (state.path.length() > 1 &&
                !state.graph.labels[state.path.front()].isContig() &&
                !state.graph.labels[state.path.back()].isContig())
                tryReduceEdge(reducedEdges, state, subGraph);

            // Recurse to find more edges.
            enumeratePathRec(reducedEdges, state, subGraph, subGraph.target(aIt), num - 1);

            // Restore state again.
            state.popEdge(aIt);
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                state.printReachedPath(std::cerr, "after popping");
        }
    }

    // Try to reduce edges that are subsumed by the path in the current state.
    template <typename TSubGraph>
    void tryReduceEdge(TConstSmartGraph::EdgeMap<bool> & reducedEdges,
                       PathEnumerationState<TSubGraph> & state,
                       TSubGraph const & subGraph)
    {
        // Collect edges connecting state.firstNodeID() and state.lastNodeID().
        std::vector<lemon::SmartGraph::Arc> edges;
        for (typename TSubGraph::OutArcIt arc(subGraph, subGraph.source(state.path.front())); arc != lemon::INVALID; ++arc)
        {
            if (reducedEdges[arc])
                continue;  // Edge is already reduced.
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "LOOKING AT\t" << state.graph.graph.id(subGraph.source(arc)) << " -> "
                          << state.graph.graph.id(subGraph.target(arc)) << "\n";
            SEQAN_ASSERT_EQ(state.graph.graph.id(subGraph.source(arc)), state.firstNodeID());
            if (!state.graph.labels[arc].isContig() && state.graph.graph.id(subGraph.target(arc)) == state.lastNodeID())
            {
                edges.push_back(arc);
                if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                    std::cerr << "    FOUND\t" << state.graph.graph.id(subGraph.source(arc)) << " -> "
                              << state.graph.graph.id(subGraph.target(arc)) << "\n";
            }
        }
        if (edges.empty())
            return;  // Nothing to reduce.
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "  ==> FOUND " << edges.size() << " CANDIDATES TO SUBSUME\n";

        // Compute length and deviation of the path.
        double pathLen = 0, pathLenSD = 0;
        for (int i = 0; i < state.path.length(); ++i)
        {
            pathLen += state.graph.labels[state.path.nth(i)].lengthMean;
            double x = state.graph.labels[state.path.nth(i)].lengthStdDev;
            pathLenSD += x * x;
        }
        pathLenSD = sqrt(pathLenSD);

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "PATH LENGTH:\t" << pathLen << "\n"
                      << "PATH SD:    \t" << pathLenSD << "\n";

        // Try to reduce each edge.
        for (auto e : edges)
        {
            double edgeLen = state.graph.labels[e].lengthMean;
            double edgeLenSD = state.graph.labels[e].lengthStdDev;
            int edgeCount = state.graph.labels[e].count;
            double edgeWeight = state.graph.labels[e].weight;

            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "EDGE LENGTH:\t" << edgeLen << "\n"
                          << "EDGE SD:    \t" << edgeLenSD << "\n"
                          << "EDGE COUNT: \t" << edgeCount << "\n"
                          << "EDGE WEIGHT:\t" << edgeWeight << "\n";

            if (abs(edgeLen - pathLen) < options.mult * std::max(pathLenSD, edgeLenSD))
            {
                if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                    std::cerr << "SUBSUMING EDGE\t" << state.graph.graph.id(state.graph.graph.source(e)) << " -> "
                              << state.graph.graph.id(state.graph.graph.target(e)) << "\n";
                // Adjust weights.
                unsigned y = state.graph.labels[e].idx;
                SEQAN_ASSERT_NEQ(y, seqan::maxValue<unsigned>());
                for (int i = 0; i < state.path.length(); i += 2 /* skip contig edges */)
                {
                    unsigned x = state.graph.labels[state.path.nth(i)].idx;
                    result[x].label.count += result[y].label.count;
                    result[x].label.weight += result[y].label.weight;
                }
                result[y].label.count = 0;
                result[y].label.weight = 0;
                // Mark as reduced.
                reducedEdges[e] = true;
            }
        }
    }

    void buildReductionGraph(ReductionGraph & graph) const
    {
        // Add vertices and edges for contigs.
        for (auto info : contigInfos)
        {
            lemon::SmartGraph::Node u = graph.graph.addNode(), v = graph.graph.addNode();
            graph.contigNodes.push_back(u);
            graph.contigNodes.push_back(v);
            lemon::SmartGraph::Edge e = graph.graph.addEdge(u, v);
            graph.labels[e] = ReductionGraphLabel(ReductionGraphLabel::CONTIG, info.length);
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "LABEL FOR EDGE\t"
                          << graph.graph.id(u) << " -> " << graph.graph.id(v)
                          << "\t" << graph.graph.id(e) << "\t IS \t" << graph.labels[e] << "\n";
        }

        // Add edges for links.
        for (auto indexedLink : enumerate(links))
        {
            auto const & link = indexedLink.second;
            unsigned idxU = 2 * link.source + 1;
            unsigned idxV = 2 * link.target;
            lemon::SmartGraph::Edge e = graph.graph.addEdge(graph.contigNodes[idxU], graph.contigNodes[idxV]);
            graph.labels[e] = ReductionGraphLabel(ReductionGraphLabel::MATE, link.label.lengthMean, link.label.lengthStdDev,
                                                  indexedLink.first, link.label.count, link.label.weight);
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "LABEL FOR EDGE\t" << idxU << " -> " << idxV << "\t"
                          << graph.graph.id(e) << "\t IS \t" << graph.labels[e] << "\n";
        }
    }
};

}  // anonymous namespace

// --------------------------------------------------------------------------
// Function transitiveReduction()
// --------------------------------------------------------------------------

void transitiveReduction(std::vector<MateLink> & reduced,
                         std::vector<MateLink> const & links,
                         std::vector<ContigEdgeLabel> const & contigInfos,
                         PathMergingOptions const & options)
{
    TransitiveReducer reducer(links, contigInfos, options);
    reducer.run(reduced);
}

} // namespace scaffolder
