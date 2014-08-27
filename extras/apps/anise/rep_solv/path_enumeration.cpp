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

// TODO(holtgrew): A number of the mappings created below are superflous and should be removed.
// TODO(holtgrew): We can represent sequences by nodes again instead of edges!

#include "path_enumeration.h"

#include <iostream>
#include <memory>

#include <lemon/connectivity.h>
#include <lemon/bfs.h>
#include <lemon/matching.h>
#include <lemon/path.h>
#include <lemon/smart_graph.h>

#include <seqan/basic.h>

#include "rep_sep/read_set.h"

#include "rep_solv/contig_graph.h"
#include "rep_solv/options.h"

#include "scaffolder/overlap_resolution.h"
#include "scaffolder/scaffolding_result.h"

namespace rep_solv {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Function appendPath
// ----------------------------------------------------------------------------

template <typename TPath>
void appendPath(TPath & out, TPath const & in)
{
    for (auto arc = in.nthIt(0); arc != lemon::INVALID; ++arc)
        out.addBack(arc);
}

// ----------------------------------------------------------------------------
// Function dagify()
// ----------------------------------------------------------------------------

class DagifyHelper
{
    enum Flag { UNMARKED, TEMPORARY, PERMANENT };

public:
    DagifyHelper(lemon::ListDigraph & g) : g(g), flags(g, UNMARKED) {}

    void run();

    // Allow access to erased arcs.
    std::vector<lemon::ListDigraph::Arc> erasedArcs() const
    { return eraseArcs; }

private:
    void visit(lemon::ListDigraph::Node u, lemon::ListDigraph::Arc arc);

    // The graph to DAGify.
    lemon::ListDigraph & g;
    // Flag for the graph's nodes.
    lemon::ListDigraph::NodeMap<Flag> flags;
    // Collect arcs to erase.
    std::vector<lemon::ListDigraph::Arc> eraseArcs;
};

void DagifyHelper::run()
{
    for (lemon::ListDigraph::NodeIt u(g); u != lemon::INVALID; ++u)
        if (flags[u] == UNMARKED)
            visit(u, lemon::INVALID);
    for (auto arc : eraseArcs)
        g.erase(arc);
}

void DagifyHelper::visit(lemon::ListDigraph::Node u, lemon::ListDigraph::Arc arc)
{
    if (flags[u] == TEMPORARY)
    {
        // arc violates DAG structure
        SEQAN_ASSERT(arc != lemon::INVALID);
        SEQAN_ASSERT(std::find(eraseArcs.begin(), eraseArcs.end(), arc) == eraseArcs.end());
        eraseArcs.push_back(arc);
    }
    else if (flags[u] == UNMARKED)
    {
        flags[u] = TEMPORARY;
        for (lemon::ListDigraph::OutArcIt outArc(g, u); outArc != lemon::INVALID; ++outArc)
            visit(g.target(outArc), outArc);
        flags[u] = PERMANENT;
    }
}

void dagify(lemon::ListDigraph & g)
{
    DagifyHelper helper(g);
    helper.run();
}

// ----------------------------------------------------------------------------
// Class WeightedBipartiteGraph
// ----------------------------------------------------------------------------

// Bundles bipartite graph with edge weights;  used for path enumeration in PathEnumerator.

class WeightedBipartiteGraph
{
public:
    // Special values for contig IDs.
    static const unsigned INVALID = (unsigned)-1;

    struct EdgeLabel
    {
        double weight { 0.0 };
        unsigned source { INVALID };
        unsigned target { INVALID };
        lemon::ListDigraph::Arc arc { lemon::INVALID };  // mapping back to arc in WeightedDAG

        EdgeLabel() = default;
        EdgeLabel(double weight, unsigned source, unsigned target, lemon::ListDigraph::Arc arc) :
                weight(weight), source(source), target(target), arc(arc)
        {}
    };

    WeightedBipartiteGraph() : label(graph), weight(graph) {}

    void print(std::ostream & out) const;

    // The underlying graph structure.
    lemon::SmartGraph graph;
    // Labels for edges.
    lemon::SmartGraph::EdgeMap<EdgeLabel> label;
    // Weights for edges, required for matching algorithm.
    lemon::SmartGraph::EdgeMap<double> weight;
};

void WeightedBipartiteGraph::print(std::ostream & out) const
{
    out << "WeightedBipartiteGraph\n"
        << "EDGES\n";
    for (lemon::SmartGraph::EdgeIt e(graph); e != lemon::INVALID; ++e)
    {
        auto const & l = label[e];
        out << "\t" << graph.id(graph.u(e)) << " -- " << graph.id(graph.v(e))
            << "\t(contigs: " << l.source << " -> " << l.target << "\tweight=" << l.weight << ")\n";
    };
}

// ----------------------------------------------------------------------------
// Class WeightedDAG
// ----------------------------------------------------------------------------

// Bundles digraph with edge weights; used for path enumeration in PathEnumerator.

class WeightedDAG
{
public:
    // Special values for contig IDs.
    static const unsigned INVALID = (unsigned)-1;
    static const unsigned SOURCE = (unsigned)-2;
    static const unsigned TARGET = (unsigned)-3;

    struct NodeLabel
    {
        unsigned contigID { INVALID };
        double weight { 0.0 };
        // Corresponding nodes in bipartite graph.
        lemon::SmartGraph::Node bgNodeS { lemon::INVALID };
        lemon::SmartGraph::Node bgNodeT { lemon::INVALID };

        NodeLabel() = default;
        explicit NodeLabel(unsigned contigID) : contigID(contigID) {}
    };

    struct EdgeLabel
    {
        // Either edge represents link or contig (STRUCT = s/t link, TRANS = transitive closure, COMP = completing for
        // s-t paths.).
        enum Kind { INVALID, LINK, STRUCT, TRANS, COMP };
        // In case of contig, the contigID can be retrieved from labels of adjacent nodes.

        double weight { 0.0 };
        Kind kind { INVALID };
        lemon::SmartGraph::Edge bgEdge { lemon::INVALID };  // edge in bipartite graph
        lemon::SmartGraph::Edge cgEdge { lemon::INVALID };  // edge in input ContigGraph

        lemon::SimplePath<lemon::ListDigraph> containedArcs;  // contained arcs in case of kind == TRANS

        EdgeLabel() = default;
        EdgeLabel(double weight, Kind kind) : weight(weight), kind(kind) {}

        void clear()
        {
            weight = 0.0;
            kind = INVALID;
            bgEdge = lemon::INVALID;
            cgEdge = lemon::INVALID;
        }
    };

    WeightedDAG() : labelU(graph), labelE(graph)
    {
        init();
    }

    char const * kindStr(EdgeLabel::Kind kind) const
    {
        switch (kind)
        {
            case EdgeLabel::LINK:
                return "LINK";
            case EdgeLabel::TRANS:
                return "TRANS";
            case EdgeLabel::STRUCT:
                return "STRUCT";
            case EdgeLabel::COMP:
                return "COMP";
            default:
                return "INVALID";
        }
    }

    void print(std::ostream & out) const;

    // The nodes for s and t.
    lemon::ListDigraph::Node s, t;
    // The graph structure.
    lemon::ListDigraph graph;
    // Labels for the graph nodes.
    lemon::ListDigraph::NodeMap<NodeLabel> labelU;
    // Labels for the graph edges.
    lemon::ListDigraph::ArcMap<EdgeLabel> labelE;
    // Mapping from contigID to node in graph.
    std::vector<lemon::ListDigraph::Node> nodes;

private:

    void init()
    {
        s = graph.addNode();
        labelU[s].contigID = SOURCE;
        t = graph.addNode();
        labelU[t].contigID = TARGET;
    }
};

void WeightedDAG::print(std::ostream & out) const
{
    out << "WeightedDAG\n"
        << "\tS\t" << graph.id(s) << "\n"
        << "\tT\t" << graph.id(t) << "\n"
        << "NODES\n";
    for (lemon::ListDigraph::NodeIt u(graph); u != lemon::INVALID; ++u)
        out << "\t" << graph.id(u) << "\tcontigID=" << labelU[u].contigID << "\tweight=" << labelU[u].weight << "\n";
    out << "EDGES\n";
    for (lemon::ListDigraph::ArcIt a(graph); a != lemon::INVALID; ++a)
        out << "\t" << graph.id(graph.source(a)) << " -> " << graph.id(graph.target(a))
            << "\t(weight=" << labelE[a].weight << ", kind=" << kindStr(labelE[a].kind) << ")\n";
}

// ----------------------------------------------------------------------------
// Class PathEnumerator
// ----------------------------------------------------------------------------


// Perform path enumeration for the given ContigGraph.
//
// The aim is to build a list of s-t-paths that cover all contigs in the ContigGraph, ideally through highly weighted
// link/overlap edges.  We do this using Dilworth's/Hall's theorem computing bipartite matchings.  The approach follows
// the description in (Eriksson et al., 2008).
//
// The graph is converted to a DAG by a DFS search starting at s only using the undirected edges in left-to-right
// direction.  Each such traversed edge is converted into a directed edge and we convert our graph into an acyclic DAG.
// We then check whether there exists an s-t-path and bail out if there is no such path.  Each edge is weighted
// according to the paired/overlap links used.
//
// We can now compute a vertex-disjoint (maximum weight) path cover for this DAG G as follows.  We create a bipartite
// graph B where we add a vertex v1 to B for each vertex v from G that has non-zero out degree and a vertex v2 for each
// vertex v from G that has non-zero in degree.  We connect vertices u1 -> v2 in B if there is a path u ->..-> v in G.
//
// A (maximum weight) cover of B corresponds to a selection of paths in G.  We then complete each path starting in u and
// ending in v by an arbitrary s-u-path and an arbitrary v-t path.
//
// References:
//
//  * Eriksson, Nicholas, et al. "Viral population estimation using pyrosequencing." PLoS computational biology 4.5
//    (2008): e1000074.

class PathEnumerator
{
public:
    PathEnumerator(seqan::StringSet<TContigSeq> & outAssembly,
                   std::vector<AssemblyInfo> & outAssemblyInfos,
                   rep_sep::FeatureReadSet & outReadSet,
                   ContigGraph const & cg,
                   seqan::StringSet<TContigSeq> const & contigs,
                   rep_sep::FeatureReadSet const & coveringSet,
                   std::vector<int> const & readBirthStep,
                   std::vector<bool> const & overlapsWithFeature,
                   std::vector<AnchorType> const & anchors,
                   Options const & options) :
            outAssembly(outAssembly), outAssemblyInfos(outAssemblyInfos), outReadSet(outReadSet), cg(cg), contigs(contigs),
            coveringSet(coveringSet), readBirthStep(readBirthStep), overlapsWithFeature(overlapsWithFeature),
            anchors(anchors), options(options)
    {}

    void run();

private:

    // Build DAG to use for path enumeration from ContigGraph cg.
    std::unique_ptr<WeightedDAG> buildDAG() const;
    // Compute transitive closure of the graph.
    void transitiveClosure(WeightedDAG & dag) const;
    // Erase transitive edges again.
    void removeTransitiveEdges(WeightedDAG & dag) const;
    // Remove arcs that are not reachable from either s or t.
    void cleanseDAG(WeightedDAG & dag) const;
    // Add completing arcs.
    void addCompletingArcs(WeightedDAG & dag) const;
    // Build weighted bipartite from WeightedDAG, bgNode and bgEdge members of labels are updated in dag.
    std::unique_ptr<WeightedBipartiteGraph> buildBipartiteGraph(WeightedDAG & dag) const;
    // Use maximal weight bipartite matching algorithm to select edges from WeightedBipartiteGraph;
    std::vector<lemon::ListDigraph::Arc> matching(WeightedBipartiteGraph const & graph) const;
    // Compute heaviest s-t-path in dag.
    lemon::SimplePath<lemon::ListDigraph> heaviestPath(WeightedDAG const & dag) const;
    // Collect covering paths from matching.
    std::vector<lemon::SimplePath<lemon::ListDigraph>> coveringPaths(
            WeightedDAG const & dag,
            std::vector<lemon::ListDigraph::Arc> const & matching) const;
    // Complete the covering paths to s-t paths.
    std::vector<lemon::SimplePath<lemon::ListDigraph>> completePaths(
            WeightedDAG const & dag,
            std::vector<lemon::SimplePath<lemon::ListDigraph>> const & coveringPaths) const;
    // Build scaffolding result from paths.
    scaffolder::ScaffoldingResult buildScaffoldingResult(
            WeightedDAG const & dag,
            std::vector<lemon::SimplePath<lemon::ListDigraph>> const & coveringPaths) const;

    // Returns whether the path is a "clean" s-t-path without any completing links and thus the left and right assembly
    // "met".
    bool isMeetingPath(WeightedDAG const & dag,
                       lemon::SimplePath<lemon::ListDigraph> const & path) const
    {
        typedef lemon::SimplePath<lemon::ListDigraph> TPath;
        for (TPath::ArcIt arc(path); arc != lemon::INVALID; ++arc)
            if (dag.labelE[arc].kind == WeightedDAG::EdgeLabel::COMP)
                return false;
        return true;
    }

    // Print paths.
    void printCoveringPaths(std::ostream & out,
                            WeightedDAG const & dag,
                            std::vector<lemon::SimplePath<lemon::ListDigraph>> const & paths) const;

    template <typename TPath1, typename TPath2, typename TPath3>
    void concatenatePaths(lemon::SimplePath<lemon::ListDigraph> & out,
                          TPath1 & pS,
                          TPath2 & path,
                          TPath3 & pT,
                          WeightedDAG const & dag) const;

    // Input / Output.

    // Reference to the output assembly sequences.
    seqan::StringSet<TContigSeq> & outAssembly;
    // Output assembly infos.
    std::vector<AssemblyInfo> & outAssemblyInfos;
    // Read coverage information.
    rep_sep::FeatureReadSet & outReadSet;
    // Contig graph to enumerate the s-t paths for.
    ContigGraph const & cg;
    // Input contigs.
    seqan::StringSet<TContigSeq> const & contigs;
    // Covering read set selected earlier for contig expansion (non-conflicting subsets were merged before generating
    // the FeatureReadSet for expansion, though).  Can be used to mark nodes as being connected to "step 0" nodes.
    rep_sep::FeatureReadSet const & coveringSet;
    // Read birth step, used for selecting "step 0" nodes.
    std::vector<int> const & readBirthStep;
    // Whether or not a read overlaps with a feature.
    std::vector<bool> const & overlapsWithFeature;
    // Marker for whether a read is an anchor or not.
    std::vector<AnchorType> const & anchors;
    // Configuration of path enumerator.
    Options const & options;
};

void PathEnumerator::printCoveringPaths(
        std::ostream & out,
        WeightedDAG const & dag,
        std::vector<lemon::SimplePath<lemon::ListDigraph>> const & paths) const
{
    typedef lemon::SimplePath<lemon::ListDigraph> TPath;

    char const * label[5] = { "i", "l", "s", "t", "c" };

    out << "Covering Path List\n";
    unsigned pathID = 0;
    for (auto const & path : paths)
    {
        out << pathID++ << "\t"
            << dag.labelU[dag.graph.source(path.front())].contigID;
        for (TPath::ArcIt arc(path); arc != lemon::INVALID; ++arc)
            out << " -" << label[dag.labelE[arc].kind]
                << "-> " << dag.labelU[dag.graph.target(arc)].contigID;
        out << "\n";
    }
}

scaffolder::ScaffoldingResult PathEnumerator::buildScaffoldingResult(
        WeightedDAG const & dag,
        std::vector<lemon::SimplePath<lemon::ListDigraph>> const & coveringPaths) const
{
    typedef lemon::SimplePath<lemon::ListDigraph> const TPath;
    scaffolder::ScaffoldingResult result;

    for (auto const & path : coveringPaths)
    {
        scaffolder::ScaffoldingResult::TScaffold scaffold;
        int pos = 0, posSD = 0;

        SEQAN_ASSERT_GEQ(path.length(), 2);
        SEQAN_ASSERT(pathSource(dag.graph, path) == dag.s);
        SEQAN_ASSERT(pathTarget(dag.graph, path) == dag.t);

        for (auto arc = path.nthIt(1); arc != lemon::INVALID; ++arc)
        {
            unsigned contigID = dag.labelU[dag.graph.source(arc)].contigID;
            unsigned len = length(contigs[contigID]);
            scaffold.push_back(scaffolder::PositionedContig(pos, posSD, contigID, len));
            auto e = dag.labelE[arc].cgEdge;
            if (e == lemon::INVALID)
            {
                unsigned const MIN_DIST = 10;  // distance to add for complementing arcs
                pos += len + MIN_DIST;
            }
            else
            {
                auto label = cg.link[e];
                auto distance = label.uniqueLinks ? label.uniqueLinkDistance : label.overlapDistance;
                pos += len + distance.mean;
            }
        }

        result.scaffolds.push_back(scaffold);
    }

    return result;
}

void PathEnumerator::run()
{
    // Build weighted DAG from input graph and the corresponding bipartite graph.
    auto wDAG = buildDAG();
    if (options.verbosity >= 2)
    {
        std::cerr << "Weighted DAG\n";
        wDAG->print(std::cerr);
    }
    cleanseDAG(*wDAG);
    if (options.verbosity >= 2)
    {
        std::cerr << "Cleansed, weighted DAG\n";
        wDAG->print(std::cerr);
    }
    // Add arcs from all nodes with not outgoing arc to all nodes with no incoming arcs.
    addCompletingArcs(*wDAG);
    if (options.verbosity >= 2)
    {
        std::cerr << "DAG with more complete arcs\n";
        wDAG->print(std::cerr);
    }
    // Compute heaviest path.
    auto bestPath = heaviestPath(*wDAG);

    // Compute transitive closure.
    transitiveClosure(*wDAG);
    if (options.verbosity >= 2)
    {
        std::cerr << "Transitive closure.\n";
        wDAG->print(std::cerr);
    }
    auto wBG = buildBipartiteGraph(*wDAG);
    if (options.verbosity >= 2)
    {
        std::cerr << "Weighted Bipartite Graph\n";
        wBG->print(std::cerr);
    }

    // Compute matching an weighted bipartite graph and convert to covering paths.  Then, complete these paths.
    auto cPaths = coveringPaths(*wDAG, matching(*wBG));
    if (options.verbosity >= 2)
    {
        std::cerr << "Covering paths\n";
        printCoveringPaths(std::cerr, *wDAG, cPaths);
    }
    cPaths = completePaths(*wDAG, cPaths);
    if (options.verbosity >= 2)
    {
        std::cerr << "Completed covering paths\n";
        printCoveringPaths(std::cerr, *wDAG, cPaths);
    }

    // Prepend best path to covering paths.
    cPaths.insert(cPaths.begin(), bestPath);

    // Check whether paths are meeting.
    for (auto const & path : cPaths)
    {
        unsigned const MIN_ANCHORS = 3;
        unsigned anchoredLeft = 0;
        unsigned anchoredRight = 0;

        rep_sep::Read superRead;
        superRead.id = outReadSet.reads.size();
        for (lemon::SimplePath<lemon::ListDigraph>::ArcIt arc(path); arc != lemon::INVALID; ++arc)
            if (wDAG->graph.target(arc) != wDAG->t)
            {
                unsigned contigID = wDAG->labelU[wDAG->graph.target(arc)].contigID;
                for (unsigned readID : coveringSet.reads[contigID].subReads)
                    superRead.subReads.push_back(readID);
            }
        std::sort(superRead.subReads.begin(), superRead.subReads.end());

        // Count left anchors of first contig in path.
        {
            auto contigID = wDAG->labelU[wDAG->graph.target(path.front())].contigID;
            for (unsigned readID : coveringSet.reads[contigID].subReads)
                anchoredLeft += (anchors[readID] == AnchorType::LEFT);
        }
        // Count right anchors of first contig in path.
        {
            auto contigID = wDAG->labelU[wDAG->graph.source(path.back())].contigID;
            for (unsigned readID : coveringSet.reads[contigID].subReads)
                anchoredRight += (anchors[readID] == AnchorType::RIGHT);
        }

        bool spanning = (anchoredLeft > MIN_ANCHORS) && (anchoredRight > MIN_ANCHORS) && isMeetingPath(*wDAG, path);
        outAssemblyInfos.push_back(AssemblyInfo((anchoredLeft > MIN_ANCHORS),
                                                (anchoredRight > MIN_ANCHORS),
                                                spanning));
        outReadSet.reads.push_back(superRead);
    }

    // Build scaffolding result
    auto sr = buildScaffoldingResult(*wDAG, cPaths);
    if (options.verbosity >= 2)
    {
        std::cerr << "Scaffolding result\n";
        sr.print(std::cerr);
    }
    // Compute resulting scaffold sequences.
    resolveOverlaps(outAssembly, sr, contigs);
}

void PathEnumerator::removeTransitiveEdges(WeightedDAG & dag) const
{
    std::vector<lemon::ListDigraph::Arc> arcs;
    for (lemon::ListDigraph::ArcIt arc(dag.graph); arc != lemon::INVALID; ++arc)
        if (dag.labelE[arc].kind == WeightedDAG::EdgeLabel::TRANS)
            arcs.push_back(arc);
    for (auto arc : arcs)
        dag.graph.erase(arc);
}

void PathEnumerator::transitiveClosure(WeightedDAG & dag) const
{
    SEQAN_ASSERT_MSG(simpleGraph(dag.graph), "Must be simple!");

    lemon::ListDigraph::NodeMap<int> order(dag.graph);
    SEQAN_CHECK(checkedTopologicalSort(dag.graph, order), "Must be DAG");

    // Obtain list of nodes in reverse topological order.
    std::vector<lemon::ListDigraph::Node> nodes;
    for (lemon::ListDigraph::NodeIt u(dag.graph); u != lemon::INVALID; ++u)
        nodes.push_back(u);
    typedef lemon::ListDigraph::Node TNode;
    std::sort(nodes.begin(), nodes.end(), [&order](TNode u, TNode v) { return (order[u] > order[v]); });

    // Compute transitive closure in this order.
    std::set<unsigned> ids;  // reachable from current node
    for (auto node : nodes)
    {
        // We maintain a list of currently reachable nodes such that the graph remains simple.
        ids.clear();
        for (lemon::ListDigraph::OutArcIt arc(dag.graph, node); arc != lemon::INVALID; ++arc)
            ids.insert(dag.graph.id(dag.graph.target(arc)));

        // Extend transitive closure from node.
        for (lemon::ListDigraph::OutArcIt arc(dag.graph, node); arc != lemon::INVALID; ++arc)
        {
            auto u = dag.graph.target(arc);
            for (lemon::ListDigraph::OutArcIt arc2(dag.graph, u); arc2 != lemon::INVALID; ++arc2)
            {
                auto v = dag.graph.target(arc2);
                if (ids.count(dag.graph.id(v)))
                    continue;

                ids.insert(dag.graph.id(v));
                auto a = dag.graph.addArc(node, v);
                auto & label = dag.labelE[a];
                label.kind = WeightedDAG::EdgeLabel::TRANS;
                label.containedArcs.addBack(arc);

                auto const & label2 = dag.labelE[arc2];
                if (label2.kind == WeightedDAG::EdgeLabel::TRANS)
                    appendPath(label.containedArcs, label2.containedArcs);
                else
                    label.containedArcs.addBack(arc2);
            }
        }
    }

    SEQAN_ASSERT_MSG(simpleGraph(dag.graph), "Must still be simple!");
}

std::unique_ptr<WeightedDAG> PathEnumerator::buildDAG() const
{
    std::unique_ptr<WeightedDAG> result(new WeightedDAG);
    auto & out = *result;

    // Add vertices.
    out.nodes.resize(cg.node.size());
    for (lemon::SmartGraph::NodeIt u(cg.graph); u != lemon::INVALID; ++u)
    {
        if (u == cg.s || u == cg.t)
            continue;  // do not copy over s and t

        auto v = out.graph.addNode();
        out.nodes[cg.contig[u].id] = v;
        out.labelU[v].contigID = cg.contig[u].id;
    }

    // Weight nodes with the number of level 0 nodes that they have.
    for (auto const & pair : cg.readToContigs)
        if (readBirthStep.at(pair.first) == 0 && overlapsWithFeature.at(pair.first))
        {
            double weight = 1.0 / pair.second.size();
            for (unsigned contigID : pair.second)
                out.labelU[out.nodes.at(contigID)].weight += weight;
        }

    // TODO(holtgrew): Forward weights with linked contigs using coveringReadSet.

    // Add edges.
    for (lemon::SmartGraph::EdgeIt edge(cg.graph); edge != lemon::INVALID; ++edge)
    {
        auto link = cg.link[edge];
        if (link.leftID == EdgeLabel::SOURCE)
        {
            auto e = out.graph.addArc(out.s, out.nodes[link.rightID]);
            out.labelE[e].kind = WeightedDAG::EdgeLabel::STRUCT;
            out.labelE[e].cgEdge = edge;
        }
        else if (link.rightID == EdgeLabel::TARGET)
        {
            auto e = out.graph.addArc(out.nodes[link.leftID], out.t);
            out.labelE[e].kind = WeightedDAG::EdgeLabel::STRUCT;
            out.labelE[e].cgEdge = edge;
        }
        else
        {
            auto e = out.graph.addArc(out.nodes[link.leftID], out.nodes[link.rightID]);
            out.labelE[e].kind = WeightedDAG::EdgeLabel::LINK;
            out.labelE[e].weight = link.uniqueLinks;
            out.labelE[e].cgEdge = edge;
        }
    }

    // At this point, result->graph might still contain directed cycles. Now, DAG-ify it.
    dagify(result->graph);

    return result;
}

void PathEnumerator::addCompletingArcs(WeightedDAG & dag) const
{
    // This is called after directed tree growing, i.e. we do not have to check for reachability.

    // Compute topological sorting of the graph.
    lemon::ListDigraph::NodeMap<int> order(dag.graph);
    SEQAN_CHECK(checkedTopologicalSort(dag.graph, order), "Must be a DAG.");
    // Compute connected components (via undirected graph).
    lemon::SmartGraph tmpGraph;
    lemon::ListDigraph::NodeMap<lemon::SmartGraph::Node> tmpNode(dag.graph);
    graphCopy(undirector(dag.graph), tmpGraph).nodeRef(tmpNode).run();
    lemon::SmartGraph::NodeMap<int> component(tmpGraph);
    connectedComponents(tmpGraph, component);
    // We will only connect nodes in the component of s to the nodes in the component of t.  If the
    // two nodes are in the same component then we will only connect them if the left end is topological
    // smaller than the right end.
    int cS = component[tmpNode[dag.s]], cT = component[tmpNode[dag.t]];

    // Compute which node is a left or right end.
    lemon::ListDigraph::NodeMap<bool> leftEnd(dag.graph, true), rightEnd(dag.graph, true);
    for (lemon::ListDigraph::ArcIt arc(dag.graph); arc != lemon::INVALID; ++arc)
    {
        leftEnd[dag.graph.source(arc)] = false;
        rightEnd[dag.graph.target(arc)] = false;
    }

    // Collect nodes to compute completing paths for.
    std::vector<lemon::ListDigraph::Node> leftEnds, rightEnds;
    for (lemon::ListDigraph::NodeIt u(dag.graph); u != lemon::INVALID; ++u)
        if (leftEnd[u] != rightEnd[u])
        {
            if (leftEnd[u] && u != dag.t)
                leftEnds.push_back(u);
            if (rightEnd[u] && u != dag.s)
                rightEnds.push_back(u);
        }

    // Connect left and right ends.
    for (auto u : leftEnds)
        for (auto v : rightEnds)
        {
            if (component[tmpNode[u]] != cS || component[tmpNode[v]] != cT)
                continue;  // Skip, wrong components.
            if (component[tmpNode[u]] == component[tmpNode[v]] && order[u] >= order[v])
                continue;  // Skip, no loops and only correct order.

            auto arc = dag.graph.addArc(u, v);
            dag.labelE[arc] = WeightedDAG::EdgeLabel();
            dag.labelE[arc].kind = WeightedDAG::EdgeLabel::COMP;
        }
    // Connect left ends to t and right ends to s.
    for (auto u : leftEnds)
    {
        auto arc = dag.graph.addArc(u, dag.t);
        dag.labelE[arc] = WeightedDAG::EdgeLabel();
        dag.labelE[arc].kind = WeightedDAG::EdgeLabel::COMP;
    }
    for (auto u : rightEnds)
    {
        auto arc = dag.graph.addArc(dag.s, u);
        dag.labelE[arc] = WeightedDAG::EdgeLabel();
        dag.labelE[arc].kind = WeightedDAG::EdgeLabel::COMP;
    }

    // Check that there is an s-t-path in the DAG.
    lemon::Bfs<lemon::ListDigraph> checkBfs(dag.graph);
    SEQAN_CHECK(checkBfs.run(dag.s, dag.t), "There must be an s-t-path!");
}

void PathEnumerator::cleanseDAG(WeightedDAG & dag) const
{
    // Perform DFS from s.
    lemon::Dfs<lemon::ListDigraph>::ReachedMap reachedS(dag.graph, false);
    bfs(dag.graph).reachedMap(reachedS).run(dag.s);
    // Perform DFS from t on reverse graph.
    lemon::Dfs<lemon::ListDigraph>::ReachedMap reachedT(dag.graph, false);
    bfs(reverseDigraph(dag.graph)).reachedMap(reachedT).run(dag.t);

    // Remove unreachable edges.
    for (lemon::ListDigraph::ArcIt arc(dag.graph); arc != lemon::INVALID; ++arc)
        if (!reachedS[dag.graph.source(arc)] && !reachedT[dag.graph.source(arc)] &&
            !reachedS[dag.graph.target(arc)] && !reachedT[dag.graph.target(arc)])
        {
            dag.labelE[arc].clear();
            dag.graph.erase(arc);
        }
}

std::unique_ptr<WeightedBipartiteGraph> PathEnumerator::buildBipartiteGraph(WeightedDAG & dag) const
{
    std::unique_ptr<WeightedBipartiteGraph> result(new WeightedBipartiteGraph);
    auto & out = *result;

    // Build vertices.
    for (lemon::ListDigraph::NodeIt u(dag.graph); u != lemon::INVALID; ++u)
    {
        dag.labelU[u].bgNodeS = out.graph.addNode();
        dag.labelU[u].bgNodeT = out.graph.addNode();
    }

    // Add edges.
    for (lemon::ListDigraph::ArcIt arc(dag.graph); arc != lemon::INVALID; ++arc)
    {
        auto u = dag.labelU[dag.graph.source(arc)].bgNodeS;
        auto v = dag.labelU[dag.graph.target(arc)].bgNodeT;
        auto e = out.graph.addEdge(u, v);
        dag.labelE[arc].bgEdge = e;
        out.weight[e] = dag.labelE[arc].weight;
        out.label[e] = WeightedBipartiteGraph::EdgeLabel(dag.labelE[arc].weight,
                                                         dag.labelU[dag.graph.source(arc)].contigID,
                                                         dag.labelU[dag.graph.target(arc)].contigID,
                                                         arc);
    }

    return result;
}

std::vector<lemon::ListDigraph::Arc> PathEnumerator::matching(
        WeightedBipartiteGraph const & graph) const
{
    std::vector<lemon::ListDigraph::Arc> result;

    // Compute maximal matching.
    typedef lemon::MaxMatching<decltype(graph.graph)> TAlgo;
    TAlgo algo(graph.graph);
    algo.run();

    if (options.verbosity >= 2)
        std::cerr << "Matching Result\n";

    // Write out matched edges.
    for (lemon::SmartGraph::EdgeIt edge(graph.graph); edge != lemon::INVALID; ++edge)
        if (algo.matching(edge))
        {
            auto const & l = graph.label[edge];
            if (options.verbosity >= 2)
                std::cerr << "\t" << graph.graph.id(graph.graph.u(edge)) << " -- "
                          << graph.graph.id(graph.graph.v(edge))
                          << "\t(" << l.source << " -> " << l.target << "\tweight=" << l.weight << ")\n";
            result.push_back(l.arc);
        }

    return result;
}

lemon::SimplePath<lemon::ListDigraph> PathEnumerator::heaviestPath(WeightedDAG const & dag) const
{
    // Compute topological order.
    lemon::ListDigraph::NodeMap<int> order(dag.graph);
    SEQAN_CHECK(checkedTopologicalSort(dag.graph, order), "Must be a DAG.");
    std::vector<lemon::ListDigraph::Node> nodes;
    for (lemon::ListDigraph::NodeIt u(dag.graph); u != lemon::INVALID; ++u)
        nodes.push_back(u);
    std::sort(nodes.begin(), nodes.end(),
              [&order](lemon::ListDigraph::Node lhs, lemon::ListDigraph::Node rhs) {
                  return (order[lhs] < order[rhs]);
              });

    // Traverse DAG in topological order and compute weight of paths up to each node and predecessor.
    lemon::ListDigraph::NodeMap<lemon::ListDigraph::Arc> predArc(dag.graph, lemon::INVALID);
    lemon::ListDigraph::NodeMap<double> weight(dag.graph, 0);
    for (auto u : nodes)
    {
        double bestWeight = -1.0;
        lemon::ListDigraph::Arc bestPred = lemon::INVALID;
        for (lemon::ListDigraph::InArcIt arc(dag.graph, u); arc != lemon::INVALID; ++arc)
        {
            double mateCount = 0;
            if (u != dag.s && u != dag.t)
                mateCount = cg.contig[cg.node[dag.labelU[u].contigID]].mateCount;
            double w = weight[dag.graph.source(arc)] + dag.labelE[arc].weight +
                    dag.labelU[dag.graph.target(arc)].weight + mateCount;
            if (w > bestWeight)
            {
                bestWeight = w;
                bestPred = arc;
            }
        }
        weight[u] = std::max(0.0, bestWeight);
        predArc[u] = bestPred;
    }

    // Reconstruct best path up to t (must start at s).
    std::vector<lemon::ListDigraph::Arc> arcs;  // will collect in reverse order
    for (auto u = dag.t; u != lemon::INVALID; /*see below*/)
    {
        if (predArc[u] == lemon::INVALID)
            break;
        arcs.push_back(predArc[u]);
        u = dag.graph.source(predArc[u]);
    }
    SEQAN_ASSERT(dag.graph.source(arcs.back()) == dag.s);
    SEQAN_ASSERT(dag.graph.target(arcs.front()) == dag.t);

    lemon::SimplePath<lemon::ListDigraph> result;
    for (auto it = arcs.rbegin(); it != arcs.rend(); ++it)
        result.addBack(*it);
    return result;
}

std::vector<lemon::SimplePath<lemon::ListDigraph>> PathEnumerator::coveringPaths(
            WeightedDAG const & dag,
            std::vector<lemon::ListDigraph::Arc> const & matching) const
{
    std::vector<lemon::SimplePath<lemon::ListDigraph>> result;

    // Build mapping from contigID to selected out arc (if any).
    std::map<unsigned, lemon::ListDigraph::Arc> outArc;
    for (auto const & arc : matching)
        outArc[dag.labelU[dag.graph.source(arc)].contigID] = arc;

    // Compute topological sorting for dag.graph.
    lemon::ListDigraph::NodeMap<int> order(dag.graph);
    SEQAN_CHECK(checkedTopologicalSort(dag.graph, order), "Must be a DAG.");
    std::vector<lemon::ListDigraph::Node> nodes;
    for (lemon::ListDigraph::NodeIt u(dag.graph); u != lemon::INVALID; ++u)
        nodes.push_back(u);
    std::sort(nodes.begin(), nodes.end(), [&order](lemon::ListDigraph::Node lhs, lemon::ListDigraph::Node rhs) {
            return (order[lhs] < order[rhs]);
        });

    // Traverse dag.graph structure in topologically sorted order and start DFS in matching (using outArc mapping) and
    // collect paths into result.
    lemon::ListDigraph::NodeMap<bool> reached(dag.graph, false);
    for (auto u : nodes)
    {
        // Skip if reached already, mark as reached, and skip if there is no out arc for this contig.
        if (reached[u])
            continue;
        reached[u] = true;
        unsigned id = dag.labelU[u].contigID;
        if (!outArc.count(id))
            continue;

        // There is an arc from this (yet unreached) node that starts a matched (and to be covered) path in dag.graph.
        // Build a path for this.
        lemon::SimplePath<lemon::ListDigraph> path;
        while (outArc.count(id))
        {
            auto const arc = outArc.find(id)->second;
            path.addBack(arc);
            reached[dag.graph.target(arc)] = true;
            id = dag.labelU[dag.graph.target(arc)].contigID;
        }
        SEQAN_CHECK(checkPath(dag.graph, path), "Must be valid!");
        result.push_back(path);
    }

    return result;
}

// Helper function for path concatenation.
template <typename TPath1, typename TPath2, typename TPath3>
void PathEnumerator::concatenatePaths(lemon::SimplePath<lemon::ListDigraph> & out,
                      TPath1 & pS,
                      TPath2 & path,
                      TPath3 & pT,
                      WeightedDAG const & dag) const
{
    typedef lemon::ListPath<lemon::ListDigraph> TListPath;

    // Build result by concatenating paths.
    TListPath p;
    for (typename TPath3::RevArcIt arc(pT); arc != lemon::INVALID; ++arc)
    {
        if (options.verbosity >= 2)
            std::cerr << "ADDING TO FRONT " << dag.labelU[dag.graph.source(arc)].contigID << " -> "
                      << dag.labelU[dag.graph.target(arc)].contigID << "\n";
        p.addFront(arc);
    }
    SEQAN_CHECK(checkPath(dag.graph, p), "Must be valid!");
    for (int i = 0; i < path.length(); ++i)
    {
        auto arc = path.nth(path.length() - i - 1);
        if (options.verbosity >= 2)
            std::cerr << "ADDING TO FRONT " << dag.labelU[dag.graph.source(arc)].contigID << " -> "
                      << dag.labelU[dag.graph.target(arc)].contigID << "\n";
        p.addFront(arc);
    }
    SEQAN_CHECK(checkPath(dag.graph, p), "Must be valid!");
    for (typename TPath1::RevArcIt arc(pS); arc != lemon::INVALID; ++arc)
    {
        if (options.verbosity >= 2)
            std::cerr << "ADDING TO FRONT " << dag.labelU[dag.graph.source(arc)].contigID << " -> "
                      << dag.labelU[dag.graph.target(arc)].contigID << "\n";
        p.addFront(arc);
    }
    SEQAN_CHECK(checkPath(dag.graph, p), "Must be valid!");

    // Expand transitive edges.
    for (TListPath::ArcIt arc(p); arc != lemon::INVALID; ++arc)
        if (dag.labelE[arc].kind == WeightedDAG::EdgeLabel::TRANS)
            appendPath(out, dag.labelE[arc].containedArcs);
        else
            out.addBack(arc);
}

std::vector<lemon::SimplePath<lemon::ListDigraph>> PathEnumerator::completePaths(
        WeightedDAG const & dag,
        std::vector<lemon::SimplePath<lemon::ListDigraph>> const & coveringPaths) const
{
    typedef lemon::SimplePath<lemon::ListDigraph> TSimplePath;

    std::vector<lemon::SimplePath<lemon::ListDigraph>> result;

    // Augment result with covering paths.
    for (auto const & path : coveringPaths)
    {
        lemon::Bfs<lemon::ListDigraph> bfsS(dag.graph);
        bfsS.run(dag.s, pathSource(dag.graph, path));
        auto pS = bfsS.path(pathSource(dag.graph, path));
        lemon::Bfs<lemon::ListDigraph> bfsT(dag.graph);
        bfsT.run(pathTarget(dag.graph, path), dag.t);
        auto pT = bfsT.path(dag.t);

        TSimplePath out;
        concatenatePaths(out, pS, path, pT, dag);
        result.push_back(out);
    }

    // There might be nodes u in dag.graph that were not matched, we can pick arbitrary s-u-t paths for covering them.
    lemon::ListDigraph::NodeMap<bool> reached(dag.graph, false);
    for (auto & path : result)
        for (TSimplePath::ArcIt arc(path); arc != lemon::INVALID; ++arc)
        {
            reached[dag.graph.source(arc)] = true;
            reached[dag.graph.target(arc)] = true;
        }
    for (lemon::ListDigraph::NodeIt u(dag.graph); u != lemon::INVALID; ++u)
        if (!reached[u])
        {
            if (countOutArcs(dag.graph, u) == 0)  // not even s-/t-links, removed?
            {
                SEQAN_ASSERT_EQ(countInArcs(dag.graph, u), 0);
                continue;
            }
            lemon::Bfs<lemon::ListDigraph> bfsS(dag.graph);
            bfsS.run(dag.s, u);
            auto pS = bfsS.path(u);
            lemon::Bfs<lemon::ListDigraph> bfsT(dag.graph);
            bfsT.run(u, dag.t);
            auto pT = bfsT.path(dag.t);

            TSimplePath empty;
            TSimplePath out;
            concatenatePaths(out, pS, empty, pT, dag);
            result.push_back(out);
        }

    return result;
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Variable MASK_QUAL
// ----------------------------------------------------------------------------

int const MASK_QUAL = 1;

// ----------------------------------------------------------------------------
// Function enumeratePaths()
// ----------------------------------------------------------------------------

// Build scaffold of ContigGraph cleaned earlier.

void enumeratePaths(seqan::StringSet<TContigSeq> & outAssembly,
                    std::vector<AssemblyInfo> & outAssemblyInfos,
                    rep_sep::FeatureReadSet & outReadSet,
                    ContigGraph const & cg,
                    seqan::StringSet<TContigSeq> const & contigs,
                    rep_sep::FeatureReadSet const & coveringSet,
                    std::vector<int> const & readBirthStep,
                    std::vector<bool> const & overlapsWithFeature,
                    std::vector<AnchorType> const & anchors,
                    Options const & options)
{
    PathEnumerator helper(outAssembly, outAssemblyInfos, outReadSet, cg, contigs, coveringSet, readBirthStep,
                          overlapsWithFeature, anchors, options);
    helper.run();
}

}  // namespace rep_solv
