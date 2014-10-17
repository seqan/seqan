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

#include "directed_tree.h"

#include <seqan/basic.h>

#include "rep_solv/contig_graph.h"
#include "rep_solv/options.h"

namespace rep_solv {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Class DirectedTreeHeuristic
// ----------------------------------------------------------------------------

// Helper class for directedTreeGrowing().

class DirectedTreeHeuristic
{
public:
    DirectedTreeHeuristic(lemon::SmartGraph::NodeMap<bool> & doRemove,
                          ContigGraph const & graph,
                          Options const & options) :
            doRemove(doRemove), graph(graph), options(options)
    {}

    unsigned run();

private:

    // Start DFS from graph.s to the right.
    void dfsRight(lemon::SmartGraph::NodeMap<bool> & reached);
    // Recursion helper for DFS from graph.s to the right.
    void dfsRightRec(lemon::SmartGraph::NodeMap<bool> & reached, lemon::SmartGraph::Node u);
    // Start DFS from graph.t to the left.
    void dfsLeft(lemon::SmartGraph::NodeMap<bool> & reached);
    // Recursion helper for DFS from graph.t to the left.
    void dfsLeftRec(lemon::SmartGraph::NodeMap<bool> & reached, lemon::SmartGraph::Node u);

    // Input / output.

    lemon::SmartGraph::NodeMap<bool> & doRemove;
    ContigGraph const & graph;
    Options const & options;
};

unsigned DirectedTreeHeuristic::run()
{
    if (options.verbosity >= 2)
    {
        std::cerr << "INPUT TO DIRECTED TREE HEURISTIC\n";
        graph.print(std::cerr);
    }

    lemon::SmartGraph::NodeMap<bool> reachedR(graph.graph, false), reachedL(graph.graph, false);
    dfsRight(reachedR);
    dfsLeft(reachedL);

    unsigned numUnreached = 0;

    // Build resulting doRemove map.
    for (lemon::SmartGraph::NodeIt u(graph.graph); u != lemon::INVALID; ++u)
    {
        doRemove[u] = (!reachedR[u] && !reachedL[u]);
        // if (doRemove[u])
        //     std::cerr << "uncovered by tree " << graph.contig[u].id << "\n";
        numUnreached += doRemove[u];
    }

    return numUnreached;
}

void DirectedTreeHeuristic::dfsRight(lemon::SmartGraph::NodeMap<bool> & reached)
{
    reached[graph.s] = true;
    if (options.verbosity >= 2)
        std::cerr << "dth reached right\t" << graph.contig[graph.s].id << "\n";
    for (lemon::SmartGraph::OutArcIt arc(graph.graph, graph.s); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(graph.s));
        if (graph.link[arc].leftID != graph.contig[graph.s].id)
            continue;  // skip, link has wrong direction
        dfsRightRec(reached, graph.graph.target(arc));
    }
}

void DirectedTreeHeuristic::dfsRightRec(lemon::SmartGraph::NodeMap<bool> & reached, lemon::SmartGraph::Node u)
{
    if (options.verbosity >= 2)
        std::cerr << "doRemove[u] == " << doRemove[u] << ", u == graph.t == " << (u == graph.t)
                  << ", reached[u] == " << reached[u] << "\n";  
    if (doRemove[u])
        return;  // break recursion if to be removed, link does not exist.
    if (u == graph.t || reached[u])
        return;  // break recursion at s and reached nodes
    reached[u] = true;
    if (options.verbosity >= 2)
        std::cerr << "dth reached right\t" << graph.contig[u].id << "\n";

    // Iterate over left-to-right edges.
    for (lemon::SmartGraph::OutArcIt arc(graph.graph, u); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(u));
        if (graph.link[arc].leftID != graph.contig[u].id)
            continue;  // skip, link has wrong direction
        if (options.verbosity >= 2)
            std::cerr << "dth right traversing " << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
        dfsRightRec(reached, graph.graph.target(arc));
    }
}

void DirectedTreeHeuristic::dfsLeft(lemon::SmartGraph::NodeMap<bool> & reached)
{
    reached[graph.t] = true;
    if (options.verbosity >= 2)
        std::cerr << "dth reached left\t" << graph.contig[graph.t].id << "\n";
    for (lemon::SmartGraph::InArcIt arc(graph.graph, graph.t); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.target(arc)), graph.graph.id(graph.t));
        if (graph.link[arc].rightID != graph.contig[graph.t].id)
            continue;  // skip, link has wrong direction
        dfsLeftRec(reached, graph.graph.source(arc));
    }
}

void DirectedTreeHeuristic::dfsLeftRec(lemon::SmartGraph::NodeMap<bool> & reached, lemon::SmartGraph::Node u)
{
    if (options.verbosity >= 2)
        std::cerr << "doRemove[u] == " << doRemove[u] << ", u == graph.t == " << (u == graph.t)
                  << ", reached[u] == " << reached[u] << "\n";  
    if (doRemove[u])
        return;  // break recursion if to be removed, link does not exist.
    if (u == graph.t || reached[u])
        return;  // break recursion at t and reached nodes
    reached[u] = true;
    if (options.verbosity >= 2)
      std::cerr << "dth reached left\t" << graph.contig[u].id << "\n";

    // Iterate over right-to-left edges.
    for (lemon::SmartGraph::InArcIt arc(graph.graph, u); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.target(arc)), graph.graph.id(u));
        if (graph.link[arc].rightID != graph.contig[u].id)
            continue;  // skip, link has wrong direction
        if (options.verbosity >= 2)
            std::cerr << "dth left traversing " << graph.link[arc].rightID << " -- " << graph.link[arc].leftID << "\n";
        dfsLeftRec(reached, graph.graph.source(arc));
    }
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function directedTreeGrowing()
// ----------------------------------------------------------------------------

unsigned directedTreeGrowing(lemon::SmartGraph::NodeMap<bool> & doRemove,
                             ContigGraph const & graph,
                             Options const & options)
{
    DirectedTreeHeuristic helper(doRemove, graph, options);
    return helper.run();
}

}  // namespace rep_solv
