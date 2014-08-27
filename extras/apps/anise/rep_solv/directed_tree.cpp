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
                          ContigGraph const & graph) :
            doRemove(doRemove), graph(graph)
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
};

unsigned DirectedTreeHeuristic::run()
{
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
    // std::cerr << "reached right\t" << graph.contig[graph.s].id << "\n";
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
    if (doRemove[u])
        return;  // break recursion if to be removed, link does not exist.
    if (u == graph.t || reached[u])
        return;  // break recursion at s and reached nodes
    reached[u] = true;
    // std::cerr << "reached right\t" << graph.contig[u].id << "\n";

    // Iterate over left-to-right edges.
    for (lemon::SmartGraph::OutArcIt arc(graph.graph, u); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(u));
        if (graph.link[arc].leftID != graph.contig[u].id)
            continue;  // skip, link has wrong direction
        // std::cerr << "traversing " << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
        dfsRightRec(reached, graph.graph.target(arc));
    }
}

void DirectedTreeHeuristic::dfsLeft(lemon::SmartGraph::NodeMap<bool> & reached)
{
    reached[graph.t] = true;
    // std::cerr << "reached left\t" << graph.contig[graph.t].id << "\n";
    for (lemon::SmartGraph::OutArcIt arc(graph.graph, graph.t); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(graph.t));
        if (graph.link[arc].rightID != graph.contig[graph.t].id)
            continue;  // skip, link has wrong direction
        dfsLeftRec(reached, graph.graph.target(arc));
    }
}

void DirectedTreeHeuristic::dfsLeftRec(lemon::SmartGraph::NodeMap<bool> & reached, lemon::SmartGraph::Node u)
{
    if (doRemove[u])
        return;  // break recursion if to be removed, link does not exist.
    if (u == graph.t || reached[u])
        return;  // break recursion at t and reached nodes
    reached[u] = true;
    // std::cerr << "reached left\t" << graph.contig[u].id << "\n";

    // Iterate over right-to-left edges.
    for (lemon::SmartGraph::OutArcIt arc(graph.graph, u); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(u));
        if (graph.link[arc].rightID != graph.contig[u].id)
            continue;  // skip, link has wrong direction
        // std::cerr << "traversing " << graph.link[arc].leftID << " -- " << graph.link[arc].rightID << "\n";
        dfsLeftRec(reached, graph.graph.target(arc));
    }
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function directedTreeGrowing()
// ----------------------------------------------------------------------------

unsigned directedTreeGrowing(lemon::SmartGraph::NodeMap<bool> & doRemove,
                             ContigGraph const & graph)
{
    DirectedTreeHeuristic helper(doRemove, graph);
    return helper.run();
}

}  // namespace rep_solv
