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

#include "dead_branch.h"

#include <seqan/basic.h>

#include "rep_solv/contig_graph.h"
#include "rep_solv/options.h"

namespace rep_solv {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Class DeadBranchRemover
// ----------------------------------------------------------------------------

// Helper class for removeDeadBranches().

class DeadBranchRemover
{
public:
    DeadBranchRemover(lemon::SmartGraph::NodeMap<bool> & doRemove,
                      ContigGraph const & graph,
                      unsigned stepNo,
                      std::vector<int> const & readBirthStep,
                      Options const & options) :
            doRemove(doRemove), graph(graph), stepNo(stepNo), readBirthStep(readBirthStep), options(options)
    {
        contigBirthStep = buildContigBirthStep();
    }

    unsigned run();

private:

    // Start DFS from graph.s to the right.
    void dfsRight(lemon::SmartGraph::NodeMap<int> & maxStep);
    // Recursion helper for DFS from graph.s to the right.
    int dfsRightRec(lemon::SmartGraph::NodeMap<bool> & reached,
                    lemon::SmartGraph::NodeMap<int> & maxStep,
                    lemon::SmartGraph::Node u);
    // Start DFS from graph.t to the left.
    void dfsLeft(lemon::SmartGraph::NodeMap<int> & maxStep);
    // Recursion helper for DFS from graph.t to the left.
    int dfsLeftRec(lemon::SmartGraph::NodeMap<bool> & reached,
                   lemon::SmartGraph::NodeMap<int> & maxStep,
                   lemon::SmartGraph::Node u);

    // Build value for contigBirthStep member.
    std::vector<int> buildContigBirthStep() const;

    // Reads for each contig.
    std::vector<std::vector<unsigned>> contigReads;
    // Maximal "birth step" for each contig.
    std::vector<int> contigBirthStep;

    // Input / output.
    lemon::SmartGraph::NodeMap<bool> & doRemove;
    ContigGraph const & graph;
    unsigned stepNo;
    std::vector<int> const & readBirthStep;
    Options const & options;
};

std::vector<int> DeadBranchRemover::buildContigBirthStep() const
{
    std::vector<int> result(graph.node.size(), 0);

    for (unsigned readID = 0; readID < readBirthStep.size(); ++readID)
        for (unsigned contigID : graph.readToContigs.at(readID))
            result[contigID] = std::max(result[contigID], readBirthStep[readID]);

    for (unsigned contigID = 0; contigID < result.size(); ++contigID)
        if (options.verbosity >= 2)
            std::cerr << "contigBirthStep[" << contigID << "] == " << result[contigID] << "\n";

    return result;
}

unsigned DeadBranchRemover::run()
{
    // Compute maximal birth step for each contig.
    lemon::SmartGraph::NodeMap<int> maxStepR(graph.graph, 0), maxStepL(graph.graph, 0);
    dfsRight(maxStepR);
    dfsLeft(maxStepL);

    // Combine information from left and right search.
    lemon::SmartGraph::NodeMap<int> latestUpdate(graph.graph, 0);
    for (lemon::SmartGraph::NodeIt u(graph.graph); u != lemon::INVALID; ++u)
    {
        latestUpdate[u] = std::max(maxStepL[u], maxStepR[u]);
        if (options.verbosity >= 2)
            std::cerr << "Latest update " << graph.contig[u].id
                      << " is " << latestUpdate[u] << "\n";
    }

    unsigned numRemoved = 0;

    // Update output bool map.
    for (lemon::SmartGraph::NodeIt u(graph.graph); u != lemon::INVALID; ++u)
        if (u != graph.s && u != graph.t)
        {
            doRemove[u] = (stepNo - latestUpdate[u] >= options.maxAge);
            if (options.verbosity >= 2 && doRemove[u])
                std::cerr << "Removing " << graph.contig[u].id << " (age == " << (stepNo - latestUpdate[u]) << ")\n";
            numRemoved += !!doRemove[u];
        }

    return numRemoved;
}

void DeadBranchRemover::dfsRight(lemon::SmartGraph::NodeMap<int> & maxStep)
{
    // For all out arcs of s.
    lemon::SmartGraph::NodeMap<bool> reached(graph.graph, false);
    reached[graph.s] = true;
    for (lemon::SmartGraph::OutArcIt arc(graph.graph, graph.s); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(graph.s));
        dfsRightRec(reached, maxStep, graph.graph.target(arc));
    }
}

int DeadBranchRemover::dfsRightRec(lemon::SmartGraph::NodeMap<bool> & reached,
                                   lemon::SmartGraph::NodeMap<int> & maxStep,
                                   lemon::SmartGraph::Node u)
{
    if (options.verbosity >= 2)
        std::cerr << "dfs right rec reached " << graph.graph.id(u) << "\n";
    if (u == graph.t)  // short-circuit when reaching t
        return 0;
    if (reached[u])
        return maxStep[u];
    reached[u] = true;

    // Start with contig birth step.
    maxStep[u] = contigBirthStep.at(graph.contig[u].id);

    // Get maximal birth step of child.
    for (lemon::SmartGraph::OutArcIt arc(graph.graph, u); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.source(arc)), graph.graph.id(u));
        if (graph.link[arc].leftID != graph.contig[u].id)
            continue;  // skip, link has wrong direction
        maxStep[u] = std::max(maxStep[u], dfsRightRec(reached, maxStep, graph.graph.target(arc)));
    }

    if (options.verbosity >= 2)
        std::cerr << "DFS RIGHT maxStep[" << graph.contig[u].id << "] == " << maxStep[u] << "\n";

    return maxStep[u];
}

void DeadBranchRemover::dfsLeft(lemon::SmartGraph::NodeMap<int> & maxStep)
{
    // For all out arcs of t.
    lemon::SmartGraph::NodeMap<bool> reached(graph.graph, false);
    reached[graph.t] = true;
    for (lemon::SmartGraph::InArcIt arc(graph.graph, graph.t); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.target(arc)), graph.graph.id(graph.t));
        dfsLeftRec(reached, maxStep, graph.graph.source(arc));
    }
}

int DeadBranchRemover::dfsLeftRec(lemon::SmartGraph::NodeMap<bool> & reached,
                                  lemon::SmartGraph::NodeMap<int> & maxStep,
                                  lemon::SmartGraph::Node u)
{
    if (options.verbosity >= 2)
        std::cerr << "dfs left rec reached " << graph.graph.id(u) << "\n";
    if (u == graph.s)  // short-circuit when reaching s
        return 0;
    if (reached[u])
        return maxStep[u];
    reached[u] = true;

    // Start with contig birth step.
    maxStep[u] = contigBirthStep.at(graph.contig[u].id);

    // Get maximal birth step of child.
    for (lemon::SmartGraph::InArcIt arc(graph.graph, u); arc != lemon::INVALID; ++arc)
    {
        SEQAN_ASSERT_EQ(graph.graph.id(graph.graph.target(arc)), graph.graph.id(u));
        if (graph.link[arc].rightID != graph.contig[u].id)
            continue;  // skip, link has wrong direction
        maxStep[u] = std::max(maxStep[u], dfsLeftRec(reached, maxStep, graph.graph.source(arc)));
    }

    if (options.verbosity >= 2)
        std::cerr << "DFS LEFT maxStep[" << graph.contig[u].id << "] == " << maxStep[u] << "\n";

    return maxStep[u];
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function removeDeadBranches()
// ----------------------------------------------------------------------------

// Perform dead branch removal.

unsigned removeDeadBranches(lemon::SmartGraph::NodeMap<bool> & doRemove,
                            ContigGraph const & graph,
                            unsigned stepNo,
                            std::vector<int> const & readBirthStep,
                            Options const & options)
{
    DeadBranchRemover helper(doRemove, graph, stepNo, readBirthStep, options);
    return helper.run();
}

}  // namespace rep_solv
