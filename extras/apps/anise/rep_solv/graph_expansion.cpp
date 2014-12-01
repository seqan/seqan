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

#include "graph_expansion.h"

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include "rep_sep/read_set.h"
#include "rep_sep/feature_map.h"

#include "rep_solv/options.h"
#include "rep_solv/contig_graph.h"

namespace rep_solv {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Class ContigGraphExpander
// ----------------------------------------------------------------------------

class ContigGraphExpander
{
public:
    ContigGraphExpander(ContigGraph & cg,
                        std::vector<unsigned> & expandContigMap,
                        ContigGraph const & inCG,
                        seqan::StringSet<seqan::Dna5String> const & inContigs,
                        rep_sep::FeatureReadSet const & readSet,
                        Options const & options) :
            cg(cg), expandContigMap(expandContigMap), inCG(inCG), inContigs(inContigs), readSet(readSet),
            options(options)
    {}

    void run();

private:

    // Perform some initialization.
    void init();

    // Build old-to-new mapping.
    void buildOldToNew();
    // Create nodes in output.
    void createOutputNodes();
    // Create edges in output.
    void createOutputEdges();
    // Build readToContigs.
    void createReadToContigs();

    // Input / Output.

    ContigGraph & cg;
    std::vector<unsigned> & expandContigMap;
    ContigGraph const & inCG;
    seqan::StringSet<seqan::Dna5String> const & inContigs;
    rep_sep::FeatureReadSet const & readSet;
    Options const & options;

    // State.

    // Mapping from old contig ID to new contig ID and vice versa, built in init().
    std::vector<std::vector<unsigned>> oldToNew;
    std::vector<unsigned> newToOld;
};

void ContigGraphExpander::init()
{
    (void)options;  // ignore for now

    SEQAN_CHECK(std::is_sorted(readSet.reads.begin(), readSet.reads.end(), rep_sep::ltContigID),
                "Must be sorted by contigID.");
}

void ContigGraphExpander::run()
{
    buildOldToNew();
    createOutputNodes();
    createOutputEdges();
    createReadToContigs();

    swap(expandContigMap, newToOld);

    if (options.verbosity >= 2)
    {
        std::cerr << "NEW TO OLD MAP\n";
        for (unsigned i = 0; i < expandContigMap.size(); ++i)
            std::cerr << i << "\t->\t" << expandContigMap[i] << "\n";
    }
}

void ContigGraphExpander::buildOldToNew()
{
    oldToNew.resize(inCG.node.size());
    newToOld.resize(readSet.size(), (unsigned)-1);

    for (unsigned newID = 0; newID < readSet.reads.size(); ++newID)
    {
        newToOld.at(newID) = readSet.reads[newID].contigID;
        oldToNew.at(readSet.reads[newID].contigID).push_back(newID);
    }
}

void ContigGraphExpander::createOutputNodes()
{
    cg.node.resize(newToOld.size(), lemon::INVALID);
    for (unsigned newID = 0; newID < newToOld.size(); ++newID)
    {
        // Add node.
        auto u = cg.graph.addNode();
        cg.node[newID] = u;

        // Copy contig.
        cg.contig[u] = inCG.contig[inCG.node[newToOld[newID]]];
        cg.contig[u].id = newID;
        cg.contig[u].length = length(inContigs[newToOld.at(newID)]);
    }
}

void ContigGraphExpander::createOutputEdges()
{
    // Note that at this point, we can only have overlap edges since this function does not know about read
    // duplication.
    for (lemon::SmartGraph::EdgeIt edge(inCG.graph); edge != lemon::INVALID; ++edge)
    {
        auto label = inCG.link[edge];
        if (label.leftID == EdgeLabel::SOURCE)
        {
            for (auto right : oldToNew.at(label.rightID))
            {
                auto l = label;
                l.rightID = right;
                cg.addEdge(cg.s, cg.node[right], l);
            }
        }
        else if (label.rightID == EdgeLabel::TARGET)
        {
            for (auto left : oldToNew.at(label.leftID))
            {
                auto l = label;
                l.leftID = left;
                cg.addEdge(cg.node[left], cg.t, l);
            }
        }
        else
        {
            for (auto left : oldToNew.at(label.leftID))
                for (auto right : oldToNew.at(label.rightID))
                {
                    auto l = label;
                    l.leftID = left;
                    l.rightID = right;
                    cg.addEdge(cg.node[left], cg.node[right], l);
                }
        }
    }
}

void ContigGraphExpander::createReadToContigs()
{
    unsigned contigID = 0;
    for (auto const & read : readSet.reads)
    {
        for (auto readID : read.subReads)
            cg.readToContigs[readID].push_back(contigID);
        contigID++;
    }
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function expandContigGraph()
// ----------------------------------------------------------------------------

void expandContigGraph(ContigGraph & cg,
                       std::vector<unsigned> & expandContigMap,
                       ContigGraph const & inCG,
                       seqan::StringSet<seqan::Dna5String> const & inContigs,
                       rep_sep::FeatureReadSet const & readSet,
                       Options const & options)
{
    ContigGraphExpander helper(cg, expandContigMap, inCG, inContigs, readSet, options);
    helper.run();
}


// ----------------------------------------------------------------------------
// Function expandMateInfos()
// ----------------------------------------------------------------------------

void expandMateInfos(MateInfos & out,
                     MateInfos const & in,
                     rep_sep::FeatureReadSet const & readSet,  // must correspond to vertices
                     ContigGraph const & finalGraph,
                     Options const & options)
{
    out.libraries = in.libraries;

    // Expand mate infos.
    for (auto const & info : in.records)
    {
        if (!finalGraph.readToContigs.count(info.readL) ||
            finalGraph.readToContigs.find(info.readL)->second.empty() ||
            !finalGraph.readToContigs.count(info.readR) ||
            finalGraph.readToContigs.find(info.readR)->second.empty())
            continue;  // skip, at least one has no match any more
        auto const & leftIDs = finalGraph.readToContigs.find(info.readL)->second;
        auto const & rightIDs = finalGraph.readToContigs.find(info.readR)->second;
        unsigned numCopies = leftIDs.size() * rightIDs.size();
        for (auto leftID : leftIDs)
            for (auto rightID : rightIDs)
            {
                if (options.verbosity >= 3)
                    std::cerr << "conflicts\t" << leftID << "\t" << rightID << "\t"
                              << readSet.reads[leftID].conflicts(readSet.reads[rightID]) << "\n";
                if (options.maxMateConflicts > -1 &&
                    readSet.reads[leftID].conflicts(readSet.reads[rightID]) > options.maxMateConflicts)
                    continue;  // skip, too many conflicts
                auto outInfo = info;
                outInfo.leftID = leftID;
                outInfo.rightID = rightID;
                outInfo.numCopies = numCopies;
                out.insert(outInfo);
            }
    }

    // Update the numCopies members.
    std::stable_sort(out.records.begin(), out.records.end(), lt);
    auto it = out.records.begin(), itEnd = out.records.begin();
    while (itEnd != out.records.end())
    {
        // Skip over records.
        while (itEnd != out.records.end() && !lt(*it, *itEnd))
            ++itEnd;

        // Update count.
        for (auto it2 = it; it2 != itEnd; ++it2)
            it2->numCopies = (itEnd - it);

        it = itEnd;  // next batch
    }
}

}  // namespace rep_solv
