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

#include "bog_unitigger.h"

#include <seqan/file.h>

#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>

namespace assembler {

// ----------------------------------------------------------------------------
// Class BogUnitigger
// ----------------------------------------------------------------------------

std::vector<unsigned> BogUnitigger::computeReadLengths(
        seqan::StringSet<seqan::Dna5String> const & seqs) const
{
    std::vector<unsigned> readLengths;
    for (unsigned i = 0; i < length(seqs); ++i)
        readLengths.push_back(length(seqs[i]));
    return readLengths;
}

std::unique_ptr<ContigGraph> BogUnitigger::augmentContigGraph(
        std::unique_ptr<ContigGraph> in,
        ContainmentFilter const & filter,
        seqan::StringSet<seqan::Dna5String> const & seqs) const
{
    std::unique_ptr<ContigGraph> result(in.release());
    ContigGraph & out = *result;

    // Build approximate read layouts for all contigs.
    std::vector<ApproximateReadLayout> layouts(out.node.size());

    // Fill layouts.
    for (unsigned contigID = 0; contigID < out.node.size(); ++contigID)
    {
        auto const & contig = out.contig[out.node[contigID]];
        auto & layout = layouts[contigID];

        for (auto posRead : contig.posReads)
            layout.add(ApproximateReadLayout::LayoutEntry(posRead.readID, posRead.beginPos, posRead.length));
        if (logging)
        {
            std::cerr << "INITIAL LAYOUT #" << contigID << "\n";
            layout.print(std::cerr);
        }
    }

    bool DEBUG_TOPOSORT = false;

    // Get contained ids and order them topologically so the anchoring/containing alignment has been inserted even for
    // stacks.
    lemon::SmartDigraph graph;
    lemon::SmartDigraph::NodeMap<unsigned> mapReadID(graph);
    std::map<unsigned, lemon::SmartDigraph::Node> nodes;
    if (DEBUG_TOPOSORT)
        std::cerr << "Starting toposort step\n";
    for (unsigned readID : filter.containedIDs())  // add vertices
    {
        auto u = graph.addNode();
        nodes[readID] = u;
        mapReadID[u] = readID;
        // std::cerr << "adding node for readID==" << readID << "\n";
    }
    for (unsigned readID : filter.containedIDs())  // add edges
    {
        auto ovl = filter.overlap(readID);  // longer -> shorter
        // std::cerr << "readID==" << readID << ", ovl==" << ovl << "\n";
        unsigned otherID = (ovl.seq0 == readID) ? ovl.seq1 : ovl.seq0;
        if (filter.isContained(otherID))
        {
            graph.addArc(nodes[ovl.seq0], nodes[ovl.seq1]);  // longer -> shorter for toposort
            if (DEBUG_TOPOSORT)
                std::cerr << "    arc " << readID << " -> " << otherID << "\n";
        }
    }
    lemon::SmartDigraph::NodeMap<int> order(graph);
    SEQAN_CHECK(checkedTopologicalSort(graph, order), "Must be a DAG!");
    std::vector<unsigned> sortedReadIDs(filter.containedIDs());  // copy
    std::sort(sortedReadIDs.begin(), sortedReadIDs.end(),
              [&](unsigned lhs, unsigned rhs) { return order[nodes[lhs]] < order[nodes[rhs]]; });
    if (DEBUG_TOPOSORT)
    {
        std::cerr << "DONE with toposort step\n";
        std::copy(sortedReadIDs.begin(), sortedReadIDs.end(), std::ostream_iterator<unsigned>(std::cerr, " "));
        std::cerr << "\n";
    }

    // Go over contained reads and put them into the approximate read layouts.
    for (unsigned readID : sortedReadIDs)
    {
        unsigned repr = filter.container(readID);  // representant!
        Overlap ovl = filter.overlap(readID);
        SEQAN_CHECK(ovl.seq0 == readID || ovl.seq1 == readID,
                    "ovl.seq0=%u, ovl.seq1=%u, readID=%u", ovl.seq0, ovl.seq1, readID);
        unsigned otherID = (ovl.seq0 == readID) ? ovl.seq1 : ovl.seq0;
        if (out.readToContig.count(readID))
            continue;  // Already placed, happens in case of contained overlap that is previously aligned.
        // SEQAN_CHECK(!out.readToContig.count(readID), "readID==%u", readID);
        SEQAN_CHECK(out.readToContig.count(repr), "representant=%u", repr);
        SEQAN_CHECK(out.readToContig.count(otherID), "otherID=%u", otherID);
        SEQAN_CHECK(out.readToContig[repr] == out.readToContig[otherID],
                    "repr=%u, otherID=%u", repr, otherID);

        auto & contig = out.contig[out.readToContig[repr]];
        auto & layout = layouts[contig.id];

        if (logging)
        {
            std::cerr << "Inserting " << readID << " into " << contig.id << " other read ID="
                      << otherID << "\n"
                      << "\tovl=" << ovl << "\n";
            layout.print(std::cerr);
        }

        int delta = (readID == ovl.seq0) ? ovl.begin1 : -(int)ovl.begin1;
        int otherPos = layout.readPos(otherID);
        layout.add(ApproximateReadLayout::LayoutEntry(
                readID, otherPos + delta, length(seqs[readID])));
        out.readToContig[readID] = out.readToContig[repr];
    }

    // Shift all layouts to 0.
    std::for_each(layouts.begin(), layouts.end(), [](ApproximateReadLayout & layout) { layout.shift(); });

    // Copy out layouts to contigs again.
    for (unsigned contigID = 0; contigID < out.node.size(); ++contigID)
    {
        out.contig[out.node[contigID]] = layouts[contigID].toContig(contigID);
        for (auto posContig : out.contig[out.node[contigID]].posReads)
        {
            auto readID = posContig.readID;
            out.readToContig[readID] = out.node[contigID];
            out.containedReadIDs.insert(readID);
        }
    }

    return result;
}

std::unique_ptr<ContigGraph> BogUnitigger::computeContigGraph(
        BogPathSet const & bogPathSet,
        ContainmentFilter const & filter,
        std::vector<Overlap> const & overlaps,
        seqan::StringSet<seqan::Dna5String> const & seqs) const
{
    std::unique_ptr<ContigGraph> result(new ContigGraph);
    ContigGraph & out = *result;
    (void)overlaps; // TODO(holtgrew): Remove!

    // Add nodes/contigs to ContigGraph.
    for (auto const & path : bogPathSet)
    {
// #if SEQAN_ENABLE_DEBUG  // NOT TRUE ANY MORE
//         for (unsigned readID : path.reads)
//             if (filter.isContained(readID))
//                 SEQAN_ASSERT_EQ_MSG(path.reads.size(), 1u, "readID=%u", readID);
// #endif  // #if SEQAN_ENABLE_DEBUG

        // Skip reads that are are the only one on the path since they are contained.  These will be added back later
        // with all other contained reads.
        if (path.reads.size() == 1u && filter.isContained(path.reads.front()))
            continue;
        
        Contig contig;
        contig.id = out.node.size();

        unsigned firstReadID = path.reads.front();
        contig.posReads.push_back(Contig::PositionedRead(firstReadID, 0, length(seqs[firstReadID])));
        int pos = 0;
        for (auto const & ovl : path.overlaps)
        {
            pos += ovl.begin1;
            contig.posReads.push_back(Contig::PositionedRead(ovl.id1, pos, length(seqs[ovl.id1])));
        }

        if (logging)
        {
            std::cerr << "CREATED CONTIG #" << contig.id << "\n";
            contig.print(std::cerr);
            std::cerr << "\n";
        }
        
        out.addContig(std::move(contig));
    }

    return result;
}

std::unique_ptr<ContigGraph> BogUnitigger::run(
        seqan::StringSet<seqan::Dna5String> const & seqs,
        std::vector<Overlap> const & overlaps) const
{
    // bool logging = this->logging || conflictStore != nullptr;
    // TODO(holtgrew): Use conflictStore for filtering out overlaps.

    // Filter out containment overlaps.
    if (logging)
        std::cerr << "Filtering out containment overlaps...\n";
    ContainmentFilter filter(logging);
    std::vector<Overlap> filteredOverlaps = filter.run(length(seqs), overlaps);
    if (logging)
    {
        std::cerr << "OVERLAPS without containment\n";
        for (auto const & ovl : filteredOverlaps)
            std::cerr << ovl << "\n";
    }
    if (logging)
    {
        std::cerr << "SEQS\n";
        for (unsigned i = 0; i < length(seqs); ++i)
            std::cerr << i << "\t" << seqs[i] << "\n";
        std::cerr << "CONTAINED READS\n";
        std::copy(filter.containedIDs().begin(), filter.containedIDs().end(),
                  std::ostream_iterator<unsigned>(std::cerr, "\n"));
        std::cerr << "/END OF CONTAINED READS\n";
    }
    
    // Create best overlap graph, remove cycles, and initial paths.
    if (logging)
        std::cerr << "Building Best Overlap Graph...\n";
    BestOverlapGraph bog = buildBestOverlapGraph(length(seqs), filteredOverlaps, logging);
    // if (bamRecords)
    //     updateBestOverlapGraphFromAlignments(bog, *bamRecords);
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
    }
    BogPathLengths pathLengths = computePathLengths(bog);
    auto readLengths = computeReadLengths(seqs);
    BogPathSet bps = buildInitialPaths(readLengths, bog, pathLengths, logging);
    // Perform path splitting, we can build our results from that.
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
    if (true || logging)
        std::cerr << "Performing Path Splitting...\n";
    if (true || logging)
    {
        std::cerr << "Paths before splitting\n";
        for (unsigned i = 0; i < bps.size(); ++i)
        {
            std::cerr << "path #" << i << "\n";
            std::cerr << bps.paths[i] << "\n";
        }
    }
    UnitiggerOptions options;  // TODO(holtgrew): Make configurable from the outside.
    performPathSplitting(bps, filter, options, true || logging);
    if (true || logging)
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
