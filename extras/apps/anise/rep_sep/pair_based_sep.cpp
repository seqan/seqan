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

#include "pair_based_sep.h"

#include <set>
#include <vector>
#include <lemon/smart_graph.h>

#include "asm/overlapper.h"

#include "rep_sep/read_set.h"
#include "rep_sep/feature_map.h"
#include "rep_sep/weighted_feature_vector.h"

#include "rep_solv/contig_graph.h"

#include "scaffolder/gpm.h"
#include "scaffolder/gpm_options.h"
#include "scaffolder/mate_link.h"

namespace rep_sep {

namespace {  // anonymous namespace

// --------------------------------------------------------------------------
// Class PairSeparatorGraph
// --------------------------------------------------------------------------

// Collect information required by the PairBasedSeparator.

struct PairSeparatorGraph
{
    struct ContigInfo
    {
        int contigID { 0 };
        int length   { 0 };
        int numReads { 0 };
    };

    struct LinkInfo
    {
        int lenMean { 0 };  // mean length
        int lenSD   { 0 };  // length standard deviation
        int weight  { 0 };  // number of pairs supporting the link
    };

    PairSeparatorGraph() : contig(graph), link(graph)
    {}

    lemon::SmartDigraph graph;
    std::vector<lemon::SmartDigraph::Node> nodes;
    lemon::SmartDigraph::NodeMap<ContigInfo> contig;
    lemon::SmartDigraph::ArcMap<LinkInfo> link;

    void print(std::ostream & out) const
    {
        out << "PAIR SEPARATOR GRAPH\n";
        out << "CONTIG INFOS\n";
        for (unsigned i = 0; i < nodes.size(); ++i)
            out << i << "\tid=" << graph.id(nodes[i]) << "\tlen=" << contig[nodes[i]].length << ", numReads="
                << contig[nodes[i]].numReads << "\n";
        out << "LINK INFOS\n";
        for (lemon::SmartDigraph::ArcIt arc(graph); arc != lemon::INVALID; ++arc)
            out << "id=" << graph.id(graph.source(arc)) << " -> " << graph.id(graph.target(arc))
                << "\tlenMen=" << link[arc].lenMean << ", lenSD=" << link[arc].lenSD
                << ", weight=" << link[arc].weight << "\n";
    }
};

// --------------------------------------------------------------------------
// Class PairBasedSeparator
// --------------------------------------------------------------------------

class PairBasedSeparator
{
public:
    // Returns one map for each new feature assigning where touched readIDs are assigned the feature value.
    typedef std::vector<std::map<unsigned, int>> TResult;

    PairBasedSeparator(TFragmentStore & store,
                       rep_solv::MateInfos const & mateInfos,
                       bool logging = false) :
            store(store), mateInfos(mateInfos), logging(logging)
    {}

    TResult run();

private:
    // Build PairSeparatorGraph from store.
    std::unique_ptr<PairSeparatorGraph> buildGraph(TFragmentStore const & store);
    // Compute begin/end position of alignment.
    template <typename TAlignedRead> std::pair<int, int> alignmentRange(TAlignedRead const & el) const;
    // Compute mate links.
    std::vector<scaffolder::MateLink> computeMateLinks(rep_solv::MateInfos const & mateInfos,
                                                       PairSeparatorGraph const & out);

    // // Search for overlaps indicated by link info that is not actually there.
    // void searchForConflictingOverlaps(PairSeparatorGraph const & psg,
    //                                   THandler handler) const;

    // Search for overlaps indicated by a branch.
    std::vector<std::map<unsigned, int>> searchForConflictingBranches(PairSeparatorGraph const & psg) const;

    // // Split because of overlap between the link between the given contigs and the distance mean.
    // void splitBecauseOfLink(PairSeparatorGraph const & psg, unsigned fromID, unsigned toID,
    //                         int lenMean, THandler handler) const;

    // Split because of a branch indicating an non-existing overlap.
    std::map<unsigned, int> splitBecauseOfBranch(unsigned fromID, unsigned toID0, unsigned toID1) const;

    // Mate indices, filled in run().
    seqan::String<unsigned> mateIdx;

    // Input.

    // Fragment store to perform the separation upon.
    TFragmentStore & store;
    // Paired-end library information.
    rep_solv::MateInfos const & mateInfos;
    // Flag to enable/disable logging.
    bool logging;
};

PairBasedSeparator::TResult PairBasedSeparator::run()
{
    calculateMateIndices(mateIdx, store);

    if (logging)
        std::cerr << "Running pair-based separator...\n";

    // Build the graph.
    auto psg = buildGraph(store);
    if (logging)
        psg->print(std::cerr);

    TResult result;

    // // Look for edges that indicate overlaps of contigs that do not exist.  Pairs that support this overlap can be
    // // clustered and separated against overlapping reads.
    // searchForConflictingOverlaps(*psg, handler);

    // Look for edges in Y-branches that do not exist.
    auto branchResult = searchForConflictingBranches(*psg);
    for (auto & res : branchResult)
        result.push_back(res);

    if (logging)
        std::cerr << "Done with pair-based separation.\n";

    return result;
}

std::unique_ptr<PairSeparatorGraph> PairBasedSeparator::buildGraph(TFragmentStore const & store)
{
    // Minimal support for each link and minimal coverage to require for contigs.
    unsigned const MIN_SUPPORT = 2;  // TODO(holtgrew): Get from configuration?

    std::unique_ptr<PairSeparatorGraph> result(new PairSeparatorGraph);
    PairSeparatorGraph & out = *result;

    // Add nodes for contigs and store their lengths.
    unsigned contigID = 0;
    for (auto const & contig : store.contigStore)
    {
        auto u = out.graph.addNode();
        out.nodes.push_back(u);
        out.contig[u].contigID = contigID++;
        out.contig[u].length = length(contig.seq);
    }
    // Count reads for each contig.
    for (auto const & el : store.alignedReadStore)
        ++out.contig[out.nodes[el.contigId]].numReads;

    // Compute links from reads aligning on the forward strand and bundle them.
    auto bundledLinks = bundleLinks(computeMateLinks(mateInfos, out),
                                    scaffolder::PathMergingOptions());
    // Remove the links with a too low support.
    auto itEnd = std::copy_if(bundledLinks.begin(), bundledLinks.end(), bundledLinks.begin(),
                              [MIN_SUPPORT](scaffolder::MateLink const & link) {
                                  return link.label.weight > MIN_SUPPORT;
                              });
    bundledLinks.resize(itEnd - bundledLinks.begin());

    // TODO(holtgrew): Discard links where the branching could be explained by overlap in reads:
    //                 Consider the graph 0->1 , 0->2.  We should not split anything in 0 if 1 and 2 overlap and thus
    //                 a consistent result could be achieved by good scaffold construction & merging.

    // TODO(holtgrew): Ignore links implying an overlap that is not there.

    // Construct links, if more than one link bundle exists for a connection then pick the one that has highest
    // support/weight.
    std::map<std::pair<unsigned, unsigned>, PairSeparatorGraph::LinkInfo> linkInfos;
    for (auto const & link : bundledLinks)
    {
        auto & linkInfo = linkInfos[std::make_pair(link.source, link.target)];
        if (link.label.weight > (unsigned)linkInfo.weight)
        {
            linkInfo.lenMean = link.label.lengthMean;
            linkInfo.lenSD = link.label.lengthStdDev;
            linkInfo.weight = link.label.weight;
        }
    }

    // Construct edges.
    for (auto const & pair : linkInfos)
    {
        auto arc = out.graph.addArc(out.nodes[pair.first.first], out.nodes[pair.first.second]);
        out.link[arc] = pair.second;
    }

    return result;
}

/*
void PairBasedSeparator::splitBecauseOfLink(PairSeparatorGraph const & psg, unsigned fromID, unsigned toID,
                                            int lenMean, THandler handler) const
{
    if (logging)
        std::cerr << "Splitting because of link " << fromID << " -> " << toID << "\n";

    auto const MULT = 3;
    auto const OVERLAP = 0.9;  // required overlap with linking reads

    TSeparatedReads clusters(2);  // build two clusters, supporting and non-supporting link

    unsigned support = 0;
    std::pair<int, int> rangeFrom, rangeTo;

    // Collect read alignments that support the link (fromID -> toID) with given distance.
    //
    // At the moment, we store the whole covered interval of all supporting links.  This makes it easier to compute the
    // reads that we want to separate against later but a better implementation would store the individual alignment
    // intervals in an interval tree.
    for (auto const & el : store.alignedReadStore)
    {
        seqan::BamAlignmentRecord const & record = bamRecords[el.readId];
        if (hasFlagRC(record))
            continue;  // Only consider "left-to-right" links.
        // Get handle to other aligned read entry.
        auto const & otherEl = store.alignedReadStore[mateIdx[el.readId]];
        if (el.contigId != fromID || otherEl.contigId != toID)
            continue;  // cannot support link
        // Compute begin and end position for both alignments.
        auto range = alignmentRange(el);
        auto otherRange = alignmentRange(otherEl);
        // Compute inferred distance.
        int delta = psg.contig[psg.nodes[el.contigId]].length - range.first + otherRange.second;
        int distance = libraryInfo.median - delta;
        int distanceSD = libraryInfo.stdDev;
        // Check whether fits lenMean.
        if (abs(lenMean - distance) > MULT * distanceSD)
            continue;  // does not support this distance

        clusters.front().insert(el.readId);
        clusters.front().insert(otherEl.readId);

        if (support == 0u)
        {
            rangeFrom = range;
            rangeTo = otherRange;
        }
        else
        {
            rangeFrom.first = std::min(rangeFrom.first, range.first);
            rangeFrom.second = std::max(rangeFrom.second, range.second);
            rangeTo.first = std::min(rangeTo.first, otherRange.first);
            rangeTo.second = std::max(rangeTo.second, otherRange.second);
        }
        ++support;
    }

    SEQAN_CHECK(support > 0, "Should have support fromID=%u, toID=%u", fromID, toID);

    // Collect reads ids that overlap 100*OVERLAP % with range and add them to the second cluster.
    for (auto const & el : store.alignedReadStore)
    {
        if (el.contigId != fromID && el.contigId != toID)
            continue;  // Ignore.
        seqan::BamAlignmentRecord const & record = bamRecords[el.readId];
        if (el.contigId == fromID && hasFlagRC(record))
            continue;  // only interested when pointing to the right
        else if (el.contigId == toID && !hasFlagRC(record))
            continue;  // only interested when pointing to the left
        auto range = alignmentRange(el);
        auto getOverlapLen = [](int begin0, int end0, int begin1, int end1) {
            return std::max(0, std::min(end0, end1) - std::min(begin0, begin1));
        };
        int overlapLen = 0;
        if (el.contigId == fromID)
            overlapLen = getOverlapLen(range.first, range.second, rangeFrom.first, rangeFrom.second);
        else if (el.contigId == toID)
            overlapLen = getOverlapLen(range.first, range.second, rangeTo.first, rangeTo.second);
        if (overlapLen > OVERLAP * (range.second - range.first))
            clusters.back().insert(el.readId);
    }

    if (logging)
    {
        std::cerr << "pair-based cluster (fromID= " << fromID << ", toID=" << toID << ")\n";
        for (unsigned i = 0; i < clusters.size(); ++i)
            for (auto readID : clusters[i])
                std::cerr << readID << " -> " << i << "\n";
    }

    handler(clusters);
}
*/

std::map<unsigned, int> PairBasedSeparator::splitBecauseOfBranch(
        unsigned fromID, unsigned toID0, unsigned toID1) const
{
    if (logging)
        std::cerr << "Splitting because of branch " << fromID << " -> {" << toID0 << ", " << toID1 << "}\n";
    std::map<unsigned, int>  result;

    for (auto const & el : store.alignedReadStore)
    {
        // Get handle to other aligned read entry.
        auto const & otherEl = store.alignedReadStore[mateIdx[el.readId]];
        if (el.contigId == fromID && otherEl.contigId == toID0)
        {
            result[el.readId] = 0;
            result[otherEl.readId] = 0;
        }
        else if (el.contigId == fromID && otherEl.contigId == toID1)
        {
            result[el.readId] = 1;
            result[otherEl.readId] = 1;
        }
    }

    if (logging)
        for (auto pair : result)
            std::cerr << "\t" << pair.first << " -> " << pair.second << "\n";

    return result;
}

/*
void PairBasedSeparator::searchForConflictingOverlaps(PairSeparatorGraph const & psg,
                                                      THandler handler) const
{
    auto const MULT = 3;
    auto const MAX_ERROR_RATE = 0.02;

    for (lemon::SmartDigraph::ArcIt arc(psg.graph); arc != lemon::INVALID; ++arc)
    {
        auto from = psg.contig[psg.graph.source(arc)];
        auto to = psg.contig[psg.graph.target(arc)];
        auto link = psg.link[arc];
        if (link.lenMean + MULT * link.lenSD > 0)
            continue;  // no overlap necessary

        int band = MULT * link.lenSD;
        if (overlapExists(store.contigStore[from.contigID].seq, store.contigStore[to.contigID].seq,
                          -link.lenMean, band, MAX_ERROR_RATE, logging))
            continue;  // overlap actually exists

        // In the case that the indicated overlap does not exist, we separate based on reads supporting the overlap.
        splitBecauseOfLink(psg, from.contigID, to.contigID, link.lenMean, handler);
    }
}
*/

std::vector<std::map<unsigned, int>> PairBasedSeparator::searchForConflictingBranches(
        PairSeparatorGraph const & psg) const
{
    unsigned const MULT = 3;

    if (logging)
        std::cerr << "Searching for conflicting branches.\n";

    // The overlapper to use for compatibility checking.
    assembler::OverlapperOptions ovlOptions;
    ovlOptions.overlapErrorRate = 0.03;
    ovlOptions.overlapMinLength = 30;
    ovlOptions.logging = logging;
    assembler::Overlapper overlapper(ovlOptions);

    // Approximately positioned contig.
    struct ContigAtApproximatePos
    {
        ContigAtApproximatePos() = default;

        ContigAtApproximatePos(int pos, int posSD, int len, unsigned contigID) :
                pos(pos), posSD(posSD), len(len), contigID(contigID)
        {}

        int pos { 0 };
        int posSD { 0 };
        int len { 0 };
        unsigned contigID { 0 };
    };

    // Result of this method.
    std::vector<std::map<unsigned, int>> result;

    // For each node: anchor connected contig at node via outgoing/incoming arcs and look for overlaps indicated by mate
    // links that cannot be verified using the sequence.

    // Helper function to handle positioned contigs.
    auto handlePositionedContigs = [&](lemon::SmartDigraph::Node u,
                                       std::vector<ContigAtApproximatePos> const & posCtgs) {
        // Pairwise overlap check.
        for (auto it = posCtgs.begin(); std::next(it) != posCtgs.end(); ++it)
            for (auto it2 = std::next(it); it2 != posCtgs.end(); ++it2)
            {
                auto left = *it, right = *it2;
                if (left.pos > right.pos)
                    std::swap(left, right);
                if (left.pos - MULT * left.posSD + left.len < right.pos + MULT * right.posSD)
                    continue;  // no overlap possible
                int overlapLen = left.pos + left.len - right.pos;
                if (overlapLen < ovlOptions.overlapMinLength)
                    continue;  // overlap not significant enough
                int diagonal = (int)(length(store.contigStore[left.contigID].seq)) - overlapLen;
                int band = MULT * (left.posSD + right.posSD);
                band = std::max(band, 20);
                if (overlapper.overlapExists(
                            store.contigStore[left.contigID].seq, store.contigStore[right.contigID].seq,
                            diagonal, band))
                    continue;

                // In the case that the indicated overlap between left and right does not exist, we separted based on
                // reads supporting the overlap.
                auto res = splitBecauseOfBranch(psg.contig[u].contigID, left.contigID, right.contigID);
                result.insert(result.end(), res);
            }
    };

    // First for outgoing arcs.
    for (lemon::SmartDigraph::NodeIt u(psg.graph); u != lemon::INVALID; ++u)
    {
        std::vector<ContigAtApproximatePos> posCtgs;
        for (lemon::SmartDigraph::OutArcIt arc(psg.graph, u); arc != lemon::INVALID; ++arc)
            posCtgs.push_back(ContigAtApproximatePos(
                    psg.link[arc].lenMean,
                    psg.link[arc].lenSD,
                    psg.contig[psg.graph.target(arc)].length,
                    psg.contig[psg.graph.target(arc)].contigID));
        if (posCtgs.size() > 1u)
            handlePositionedContigs(u, posCtgs);
    }
    // Then for incoming arcs.
    for (lemon::SmartDigraph::NodeIt u(psg.graph); u != lemon::INVALID; ++u)
    {
        std::vector<ContigAtApproximatePos> posCtgs;
        for (lemon::SmartDigraph::InArcIt arc(psg.graph, u); arc != lemon::INVALID; ++arc)
            posCtgs.push_back(ContigAtApproximatePos(
                    psg.link[arc].lenMean,
                    psg.link[arc].lenSD,
                    psg.contig[psg.graph.source(arc)].length,
                    psg.contig[psg.graph.source(arc)].contigID));
        if (posCtgs.size() > 1u)
            handlePositionedContigs(u, posCtgs);
    }
                                                                                             
    return result;
}

std::vector<scaffolder::MateLink> PairBasedSeparator::computeMateLinks(
        rep_solv::MateInfos const & mateInfos, PairSeparatorGraph const & out)
{
    // Minimal support for each link and minimal coverage to require for contigs.
    unsigned const MIN_SUPPORT = 4;  // TODO(holtgrew): Get from configuration?
    unsigned const MIN_COVERAGE = 2 * MIN_SUPPORT;

    std::vector<scaffolder::MateLink> result;

    for (auto const & mateInfo : mateInfos.records)
    {
        if (out.contig[out.nodes[mateInfo.leftID]].numReads < (int)MIN_COVERAGE ||
            out.contig[out.nodes[mateInfo.rightID]].numReads < (int)MIN_COVERAGE)
            continue;  // skip with too low support
        int delta = out.contig[out.nodes[mateInfo.leftID]].length - mateInfo.beginPos + mateInfo.endPos;
        int distance = mateInfos.libraries[mateInfo.libraryID].mean - delta;
        int distanceSD = mateInfos.libraries[mateInfo.libraryID].sd;
        result.push_back(scaffolder::MateLink(
                mateInfo.leftID, mateInfo.rightID,
                scaffolder::MateEdgeLabel(distance, distanceSD)));
    }

    return result;
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function performPairBasedSeparation()
// ----------------------------------------------------------------------------

// Perform separation based on pair information.  FeatureMap and the two read sets.

void performPairBasedSeparation(FeatureMap & featureMap,
                                FeatureReadSet & readSet,
                                FeatureReadSet & atomicReadSet,
                                TFragmentStore & fragStore,
                                rep_solv::MateInfos const & mateInfos)
{
    PairBasedSeparator sep(fragStore, mateInfos, /*logging=*/false);
    PairBasedSeparator::TResult result = sep.run();

    // Process single read-to-feature assignments into global "new features for read" map.
    std::map<unsigned, WeightedFeatureVector> newFeatures;  // readID -> [(featureID, value)]
    for (unsigned i = 0; i < result.size(); ++i)
    {
        unsigned featureID = featureMap.size();
        featureMap.insert(FeatureDescription(featureID, FeatureDescription::LINK));
        for (auto pair : result[i])
            newFeatures[pair.first].insert(WeightedFeature(featureID, pair.second));
    }

    for (auto & read : readSet)
        for (auto readID : read.subReads)
            if (newFeatures.count(readID))
                read.features.mergeWithThis(newFeatures.find(readID)->second);
    for (auto & read : atomicReadSet)
        for (auto readID : read.subReads)
            if (newFeatures.count(readID))
                read.features.mergeWithThis(newFeatures.find(readID)->second);

    featureMap.refresh();
}

}  // namespace rep_sep
