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

#include "consensus_builder.h"

#include <random>
#include <algorithm>

#include <seqan/graph_types.h>
#include <seqan/graph_msa.h>
#include <seqan/realign.h>

#include "asm/contig_graph.h"
#include "asm/overlapper.h"
#include "asm/frag_store_builder.h"
#include "asm/ag_to_frag_store.h"

namespace assembler {

namespace {  // anonymous namespace
}  // anonymous namespace

// --------------------------------------------------------------------------
// Class ConsensusBuilderImpl
// --------------------------------------------------------------------------

class ConsensusBuilderImpl
{
public:
    ConsensusBuilderImpl(Options const & options) : options(options)
    {}

    std::unique_ptr<TFragmentStore> run(
            ContigGraph const & cg, seqan::StringSet<seqan::Dna5String> const & seqs) const;

private:

    // Glorified pair<unsigned, int> with more verbose name for readID and offset.
    struct AlignmentTarget
    {
        AlignmentTarget() = default;
        AlignmentTarget(unsigned readID, int offset) : readID(readID), offset(offset) {}

        unsigned readID { (unsigned)-1 };
        int offset { 0 };
    };

    // Glorified triple for overlap computation.
    struct PositionedRead
    {
        PositionedRead() = default;
        PositionedRead(int beginPos, int endPos, unsigned readID) :
            beginPos(beginPos), endPos(endPos), readID(readID)
        {}

        int beginPos { 0 };
        int endPos { 0 };
        unsigned readID { (unsigned)-1 };

        bool operator<(PositionedRead const & other) const
        {
            return (std::make_tuple(beginPos, endPos, readID) <
                    std::make_tuple(other.beginPos, other.endPos, other.readID));
        }
    };

    // Compute pairwise overlaps for the reads of each contig in cg.
    //
    // In case that there are regions higher than options.startConsensusCompression, we add subsample these reads and
    // write out an entry (other read id, offset) to use as a position marker in the resulting MSA.
    void computeOverlaps(
            TFragments & fragments,
            seqan::String<int> & scores,
            seqan::Graph<seqan::Undirected<double>> & distances,
            std::map<unsigned, AlignmentTarget> & alignTo,
            seqan::StringSet<seqan::Dna5String> const & seqs,
            ContigGraph const & cg) const;

    // Integrate reads from alignTo into store using alignTo.
    void integrateReads(TFragmentStore & store,
                        std::map<unsigned, AlignmentTarget> const & alignTo) const;

    // Sample pos reads, writing out non-copied reads to alignTo, return filtered.
    std::vector<PositionedRead> samplePosReads(std::map<unsigned, AlignmentTarget> & alignTo,
                                               std::vector<PositionedRead> posReads) const;

    // Configuration to use for consensus building.
    Options options;
};


std::vector<ConsensusBuilderImpl::PositionedRead>
ConsensusBuilderImpl::samplePosReads(std::map<unsigned, ConsensusBuilderImpl::AlignmentTarget> & alignTo,
                                     std::vector<ConsensusBuilderImpl::PositionedRead> posReads) const
{
    int seed = 0;  // initialize from data
    // Ids of reads to keep, must keep first as anchor.
    std::set<unsigned> keep;

    // Build array of coverages.
    std::vector<int> counts;
    for (auto const & posRead : posReads)
    {
        seed += posRead.beginPos;
        counts.resize(std::max((int)counts.size(), posRead.endPos), 0);
    }

    // Shuffle input.
    std::mt19937 mt(seed);
    std::shuffle(posReads.begin(), posReads.end(), mt);

    // Now, place reads in order of shuffling in counts.  Build set of read ids to keep.
    for (auto const & posRead : posReads)
    {
        bool keepThis = false;
        for (int pos = posRead.beginPos; pos < posRead.endPos; ++pos)
            keepThis = (counts.at(pos)++ < 50/*options.startConsensusCompression*/) && keepThis;
        if (keepThis)
            keep.insert(posRead.readID);
    }

    // Build result and fill alignTo, iterate over posReads sorted by coordinate again.
    std::sort(posReads.begin(), posReads.end());
    keep.insert(posReads.front().readID);  // must keep as anchor
    std::vector<ConsensusBuilderImpl::PositionedRead> result;
    auto anchor = posReads.front();
    for (auto posRead : posReads)
        if (keep.count(posRead.readID))
        {
            anchor = posRead;
            result.push_back(posRead);
        }
        else
        {
            alignTo[posRead.readID] = AlignmentTarget(anchor.readID, posRead.beginPos - anchor.beginPos);
        }

    return result;
}

void ConsensusBuilderImpl::integrateReads(
    TFragmentStore & store,
    std::map<unsigned, ConsensusBuilderImpl::AlignmentTarget> const & alignTo) const
{
    // TODO(holtgrew): Perform pairwise alignment against consensus instead of relying on subsequence realignment.

    // IDs of the contigs that are removed since there are no more singletons on them.
    std::set<unsigned> removeContigs;

    // Sort aligned reads by read id.
    sortAlignedReads(store.alignedReadStore, seqan::SortReadId());

    // Integrate all reads from alignTo.
    for (auto pair : alignTo)
    {
        // Get shortcuts to relevant store entries and "pair" members.
        unsigned readID = pair.first;
        unsigned targetID = pair.second.readID;
        auto & alignedRead = store.alignedReadStore[readID];
        auto const & otherAlignedRead = store.alignedReadStore[targetID];
        SEQAN_CHECK(alignedRead.readId == readID, "Must correspond.");

        // Update aligned read store element.
        removeContigs.insert(alignedRead.contigId);
        alignedRead.contigId = otherAlignedRead.contigId;
        alignedRead.beginPos = otherAlignedRead.beginPos + pair.second.offset;
        alignedRead.endPos = alignedRead.beginPos + length(store.readSeqStore[readID]);
        clear(alignedRead.gaps);
    }

    // Update contig store and aligned read store contig ids.
    std::map<unsigned, unsigned> oldToNew;
    unsigned newCount = 0;
    for (unsigned contigID = 0; contigID < length(store.contigStore); ++contigID)
        if (!removeContigs.count(contigID))
        {
            using std::swap;
            if (contigID != newCount)
                swap(store.contigStore[contigID], store.contigStore[newCount]);
            oldToNew[contigID] = newCount++;
        }
    resize(store.contigStore, newCount);
    for (auto & el : store.alignedReadStore)
    {
        SEQAN_CHECK(oldToNew.count(el.contigId), "Must not be on removed contig.");
        el.contigId = oldToNew.find(el.contigId)->second;
    }
}

void ConsensusBuilderImpl::computeOverlaps(
            TFragments & fragments,
            seqan::String<int> & scores,
            seqan::Graph<seqan::Undirected<double>> & distances,
            std::map<unsigned, ConsensusBuilderImpl::AlignmentTarget> & alignTo,
            seqan::StringSet<seqan::Dna5String> const & seqs,
            ContigGraph const & cg) const
{
    // Configure overlapper, no limit in error rate or minimum length.
    OverlapperOptions ovlOptions;
    ovlOptions.logging = (options.verbosity >= 3);
    ovlOptions.overlapErrorRate = options.readOverlapErrorCount;
    ovlOptions.overlapMinLength = 0;  // accept all overlaps
    Overlapper overlapper(ovlOptions);

    // TODO(holtgrew): Build graph contig graph node-wise.

    for (auto const & node : cg.node)
    {
        auto const & contig = cg.contig[node];
        if (contig.posReads.empty())
            continue;  // no positioned reads
        // Build begin/end intervals.
        std::vector<PositionedRead> posReads;
        std::transform(contig.posReads.begin(), contig.posReads.end(), std::back_inserter(posReads),
                       [&](Contig::PositionedRead posRead) {
                           return PositionedRead(posRead.beginPos,
                                                 posRead.beginPos + posRead.length,
                                                 posRead.readID);
                       });
        posReads = samplePosReads(alignTo, posReads);

        if (options.verbosity >= 3)
        {
            std::cerr << "COORDS\n";
            for (auto const & posRead : posReads)
                std::cerr << posRead.readID << ":\t" << posRead.beginPos << '-' << posRead.endPos << '\n';
        }
        // Add all pairwise overlaps.
        Overlap ovl;
        TFragments frags;
        for (auto it0 = posReads.begin(); std::next(it0) != posReads.end(); ++it0)
            for (auto it1 = std::next(it0); it1 != posReads.end(); ++it1)
            {
                if (it1->beginPos >= it0->endPos)
                    break;  // no more overlap possible, sorted by begin position
                int d0 = it1->beginPos - it0->beginPos;
                int d1 = it0->endPos - it1->beginPos;
                int band = 0;
                if (options.readOverlapErrorCount >= 0)
                    band = options.readOverlapErrorCount;
                if (options.readOverlapErrorRate > 0)
                    band = std::max(band, (int)(ceil(d1 * options.readOverlapErrorRate / 100.0)));
                int uDiag = d0 + band;
                int lDiag = d0 - band;
                unsigned id0 = it0->readID;
                unsigned id1 = it1->readID;
                OverlapCandidate candidate(id0, id1, lDiag, uDiag);
                if (options.verbosity >= 3)
                    std::cerr << "VERIFYING CANDIDATE\t" << candidate << "\n";
                // TODO(holtgrew): Use OverlapStore as cache.
                if (overlapper.computeOverlap(ovl, frags, seqs[id0], seqs[id1], candidate))
                {
                    append(fragments, frags);
                    int ovlScore = ovl.length() - ovl.errors;  // pseudo-score
                    resize(scores, length(fragments), ovlScore);
                    int quality = (ovl.length() > 0) ? (ovl.length() - ovl.errors) * 100 / ovl.length() : 0;
                    addEdge(distances, idToPosition(seqs, ovl.seq0), idToPosition(seqs, ovl.seq1), quality);
                    if (options.verbosity >= 3)
                        std::cerr << "VERIFIED OVERLAP " << ovl << "\n";
                }
                else
                {
                    if (options.verbosity >= 3)
                        std::cerr << "DISCARDING OVERLAP CANDIDATE " << candidate << "\n";
                }
            }
    }
}

std::unique_ptr<TFragmentStore> ConsensusBuilderImpl::run(
        ContigGraph const & cg, seqan::StringSet<seqan::Dna5String> const & seqsC) const
{
    // The alignments as required by consensus MSA module.
    TFragments fragments;
    seqan::String<int> scores;
    seqan::Graph<seqan::Undirected<double>> distances;
    _resizeWithRespectToDistance(distances, length(seqsC));

    // Obtain alignments according to contig graph.
    std::map<unsigned, AlignmentTarget> alignTo;
    computeOverlaps(fragments, scores, distances, alignTo, seqsC, cg);

    // Build the alignment graph.
    // TODO(holtgrew): Fix const-correctness!
    seqan::StringSet<seqan::Dna5String> & seqs = const_cast<seqan::StringSet<seqan::Dna5String> &>(seqsC);
    seqan::Score<int, seqan::Simple> msaScoringScheme(2, -6, -4, -9);  // TODO(holtgrew): Get from options!
    typedef seqan::StringSet<seqan::Dna5String, seqan::Dependent<> > TDepReadSet;
    typedef seqan::Graph<seqan::Alignment<TDepReadSet, unsigned> > TInGraph;
    TDepReadSet depSeqs(seqs);
    TInGraph inGraph(depSeqs);
    buildAlignmentGraph(fragments, scores, inGraph, msaScoringScheme, seqan::ReScore());

    if (options.verbosity >= 2)
        std::cerr << "Building alignment graph\n"
                  << "  # fragments: " << length(fragments) << "\n"
                  << "  # seqs: " << length(seqs) << "\n";

    // Perform triplet library extension.
    if (length(seqs) < (unsigned)options.texStopCount)
        tripletLibraryExtension(inGraph);

    seqan::Graph<seqan::Tree<double> > guideTree;
    seqan::Graph<seqan::Undirected<double> > dCopy(distances);
    upgmaTree(dCopy, guideTree);

    // Perform progressive alignment.
    if (options.verbosity >= 2)
        std::cerr << "progressive alignment...\n";
    TInGraph graph(depSeqs);
    assignStringSet(graph, stringSet(inGraph));
    seqan::progressiveAlignment(inGraph, guideTree, graph);

    // Build alignment graph.
    if (options.verbosity >= 2)
        std::cerr << "build AG..\n";
    std::unique_ptr<TFragmentStore> result(new TFragmentStore);
    seqan::String<unsigned> component;
    seqan::String<unsigned> order;
    std::map<unsigned, unsigned> componentLength;
    SEQAN_ASSERT_NOT(empty(graph));
    bool b = convertAlignment(graph, component, order, componentLength);
    (void) b;
    SEQAN_ASSERT(b);
    unsigned numComponents = length(order);
    if (options.verbosity >= 2)
        std::cerr << "AG\n"
                  << "  numVertices = " << numVertices(graph) << "\n"
                  << "  numEdges = " << numEdges(graph) << "\n"
                  << "  componentLength.size() = " << componentLength.size() << "\n"
                  << "  numComponents = " << numComponents << "\n";

    if (options.verbosity >= 2)
        std::cerr << "build store..\n";
    alignmentGraphToFragmentStore(*result, graph, distances, component, order, numComponents,
                                  /*logging=*/(options.verbosity >= 3));

    if (options.verbosity >= 2)
    {
        std::cerr << "CONSENSUS RESULT\n";
        for (auto const & el : result->alignedReadStore)
            std::cerr << "contigID=" << el.contigId
                      << ", beginPos=" << std::min(el.beginPos, el.endPos)
                      << ", readID=" << el.readId
                      << ", seq=" << result->readSeqStore[el.readId]
                      << "\n";
        std::cerr << "\n";
    }

    if (options.verbosity >= 2)
        std::cerr << "realigning...\n";
    for (unsigned i = 0; i < length(result->contigStore); ++i)
        reAlignment(*result, i, 1, 10, false);

    if (!alignTo.empty())
    {
        integrateReads(*result, alignTo);
        for (unsigned i = 0; i < length(result->contigStore); ++i)
            reAlignment(*result, i, 1, 10, false);
    }

    return result;
}

// --------------------------------------------------------------------------
// Class ConsensusBuilder
// --------------------------------------------------------------------------

ConsensusBuilder::ConsensusBuilder(Options const & options) :
        impl(new ConsensusBuilderImpl(options))
{}

ConsensusBuilder::~ConsensusBuilder()
{}

std::unique_ptr<TFragmentStore> ConsensusBuilder::run(
        ContigGraph const & cg,
        seqan::StringSet<seqan::Dna5String> const & seqs) const
{
    return impl->run(cg, seqs);
}

}  // namespace assembler
