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

#include <iostream>
#include <stdexcept>

#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>

#include "asm/consensus_builder.h"
#include "asm/contig_graph.h"
#include "asm/frag_store.h"
#include "asm/overlapper.h"
#include "asm/store_utils.h"
#include "asm/unitigger.h"

#include "rep_sep/feature_map.h"
#include "rep_sep/read_set.h"
#include "rep_sep/rep_sep.h"
#include "rep_sep/read_correction.h"
#include "rep_sep/rep_sep_options.h"

#include "rep_solv/contig_graph.h"
#include "rep_solv/dead_branch.h"
#include "rep_solv/directed_tree.h"
#include "rep_solv/graph_expansion.h"
#include "rep_solv/options.h"
#include "rep_solv/path_enumeration.h"
#include "rep_solv/scaffolding.h"

// TODO(holtgrew): Move these into their own module.

static const char * KEY_ORIG_REF = "oR";
static const char * KEY_ORIG_CIGAR = "oC";
static const char * KEY_BIRTH_STEP = "mS";

// ----------------------------------------------------------------------------
// Function eraseIf()
// ----------------------------------------------------------------------------

template <typename TValue, typename TPred>
void eraseIf(std::vector<TValue> & v, TPred pred)
{
    v.erase(std::remove_if(v.begin(), v.end(), pred), v.end());
}

template <typename TValue, typename TPred>
void eraseIf(seqan::String<TValue> & v, TPred pred)
{
    resize(v, std::remove_if(begin(v, seqan::Standard()), end(v, seqan::Standard()), pred) - begin(v, seqan::Standard()));
}


// TODO(holtgrew): Implement pair-based separation again.

void loadSeqsFromSam(seqan::StringSet<seqan::Dna5String> & seqs,
                     std::vector<seqan::BamAlignmentRecord> & records,
                     char const * path)
{
    seqan::BamStream inSam(path);
    if (!isGood(inSam))
        throw std::runtime_error("Problem reading from SAM file\n");

    while (!atEnd(inSam))
    {
        seqan::BamAlignmentRecord record0, record1;
        if (readRecord(record0, inSam) != 0 || readRecord(record1, inSam) != 0)
            throw std::runtime_error("Problem reading from SAM file\n");

        if (record0.qName != record1.qName)
            throw std::runtime_error("SAM file should contain read0/1,read0/2, read1/1,read1/2, ...");

        if (hasFlagUnmapped(record0) && hasFlagUnmapped(record1))
            continue;  // Skip orphans.

        if (!hasFlagUnmapped(record0) && hasFlagUnmapped(record1) &&(hasFlagRC(record0) == hasFlagRC(record1)))
        {
            reverseComplement(record1.seq);
            record1.flag ^= seqan::BAM_FLAG_RC;
        }
        else if (!hasFlagUnmapped(record1) && hasFlagUnmapped(record0) && (hasFlagRC(record0) == hasFlagRC(record1)))
        {
            reverseComplement(record0.seq);
            record0.flag ^= seqan::BAM_FLAG_RC;
        }

        if (!hasFlagRC(record0))
        {
            records.push_back(record0);  // --->
            records.push_back(record1);  // <---
            appendValue(seqs, record0.seq);
            appendValue(seqs, record1.seq);
        }
        else
        {
            records.push_back(record1);  // --->
            records.push_back(record0);  // <---
            appendValue(seqs, record1.seq);
            appendValue(seqs, record0.seq);
        }
    }
}

// Fill readNameStore of TFragmentStore etc.
void setUpStore(TFragmentStore & store)
{
    for (unsigned readID = length(store.readNameStore); readID < length(store.readSeqStore); ++readID)
    {
        std::stringstream ss;
        ss << "read_" << readID;
        appendValue(store.readNameStore, ss.str().c_str());
    }
}

// Augment ContigGraph with edges from s/to t.
//
// records[id] contain the BamAlignmentRecord for the read with the given id.

void addSTLinks(rep_solv::ContigGraph & graph,
                std::vector<seqan::BamAlignmentRecord> const & records,
                bool debug)
{
    std::set<unsigned> s, t;  // connected to s/t

    for (unsigned readID = 0; readID < records.size(); ++readID)
    {
        seqan::BamAlignmentRecord & record = const_cast<seqan::BamAlignmentRecord &>(records[readID]);
        seqan::BamTagsDict tagsDict(record.tags);
        unsigned idx = 0;
        if (!findTagKey(idx, tagsDict, KEY_ORIG_REF))
            continue;  // ignore, no original alignment
        // Ignore if has clipping.
        seqan::CharString cigar;
        if (findTagKey(idx, tagsDict, KEY_ORIG_CIGAR) && extractTagValue(cigar, tagsDict, idx) &&
            (std::find(begin(cigar, seqan::Standard()), end(cigar, seqan::Standard()), 'S') !=
             end(cigar, seqan::Standard())))
                continue;  // Skip, used to be clipped and is shadow now.

        if (debug)
            std::cerr << "s-t link\t" << record.qName << "\t" << cigar << "\t"
                      << "rc=" << hasFlagRC(record) << "\n";

        if (!hasFlagRC(record))
        {
            for (auto contigID : graph.readToContigs[readID])
            {
                if (s.count(contigID))
                    continue;  // already connected
                s.insert(contigID);
                graph.addEdge(graph.s, graph.node[contigID], rep_solv::EdgeLabel(rep_solv::EdgeLabel::SOURCE, contigID));
            }
        }
        else
        {
            for (auto contigID : graph.readToContigs[readID])
            {
                if (t.count(contigID))
                    continue;  // already connected
                t.insert(contigID);
                graph.addEdge(graph.node[contigID], graph.t, rep_solv::EdgeLabel(contigID, rep_solv::EdgeLabel::TARGET));
            }
        }
    }
}

// Build vector with the steps numbers that the reads were added to the graph.

std::vector<int> buildReadBirthVector(std::vector<seqan::BamAlignmentRecord> const & records)
{
    std::vector<int> result(records.size(), 0);

    for (unsigned recordID = 0; recordID < records.size(); ++recordID)
    {
        seqan::BamAlignmentRecord & record = const_cast<seqan::BamAlignmentRecord &>(records[recordID]);

        seqan::BamTagsDict tagsDict(record.tags);
        unsigned idx = 0;
        if (!findTagKey(idx, tagsDict, KEY_BIRTH_STEP))
            continue;  // ignore
        if (!extractTagValue(result[recordID], tagsDict, idx))
            continue;  // could not extract but ignore this
    }

    return result;
}

// Print the store haplotype by haplotype.

void printHaplotypes(std::ostream & out,
                     TFragmentStore const & store,
                     rep_sep::FeatureReadSet const & readSet)
{
    unsigned lastPrinted = (unsigned)-1;
    for (unsigned idx = 0; idx < readSet.reads.size(); ++idx)
    {
        auto const & read = readSet.reads[idx];
        unsigned contigID = read.contigID;
        std::set<unsigned> readIDs(read.subReads.begin(), read.subReads.end());

        auto tmp = store;

        if (lastPrinted != contigID)
        {
            out << "ORIGINAL\n";
            printStore(out, tmp, contigID);
            lastPrinted = contigID;
        }

        eraseIf(tmp.alignedReadStore, [&](decltype(tmp.alignedReadStore[0]) el) {
                return (el.contigId != contigID || !readIDs.count(el.readId));
            });

        out << "HAPLOTYPE\n";
        printStore(out, tmp, contigID);
    }
}

void expandContigs(seqan::StringSet<seqan::Dna5String> & expandedContigs,
                   TFragmentStore const & fragStore,
                   rep_sep::FeatureReadSet const & readSet,
                   rep_sep::FeatureMap const & featureMap,
                   std::vector<unsigned> const & oldToNew)
{
    // Build mapping from contig IDs without removed to contigIDs before removal.
    std::vector<unsigned> newToOld(length(fragStore.contigStore), (unsigned)-1);
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
        if (oldToNew.at(i) != (unsigned)-1)
            newToOld.at(oldToNew[i]) = i;

    // For each super read in readSet, compute the corresponding haplotype.
    clear(expandedContigs);
    resize(expandedContigs, readSet.size());
    for (unsigned contigID = 0; contigID < readSet.size(); ++contigID)  // output contigID
    {
        // Shortcuts to required values.
        auto & out = expandedContigs[contigID];
        auto const & read = readSet.reads[contigID];
        auto const & cEl = fragStore.contigStore[newToOld[read.contigID]];
        auto const & inSeq = cEl.seq;
        auto const & inAnchors = cEl.gaps;

        // Obtain range of features from featureMap.
        auto featureRange = featureMap.contigFeatures(read.contigID);
        auto itF = featureRange.first;
        auto const & readFeatures = read.features;

        // Typedefs for output sequence construction.
        typedef typename TFragmentStore::TContigStore const TContigStore;
        typedef typename seqan::Value<TContigStore>::Type TContig;
        typedef typename TFragmentStore::TContigSeq const TContigSeq;
        typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<TContig::TGapAnchors>> TContigGaps;

        // Build output sequence from FragmentStore's contig gaps and contig features.
        TContigGaps contigGaps(inSeq, inAnchors);
        auto itBegin = begin(contigGaps, seqan::Standard());
        auto itEnd = end(contigGaps, seqan::Standard());
        int pos = 0;
        for (auto it = itBegin; it != itEnd; ++it, ++pos)
        {
            // Not on a feature, copy over consensus char.
            if (itF == featureRange.second || itF->pos != pos)
            {
                if (itF != featureRange.second)
                    SEQAN_ASSERT_GT(itF->pos, pos);
                if (!isGap(it))
                    appendValue(out, *it);
                continue;
            }

            // If we reach here then we hit a feature.  Apply feature value of read if present.
            int const GAP = 5;  // feature value for a gap
            auto itV = readFeatures.find(itF->id);
            if (itV != readFeatures.end())
            {
                // Feature found in read, append feature value, ignoring gaps.
                if (itV->value != GAP)
                    appendValue(out, seqan::Dna5(itV->value));
            }
            else
            {
                // Feature not found in read, append consensus value if not gap.
                if (!isGap(it))
                    appendValue(out, *it);
            }

            // Advance over feature.
            ++itF;
            if (itF != featureRange.second)
                SEQAN_ASSERT_GT(itF->pos, pos);
        }
    }
}

int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: asm IN.sam [DEBUG]\n";
        return 1;
    }
    bool debug = (argc >= 3);

    // Load sequences.
    std::cerr << "Loading sequences...\n";
    seqan::StringSet<seqan::Dna5String> seqs;
    std::vector<seqan::BamAlignmentRecord> records;
    loadSeqsFromSam(seqs, records, argv[1]);
    if (true || debug)
    {
        std::cerr << "SEQUENCES\n";
        for (unsigned i = 0; i < length(seqs); ++i)
            std::cerr << i << '\t' << seqs[i] << '\n';
    }

    // ------------------------------------------------------------------------
    // Compute overlaps.
    // ------------------------------------------------------------------------

    std::cerr << "Computing overlaps...\n";
    assembler::OverlapCandidateGeneratorOptions candOptions;
    // candOptions.logging = debug;
    assembler::OverlapCandidateGenerator candGen(candOptions);
    assembler::OverlapStore ovlStore;
    assembler::OverlapperOptions ovlOptions;
    ovlOptions.logging = false && debug;
    assembler::Overlapper overlapper(ovlOptions);
    auto verify = [&](assembler::OverlapCandidate const & cand) {
        overlapper.computeOverlap(ovlStore, seqs[cand.seq0], seqs[cand.seq1], cand);
    };
    candGen.run(verify, seqs);
    ovlStore.refresh();
    std::cerr << "  Computed " << ovlStore.size() << " overlaps\n";

    // ------------------------------------------------------------------------
    // Unitigging.
    // ------------------------------------------------------------------------

    std::cerr << "Performing Unitigging...\n";
    auto cgPtr = assembler::Unitigger(assembler::Unitigger::MIRA_LIKE, true || debug).run(seqs, ovlStore.overlaps());
    if (debug)
        cgPtr->print(std::cerr);
    std::cerr << "Building Consensus...\n";
    assembler::Options asmOptions;
    // asmOptions.verbosity = 2;
    assembler::ConsensusBuilder consensusBuilder(asmOptions);
    auto fsPtr = consensusBuilder.run(*cgPtr, seqs);

    // There must be a 1:1 correspondence between contig in ContigGraph and FragmentStore.
    SEQAN_CHECK(length(fsPtr->contigStore) == cgPtr->node.size(), "Must correspond!");
    // Reconstruct this correspondence and permute the contigStore.
    std::vector<unsigned> contigPermutation(length(fsPtr->contigStore), (unsigned)-1);
    for (auto const & el : fsPtr->alignedReadStore)
    {
        unsigned cgContigID = cgPtr->contig[cgPtr->readToContig[el.readId]].id;
        if (contigPermutation.at(el.contigId) == (unsigned)-1)
            contigPermutation[el.contigId] = cgContigID;
        else
            SEQAN_ASSERT_EQ(contigPermutation[el.contigId], cgContigID);
    }
    permuteStoreContigs(*fsPtr, contigPermutation);

    if (true || debug)
    {
        std::cerr << "ContigGraph from Assembly\n";
        cgPtr->print(std::cerr);

        std::cerr << "Resulting fragment store\n";
        printStore(std::cerr, *fsPtr);

        // Print store to output.
        {
            for (unsigned i = 0; i < length(fsPtr->contigStore); ++i)
            {
                std::cout << ">contig_" << i << "\n";
                std::cout << fsPtr->contigStore[i].seq << "\n";
            }
        }
    }

    // Fill read name store for FragmentStore, ConsensusBuilder only fills in contig and read seqs.
    setUpStore(*fsPtr);

    // ------------------------------------------------------------------------
    // Perform Repeat Separation.
    // ------------------------------------------------------------------------

    // Perform repeat separation (enumerates separating columns with Tammi's test and then uses clique enumeration to
    // separate the repeats.
    rep_sep::ReadSeparatorOptions rsOptions;
    // rsOptions.verbosity = 2;
    rep_sep::FeatureReadSet preReadSet;  // pre since we will remove some contigs below
    rep_sep::FeatureReadSet atomicReadSet;  // features for each atomic read, can be used to fix counts
    rep_sep::FeatureMap featureMap;
    separateRepeats(preReadSet, atomicReadSet, featureMap, *fsPtr, rsOptions);

    if (true || debug)
    {
        std::cerr << "Atomic read set\n";
        atomicReadSet.print(std::cerr);
        std::cerr << "(pre) Read set after Separation\n";
        preReadSet.print(std::cerr);
        std::cerr << "Feature Map\n";
        featureMap.print(std::cerr);

        printHaplotypes(std::cerr, *fsPtr, preReadSet);
    }

    // ------------------------------------------------------------------------
    // Perform Error Correction
    // ------------------------------------------------------------------------

    auto returnTrue = [&](seqan::BamAlignmentRecord const &) { return true; };
    performReadCorrection(records, *fsPtr, featureMap, returnTrue);

    std::cerr << "Store after error correction\n";
    printStore(std::cerr, *fsPtr);

    // ------------------------------------------------------------------------
    // Perform Repeat Resolution.
    // ------------------------------------------------------------------------

    // Build mate-pair infos.
    sortAlignedReads(fsPtr->alignedReadStore, seqan::SortReadId());
    rep_solv::MateInfos mateInfos;
    SEQAN_ASSERT_EQ(records.size() % 2u, 0u);
    // TODO(holtgrew): Currently strong assumption that mate pair are consecutive in records.
    if (debug)
        std::cerr << "ADDING PE LINKS\n";
    for (unsigned idx = 0; idx < records.size(); idx += 2)
    {
        auto const & recordL = records[idx];
        auto const & recordR = records[idx + 1];
        SEQAN_ASSERT_EQ(fsPtr->alignedReadStore[idx].readId, idx);
        SEQAN_ASSERT_EQ(fsPtr->alignedReadStore[idx + 1].readId, idx + 1);
        auto const & elL = fsPtr->alignedReadStore[idx];
        auto const & elR = fsPtr->alignedReadStore[idx + 1];
        if (elL.contigId == elR.contigId)
            continue;  // skip links on the same contig

        if (debug)
            std::cerr << "LINK FROM\n"
                      << "\trecordL.qName=" << recordL.qName << "\trecordL.flag=" << recordL.flag
                      << "\telL.contigId=" << elL.contigId << "\telL.beginPos=" << elL.beginPos
                      << "\telL.endPos=" << elL.endPos
                      << "\n"
                      << "\trecordR.qName=" << recordR.qName << "\trecordR.flag=" << recordR.flag
                      << "\telR.contigId=" << elR.contigId << "\telR.beginPos=" << elR.beginPos
                      << "\telR.endPos=" << elR.endPos
                      << "\n";

        typedef typename TFragmentStore::TReadPos         TReadPos;
        typedef typename TFragmentStore::TContigStore     TContigStore;
        typedef typename seqan::Value<TContigStore>::Type TContig;
        typedef typename TFragmentStore::TContigSeq       TContigSeq;
        typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
        TContigGaps contigGapsL(fsPtr->contigStore[elL.contigId].seq, fsPtr->contigStore[elL.contigId].gaps);
        TContigGaps contigGapsR(fsPtr->contigStore[elR.contigId].seq, fsPtr->contigStore[elR.contigId].gaps);

        SEQAN_ASSERT_MSG(!hasFlagRC(recordL), "QNAME=%s", toCString(recordL.qName));
        SEQAN_ASSERT_MSG(hasFlagRC(recordR), "QNAME=%s", toCString(recordR.qName));
        rep_solv::MateInfo info(idx, idx + 1,
                                elL.contigId, elR.contigId,
                                toSourcePosition(contigGapsL, std::min(elL.beginPos, elL.endPos)),
                                toSourcePosition(contigGapsR, std::max(elR.beginPos, elR.endPos)),
                                /*libraryID=*/0);
        if (debug)
            std::cerr << "\t" << info << "\n";
        mateInfos.insert(info);
    }
    // Hard-coded library info for now.
    mateInfos.libraries = { rep_solv::LibraryInfo(387.0, 40.55) };

    // Build contig graph for repeat resolution from assembled contig graph.  We call it pre contig graph since we will
    // remove some contigs in our heuristics.
    std::vector<unsigned> contigLengths;
    for (unsigned i = 0; i < length(fsPtr->contigStore); ++i)
        contigLengths.push_back(length(fsPtr->contigStore[i].seq));
    rep_solv::ContigGraph preContigGraph;
    buildContigGraph(preContigGraph, *cgPtr, contigLengths);
    // Augment the graph with s- and t-links.
    addSTLinks(preContigGraph, records, debug);
    if (debug)
    {
        std::cerr << "PRE ContigGraph (with s-/t-links)\n";
        preContigGraph.print(std::cerr);
    }
    // Execute dead branch and directed tree reachability heuristic.
    //
    // For this, we need mate links that we do not want for the graph expansion, however.  Thus, we create a copy here.
    rep_solv::Options repSolvOptions;
    // repSolvOptions.verbosity = 2;
    repSolvOptions.maxMateConflicts = -1;  // infinity
    unsigned stepNo = 1;
    // We build a mapping from contigID to flag that indicates a deletion wish.
    std::vector<bool> removeContigMap(preContigGraph.node.size(), false);
    {
        rep_solv::ContigGraph graphCopy(preContigGraph);
        addPairedEndLinks(graphCopy, mateInfos);
        if (debug)
        {
            std::cerr << "PRE ContigGraph (with s-/t-links, PE links) [" << __LINE__ << "]\n";
            graphCopy.print(std::cerr);
        }

        auto readBirthStep = buildReadBirthVector(records);
        lemon::SmartGraph::NodeMap<bool> doRemove(graphCopy.graph, false);
        unsigned numRemoved = removeDeadBranches(doRemove, graphCopy, stepNo, readBirthStep, repSolvOptions);
        std::cerr << "Dead branch removal: removed " << numRemoved << " vertices\n";
        // Execute directed tree growing heuristic.
        lemon::SmartGraph::NodeMap<bool> reachable(graphCopy.graph, false);
        unsigned numUnreached = directedTreeGrowing(reachable, graphCopy);
        std::cerr << "Directed tree heuristic: " << numUnreached << " unreachable nodes\n";
        // Combine doRemove and reachable into removeContigMap.
        for (lemon::SmartGraph::NodeIt u(graphCopy.graph); u != lemon::INVALID; ++u)
            if (!reachable[u] || doRemove[u])
            {
                removeContigMap.at(graphCopy.contig[u].id) = true;
                if (true || debug)
                    std::cerr << "removing contig " << graphCopy.contig[u].id << "\n";
            }
    }

    // Copy over removeContigMap into map for preContigGraph;
    lemon::SmartGraph::NodeMap<bool> doRemove(preContigGraph.graph, false);
    for (unsigned i = 0; i < removeContigMap.size(); ++i)
        if (removeContigMap[i])
            doRemove[preContigGraph.node[i]] = true;
    // Build a new contig graph without dead branches and branches into the "wrong" direction.
    rep_solv::ContigGraph contigGraph;
    removeContigs(contigGraph, preContigGraph, doRemove);
    // TODO(holtgrew): Also remove contigs with too low support against noise!
    std::cerr << "TODO: remove contigs with too low support against noise!\n";
    rep_sep::FeatureReadSet readSet;
    removeReads(readSet, preReadSet, removeContigMap);
    // Fix feature counts in readSet.
    fixFeatureCounts(readSet, atomicReadSet);
    // Copy out consensus seqs, ignore removed ones.
    seqan::StringSet<seqan::Dna5String> consensusSeqs;
    for (unsigned i = 0; i < length(fsPtr->contigStore); ++i)
        if (!removeContigMap[i])
            appendValue(consensusSeqs, fsPtr->contigStore[i].seq);

    if (debug)
    {
        std::cerr << "ContigGraph (after removing) [" << __LINE__ << "]\n";
        contigGraph.print(std::cerr);
        std::cerr << "Consensus seqs (after removing)\n";
        for (unsigned i = 0; i < length(consensusSeqs); ++i)
            std::cerr << ">" << i << "\n"
                      << consensusSeqs[i] << "\n";
    }
    if (true || debug)
    {
        std::cerr << "(pre) Read set after removing\n";
        readSet.print(std::cerr);
    }

    // Build renumber mapping.
    std::vector<unsigned> oldToNew;
    oldToNew.resize(preContigGraph.node.size(), (unsigned)-1);
    unsigned id = 0;
    for (unsigned i = 0; i < oldToNew.size(); ++i)
        if (!removeContigMap[i])
        {
            if (debug)
                std::cerr << "oldToNew[" << i << "] = " << id << "\n";
            oldToNew[i] = id++;
        }
    // Update FeatureMap.
    featureMap.applyContigIDMapping(oldToNew);
    if (true || debug)
    {
        std::cerr << "Feature Map after updating contig ID\n";
        featureMap.print(std::cerr);
    }

    // Now, expand the contigs using the readSet, featureMap, and consensus sequences from *fsPtr.
    rep_solv::ContigGraph finalGraph;
    std::vector<unsigned> expandContigMap;  // contig id new to old mapping, unused
    expandContigGraph(finalGraph, expandContigMap,  // output
                      contigGraph, consensusSeqs, readSet, repSolvOptions);  // input
    seqan::StringSet<seqan::Dna5String> expandedContigs;
    expandContigs(expandedContigs, *fsPtr, readSet, featureMap, oldToNew);
    if (debug)
    {
        std::cerr << "ORIGINAL MATE INFOS [" << __LINE__ << "]\n";
        mateInfos.print(std::cerr);
        std::cerr << "FINAL GRAPH before PE links [" << __LINE__ << "]\n";
        finalGraph.print(std::cerr);
    }
    if (true || debug)
    {
        std::cerr << "EXPANDED CONTIGS\n";
        for (unsigned i = 0; i < length(expandedContigs); ++i)
            std::cerr << ">" << i << "\n"
                      << expandedContigs[i] << "\n";
    }
    // Remove MateInfo objects that map to objects that were removed.  Then, renumber contigs.
    {
        // Erase links that touch a removed contig.
        eraseIf(mateInfos.records, [&oldToNew](rep_solv::MateInfo const & info) {
                return (oldToNew.at(info.leftID) == (unsigned)-1 || oldToNew.at(info.rightID) == (unsigned)-1);
            });
        // Renumber links.
        std::for_each(mateInfos.records.begin(), mateInfos.records.end(),
                      [&oldToNew](rep_solv::MateInfo & info) {
                          info.leftID = oldToNew.at(info.leftID);
                          info.rightID = oldToNew.at(info.rightID);
                      });
    }
    if (debug)
    {
        std::cerr << "CLEANSED MATE INFOS [" << __LINE__ << "]\n";
        mateInfos.print(std::cerr);
    }
    // Expand mate infos.
    rep_solv::MateInfos expandedMateInfos;
    expandMateInfos(expandedMateInfos, mateInfos, readSet, finalGraph, repSolvOptions);
    if (debug)
    {
        std::cerr << "EXPANDED MATE INFOS [" << __LINE__ << "]\n";
        expandedMateInfos.print(std::cerr);
    }
    addPairedEndLinks(finalGraph, expandedMateInfos);

    if (true || debug)
    {
        std::cerr << "FINAL GRAPH [" << __LINE__ << "]\n";
        finalGraph.print(std::cerr);
    }

    // Enumerate explaining haplotypes in case of existing s-t path.
    seqan::StringSet<seqan::Dna5String> scaffolds;
    if (stPathExists(finalGraph))
        enumeratePaths(scaffolds, finalGraph, expandedContigs, repSolvOptions);
    else
        buildScaffold(scaffolds, finalGraph, expandedContigs, repSolvOptions);

    std::cerr << "RESULTING SCAFFOLDS\n";
    for (unsigned i = 0; i < length(scaffolds); ++i)
        std::cout << ">" << i << "\n"
                  << scaffolds[i] << "\n";

    return 0;
}
