// ==========================================================================
//                                      ANISE
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

// TODO(holtgrew): In the first step, assembly is not really really requires, mapped reads should not become orphans. Maybe fill gaps with reference sequence?
// TODO(holtgrew): Handle singletons by purging the unitig/contig graph.

#include "assembler.h"

#include <memory>
#include <set>
#include <map>

#include <seqan/basic.h>

#include "asm/consensus_builder.h"
#include "asm/contig_graph.h"
#include "asm/frag_store.h"
#include "asm/store_utils.h"
#include "asm/unitigger.h"

#include "rep_sep/feature_map.h"
#include "rep_sep/read_set.h"
#include "rep_sep/read_correction.h"
#include "rep_sep/rep_sep.h"
#include "rep_sep/rep_sep_options.h"

#include "rep_solv/contig_graph.h"
#include "rep_solv/dead_branch.h"
#include "rep_solv/directed_tree.h"
#include "rep_solv/graph_expansion.h"
#include "rep_solv/options.h"
#include "rep_solv/path_enumeration.h"
#include "rep_solv/scaffolding.h"

#include "anise/bam_tag_names.h"
#include "anise/anise_options.h"
#include "anise/app_state.h"
#include "anise/site_data.h"
#include "anise/site_data_updater.h"
#include "anise/temporary_file_manager.h"

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Function buildReadBirthVector()
// ----------------------------------------------------------------------------

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

    // if (options.verbosity >= 2)
    // {
    //     std::cerr << "READ BIRTH STEP\n";
    //     for (unsigned recordID = 0; recordID < records.size(); ++recordID)
    //         std::cerr << recordID << "\t->\t" << result[recordID] << "\n";
    // }

    return result;
}

// ----------------------------------------------------------------------------
// Function computeMateInfos()
// ----------------------------------------------------------------------------

// Compute mate information from library Info, fragmentStore and records (latter for orientation).

rep_solv::MateInfos computeMateInfos(
        std::vector<rep_solv::LibraryInfo> const & libInfos,
        TFragmentStore const & fragStore,
        std::vector<seqan::BamAlignmentRecord> const & records,
        bool logging = false)
{
    sortAlignedReads(fragStore.alignedReadStore, seqan::SortReadId());
    SEQAN_ASSERT_EQ_MSG(records.size() % 2u, 0u, "Only paired-end reads at the moment.");

    rep_solv::MateInfos mateInfos;
    // Copy over library infos.
    mateInfos.libraries = libInfos;

    if (logging)
        std::cerr << "ADDING PE LINKS\n";
    for (unsigned idx = 0; idx < records.size(); idx += 2)
    {
        auto const & recordL = records[idx];
        auto const & recordR = records[idx + 1];
        SEQAN_CHECK(fragStore.alignedReadStore[idx].readId == idx, "Must be the same.");
        SEQAN_CHECK(fragStore.alignedReadStore[idx + 1].readId == idx + 1, "Must be the same.");
        auto const & elL = fragStore.alignedReadStore[idx];
        auto const & elR = fragStore.alignedReadStore[idx + 1];
        if (elL.contigId == elR.contigId)
            continue;  // skip links on the same contig

        SEQAN_CHECK(recordL.qName == recordR.qName, "Must be the same.");

        if (logging)
            std::cerr << "LINK FROM\n"
                      << "\trecordL.qName=" << recordL.qName << "\trecordL.flag=" << recordL.flag
                      << "\telL.contigId=" << elL.contigId << "\telL.beginPos=" << elL.beginPos
                      << "\telL.endPos=" << elL.endPos
                      << "\trecordL.seq=" << recordL.seq
                      << "\tfragStore.readSeqStore[" << elL.readId << "]=" << fragStore.readSeqStore[elL.readId]
                      << "\n"
                      << "\trecordR.qName=" << recordR.qName << "\trecordR.flag=" << recordR.flag
                      << "\telR.contigId=" << elR.contigId << "\telR.beginPos=" << elR.beginPos
                      << "\telR.endPos=" << elR.endPos
                      << "\trecordR.seq=" << recordR.seq
                      << "\tfragStore.readSeqStore[" << elR.readId << "]=" << fragStore.readSeqStore[elR.readId]
                      << "\n";

        typedef typename TFragmentStore::TReadPos         TReadPos;
        typedef typename TFragmentStore::TContigStore     TContigStore;
        typedef typename seqan::Value<TContigStore>::Type TContig;
        typedef typename TFragmentStore::TContigSeq       TContigSeq;
        typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
        TContigGaps contigGapsL(fragStore.contigStore[elL.contigId].seq, fragStore.contigStore[elL.contigId].gaps);
        TContigGaps contigGapsR(fragStore.contigStore[elR.contigId].seq, fragStore.contigStore[elR.contigId].gaps);

        SEQAN_ASSERT_NEQ(hasFlagRC(recordL), hasFlagRC(recordR));
        rep_solv::MateInfo info;
        auto leftBegin = toSourcePosition(contigGapsL, std::min(elL.beginPos, elL.endPos));
        auto leftEnd = toSourcePosition(contigGapsL, std::max(elL.beginPos, elL.endPos));
        auto rightBegin = toSourcePosition(contigGapsR, std::min(elR.beginPos, elR.endPos));
        auto rightEnd = toSourcePosition(contigGapsR, std::max(elR.beginPos, elR.endPos));
        if (!hasFlagRC(recordL))
            info = rep_solv::MateInfo(idx,
                                      idx + 1,
                                      elL.contigId,
                                      elR.contigId,
                                      leftBegin,
                                      rightEnd,
                                      leftEnd - leftBegin,
                                      rightEnd - rightBegin,
                                      /*libraryID=*/0);
        else
            info = rep_solv::MateInfo(idx + 1,
                                      idx,
                                      elR.contigId,
                                      elL.contigId,
                                      rightBegin,
                                      leftEnd,
                                      rightEnd - rightBegin,
                                      leftEnd - leftBegin,
                                      /*libraryID=*/0);
        if (logging)
            std::cerr << "\t" << info << "\n";
        mateInfos.insert(info);
    }

    return mateInfos;
}

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
    resize(v, std::remove_if(begin(v, seqan::Standard()),
                             end(v, seqan::Standard()), pred) - begin(v, seqan::Standard()));
}

// ----------------------------------------------------------------------------
// Function printHaplotypes()
// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------
// Function setUpStore()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Really required? Refine!

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

// ----------------------------------------------------------------------------
// UnitiggerResult
// ----------------------------------------------------------------------------

// Result from unitigging step, returned by UnitiggingDriver.

struct UnitiggerResult
{
    // Pointer to ContigGraph with the resulting contigs.
    std::unique_ptr<assembler::ContigGraph> cgPtr;
    // FragmentStore with the consensus of cgPtr.
    std::unique_ptr<TFragmentStore> fsPtr;
};

// ----------------------------------------------------------------------------
// UnitiggingDriver
// ----------------------------------------------------------------------------

// Driver for unitigging.

class UnitiggingDriver
{
public:
    UnitiggingDriver(seqan::StringSet<seqan::Dna5String> const & seqs,
                     AniseOptions const & options,
                     int siteID) :
            seqs(seqs), options(options), siteID(siteID)
    {}

    UnitiggerResult run();

private:

    // Compute overlaps of input sequences and write to ovlStore.
    void computeOverlaps(assembler::OverlapStore & ovlStore);
    // Perform the unitigging.
    void performUnitigging(assembler::OverlapStore const & ovlStore);
    // Compute the consensus sequence.
    void buildConsensus();

    // Input / State

    // Read sequences.
    seqan::StringSet<seqan::Dna5String> const & seqs;
    // Configuration.
    AniseOptions const & options;
    // Current site's ID.
    int siteID;

    // Output

    UnitiggerResult result;
};

void UnitiggingDriver::computeOverlaps(assembler::OverlapStore & ovlStore)
{
    if (options.verbosity >= 2)
        std::cerr << "Computing overlaps...\n";
    assembler::OverlapCandidateGeneratorOptions candOptions;
    candOptions.logging = (options.verbosity >= 3);
    assembler::OverlapCandidateGenerator candGen(candOptions);

    assembler::OverlapperOptions ovlOptions;
    ovlOptions.logging = (options.verbosity >= 3);
    assembler::Overlapper overlapper(ovlOptions);

    auto verify = [&](assembler::OverlapCandidate const & cand) {
        overlapper.computeOverlap(ovlStore, seqs[cand.seq0], seqs[cand.seq1], cand);
    };

    candGen.run(verify, seqs);

    ovlStore.refresh();
    if (options.verbosity >= 2)
        std::cerr << "  Computed " << ovlStore.size() << " overlaps\n";

}

void UnitiggingDriver::performUnitigging(assembler::OverlapStore const & ovlStore)
{
    if (options.verbosity >= 2)
        std::cerr << "Performing Unitigging...\n";
    bool debug = (options.verbosity >= 3);
    result.cgPtr = assembler::Unitigger(assembler::Unitigger::MIRA_LIKE,
                                        debug).run(seqs, ovlStore.overlaps());
    if (options.verbosity >= 2)
        result.cgPtr->print(std::cerr);
}

void UnitiggingDriver::buildConsensus()
{
    // Compute consensus into FragmentStore.
    if (options.verbosity >= 2)
        std::cerr << "Building Consensus...\n";
    assembler::Options asmOptions;
    asmOptions.verbosity = options.verbosity;
    assembler::ConsensusBuilder consensusBuilder(asmOptions);
    result.fsPtr = consensusBuilder.run(*result.cgPtr, seqs);

    // There must be a 1:1 correspondence between contig in ContigGraph and FragmentStore.
    SEQAN_CHECK(length(result.fsPtr->contigStore) == result.cgPtr->node.size(),
                "Must correspond! siteID=%d, contigStore contigs: %u, num nodes: %u", siteID,
                (unsigned)length(result.fsPtr->contigStore), (unsigned)result.cgPtr->node.size());
    // Reconstruct this correspondence and permute the contigStore.
    std::vector<unsigned> contigPermutation(length(result.fsPtr->contigStore), (unsigned)-1);
    for (auto const & el : result.fsPtr->alignedReadStore)
    {
        unsigned cgContigID = result.cgPtr->contig[result.cgPtr->readToContig[el.readId]].id;
        if (contigPermutation.at(el.contigId) == (unsigned)-1)
            contigPermutation[el.contigId] = cgContigID;
        else
            SEQAN_ASSERT_EQ(contigPermutation[el.contigId], cgContigID);
    }
    permuteStoreContigs(*result.fsPtr, contigPermutation);
}

UnitiggerResult UnitiggingDriver::run()
{
    // Compute the overlaps.
    assembler::OverlapStore ovlStore;
    computeOverlaps(ovlStore);

    // Perform Unitigging and write result to *result.cgPtr.
    performUnitigging(ovlStore);

    // Build consensus and write result to *result.fsPtr.
    buildConsensus();

    if (options.verbosity >= 2)
    {
        std::cerr << "ContigGraph from Assembly\n";
        result.cgPtr->print(std::cerr);

        std::cerr << "Resulting fragment store\n";
        printStore(std::cerr, *result.fsPtr);

        // Print store to output.
        for (unsigned i = 0; i < length(result.fsPtr->contigStore); ++i)
        {
            std::cout << ">contig_" << i << "\n";
            std::cout << result.fsPtr->contigStore[i].seq << "\n";
        }
    }

    // Fill read name store for FragmentStore, ConsensusBuilder only fills in contig and read seqs.
    setUpStore(*result.fsPtr);  // TODO(holtgrew): Really required?

    return std::move(result);
}

// ----------------------------------------------------------------------------
// Class RepeatSeparationResult
// ----------------------------------------------------------------------------

// Bundles the result of the repeat separation step.

struct RepeatSeparationResult
{
    // Read set as generated from the repeat separation.
    rep_sep::FeatureReadSet readSet;
    // Read set with covering reads (before independent read merging).  Can be used to select compatible contigs in the
    // contig graph.
    rep_sep::FeatureReadSet coveringReadSet;
    // Read set with atomic reads, i.e. each entry has exactly one atomic read and the id of the read is the same as the
    // id of the contained read.
    rep_sep::FeatureReadSet atomicReadSet;
    // Map of the features detected during repeat separation.
    rep_sep::FeatureMap featureMap;
};

// ----------------------------------------------------------------------------
// Class RepeatSeparationDriver
// ----------------------------------------------------------------------------

// Driver class for repeat separation module.

class RepeatSeparationDriver
{
public:
    RepeatSeparationDriver(TFragmentStore & fragStore,
                           rep_solv::MateInfos const & mateInfos,
                           std::vector<seqan::BamAlignmentRecord> const & records,
                           int stepNo,
                           AniseOptions const & options) :
            fragStore(fragStore), mateInfos(mateInfos), records(records), stepNo(stepNo), options(options)
    {
        (void)this->stepNo;  // only used during debugging
    }

    RepeatSeparationResult run();

private:

    // // Filters the read set by tyring to select the "right" haplotype for each contig.
    // //
    // // For each input contig with features, selects the haplotypes that contain at least 5% of this contig's reads with
    // // features from the previous round.
    // rep_sep::FeatureReadSet selectGoodHaplotypes(rep_sep::FeatureReadSet const & readSet,
    //                                              rep_sep::FeatureReadSet const & atomicReadSet);

    // Input

    // Pointer to FragmentStore from consensus building.  Non-const since we need to reorder reads for realignment.
    TFragmentStore & fragStore;
    // Information about the links.
    rep_solv::MateInfos const & mateInfos;
    // The BAM records.
    std::vector<seqan::BamAlignmentRecord> const & records;
    // The current step number.
    int stepNo;
    // Configuration.
    AniseOptions const & options;

    // Result.

    RepeatSeparationResult result;
};

// rep_sep::FeatureReadSet RepeatSeparationDriver::selectGoodHaplotypes(
//         rep_sep::FeatureReadSet const & readSet,
//         rep_sep::FeatureReadSet const & atomicReadSet)
// {
//     rep_sep::FeatureReadSet result = readSet;
//     if (readSet.empty())
//         return result;

//     // Compute highest contigID.
//     unsigned numContigs = std::max_element(readSet.begin(), readSet.end(),
//                                            [](rep_sep::Read const & lhs, rep_sep::Read const & rhs) {
//                                                return (lhs.contigID < rhs.contigID);
//                                            })->contigID + 1;

//     std::cerr << "SELECTING GOOD HAPLOTYPES (step #" << stepNo << ")\n";
//     std::cerr << "Input Read Set\n";
//     readSet.print(std::cerr);

//     // Select all contigs that have features on them and factorize reads by their contig.
//     std::vector<std::set<unsigned>> contigReads(numContigs);  // reads with feature by contig
//     for (auto const & read : readSet)
//         for (auto readID : read.subReads)
//         {
//             if (atomicReadSet.reads.at(readID).features.empty())
//                 continue;  // skip, has no feature
//             int readStepNo = 0;
//             unsigned idx = 0;
//             seqan::BamTagsDict tags(const_cast<seqan::CharString &>(records[readID].tags));
//             if (findTagKey(idx, tags, KEY_BIRTH_STEP))
//                 extractTagValue(readStepNo, tags, idx);

//             if (readStepNo + 1 < stepNo)
//                 contigReads.at(read.contigID).insert(readID);
//         }

//     std::cerr << "CONTIG READS (with features and from previous round)\n";
//     for (unsigned contigID = 0; contigID < contigReads.size(); ++contigID)
//     {
//         std::cerr << contigID << "\t";
//         for (auto readID : contigReads[contigID])
//             std::cerr << readID << " ";
//         std::cerr << "\n";
//     }

//     // Remove all elements in readSet that have less than MIN_CONTAINED reads with features.
//     double const MIN_CONTAINED = 0.05;  // smallest number of feature reads to contain on contig
//     unsigned const MIN_READS = 5u;  // smallest number of reads with features to require for deletion
//     eraseIf(result.reads, [&](rep_sep::Read const & read) {
//             unsigned numTotal = contigReads.at(read.contigID).size();
//             if (numTotal < MIN_READS)
//                 return false;

//             std::set<unsigned> onRead(read.subReads.begin(), read.subReads.end());
//             std::vector<unsigned> tmp;
//             std::set_intersection(onRead.begin(), onRead.end(),
//                                   contigReads.at(read.contigID).begin(), contigReads.at(read.contigID).end(),
//                                   std::back_inserter(tmp));
//             unsigned numOnRead = tmp.size();
//             return (numOnRead <= std::ceil(MIN_CONTAINED * numTotal));
//         });

//     std::cerr << "Resulting read set\n";
//     readSet.print(std::cerr);

//     return result;
// }

RepeatSeparationResult RepeatSeparationDriver::run()
{
    if (options.verbosity >= 2)
        std::cerr << "Performing repeat separation...\n";

    // Perform repeat separation (enumerates separating columns with Tammi's test and then cluster linking for
    // globalizing haplotypes).
    rep_sep::ReadSeparatorOptions rsOptions;
    auto rsOptions2 = rsOptions;
    rsOptions.verbosity = options.verbosity;
    auto readBirthStep = buildReadBirthVector(records);
    separateRepeats(result.readSet, result.coveringReadSet, result.atomicReadSet, result.featureMap, fragStore,
                    mateInfos, readBirthStep, rsOptions);

    if (options.verbosity >= 2)
    {
        std::cerr << "Atomic read set\n";
        result.atomicReadSet.print(std::cerr);
        std::cerr << "Read set after Separation\n";
        result.readSet.print(std::cerr);
        std::cerr << "Feature Map\n";
        result.featureMap.print(std::cerr);
        std::cerr << "Haplotypes\n";
        printHaplotypes(std::cerr, fragStore, result.readSet);
    }

    return std::move(result);
}

// ----------------------------------------------------------------------------
// Class RepeatResolutionResult
// ----------------------------------------------------------------------------

struct RepeatResolutionResult
{
    // Resulting super read set, corresponding to scaffoldSeqs.
    rep_sep::FeatureReadSet superReads;
    // Scaffold information.
    std::vector<rep_solv::AssemblyInfo> scaffoldInfos;
    // Resulting scaffold sequences.
    seqan::StringSet<seqan::Dna5String> scaffoldSeqs;
    // Masks for scaffold sequences, if false, corresponding char from scaffoldSeqs is written out as lower case.
    seqan::StringSet<seqan::String<bool>> masks;
};

// ----------------------------------------------------------------------------
// Class RepeatResolutionDriver
// ----------------------------------------------------------------------------

// Driver for the repeat resolution process.

class RepeatResolutionDriver
{
public:

    RepeatResolutionDriver(UnitiggerResult & unitiggerResult,
                           RepeatSeparationResult & repSepResult,
                           std::vector<seqan::BamAlignmentRecord> const & records,
                           std::vector<rep_solv::LibraryInfo> const & libInfos,
                           rep_solv::MateInfos const & mateInfos,
                           TFragmentStore & fragStore,
                           AniseOptions const & options,
                           unsigned stepNo) :
            unitiggerResult(unitiggerResult), repSepResult(repSepResult), records(records), libInfos(libInfos),
            mateInfos(mateInfos), fragStore(fragStore), options(options), stepNo(stepNo), repSolvOptions(buildRepSolvOptions())
    {}

    RepeatResolutionResult run();

private:

    // Fill repeat resolution options.
    rep_solv::Options buildRepSolvOptions() const
    {
        rep_solv::Options result;
        result.verbosity = options.verbosity;
        result.maxMateConflicts = -1;  // infinity
        return result;
    }

    // Use pair information in records to update the matePairId of the ReadStoreElement, add appropriate
    // MatePairStorElement objects.
    void updateStore();
    // Build contig graph for repeat resolution from assembled contig graph.  We call it pre contig graph since we will
    // remove some contigs in our heuristics.
    std::unique_ptr<rep_solv::ContigGraph> buildInitialContigGraph() const;
    // Augment ContigGraph with edges from s/to t.
    void addSTLinks(rep_solv::ContigGraph & graph) const;
    // Perform dead branch removal and directed tree reachability heuristic.  Also remove low-coverage contigs against
    // noise.  We do not directly modify the graph but work on a copy and return a bool vector with marks for to be
    // removed contigs.
    std::vector<bool> computePurgeMap(rep_solv::ContigGraph const & graph) const;
    // Compute old-to-new map for contig IDs.
    std::vector<unsigned> computeOldToNew(rep_solv::ContigGraph const & graph,
                                          std::vector<bool> const & toRemove) const;
    // Perform contig removal and old-to-new mapping.
    void performContigRemoval(rep_solv::ContigGraph & contigGraph,
                              rep_sep::FeatureReadSet & readSet,
                              seqan::StringSet<seqan::Dna5String> & consensusSeqs,
                              rep_solv::ContigGraph const & inGraph,
                              std::vector<bool> const & toRemove,
                              std::vector<unsigned> const & oldToNew);
    // Perform the contig and link expansion.
    void performExpansion(rep_solv::ContigGraph & finalGraph,
                          seqan::StringSet<seqan::Dna5String> & expandedContigs,
                          rep_solv::MateInfos &  mateInfos,  // we modify
                          rep_solv::ContigGraph const & contigGraph,
                          seqan::StringSet<seqan::Dna5String> const & consensusSeqs,
                          rep_sep::FeatureReadSet const & readSet,
                          std::vector<unsigned> const & oldToNew);
    // Apply masking to expanded contigs using the information about the aligned reads in unitiggerResult.
    void maskExpandedContigs(seqan::StringSet<rep_solv::TContigSeq> & maskedExpandedContigs,
                             seqan::StringSet<seqan::Dna5String> const & expandedContigs,
                             rep_sep::FeatureReadSet const & superReads,
                             std::vector<unsigned> const & oldToNew);

    // Perform scaffolding or path enumeration.
    void performScaffolding(seqan::StringSet<rep_solv::TContigSeq> & scaffolds,
                            std::vector<rep_solv::AssemblyInfo> & scaffoldInfos,
                            rep_sep::FeatureReadSet & superReads,
                            rep_solv::ContigGraph const & finalGraph,
                            seqan::StringSet<rep_solv::TContigSeq> const & expandedContigs) const;
    // Expand the contigs using the given FeatureReadSet, FeatureMap, and FragmentStore.  oldToNew maps fragStore
    // contig IDs to contigIDs after contig removal.
    void expandContigs(seqan::StringSet<seqan::Dna5String> & expandedContigs,
                       TFragmentStore const & fragStore,
                       rep_sep::FeatureReadSet const & readSet,
                       rep_sep::FeatureMap const & featureMap,
                       std::vector<unsigned> const & oldToNew) const;
    // Copy sequence in scaffold seqs that is masked in Dna5Q Strings to repeat resolution result with explicit masks.
    void writeMaskedScaffolds(RepeatResolutionResult & result,
                              seqan::StringSet<rep_solv::TContigSeq> const & scaffolds) const;

    // Input.

    // Result of the unitigging step, required for contig overlap information.
    UnitiggerResult & unitiggerResult;
    // Result of the repeat separation.
    RepeatSeparationResult & repSepResult;
    // BAM alignment records, used for the sequence and and names.
    std::vector<seqan::BamAlignmentRecord> const & records;
    // Library information.
    std::vector<rep_solv::LibraryInfo> libInfos;
    // Information about links; store copy, we modify.
    rep_solv::MateInfos mateInfos;
    // Fragment store with the alignment of the records to the assembled contigs, non-const since we need to be able to
    // sort the alignments.
    TFragmentStore & fragStore;
    // Configuration.
    AniseOptions const & options;
    // Current step number.
    unsigned stepNo;

    // State

    // Options for repeat resolution.
    rep_solv::Options repSolvOptions;
};

void RepeatResolutionDriver::updateStore()
{
    // Assumption: records is clustered by query name, i.e. if there is a pair then the two records occur after each
    // other.

    clear(fragStore.matePairStore);
    SEQAN_ASSERT_EQ(length(fragStore.readStore), records.size());
    for (unsigned idx = 0; idx < records.size();)
        if (!hasFlagMultiple(records[idx]))  // unpaired
        {
            fragStore.readStore[idx].matePairId = seqan::ReadStoreElement<>::INVALID_ID;
            idx += 1;  // advance
        }
        else  // is paired
        {
            // Check that assumption from above is correct.
            SEQAN_ASSERT_LT(idx + 1, records.size());
            SEQAN_ASSERT_EQ(records[idx].qName, records[idx + 1].qName);

            // Consistency checks.
            unsigned idF = idx, idL = idx + 1;
            SEQAN_ASSERT_NEQ(hasFlagFirst(records[idx]), hasFlagFirst(records[idx + 1]));
            SEQAN_ASSERT_NEQ(hasFlagLast(records[idx]), hasFlagLast(records[idx + 1]));
            SEQAN_ASSERT_EQ(hasFlagFirst(records[idx]), hasFlagLast(records[idx + 1]));
            SEQAN_ASSERT_EQ(hasFlagLast(records[idx]), hasFlagFirst(records[idx + 1]));
            if (!hasFlagFirst(records[idF]))
                std::swap(idF, idL);

            // Create new entry in matePairStore.
            unsigned matePairID = length(fragStore.matePairStore);
            resize(fragStore.matePairStore, matePairID + 1);
            back(fragStore.matePairStore).readId[0] = idF;
            back(fragStore.matePairStore).readId[1] = idL;

            fragStore.readStore[idF].matePairId = matePairID;
            fragStore.readStore[idL].matePairId = matePairID;

            idx += 2;  // advance
        }
}

RepeatResolutionResult RepeatResolutionDriver::run()
{
    RepeatResolutionResult result;

    if (options.verbosity >= 2)
        std::cerr << "Performing repeat resolution...\n";

    // Update mate pair info in fragment store.
    updateStore();

    // Compute initial contig graph from unitigger contig graph (only contains links from overlap and s-/t- links).
    std::unique_ptr<rep_solv::ContigGraph> initialContigGraph = buildInitialContigGraph();

    // Build string set of initial consensus seqs, so containment links in addPairedEndLinks can be resolved.
    seqan::StringSet<seqan::Dna5String> rawConsensusSeqs;
    for (unsigned i = 0; i < length(unitiggerResult.fsPtr->contigStore); ++i)
        appendValue(rawConsensusSeqs, unitiggerResult.fsPtr->contigStore[i].seq);
    // Compute which contigs to remove by dead branch removal and directed tree growing heuristic.
    addPairedEndLinks(*initialContigGraph, mateInfos, &rawConsensusSeqs);
    std::vector<bool> toRemove = computePurgeMap(*initialContigGraph);
    // Compute old-to-new map for contigs from toRemove.
    std::vector<unsigned> oldToNew = computeOldToNew(*initialContigGraph, toRemove);

    // Remove flagged contigs from *initialContigGraph, repSepResult.readSet, and get non-removed consensus sequences.
    // Also update repSepResult.featureMap.
    rep_solv::ContigGraph contigGraph;  // contig graph without removed contigs
    seqan::StringSet<seqan::Dna5String> consensusSeqs;  // non-removed consensus sequences
    if (options.verbosity >= 2)
    {
        std::cerr << "INPUT TO CONTIG REMOVAL\n";
        initialContigGraph->print(std::cerr);
    }
    performContigRemoval(contigGraph, result.superReads, consensusSeqs, *initialContigGraph, toRemove, oldToNew);

    // Return if there are no contigs left.
    if (contigGraph.node.empty())
        return result;

    if (options.verbosity >= 2)
    {
        std::cerr << "READ SET AFTER CONTIG REMOVAL\n";
        result.superReads.print(std::cerr);
        std::cerr << "ORIGINAL MATE INFOS [" << __LINE__ << "]\n";
        mateInfos.print(std::cerr);
    }
    // Erase mateInfos that touch a removed contig.
    eraseIf(mateInfos.records, [&oldToNew](rep_solv::MateInfo const & info) {
            return (oldToNew.at(info.leftID) == (unsigned)-1 || oldToNew.at(info.rightID) == (unsigned)-1);
        });
    // Renumber links.
    std::for_each(mateInfos.records.begin(), mateInfos.records.end(),
                  [&oldToNew](rep_solv::MateInfo & info) {
                      info.leftID = oldToNew.at(info.leftID);
                      info.rightID = oldToNew.at(info.rightID);
                  });
    if (options.verbosity >= 2)
    {
        std::cerr << "CLEANSED MATE INFOS [" << __LINE__ << "]\n";
        mateInfos.print(std::cerr);
    }

    // Now, expand the contigs using the readSet, (updated) repSepResult.featureMap, and consensus sequences from
    // *unitiggerResult.fsPtr.  During this, we can build the final graph that we will subject to path enumeration or
    // scaffolding.
    rep_solv::ContigGraph finalGraph;
    seqan::StringSet<seqan::Dna5String> expandedContigs;
    performExpansion(finalGraph, expandedContigs, mateInfos, contigGraph, consensusSeqs, result.superReads, oldToNew);

    // Apply masks to expanded contigs
    seqan::StringSet<rep_solv::TContigSeq> maskedExpandedContigs;
    maskExpandedContigs(maskedExpandedContigs, expandedContigs, result.superReads, oldToNew);

    // Enumerate explaining global haplotypes.
    seqan::StringSet<rep_solv::TContigSeq> scaffolds;
    performScaffolding(scaffolds, result.scaffoldInfos, result.superReads,
                       finalGraph, maskedExpandedContigs);

    // Copy out scaffold seqs and masks.
    writeMaskedScaffolds(result, scaffolds);

    // Remove mask in first step.
    if (stepNo <= 1u)
        clear(result.masks);

    return result;
}

void RepeatResolutionDriver::maskExpandedContigs(
        seqan::StringSet<rep_solv::TContigSeq> & maskedExpandedContigs,
        seqan::StringSet<seqan::Dna5String> const & expandedContigs,
        rep_sep::FeatureReadSet const & superReads,
        std::vector<unsigned> const & oldToNew)
{
    // Build mapping from contig IDs without removed to contigIDs before removal.
    std::vector<unsigned> newToOld(length(fragStore.contigStore), (unsigned)-1);
    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
        if (oldToNew.at(i) != (unsigned)-1)
            newToOld.at(oldToNew[i]) = i;

    sortAlignedReads(unitiggerResult.fsPtr->alignedReadStore, seqan::SortContigId());

    if (options.verbosity >= 2)
        std::cerr << "MASKING CONTIGS\n";

    SEQAN_ASSERT_EQ(length(expandedContigs), superReads.size());
    for (unsigned i = 0; i < length(expandedContigs); ++i)
    {
        unsigned contigID = superReads.reads[i].contigID;  // contigID in output
        unsigned storeContigID = newToOld[contigID];

        if (options.verbosity >= 2)
            std::cerr << "i=" << i << ", contigID=" << contigID << "\n"
                      << "seq=     " << unitiggerResult.fsPtr->contigStore[contigID].seq << "\n"
                      << "expanded=" << expandedContigs[i] << "\n";
        auto const & contig = expandedContigs[i];
        appendValue(maskedExpandedContigs, contig);
        auto & maskedContig = back(maskedExpandedContigs);


        auto itBegin = lowerBoundAlignedReads(unitiggerResult.fsPtr->alignedReadStore,
                                              storeContigID, seqan::SortContigId());
        auto itEnd = upperBoundAlignedReads(unitiggerResult.fsPtr->alignedReadStore,
                                            storeContigID, seqan::SortContigId());

        // Obtain gaps for coordinate translation.
        typedef typename TFragmentStore::TContigStore const TContigStore;
        typedef typename seqan::Value<TContigStore>::Type TContig;
        typedef typename TFragmentStore::TContigSeq const TContigSeq;
        typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<TContig::TGapAnchors>> TContigGaps;
        TContigGaps contigGaps(unitiggerResult.fsPtr->contigStore[storeContigID].seq,
                               unitiggerResult.fsPtr->contigStore[storeContigID].gaps);

        // Mask all as invalid, will mask new parts as valid again.
        for (auto itC = begin(maskedContig, seqan::Standard()); itC != end(maskedContig, seqan::Standard()); ++itC)
            assignQualityValue(*itC, rep_solv::MASK_QUAL);

        for (auto it = itBegin; it != itEnd; ++it)
        {
            // TODO(holtgrew): Which one is the better approach here?
            int readStepNo = 0;
            unsigned idx = 0;
            seqan::BamTagsDict tags(const_cast<seqan::CharString &>(records[it->readId].tags));
            if (findTagKey(idx, tags, KEY_BIRTH_STEP))
                extractTagValue(readStepNo, tags, idx);
            if ((unsigned)readStepNo + 1 < stepNo)
            // if (!hasFlagUnmapped(records[it->readId]))
            {
                if (options.verbosity >= 2)
                    std::cerr << "NOT UNMASKING WITH\t" << records[it->readId].qName << "\n";
                continue;  // mapped to early
            }

            auto beginPos = it->beginPos, endPos = it->endPos;
            if (beginPos > endPos)
                std::swap(beginPos, endPos);

            beginPos = toSourcePosition(contigGaps, beginPos);
            endPos = toSourcePosition(contigGaps, endPos);

            if ((unsigned)beginPos > length(maskedContig))
                beginPos = length(maskedContig);
            if ((unsigned)endPos > length(maskedContig))
                endPos = length(maskedContig);

            if (options.verbosity >= 2)
                std::cerr << "UNMASKING\t" << contigID << "\t" << beginPos << "\t" << endPos
                          << "\t" << records[it->readId].flag
                          << "\t" << records[it->readId].qName << "\t" << records[it->readId].seq
                          << "\t" << unitiggerResult.fsPtr->readSeqStore[it->readId] << "\n";

            // Mask aligned read by assigning the mask quality.
            for (auto itC = iter(maskedContig, beginPos); itC != iter(maskedContig, endPos); ++itC)
                assignQualityValue(*itC, 40);
        }
    }

    if (options.verbosity >= 2)
    {
        std::cerr << "EXPANDED+MASKED CONTIGS\n";
        for (auto const & eCtg : maskedExpandedContigs)
        {
            for (unsigned i = 0; i < length(eCtg); ++i)
                std::cerr << eCtg[i];
            std::cerr << "\n";
            for (unsigned i = 0; i < length(eCtg); ++i)
                std::cerr << (getQualityValue(eCtg[i]) != rep_solv::MASK_QUAL);
            std::cerr << "\n";
            std::cerr << "--\n";
        }
    }
}

void RepeatResolutionDriver::writeMaskedScaffolds(RepeatResolutionResult & result,
                                                  seqan::StringSet<rep_solv::TContigSeq> const & scaffolds) const
{
    clear(result.scaffoldSeqs);
    clear(result.masks);
    resize(result.masks, length(scaffolds));
    for (unsigned i = 0; i < length(scaffolds); ++i)
    {
        appendValue(result.scaffoldSeqs, scaffolds[i]);
        resize(result.masks[i], length(scaffolds[i]), true);
        unsigned j = 0;
        for (auto it = begin(scaffolds[i], seqan::Standard()); it != end(scaffolds[i], seqan::Standard()); ++it, ++j)
            result.masks[i][j] = (getQualityValue(*it) != rep_solv::MASK_QUAL);
    }
}

void RepeatResolutionDriver::performExpansion(rep_solv::ContigGraph & finalGraph,
                                              seqan::StringSet<seqan::Dna5String> & expandedContigs,
                                              rep_solv::MateInfos & mateInfos,  // we modify
                                              rep_solv::ContigGraph const & contigGraph,
                                              seqan::StringSet<seqan::Dna5String> const & consensusSeqs,
                                              rep_sep::FeatureReadSet const & readSet,
                                              std::vector<unsigned> const & oldToNew)
{
    std::vector<unsigned> expandContigMap;  // contig id new to old mapping, unused, can use readSet instead
    expandContigGraph(finalGraph,  expandContigMap,                          // output
                      contigGraph, consensusSeqs, readSet, repSolvOptions);  // input
    if (options.verbosity >= 2)
    {
        std::cerr << "FINAL GRAPH [" << __LINE__ << "]\n";
        finalGraph.print(std::cerr);
        std::cerr << "READ SET\n";
        readSet.print(std::cerr);
    }
    expandContigs(expandedContigs, *unitiggerResult.fsPtr, readSet, repSepResult.featureMap, oldToNew);
    for (auto const & node : finalGraph.node)
        finalGraph.contig[node].length = length(expandedContigs[finalGraph.contig[node].id]);
    if (options.verbosity >= 2)
    {
        std::cerr << "EXPANDED CONTIGS\n";
        for (unsigned i = 0; i < length(expandedContigs); ++i)
            std::cerr << ">" << i << "\n"
                      << expandedContigs[i] << "\n";
    }
    // Expand mate infos.
    rep_solv::MateInfos expandedMateInfos;
    expandMateInfos(expandedMateInfos, mateInfos, readSet, finalGraph, repSolvOptions);
    if (options.verbosity >= 2)
    {
        std::cerr << "EXPANDED MATE INFOS [" << __LINE__ << "]\n";
        expandedMateInfos.print(std::cerr);
    }
    // addPairedEndLinks(finalGraph, expandedMateInfos, &expandedContigs);

    if (options.verbosity >= 2)
    {
        std::cerr << "FINAL GRAPH [" << __LINE__ << "]\n";
        finalGraph.print(std::cerr);
    }
}

void RepeatResolutionDriver::performScaffolding(seqan::StringSet<rep_solv::TContigSeq> & scaffolds,
                                                std::vector<rep_solv::AssemblyInfo> & scaffoldInfos,
                                                rep_sep::FeatureReadSet & superReads,
                                                rep_solv::ContigGraph const & finalGraph,
                                                seqan::StringSet<rep_solv::TContigSeq> const & expandedContigs) const
{
    if (options.verbosity >= 2)
        std::cerr << "Running path enumeration.\n";
    std::vector<bool> overlapsWithFeature(records.size(), false);
    for (unsigned readID = 0; readID < overlapsWithFeature.size(); ++readID)
        overlapsWithFeature[readID] = !repSepResult.atomicReadSet.reads.at(readID).features.empty();

    std::vector<rep_solv::AnchorType> anchors;
    seqan::CharString str;
    for (auto const & record : records)
    {
        clear(str);
        seqan::BamTagsDict tagsDict(const_cast<seqan::CharString &>(record.tags));
        unsigned idx = 0;
        rep_solv::AnchorType anchorType = rep_solv::AnchorType::NONE;
        if (findTagKey(idx, tagsDict, KEY_ORIG_STRAND) && extractTagValue(str, tagsDict, idx))
            anchorType = ((str == "fwd") ? rep_solv::AnchorType::LEFT : rep_solv::AnchorType::RIGHT);
        anchors.push_back(anchorType);
    }

    auto readBirthStep = buildReadBirthVector(records);
    superReads.reads.clear();  // TODO(holtgrew): Do not write into this earlier!
    enumeratePaths(scaffolds, scaffoldInfos, superReads, finalGraph, expandedContigs,
                   repSepResult.readSet, readBirthStep, overlapsWithFeature, anchors, repSolvOptions);

    if (options.verbosity >= 2)
    {
        std::cerr << "RESULTING SCAFFOLDS\n";
        for (unsigned i = 0; i < length(scaffolds); ++i)
            std::cout << ">" << i << "\n"
                      << scaffolds[i] << "\n";
    }
}

void RepeatResolutionDriver::expandContigs(seqan::StringSet<seqan::Dna5String> & expandedContigs,
                                           TFragmentStore const & fragStore,
                                           rep_sep::FeatureReadSet const & readSet,
                                           rep_sep::FeatureMap const & featureMap,
                                           std::vector<unsigned> const & oldToNew) const
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

void RepeatResolutionDriver::performContigRemoval(rep_solv::ContigGraph & contigGraph,
                                                  rep_sep::FeatureReadSet & readSet,
                                                  seqan::StringSet<seqan::Dna5String> & consensusSeqs,
                                                  rep_solv::ContigGraph const & inGraph,
                                                  std::vector<bool> const & removeContigMap,
                                                  std::vector<unsigned> const & oldToNew)
{
    // Copy over removeContigMap into map for inGraph.
    lemon::SmartGraph::NodeMap<bool> doRemove(inGraph.graph, false);
    for (unsigned i = 0; i < removeContigMap.size(); ++i)
        if (removeContigMap[i])
        {
            if (options.verbosity >= 2)
                std::cerr << "doRemove[" << i << "] = true\n";
            doRemove[inGraph.node[i]] = true;
        }

    // Build a new contig graph without dead branches and branches into the "wrong" direction.
    if (options.verbosity >= 2)
    {
        std::cerr << "(pre) Read set before removing\n";
        readSet.print(std::cerr);
    }
    removeContigs(contigGraph, inGraph, doRemove);
    removeReads(readSet, repSepResult.readSet, removeContigMap);
    std::for_each(readSet.reads.begin(), readSet.reads.end(), [&oldToNew](rep_sep::Read & read) {
            read.contigID = oldToNew.at(read.contigID);
        });
    if (options.verbosity >= 2)
    {
        std::cerr << "(pre) Read set after removing\n";
        readSet.print(std::cerr);
    }

    // Fix feature counts in readSet.
    fixFeatureCounts(readSet, repSepResult.atomicReadSet);

    // Copy out consensus seqs, ignore removed ones.
    for (unsigned i = 0; i < length(unitiggerResult.fsPtr->contigStore); ++i)
        if (!removeContigMap[i])
            appendValue(consensusSeqs, unitiggerResult.fsPtr->contigStore[i].seq);

    if (options.verbosity >= 2)
    {
        std::cerr << "ContigGraph (after removing) [" << __LINE__ << "]\n";
        contigGraph.print(std::cerr);
        std::cerr << "Consensus seqs (after removing)\n";
        for (unsigned i = 0; i < length(consensusSeqs); ++i)
            std::cerr << ">" << i << "\n"
                      << consensusSeqs[i] << "\n";
    }
    if (options.verbosity >= 2)
    {
        std::cerr << "(pre) Read set after count fixing.\n";
        readSet.print(std::cerr);
    }

    // Fix / update feature map.
    if (options.verbosity >= 2)
    {
        std::cerr << "Feature Map before updating contig ID\n";
        repSepResult.featureMap.print(std::cerr);
    }
    repSepResult.featureMap.applyContigIDMapping(oldToNew);
    if (options.verbosity >= 2)
    {
        std::cerr << "Feature Map after updating contig ID\n";
        repSepResult.featureMap.print(std::cerr);
    }
}

std::vector<unsigned> RepeatResolutionDriver::computeOldToNew(rep_solv::ContigGraph const & graph,
                                                              std::vector<bool> const & toRemove) const
{
    std::vector<unsigned> oldToNew;
    oldToNew.resize(graph.node.size(), (unsigned)-1);
    unsigned id = 0;
    for (unsigned i = 0; i < oldToNew.size(); ++i)
        if (!toRemove[i])
        {
            if (options.verbosity >= 2)
                std::cerr << "oldToNew[" << i << "] = " << id << "\n";
            oldToNew[i] = id++;
        }
    return oldToNew;
}

std::vector<bool> RepeatResolutionDriver::computePurgeMap(rep_solv::ContigGraph const & graph) const
{
    std::vector<bool> removeContigMap(graph.node.size(), false);

    // if (options.verbosity >= 2)
    // {
    //     std::cerr << "Added PE mateInfos [" << __LINE__ << "]\n";
    //     for (auto const & info : mateInfos.records)
    //         std::cerr << "\t" << info << "\n";

    //     std::cerr << "PRE ContigGraph (with s-/t-links, PE links) [" << __LINE__ << "]\n";
    //     graph.print(std::cerr);
    // }

    auto readBirthStep = buildReadBirthVector(records);
    lemon::SmartGraph::NodeMap<bool> doRemove(graph.graph, false);
    unsigned numRemoved = removeDeadBranches(doRemove, graph, stepNo, readBirthStep, repSolvOptions);
    if (options.verbosity >= 2)
        std::cerr << "Dead branch removal: removed " << numRemoved << " vertices\n";
    // Remove contigs with too low support and too small length.
    std::vector<unsigned> support(graph.node.size(), 0);  // first collect support counts
    for (auto const & pair : graph.readToContigs)  // <(readID, <contigID0, contigID1, ..>)>
        for (auto const contigID : pair.second)
            support.at(contigID) += 1;
    for (lemon::SmartGraph::NodeIt u(graph.graph); u != lemon::INVALID; ++u)
        if (u != graph.s && u != graph.t)
            if (graph.contig[u].length < (unsigned)options.assemblyMinLength &&
                support.at(graph.contig[u].id) < (unsigned)options.assemblyMinSupport)
            {
                if (options.verbosity >= 2)
                    std::cerr << "support doRemove[" << graph.contig[u].id << "] = true\n";
                doRemove[u] = true;
            }
    // Execute directed tree growing heuristic.
    unsigned numUnreached = directedTreeGrowing(doRemove, graph, repSolvOptions);
    if (options.verbosity >= 2)
        std::cerr << "Directed tree heuristic: " << numUnreached << " unreachable nodes\n";
    // Write doRemove into removeContigMap.
    for (lemon::SmartGraph::NodeIt u(graph.graph); u != lemon::INVALID; ++u)
        if (doRemove[u])
        {
            removeContigMap.at(graph.contig[u].id) = true;
            if (options.verbosity >= 2)
                std::cerr << "removing contig " << graph.contig[u].id << "\n";
        }

    return removeContigMap;
}

std::unique_ptr<rep_solv::ContigGraph> RepeatResolutionDriver::buildInitialContigGraph() const
{
    std::unique_ptr<rep_solv::ContigGraph> result(new rep_solv::ContigGraph);
    auto & preContigGraph = *result;
    buildContigGraph(preContigGraph, *unitiggerResult.cgPtr, *unitiggerResult.fsPtr);
    // Augment the graph with s- and t-links.
    addSTLinks(preContigGraph);
    if (options.verbosity >= 2)
    {
        std::cerr << "PRE ContigGraph (with s-/t-links)\n";
        preContigGraph.print(std::cerr);
    }
    return result;
}

void RepeatResolutionDriver::addSTLinks(rep_solv::ContigGraph & graph) const
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

        if (options.verbosity >= 2)
            std::cerr << "s-t link\t" << record.qName << "\t" << cigar << "\t"
                      << "rc=" << hasFlagRC(record) << "\n";

        if (!hasFlagRC(record))
        {
            for (auto contigID : graph.readToContigs[readID])
            {
                if (s.count(contigID))
                    continue;  // already connected
                s.insert(contigID);
                graph.addEdge(graph.s, graph.node[contigID],
                              rep_solv::EdgeLabel(rep_solv::EdgeLabel::SOURCE, contigID));
            }
        }
        else
        {
            for (auto contigID : graph.readToContigs[readID])
            {
                if (t.count(contigID))
                    continue;  // already connected
                t.insert(contigID);
                graph.addEdge(graph.node[contigID], graph.t,
                              rep_solv::EdgeLabel(contigID, rep_solv::EdgeLabel::TARGET));
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Class StepAssembler
// ----------------------------------------------------------------------------

// Helper class for better structuring of the algorithm with state.

class StepAssembler
{
public:
    StepAssembler(SiteData & result,
                  SiteData const & in,
                  AniseOptions const & options,
                  AppState const & appState,
                  TemporaryFileManager & tmpMgr) :
            result(result), in(in), options(options), appState(appState), tmpMgr(tmpMgr)
    {}

    void run();

private:

    // Get library infos from application state.
    std::vector<rep_solv::LibraryInfo> libInfos() const
    {
        std::vector<rep_solv::LibraryInfo> result { rep_solv::LibraryInfo(appState.libraryInfos[0].median,
                                                                          appState.libraryInfos[0].stdDev) };
        return result;
    };

    // Input / Output

    // Output SiteData.
    SiteData & result;
    // Input SiteData.
    SiteData const & in;
    // Configuration of ANISE.
    AniseOptions const & options;
    // Current application state.
    AppState const & appState;
    // Temporary file manager to use for creating/opening temporary files.
    TemporaryFileManager & tmpMgr;
};

void StepAssembler::run()
{
    if (options.verbosity >= 2)
    {
        std::cerr << "SEQUENCES\n";
        for (unsigned i = 0; i < length(in.readSet.seqs); ++i)
            std::cerr << i << "\t" << in.readSet.seqs[i] << "\n";
    }

    // Perform unitigging.
    UnitiggingDriver unitigger(in.readSet.seqs, options, in.state.siteID);
    auto unitiggerResult = unitigger.run();

    // Write out uncorrected store.
    {
        std::fstream f;
        tmpMgr.open(f, std::ios::binary | std::ios::out, "raw_msa", ".txt", in.state.siteID, appState.stepNo);
        if (!f.good())
            throw std::runtime_error("Could not open raw_msa file for writing.");
        f << "# Read MSA before correction step " << appState.stepNo << " for site "
                  << in.state.siteID << "\n";
        printStore(f, *unitiggerResult.fsPtr);
    }

    // Compute mate infos, required in both repeat separation and resolution.
    rep_solv::MateInfos mateInfos = computeMateInfos(libInfos(), *unitiggerResult.fsPtr, in.readSet.bamRecords,
                                                     (options.verbosity >= 2));

    // Perform repeat separation.
    RepeatSeparationDriver repSeptor(*unitiggerResult.fsPtr, mateInfos, in.readSet.bamRecords, appState.stepNo,
                                     options);
    auto repSepResult = repSeptor.run();

    // Write out haplotypes.
    {
        std::fstream f;
        tmpMgr.open(f, std::ios::binary | std::ios::out, "haplos", ".txt", in.state.siteID, appState.stepNo);
        if (!f.good())
            throw std::runtime_error("Could not open haplos file for writing.");
        f << "# Read Sets\n";
        repSepResult.readSet.print(f);
        f << "# Feature Map\n";
        repSepResult.featureMap.print(f);
        f << "# Haplotypes for step " << appState.stepNo << " for site " << in.state.siteID << "\n";
        printHaplotypes(f, *unitiggerResult.fsPtr, repSepResult.readSet);
    }

    // Obtain a copy of the feature map before repeat resolution and contig shuffling.  We need it for the read
    // correction.
    decltype(repSepResult.featureMap) featureMapCopy;
    if (options.readCorrection)
        featureMapCopy = repSepResult.featureMap;

    // Perform repeat resolution.
    RepeatResolutionDriver repSolver(unitiggerResult, repSepResult, in.readSet.bamRecords, libInfos(),
                                     mateInfos, *unitiggerResult.fsPtr, options, appState.stepNo);
    auto repSolvResult = repSolver.run();
    // Check whether repeat resolution purged all contigs.
    if (empty(repSolvResult.scaffoldSeqs))
    {
        result.state = in.state;
        result.state.comment = "deactivating: empty scaffold seqs";
        result.state.active = false;
        return;
    }

    // TODO(holtgrew): Flag reads on removed contigs as removed.

    // Build result, starting with input results.
    result.state = in.state;
    result.readSet = in.readSet;

    if (options.readCorrection)
    {
        // Correct read sequences in result.readSet using the alignments from *unitiggerResult.fsPtr (for the consensus
        // sequence and alignment of reads against the same) and the features (for their location) in
        // repSepResult.featureMap.
        auto isReadCorrectable =[&](seqan::BamAlignmentRecord const & recordC) {
            // Get birth step of read for el.
            auto & record = const_cast<seqan::BamAlignmentRecord &>(recordC);
            seqan::BamTagsDict tagsDict(record.tags);
            int birthStep = 0;
            unsigned idx = 0;
            if (findTagKey(idx, tagsDict, KEY_BIRTH_STEP))
                extractTagValue(birthStep, tagsDict, idx);
            // Correct read if read old enough.
            return (birthStep < (int)appState.stepNo - 1);
        };
        performReadCorrection(result.readSet.bamRecords, *unitiggerResult.fsPtr,
                              featureMapCopy, isReadCorrectable);
    }

    // Write out corrected MSA.
    {
        std::fstream f;
        tmpMgr.open(f, std::ios::binary | std::ios::out, "corr_msa", ".txt", in.state.siteID, appState.stepNo);
        if (!f.good())
            throw std::runtime_error("Could not open corr_msa file for writing.");
        f << "# Read MSA after correction in step " << appState.stepNo << " for site "
                  << in.state.siteID << "\n";
        printStore(f, *unitiggerResult.fsPtr);
    }

    // Write out spanning state.
    for (auto info : repSolvResult.scaffoldInfos)
        appendValue(result.scaffold.scaffoldInfos,
                    Scaffold::ScaffoldInfo(info.anchoredLeft, info.anchoredRight, info.spanning));

    // Need to re-sort store for site data updating.
    sortAlignedReads(unitiggerResult.fsPtr->alignedReadStore, seqan::SortReadId());
    // Rebuild site data, scaffold seqs come directly from repeat separation.
    result.scaffold.siteName = in.scaffold.siteName;
    result.scaffold.seqs = repSolvResult.scaffoldSeqs;
    result.scaffold.masks = repSolvResult.masks;
    updateSiteData(result, *unitiggerResult.fsPtr, repSolvResult.superReads, options);

    if (empty(result.scaffold.seqs))
    {
        result.state.comment = "deactivating: empty scaffold seqs";
        result.state.active = false;
        return;
    }
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function performAssemblyStep()
// ----------------------------------------------------------------------------

void performAssemblyStep(SiteData & result, SiteData const & in, AniseOptions const & options,
                         AppState const & appState, TemporaryFileManager & tmpMgr)
{
    AniseOptions options2(options);
    // options2.verbosity = 2;
    StepAssembler helper(result, in, options2, appState, tmpMgr);
    helper.run();
}
