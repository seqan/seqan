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

#include "rep_sep.h"

#include <map>

#include "rep_sep/clique_enum.h"
#include "rep_sep/clique_set.h"
#include "rep_sep/cluster_linking.h"
#include "rep_sep/feature_map.h"
#include "rep_sep/local_variation_store.h"
#include "rep_sep/merge_read_set.h"
#include "rep_sep/pair_based_sep.h"
#include "rep_sep/read.h"
#include "rep_sep/read_set.h"
#include "rep_sep/read_cover.h"
#include "rep_sep/separating_columns.h"
#include "rep_sep/string_packed_pop_count.h"

namespace rep_sep {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Function splitReadSetByContig()
// ----------------------------------------------------------------------------

FeatureReadSet splitReadSetByContig(FeatureReadSet const & in,
                                    FeatureReadSet const & atomicReadSet)
{
    FeatureReadSet result;

    for (auto const & inRead : in)
    {
        // Factorize inRead by contigID, subClusters stores the read ids by cluster.
        std::map<unsigned, std::set<unsigned>> subClusters;
        for (auto readID : inRead.subReads)
            subClusters[atomicReadSet.reads.at(readID).contigID].insert(readID);

        // Create a Read in the output for each entry in subClusters.
        for (auto pair : subClusters)
        {
            Read outRead;
            // unsigned contigID = pair.first;
            auto const & subCluster = pair.second;
            bool first = true;
            for (auto readID : subCluster)
                if (first)
                {
                    outRead = atomicReadSet.reads.at(readID);
                    SEQAN_ASSERT_EQ(outRead.id, readID);
                    outRead.id = result.reads.size();
                    first = false;
                }
                else
                {
                    outRead.mergeWithThis(atomicReadSet.reads.at(readID), -1, true);
                }
            result.reads.emplace_back(std::move(outRead));
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
// Function eraseIf()
// ----------------------------------------------------------------------------

template <typename TValue, typename TPred>
void eraseIf(std::vector<TValue> & v, TPred pred)
{
    v.erase(std::remove_if(v.begin(), v.end(), pred), v.end());
}

// ----------------------------------------------------------------------------
// Function updateCounts()
// ----------------------------------------------------------------------------

// Update read set counts from atomic read set.

void updateCounts(FeatureReadSet & readSet, FeatureReadSet const & atomicReadSet)
{
    // Rebuild read set from reads in atomicReadSet, this will give us valid counts.
    FeatureReadSet result;

    for (auto const & inRead : readSet.reads)
    {
        // std::cerr << "RECREATING " << inRead << "\n";

        Read read(inRead.id, inRead.contigID, inRead.beginPos, inRead.endPos);
        for (auto readID : inRead.subReads)
        {
            SEQAN_CHECK(atomicReadSet.reads[readID].id == readID, "Must be equal readID=%u", readID);
            if (readID != atomicReadSet.reads[readID].subReads.front())
                continue;  // only once for pairs
            SEQAN_ASSERT_LEQ(atomicReadSet.reads[readID].subReads.size(), 2u);
            SEQAN_ASSERT_EQ(readID, atomicReadSet.reads[readID].subReads.front());
            // std::cerr << "  MERGING " << atomicReadSet.reads[readID] << " INTO " << read << "\n";
            int const IGNORE_COUNT = -1;
            read.mergeWithThis(atomicReadSet.reads[readID], IGNORE_COUNT);
        }

        result.reads.push_back(read);
    }

    swap(readSet, result);
}

// ----------------------------------------------------------------------------
// Function buildFeatureMap()
// ----------------------------------------------------------------------------

// Build feature map from local variation store.

FeatureMap buildFeatureMap(LocalVariationStore const & varStore)
{
    FeatureMap featureMap;

    for (unsigned columnID = 0; columnID < length(varStore.columns); ++columnID)
        featureMap.insert(FeatureDescription(/*id=*/columnID,
                                             /*kind=*/FeatureDescription::COLUMN,
                                             /*contigID=*/varStore.positions[columnID].first,
                                             /*pos=*/varStore.positions[columnID].second,
                                             /*coverage=*/length(varStore.columns[columnID])));

    featureMap.refresh();
    return featureMap;
}

// ----------------------------------------------------------------------------
// Class ClusterLinkingDriver
// ----------------------------------------------------------------------------

// Helper for running cluster linking.
//
// The result is a FeatureReadSet, each Read representing a cluster.  The read can span multiple contigs.  We have to
// distribute the reads to their contig downstream.

class ClusterLinkingDriver
{
    struct ReadFeature
    {
        unsigned readID { (unsigned)-1 };
        int value { -1 };

        bool operator<(ReadFeature const & other) const
        { return std::make_pair(readID, value) < std::make_pair(other.readID, other.value); }

        ReadFeature() = default;
        ReadFeature(unsigned readID, int value) : readID(readID), value(value) {}
    };

public:
    ClusterLinkingDriver(FeatureReadSet & readSet,
                         FeatureReadSet const & atomicReadSet,
                         FeatureMap const & featureMap,
                         TFragmentStore /*const*/ & inStore,
                         ReadSeparatorOptions const & options) :
            readSet(readSet), atomicReadSet(atomicReadSet), featureMap(featureMap), /*inStore(inStore),*/
            options(options),
            contigID(atomicReadSet.numActualReads, (unsigned)-1),
            globalizer(/*minOverlap=*/3, (options.verbosity >= 2))
    {
        globalizer.init(inStore);  // for mate ids
    }

    void run();

private:

    // Build map from feature id to read.
    void buildMappings();
    // Perform cluster linking step.
    void performLinking();
    // Build resulting FeatureReadSet.
    void buildResult();

    // Output / Input

    FeatureReadSet & readSet;
    FeatureReadSet const & atomicReadSet;
    FeatureMap const & featureMap;
    // Input fragment store.
    // TFragmentStore /*const*/ & inStore;
    // Configuration for repeat separation.
    ReadSeparatorOptions const & options;

    // State

    std::map<unsigned, std::set<ReadFeature>> readsForFeature;
    std::vector<unsigned> contigID;  // read to contig
    ClusterGlobalizer globalizer;
};

void ClusterLinkingDriver::buildMappings()
{
    // Build map from feature id to read.
    unsigned readSetEntryID = 0;
    for (auto const & read : atomicReadSet)
    {
        for (unsigned readID : read.subReads)
            contigID.at(readID) = read.contigID;

        for (auto feature : read.features)
        {
            auto it = featureMap.find(feature.id);
            if (it == featureMap.end() || it->kind != FeatureDescription::COLUMN)
                continue;
            for (unsigned readID : read.subReads)
                readsForFeature[feature.id].insert(ReadFeature(readID, feature.value));
        }

        ++readSetEntryID;
    }
}

void ClusterLinkingDriver::performLinking()
{
    // Perform cluster linking steps.
    //
    // The mapping values is for mapping the actual values to 0..(k-1).
    for (auto const & featureReads : readsForFeature)
    {
        auto const & featureValues = featureReads.second;
        std::map<unsigned, unsigned> mapping;  // read -> value (first actual, later 0..(k-1)).
        for (auto const & pair : featureValues)
            mapping[pair.readID] = pair.value;
        std::map<int, int> values;  // actual -> 0..(k-1)
        for (auto const & pair : mapping)
            if (!values.count(pair.second))
            {
                unsigned k = values.size();
                values[pair.second] = k;
            }
        if (options.verbosity >= 3)
        {
            std::cerr << "VALUES\n";
            for (auto const & pair : values)
                std::cerr << pair.first << " -> " << pair.second << "\n";
        }
        for (auto & pair : mapping)
            pair.second = values[pair.second];
        globalizer.run(values.size(), mapping);
    }
}

void ClusterLinkingDriver::buildResult()
{
    // Finalize and build resulting FeatureReadSet.
    readSet.reads.clear();
    globalizer.stop();

    for (auto const & cluster : globalizer.partition())
    {
        if (cluster.empty())
            continue;
        Read read;
        bool first = true;
        for (auto readID : cluster)
            if (first)
            {
                read = atomicReadSet.reads.at(readID);
                SEQAN_ASSERT_EQ(read.id, readID);
                read.id = readSet.reads.size();
                first = false;
            }
            else
            {
                read.mergeWithThis(atomicReadSet.reads.at(readID), -1, true, true);
            }
        readSet.reads.emplace_back(std::move(read));
    }
}

void ClusterLinkingDriver::run()
{
    buildMappings();
    performLinking();
    buildResult();

    if (options.verbosity >= 2)
        std::cerr << "Done Cluster Linking.\nComputing read cover.\n";
}

// ----------------------------------------------------------------------------
// Class RepeatSeparator
// ----------------------------------------------------------------------------

// Helper for separateRepeats().

class RepeatSeparator
{
public:
    RepeatSeparator(TFragmentStore /*const*/ & inStore,
                    rep_solv::MateInfos const & mateInfos,
                    std::vector<int> const & readBirthStep,
                    ReadSeparatorOptions const & options) :
            inStore(inStore), inMateInfos(mateInfos), readBirthStep(readBirthStep), options(options)
    {}

    void run(FeatureReadSet & readSet,
             FeatureReadSet & outCoveringSet,
             FeatureReadSet & atomicReadSet,
             FeatureMap & featureMap);

private:

    // // Run clique enumeration
    // void runCliqueEnumeration(FeatureReadSet & readSet, FeatureMap const & featureMap);
    // Run Kuchenbecker's cluster linking.
    //
    // Read set must be "atomic", each read has one sub read with the id of the read itself.
    void runClusterLinking(FeatureReadSet & readSet,
                           FeatureReadSet const & atomicReadSet,
                           FeatureMap const & featureMap);

    // Build local variation store.
    void buildLocalVariationStore(LocalVariationStore & varStore)
    {
        sortAlignedReads(inStore.alignedReadStore, seqan::SortBeginPos());
        sortAlignedReads(inStore.alignedReadStore, seqan::SortContigId());
        fill(varStore, inStore, options);
        // Remove columns with too high coverage.
        seqan::String<unsigned> columnIds;
        for (unsigned columnID = 0; columnID < length(varStore.columns); ++columnID)
            if (length(varStore.columns[columnID]) < options.maxColumnSize)
                appendValue(columnIds, columnID);
        filter(varStore, columnIds);
        // Enumerate separating columns.
        clear(columnIds);
        auto options2 = options;
        // options2.verbosity = 3;
        enumerateSeparatingColumns(columnIds, varStore, options2);
        filter(varStore, columnIds);
    }

    // Input fragment store.  We copy since we have to add pseudo-read entries for the previous round.
    TFragmentStore /*const*/ inStore;
    // Information about mate links between the contigs.
    rep_solv::MateInfos const & inMateInfos;
    // Birth step of each read, if 0 then was part of the original mapping.
    std::vector<int> const & readBirthStep;
    // Configuration for repeat separation.
    ReadSeparatorOptions const & options;
};

void RepeatSeparator::runClusterLinking(FeatureReadSet & readSet,
                                        FeatureReadSet const & atomicReadSet,
                                        FeatureMap const & featureMap)
{
    ClusterLinkingDriver helper(readSet, atomicReadSet, featureMap, inStore, options);
    helper.run();
}

// void RepeatSeparator::runCliqueEnumeration(FeatureReadSet & readSet, FeatureMap const & /*featureMap*/)
// {
//     // Iterative cluster enumeration and read merging.
//     for (unsigned round = 0; true; ++round)
//     {
//         if (options.verbosity >= 3)
//             std::cerr << "Starting enumeration round (#" << round << "), readSet.size() == " << readSet.size() << ".\n";
//         // Perform clique enumeration / clustering.
//         CliqueSet cliqueSet(readSet.size());
//         if (!performCliqueEnumeration(cliqueSet, readSet, options))
//             break;

//         if (options.verbosity >= 3)
//         {
//             std::cerr << "Resulting clique set (#" << round << ")\n";
//             for (auto const & clique : cliqueSet.cliques())
//                 std::cerr << clique << "\n";
//         }

//         // Merge reads to super reads to yield input for next iteration.
//         FeatureReadSet mergedSet;
//         auto ltId = [](Read const & lhs, Read const & rhs) { return (lhs.id < rhs.id); };
//         std::sort(readSet.reads.begin(), readSet.reads.end(), ltId);
//         mergeToSuperReads(mergedSet, cliqueSet, readSet, (options.verbosity >= 3));
//         swap(readSet, mergedSet);

//         if (options.verbosity >= 3)
//         {
//             std::cerr << "Round #" << round << " -- read set\n";
//             for (auto const & read : readSet.reads)
//                 std::cerr << read << "\n";
//         }
//     }

//     if (options.verbosity >= 2)
//         std::cerr << "Done enumerating.\nComputing read cover.\n";
// }

void RepeatSeparator::run(FeatureReadSet & outReadSet,
                          FeatureReadSet & outCoveringSet,
                          FeatureReadSet & atomicReadSet,
                          FeatureMap & outFeatureMap)
{
    // Build LocalVariationStore, call sepearating columns, and filter down to these columns.
    LocalVariationStore varStore;
    buildLocalVariationStore(varStore);

    // Build a FeatureMap from LocalVariationStore;
    outFeatureMap = buildFeatureMap(varStore);
    if (options.verbosity >= 3)
        outFeatureMap.print(std::cerr);

    // Build the FeatureReadSets.  We build atomicReadSet which basically is an annotation of reads with their covered
    // features and readSet in which pairs mapping on the same contig are joint into the same Read.  Reads without any
    // features are not included in readSet for running time reasons.
    FeatureReadSet readSet;
    buildAtomicReadSet(atomicReadSet, inStore, varStore);
    buildFeatureReadSet(readSet, inStore, varStore);

    // Augment the FeatureMap and FeatureReadSets with pair-based information.
    performPairBasedSeparation(outFeatureMap, readSet, atomicReadSet, inStore, inMateInfos);

    if (options.verbosity >= 2)
    {
        std::cerr << "Augmented Feature Map\n";
        outFeatureMap.print(std::cerr);
        std::cerr << "Atomic read set\n";
        atomicReadSet.print(std::cerr);
        std::cerr << "Initial read set\n";
        readSet.print(std::cerr);
    }

    // Store copy of initial read set.
    FeatureReadSet initialReadSet(readSet);
    // Run the cluster linking algorithm.
    runClusterLinking(readSet, atomicReadSet, outFeatureMap);

    if (options.verbosity >= 2)
    {
        std::cerr << "Selecting read cover from the following reads\n";
        for (auto const & read : readSet.reads)
            std::cerr << read << "\n";
    }

    // Select an approximately minimal covering subset of reads.
    selectReadCover(outCoveringSet, readSet, atomicReadSet, outFeatureMap, inStore, options);

    if (options.verbosity >= 2)
    {
        std::cerr << "Selected read cover is\n";
        for (auto const & read : outCoveringSet.reads)
            std::cerr << read << "\n";
    }

    FeatureReadSet coveringSetByContig = splitReadSetByContig(outCoveringSet, atomicReadSet);

    if (options.verbosity >= 2)
    {
        std::cerr << "Selected read cover, factorized by contigID, is\n";
        for (auto const & read : coveringSetByContig.reads)
            std::cerr << read << "\n";
    }

    // Update counts in covering read set.
    updateCounts(coveringSetByContig, atomicReadSet);

    if (options.verbosity >= 2)
    {
        std::cerr << "Covering read set (this round), before aux removal\n";
        coveringSetByContig.print(std::cerr);
    }

    // Merge non-conflicting reads on contig.
    outReadSet.numActualReads = initialReadSet.numActualReads;
    outReadSet.numFeatures = initialReadSet.numFeatures;
    mergeNonConflictingReadsOnContig(outReadSet, coveringSetByContig, initialReadSet, readBirthStep, options);
    std::stable_sort(outReadSet.reads.begin(), outReadSet.reads.end(), ltContigID);
    // TODO(holtgrew): Currently makes count incorrect again. updateCounts() again?

    if (options.verbosity >= 2)
    {
        std::cerr << "After merging non-conflicting reads\n";
        outReadSet.print(std::cerr);
    }

    // Distribute yet unplaced reads from atomicReads to outReadSet.
    distributeUnplacedReads(outReadSet, atomicReadSet);

    // Renumber outReadSet.
    for (unsigned i = 0; i < outReadSet.size(); ++i)
        outReadSet.reads[i].id = i;

    if (options.verbosity >= 2)
    {
        std::cerr << "After distributing unplaced reads\n";
        outReadSet.print(std::cerr);
    }

    if (options.verbosity >= 2)
        std::cerr << "Done computing read cover.\n";
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function separateRepeats()
// ----------------------------------------------------------------------------

void separateRepeats(FeatureReadSet & readSet,
                     FeatureReadSet & outCoveringSet,
                     FeatureReadSet & atomicReadSet,
                     FeatureMap & featureMap,
                     TFragmentStore /*const*/ & inStore,
                     rep_solv::MateInfos const & mateInfos,
                     std::vector<int> const & readBirthStep,
                     ReadSeparatorOptions const & options)
{
    RepeatSeparator sep(inStore, mateInfos, readBirthStep, options);
    sep.run(readSet, outCoveringSet, atomicReadSet, featureMap);
}

}  // namespace rep_sep
