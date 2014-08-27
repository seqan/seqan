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

#include "read_cover.h"

#include <boost/heap/fibonacci_heap.hpp>

#include <utility>
#include <vector>

namespace rep_sep {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Class FeatureCoverageInformation
// ----------------------------------------------------------------------------

// Helper struct with information about the reads per feature.

class FeatureCoverageInformation
{
public:
    FeatureCoverageInformation(unsigned numFeatures, unsigned numReads, unsigned numSuperReads) :
            numFeatures(numFeatures), numReads(numReads), numSuperReads(numSuperReads)
    {
        readsPerFeature.resize(numSuperReads);
        for (auto & vec : readsPerFeature)
            vec.resize(numFeatures, 0);
        totalReadsPerFeature.resize(numFeatures, 0);
    }

    // Decrease counters for super reads for given feature.
    void decrease(unsigned featureID, std::vector<unsigned> const & superReadIDs)
    {
        bool LOGGING = false;
        SEQAN_CHECK(featureID < numFeatures, "Must be valid feature.");

        for (auto readID : superReadIDs)
        {
            SEQAN_ASSERT_GT(readsPerFeature[readID][featureID], 0u);
            if (LOGGING)
                std::cerr << "readsPerFeature[" << readID << "][" << featureID << "] ( == "
                          << readsPerFeature[readID][featureID] << ") -= 1\n";
            readsPerFeature[readID][featureID] -= 1;
        }
    }

    // Number of COLUMN features, we do not handle LINK features here.
    unsigned numFeatures;
    // Number of actual and super reads.
    unsigned numReads;
    unsigned numSuperReads;

    // readsPerFeature[i][j] is the number of actual reads for super read "i" on feature "j".
    std::vector<std::vector<unsigned>> readsPerFeature;
    // Total actual reads on each feature.
    std::vector<unsigned> totalReadsPerFeature;
};

// ----------------------------------------------------------------------------
// Class ReadInfo
// ----------------------------------------------------------------------------

// Information about the uncovered part of the actual read set.

class ReadInfo
{
public:
    ReadInfo(unsigned numFeatures, unsigned numReads) :
            numFeatures(numFeatures), numReads(numReads), superReadsForRead(numReads),
            featuresForRead(numReads), readsForFeature(numFeatures)
    {}

    static const unsigned INVALID = (unsigned)-1;

    // Mark read as being removed.
    void remove(unsigned readID)
    {
        SEQAN_CHECK(!featuresForRead.at(readID).empty(), "Must be valid readID=%u", readID);
        std::set<unsigned> toErase;
        for (auto featureID : featuresForRead.at(readID))
        {
            readsForFeature.at(featureID).erase(readID);
            toErase.insert(featureID);
        }
        for (auto featureID : toErase)
            featuresForRead[readID].erase(featureID);
    }

    // Returns whether was previously removed.
    bool isRemoved(unsigned readID) const
    {
        return featuresForRead[readID].empty();
    }

    // Number of COLUMN features, we do not handle LINK features here.
    unsigned numFeatures;
    // Number of reads.
    unsigned numReads;

    // The super reads for each read.
    std::vector<std::vector<unsigned>> superReadsForRead;
    // The feature for each read.
    std::vector<std::set<unsigned>> featuresForRead;
    // The reads covering each feature.
    std::vector<std::set<unsigned>> readsForFeature;
};

// ----------------------------------------------------------------------------
// Class ReadCoverAlgorithm
// ----------------------------------------------------------------------------

// Implementation of read cover algorithm.

class ReadCoverAlgorithm
{
public:
    ReadCoverAlgorithm(FeatureReadSet const & readSet,
                       FeatureReadSet const & atomicReadSet,
                       FeatureMap const & featureMap,
                       TFragmentStore const & fragStore,
                       ReadSeparatorOptions const & options) :
            featureCoverageInfo(featureMap.numFeatures(FeatureDescription::COLUMN),
                                length(fragStore.readSeqStore),
                                readSet.size()),
            readInfo(featureMap.numFeatures(FeatureDescription::COLUMN), length(fragStore.readSeqStore)),
            featureMap(featureMap),
            fragStore(fragStore),
            readSet(readSet),
            atomicReadSet(atomicReadSet),
            options(options)
    {
        fillReadInfo();
    }

    void run(FeatureReadSet & out);

private:

    // Helper struct, to be used in priority queue.

    struct WeightedReadID
    {
        static const unsigned INVALID = (unsigned)-1;

        WeightedReadID() = default;
        WeightedReadID(unsigned weight, unsigned readID) : weight(weight), readID(readID) {}

        bool operator<(WeightedReadID const & other) const
        { return (std::make_pair(weight, readID) < std::make_pair(other.weight, other.readID)); }

        // Weight of the super read.
        unsigned weight { 0 };
        // ID of read.
        unsigned readID { (unsigned)-1 };
    };

    typedef boost::heap::fibonacci_heap<WeightedReadID> TPriorityQueue;
    typedef TPriorityQueue::handle_type THandle;
    TPriorityQueue pQueue;
    std::vector<THandle> pqHandles;
    std::vector<WeightedReadID> weightedReadIDs;

    // Book keeping about feature coverage.
    FeatureCoverageInformation featureCoverageInfo;
    // Book keeping about the reads covering their features.
    ReadInfo readInfo;

    // FeatureMap used.
    FeatureMap const & featureMap;
    // Input fragment store.
    TFragmentStore const & fragStore;
    // Input read set
    FeatureReadSet const & readSet;
    // Atomic read set of original reads.
    FeatureReadSet const & atomicReadSet;
    // Configuration to use.
    ReadSeparatorOptions const & options;

    // Returns true in case the read has sufficient support.
    bool significant(Read const & read) const;

    // Fill the readInfo field.
    void fillReadInfo();

    // Mark read and subreads as selected.
    void select(Read const & read);

};

bool ReadCoverAlgorithm::significant(Read const & read) const
{
    for (unsigned featureID = 0; featureID < featureCoverageInfo.numFeatures; ++featureID)
    {
        auto readsOnFeature = featureCoverageInfo.totalReadsPerFeature[featureID];
        if (!readsOnFeature)
            continue;
        auto coveredReadsOnFeature = featureCoverageInfo.readsPerFeature[read.id][featureID];
        if (options.verbosity >= 3)
            std::cerr << "Read " << read << " covering " << coveredReadsOnFeature << " / " << readsOnFeature << "\n";
        if ((int)(100 * coveredReadsOnFeature / readsOnFeature) >= options.minCoverPercentage)
        {
            if (options.verbosity >= 3)
                std::cerr << "YES, significant: " << read << "\n";
            return true;
        }
    }
    if (options.verbosity >= 3)
        std::cerr << "NO, not significant: " << read << "\n";
    return false;
}

void ReadCoverAlgorithm::fillReadInfo()
{
    for (auto const & el : fragStore.alignedReadStore)
    {
        SEQAN_ASSERT_EQ(atomicReadSet.reads.at(el.readId).id, el.readId);
        SEQAN_ASSERT_EQ(atomicReadSet.reads.at(el.readId).subReads.size(), 1u);
        SEQAN_ASSERT_EQ(atomicReadSet.reads.at(el.readId).subReads[0], el.readId);
        auto const & read = atomicReadSet.reads[el.readId];
        for (auto const & feature : read.features)
        {
            unsigned featureID = feature.id;
            if (featureID >= readInfo.numFeatures)
                continue;  // is not a COLUMN feature
            if (options.verbosity >= 3)
                std::cerr << "readInfo.featuresForRead[" << el.readId << "].append(" << featureID << ")\n"
                          << "readInfo.readsForFeature[" << featureID << "].insert(" << el.readId << ")\n";
            readInfo.featuresForRead.at(el.readId).insert(featureID);
            readInfo.readsForFeature.at(featureID).insert(el.readId);
            featureCoverageInfo.totalReadsPerFeature.at(featureID) += 1;
        }
    }
    for (auto const & read : readSet.reads)
        for (auto readID : read.subReads)
        {
            if (options.verbosity >= 3)
                std::cerr << "readInfo.superReadsForRead[" << readID << "].push_back(" << read.id << ")\n";
            SEQAN_ASSERT_LT(read.id, featureCoverageInfo.readsPerFeature.size());
            SEQAN_ASSERT_LT(readID, readInfo.featuresForRead.size());
            for (auto featureID : readInfo.featuresForRead[readID])
            {
                if (featureID >= readInfo.numFeatures)
                    continue;  // is not a COLUMN feature
                SEQAN_ASSERT_LT(featureID, featureCoverageInfo.readsPerFeature[read.id].size());
                featureCoverageInfo.readsPerFeature[read.id][featureID] += 1;
            }
            readInfo.superReadsForRead[readID].push_back(read.id);
        }
}

void ReadCoverAlgorithm::select(Read const & read)
{
    if (options.verbosity >= 3)
        std::cerr << "SELECTING\t" << read.id << "\tWEIGHT\t" << weightedReadIDs[read.id].weight << "\n";

    for (auto readID : read.subReads)
        if (!readInfo.isRemoved(readID))
        {
            for (auto superReadID : readInfo.superReadsForRead[readID])
                if (weightedReadIDs[superReadID].readID != WeightedReadID::INVALID)
                {
                    weightedReadIDs[superReadID].weight -= 1;
                    pQueue.decrease(pqHandles[superReadID], weightedReadIDs[superReadID]);
                }
            for (auto featureID : readInfo.featuresForRead[readID])
            {
                if (featureID >= readInfo.numFeatures)
                    continue;  // is not a COLUMN feature
                featureCoverageInfo.decrease(featureID, readInfo.superReadsForRead[readID]);
            }
            readInfo.remove(readID);
        }
}

void ReadCoverAlgorithm::run(FeatureReadSet & out)
{
    if (options.verbosity >= 2)
        std::cerr << "Running read cover algorithm.\n";

    // Copy basic properties.
    std::set<unsigned> actualReadIDs;  // actual read ids selected
    out.numFeatures = readSet.numFeatures;

    // Priority queues and handles into its elements for decrease-key.
    pqHandles.resize(readSet.size());
    weightedReadIDs.resize(readSet.size());

    // Fill priority queue and handles.
    for (auto const & read : readSet.reads)
    {
        weightedReadIDs[read.id] = WeightedReadID(read.subReads.size(), read.id);
        if (options.verbosity >= 3)
            std::cerr << "weightedReadIDs[" << read.id << "] == (WeightedReadID("
                      << weightedReadIDs[read.id].weight << ", " << weightedReadIDs[read.id].readID << ")\n";
        pqHandles[read.id] = pQueue.push(weightedReadIDs[read.id]);
    }

    // Perform covering set selection.
    while (!pQueue.empty())
    {
        // Get shortcut to current best-looking read.
        Read const & read = readSet.reads[pQueue.top().readID];
        // Pop from queue and mark as removed so it is not considered further in select()/decrease().
        pQueue.pop();
        weightedReadIDs[read.id].readID = WeightedReadID::INVALID;

        // Check that the read is important enough.
        if (significant(read))
        {
            // Write to output.
            out.reads.push_back(read);
            out.reads.back().id = out.reads.size() - 1;  // assign new read id
            actualReadIDs.insert(read.subReads.begin(), read.subReads.end());

            // Mark as selected, will call the decreaseKey() and modify queue.
            select(read);
        }
    }

    // Write out actual reads.
    out.numActualReads = actualReadIDs.size();
    if (options.verbosity >= 2)
        std::cerr << "Done selecting covering reads\n";
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function selectReadCover()
// ----------------------------------------------------------------------------

void selectReadCover(FeatureReadSet & out,
                     FeatureReadSet const & in,
                     FeatureReadSet const & atomicReadSet,
                     FeatureMap const & featureMap,
                     TFragmentStore const & fragStore,
                     ReadSeparatorOptions const & options)
{
    ReadCoverAlgorithm algo(in, atomicReadSet, featureMap, fragStore, options);
    algo.run(out);
}

}  // namespace rep_sep
