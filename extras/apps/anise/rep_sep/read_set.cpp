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

#include "read_set.h"

#include <algorithm>
#include <ostream>
#include <map>
#include <numeric>

#include <seqan/basic.h>

#include "rep_sep/read.h"

namespace rep_sep {

// --------------------------------------------------------------------------
// Class FeatureReadSet
// --------------------------------------------------------------------------

void FeatureReadSet::print(std::ostream & out) const
{
    out << "Read Set (#reads= " << reads.size() << ")\n";
    for (auto const & read : reads)
        out << read << "\n";
}

// ----------------------------------------------------------------------------
// Function removeReads()
// ----------------------------------------------------------------------------

void removeReads(FeatureReadSet & out,
                 FeatureReadSet const & in,
                 std::vector<bool> removeContig)
{
    out.numActualReads = in.numActualReads;
    out.numFeatures = in.numFeatures;

    // Remove reads not written out.
    std::vector<unsigned> idMapping;
    for (auto const & read : in.reads)
        if (!removeContig[read.contigID])
            out.reads.push_back(read);

    // // Update contigID members.
    // std::sort(out.reads.begin(), out.reads.end(), ltContigID);
    // unsigned id = 0;
    // std::map<unsigned, unsigned> ids;
    // for (auto & read : out.reads)
    //     if (ids.count(read.contigID))
    //     {
    //         read.contigID = ids[read.contigID];
    //     }
    //     else
    //     {
    //         ids[read.contigID] = id;
    //         read.contigID = id++;
    //     }

    // Update IDs.
    std::sort(out.reads.begin(), out.reads.end(), ltID);
    for (unsigned i = 0; i < out.reads.size(); ++i)
        out.reads[i].id = i;
}

// ----------------------------------------------------------------------------
// Function fixFeatureCounts()
// ----------------------------------------------------------------------------

void fixFeatureCounts(FeatureReadSet & readSet, FeatureReadSet const & atomicReadSet)
{
#if SEQAN_ENABLE_DEBUG
    // Check that we can directly access reads by their ID.
    for (unsigned i = 0; i < atomicReadSet.reads.size(); ++i)
        SEQAN_ASSERT_EQ(atomicReadSet.reads[i].id, i);
#endif  // #if SEQAN_ENABLE_DEBUG
    
    for (auto & superRead : readSet.reads)
    {
        // Rebuild WeightedFeatureVector for superRead using atomicReadSet.

        // First, count number of occurences for each value of each feature.
        std::map<unsigned, std::map<int, int>> count;  // count[feature_id][value]
        for (auto readID : superRead.subReads)
            for (auto feature : atomicReadSet.reads[readID].features)
                ++count[feature.id][feature.value];

        // Then, rebuild features and swap out superRead.features for rebuilt one.
        WeightedFeatureVector vec;
        for (auto const & pair : count)
        {
            WeightedFeature ftr(pair.first, 0, 0);
            for (auto const valueCount : pair.second)
                if (valueCount.second > ftr.count)
                {
                    ftr.value = valueCount.first;
                    ftr.count = valueCount.second;
                }
            vec.insert(ftr);
        }
        swap(superRead.features, vec);
    }
}

// ----------------------------------------------------------------------------
// Function distributeUnplacedReads()
// ----------------------------------------------------------------------------

void distributeUnplacedReads(FeatureReadSet & readSet,
                             FeatureReadSet const & atomicReadSet)
{
    // Contigs occuring in atomicReadSet, such that we can create reads for contigs missing in readSet.
    std::set<unsigned> contigIDs;
    for (auto const & read : atomicReadSet.reads)
        contigIDs.insert(read.contigID);

    // Obtain set of unplaced reads and mapping from contigID to super reads in ReadSet that map to the contig.
    std::set<unsigned> unplaced;
    for (unsigned i = 0; i < atomicReadSet.size(); ++i)
        unplaced.insert(unplaced.end(), i);
    std::map<unsigned, std::set<unsigned> > contigSuperReads;
    unsigned superReadID = 0;
    for (auto const & read : readSet.reads)
    {
        SEQAN_ASSERT_EQ(read.id, superReadID);
        contigSuperReads[read.contigID].insert(read.id);
        for (auto readID : read.subReads)
            unplaced.erase(readID);
        ++superReadID;
    }

    // Now, distribute all yet unplaced atomic reads that have no feature to all super reads that cover this contig.
    for (auto readID : unplaced)
    {
        auto const & atomicRead = atomicReadSet.reads[readID];
        SEQAN_ASSERT_EQ(atomicRead.id, readID);
        if (atomicRead.features.empty())
            for (auto superReadID : contigSuperReads[atomicRead.contigID])
            {
                contigIDs.erase(atomicRead.contigID);
                readSet.reads[superReadID].mergeWithThis(atomicRead, 0, true);
            }
    }

    // Now, create one read in readSet for each yet unwritten contigID.
    for (auto contigID : contigIDs)
    {
        Read newRead;
        newRead.id = readSet.size();
        bool first = true;
        for (auto readID : unplaced)
        {
            auto const & atomicRead = atomicReadSet.reads[readID];
            if (!atomicRead.features.empty() || atomicRead.contigID != contigID)
                continue;  // ignore
            if (first)
            {
                newRead = atomicRead;
                first = false;
            }
            else
            {
                newRead.mergeWithThis(atomicRead, 0, true);
            }
        }
        if (!first)
            readSet.reads.push_back(newRead);
    }
}

}  // namespace rep_sep
