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
// Haplotype reconstruction / repeat separation globalization using clique
// enumeration.
// ==========================================================================

// TODO(holtgrew): Use other significance in read cover, cover certain percentage of each feature.

#include "clique_enum.h"

#include <ostream>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include <seqan/random.h>
#include <seqan/realign.h>
#include <seqan/misc/misc_union_find.h>

#include <boost/heap/fibonacci_heap.hpp>

#include "asm/store_utils.h"

#include "rep_sep/feature_map.h"
#include "rep_sep/rep_sep_options.h"
#include "rep_sep/local_variation_store.h"
#include "rep_sep/read.h"
#include "rep_sep/read_set.h"
#include "rep_sep/clique_set.h"

namespace rep_sep {

namespace {  // anonymous namespace

// --------------------------------------------------------------------------
// FragmentStoreSplitter
// --------------------------------------------------------------------------

// Helper struct for "split()" function.
//
// Constructs a new FragmentStore and vector of BamAlignmentRecord objects from the input store and records based on the
// haplotypes in the readSet.
//
// We consider each super read R in turn.  Each read pair where both reads are part of R are copied to a new contig.  We
// "clone" the reads with new names in case that we already copied out these reads before and the original names are
// copied to a tag.  The case where the reads from one pair are placed into different super reads (and thus also
// contigs) is handled separately.
//
// At the moment, we create a copy of this pair for each combination.  This might lead to larger stack but should offer
// few problems and can be fixed by subsampling later.

class FragmentStoreSplitter
{
public:
    FragmentStoreSplitter(
            TFragmentStore & outStore,
            FeatureReadSet const & readSet,
            TFragmentStore const & inStore,
            bool logging = false) :
            rng(length(inStore.alignedReadStore)), outStore(outStore), readSet(readSet),
            inStore(inStore), logging(logging)
    {
        buildAverageCoverages();
    }

    void run();

private:

    // Check that inStore.alignedReadStore is sorted by read ID.
    void checkAlignedReadStore() const;

    // Copy out pairs where reads align to more than one contig.
    void writeOutLinks(seqan::Rng<> & rng);

    // Copy out alignments for inReadID0 and inReadID1 to output to contig outContigID0, outContigID1.
    void copyOut(unsigned inReadID0, unsigned inReadID1, unsigned outContigID0, unsigned outContigID1);

    // State

    // The average coverage for each super read in readSet.
    std::vector<double> avgCoverages;
    // Assigned super reads / contigs for each read.
    std::vector<std::vector<unsigned>> assignments;
    // Return vector of average coveragges for each input super read / output contig.  Multi-assigned reads contribute
    // only a fraction.  Also builds assignments.
    void buildAverageCoverages()
    {
        assignments.resize(length(inStore.readSeqStore));
        for (auto const & superRead : readSet.reads)
            for (auto const readID : superRead.subReads)
                assignments[readID].push_back(superRead.id);

        std::vector<double> bases(readSet.size(), 0);
        for (auto const & superRead : readSet.reads)
            for (auto const readID : superRead.subReads)
                bases[superRead.id] += 1.0 * length(inStore.readSeqStore[readID]) / assignments[readID].size();

        for (unsigned i = 0; i < readSet.size(); ++i)
            avgCoverages.push_back(1.0 * bases[i] / (readSet.reads[i].endPos - readSet.reads[i].beginPos));
    }

    // Set of read ids written out once already.
    std::set<unsigned> seenReadIDs;
    // Counter for next read id to generate.
    unsigned nextID { 0 };

    // Pseudo RNG.
    seqan::Rng<> rng;

    // Output.

    TFragmentStore & outStore;

    // Input

    // Input FeatureReadSet.
    FeatureReadSet const & readSet;
    // Input FragmentStore.
    TFragmentStore const & inStore;

    bool logging;
};

void FragmentStoreSplitter::run()
{
    checkAlignedReadStore();  // Check that the aligned read store is sorted by readID.

    // Begin out store with a copy of inStore.
    outStore = inStore;
    // Determine which output contig maps to which input contig, required for copying out contigs.
    std::vector<unsigned> newToOld;
    for (auto const & read : readSet.reads)
    {
        SEQAN_ASSERT_NOT(read.subReads.empty());
        newToOld.push_back(inStore.alignedReadStore[read.subReads.front()].contigId);
    }
    // Create copies of the contigs.
    clear(outStore.contigStore);
    clear(outStore.contigNameStore);
    // resize(outStore.contigStore, newToOld.size());
    for (unsigned oldID : newToOld)
    {
        appendValue(outStore.contigStore, inStore.contigStore[oldID]);
        std::stringstream ss;
        ss << "contig_" << length(outStore.contigNameStore);
        appendValue(outStore.contigNameStore, ss.str().c_str());
    };

    // Collect all reads where both pairs align on the same contig.
    std::set<std::pair<unsigned, unsigned>> contigPairedReads;
    for (auto const & read : readSet.reads)
    {
        SEQAN_ASSERT(std::is_sorted(read.subReads.begin(), read.subReads.end()));
        for (unsigned readID : read.subReads)
        {
            auto const & matePair = inStore.matePairStore[inStore.readStore[readID].matePairId];
            unsigned otherID = (matePair.readId[0] == readID) ? matePair.readId[1] : matePair.readId[0];
            if (std::binary_search(read.subReads.begin(), read.subReads.end(), otherID))
                contigPairedReads.insert(std::make_pair(matePair.readId[0], matePair.readId[1]));
        }
    }

    // Distribute reads proportionally to coverage on contigs.
    std::vector<double> borders;
    std::vector<unsigned> readToContig(length(inStore.readStore), (unsigned)-1);
    for (auto const & pair : contigPairedReads)
    {
        double sum = 0;
        for (auto const & cID : assignments.at(pair.first))
            sum += avgCoverages[cID];

        unsigned dest = assignments[pair.first].back();
        double y = 0;
        double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double>>(0, sum));
        for (auto const & cID : assignments[pair.first])
        {
            y += avgCoverages[cID];
            if (x < y)
            {
                dest = cID;
                break;
            }
        }

        if (logging)
            std::cerr << "READ0\t" << outStore.readNameStore[pair.first] << "\t" << outStore.readSeqStore[pair.first] << " => " << dest << "\n"
                      << "READ1\t" << outStore.readNameStore[pair.second] << "\t" << outStore.readSeqStore[pair.second] << " => " << dest << "\n";
        outStore.alignedReadStore[pair.first].contigId = dest;
        outStore.alignedReadStore[pair.second].contigId = dest;
    }

    writeOutLinks(rng);
}

void FragmentStoreSplitter::writeOutLinks(seqan::Rng<> & rng)
{
    unsigned numReads = length(inStore.readSeqStore);

    // Collect all reads where both pairs align on different contigs.
    std::set<std::pair<unsigned, unsigned>> contigPairedReads;  // TODO(holtgrew): rename?
    for (auto const & read : readSet.reads)
    {
        SEQAN_ASSERT(std::is_sorted(read.subReads.begin(), read.subReads.end()));
        for (unsigned readID : read.subReads)
        {
            auto const & matePair = inStore.matePairStore[inStore.readStore[readID].matePairId];
            unsigned otherID = (matePair.readId[0] == readID) ? matePair.readId[1] : matePair.readId[0];
            if (!std::binary_search(read.subReads.begin(), read.subReads.end(), otherID))
                contigPairedReads.insert(std::make_pair(matePair.readId[0], matePair.readId[1]));
        }
    }

    // Assign reads proportionally to coverage on output contigs.
    SEQAN_FAIL("Continue here!");


    // Enumerate super reads for each read.
    std::vector<std::set<unsigned>> superReads(numReads);
    for (auto const & read : readSet.reads)
        for (unsigned readID : read.subReads)
        {
            superReads[readID].insert(read.id);
            if (logging)
                std::cerr << "superReads[" << readID << "].insert(" << read.id << ")\n";
        }

    // Determine which read pairs occur in different super reads.
    for (unsigned readID = 0; readID < numReads; ++readID)
    {
        auto const & matePair = inStore.matePairStore[inStore.readStore[readID].matePairId];
        if (matePair.readId[0] == readID)
            continue;  // process each mate pair only once

        unsigned readID0 = matePair.readId[0], readID1 = matePair.readId[1];
        if (superReads.at(readID0).empty() || superReads.at(readID1).empty())
            continue;  // Skip if removed since not significant.

        if (superReads.at(readID0) == superReads.at(readID1))
            continue;  // already written out previously

        // Pick a super read ID for each read for the output.
        unsigned superReadID0 = *superReads.at(readID0).begin();
        {
            double sum = 0;
            for (unsigned superReadID0 : superReads.at(readID0))
                sum += avgCoverages[superReadID0];
            double y = 0;
            double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double>>(0, sum));
            for (auto const & cID : superReads[readID0])
            {
                y += avgCoverages[cID];
                if (x < y)
                {
                    superReadID0 = cID;
                    break;
                }
            }
        }
        unsigned superReadID1 = *superReads.at(readID1).begin();
        {
            double sum = 0;
            for (unsigned superReadID1 : superReads.at(readID1))
                sum += avgCoverages[superReadID1];
            double y = 0;
            double x = pickRandomNumber(rng, seqan::Pdf<seqan::Uniform<double>>(0, sum));
            for (auto const & cID : superReads[readID1])
            {
                y += avgCoverages[cID];
                if (x < y)
                {
                    superReadID1 = cID;
                    break;
                }
            }
        }

        // std::cerr << "READ0*\t" << outStore.readNameStore[readID0] << "\t" << outStore.readSeqStore[readID0] << " => " << dest << "\n";
        // std::cerr << "READ1*\t" << outStore.readNameStore[readID1] << "\t" << outStore.readSeqStore[readID1] << " => " << dest << "\n";

        // Write out to selected contig pair.
        outStore.alignedReadStore[readID0].contigId = superReadID0;
        outStore.alignedReadStore[readID1].contigId = superReadID1;
    }
}

void FragmentStoreSplitter::checkAlignedReadStore() const
{
    typedef decltype(inStore.alignedReadStore[0]) TElemRef;
    auto ltReadID = [](TElemRef lhs, TElemRef rhs) { return (lhs.readId < rhs.readId); };
    SEQAN_CHECK(std::is_sorted(begin(inStore.alignedReadStore, seqan::Standard()),
                               end(inStore.alignedReadStore, seqan::Standard()),
                               ltReadID),
                "Must be sorted by read ID");

    for (unsigned i = 0; i < length(inStore.alignedReadStore); ++i)
        SEQAN_CHECK(inStore.alignedReadStore[i].readId == i, "i = %u, readID = %u",
                    i, inStore.alignedReadStore[i].readId);
}

}  // anonymous namespace

// --------------------------------------------------------------------------
// Function mergeToSuperReads()
// --------------------------------------------------------------------------

void mergeToSuperReads(FeatureReadSet & out, CliqueSet const & cliqueSet, FeatureReadSet const & in, bool logging)
{
    out.numActualReads = in.numActualReads;
    out.numFeatures = in.numFeatures;

    if (logging)
    {
        std::cerr << "INPUT READ SET\n"
                  << "READS\n";
        std::copy(in.reads.begin(), in.reads.end(), std::ostream_iterator<Read>(std::cerr, "\n"));
    }

    // Create super reads.
    for (auto const & clique : cliqueSet.cliques())
    {
        if (logging)
            std::cerr << "processing next clique\n";
        Read outRead(out.size(), clique.contigID, clique.beginPos, clique.endPos);
        std::set<unsigned> ids;
        for (unsigned readID = 0; readID < clique.bitSet.size(); ++readID)
            if (clique.bitSet[readID])
            {
                if (logging)
                    std::cerr << "processing read id: " << readID << "\n";
                
                // Merge read into out read.
                outRead.mergeWithThis(in.reads.at(readID));
            }

        out.reads.push_back(outRead);
    }

    if (logging)
    {
        std::cerr << "OUTPUT READ SET\n"
                  << "READS\n";
        std::copy(out.reads.begin(), out.reads.end(), std::ostream_iterator<Read>(std::cerr, "\n"));
    }
}

// --------------------------------------------------------------------------
// Function performCliqueEnumeration()
// --------------------------------------------------------------------------

bool performCliqueEnumeration(CliqueSet & cliqueSet, FeatureReadSet const & in,
                              ReadSeparatorOptions const & options)
{
    // neighbourhood in [it, itBegin).
    std::vector<Read>::const_iterator itBegin = in.reads.begin(), it = in.reads.begin();

    // Returns true if lhs is truly left of rhs.
    auto leftOf = [](Read const & lhs, Read const & rhs) {
        return std::make_pair(lhs.contigID, lhs.endPos) < std::make_pair(rhs.contigID, rhs.beginPos);
    };

    // Forward as long as no overlaps with *itEnd are possible, at most ot itEnd.
    auto fastForward = [&](std::vector<Read>::const_iterator & it,
                           std::vector<Read>::const_iterator const itEnd) {
        while (leftOf(*it, *itEnd) && (it != itEnd))
            ++it;
    };

    boost::dynamic_bitset<> nbh;
    for (; it != in.reads.end(); ++it)
    {
        nbh.clear();  // clear neighbourhood scratch memory
        nbh.resize(in.size());
        fastForward(itBegin, it);  // forward iterator

        // Build neighbourhood.
        std::for_each(itBegin, it, [&](Read const & read) {
                if (it->overlaps(read, options.minOverlapRate) && !it->conflicts(read))
                    nbh.set(read.id);
            });

        if (options.verbosity >= 3)
        {
            std::cerr << "PROCESSING    " << *it << "\n"
                      << "NEIGHBOURHOOD ";
            for (unsigned i = 0; i < nbh.size(); ++i)
                if (nbh[i])
                    std::cerr << i << ", ";
            std::cerr << "\n";
        }

        cliqueSet.processRead(*it, nbh, in, (options.verbosity >= 3));
    }

    if (options.verbosity >= 3)
    {
        std::cerr << "RESULTING CLIQUES\n";
        for (auto const & clique : cliqueSet.cliques())
            std::cerr << clique << "\n";
    }

    return (cliqueSet.cliques().size() != in.size());
}

// ----------------------------------------------------------------------------
// Function buildAtomicReadSet()
// ----------------------------------------------------------------------------

// Builds a FeatureReadSet where each Read represents the corresponding atomic read from the store.

void buildAtomicReadSet(FeatureReadSet & out,
                        TFragmentStore const & fragStore,
                        LocalVariationStore const & localVarStore)
{
    out.numActualReads = length(fragStore.readSeqStore);
    out.numFeatures = length(localVarStore.columns);

    // Build WeightedFeatureVector for each actual read in fragStore from localVarStore.
    out.reads.resize(out.numActualReads);

    for (auto const & el : fragStore.alignedReadStore)
    {
        auto & read = out.reads[el.readId];
        read.id = el.readId;
        read.contigID = el.contigId;
        read.beginPos = std::min(el.beginPos, el.endPos);
        read.endPos = std::max(el.beginPos, el.endPos);
        read.subReads.push_back(el.readId);
    }

    for (unsigned columnID = 0; columnID < length(localVarStore.columns); ++columnID)  // columnID == featureID
    {
        auto const & column = localVarStore.columns[columnID];
        for (unsigned entryID = 0; entryID < length(column); ++entryID)
        {
            auto const & entry = column[entryID];
            unsigned readID = entry.i1;
            int value = ordValue(entry.i2);
            out.reads[readID].features.insert(WeightedFeature(columnID, value, 1));
        }
    }
}

// ----------------------------------------------------------------------------
// Function buildFeatureReadSet()
// ----------------------------------------------------------------------------

void buildFeatureReadSet(FeatureReadSet & out,
                         TFragmentStore const & fragStore,
                         LocalVariationStore const & localVarStore)
{
    // Set basic properties.
    out.numActualReads = length(fragStore.readSeqStore);
    out.numFeatures = length(localVarStore.columns);

    // Create read entries; we will create one entry for each pair that aligns completely on a contig and one for reach
    // read where the pair is not completely on the same contig.

    // Build atomic read set.
    FeatureReadSet atomicReadSet;
    buildAtomicReadSet(atomicReadSet, fragStore, localVarStore);

    // Compute mate indices.
    seqan::String<unsigned> mateIndices;
    calculateMateIndices(mateIndices, fragStore);

    // Next, build a vector of pairs with indices into the alignedReadStore.  If the second value is (unsigned)-1 then
    // this read is not part of a pair that aligns on the same contig.
    unsigned const INVALID_ID = std::remove_reference<decltype(fragStore.alignedReadStore[0])>::type::INVALID_ID;
    std::vector<std::pair<unsigned, unsigned>> pairs;
    for (unsigned idx0 = 0; idx0 < length(fragStore.alignedReadStore); ++idx0)
    {
        auto const & el0 = fragStore.alignedReadStore[idx0];
        unsigned id0 = el0.readId;

        unsigned idx1 = mateIndices[id0];
        SEQAN_CHECK(idx1 != INVALID_ID, "Must be part of a pair!");
        auto const & el1 = fragStore.alignedReadStore[idx1];

        if (el0.contigId == el1.contigId)  // Case of mapping to same contig.
        {
            if (el0.readId > el1.readId)
                continue;  // only handle once
            pairs.push_back(std::make_pair(idx0, idx1));
        }
        else
        {
            pairs.push_back(std::make_pair(idx0, (unsigned)-1));
        }
    }

    auto const & reads = atomicReadSet.reads;

    // Build read set from this.
    for (auto const & pair : pairs)
    {
        auto const & el = fragStore.alignedReadStore[pair.first];
        Read read(out.reads.size(), el.contigId, std::min(el.beginPos, el.endPos), std::max(el.beginPos, el.endPos));
        read.subReads.push_back(el.readId);
        read.features.mergeWithThis(reads.at(el.readId).features);

        if (pair.second != (unsigned)-1)
        {
            auto const & otherEl = fragStore.alignedReadStore[pair.second];
            SEQAN_ASSERT_EQ(el.contigId, otherEl.contigId);
            read.subReads.push_back(otherEl.readId);
            read.features.mergeWithThis(reads.at(otherEl.readId).features, -1/*=> ignore conflicts*/);
            read.beginPos = std::min(read.beginPos, (int)std::min(otherEl.beginPos, otherEl.endPos));
            read.endPos = std::max(read.endPos, (int)std::max(otherEl.beginPos, otherEl.endPos));
            std::sort(read.subReads.begin(), read.subReads.end());
        }
        else
        {
            auto const & otherEl = fragStore.alignedReadStore[mateIndices[el.readId]];
            read.features.mergeWithThis(reads.at(otherEl.readId).features, -1/*=> ignore conflicts*/);
        }

        // // Always merge with other mate's features.
        // auto const & otherEl = fragStore.alignedReadStore[mateIndices[el.readId]];
        // read.features.mergeWithThis(reads.at(otherEl.readId).features);

        if (!read.features.empty())  // only include reads overlapping with features
            out.reads.push_back(read);
    }

    // Sort reads.
    auto ltCoordinate = [](Read const & lhs, Read const & rhs) {
        return (std::make_tuple(lhs.contigID, lhs.beginPos, lhs.endPos, lhs.id) <
                std::make_tuple(rhs.contigID, rhs.beginPos, rhs.endPos, rhs.id));
    };
    std::sort(out.reads.begin(), out.reads.end(), ltCoordinate);
}

// ----------------------------------------------------------------------------
// Function split()
// ----------------------------------------------------------------------------

// // Removes contig pseudo read again from store that was added in reAlignment().
// void deleteFakeContigRead(TFragmentStore & store)
// {
//     eraseBack(store.alignedReadStore);
//     eraseBack(store.readNameStore);
//     eraseBack(store.readSeqStore);
//     eraseBack(store.readStore);
// }

void split(TFragmentStore & outStore,
           FeatureReadSet const & readSet,
           TFragmentStore const & inStore,
           ReadSeparatorOptions const & options)
{
    bool LOGGING = (options.verbosity >= 3);

    FragmentStoreSplitter splitter(outStore, readSet, inStore, (options.verbosity >= 3));
    splitter.run();

    if (LOGGING)
    {
        std::cerr << "STORE AFTER REP SEP BEFORE REALIGNMENT\n";
        printStore(std::cerr, outStore);
    }

    // Remove differing reads after splitting, also splits a

    for (unsigned i = 0; i < length(outStore.contigStore); ++i)
        // Realign and build consensus.
        reAlignment(outStore, i, 1, 20, false);

    if (LOGGING)
    {
        std::cerr << "STORE AFTER SPLITTING\n";
        printStore(std::cerr, outStore);
    }
}

}  // namespace rep_sep
