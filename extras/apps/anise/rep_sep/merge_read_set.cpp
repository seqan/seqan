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

#include "merge_read_set.h"

#include <algorithm>
#include <list>
#include <map>

#include <seqan/basic.h>

namespace rep_sep {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Class IndependentReadMerger
// ----------------------------------------------------------------------------

// Helper for mergeNonConflictingReadsOnContig(), simple heuristic via linear scan.

class IndependentReadMerger
{
public:
    typedef std::vector<Read>::const_iterator TIterator;

    IndependentReadMerger(FeatureReadSet const & in,
                          FeatureReadSet const & initialReadSet,
                          std::vector<int> const & readBirthStep,
                          ReadSeparatorOptions const & options) :
            input(in), initialReadSet(initialReadSet), readBirthStep(readBirthStep), options(options)
    {
        sortInput();
        buildInitialReadID();
    }

    // Perform the merging.
    FeatureReadSet run();

private:
    // Sort the input read set.
    void sortInput();
    // Fill initialReadID();
    void buildInitialReadID();

    // Merge the reads for one contig.
    void mergeContigReads(FeatureReadSet & out, TIterator itBegin, TIterator itEnd) const;
    // For all super reads on one contig after mergeContigReads(), distribute reads arbitrarily as long as there are 0
    // conflicts.
    void arbitraryReadDistribution(std::vector<Read>::iterator itBegin,
                                   std::vector<Read>::iterator itEnd) const;

    // Copy of the input, we will need to sort by coordinate.
    FeatureReadSet input;
    // Initial read set, used for conflict check.
    FeatureReadSet const & initialReadSet;
    // Birth step for each read.
    std::vector<int> const & readBirthStep;
    // Configuration.
    ReadSeparatorOptions const & options;

    // Mapping from readID to initialReadSet.readSet.
    std::vector<unsigned> initialReadID;
};

void IndependentReadMerger::buildInitialReadID()
{
    initialReadID.resize(initialReadSet.numActualReads, (unsigned)-1);
    unsigned idx = 0;
    for (auto const & read : initialReadSet.reads)
    {
        for (unsigned readID : read.subReads)
            initialReadID[readID] = idx;
        ++idx;
    }
}

void IndependentReadMerger::mergeContigReads(FeatureReadSet & out,
                                             IndependentReadMerger::TIterator itBegin,
                                             IndependentReadMerger::TIterator itEnd) const
{
    if (itBegin != itEnd)
        SEQAN_CHECK(std::all_of(itBegin, itEnd, [itBegin](Read const & read) {
                    return (read.contigID == itBegin->contigID); }),
            "All must be on same contig.");

    // Add all iterators as uncovered.
    std::list<TIterator> uncovered;
    for (auto it = itBegin; it != itEnd; ++it)
        uncovered.push_back(it);

    // Count step 0 generation reads.
    std::map<TIterator, int> numFirstGeneration;
    for (auto it : uncovered)
        for (auto readID : it->subReads)
            if (readBirthStep.at(readID) == 0)
                numFirstGeneration[it] += 1;

    // Sort iterators by (num first gen reads, begin pos).
    uncovered.sort([&numFirstGeneration](TIterator const lhs, TIterator const rhs) {
            return (std::make_pair(-numFirstGeneration[lhs], lhs->beginPos) <
                    std::make_pair(-numFirstGeneration[rhs], rhs->beginPos));
        });

    // Scan over uncovered reads (sorted by (num first gen reads, begin positions)) and merge into final read set.
    while (!uncovered.empty())
    {
        Read read(**uncovered.begin());  // always take first
        read.id = out.reads.size();
        uncovered.pop_front();

        for (auto it = uncovered.begin(); it != uncovered.end(); /*see below*/)
            if (read.conflicts(**it) <= options.readMergeIgnoreConflicts)
            {
                read.mergeWithThis(**it, options.readMergeIgnoreConflicts, true);
                it = uncovered.erase(it);
            }
            else
            {
                ++it;
            }
        out.reads.push_back(read);
    }
}

void IndependentReadMerger::arbitraryReadDistribution(
        std::vector<Read>::iterator itBegin,
        std::vector<Read>::iterator itEnd) const
{
    // Distribute the elementary reads from *it to *it2 that do not conflict with *it2.
    for (auto it = itBegin; it != itEnd; ++it)
        for (auto it2 = itBegin; it2 != itEnd; ++it2)
            if (it != it2)
                for (auto readID : it->subReads)
                {
                    auto const & actualRead = initialReadSet.reads.at(initialReadID.at(readID));
                    // std::cerr << "CONFLICTS BETWEEN (#=" << it2->conflicts(actualRead) << ")\n"
                    //           << "\t" << actualRead << "\n"
                    //           << "\t" << *it2 << "\n";
                    if (!it2->conflicts(actualRead))
                        it2->mergeWithThis(actualRead, 0, true);
                }
}

FeatureReadSet IndependentReadMerger::run()
{
    FeatureReadSet result;

    auto itBegin = input.reads.begin();
    while (itBegin != input.reads.end())
    {
        // Get end of contig subrange.
        auto itEnd = itBegin;
        for (; itEnd != input.reads.end() && itEnd->contigID == itBegin->contigID; ++itEnd)
            continue;

        unsigned oldSize = result.size();
        mergeContigReads(result, itBegin, itEnd);
        arbitraryReadDistribution(result.reads.begin() + oldSize, result.reads.end());

        // Start anew for next contig.
        itBegin = itEnd;
    }

    return result;
}

void IndependentReadMerger::sortInput()
{
    std::sort(input.reads.begin(), input.reads.end(), ltCoordinate);
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function mergeNonConflictingReadsOnContig()
// ----------------------------------------------------------------------------

void mergeNonConflictingReadsOnContig(FeatureReadSet & out,
                                      FeatureReadSet const & in,
                                      FeatureReadSet const & initialReadSet,
                                      std::vector<int> const & readBirthStep,
                                      ReadSeparatorOptions const & options)
{
    IndependentReadMerger helper(in, initialReadSet, readBirthStep, options);
    out = helper.run();
}

}  // namespace rep_sep
