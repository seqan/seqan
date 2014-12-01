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

#include "local_variation_store.h"

#include "rep_sep/profile_computation.h"
#include "rep_sep/shared_typedefs.h"
#include "rep_sep/rep_sep_options.h"

namespace rep_sep {

namespace {  // anonymous namespace

// --------------------------------------------------------------------------
// Class LocalVariationStoreFiller
// --------------------------------------------------------------------------

// Helper for fill()

class LocalVariationStoreFiller
{
public:
    LocalVariationStoreFiller(LocalVariationStore & varStore,
                              TFragmentStore const & fragStore,
                              ReadSeparatorOptions const & options) :
            varStore(varStore), fragStore(fragStore), options(options)
    {}

    void runForContig(unsigned contigID);

private:
    // Check sorting of aligned reads.
    void checkSortedness() const;

    // Compute consensus, positions, and profile for varStore.
    void initialFill(unsigned contigID, unsigned offset);

    // Collect remaining column information for all deviating columns.
    void fillColumns(unsigned oldSize);

    // State

    // The profile of the current contig.
    seqan::String<TProfileChar> profile;
    
    // Input / Output
    
    // Local variation store.
    LocalVariationStore & varStore;
    // Fragment store.
    TFragmentStore const & fragStore;
    // Configuration for the read separator.
    ReadSeparatorOptions const & options;
};

void LocalVariationStoreFiller::checkSortedness() const
{
    // Check that the aligned reads are sorted by contig id.
    typedef decltype(fragStore.alignedReadStore[0]) TElemRef;
    auto lessThanCoord = [](TElemRef lhs, TElemRef rhs) {
        return (std::make_pair(lhs.contigId, std::min(lhs.beginPos, lhs.endPos)) <
                std::make_pair(rhs.contigId, std::min(rhs.beginPos, rhs.endPos)));
    };
    SEQAN_CHECK(std::is_sorted(begin(fragStore.alignedReadStore, seqan::Standard()),
                               end(fragStore.alignedReadStore, seqan::Standard()),
                               lessThanCoord),
                "Aligned reads must be sorted ascendingly by coordinate.");
}

void LocalVariationStoreFiller::runForContig(unsigned contigID)
{
    // Check for fragment store being sorted properly.
    checkSortedness();

    // Compute the profile.  We will use it to compute the consensus and deviations thereof.
    clear(profile);
    unsigned offset = 0;  // first aligned base on current contig
    computeProfile(profile, offset, fragStore, contigID);
    
    // Store old size of the local variation store to know where to start updating the store.
    unsigned oldSize = length(varStore);

    // Fill store with consensus, positions, and profile of deviating columns.
    initialFill(contigID, offset);

    // Fill remaining info.
    fillColumns(oldSize);
}

void LocalVariationStoreFiller::initialFill(unsigned contigID, unsigned offset)
{
    using namespace seqan;

    typedef typename Iterator<String<TProfileChar>, Standard>::Type TProfileIter;
    unsigned pos = offset;
    for (TProfileIter it = begin(profile, Standard()); it != end(profile, Standard()); ++it, ++pos)
        for (unsigned i = 0; i < 6u; ++i)  // Note: N does not count as deviation, see below.
            if ((i != 4) && (i != _getMaxIndex(*it)) && (int)it->count[i] >= options.minDeviations)
            {
                typedef typename LocalVariationStore::TConsensusAlphabet TConsensusAlphabet;
                appendValue(varStore.consensus, TConsensusAlphabet(_getMaxIndex(*it)));
                appendValue(varStore.positions, std::make_pair(contigID, pos));
                appendValue(varStore.profile, *it);
                break;
            }
    resize(varStore.columns, length(varStore.positions));
    resize(varStore.coveringReads, length(varStore.positions));
}

void LocalVariationStoreFiller::fillColumns(unsigned oldSize)
{
    using namespace seqan;

    // Collect column information for all positions.
    typedef typename LocalVariationStore::TColumnEntry TColumnEntry;
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore const, Standard>::Type TAlignedReadIter;
    typedef typename Iterator<String<std::pair<int, int> > >::Type TPositionIter;

    // We iterate over the aligned reads (sorted by begin position) and the positions array at the same time.
    // TODO(holtgrew): Limit to current contig?
    std::pair<TAlignedReadIter, TAlignedReadIter> aliIt = std::make_pair(begin(fragStore.alignedReadStore, Standard()),
                                                                         begin(fragStore.alignedReadStore, Standard()));
    TAlignedReadIter itAEnd = end(fragStore.alignedReadStore, Standard());
    TPositionIter itP = iter(varStore.positions, oldSize, Standard());
    TPositionIter itPEnd = end(varStore.positions, Standard());

    for (unsigned idx = oldSize; itP != itPEnd; ++idx, ++itP)  // idx is the column in positions/consensus/columns
    {
        // Compute interval of alignments that might overlap with *itP.
        while (aliIt.first != itAEnd && std::make_pair((int)aliIt.first->contigId, (int)std::max(aliIt.first->beginPos, aliIt.first->endPos)) < *itP)
            ++aliIt.first;
        if (aliIt.first == itAEnd)
            break;  // Done, no overlap.
        aliIt.second = aliIt.first;
        while (aliIt.second != itAEnd && std::make_pair((int)aliIt.second->contigId, (int)std::min(aliIt.second->beginPos, aliIt.second->endPos)) <= *itP)
            ++aliIt.second;
        if (aliIt.first != itAEnd)
            SEQAN_ASSERT(std::make_pair((int)aliIt.first->contigId, (int)std::min(aliIt.first->beginPos, aliIt.first->endPos)) <= *itP);
        if (aliIt.second != itAEnd)
            SEQAN_ASSERT(std::make_pair((int)aliIt.second->contigId, (int)std::min(aliIt.second->beginPos, aliIt.second->endPos)) > *itP);

        // For each of these alignments, verify that they actually overlap.  If this is the case then compute the value
        // of the column on this alignment in the row.
        for (TAlignedReadIter it = aliIt.first; it != aliIt.second; ++it)
        {
            if (std::make_pair((int)it->contigId, (int)std::min(it->beginPos, it->endPos)) > *itP ||
                std::make_pair((int)it->contigId, (int)std::max(it->beginPos, it->endPos)) <= *itP)
                continue;  // Skip, does not overlap.

            typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
            // typedef typename Iterator<TAlignedReadStore const, Rooted>::Type TAlignedReadIter;
            // typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
            typedef typename TFragmentStore::TReadSeq TReadSeq;
            // typedef typename Position<TReadSeq>::Type TReadPos;

            // typedef typename TFragmentStore::TContigStore TContigStore;
            // typedef typename Value<TContigStore>::Type TContig;
            typedef typename Value<TAlignedReadStore>::Type TAlignedRead;

            typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> > TReadGaps;
            typedef typename Iterator<TReadGaps, Standard>::Type TReadGapsIter;

            TReadSeq readSeq = fragStore.readSeqStore[it->readId];
            if (it->endPos < it->endPos)
                reverseComplement(readSeq);
            TReadGaps readGaps(readSeq, it->gaps);

            unsigned pos = itP->second - std::min(it->beginPos, it->endPos);  // position in gaps
            TReadGapsIter itRG = iter(readGaps, pos, Standard());
            TColumnEntry val(it->readId, '-', 'I');  // TODO(holtgrew): from neighbourhood?
            if (!isGap(itRG))
            {
                val.i2 = readSeq[toSourcePosition(readGaps, pos)];
                val.i3 = '!' + getQualityValue(readSeq[toSourcePosition(readGaps, pos)]);
            }
            appendValue(varStore.columns[idx], val);

            appendValue(varStore.coveringReads[idx], it->readId);
            // Build compressionMap.
            resize(varStore.compressionMap[it->readId], 1);
            varStore.compressionMap[it->readId][0] = it->readId;
        }
        std::sort(begin(varStore.columns[idx], Standard()), end(varStore.columns[idx], Standard()));
        std::sort(begin(varStore.coveringReads[idx], Standard()), end(varStore.coveringReads[idx], Standard()));
    }
}

}  // anonymous namespace

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

inline void clear(LocalVariationStore & store)
{
    clear(store.positions);
    clear(store.consensus);
    clear(store.profile);
    clear(store.columns);
    clear(store.coveringReads);
    store.compressionMap.clear();
}

// --------------------------------------------------------------------------
// Function fill()
// --------------------------------------------------------------------------

void fill(LocalVariationStore & varStore,
          TFragmentStore const & fragStore,
          ReadSeparatorOptions const & options)
{
    LocalVariationStoreFiller helper(varStore, fragStore, options);
    clear(varStore);
    for (unsigned contigID = 0; contigID < length(fragStore.contigStore); ++contigID)
        helper.runForContig(contigID);
}

// ----------------------------------------------------------------------------
// Function filter()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Reduce size of compression map?

void filter(LocalVariationStore & varStore,
            seqan::String<unsigned> const & columnIds)
{
    using namespace seqan;

    std::set<unsigned> readIDs;

    for (unsigned i = 0; i < length(columnIds); ++i)
    {
        unsigned cID = columnIds[i];

        varStore.positions[i] = varStore.positions[cID];
        varStore.consensus[i] = varStore.consensus[cID];
        varStore.profile[i] = varStore.profile[cID];
        varStore.columns[i] = varStore.columns[cID];
        varStore.coveringReads[i] = varStore.coveringReads[cID];

        readIDs.insert(begin(varStore.coveringReads[i], Standard()),
                       end(varStore.coveringReads[i], Standard()));
    }

    resize(varStore.positions, length(columnIds));
    resize(varStore.consensus, length(columnIds));
    resize(varStore.profile, length(columnIds));
    resize(varStore.columns, length(columnIds));
    resize(varStore.coveringReads, length(columnIds));

    // Update compression map.
    for (std::map<unsigned, String<unsigned> >::iterator it = varStore.compressionMap.begin();
         it != varStore.compressionMap.end(); /* see below */)
    {
        if (!readIDs.count(it->first))
        {
            varStore.compressionMap.erase(it++);
            continue;
        }

        Iterator<String<unsigned>, Standard>::Type it2;
        it2 = std::set_intersection(begin(it->second, Standard()), end(it->second, Standard()),
                                    readIDs.begin(), readIDs.end(),
                                    begin(it->second, Standard()));
        resize(it->second, it2 - begin(it->second, Standard()));

        // Advance to next, required here since we erase during iteration above.
        ++it;
    }
}

}  // namespace rep_sep
