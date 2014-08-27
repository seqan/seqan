// ==========================================================================
//                                      ANISE
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

#include "site_data_updater.h"

#include "rep_sep/read.h"
#include "rep_sep/read_set.h"

#include "anise/anise_options.h"
#include "anise/site_data.h"

#include "shared/bam_io_utils.h"

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Function getCigar()
// ----------------------------------------------------------------------------

// Get CIGAR string from alignment (via iterators).

template <typename TContigGapsIter, typename TReadGapsIter>
void getCigar(seqan::String<seqan::CigarElement<>> & cigar,
              TContigGapsIter itCBegin, TContigGapsIter itCEnd,
              TReadGapsIter itRBegin, TReadGapsIter itREnd)
{
    clear(cigar);

    while (itCBegin != itCEnd)
    {
        int len = 0;
        if (!isGap(itCBegin) && !isGap(itRBegin))
        {
            len = std::min(countCharacters(itCBegin), countCharacters(itRBegin));
            appendValue(cigar, seqan::CigarElement<>('M', len));
        }
        else if (!isGap(itCBegin) && isGap(itRBegin))
        {
            len = countGaps(itRBegin);
            appendValue(cigar, seqan::CigarElement<>('D', len));
        }
        else if (isGap(itCBegin) && !isGap(itRBegin))
        {
            len = countGaps(itCBegin);
            appendValue(cigar, seqan::CigarElement<>('I', len));
        }
        else if (isGap(itCBegin) && isGap(itRBegin))
        {
            len = std::min(countGaps(itCBegin), countGaps(itRBegin));
        }
        goFurther(itCBegin, len);
        goFurther(itRBegin, len);
    }

    SEQAN_ASSERT(itRBegin == itREnd);
    (void)itREnd;  // only used in assertion
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
// Function SiteDataUpdater()
// ----------------------------------------------------------------------------

class SiteDataUpdater
{
public:
    SiteDataUpdater(SiteData & siteData,
                    TFragmentStore const & fragStore,
                    rep_sep::FeatureReadSet const & superReads,
                    AniseOptions const & options) :
            siteData(siteData), fragStore(fragStore), superReads(superReads), options(options)
    {}


    void run();

private:

    // Mark records as invalid that are not part of a superRead and older than 2 steps.
    void removeUnplacedReads();
    // Filter out elements from superReads if are not supported sufficiently or too short.
    void filterSuperReads();
    // Update single-end reads.
    void updateSingleEnds();
    // Update paired-end reads that align on the same contig.
    void updatePairedEndSameContig();
    // Update paired-end reads that align on different contigs.
    void updatePairedEndDifferentContig();
    // Update all reads whose ID does not occur in processedIDs.
    void updateRemainingReads();
    // Rebuild siteData.scaffold.refNames.
    void buildRefNames();

    // Input / Output

    // Site data with scaffold and read alignmetns to update.
    SiteData & siteData;
    // FragmentStore with the alignments.
    TFragmentStore const & fragStore;
    // The haplotypes as super read set.
    rep_sep::FeatureReadSet superReads;
    // Configuration.
    AniseOptions const & options;

    // State

    // Set of read ids that have already been written out as primary alignments.
    std::set<unsigned> processedIDs;
    // Collection of BamAlignmentRecord objects added for duplicate copies.
    std::vector<seqan::BamAlignmentRecord> dupRecords;
};

void SiteDataUpdater::run()
{
    removeUnplacedReads();

    // filterSuperReads();  // TODO(holtgrew): Done previously already, remove here!

    // TODO(holtgrew): Need to incorporate use ScaffoldingResult here!
    // updateSingleEnds();
    // updatePairedEndSameContig();
    // updatePairedEndDifferentContig();
    updateRemainingReads();

    // Copy out the duplicate records.
    std::copy(dupRecords.begin(), dupRecords.end(), std::back_inserter(siteData.readSet.bamRecords));

    buildRefNames();
}

void SiteDataUpdater::removeUnplacedReads()
{
    std::set<unsigned> seenIDs;  // seen read ids
    for (auto const & superRead : superReads)
        std::copy(superRead.subReads.begin(), superRead.subReads.end(),
                  std::inserter(seenIDs, seenIDs.end()));

    int stepNo = siteData.state.stepNo;

    for (unsigned readID = 0; readID < siteData.readSet.bamRecords.size(); ++readID)
        if (!seenIDs.count(readID))
        {
            // Check age.
            auto & record = siteData.readSet.bamRecords[readID];
            seqan::BamTagsDict tagsDict(record.tags);
            unsigned idx = 0;
            int mappedStep = 0;
            if (findTagKey(idx, tagsDict, "mS"))
                extractTagValue(mappedStep, tagsDict, idx);
            int age = mappedStep - stepNo;
            int MIN_AGE = 2;
            if (age > MIN_AGE)
                siteData.readSet.disableReadAndMate(readID);
        }
}

void SiteDataUpdater::updateRemainingReads()
{
    // Collect IDs of reads that are assigned to the first contig (both of a pair must be)!
    std::set<unsigned> firstIDs;
    for (unsigned readID : superReads.reads[0].subReads)
        firstIDs.insert(readID);
    std::set<unsigned> toErase;
    for (auto readID : firstIDs)
        if (!firstIDs.count((readID / 2) * 2) || !firstIDs.count((readID / 2) * 2 + 1))
        {
            toErase.insert((readID / 2) * 2);
            toErase.insert((readID / 2) * 2 + 1);
        }
    for (auto readID : toErase)
        firstIDs.erase(readID);
    // Otherwise, assign pair to contig they are first seen on.
    std::map<unsigned, unsigned> rID;
    for (auto const & superRead : superReads)
        for (unsigned readID : superRead.subReads)
        {
            if (firstIDs.count(readID))
            {
                rID[readID] = 0;
            }
            else if (!rID.count((readID / 2) * 2) && !rID.count((readID / 2) * 2 + 1))
            {
                rID[(readID / 2) * 2] = superRead.id;
                rID[(readID / 2) * 2 + 1] = superRead.id;
            }
        }

    // Simply move all reads to a valid contig at a valid position such that there are no more conflicts on
    // reading/writing.
    for (auto const & el : fragStore.alignedReadStore)
        if (!processedIDs.count(el.readId))
        {
            auto & record = siteData.readSet.bamRecords[el.readId];
            clearFlag(record, seqan::BAM_FLAG_UNMAPPED);
            clearFlag(record, seqan::BAM_FLAG_NEXT_UNMAPPED);
            setFlag(record, seqan::BAM_FLAG_ALL_PROPER);
            resize(record.cigar, 1);
            record.cigar[0].operation = 'M';
            record.cigar[0].count = length(record.seq);
            record.rID = record.rNextId = rID.find(el.readId)->second;
            record.beginPos = record.pNext = 0;
            record.tLen = 0;
        }
}

void SiteDataUpdater::filterSuperReads()
{
    eraseIf(superReads.reads, [&](rep_sep::Read const & read) {
            return ((int)read.subReads.size() < options.assemblyMinSupport ||
                    (int)length(fragStore.contigStore[read.contigID]) < options.assemblyMinLength);
        });
}

void SiteDataUpdater::updateSingleEnds()
{
    // TODO(holtgrew): Only required later after pair gap closing.
    for (auto const & superRead : superReads)
        for (auto readID : superRead.subReads)
            if (!hasFlagMultiple(siteData.readSet.bamRecords[readID]))
                SEQAN_FAIL("Cannot process single-end reads at the moment!");
}

// Helper functions for obtaining begin/end positions.
int beginPos(TAlignedReadStoreElement const & el) { return std::min(el.beginPos, el.endPos); }
int endPos(TAlignedReadStoreElement const & el) { return std::max(el.beginPos, el.endPos); }

void SiteDataUpdater::updatePairedEndSameContig()
{
    int rID = 0;  // output rID
    for (auto const & superRead : superReads)
    {
        // Obtain Gaps object for contig.
        typedef typename TFragmentStore::TContigStore const TContigStore;
        typedef typename seqan::Value<TContigStore>::Type TContig;
        typedef typename TFragmentStore::TContigSeq const TContigSeq;
        typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<TContig::TGapAnchors>> TContigGaps;

        // Build output sequence from FragmentStore's contig gaps and contig features.
        unsigned contigID = superRead.contigID;
        TContigGaps contigGaps(fragStore.contigStore[contigID].seq, fragStore.contigStore[contigID].gaps);

        // Get shortcut to BAM records.
        auto & records = siteData.readSet.bamRecords;

        // Collect ids of reads on the current haplotype.
        std::set<unsigned> contigReadIDs;
        for (auto readID : superRead.subReads)
            if (hasFlagMultiple(records[readID]))
                contigReadIDs.insert(readID);

        // Process all reads on this contig.
        for (auto readID : contigReadIDs)
        {
            // Get id of other read in pair, only process if otherID > current, to process each pair only once.
            unsigned matePairID = fragStore.readStore[readID].matePairId;
            int mateNo = getMateNo(fragStore, readID);
            SEQAN_ASSERT_NEQ(mateNo, -1);
            unsigned otherID = fragStore.matePairStore[matePairID].readId[1 - mateNo];
            if (!contigReadIDs.count(otherID) || otherID < readID)
                continue;
            SEQAN_ASSERT_EQ(records[readID].qName, records[otherID].qName);
            // Update readID and otherID such that readID is the first record.
            if (!hasFlagFirst(records[readID]))
                std::swap(readID, otherID);

            // Obtain aligned read element for both reads.
            auto const & elF = fragStore.alignedReadStore[readID];
            SEQAN_ASSERT_EQ(elF.readId, readID);
            auto const & elL = fragStore.alignedReadStore[otherID];
            SEQAN_ASSERT_EQ(elL.readId, otherID);

            // Check whether we process a duplicate, create duplicate records in this case, otherwise mark ids as
            // processed.
            bool processed = processedIDs.count(readID);
            if (processed)
            {
                dupRecords.push_back(records[readID]);
                dupRecords.back().flag |= seqan::BAM_FLAG_SECONDARY;
                dupRecords.push_back(records[otherID]);
                dupRecords.back().flag |= seqan::BAM_FLAG_SECONDARY;
            }
            else
            {
                // Flag read ids as processed.
                processedIDs.insert(readID);
                processedIDs.insert(otherID);
            }

            // Obtain references to the first/second record in the pair.
            auto & recordF = processed ? dupRecords[dupRecords.size() - 2] : records[readID];
            auto & recordL = processed ? dupRecords[dupRecords.size() - 1] : records[otherID];

            // Update coordinate.
            recordF.beginPos = recordL.pNext = toSourcePosition(contigGaps, beginPos(elF));
            recordF.rID = recordF.rNextId = rID;
            recordL.beginPos = recordF.pNext = toSourcePosition(contigGaps, beginPos(elL));
            recordL.rID = recordL.rNextId = rID;
            // Update TLEN value.
            if (recordF.beginPos < recordL.beginPos)
                recordF.tLen = (int)(toSourcePosition(contigGaps, endPos(elL)) -
                                     toSourcePosition(contigGaps, beginPos(elF)));
            else
                recordF.tLen = -(int)(toSourcePosition(contigGaps, endPos(elF)) -
                                     toSourcePosition(contigGaps, beginPos(elL)));
            recordL.tLen = -recordF.tLen;

            // Update alignment in CIGAR string.
            typedef typename std::remove_reference<decltype(fragStore.readSeqStore[0])>::type TReadSeq;
            typedef typename std::remove_reference<decltype(elF.gaps)>::type TReadGapAnchors;
            typedef seqan::Gaps<TReadSeq, seqan::AnchorGaps<TReadGapAnchors> > TReadGaps;
            TReadGaps readGapsF(fragStore.readSeqStore[elF.readId], elF.gaps);
            getCigar(recordF.cigar,
                     iter(contigGaps, beginPos(elF), seqan::Standard()), iter(contigGaps, endPos(elF), seqan::Standard()),
                     begin(readGapsF, seqan::Standard()), end(readGapsF, seqan::Standard()));
            TReadGaps readGapsL(fragStore.readSeqStore[elL.readId], elL.gaps);
            getCigar(recordL.cigar,
                     iter(contigGaps, beginPos(elL), seqan::Standard()), iter(contigGaps, endPos(elL), seqan::Standard()),
                     begin(readGapsL, seqan::Standard()), end(readGapsL, seqan::Standard()));
        }

        // Advance to next contig/haplotype.
        rID += 1;
    }
}

void SiteDataUpdater::updatePairedEndDifferentContig()
{
    std::cerr << "WARNING: Write me!\n";
}

void SiteDataUpdater::buildRefNames()
{
    clear(siteData.scaffold.refNames);

    for (unsigned contigID = 0; contigID < length(siteData.scaffold.seqs); ++contigID)
    {
        std::stringstream ss;
        ss << siteData.scaffold.siteName << "_contig_" << contigID;
        appendValue(siteData.scaffold.refNames, ss.str().c_str());
    }
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function updateSiteData()
// ----------------------------------------------------------------------------

// Update SiteData object (scaffold sequences and aligned reads) from the contigs in the given TFragmentStore and the
// haplotypes stored in superReads.  The names of the scaffold sequences are computed from refNamePrefix (prefix is
// concatenated with the numeric index of the contig).

void updateSiteData(SiteData & siteData,
                    TFragmentStore const & fragStore,
                    rep_sep::FeatureReadSet const & superReads,
                    AniseOptions const & options)
{
    SiteDataUpdater helper(siteData, fragStore, superReads, options);
    helper.run();
}
