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

#include "read_correction.h"

#include "anise/store_utils.h"
#include "anise/bam_tag_names.h"

#include "rep_sep/feature_map.h"
#include "rep_sep/rep_sep_options.h"

namespace rep_sep {

namespace {  // anonymous namespace

// ----------------------------------------------------------------------------
// Class ReadCorrector
// ----------------------------------------------------------------------------

// Helper class for read correction.

class ReadCorrector
{
public:
    ReadCorrector(std::vector<seqan::BamAlignmentRecord> & records,
                  TFragmentStore & fragStore,
                  FeatureMap const & featureMap,
                  TIsReadCorrectable isReadCorrectable) :
            records(records), fragStore(fragStore), featureMap(featureMap), isReadCorrectable(isReadCorrectable)
    {}

    void run();

private:
    // Correct read pointed to by el.
    void correctRead(TAlignedReadStoreElement & el);
    // Return range of features within the given range.
    std::pair<FeatureMap::TIterator, FeatureMap::TIterator> getFeatureRange(
            unsigned contigID, unsigned beginPos, unsigned endPos) const;
    // Save original sequence for given readID in tag "oS".
    void saveOriginalSequence(unsigned readID);

    // The to be corrected read set.
    std::vector<seqan::BamAlignmentRecord> & records;
    // The FragmentStore to use for consensus and read alignment in assembly.x
    TFragmentStore & fragStore;
    // The FeatureMap to use for feature locations.
    FeatureMap const & featureMap;
    // Function returning whether the record holds a correctable read.
    TIsReadCorrectable isReadCorrectable;
};

void ReadCorrector::run()
{
    double const MAX_ERROR_RATE = 0.08;  // maximal distance to consensus
    for (auto & el : fragStore.alignedReadStore)
    {
        int maxErrors = ceil(MAX_ERROR_RATE * length(fragStore.readSeqStore[el.readId]));
        if (isReadCorrectable(records.at(el.readId)) && countErrors(fragStore, el) < maxErrors)
            correctRead(el);
    }
}

void ReadCorrector::correctRead(TAlignedReadStoreElement & el)
{
    // Shortcuts to begin and end position of read alignment in store.
    auto beginPos = std::min(el.beginPos, el.endPos);
    auto endPos = std::max(el.beginPos, el.endPos);
    // Obtain range of features.
    auto featureRange = getFeatureRange(el.contigId, beginPos, endPos);

    // Obtain Gaps object for the contig and read.
    typedef typename TFragmentStore::TContigStore const TContigStore;
    typedef typename seqan::Value<TContigStore>::Type TContig;
    typedef typename TFragmentStore::TContigSeq const TContigSeq;
    typedef seqan::Gaps<TContigSeq, seqan::AnchorGaps<TContig::TGapAnchors>> TContigGaps;

	typedef typename TFragmentStore::TReadSeqStore const TReadSeqStore;
	typedef typename seqan::Value<TReadSeqStore>::Type TReadSeq;
    typedef seqan::Gaps<TReadSeq, seqan::AnchorGaps<seqan::String<TFragmentStore::TReadGapAnchor>>> TReadGaps;

    TContigGaps contigGaps(fragStore.contigStore[el.contigId].seq, fragStore.contigStore[el.contigId].gaps);
    TReadGaps readGaps(fragStore.readSeqStore[el.readId], el.gaps);

    seqan::Dna5String seq;  // rebuild sequence

    auto pos = beginPos;
    auto itC = iter(contigGaps, beginPos);
    auto itR = begin(readGaps, seqan::Standard());
    auto itF = featureRange.first;
    std::vector<decltype(pos)> gapPositions;  // non-consensus gap positions
    for (; pos != endPos; ++pos, ++itC, ++itR)
    {
        SEQAN_ASSERT_NOT(itC == end(contigGaps, seqan::Standard()));
        SEQAN_ASSERT_NOT(itR == end(readGaps, seqan::Standard()));
        if (itF != featureRange.second)
            SEQAN_ASSERT_GEQ(itF->pos, (int)pos);

        if (itF == featureRange.second || itF->pos != (int)pos)
        {
            // We are not on a feature, use consensus.
            if (!isGap(itC))
                appendValue(seq, *itC);
        }
        else
        {
            // We are on a feature, use read value.
            if (!isGap(itR))
                appendValue(seq, *itR);
            else
                gapPositions.push_back(pos);
            ++itF;
        }
    }

    if (records[el.readId].seq == seq)
        return;  // nothing to update

    // Update BAM records.
    saveOriginalSequence(el.readId);
    records[el.readId].seq = seq;
    resize(records[el.readId].qual, length(seq), 'I');  // TODO(holtgrew): Use Dna5Q above instead of this hack.

    // Update store seq and rebuild gaps.
    fragStore.readSeqStore[el.readId] = seq;
    clearGaps(readGaps);
    pos = beginPos;
    for (auto itC = iter(contigGaps, beginPos, seqan::Rooted()),
                 itR = begin(readGaps, seqan::Rooted()); pos < endPos;
         ++pos, ++itC, ++itR)
        if (isGap(itC) || std::binary_search(gapPositions.begin(), gapPositions.end(), pos))
            insertGap(itR);

    // Remove leading and trailing gaps.
    auto itFirst = begin(readGaps, seqan::Standard());
    while (isGap(itFirst))
    {
        removeGap(itFirst);
        beginPos += 1;
        itFirst = begin(readGaps, seqan::Standard());
    }
    auto itLast = iter(readGaps, length(readGaps) - 1, seqan::Standard());
    while (isGap(itLast))
    {
        auto tmp = itLast;
        --itLast;
        removeGap(tmp);
        endPos -= 1;
    }
    bool flip = (el.beginPos > el.endPos);
    el.beginPos = beginPos;
    el.endPos = endPos;
    if (flip)
        std::swap(beginPos, endPos);
}

void ReadCorrector::saveOriginalSequence(unsigned readID)
{
    // TODO(holtgrew): "oS" should come from configuration.
    seqan::BamAlignmentRecord & record = records.at(readID);
    seqan::BamTagsDict dict(record.tags);
    unsigned idx = 0;
    if (findTagKey(idx, dict, KEY_ORIG_SEQ))
        return;  // already stored original sequence
    (void)idx;  // unused
    setTagValue(dict, KEY_ORIG_SEQ, toCString(record.seq));
    idx = 0;
    if (findTagKey(idx, dict, KEY_ORIG_QUAL))
        return;  // already stored original qual string 
    setTagValue(dict, KEY_ORIG_QUAL, toCString(record.qual));
}

std::pair<FeatureMap::TIterator, FeatureMap::TIterator> ReadCorrector::getFeatureRange(
            unsigned contigID, unsigned beginPos, unsigned endPos) const
{
    auto contigRange = featureMap.contigFeatures(contigID);
    FeatureDescription valL(/*id=*/0, FeatureDescription::COLUMN,
                            contigID, beginPos, 0);
    FeatureDescription valU(/*id=*/0, FeatureDescription::COLUMN,
                            contigID, endPos, 0);
    return std::make_pair(std::lower_bound(contigRange.first, contigRange.second, valL),
                          std::upper_bound(contigRange.first, contigRange.second, valU));
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function performReadCorrection()
// ----------------------------------------------------------------------------

void performReadCorrection(std::vector<seqan::BamAlignmentRecord> & records,
                           TFragmentStore & store,
                           FeatureMap const & featureMap,
                           TIsReadCorrectable isReadCorrectable)
{
    // std::cerr << "CORRRECTING\n";
    // featureMap.print(std::cerr);
    // printStore(std::cerr, store);

    ReadCorrector helper(records, store, featureMap, isReadCorrectable);
    helper.run();

    // TODO(holtgrew): Perform realignment to get rid of all-gaps columns?
}

}  // namespace rep_sep
