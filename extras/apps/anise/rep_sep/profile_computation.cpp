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

#include "profile_computation.h"

namespace rep_sep {

// ----------------------------------------------------------------------------
// Function computeProfile()
// ----------------------------------------------------------------------------

void computeProfile(seqan::String<TProfileChar> & profile,
                    unsigned & offset,
                    TFragmentStore const & store,
                    unsigned contigID)
{
    using namespace seqan;

    // Check that the aligned reads are sorted by contig id.
    typedef decltype(store.alignedReadStore[0]) TElemRef;
    auto lessThanContigID = [](TElemRef lhs, TElemRef rhs) { return (lhs.contigId < rhs.contigId); };
    SEQAN_CHECK(std::is_sorted(begin(store.alignedReadStore, seqan::Standard()),
                               end(store.alignedReadStore, seqan::Standard()),
                               lessThanContigID),
                "Aligned reads must be sorted ascendingly by contigID.");

    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore const, Standard>::Type TAlignedReadIter;
    // typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename Position<TReadSeq>::Type TReadPos;

    // typedef typename TFragmentStore::TContigStore TContigStore;
    // typedef typename Value<TContigStore>::Type TContig;
    typedef typename Value<TAlignedReadStore>::Type TAlignedRead;

    clear(profile);

    // Get begin and end position of aligned reads on contig.
    TAlignedReadIter itBegin = lowerBoundAlignedReads(store.alignedReadStore, contigID, SortContigId());
    TAlignedReadIter itEnd = upperBoundAlignedReads(store.alignedReadStore, contigID, SortContigId());

    // Compute smallest and largest position of the alignments in the fragment store.
    TReadPos minPos = maxValue<TReadPos>();
    TReadPos maxPos = 0;
    for (TAlignedReadIter it = itBegin; it != itEnd; ++it)
    {
        minPos = std::min(minPos, (TReadPos)std::min(it->beginPos, it->endPos));
        maxPos = std::max(maxPos, (TReadPos)std::max(it->beginPos, it->endPos));
    }
    if (minPos == maxValue<TReadPos>())
        return;  // Empty.
    offset = minPos;  // First part of result.

    TProfileChar proto;
    resize(profile, maxPos - minPos, proto);

    // Iterate over each read alignment and fill the profile.
    for (TAlignedReadIter it = itBegin; it != itEnd; ++it)
    {
        typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> > TReadGaps;
        typedef typename Iterator<TReadGaps, Standard>::Type TReadGapsIter;

        TReadSeq readSeq = store.readSeqStore[it->readId];
        if (it->endPos < it->endPos)
            reverseComplement(readSeq);
        TReadGaps readGaps(readSeq, it->gaps);

        unsigned pos = std::min(it->beginPos, it->endPos) - offset;
        // std::cerr << it->beginPos << "\t" << it->endPos << "\t" << pos << "\t" << readSeq << "\n";
        for (TReadGapsIter it = begin(readGaps, Standard()), itEnd = end(readGaps, Standard());
             it != itEnd; ++it, ++pos)
        {
            if (isGap(it))
                profile[pos].count[ValueSize<Dna5>::VALUE] += 1;
            else
                profile[pos].count[ordValue(*it)] += 1;
        }
    }
}

void computeProfile(seqan::String<TProfileChar> & profile,
                    TFragmentStore const & store,
                    unsigned contigID)
{
    unsigned ignored = 0;
    (void)ignored;
    computeProfile(profile, ignored, store, contigID);
}

}  // namespace rep_sep
