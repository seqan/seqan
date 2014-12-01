// ==========================================================================
//                                 BASIL
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

#include "filter_remove_clipping.h"

#include "utils.h"

#include <unordered_map>

// #define BASIL_DEBUG

// ----------------------------------------------------------------------------
// Class RemoveClippingFilterImpl
// ----------------------------------------------------------------------------

class RemoveClippingFilterImpl
{
public:
    RemoveClippingFilterImpl(int maxBeginPosDistance, int minClippingLength) :
            maxBeginPosDistance(maxBeginPosDistance), minClippingLength(minClippingLength),
            maxBufferSize(0)
    {}

    void filter(std::vector<seqan::BamAlignmentRecord *> & out,
                std::vector<seqan::BamAlignmentRecord *> const & in);

    void finish(std::vector<seqan::BamAlignmentRecord *> & out);

private:
    typedef std::list<seqan::BamAlignmentRecord *> TBuffer;
    typedef TBuffer::iterator TBufferIter;
    // List of ExtraInfo/BamAlignmentRecord pairs.  Sorted by coordinate.
    TBuffer buffer;
    // Get position in output buffer from name of the read and is second flag.
    typedef std::unordered_map<std::pair<seqan::CharString, bool>, TBufferIter, HashPair> TBufferMap;
    typedef TBufferMap::iterator TBufferMapIter;
    TBufferMap bufferMap;
    // Maximal distance between the begin positions of a valid pair.  This can be the sum of the maximal fragment size
    // and the maximal alignment length.
    int maxBeginPosDistance;
    int minClippingLength;

    // Process from the beginning of buffer until the distance is greater than maxBeginPosDistance.
    void advance(std::vector<seqan::BamAlignmentRecord *> & out,
                 bool isFinish = false);

    bool spanAboveFragmentSize(seqan::BamAlignmentRecord const & lhs,
                               seqan::BamAlignmentRecord const & rhs,
                               int maxBeginPosDistance) const;

    unsigned maxBufferSize;
};

void RemoveClippingFilterImpl::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                                      std::vector<seqan::BamAlignmentRecord *> const & in)
{
    // Copy in all reads, adding them to the map.
    for (auto ptr : in)
    {
        std::pair<seqan::CharString, bool> key(ptr->qName, hasFlagLast(*ptr));
        // We can only detect the duplicate alignments here if they are close.  They should have been filtered
        // out by the CleanPairsFilter.
        if (bufferMap.find(key) != bufferMap.end())
            SEQAN_FAIL("There are two alignments in the BAM file of the same read %s/%d.",
                       toCString(ptr->qName), (hasFlagLast(*ptr) + 1));
        // Register with key in map.
#ifdef BASIL_DEBUG
        std::cerr << "INTO ANTI CLIPPING\t" << ptr->qName << "/" << (hasFlagLast(*ptr) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
        buffer.push_back(ptr);
        bufferMap[key] = buffer.end();
        --bufferMap[key];

        // Advance in buffers.
        advance(out);
    }

    maxBufferSize = std::max(maxBufferSize, (unsigned)buffer.size());

}

void RemoveClippingFilterImpl::finish(std::vector<seqan::BamAlignmentRecord *> & out)
{
    advance(out, true);
    SEQAN_CHECK(buffer.empty(), "size = %d", (int)buffer.size());
    //fprintf(stderr, "REMOVE CLIPPING MAX BUFFER SIZE\t%d\n", maxBufferSize);
}

void RemoveClippingFilterImpl::advance(std::vector<seqan::BamAlignmentRecord *> & out,
                                       bool isFinish)
{
    while (!buffer.empty() && (isFinish || spanAboveFragmentSize(*buffer.front(), *buffer.back(), maxBeginPosDistance)))
    {
        seqan::BamAlignmentRecord const & record = *buffer.front();  // shortcut

        // Find both mates in the buffer.
        TBufferMapIter itRecord = bufferMap.find(std::make_pair(record.qName, hasFlagLast(record)));
        TBufferMapIter itOther = bufferMap.find(std::make_pair(record.qName, !hasFlagLast(record)));

#ifdef BASIL_DEBUG
        std::cerr << "ANTI CLIPPING HANDLING (" << __LINE__ << ")\t" << record.qName << "/" << (hasFlagLast(record) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG

        SEQAN_CHECK(itRecord != bufferMap.end(), "Should find this read: qname = %s/%d", toCString(record.qName), (hasFlagLast(record) + 1));
        SEQAN_CHECK(itOther != bufferMap.end(), "Should find this read: qname = %s/%d", toCString(record.qName), (!hasFlagLast(record) + 1));

        if (hasFlagUnmapped(**itRecord->second) != hasFlagNextUnmapped(**itRecord->second))
        {
#ifdef BASIL_DEBUG
            std::cerr << "  ANCHOR/SHADOW PAIR\n";
#endif  // #ifdef BASIL_DEBUG

            // If this is a read/OEA pair then we write out the two records directly to the output, with the anchor one
            // leading.
            if (hasFlagUnmapped(**itRecord->second))
                std::swap(itRecord, itOther);
            // *itRecord->second is aligned now, mark as anchor
            seqan::BamTagsDict anchorTags((*itRecord->second)->tags);
            setTagValue(anchorTags, "an", (int)1);
            // Write out and erase from map.
            out.push_back(*itRecord->second);
            out.push_back(*itOther->second);
#ifdef BASIL_DEBUG
            std::cerr << "TO OUT (" << __LINE__ << ")\titRecord: " << (*itRecord->second)->qName << "/" << (hasFlagLast(**itRecord->second) + 1) << "\n";
            std::cerr << "TO OUT (" << __LINE__ << ")\titOther: " << (*itOther->second)->qName << "/" << (hasFlagLast(**itOther->second) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
            buffer.erase(itRecord->second);
            buffer.erase(itOther->second);
            bufferMap.erase(itRecord);
            bufferMap.erase(itOther);
        }
        else if (hasClipping(**itRecord->second, minClippingLength) != hasClipping(**itOther->second, minClippingLength))
        {
#ifdef BASIL_DEBUG
            std::cerr << "  CLIPPED\n";
#endif  // #ifdef BASIL_DEBUG

            // If exactly one of the records is clipped then we can transform it into a shadow record if the clipping is
            // at the outer end of the fragment.
            if (hasClipping(**itRecord->second, minClippingLength))
                std::swap(itRecord, itOther);
            // *itRecord->second is unclipped now.
            // Get shortcuts with better labels (anchor and shadow to be).
            seqan::BamAlignmentRecord & anchor = **itRecord->second;
            seqan::BamAlignmentRecord & shadow = **itOther->second;
            // TODO(holtgrew): The following only works in R+ (--> <--) orientation.
            if ((hasFlagRC(shadow) && startsWithClipping(shadow.cigar, minClippingLength)) ||
                (!hasFlagRC(shadow) && endsWithClipping(shadow.cigar, minClippingLength)))
            {
                // Clipped record is clipped on the wrong side, discard it.
#ifdef BASIL_DEBUG
                std::cerr << "DELETING (" << __LINE__ << ")\titRecord: " << (*itRecord->second)->qName << "/" << (hasFlagLast(**itRecord->second) + 1) << "\n";
                std::cerr << "DELETING (" << __LINE__ << ")\titOther: " << (*itOther->second)->qName << "/" << (hasFlagLast(**itOther->second) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
                delete *itRecord->second;
                delete *itOther->second;
                buffer.erase(itRecord->second);
                buffer.erase(itOther->second);
                bufferMap.erase(itRecord);
                bufferMap.erase(itOther);
                continue;
            }
            // Make shadow unaligned.
            anchor.flag = anchor.flag | seqan::BAM_FLAG_NEXT_UNMAPPED;
            anchor.flag = anchor.flag ^ seqan::BAM_FLAG_ALL_PROPER;
            anchor.pNext = anchor.beginPos;
            anchor.tLen = 0;
            seqan::BamTagsDict anchorTags(anchor.tags);
            setTagValue(anchorTags, "an", (int)1);
            shadow.flag = shadow.flag | seqan::BAM_FLAG_UNMAPPED;
            shadow.flag = shadow.flag ^ seqan::BAM_FLAG_ALL_PROPER;
            shadow.beginPos = anchor.beginPos;
            clear(shadow.cigar);
            shadow.mapQ = 0;
            shadow.tLen = 0;
            // Write out and erase from map.
            out.push_back(*itRecord->second);
            out.push_back(*itOther->second);
#ifdef BASIL_DEBUG
            std::cerr << "TO OUT (" << __LINE__ << ")\titRecord: " << (*itRecord->second)->qName << "/" << (hasFlagLast(**itRecord->second) + 1) << "\n";
            std::cerr << "TO OUT (" << __LINE__ << ")\titOther: " << (*itOther->second)->qName << "/" << (hasFlagLast(**itOther->second) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
            buffer.erase(itRecord->second);
            buffer.erase(itOther->second);
            bufferMap.erase(itRecord);
            bufferMap.erase(itOther);
        }
        else
        {
#ifdef BASIL_DEBUG
            std::cerr << "  DISCARDING (" << __LINE__ << ")\n";
            std::cerr << "  " << record.qName << "/" << (hasFlagLast(record) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG

            // If this is not an anchor/shadow pair and it is no clipped/unclipped pair then discard it.
#ifdef BASIL_DEBUG
            std::cerr << "DELETING (" << __LINE__ << ")\titRecord: " << (*itRecord->second)->qName << "/" << (hasFlagLast(**itRecord->second) + 1) << "\n";
            std::cerr << "DELETING (" << __LINE__ << ")\titOther: " << (*itOther->second)->qName << "/" << (hasFlagLast(**itOther->second) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
            delete *itRecord->second;
            delete *itOther->second;
            buffer.erase(itRecord->second);
            buffer.erase(itOther->second);
            bufferMap.erase(itRecord);
            bufferMap.erase(itOther);
        }
    }
}

bool RemoveClippingFilterImpl::spanAboveFragmentSize(seqan::BamAlignmentRecord const & lhs,
                                                     seqan::BamAlignmentRecord const & rhs,
                                                     int maxBeginPosDistance) const
{
    if (lhs.rID != rhs.rID)
        return true;
    return abs(lhs.beginPos - rhs.beginPos) > maxBeginPosDistance;
}

// ----------------------------------------------------------------------------
// Class RemoveClippingFilter
// ----------------------------------------------------------------------------

RemoveClippingFilter::RemoveClippingFilter(int maxBeginPosDistance, int minClippingDistance) :
        impl(new RemoveClippingFilterImpl(maxBeginPosDistance, minClippingDistance))
{}

RemoveClippingFilter::~RemoveClippingFilter()
{}

void RemoveClippingFilter::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                                  std::vector<seqan::BamAlignmentRecord *> const & in)
{
    impl->filter(out, in);
}

void RemoveClippingFilter::finish(std::vector<seqan::BamAlignmentRecord *> & out)
{
    impl->finish(out);
}
