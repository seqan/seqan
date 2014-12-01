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

#include "filter_clean_pairs.h"

#include "utils.h"

#include <queue>
#include <unordered_map>

// #define BASIL_DEBUG

// ----------------------------------------------------------------------------
// Class CleanPairsFilterImpl
// ----------------------------------------------------------------------------

class CleanPairsFilterImpl
{
public:
    CleanPairsFilterImpl(int maxBeginPosDistance) : maxBeginPosDistance(maxBeginPosDistance), maxBufferSize(0)
    {}

    // Read from input list writing to output list.
    //
    // We will buffer internally, so the output might contain an arbitrary number (even 0) of records.
    void filter(std::vector<seqan::BamAlignmentRecord *> & out,
                std::vector<seqan::BamAlignmentRecord *> const & in);

    // Flush everything into the out buffer.
    void finish(std::vector<seqan::BamAlignmentRecord *> & out);

private:
    typedef std::deque<seqan::BamAlignmentRecord *> TBuffer;
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

    // Process from the beginning of buffer until the distance is greater than maxBeginPosDistance.
    void advance(std::vector<seqan::BamAlignmentRecord *> & out,
                 bool isFinish = false);

    // Returns whether the begin position of the two records is above the given maximal begin pos distance.
    bool spanAboveFragmentSize(seqan::BamAlignmentRecord const & lhs,
                               seqan::BamAlignmentRecord const & rhs,
                               int maxBeginPosDistance) const;

    unsigned maxBufferSize;
};

void CleanPairsFilterImpl::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                                  std::vector<seqan::BamAlignmentRecord *> const & in)
{
    for (auto ptr : in)
    {
        // std::cerr << "INTO CLEAN PAIRS\t" << ptr->qName << "/" << (hasFlagLast(*ptr) + 1) << "\n";
        std::pair<seqan::CharString, bool> key(ptr->qName, hasFlagLast(*ptr));
        if (bufferMap.find(key) != bufferMap.end())
        {
            std::cerr << "\nWARNING: There are two alignments in the BAM file of the same read "
                      << ptr->qName << "/" << (hasFlagLast(*ptr) + 1) << " pos: "
                      << (ptr->beginPos + 1) << ".  Ignoring second.\n";
            delete ptr;
        }
        else
        {
            buffer.push_back(ptr);
            bufferMap[key] = buffer.end();
            --bufferMap[key];
        }
    }

    maxBufferSize = std::max(maxBufferSize, (unsigned)buffer.size());

    advance(out);
}

void CleanPairsFilterImpl::finish(std::vector<seqan::BamAlignmentRecord *> & out)
{
    advance(out, true);
    SEQAN_CHECK(buffer.empty(), "Must be empty at end!");
    //fprintf(stderr, "CLEAN PAIRS MAX BUFFER SIZE\t%d\n", maxBufferSize);
}

void CleanPairsFilterImpl::advance(std::vector<seqan::BamAlignmentRecord *> & out,
                                   bool isFinish)
{
    // Process all records where their mate must be in the current buffer.
    //
    // If the current record is the left one then we discard it if its mate is not in the buffer.  If its mate is in the
    // buffer then we write it out and remove the mate from the map.  If the current record is not in the map then its
    // mate must have been written out and thus is automatically written out.
    while (!buffer.empty() && (isFinish || spanAboveFragmentSize(*buffer.front(), *buffer.back(), maxBeginPosDistance)))
    {
        seqan::BamAlignmentRecord const & record = *buffer.front();  // shortcut
        // Stop on orphans.
        if (record.rID == seqan::BamAlignmentRecord::INVALID_REFID)
            break;

        // Handle first case of mate being written out already.
        TBufferMapIter itRecord = bufferMap.find(std::make_pair(record.qName, hasFlagLast(record)));
        if (itRecord == bufferMap.end())
        {
            // If the current record's mate aligns to the right then it cannot have been written out previously:
            // the records are processed in genomic order.
            SEQAN_CHECK(record.beginPos >= record.pNext, "Mate cannot have been written out previously: qname = %s/%d)",
                        toCString(record.qName), (hasFlagLast(record) + 1));

            // Mate has already been written out, splice first element of buffer into outBuffer.
#ifdef BASIL_DEBUG
            std::cerr << "CLEAN PAIRS\tWRITING OUT\t" << buffer.front()->qName << "/" << (hasFlagLast(*buffer.front()) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
            out.push_back(buffer.front());
            buffer.pop_front();
            continue;
        }

        // Otherwise, look whether we can find the other mate.
        TBufferMapIter itOther = bufferMap.find(std::make_pair(record.qName, !hasFlagLast(record)));
        if (itOther == bufferMap.end())
        {
            // Discard this mate.
#ifdef BASIL_DEBUG
            std::cerr << "CLEAN PAIRS\tDELETING\t" << buffer.front()->qName << "/" << (hasFlagLast(*buffer.front()) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
            delete buffer.front();
            buffer.pop_front();
            bufferMap.erase(itRecord);
            continue;
        }

        // If we can find a map entry for the other mate then we write out this record and remove both records from the
        // record map.
#ifdef BASIL_DEBUG
        std::cerr << "CLEAN PAIRS\tWRITING OUT\t" << buffer.front()->qName << "/" << (hasFlagLast(*buffer.front()) + 1) << "\n";
#endif  // #ifdef BASIL_DEBUG
        out.push_back(buffer.front());
        buffer.pop_front();
        bufferMap.erase(itRecord);
        bufferMap.erase(itOther);
    }
}

bool CleanPairsFilterImpl::spanAboveFragmentSize(seqan::BamAlignmentRecord const & lhs,
                                                 seqan::BamAlignmentRecord const & rhs,
                                                 int maxBeginPosDistance) const
{
    if (lhs.rID != rhs.rID)
        return true;
    return abs(lhs.beginPos - rhs.beginPos) > maxBeginPosDistance;
}

// ----------------------------------------------------------------------------
// Class CleanPairsFilter
// ----------------------------------------------------------------------------

CleanPairsFilter::CleanPairsFilter(int maxBeginPosDistance) : impl(new CleanPairsFilterImpl(maxBeginPosDistance))
{}

CleanPairsFilter::~CleanPairsFilter()
{}

void CleanPairsFilter::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                              std::vector<seqan::BamAlignmentRecord *> const & in)
{
    impl->filter(out, in);
}

void CleanPairsFilter::finish(std::vector<seqan::BamAlignmentRecord *> & out)
{
    impl->finish(out);
}
