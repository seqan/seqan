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

#include "filter_realignment.h"

#include <thread>
#include <atomic>

#include <seqan/align.h>
#include <seqan/parallel.h>
#include <seqan/seq_io.h>

#include "filter_shared_types.h"

namespace  // anonymous
{

// ----------------------------------------------------------------------------
// Function cigarLength()
// ----------------------------------------------------------------------------

int cigarLength(seqan::String<seqan::CigarElement<> > const & str)
{
    int res = 0;
    for (unsigned i = 0; i < length(str); ++i)
        if (str[i].operation == 'M' || str[i].operation == 'I' || str[i].operation == 'X')
            res += str[i].count;
    return res;
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class RealignmentFilterImpl()
// ----------------------------------------------------------------------------

class RealignmentFilterImpl
{
public:
    RealignmentFilterImpl(int maxFragmentSize, BamLibraryInfo::Orientation orientation,
                          seqan::CharString const & refFilename, TBamIOContext & bamIOContext,
                          int numThreads, int chunkSize) :
            maxFragmentSize(maxFragmentSize), orientation(orientation), refFilename(refFilename), bamIOContext(bamIOContext),
            numThreads(numThreads), chunkSize(chunkSize), maxBufferSize(0)
    {
        (void)this->bamIOContext;  // used only for debugging
        (void)this->orientation;  // TODO(holtgrew): Remove, unused.
        init();
    }

    void filter(std::vector<seqan::BamAlignmentRecord *> & out,
                std::vector<seqan::BamAlignmentRecord *> const & in);

    void finish(std::vector<seqan::BamAlignmentRecord *> & out)
    {
        processBuffer(out);
        //fprintf(stderr, "REALIGNMENT MAX BUFFER SIZE\t%d\n", maxBufferSize);
    }

private:
    void init();
    void processBuffer(std::vector<seqan::BamAlignmentRecord *> & out);

    // Write pairs of records from buffer to out with some filtering, deleting unused objects.
    void flushBuffer(std::vector<seqan::BamAlignmentRecord *> & out);

    // Grab elements from buffer using the given atomic index into the buffer until there are none left.
    void realignRecords(std::atomic<unsigned> & atomic);

    typedef std::vector<seqan::BamAlignmentRecord *> TBuffer;

    int maxFragmentSize;
    BamLibraryInfo::Orientation orientation;

    seqan::FaiIndex faiIndex;
    seqan::CharString refFilename;

    TBamIOContext & bamIOContext;

    int numThreads;
    int chunkSize;
    TBuffer buffer;

    unsigned maxBufferSize;

    static const int MIN_MATCH_LEN = 20;  // Taken from BWA
    static const int IGNORE_CLIPPING = 2;  // Ignore that many clipped bases.
};

void RealignmentFilterImpl::realignRecords(std::atomic<unsigned> & atomic)
{
    // Buffers used in the loop below.
    seqan::Dna5String refSeq, readSeq;
    seqan::Align<seqan::Dna5String> align;
    resize(rows(align), 2);

    while (true)
    {
        unsigned idx = atomic_fetch_add(&atomic, 2u);
        if (idx >= buffer.size())
            break;  // done, processed all

        seqan::BamAlignmentRecord & anchor = *buffer[idx];
        seqan::BamAlignmentRecord & shadow = *buffer[idx + 1];

        if (!hasFlagUnmapped(shadow))
            continue;
        SEQAN_ASSERT(!hasFlagUnmapped(anchor));

        // Get reference and read sequence.
        SEQAN_ASSERT_EQ(orientation, BamLibraryInfo::R_PLUS);
        bool directionLeft = hasFlagNextRC(shadow);
        int beginPos = 0;  // TODO(holtgrew): Fix this.
        int endPos = 0;
        if (directionLeft)
        {
            endPos = shadow.beginPos;
            beginPos = std::max(0, endPos - maxFragmentSize + (int)length(shadow.seq));
        }
        else
        {
            beginPos = shadow.beginPos;
            endPos = std::min((int)sequenceLength(faiIndex, shadow.rID), beginPos + maxFragmentSize);
        }
#ifdef BASIL_DEBUG
        std::cerr << "Reading region " << sequenceName(faiIndex, shadow.rID) << ":" << (beginPos + 1) << "-" << endPos << "\n";
#endif  // #ifdef BASIL_DEBUG
        if (readRegion(refSeq, faiIndex, shadow.rID, beginPos, endPos) != 0)
        {
            std::cerr << "ERROR: Could not read region from FASTA file/FAI Index\n";
            exit(1);
        }

        readSeq = shadow.seq;
        SEQAN_ASSERT_EQ(orientation, BamLibraryInfo::R_PLUS);
        if (hasFlagRC(shadow) == hasFlagNextRC(shadow))
            reverseComplement(readSeq);

        // Assign reference and read into Align object.
        clear(row(align, 0));
        setSource(row(align, 0), refSeq);
        clear(row(align, 1));
        setSource(row(align, 1), readSeq);

#ifdef BASIL_DEBUG
        std::cerr << "Local alignment filter\n"
                  << "READ\t" << (*it)->qName << "/" << (hasFlagLast(**it) + 1) << "\n";
        std::cerr << "ref:  " << refSeq << "\n"
                  << "read: " << readSeq << "\n";
#endif  // #ifdef BASIL_DEBUG

        // Perform local alignment.
        seqan::SimpleScore scoringScheme(11, -19, -9, -37);  // taken from BWA, except not penalizing N alignments less
        int score = localAlignment(align, scoringScheme);
        (void)score;
        int leadingBases = toSourcePosition(row(align, 1), 0);
        int trailingBases = length(readSeq) - toSourcePosition(row(align, 1), length(row(align, 1)));
#ifdef BASIL_DEBUG
        std::cerr << "resulting align is\n"
                  << align;
#endif  // #ifdef BASIL_DEBUG
        if (leadingBases > IGNORE_CLIPPING && trailingBases > IGNORE_CLIPPING)
            continue;  // Skip if clipped on both sides.
        // Count number of matching bases.
        int alignedBases = 0;
        seqan::Iterator<seqan::Gaps<seqan::Dna5String> >::Type readGapsIt = begin(row(align, 1));
        seqan::Iterator<seqan::Gaps<seqan::Dna5String> >::Type refGapsIt = begin(row(align, 0));
        for (; readGapsIt != end(row(align, 1)); ++readGapsIt, ++refGapsIt)
            if (!isGap(readGapsIt) && !isGap(refGapsIt))
                alignedBases += (*readGapsIt == *refGapsIt);
        // If we could align enough bases and the clipping is on the right side then we update the record such that both
        // records are aligned again.
#ifdef BASIL_DEBUG
        std::cerr << alignedBases << "\t" << directionLeft << "\t" << leadingBases << "\t" << trailingBases << "\n";
        std::cerr << (alignedBases >= MIN_MATCH_LEN) << "\t"
                  << (directionLeft && trailingBases <= IGNORE_CLIPPING) << "\t"
                  << (!directionLeft && leadingBases <= IGNORE_CLIPPING) << "\n";
#endif  // #ifdef BASIL_DEBUG
        if (alignedBases >= MIN_MATCH_LEN &&
            ((!directionLeft && leadingBases <= IGNORE_CLIPPING) || (directionLeft && trailingBases <= IGNORE_CLIPPING)))
        {
            int pos = beginPos + toSourcePosition(row(align, 0), 0);
            seqan::String<seqan::CigarElement<> > cigar;
            if (leadingBases > 0)
                appendValue(cigar, seqan::CigarElement<>('S', leadingBases));
            typedef seqan::Iterator<seqan::Gaps<seqan::Dna5String>, seqan::Rooted>::Type TGapsIter;
            TGapsIter it0 = begin(row(align, 0), seqan::Rooted());
            TGapsIter it1 = begin(row(align, 1), seqan::Rooted());
            while (!atEnd(it0) && !atEnd(it1))
            {
                int len = 0;
                if (!isGap(it0) && !isGap(it1))
                {
                    len = std::min(countCharacters(it0), countCharacters(it1));
                    appendValue(cigar, seqan::CigarElement<>('M', len));
                }
                else if (!isGap(it0) && isGap(it1))
                {
                    len = countGaps(it1);
                    appendValue(cigar, seqan::CigarElement<>('D', len));
                }
                else if (isGap(it0) && !isGap(it1))
                {
                    len = countGaps(it0);
                    appendValue(cigar, seqan::CigarElement<>('I', len));
                }
                else if (isGap(it0) && isGap(it1))
                {
                    len = std::min(countGaps(it0), countGaps(it1));
                }
                goFurther(it0, len);
                goFurther(it1, len);
            }
            if (trailingBases > 0)
                appendValue(cigar, seqan::CigarElement<>('S', trailingBases));

            // Update record and mate.
            // std::cerr << "ref\t" << beginPos << "\t" << endPos << "\t" << refSeq << "\n";
            // std::cerr << "read\t" << readSeq << "\n";
            // std::cerr << "alignment\n" << align;
            // std::cerr << "leading, aligned, trailing\t" << leadingBases << ", " << alignedBases << ", " << trailingBases << "\n";
            // std::cerr << "pos == " << pos << "\n";
            // std::cerr << "score == " << score << "\n";

#ifdef BASIL_DEBUG
            std::cerr << "before update\n";
            write2(std::cerr, anchor, bamIOContext, seqan::Sam());
            write2(std::cerr, shadow, bamIOContext, seqan::Sam());
#endif  // #ifdef BASIL_DEBUG

            anchor.pNext = pos;
            anchor.flag ^= seqan::BAM_FLAG_NEXT_UNMAPPED;
            anchor.flag |= seqan::BAM_FLAG_ALL_PROPER;
            shadow.flag ^= seqan::BAM_FLAG_UNMAPPED;
            shadow.flag |= seqan::BAM_FLAG_ALL_PROPER;
            shadow.beginPos = pos;
            shadow.cigar = cigar;
            shadow.seq = readSeq;

            if (anchor.beginPos < shadow.beginPos)
                anchor.tLen = shadow.beginPos + cigarLength(shadow.cigar) - anchor.beginPos;
            else
                anchor.tLen = shadow.beginPos - cigarLength(anchor.cigar) - anchor.beginPos;
            shadow.tLen = -anchor.tLen;

#ifdef BASIL_DEBUG
            std::cerr << "after update\n";
            write2(std::cerr, anchor, bamIOContext, seqan::Sam());
            write2(std::cerr, shadow, bamIOContext, seqan::Sam());
#endif  // #ifdef BASIL_DEBUG
        }
        else
        {
#ifdef BASIL_DEBUG
            std::cerr << "not updating\n";
            write2(std::cerr, anchor, bamIOContext, seqan::Sam());
            write2(std::cerr, shadow, bamIOContext, seqan::Sam());
#endif  // #ifdef BASIL_DEBUG
        }
    }
}

void RealignmentFilterImpl::processBuffer(std::vector<seqan::BamAlignmentRecord *> & out)
{

    SEQAN_CHECK(buffer.size() % 2 == 0u, "Pairs, length must be even.");

    // Align all records with multiple threads.
    std::atomic<unsigned> nextIdx(0);
    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; ++i)
        threads.push_back(std::thread([&,this]{ realignRecords(nextIdx); }));

    // Wait for all threads.
    for (auto & t : threads)
        t.join();

    // Remove pairs where both mates align now, i.e. the previous shadow now has less than IGNORE_CLIPPING clipped bases
    // on either side.  Note that we ignored the pairs where the shadow is aligned with clipping on both sides.
    flushBuffer(out);
}

void RealignmentFilterImpl::init()
{
    if (read(faiIndex, toCString(refFilename)) != 0)
        throw BamFilterException("Problem reading FAI index.");
}

void RealignmentFilterImpl::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                                   std::vector<seqan::BamAlignmentRecord *> const & in)
{
    std::copy(in.begin(), in.end(), std::back_inserter(buffer));
    maxBufferSize = std::max(maxBufferSize, (unsigned)buffer.size());

    // If there is a sufficient number of elements in the buffer, split it and process it in parallel.
    if (buffer.size() > (unsigned)(numThreads * chunkSize))
        processBuffer(out);
}

void RealignmentFilterImpl::flushBuffer(std::vector<seqan::BamAlignmentRecord *> & out)
{
    for (TBuffer::iterator it = buffer.begin(); it != buffer.end(); )
    {
        TBuffer::iterator itAnchor = it++;
        TBuffer::iterator itShadow = it++;

        if (empty((*itShadow)->cigar))  // remains full shadow, write to out
        {
            out.push_back(*itAnchor);
            out.push_back(*itShadow);
            continue;
        }

        bool frontClipping = (front((*itShadow)->cigar).operation == 'S' && (int)front((*itShadow)->cigar).count > IGNORE_CLIPPING);
        bool backClipping = (back((*itShadow)->cigar).operation == 'S' && (int)back((*itShadow)->cigar).count > IGNORE_CLIPPING);
        if (frontClipping && backClipping)
        {
            delete *itAnchor;
            delete *itShadow;
        }
        else
        {
            out.push_back(*itAnchor);
            out.push_back(*itShadow);
        }
    }

    // Clear buffer.
    buffer.clear();
}

// ----------------------------------------------------------------------------
// Class RealignmentFilter()
// ----------------------------------------------------------------------------

RealignmentFilter::RealignmentFilter(int maxFragmentSize, BamLibraryInfo::Orientation orientation,
                                     seqan::CharString const & refFilename, TBamIOContext & bamIOContext,
                                     int numThreads, int chunkSize) :
        impl(new RealignmentFilterImpl(maxFragmentSize, orientation, refFilename, bamIOContext, numThreads, chunkSize))
{}

RealignmentFilter::~RealignmentFilter()
{}

void RealignmentFilter::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                               std::vector<seqan::BamAlignmentRecord *> const & in)
{
    impl->filter(out, in);
}

void RealignmentFilter::finish(std::vector<seqan::BamAlignmentRecord *> & out)
{
    impl->finish(out);
}
