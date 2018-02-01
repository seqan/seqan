// ==========================================================================
//                   NGS: Regions of Interest Analysis
// ==========================================================================
// Copyright (c) 2012-2018, Bernd Jagla, Institut Pasteur
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

#include "roi_builder.h"

// ---------------------------------------------------------------------------
// Static Members
// ---------------------------------------------------------------------------

int RoiBuilder::nextId = 0;

// ---------------------------------------------------------------------------
// Member Function RoiBuilder::createRoi()
// ---------------------------------------------------------------------------

// Create a new ROI.
void RoiBuilder::createRoi(seqan::BamAlignmentRecord const & record)
{
    clear(currentProfile);
    clear(connective);
    clear(currentRoi);

    readsInCurrentRoi = 0;

    currentRoi.ref = refNames[record.rID];
    currentRoi.rID = record.rID;
    currentRoi.beginPos = record.beginPos;
    currentRoi.endPos = record.beginPos + getAlignmentLengthInRef(record);
    if (options.strandSpecific)
        currentRoi.strand = hasFlagRC(record) ? '-' : '+';
    else
        currentRoi.strand = '+';  // BED does not allow '.' here.
    currentRoi.len = currentRoi.endPos - currentRoi.beginPos;
    currentRoi.countMax = 0;

    extendRoi(record);
}

// ---------------------------------------------------------------------------
// Member Function RoiBuilder::extendRoi()
// ---------------------------------------------------------------------------

// Extend the current ROI.
void RoiBuilder::extendRoi(seqan::BamAlignmentRecord const & record)
{
    // We will use this for extending the connective array in case of using pairing.
    bool fillConnective = false;

    readsInCurrentRoi += 1;

    // Apply pairing information if configured and there is any and we have the left mate.
    if (options.usePairing && hasFlagMultiple(record) && record.tLen > 0)
    {
        if ((int)currentRoi.len < record.tLen + (record.beginPos - currentRoi.beginPos))
            currentRoi.len = record.tLen + (record.beginPos - currentRoi.beginPos);
        currentRoi.endPos = currentRoi.beginPos + currentRoi.len;
        fillConnective = true;
    }

    // Apply alignment information from record.
    if (currentRoi.endPos < record.beginPos + (int)getAlignmentLengthInRef(record))
    {
        currentRoi.endPos = record.beginPos + getAlignmentLengthInRef(record);
        currentRoi.len = currentRoi.endPos - currentRoi.beginPos;
    }

    // Resize coverage count array and profile, case where mates are not considered.
    if (length(currentRoi.count) < (unsigned)currentRoi.len)
    {
        resize(currentRoi.count, currentRoi.len, 0);
        resize(currentProfile, currentRoi.len, TProfileChar());
        resize(connective, currentRoi.len, fillConnective);
    }

    // Increase counts.
    typedef seqan::Iterator<seqan::String<seqan::CigarElement<> > const, seqan::Rooted>::Type TCigarIt;
    unsigned beginPos = record.beginPos - currentRoi.beginPos;
    unsigned k = 0;
    for (TCigarIt it = begin(record.cigar, seqan::Rooted()); !atEnd(it); goNext(it))
    {
        switch (it->operation)
        {
            case 'M':
            case 'X':
            case '=':
                for (unsigned i = 0; i < it->count; ++i, ++beginPos, ++k)
                {
                    currentProfile[beginPos].count[ordValue(seqan::Dna5(record.seq[k]))] += 1;
                    currentRoi.count[beginPos] += 1;
                    connective[beginPos] = true;
                }
                break;
            case 'D':
                for (unsigned i = 0; i < it->count; ++i, ++beginPos)
                {
                    currentRoi.count[beginPos] += 1;
                    connective[beginPos] = true;
                }
                break;
            case 'I':
                k += 1;
                break;
            case 'N':
                if (options.linkOverSkipped)
                {
                    for (unsigned i = 0; i < it->count; ++i, ++beginPos)
                        connective[beginPos] = true;
                }
                else
                {
                    beginPos += it->count;
                }
                break;
            default:
                break;  // NOP
        }
    }
}

// ---------------------------------------------------------------------------
// Member Function RoiBuilder::pushRecord()
// ---------------------------------------------------------------------------

void RoiBuilder::pushRecord(seqan::BamAlignmentRecord const & record)
{
    if (options.verbosity >= 2 && currentRoi.rID != record.rID)
        std::cerr << "switching to rId == " << record.rID << "\n";
    // static int i = 0;
    // if (options.verbosity >= 2 && ++i)
    // {
    //     if (i % 10000 == 0)
    //         std::cerr << "i == " << i << "\n";
    // }

    if (currentRoi.rID != record.rID || currentRoi.endPos < record.beginPos)
    {
        // The record does not overlap with the current ROI.  Write out current record (if any) and start a new ROI.
        writeCurrentRecord();
        createRoi(record);
    }
    else
    {
        // The record overlaps with the current ROI.  Extend it.
        extendRoi(record);
    }
}

// ---------------------------------------------------------------------------
// Member Function RoiBuilder::writeHeader()
// ---------------------------------------------------------------------------

void RoiBuilder::writeHeader()
{
    seqan::RoiHeader header;
    appendValue(header.extraColumns, "num_reads");
    appendValue(header.extraColumns, "gc_content");
    seqan::writeHeader(roiFileOut, header);
}

// ---------------------------------------------------------------------------
// Member Function RoiBuilder::writeRecord()
// ---------------------------------------------------------------------------

int RoiBuilder::writeRecord(MyRoiRecord & record,
                            seqan::String<TProfileChar> const & profile)
{
    // Compute maximal count.
    for (unsigned i = 1; i < record.len; i++)
        record.countMax = std::max(record.countMax, record.count[i]);

    // Get ROI id.
    std::stringstream roiIdStr;
    roiIdStr << "region_" << nextId++;
    record.name = roiIdStr.str();

    // Add information on number of reads in region.
    char buffer[100];
    snprintf(buffer, 99, "%d", readsInCurrentRoi);
    appendValue(record.data, buffer);

    // Compute C+G content from profile.
    unsigned cg = 0, at = 0;
    for (unsigned i = 0; i < length(profile); ++i)
    {
        seqan::Dna5 x = profile[i];
        cg += (x == 'C' || x == 'G');
        at += (x == 'A' || x == 'T');
    }
    double cgContent = (cg + at == 0) ? 0 : 1.0 * cg / (cg + at);
	cgContent = floor(cgContent * 10000.0) / 10000.0;  // flooring to 4 digits
    snprintf(buffer, 99, "%4.4f", cgContent);
    appendValue(record.data, buffer);

    seqan::writeRecord(roiFileOut, record);
    return 0;
}

// ---------------------------------------------------------------------------
// Member Function RoiBuilder::writeCurrentRecord()
// ---------------------------------------------------------------------------

int RoiBuilder::writeCurrentRecord()
{
    if (currentRoi.beginPos == seqan::RoiRecord::INVALID_POS)
        return 0;  // Do not write out if not set.

    // TODO(holtgrew): Cache that there was no N in cigar string and short circuit to writeRecord()?

    // Try to find a value in connective that is false and extract infixes with consecutive stretches of 1s.
    seqan::String<seqan::Pair<int, int> > pairs;
    SEQAN_ASSERT_NOT_MSG(empty(connective), "ROI cannot be empty.");
    SEQAN_ASSERT_MSG(connective[0], "First element of ROI must be match.");
    int beginPos = 0;
    bool state = true;
    for (unsigned i = 1; i < length(connective); ++i)
    {
        if (state)  // looking for 1/0 change
        {
            if (!connective[i])
            {
                appendValue(pairs, seqan::Pair<int, int>(beginPos, i));
                state = false;
                beginPos = i;
            }
        }
        else  // looking for 0/1 change
        {
            if (connective[i])
            {
                state = true;
                beginPos = i;
            }
        }
    }
    if (state && beginPos != (int)length(connective))
        appendValue(pairs, seqan::Pair<int, int>(beginPos, length(connective)));

    // Easy case of all connective: Early exit.
    if (length(pairs) == 1u && front(pairs).i1 == 0 && back(pairs).i2 == (int)length(connective))
        return writeRecord(currentRoi, currentProfile);

    // Write out connected components.
    MyRoiRecord infixRecord;
    seqan::String<TProfileChar> infixProfile;
    for (unsigned i = 0; i < length(pairs); ++i)
    {
        clear(infixRecord);
        clear(infixProfile);

        infixRecord = currentRoi;
        infixRecord.beginPos += pairs[i].i1;
        infixRecord.endPos = infixRecord.beginPos + pairs[i].i2 - pairs[i].i1;
        infixRecord.len = infixRecord.endPos - infixRecord.beginPos;
        infixRecord.count = infix(currentRoi.count, pairs[i].i1, pairs[i].i2);

        infixProfile = infix(currentProfile, pairs[i].i1, pairs[i].i2);

        int res = 0;
        if ((res = writeRecord(infixRecord, infixProfile)) != 0)
            return res;
    }

    return 0;
}
