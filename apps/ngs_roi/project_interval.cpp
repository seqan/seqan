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
// Author: Bernd Jagla <bernd.jagla@pasteur.fr>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "project_interval.h"

// ----------------------------------------------------------------------------
// Member Functions IntersectBed::pushBed()
// ----------------------------------------------------------------------------

void IntersectBed::pushBed(TBedRecord const & bedRecord)
{
    if (options.verbosity >= 2)
    {
        std::cerr << "pushBed()\n"
                  << "  ";
        writeRecord(std::cerr, bedRecord, seqan::Bed());
    }

    // Handle contig transition, we can unconditionally update rId.
    if (bedRecord.ref != ref)
        finishContig();
    ref = bedRecord.ref;

    // Register new BED record.
    bedRecords.push_back(bedRecord);
}

// ----------------------------------------------------------------------------
// Member Function IntersectBed::writeEmptyBed()
// ----------------------------------------------------------------------------

void IntersectBed::writeEmptyBed(seqan::BedRecord<seqan::Bed6> const & bedRecord)
{
    seqan::RoiRecord roiRecord;
    roiRecord.ref = bedRecord.ref;
    roiRecord.beginPos = bedRecord.beginPos;
    roiRecord.endPos = bedRecord.endPos;
    roiRecord.strand = bedRecord.strand;
    roiRecord.name = bedRecord.name;
    roiRecord.len = roiRecord.endPos - roiRecord.beginPos;
    roiRecord.countMax = 0;
    for (int i = 0; i < (int)roiRecord.len; ++i)
        appendValue(roiRecord.count, 0);

    writeRecord(roiFileOut, roiRecord);
}

// ----------------------------------------------------------------------------
// Member Functions IntersectBed::processFirstBedRecord()
// ----------------------------------------------------------------------------

// Return true iff lhs and rhs overlap (both are BED or ROI records).
template <typename TLeft, typename TRight>
bool overlap(TLeft const & lhs, TRight const & rhs)
{
    return (rhs.beginPos < lhs.endPos) && (lhs.beginPos < rhs.endPos);
}

void IntersectBed::processFirstBedRecord()
{
    if (options.verbosity >= 2)
    {
        std::cerr << "processFirstBedRecord()\n";
        if (!bedRecords.empty())
        {
            std::cerr << "  ";
            writeRecord(std::cerr, bedRecords.front(), seqan::Bed());
        }
    }

    // When computing symmetric difference, we write out the BED record as ROI if there is no overlapping ROI record.
    if (options.mode == IntersectBedOptions::DIFF)
    {
        if (roiRecords.empty() || roiRecords.front().beginPos >= bedRecords.front().endPos)
        {
            writeEmptyBed(bedRecords.front());
            return;
        }
    }

    if (roiRecords.empty())
        return;

    // When in DIFF mode:  Mark all ROI records overlapping with current front BED record as such and stop.
    typedef std::list<TRoiRecord>::/*const_*/iterator TRoiIter;
    if (options.mode == IntersectBedOptions::DIFF)
    {
        for (TRoiIter it = roiRecords.begin(); it != roiRecords.end(); ++it)
            if (overlap(*it, bedRecords.front()))
                back(it->data)[0] = '*';  // Mark as overlapping with BED.
        return;
    }

    // Get range of ROI records overlapping with BED record.
    TRoiIter rangeBegin = roiRecords.begin();
    TRoiIter rangeEnd = rangeBegin;
    for (; rangeEnd != roiRecords.end() && overlap(*rangeEnd, bedRecords.front()); ++rangeEnd)
    {
        if (options.verbosity >= 2)
        {
            std::cerr << "ROI RECORD\t";
            writeRecord(std::cerr, *rangeEnd, seqan::Roi());
        }
        continue;
    }

    // ------------------------------------------------------------------------
    // Create a result ROI from that depending on the configuration.
    // ------------------------------------------------------------------------

    // Compute the smallest begin position and the largest end position of BED and all ROI objects.
    TBedRecord const & bedRecord = bedRecords.front();
    int beginPos = bedRecord.beginPos;
    int endPos = bedRecord.endPos;
    if (options.verbosity >= 2)
        std::cerr << "beginPos, endPos == " << beginPos << ", " << endPos << "\n";
    for (TRoiIter it = rangeBegin; it != rangeEnd; ++it)
    {
        if (options.verbosity >= 2)
            std::cerr << "it->beginPos, it->endPos == " << it->beginPos << ", " << it->endPos << "\n";
        beginPos = std::min(it->beginPos, beginPos);
        endPos = std::max(it->endPos, endPos);
        if (options.verbosity >= 2)
            std::cerr << "beginPos, endPos == " << beginPos << ", " << endPos << "\n";
    }

    // Create bitmap marking positions as connected/covered.
    seqan::String<bool> connected;
    resize(connected, endPos - beginPos, false);
    updateConnectivity(connected, beginPos, bedRecord, rangeBegin, rangeEnd);

    // Create array with counts.
    seqan::String<int> counts;
    resize(counts, endPos - beginPos, 0);
    for (TRoiIter it = rangeBegin; it != rangeEnd; ++it)
        for (unsigned j = it->beginPos - beginPos, i = 0; i < length(it->count); ++i, ++j)
            counts[j] += it->count[i];  // TODO(holtgrew): Increment or simply set?

    // Write resulting intervals to output file.
    writeRois(connected, counts, beginPos, bedRecord);
}

// ----------------------------------------------------------------------------
// Member Functions IntersectBed::cleanRoiRecords()
// ----------------------------------------------------------------------------

void IntersectBed::cleanRoiRecords()
{
    if (options.verbosity >= 2)
    {
        std::cerr << "cleanRoiRecords()\n";
    }

    while (!roiRecords.empty() && (bedRecords.empty() || roiRecords.front().endPos <= bedRecords.front().beginPos))
    {
        // Remove data from front ROI, write it out, remove it.
        if (options.mode == IntersectBedOptions::DIFF && back(roiRecords.front().data)[0] == '.')
        {
            clear(roiRecords.front().data);
            writeRecord(roiFileOut, roiRecords.front());
        }
        roiRecords.pop_front();
    }
}

// ----------------------------------------------------------------------------
// Member Functions IntersectBed::pushRoi()
// ----------------------------------------------------------------------------

void IntersectBed::pushRoi(TRoiRecord const & roiRecord)
{
    if (options.verbosity >= 2)
    {
        std::cerr << "pushRoi()\n"
                  << "  ";
        writeRecord(std::cerr, roiRecord, seqan::Roi());
    }

    // Handle contig transition, we can unconditionally update ref.
    if (roiRecord.ref != ref)
        finishContig();
    ref = roiRecord.ref;

    // Complete processing all BED records that cannot overlap with roiRecord.  After removing each BED record, remove
    // all ROI records that do not overlap with first BED record any more.
    while (!bedRecords.empty() && bedRecords.front().endPos <= roiRecord.beginPos)
    {
        processFirstBedRecord();
        bedRecords.pop_front();
        cleanRoiRecords();
    }

    // Finally, register ROI records.
    roiRecords.push_back(roiRecord);
    // Mark ROI record as not touched yet.  We do this using the data entry.
    appendValue(roiRecords.back().data, ".");
}

// ----------------------------------------------------------------------------
// Member Functions IntersectBed::finishContig()
// ----------------------------------------------------------------------------

void IntersectBed::finishContig()
{
    if (options.verbosity >= 2)
        std::cerr << "finishContig()\n";

    while (!bedRecords.empty())
    {
        processFirstBedRecord();
        bedRecords.pop_front();
    }

    cleanRoiRecords();
}

// ----------------------------------------------------------------------------
// Member Functions IntersectBed::updateConnectivity()
// ----------------------------------------------------------------------------

void IntersectBed::updateConnectivity(seqan::String<bool> & bitmap,
                                      int beginPos,  // position of bitmap[0]
                                      seqan::BedRecord<seqan::Bed6> const & bedRecord,
                                      std::list<TRoiRecord>::const_iterator const & rangeBegin,
                                      std::list<TRoiRecord>::const_iterator const & rangeEnd)
{
    if (options.verbosity >= 2)
        std::cerr << "updateConnectivity()\n";

    typedef std::list<TRoiRecord>::const_iterator TRoiIter;

    switch (options.mode)
    {
        case IntersectBedOptions::PROJECTION:
            // Mark whole BED record range as covered.
            for (int pos = bedRecord.beginPos; pos != bedRecord.endPos; ++pos)
                bitmap[pos - beginPos] = true;
            break;
        case IntersectBedOptions::INTERSECTION:
            // Mark intersection of BED record and ROI records as covered.
            for (TRoiIter it = rangeBegin; it != rangeEnd; ++it)
            {
                for (int pos = it->beginPos; pos != it->endPos; ++pos)
                    if (pos >= bedRecord.beginPos && pos < bedRecord.endPos)
                        bitmap[pos - beginPos] = true;
            }
            break;
        case IntersectBedOptions::UNION:
            // Mark union of BED and ROI intervals as being covered.
            for (unsigned i = 0; i < length(bitmap); ++i)
                bitmap[i] = true;
            break;
        case IntersectBedOptions::DIFF:
        default:
            SEQAN_FAIL("Cannot reach here!");
    }
}

// ----------------------------------------------------------------------------
// Member Functions IntersectBed::writeRois()
// ----------------------------------------------------------------------------

void IntersectBed::writeRois(seqan::String<bool> const & bitmap,
                             seqan::String<int> const & counts,
                             int beginPos,  // position of bitmap[0] and counts[0]
                             seqan::BedRecord<seqan::Bed6> const & bedRecord)
{
    // We will build a string of ranges in bitmap/counts for writing out.  The bitmap is used depending on the
    // combination mode.
    typedef seqan::Pair<int, int> TIntPair;
    seqan::String<TIntPair> pairs;
    // We also generate the names of the new regions.
    seqan::StringSet<seqan::CharString> names;

    switch (options.mode)
    {
        case IntersectBedOptions::PROJECTION:
            // Use BED range only.
            appendValue(pairs, TIntPair(bedRecord.beginPos - beginPos, bedRecord.endPos - beginPos));
            break;
        case IntersectBedOptions::INTERSECTION:
            // Search for 0/1 and 1/0 flanks.
            {
                bool state = bitmap[0];
                int lastFlank = 0;
                for (unsigned i = 1; i < length(bitmap); ++i)
                {
                    if (state)
                    {
                        if (!bitmap[i])  // Found a 1/0 flank.
                        {
                            appendValue(pairs, TIntPair(lastFlank, i));
                            state = false;
                            lastFlank = i;
                        }
                    }
                    else
                    {
                        if (bitmap[i])  // Found a 0/1 flank.
                        {
                            state = true;
                            lastFlank = i;
                        }
                    }
                }
                if (state)
                    appendValue(pairs, TIntPair(lastFlank, length(bitmap)));
            }
            break;
        case IntersectBedOptions::UNION:
            // Use whole range, including overlaps to the left of right.
            appendValue(pairs, TIntPair(0, length(counts)));
            break;
        case IntersectBedOptions::DIFF:
            SEQAN_FAIL("Difference does not implemented yet!");
            break;
        default:
            SEQAN_FAIL("Cannot reach here!");
    }

    if (options.verbosity >= 2)
        std::cerr << "Adding " << length(pairs) << " for BED " << bedRecord.name << "\n";

    for (unsigned i = 0; i < length(pairs); ++i)
    {
        std::stringstream ss;
        ss << bedRecord.name << "_" << i;
        appendValue(names, ss.str());
    }

    seqan::RoiRecord record;
    for (unsigned i = 0; i < length(pairs); ++i)
    {
        clear(record);

        record.ref = bedRecord.ref;
        record.beginPos = beginPos + pairs[i].i1;
        record.endPos = beginPos + pairs[i].i2;
        record.len = record.endPos - record.beginPos;
        record.name = names[i];
        record.strand = bedRecord.strand;
        record.countMax = 0;
        resize(record.count, pairs[i].i2 - pairs[i].i1, 0);
        for (int k = 0, j = pairs[i].i1; j != pairs[i].i2; ++j, ++k)
        {
            record.countMax = std::max(record.countMax, (unsigned)counts[j]);
            record.count[k] = counts[j];
        }

        writeRecord(roiFileOut, record);
    }
}

