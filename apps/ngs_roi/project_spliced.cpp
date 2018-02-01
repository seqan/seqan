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

#include "project_spliced.h"

// ---------------------------------------------------------------------------
// Member Function ProjectSplicedRoi::~ProjectSplicedRoi()
// ---------------------------------------------------------------------------

ProjectSplicedRoi::~ProjectSplicedRoi()
{
    while (currentGroup < length(ranges))
    {
        writeCurrentGroup();
        purgeGffAndRoiRecords();
        currentGroup += 1;
    }
}

// ---------------------------------------------------------------------------
// Member Function ProjectedSplicedRoi::beginContig()
// ---------------------------------------------------------------------------

void ProjectSplicedRoi::beginContig()
{
    // Write out all remaining records for this contig.
    while (currentGroup < length(ranges))
    {
        writeCurrentGroup();
        purgeGffAndRoiRecords();
        currentGroup += 1;
    }
    
    // Remove pending ROI records.
    roiRecords.clear();

    // Clear open GFF records and names for the records.
    // TODO(holtgrew): Can we add some checks here that confirm that all has be handled?
    clear(gffRecords);
    gffGroups.clear();

    // Clear group name store and cache and ranges.
    clear(groupNames);
    refresh(groupNamesCache);
    clear(ranges);

    // Set current group to invalid.
    currentGroup = std::numeric_limits<unsigned>::max();
}

// ---------------------------------------------------------------------------
// Member Function ProjectSplicedRoi::beginSecondPass()
// ---------------------------------------------------------------------------

void ProjectSplicedRoi::beginSecondPass()
{
    if (verbosity >= 3)
        for (unsigned i = 0; i < length(groupNames); ++i)
            std::cerr << groupNames[i] << "\t" << ranges[i].i1 << "\t" << ranges[i].i2 << "\n";

    // Sort the group names by begin position of the group (ties are broken by endPosition).
    seqan::String<seqan::Pair<TIntPair, seqan::CharString> > buffer;
    resize(buffer, length(groupNames), seqan::Exact());
    for (unsigned i = 0; i < length(groupNames); ++i)
    {
        buffer[i].i1 = ranges[i];
        buffer[i].i2 = groupNames[i];
    }
    std::sort(begin(buffer, seqan::Standard()), end(buffer, seqan::Standard()));
    // Write back in correct order and refresh group names cache.
    for (unsigned i = 0; i < length(buffer); ++i)
    {
        ranges[i] = buffer[i].i1;
        groupNames[i] = buffer[i].i2;
    }
    refresh(groupNamesCache);

    // Set first group as active.
    currentGroup = 0;
}

// ---------------------------------------------------------------------------
// Member Function ProjectSplicedRoi::_updateRanges()
// ---------------------------------------------------------------------------

void ProjectSplicedRoi::_updateRanges(seqan::GffRecord const & record,
                                      seqan::Segment<seqan::CharString, seqan::InfixSegment> const & name)
{
    if (verbosity >= 3)
        std::cerr << "Updating " << name << "\t" << record.beginPos << "\t" << record.endPos << "\n";

    unsigned idx = 0;
    if (getIdByName(idx, groupNamesCache, name))
    {
        ranges[idx].i1 = std::min(ranges[idx].i1, (int)record.beginPos);
        ranges[idx].i2 = std::max(ranges[idx].i2, (int)record.endPos);
    }
    else
    {
        idx = length(groupNames);
        appendName(groupNamesCache, name);
        appendValue(ranges, TIntPair(record.beginPos, record.endPos));
    }
}

// ---------------------------------------------------------------------------
// Member Function ProjectSplicedRoi::updateRanges()
// ---------------------------------------------------------------------------

void ProjectSplicedRoi::updateRanges(seqan::GffRecord const & record)
{
    // Get group name (possibly comma separated list).
    seqan::CharString groupNames;
    for (unsigned i = 0; i < length(record.tagNames); ++i)
        if (record.tagNames[i] == groupBy)
            groupNames = record.tagValues[i];
    if (empty(groupNames))
        return;  // Record has no group names.

    // Parse out the group names.
    unsigned beginPos = 0;
    unsigned endPos = 0;
    for (; endPos <= length(groupNames); ++endPos)
    {
        if (endPos == length(groupNames) && beginPos < endPos)  // Ignore empty keys.
        {
            _updateRanges(record, infix(groupNames, beginPos, endPos));
        }
        else if (groupNames[endPos] == ',')
        {
            if (beginPos < endPos)  // Ignore empty keys.
                _updateRanges(record, infix(groupNames, beginPos, endPos));
            beginPos = endPos + 1;
        }
    }
}

// ---------------------------------------------------------------------------
// Member Function ProjectSplicedRoi::pushGff()
// ---------------------------------------------------------------------------

void ProjectSplicedRoi::pushGff(seqan::GffRecord const & record)
{
    if (verbosity >= 3)
    {
        std::cerr << "Pushing GFF record.\n  ";
        writeRecord(std::cerr, const_cast<seqan::GffRecord &>(record), seqan::Gff());
    }

    gffRecords.push_back(record);
    // Get string set of group names for the record.
    seqan::StringSet<seqan::CharString> groups;
    for (unsigned i = 0; i < length(record.tagNames); ++i)
        if (record.tagNames[i] == groupBy)
        {
            strSplit(groups, record.tagValues[i], seqan::EqualsChar<','>());
            break;
        }
    gffGroups.push_back(groups);
    if (verbosity >= 3)
    {
        std::cerr << "Groups:";
        for (unsigned i = 0; i < length(groups); ++i)
            std::cerr << " " << groups[i];
        std::cerr << "\n";
    }
}

// ---------------------------------------------------------------------------
// Member Function ProjectSplicedRoi::pushRoi()
// ---------------------------------------------------------------------------

void ProjectSplicedRoi::pushRoi(seqan::RoiRecord const & record)
{
    if (verbosity >= 3)
    {
        std::cerr << "Pushing ROI record.\n";
        std::cerr << "  ";
        writeRecord(std::cerr, record, seqan::Roi());
    }

    roiRecords.push_back(record);

    // If the current group cannot overlap with any ROIs on from record then we can write it out.  This is followed by
    // removing any GFF and ROI records that are used up.
    if (currentGroup < length(ranges) && record.beginPos >= ranges[currentGroup].i2)
    {
        writeCurrentGroup();
        purgeGffAndRoiRecords();
        currentGroup += 1;
    }
}

// ---------------------------------------------------------------------------
// Member Function ProjectSplicedRoi::writeCurrentGroup()
// ---------------------------------------------------------------------------

// Return true iff lhs and rhs overlap (both are BED or ROI records).
template <typename TLeft, typename TRight>
bool overlap(TLeft const & lhs, TRight const & rhs)
{
    return (rhs.beginPos < lhs.endPos) && (lhs.beginPos < rhs.endPos);
}

void ProjectSplicedRoi::writeCurrentGroup()
{
    if (verbosity >= 3)
        std::cerr << "Writing current group (" << currentGroup << ": " << groupNames[currentGroup] << ")\n";

    // Collect all begin/end position pairs for GFF records with the grouping key set to currentName.
    seqan::Reference<TNameStore>::Type currentName = groupNames[currentGroup];
    seqan::String<seqan::Pair<int, int> > pairs;
    typedef std::list<seqan::StringSet<seqan::CharString> >::iterator TGroupsIter;
    TGroupsIter itG = gffGroups.begin();
    typedef std::list<seqan::GffRecord>::iterator TGffIter;
    for (TGffIter it = gffRecords.begin(); it != gffRecords.end(); ++it, ++itG)
    {
        if (std::find(begin(*itG, seqan::Standard()), end(*itG, seqan::Standard()), currentName) == end(*itG))
            continue;  // No match.
        appendValue(pairs, seqan::Pair<int, int>(it->beginPos, it->endPos));
    }

    if (verbosity >= 3)
    {
        std::cerr << "Has the following " << length(pairs) << " GFF intervals\n";
        for (unsigned i = 0; i < length(pairs); ++i)
            std::cerr << "  " << pairs[i].i1 << "\t" << pairs[i].i2 << "\n";
    }

    if (empty(pairs))
        return;  // Nothing to do.

    // Sort the begin/end positions by begin position.
    std::sort(begin(pairs, seqan::Standard()), end(pairs, seqan::Standard()));
    // Compute prefix sums of positions in result.
    seqan::String<int> beginPositions;
    appendValue(beginPositions, 0);
    for (unsigned i = 0; i < length(pairs); ++i)
        appendValue(beginPositions, back(beginPositions) + (pairs[i].i2 - pairs[i].i1));

    if (verbosity >= 3)
    {
        std::cerr << "Begin positions:";
        for (unsigned i = 0; i < length(beginPositions); ++i)
            std::cerr << " " << beginPositions[i];
        std::cerr << "\n";
    }

    // TODO(holtgrew): Check that the intervals in pairs don't overlap?

    // Create resulting ROI.
    seqan::RoiRecord record;
    record.ref = gffRecords.front().ref;
    record.beginPos = front(pairs).i1;
    record.endPos = back(pairs).i2;
    record.strand = gffRecords.front().strand;
    record.name = currentName; // TODO(holtgrew): Make unique.
    record.len = back(beginPositions);
    record.countMax = 0;
    resize(record.count, record.len, 0);

    // Project the ROI counts on the GFF intervals.
    typedef std::list<seqan::RoiRecord>::iterator TRoiIter;
    for (TRoiIter it = roiRecords.begin(); it != roiRecords.end(); ++it)
    {
        for (unsigned i = 0; i < length(pairs); ++i)
        {
            if (verbosity >= 3)
            {
                std::cerr << "ROI record\n  ";
                writeRecord(std::cerr, *it, seqan::Roi());
            }
            if (!(pairs[i].i1 < it->endPos && it->beginPos < pairs[i].i2))
                continue;
            if (verbosity >= 3)
                std::cerr << "=> overlapping\n";

            // Begin and end position of projected interval on contig.
            int beginPosI = std::max(it->beginPos, pairs[i].i1);
            int endPosI = std::min(it->endPos, pairs[i].i2);
            if (beginPosI >= endPosI)
                continue;  // Skip
            // Begin position in record.count.
            int offsetR = beginPositions[i] + beginPosI - pairs[i].i1;
            SEQAN_ASSERT_GEQ(offsetR, 0);
            // Begin position in it->count.
            int offsetC = beginPosI - it->beginPos;
            SEQAN_ASSERT_GEQ(offsetC, 0);

            if (verbosity >= 3)
                std::cerr << ">>> beginPosI = " << beginPosI << "\n"
                          << ">>> endPosI   = " << endPosI << "\n"
                          << ">>> offsetR   = " << offsetR << "\n"
                          << ">>> offsetC   = " << offsetC << "\n"
                          << ">>> beginPositions[i] = " << beginPositions[i] << "\n"
                          << ">>> it->beginPos      = " << it->beginPos << "\n";

            SEQAN_ASSERT_LEQ(offsetR + endPosI - beginPosI, (int)record.len);
            SEQAN_ASSERT_LEQ(offsetC + endPosI - beginPosI, (int)length(it->count));

            if (verbosity >= 3)
            {
                std::cerr << "PROJECTING [" << offsetC << ", " << (offsetC + endPosI - beginPosI) << ") TO "
                          << "[" << offsetR << ", " << offsetR + endPosI - beginPosI << ")\n";
            }

            for (int i = 0, len = endPosI - beginPosI; i < len; ++i)
                record.count[offsetR + i] = it->count[offsetC + i];
        }
    }
    for (unsigned i = 0; i < length(record.count); ++i)
        record.countMax = std::max(record.countMax, record.count[i]);

    writeRecord(roiFileOut, record);
    if (verbosity >= 3)
        std::cerr << "RESULT\t record.ref == " << record.ref << "\n";
}

// ---------------------------------------------------------------------------
// Member Function ProjectSplicedRoi::purgeGffAndRoiRecords()
// ---------------------------------------------------------------------------

void ProjectSplicedRoi::purgeGffAndRoiRecords()
{
    if (verbosity >= 3)
        std::cerr << "Purging GFF and ROI records.\n"
                  << "  current group: " << ranges[currentGroup].i1 << ", " << ranges[currentGroup].i2 << "\n";

    typedef std::list<seqan::GffRecord>::iterator TGffIter;
    typedef std::list<seqan::StringSet<seqan::CharString> >::iterator TGroupsIter;
    TGroupsIter itG = gffGroups.begin();
    for (TGffIter it = gffRecords.begin(); it != gffRecords.end();)
        if ((int)it->endPos <= ranges[currentGroup].i1)
        {
            if (verbosity >= 3)
            {
                std::cerr << "Purging\t";
                writeRecord(std::cerr, *it, seqan::Gff());
            }
            TGffIter itE = it++;
            gffRecords.erase(itE);
            TGroupsIter itEG = itG++;
            gffGroups.erase(itEG);
        }
        else
        {
            ++it;
            ++itG;
        }
    for (TGffIter it = gffRecords.begin(); it != gffRecords.end(); ++it)
        SEQAN_ASSERT_GT((int)it->endPos, ranges[currentGroup].i1);

    typedef std::list<seqan::RoiRecord>::iterator TRoiIter;
    for (TRoiIter it = roiRecords.begin(); it != roiRecords.end();)
        if (it->endPos <= ranges[currentGroup].i1)
        {
            if (verbosity >= 3)
            {
                std::cerr << "Purging\t";
                writeRecord(std::cerr, *it, seqan::Roi());
            }
            TRoiIter itE = it++;
            roiRecords.erase(itE);
        }
        else
        {
            ++it;
        }
}

