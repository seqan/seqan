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
// The ProjectSplicedRoi class allows one to do the projection of ROIs to grouped
// GTF/GFF intervals.
//
// This mode is enabled when the --gff-group-by KEY option is used.  In this
// case, the GFF/GTF records with the same value for KEY are joined together
// in the order of their start position after projecting the ROIs on them.
// This implies that the GFF/GTF records have to be disjoint which is
// generally the case for exons for one transcript in a GFF/GTF file.  The
// type of the records to use can be selected with the --gff-type option.
// ==========================================================================

#ifndef SANDBOX_JAGLA_APPS_ROI_INTERSECT_PROJECT_SPLICED_H_
#define SANDBOX_JAGLA_APPS_ROI_INTERSECT_PROJECT_SPLICED_H_

#include <fstream>

#include <seqan/basic.h>
#include <seqan/gff_io.h>
#include <seqan/roi_io.h>
#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/stream.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Projection of ROI to grouped GFF/GTF intervals.  Used by the
// RoiIntersectApp when --gff-group-by is given to the command line.
//
// The algorithm implied here works reference by reference sequence (the GFF
// records have to be clustered by reference and then sorted by begin
// position.
//
// First, all GFF/GTF records for the contig are read and we compute the
// smallest begin and largest end position for each group.  In a second pass,
// we then do the overlap computation (with the groups, not the GFF/GTF
// records) and projection.  The projections of the ROI to the GFF/GTF records
// is then joined for each group and written as one ROI record to the output.

class ProjectSplicedRoi
{
public:
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::Pair<int, int> TIntPair;

    // Resulting ROI file is written to out.
    seqan::RoiFileOut & roiFileOut;

    // The names of the groups that we have seen so far.
    TNameStore groupNames;
    TNameStoreCache groupNamesCache;
    // The begin/end positions for each of the groups.
    seqan::String<TIntPair> ranges;

    // The group by key.
    seqan::CharString groupBy;

    // Index of the currently active group.
    unsigned currentGroup;

    // The currently active GFF records.
    std::list<seqan::GffRecord> gffRecords;
    // Groups that the GFF records are part of.
    std::list<seqan::StringSet<seqan::CharString> > gffGroups;
    // The currently active ROI records.
    std::list<seqan::RoiRecord> roiRecords;

    int verbosity;

    ProjectSplicedRoi(seqan::RoiFileOut & roiFileOut, seqan::CharString const & groupBy, int verbosity) :
            roiFileOut(roiFileOut), groupNamesCache(groupNames), groupBy(groupBy), currentGroup(0),
            verbosity(verbosity)
    {}

    // Write remaining to output file.
    ~ProjectSplicedRoi();
    // TODO(holtgrew): Write destructor that writes out remaining data.

    // Start a new contig, reset the group names, caches, and ranges for first pass.
    void beginContig();

    // Start second pass, sort group names by their begin position.
    void beginSecondPass();

    // Register a GFF record for updating the ranges.
    void updateRanges(seqan::GffRecord const & record);

    // Helper function that updates a range given a record and a group name.
    void _updateRanges(seqan::GffRecord const & record,
                       seqan::Segment<seqan::CharString, seqan::InfixSegment> const & name);

    // Push a GFF record for overlapping/projection.
    void pushGff(seqan::GffRecord const & record);

    // Push a ROI record for overlapping/projection.
    void pushRoi(seqan::RoiRecord const & record);

    // Write out the ROI for the currently active group.
    void writeCurrentGroup();

    // Remove GFF and ROI records that are to the left of the active group.
    void purgeGffAndRoiRecords();
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_JAGLA_APPS_ROI_INTERSECT_PROJECT_SPLICED_H_
