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

#ifndef SANDBOX_JAGLA_APPS_ROI_FEATURE_PROJECTION_PROJECT_INTERVAL_H_
#define SANDBOX_JAGLA_APPS_ROI_FEATURE_PROJECTION_PROJECT_INTERVAL_H_

#include <fstream>
#include <list>

#include <seqan/basic.h>
#include <seqan/bed_io.h>
#include <seqan/roi_io.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class IntersectBedOptions
// ----------------------------------------------------------------------------

struct IntersectBedOptions
{
    enum CombinationMode
    {
        PROJECTION,
        INTERSECTION,
        UNION,
        DIFF
    };

    CombinationMode mode;

    int verbosity;

    IntersectBedOptions() : mode(PROJECTION), verbosity(1)
    {}
};

// ----------------------------------------------------------------------------
// Class IntersectBed
// ----------------------------------------------------------------------------

// Intersection of ROI with BED.

class IntersectBed
{
public:
    typedef seqan::BedRecord<seqan::Bed6> TBedRecord;
    typedef seqan::RoiRecord              TRoiRecord;

    // Resulting ROI file is written to out.
    seqan::RoiFileOut & roiFileOut;

    // Configuration for BED intersection.
    IntersectBedOptions options;

    // Currently active BED records, sorted by coordinate (push order).
    std::list<TBedRecord> bedRecords;
    // Currently active ROI records, sorted by coordinate (push order).
    std::list<TRoiRecord> roiRecords;

    // name of current reference.
    seqan::CharString ref;

    // Construct BED intersection object with out stream.
    IntersectBed(seqan::RoiFileOut & roiFileOut, IntersectBedOptions const & options) :
            roiFileOut(roiFileOut), options(options)
    {}

    // Destructor; simply write out remaining data.
    ~IntersectBed()
    {
        finishContig();
    }

    // Push BED record into intersection tool.
    void pushBed(TBedRecord const & bedRecord);

    // Push ROI record into intersection tool.
    void pushRoi(TRoiRecord const & roiRecord);

    // Write out everything from the current contig.
    void finishContig();

    // Write out first BED record with all overlapping ROI records combined with it.
    void processFirstBedRecord();

    // Clear all ROI records that cannot overlap with first BED record any more.
    void cleanRoiRecords();

    // Update connectivity from range in ROI records, depending on the configured combination mode.
    void updateConnectivity(seqan::String<bool> & bitmap,
                            int beginPos,  // position of bitmap[0]
                            seqan::BedRecord<seqan::Bed6> const & bedRecord,
                            std::list<TRoiRecord>::const_iterator const & rangeBegin,
                            std::list<TRoiRecord>::const_iterator const & rangeEnd);

    // Create one ROI from the covered bitmap, counts, BED record and range of ROI records.
    void writeRois(seqan::String<bool> const & bitmap,
                   seqan::String<int> const & counts,
                   int beginPos,  // position of bitmap[0] and counts[0]
                   seqan::BedRecord<seqan::Bed6> const & bedRecord);

    // Write a BED record without overlapping ROI record as ROI record.
    void writeEmptyBed(seqan::BedRecord<seqan::Bed6> const & bedRecord);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_JAGLA_APPS_ROI_FEATURE_PROJECTION_PROJECT_INTERVAL_H_
