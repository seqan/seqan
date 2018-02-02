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

#ifndef SANDBOX_JAGLA_APPS_BAM2ROI_ROI_BUILDER_H_
#define SANDBOX_JAGLA_APPS_BAM2ROI_ROI_BUILDER_H_

#include <fstream>

#include <seqan/roi_io.h>
#include <seqan/consensus.h>

#include "record_ext.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class RoiBuilderOptions
// ----------------------------------------------------------------------------

struct RoiBuilderOptions
{
    // Verbosity level.
    int verbosity;

    // Strand specific?
    bool strandSpecific;

    // Whether or not to use pairing information.
    bool usePairing;

    // Whether or not to consider "N" entries in CIGAR strings to link two ROIs.
    bool linkOverSkipped;

    RoiBuilderOptions() : verbosity(0), strandSpecific(false), usePairing(false), linkOverSkipped(false)
    {}

    RoiBuilderOptions(int verbosity, bool strandSpecific, bool usePairing, bool linkOverSkipped) :
            verbosity(verbosity), strandSpecific(strandSpecific), usePairing(usePairing),
            linkOverSkipped(linkOverSkipped)
    {}
};

// ----------------------------------------------------------------------------
// Class RoiBuilder
// ----------------------------------------------------------------------------

class RoiBuilder
{
public:
    // ID of next ROI to be generated.
    // TODO(holtgrew): Is this the way to go?
    static int nextId;

    // The reference sequence names.
    seqan::StringSet<seqan::CharString> refNames;

    // The RoiFileOut to use for writing.
    seqan::RoiFileOut & roiFileOut;

    // The options to use for building the ROI file.
    RoiBuilderOptions options;

    // The next ROI record to be written.
    MyRoiRecord currentRoi;

    // The current profile of the ROI.  Required for C+G content computation.
    typedef seqan::ProfileChar<seqan::Dna5, int> TProfileChar;
    seqan::String<TProfileChar > currentProfile;

    // A bit mask that defines whether a position in the current region is considered to be connective (coverage > 1 in
    // case of single end or !options.linkedOverSkipped and a skipped-bases link/link between two mate pairs in case of
    // options.linkOverSkipped/options.usePairing.
    seqan::String<bool> connective;

    // The number of reads in the current ROI.
    int readsInCurrentRoi;

    RoiBuilder(seqan::RoiFileOut & roiFileOut, RoiBuilderOptions const & options) :
        roiFileOut(roiFileOut), options(options), readsInCurrentRoi(0)
    {}

    ~RoiBuilder()
    {
        writeCurrentRecord();
    }

	// Add a BAM alignment record to the current ROI, create a new one if no overlap with the current ROI
    void pushRecord(seqan::BamAlignmentRecord const & record);

    // write ROI file header
    void writeHeader();

private:
	// Write current roi record to text file
    int writeCurrentRecord();

    // Write the given roi record to text file with the given profile, doing some computation.
    int writeRecord(MyRoiRecord & record,
                    seqan::String<TProfileChar> const & profile);

    // Extend the current ROI.
    void extendRoi(seqan::BamAlignmentRecord const & record);

    // Create a new ROI.
    void createRoi(seqan::BamAlignmentRecord const & record);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_JAGLA_APPS_BAM2ROI_ROI_BUILDER_H_
