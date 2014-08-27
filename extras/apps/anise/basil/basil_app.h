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

#ifndef SANDBOX_HERBARIUM_APPS_BASIL_BASIL_APP_H_
#define SANDBOX_HERBARIUM_APPS_BASIL_BASIL_APP_H_

#include <stdexcept>
#include <ostream>

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

#include "library_info.h"
#include "bam_filter_pipeline.h"
#include "oea_clustering.h"
#include "clipping_clustering.h"

// ============================================================================
// Forwards
// ============================================================================

// Forward for the implementaiton of BasilApp.
class BasilAppImpl;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BasilOptions
// ----------------------------------------------------------------------------

// Top level configuration of BASIL.

struct BasilOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    static const int QUIET = 0;
    static const int NORMAL = 1;
    static const int VERBOSE = 2;
    static const int VERY_VERBOSE = 3;

    // -----------------------------------------------------------------------
    // Input / Output
    // -----------------------------------------------------------------------

    // Path to input BAM/SAM files.
    seqan::CharString inputFile;
    // Path to reference FASTA file.
    seqan::CharString referenceFile;
    // Path to output file. "-" for stdout.
    seqan::CharString outVcfFile;
    // Path to directory with debug output files.
    seqan::CharString outputDebugDir;

    // -----------------------------------------------------------------------
    // Library Info Detection
    // -----------------------------------------------------------------------

    // Library info as given from command line.  Automatically determined if not given.
    BamLibraryInfo libraryInfo;
    // Factor to get from std dev to allowed fragment size.
    double fragmentSizeFactor;
    // This flag is used to tell that the library information was given on the command line.
    bool autoLibraryInfo;
    // The number of records to read from the mapping file for automatically determining the.
    int autoLibraryNumRecords;

    // Computes the maximal fragment size from the configured options.
    int maxFragmentSize() const
    {
        return libraryInfo.median + fragmentSizeFactor * libraryInfo.stdDev;
    }

    // -----------------------------------------------------------------------
    // Coverage Filter
    // -----------------------------------------------------------------------

    // Variants calls in a region with a coverage higher than this number are ignored.  Set to 0 to disable ignoring.
    int filterMaxCoverage;
    // Ignore anchor alignments with a coverage below this value.
    int filterMinAlignmentQuality;

    // -----------------------------------------------------------------------
    // Output Information / Output Filter
    // -----------------------------------------------------------------------

    // Identifier of the individual.
    seqan::CharString individualName;

    // -----------------------------------------------------------------------
    // Clipping Cluster BasilOptions
    // -----------------------------------------------------------------------

    // Maximal alignment length, used for reordering positions.
    int maxAlignmentLength;
    // Radius around clipping position to use for cluster computation.
    int clippingWindowRadius;
    // Smallest length of clipped sequence to consider.
    int clippingMinLength;
    // Minimal number of clipping positions to fall into a window around a clipping position for it to be called as a
    // candidate.
    int clippingMinCoverage;

    // -----------------------------------------------------------------------
    // Realignment Options
    // -----------------------------------------------------------------------

    // Number of threads to use for the realignment.
    int realignmentNumThreads;

    // -----------------------------------------------------------------------
    // Breakpoint Labeling Options
    // -----------------------------------------------------------------------

    // Radius for windows around breakpoints called by clipping scanner.
    int breakpointWindowRadius;

    // -----------------------------------------------------------------------
    // OEA-Clustering Configuration
    // -----------------------------------------------------------------------

    // Which approach to use for the OEA cluster selection.
    OeaClusterSelection oeaClusterSelection;
    // Minimal number of OEA anchors for a tentative site to be considered.
    int oeaMinSupport;
    // Minimal number of OEA anchors on each side for a tentative site to be considered.
    int oeaMinSupportEachSide;

    BasilOptions() :
            verbosity(NORMAL), fragmentSizeFactor(3), autoLibraryInfo(true), autoLibraryNumRecords(1000), filterMaxCoverage(0),
            filterMinAlignmentQuality(0), individualName("INDIVIDUAL"), maxAlignmentLength(0), clippingWindowRadius(0),
            clippingMinLength(0), clippingMinCoverage(0), realignmentNumThreads(0), breakpointWindowRadius(0),
            oeaClusterSelection(OeaClusterSelection::CHAINING),
            oeaMinSupport(0), oeaMinSupportEachSide(0)
    {}

    // Write to given out stream.
    void print(std::ostream & out);

    // Parse options from command line.  Result can be one of OK, ERROR, OTHER.
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const ** argv);
};

// ----------------------------------------------------------------------------
// Typedef BasilAppException
// ----------------------------------------------------------------------------

typedef std::runtime_error BasilAppException;

// ----------------------------------------------------------------------------
// Class BasilApp
// ----------------------------------------------------------------------------

class BasilApp
{
public:
    BasilOptions options;

    BasilApp();
    ~BasilApp();  // for pimpl

    // Entry point for the applications.
    int run(int argc, char const ** argv);

private:
    std::unique_ptr<BasilAppImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_BASIL_BASIL_APP_H_
