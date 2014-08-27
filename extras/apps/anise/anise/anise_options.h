// ==========================================================================
//                                  ANISE
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_ANISE_OPTIONS_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_ANISE_OPTIONS_H_

#include <ostream>
#include <seqan/sequence.h>

#include "rep_sep/rep_sep_options.h"
#include "library_info.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class AniseOptions
// ----------------------------------------------------------------------------

// Top level configuration of ANISE.

struct AniseOptions
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    static const int QUIET = 0;
    static const int NORMAL = 1;
    static const int VERBOSE = 2;
    static const int VERY_VERBOSE = 3;

    // Number of threads to use.
    int numThreads;

    // Verbatim command line.
    std::string commandLine;

    // id of site to debug, assembly is only run for this site, -1 to disable.
    int debugSiteID;
    // step no to debug for
    int debugStepNo;

    // Enable/disable auto-tuning of overlap length/error rate and read mapping error rate.  ANISE will try to set these
    // values in case of too small values (i.e. too small overlap, number of allowed errors in case of short Illumina
    // reads).
    bool autoTuning;

    // -------------------------------------------------------------------------
    // I/O Related Options
    // -------------------------------------------------------------------------

    // Input FASTA file with the reference.
    seqan::CharString inputReference;
    // Input VCF file with insert sites.
    seqan::CharString inputVcf;
    // Input read mapping SAM/BAM file.
    seqan::CharString inputMapping;
    // Output FASTA file with contigs.
    seqan::CharString outputFasta;
    // Output SAM file with mapping of reads to contigs.
    seqan::CharString outputMapping;

    // Directory for temporary files such as orphans or the step BAM files.
    seqan::CharString tmpDir;

    // Directory for debug output.
    seqan::CharString outputDebugDir;

    // Whether or not to remove files iteratively and clean up after program termination.
    bool cleanUpTemporaryFiles;

    // -------------------------------------------------------------------------
    // Algorithm Options
    // -------------------------------------------------------------------------

    // Maximum number of recursion steps.
    int recursionMaxSteps;
    // Enable/disable realignment after assembly.
    bool realignAssembly;  // TODO(holtgrew): Remove, always realigning.
    // If there are more than this number of reads for a site in the initial round then the site is deactivated.
    int stopInitialReadCount;
    // If there are more than this number of reads for a site then the site is deactivated.
    int stopReadCount;
    // If there are more than this number of reads for a site then the triplet extension is not performed.
    int stopTripletExtensionReadCount;
    // If the length sum of all reads for a site divided by the length sum of its contigs is higher than this value
    // before assembly then the site is deactivated.
    int stopCoverage;
    // Bandwidth to use for realignment.
    int realignmentBandwidth;
    // Border around alignments to extract from the profile.
    int realignmentBorder;

    // -----------------------------------------------------------------------
    // Read SeparationOptions
    // -----------------------------------------------------------------------
    // Whether or not to use repeat separation algorithm.
    bool separateRepeats;
    // Read separation options.
    rep_sep::ReadSeparatorOptions readSepOptions;

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

    // -------------------------------------------------------------------------
    // Assembly Options
    // -------------------------------------------------------------------------

    // Radius around the insertion sites to use for collecting BAM records.
    int assemblySiteWindowRadius;
    // Radius around the insertion site to collect clippings in.
    int assemblySiteFringeRadius;

    // -------------------------------------------------------------------------
    // Read Mapping Options
    // -------------------------------------------------------------------------

    // The maximal error rate of alignments to find in read mapping step.
    double readMappingErrorRate;
    // The number of reads to map in one batch.
    int readMappingBatchSize;

    // -------------------------------------------------------------------------
    // Overlapper Options
    // -------------------------------------------------------------------------

    // The minimal overlap of two reads to have, in percent of the longer read length.
    int overlapperMinOverlapRatio;
    // The maximal error rate that the overlap must have.
    int overlapperMaxErrorRate;

    // -------------------------------------------------------------------------
    // MSA Options
    // -------------------------------------------------------------------------

    // The affine scoring scheme to use for the multiple sequence alignment when constructing the alignment graph.
    int msaScoreMatch;
    int msaScoreMismatch;
    int msaScoreGapOpen;
    int msaScoreGapExtend;

    // -------------------------------------------------------------------------
    // Consensus Calling Options
    // -------------------------------------------------------------------------

    // Whether or not to perform read correction during assembly.
    bool readCorrection { false };

    // Minimal number of bases in consensus calling to call a non-N base.
    int consensusMinBaseSupport;
    // Minimal length of a contig such that it is not ignored in consensus calling, expressed in percent of average read
    // length.
    int consensusMinContigLengthRate;
    // Minimal support for a scaffold in terms of reads.
    int minScaffoldSupport;  // TODO(holtgrew): Add command line parameter? Print option!
    // Remove reads with more than this rate of errors.
    double minKeptReadsErrorRate;  // TOOD(holtgrew): Add command line parameter? Print option!
    // Only write out best contig for each site.
    bool onlyWriteOutBest { true };  // TODO(holtgrew): Add command line parameter? Print option!

    // -------------------------------------------------------------------------
    // Configuration for repeat resolution.
    // -------------------------------------------------------------------------

    // Minimal number of reads on contig after unitigging.
    int assemblyMinSupport;
    // Minimal length of contig after unitigging.
    int assemblyMinLength;

    // -------------------------------------------------------------------------
    // Scaffolder Options
    // -------------------------------------------------------------------------

    int minOverlapMergeLen;
    int maxOverlapMergeError;

    // We initialize all members here but defaults come from command line parsing.
    AniseOptions() :
        verbosity(0),
        numThreads(1),
        debugSiteID(-1),
        debugStepNo(-1),
        autoTuning(true),
        cleanUpTemporaryFiles(false),
        recursionMaxSteps(0),
        realignAssembly(true),
        stopInitialReadCount(0),
        stopReadCount(0),
        stopTripletExtensionReadCount(0),
        stopCoverage(0),
        realignmentBandwidth(30),
        realignmentBorder(10),
        separateRepeats(true),
        fragmentSizeFactor(3),
        autoLibraryInfo(true),
        autoLibraryNumRecords(1000),
        assemblySiteWindowRadius(0),
        assemblySiteFringeRadius(0),
        readMappingErrorRate(5.0),
        readMappingBatchSize(0),
        overlapperMinOverlapRatio(0),
        overlapperMaxErrorRate(0),
        msaScoreMatch(0),
        msaScoreMismatch(0),
        msaScoreGapOpen(0),
        msaScoreGapExtend(0),
        consensusMinBaseSupport(0),
        consensusMinContigLengthRate(0),
        minScaffoldSupport(2),
        minKeptReadsErrorRate(0.05),
        assemblyMinSupport(3),
        assemblyMinLength(150),
        minOverlapMergeLen(10),
        maxOverlapMergeError(20)
    {}

    // Write to given out stream.
    void print(std::ostream & out) const;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_ANISE_OPTIONS_H_
