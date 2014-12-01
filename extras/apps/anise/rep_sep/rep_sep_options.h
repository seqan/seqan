// ==========================================================================
//                                   ANISE
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_REP_SEP_OPTIONS_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_REP_SEP_OPTIONS_H_

#include <string>

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadSeparatorOptions
// ----------------------------------------------------------------------------

struct ReadSeparatorOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity { 1 };

    enum Method
    {
        TAMMI_SIMPLE,  // simple Tammi version using the same error probability
        TAMMI_PHRED    // Tammi version with per-base qualities
    };

    // The method to use for separating column detection.
    Method tammiMethod { TAMMI_SIMPLE };

    // Smallest number of deviations
    int minDeviations { 2 };

    // Average per-base error for simple Tammi method.
    double pErr { 0.01 };

    // Lower bound on number of common reads that a column must have to be considered for separating column test
    // (Kuchenbecker's d_min).
    unsigned minCommonReads { 2 };

    // The largest probability that a correlation occurs by changes in Tammi's model (p^{tot}_{max}).
    double maxRandomCorrelation { 0.001 };

    // The constants \tau_\min and r_\min in Algorithm 4 from Kuchenbecker, 2011.
    int tauMin { 10000 };
    int rMin { 10000 };

    // The minimal percentage of reads to cover per contig to be considered in set selection.
    int minCoverPercentage { 5 };

    // The minimal overlap between local and global classes to consider a global class (\Delta_\max from p. 39 of
    // Kuchenbecker, 2011).
    unsigned minOverlap { 2 };

    // The number of reads for one site to start compressing at.
    unsigned startCompressionAt { 100 };

    // Maximal column size (against super high read stacks).
    unsigned maxColumnSize { 200 };

    // Split at d_min deviations instead of using separating column method.
    bool splitDMin { false };

    // Prefix to use for cloned read.
    std::string readNamePrefix { "cloned_read_" };

    // Percentage of shorter read length for overlap detection.
    double minOverlapRate { 0.4 };

    // The method to split clusters further
    enum ClusterSplitMethod
    {
        SEP_COLUMNS,  // search for separating columns
        BOUND         // search for columns with >= d_min columns.
    };
    ClusterSplitMethod splitMethod { SEP_COLUMNS };

    // Number of conflicts to ignore in conflict free read merging at the end.
    int readMergeIgnoreConflicts { 2 };

    ReadSeparatorOptions() = default;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_REP_SEP_OPTIONS_H_
