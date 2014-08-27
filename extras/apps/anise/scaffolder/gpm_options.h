// ==========================================================================
//                               gpm_options.h
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_GPM_OPTIONS_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_GPM_OPTIONS_H_

namespace scaffolder {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class PathMergingOptions
// ----------------------------------------------------------------------------

// Options for the path merging algorithm.

struct PathMergingOptions
{
    enum
    {
        QUIET = 0,
        NORMAL = 1,
        VERBOSE = 2,
        VERY_VERBOSE = 3
    };

    // Verbosity, also se anonymous enum below.
    int verbosity { QUIET };
    // Multiplicator for sigma.
    int mult { 3 };
    // Path length to use for the the mate edge enumeration in transitive reduction.
    unsigned reductionPathLength { 5 };
    // Number of bases to consider a significant overlap.
    int significantOverlap { 5 };
    // Number of bases to ignore a predicted overlap for.
    int ignoreOverlap { 5 };
    // Ignore links with weight smaller than this.
    double minWeight { 3.0 };

    PathMergingOptions() = default;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace scaffolder

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_GPM_OPTIONS_H_
