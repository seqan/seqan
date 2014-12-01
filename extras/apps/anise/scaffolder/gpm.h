// ==========================================================================
//                                  ANISE
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

// TODO(holtgrew): Ordered pair for overlaps not suitable.

// TODO(holtgrew): Overlap checking, best site selection, direction, see below.

// (1) Overlap Checking
//
// The "fitting" of active contigs in the gaps is pretty stupid as implemented below: We assume that the length mean of
// the guiding edge is correct and compute the overlap lengths solely depending on this.  It would be much smarter to
// handle this the same way as the Celera assembler.  The Celera assembler uses statistical tests based on the mean and
// standard deviation to check for the overlap.  Also, we should also add in the overlap lengths instead of just binary
// overlap flags.
//
// Furthermore, overlaps are currently considered unorderedly, only.
//
// (2) Best Site Selection
//
// In case that multiple gaps are possible (large standard deviation), we should take the site that is closest to the
// mean.  This is currently not implemented.
//
// (3) Direction Checking
//
// This is a biggie.  Currently, we assume that the contigs are alredy correctly oriented and all links are conforming
// with this.  The implementation here was guided by (Huson et al., 2002) that uses edges for modeling contigs.  These
// implementation could represent ordering by connecting to the begin/end vertex of the contig edges.
//
// The Celera assembler represents contigs by vertices, however.  It remains to be seen how much more complex the code
// becomes when introducing directions.  Possibly, this could allow putting zippering left and right into one function
// with reverse iterators on lists.

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_GPM_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_GPM_H_

#include <functional>
#include <set>
#include <utility>
#include <vector>

// ============================================================================
// External Forwards
// ============================================================================

namespace assembler { class Overlap; }

namespace scaffolder {

// ============================================================================
// Forwards
// ============================================================================

struct MateLink;
struct PathMergingOptions;
struct ContigEdgeLabel;
struct ScaffoldingResult;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function bundleLinks()
// --------------------------------------------------------------------------

void bundleLinks(std::vector<MateLink> & bundled,
                 std::vector<MateLink> const & links,
                 PathMergingOptions const & options);

std::vector<MateLink> bundleLinks(std::vector<MateLink> const & links,
                                  PathMergingOptions const & options);

// --------------------------------------------------------------------------
// Function transitiveReduction()
// --------------------------------------------------------------------------

void transitiveReduction(std::vector<MateLink> & reduced,
                         std::vector<MateLink> const & links,
                         std::vector<ContigEdgeLabel> const & contigInfos,
                         PathMergingOptions const & options);

// --------------------------------------------------------------------------
// Function greedyPathMerging()
// --------------------------------------------------------------------------

// Function types for (1) checking for existing overlap and (2) computation of the overlap (errors == Overlap::INVALID
// in case of non-existing overlap.
typedef std::function<bool(unsigned seq0,
                           unsigned seq1,
                           unsigned overlapLen,
                           unsigned band)> TOverlapExistsFunc;
typedef std::function<assembler::Overlap(unsigned seq0,
                                         unsigned seq1,
                                         unsigned overlapLen,
                                         unsigned band)> TComputeOverlapFunc;

void greedyPathMerging(ScaffoldingResult & result,
                       std::vector<MateLink> const & links,
                       std::vector<ContigEdgeLabel> const & contigInfos,
                       TOverlapExistsFunc overlapExists,
                       TComputeOverlapFunc computeOverlap,
                       PathMergingOptions const & options);

}  // namespace scaffolder

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_GPM_H_
