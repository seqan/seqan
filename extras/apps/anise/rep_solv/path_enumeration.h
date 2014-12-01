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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SOLV_PATH_ENUMERATION_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SOLV_PATH_ENUMERATION_H_

#include <seqan/sequence.h>

// ============================================================================
// External Forwards
// ============================================================================

namespace rep_sep { struct FeatureReadSet; }

namespace rep_solv {

// ============================================================================
// Forwards
// ============================================================================

class ContigGraph;
struct Options;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// The type to use for contig sequences.
//
// We use Dna5Q as an alphabet so we can soft-mask sequences with PHRED quality 1 (see MASK_QUAL).

typedef seqan::String<seqan::Dna5Q> TContigSeq;

extern int const MASK_QUAL;

// Information about assembled sequence.

struct AssemblyInfo
{
    // Whether the path is anchored to the left.
    bool anchoredLeft { false };
    // Whether the path is anchored to the right.
    bool anchoredRight { false };
    // Whether or not we could create an s-t path with overlap and link info only.
    bool spanning { false };

    AssemblyInfo() = default;
    explicit AssemblyInfo(bool anchoredLeft, bool anchoredRight, bool spanning) :
            anchoredLeft(anchoredLeft), anchoredRight(anchoredRight), spanning(spanning)
    {}
};

// Annotate read for being no / left / right anchor.
enum class AnchorType { NONE, LEFT, RIGHT };

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function enumeratePaths()
// ----------------------------------------------------------------------------

// Build scaffold of ContigGraph cleaned earlier.
//
// The FeatureReadSet coveringSet was used to generate the read sets for contig expansion earlier.  It can be used for
// selecting contigs for "step 0" nodes.

void enumeratePaths(seqan::StringSet<TContigSeq> & outAssembly,
                    std::vector<AssemblyInfo> & outAssemblyInfos,
                    rep_sep::FeatureReadSet & outReadSet,
                    ContigGraph const & cg,
                    seqan::StringSet<TContigSeq> const & contigs,
                    rep_sep::FeatureReadSet const & coveringSet,
                    std::vector<int> const & readBirthStep,
                    std::vector<bool> const & overlapsWithFeature,
                    std::vector<AnchorType> const & anchors,
                    Options const & options);

}  // namespace rep_solv

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SOLV_PATH_ENUMERATION_H_
