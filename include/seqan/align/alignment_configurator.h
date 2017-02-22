// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Runtime configurable alignment settings.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_ALIGNMENT_CONFIGURATOR_H_
#define INCLUDE_SEQAN_ALIGN_ALIGNMENT_CONFIGURATOR_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
//  Class AlignmentConfigurator
// ----------------------------------------------------------------------------

template <typename TScoringScheme, typename TSpec>
class AlignmentConfigurator;

// ----------------------------------------------------------------------------
//  Class DPPolicies
// ----------------------------------------------------------------------------

struct DPSelect
{
    // ----------------------------------------------------------------------------
    //  Enum GapPenalty
    // ----------------------------------------------------------------------------

    enum GapPenalty : uint8_t
    {
        LINEAR,
        AFFINE,
        DYNAMIC,
        AUTOMATIC
    };

    // ----------------------------------------------------------------------------
    //  Class FreeEndGaps
    // ----------------------------------------------------------------------------

    struct FreeEndGaps
    {
        static constexpr uint8_t TOP{1};
        static constexpr uint8_t BOTTOM{2};
        static constexpr uint8_t LEFT{4};
        static constexpr uint8_t RIGHT{8};
    };

    // ----------------------------------------------------------------------------
    //  Enum Traceback
    // ----------------------------------------------------------------------------

    enum Traceback : uint8_t
    {
        SCORE_ONLY,
        BEST_ONE,
        BEST_ALL
    };

    // ----------------------------------------------------------------------------
    //  Enum GapsPlacement
    // ----------------------------------------------------------------------------

    enum GapsPlacement : uint8_t
    {
        LEFT,
        RIGHT
    };

    enum OutputFormat : uint8_t
    {
        ANCHOR_GAPS,
        ARRAY_GAPS,
        ALIGNMENT_GRAPH,
        FRAGMENTS
    };
};

struct DPTraits
{
    // Gocal alignment with linear gap costs.
    struct GlobalLinear
    {
        // The algorithm to choose.
        using TAlgorithmType    = GlobalAlignment_<>;
        // The Gaps to choos
        using TGapType          = LinearGaps;
        // The Band to choose.
        using TBandType         = BandOff;
        // The traceback.
        using TTracebackType    = TracebackOn<TracebackConfig_<SingleTrace, GapsLeft>>;
        // The output to choose.
        using TFormat           = ArrayGaps;
    };

    // Global alignment with affine gap costs.
    struct GlobalAffine : public GlobalLinear
    {
        using TGapType          = AffineGaps;
    };

    // Banded global alignment with linear gap costs.
    struct BandedGlobalLinear : public GlobalLinear
    {
        using TBandType         = BandOn;
    };

    // Banded global alignment with affine gap costs.
    struct BandedGlobalAffine : public BandedGlobalLinear
    {
        using TGapType          = AffineGaps;
    };

    // Local alignment with linear gap costs.
    struct LocalLinear : public GlobalLinear
    {
        using TAlgorithmType    = LocalAlignment_<>;
    };

    // Local alignment with affine gap costs.
    struct LocalAffine : public LocalLinear
    {
        using TGapType          = AffineGaps;
    };
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================
    
    
}  // namespace seqan

#endif  // INCLUDE_SEQAN_ALIGN_ALIGNMENT_CONFIGURATOR_H_
