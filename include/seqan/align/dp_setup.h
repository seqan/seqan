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

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_SETUP_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_SETUP_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// SubstituteAlignConfig_
// ----------------------------------------------------------------------------

template <typename TAlignConfig>
struct SubstituteAlignConfig_;

// 0000
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, false, false, false, TSpec> >
{
    typedef FreeEndGaps_<False, False, False, False> Type;
};

// 0001
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, false, false, true, TSpec> >
{
    typedef FreeEndGaps_<False, False, True, False> Type;
};

// 0010
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, false, true, false, TSpec> >
{
    typedef FreeEndGaps_<False, False, False, True> Type;
};


// 0011
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, false, true, true, TSpec> >
{
    typedef FreeEndGaps_<False, False, True, True> Type;
};


// 0100
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, true, false, false, TSpec> >
{
    typedef FreeEndGaps_<False, True, False, False> Type;
};


// 0101
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, true, false, true, TSpec> >
{
    typedef FreeEndGaps_<False, True, True, False> Type;
};


// 0110
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, true, true, false, TSpec> >
{
    typedef FreeEndGaps_<False, True, False, True> Type;
};


// 0111
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, true, true, true, TSpec> >
{
    typedef FreeEndGaps_<False, True, True, True> Type;
};


// 1000
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, false, false, false, TSpec> >
{
    typedef FreeEndGaps_<True, False, False, False> Type;
};


// 1001
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, false, false, true, TSpec> >
{
    typedef FreeEndGaps_<True, False, True, False> Type;
};


// 1010
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, false, true, false, TSpec> >
{
    typedef FreeEndGaps_<True, False, False, True> Type;
};


// 1011
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, false, true, true, TSpec> >
{
    typedef FreeEndGaps_<True, False, True, True> Type;
};


// 1100
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, true, false, false, TSpec> >
{
    typedef FreeEndGaps_<True, True, False, False> Type;
};


// 1101
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, true, false, true, TSpec> >
{
    typedef FreeEndGaps_<True, True, True, False> Type;
};


// 1110
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, true, true, false, TSpec> >
{
    typedef FreeEndGaps_<True, True, False, True> Type;
};

// 1111
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, true, true, true, TSpec> >
{
    typedef FreeEndGaps_<True, True, True, True> Type;
};

// ----------------------------------------------------------------------------
// Metafunction SubstituteAlgoTag_
// ----------------------------------------------------------------------------

// NOTE(rmaerker): Needed to substitute the global alingment algo tags to the correct gap model.
template <typename TTag>
struct SubstituteAlgoTag_
{
    typedef TTag Type;
};

template <>
struct SubstituteAlgoTag_<NeedlemanWunsch>
{
    typedef LinearGaps Type;
};

template <>
struct SubstituteAlgoTag_<Gotoh>
{
    typedef AffineGaps Type;
};

// ----------------------------------------------------------------------------
// SetUpAlignmentProfile
// ----------------------------------------------------------------------------

template <typename TDPType, typename TAlignConfig, typename TTraceSwitch, typename TGapCosts>
struct SetupAlignmentProfile_;

// Profile for Needleman-Wunsch algorithm.
template <typename TFreeEndGaps, typename TGapCosts, typename TTraceSwitch>
struct SetupAlignmentProfile_<DPGlobal, TFreeEndGaps, TGapCosts, TTraceSwitch>
{
    typedef DPProfile_<GlobalAlignment_<TFreeEndGaps>, TGapCosts, TTraceSwitch> Type;
};

// Profile for Smith-Waterman algorithm.
template <typename TFreeEndGaps, typename TGapCosts, typename TTraceSwitch>
struct SetupAlignmentProfile_<DPLocal, TFreeEndGaps, TGapCosts, TTraceSwitch>
{
    typedef DPProfile_<LocalAlignment_<>, TGapCosts, TTraceSwitch> Type;
};

// Profile for Waterman-Eggert algorithm
template <typename TFreeEndGaps, typename TGapCosts, typename TTraceSwitch>
struct SetupAlignmentProfile_<DPLocalEnumerate, TFreeEndGaps, TGapCosts, TTraceSwitch>
{
    typedef DPProfile_<LocalAlignment_<SuboptimalAlignment>, TGapCosts, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _usesAffineGaps()
// ----------------------------------------------------------------------------

template <typename TScoringScheme, typename TSeqH, typename TSeqV>
inline bool
_usesAffineGaps(TScoringScheme const & scoringScheme,
                TSeqH const & seqH,
                TSeqV const & seqV)
{
    typedef typename SequenceEntryForScore<TScoringScheme, TSeqH>::Type TSequenceHEntry;
    typedef typename SequenceEntryForScore<TScoringScheme, TSeqV>::Type TSequenceVEntry;

    TSequenceHEntry seqHEntry = sequenceEntryForScore(scoringScheme, seqH, 0);
    TSequenceVEntry seqVEntry = sequenceEntryForScore(scoringScheme, seqV, 0);

    return (scoreGapExtendHorizontal(scoringScheme, seqHEntry, seqVEntry) !=
            scoreGapOpenHorizontal(scoringScheme, seqHEntry, seqVEntry)) ||
           (scoreGapExtendVertical(scoringScheme, seqHEntry, seqVEntry) !=
               scoreGapOpenVertical(scoringScheme, seqHEntry, seqVEntry));
}

// ----------------------------------------------------------------------------
// Function _setUpAndRunAlignment()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapModel,
          typename TDPScoutStateSpec,
          typename TTrace,
          typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue2, typename TScoreSpec,
          typename TDPType, typename TBand, typename TFreeEndGaps, typename TTraceConfig>
typename Value<Score<TScoreValue2, TScoreSpec> >::Type
_setUpAndRunAlignment(DPContext<TScoreValue, TGapModel> & dpContext,
                      TTrace & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue2, TScoreSpec> const & scoringScheme,
                      AlignConfig2<TDPType, TBand, TFreeEndGaps, TTraceConfig> const & alignConfig)
{
    SEQAN_ASSERT_GEQ(length(seqH), 1u);
    SEQAN_ASSERT_GEQ(length(seqV), 1u);

    typedef typename SetupAlignmentProfile_<TDPType, TFreeEndGaps, TGapModel, TTraceConfig>::Type TDPProfile;
    return _computeAlignment(dpContext, traceSegments, dpScoutState, seqH, seqV, scoringScheme, alignConfig._band,
                             TDPProfile());
}

template <typename TTrace,
          typename TDPScoutStateSpec,
          typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue2, typename TScoreSpec,
          typename TDPType, typename TBand, typename TFreeEndGaps, typename TTraceConfig,
          typename TGapModel>
typename Value<Score<TScoreValue2, TScoreSpec> >::Type
_setUpAndRunAlignment(TTrace & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue2, TScoreSpec> const & scoringScheme,
                      AlignConfig2<TDPType, TBand, TFreeEndGaps, TTraceConfig> const & alignConfig,
                      TGapModel const & /**/)
{
    DPContext<TScoreValue2, TGapModel> dpContext;
    return _setUpAndRunAlignment(dpContext, traceSegments, dpScoutState, seqH, seqV, scoringScheme, alignConfig);
}

template <typename TTrace,
          typename TDPScoutStateSpec,
          typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue2, typename TScoreSpec,
          typename TDPType, typename TBand, typename TFreeEndGaps, typename TTraceConfig>
typename Value<Score<TScoreValue2, TScoreSpec> >::Type
_setUpAndRunAlignment(TTrace & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue2, TScoreSpec> const & scoringScheme,
                      AlignConfig2<TDPType, TBand, TFreeEndGaps, TTraceConfig> const & alignConfig)
{
    if (_usesAffineGaps(scoringScheme, seqH, seqV))
        return _setUpAndRunAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, alignConfig, AffineGaps());
    else
        return _setUpAndRunAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, alignConfig, LinearGaps());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_SETUP_H_
