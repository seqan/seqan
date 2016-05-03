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
// This header contains all tags, structures and meta-functions that are
// used to define the meta-profile of an alignment algorithm.
// With the meta-profile the sort of alignment can be selected such as
// a global or a local alignment. It further structures the different
// specializations of global and local alignments or selects the gap cost
// function, or enables or disables the trace-back function.
// ==========================================================================

// TODO(holtgrew): Documentation in this header necessary or internal only?

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_PROFILE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_PROFILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FreeEndGaps_
// ----------------------------------------------------------------------------

// Used to determine which end-gaps are free.
template <typename TInitGapsHorizontal = False, typename TInitGapsVertical = False, typename TEndGapsHorizontal = False,
          typename TEndGapsVertical = False>
struct FreeEndGaps_ {};

// ----------------------------------------------------------------------------
// Class SplitBreakpointAlignment
// ----------------------------------------------------------------------------
// TODO(rmaerker): maybe in a different header
// Used to specify the global alignment for split breakpoint computation.
struct AlignmentSplitBreakpoint_;
typedef Tag<AlignmentSplitBreakpoint_> SplitBreakpointAlignment;

// ----------------------------------------------------------------------------
// Class GlobalAlignment_
// ----------------------------------------------------------------------------

// This is used to select global alignments. The default is the standard global
// dp-algorithm.
//
// Note, all global alignments have to be specialized versions of GlobalAlignment_<>
template <typename TSpec = FreeEndGaps_<> >
struct GlobalAlignment_;

typedef GlobalAlignment_<> DPGlobal;

// ----------------------------------------------------------------------------
// Class SuboptimalAlignment
// ----------------------------------------------------------------------------

// TODO(rmaerker): maybe in a different header
// Used to specify the WatermanEggert algorithm.
struct AlignmentSuboptimal_;
typedef Tag<AlignmentSuboptimal_> SuboptimalAlignment;

// ----------------------------------------------------------------------------
// Class LocalAlignment_
// ----------------------------------------------------------------------------

// This is used to select local alignments. The default is the standard local
// dp-algorithm.
//
// Note, all local alignments have to be specialized versions of LocalAlignment_<>

template <typename TSpec = Default>
struct LocalAlignment_;

typedef LocalAlignment_<> DPLocal;
typedef LocalAlignment_<SuboptimalAlignment> DPLocalEnumerate;

// Use macro expansion to define all possible SIMD initialization types.

template <typename TVector, __uint8 FILL_VALUE, unsigned SIZE>
struct InitSimdTrace_;

#define SEQAN_SIMD_INIT_FILL_VALUE_2_ FILL_VALUE, FILL_VALUE
#define SEQAN_SIMD_INIT_FILL_VALUE_4_ SEQAN_SIMD_INIT_FILL_VALUE_2_, SEQAN_SIMD_INIT_FILL_VALUE_2_
#define SEQAN_SIMD_INIT_FILL_VALUE_8_ SEQAN_SIMD_INIT_FILL_VALUE_4_, SEQAN_SIMD_INIT_FILL_VALUE_4_
#define SEQAN_SIMD_INIT_FILL_VALUE_16_ SEQAN_SIMD_INIT_FILL_VALUE_8_, SEQAN_SIMD_INIT_FILL_VALUE_8_
#define SEQAN_SIMD_INIT_FILL_VALUE_32_ SEQAN_SIMD_INIT_FILL_VALUE_16_, SEQAN_SIMD_INIT_FILL_VALUE_16_

#define SEQAN_SIMD_TRACE_SETUP_2_(SIZE, ...)                                                        \
template <typename TVector, __uint8 FILL_VALUE>                                                     \
struct InitSimdTrace_<TVector, FILL_VALUE, SIZE>                                                    \
{                                                                                                   \
    static const TVector VALUE;                                                                     \
};                                                                                                  \
                                                                                                    \
template <typename TVector, __uint8 FILL_VALUE>                                                     \
const TVector InitSimdTrace_<TVector, FILL_VALUE, SIZE>::VALUE = TVector{__VA_ARGS__};

#define SEQAN_SIMD_TRACE_SETUP_1_(SIZE, MACRO) SEQAN_SIMD_TRACE_SETUP_2_(SIZE, MACRO)
#define SEQAN_SIMD_TRACE_SETUP_(SIZE) SEQAN_SIMD_TRACE_SETUP_1_(SIZE, SEQAN_SIMD_INIT_FILL_VALUE_ ## SIZE ## _)

SEQAN_SIMD_TRACE_SETUP_(2)
SEQAN_SIMD_TRACE_SETUP_(4)
SEQAN_SIMD_TRACE_SETUP_(8)
SEQAN_SIMD_TRACE_SETUP_(16)
SEQAN_SIMD_TRACE_SETUP_(32)

// Scalar version.
template <typename TValue, typename TIsSimdVector>
struct TraceValue_
{
    typedef uint8_t Type;
    static const Type NONE = 0u;                         //0000000
    static const Type DIAGONAL = 1u;                     //0000001
    static const Type HORIZONTAL = 2u;                   //0000010
    static const Type VERTICAL = 4u;                     //0000100
    static const Type HORIZONTAL_OPEN = 8u;              //0001000
    static const Type VERTICAL_OPEN = 16u;               //0010000
    static const Type MAX_FROM_HORIZONTAL_MATRIX = 32u;  //0100000
    static const Type MAX_FROM_VERTICAL_MATRIX = 64u;    //1000000
    static const Type NO_VERTICAL_TRACEBACK = ~(VERTICAL | VERTICAL_OPEN);
    static const Type NO_HORIZONTAL_TRACEBACK = ~(HORIZONTAL | HORIZONTAL_OPEN);
};

// SIMD Vector version.
template <typename TVector>
struct TraceValue_<TVector, True>
{
    typedef TVector Type;
    static const Type NONE;
    static const Type DIAGONAL;
    static const Type HORIZONTAL;
    static const Type VERTICAL;
    static const Type HORIZONTAL_OPEN;
    static const Type VERTICAL_OPEN;
    static const Type MAX_FROM_HORIZONTAL_MATRIX;
    static const Type MAX_FROM_VERTICAL_MATRIX;
    static const Type NO_HORIZONTAL_TRACEBACK;
    static const Type NO_VERTICAL_TRACEBACK;
};

// Macro expansion to define out-of-class initialization of static members.

#define SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(TRACE_VALUE)                             \
    template <typename TVector>                                                      \
    const TVector TraceValue_<TVector, True>::TRACE_VALUE = InitSimdTrace_<TVector, TraceValue_<__uint8, False>::TRACE_VALUE, LENGTH<TVector>::VALUE>::VALUE;

SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(NONE)
SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(DIAGONAL)
SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(HORIZONTAL)
SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(VERTICAL)
SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(HORIZONTAL_OPEN)
SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(VERTICAL_OPEN)
SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(MAX_FROM_HORIZONTAL_MATRIX)
SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(MAX_FROM_VERTICAL_MATRIX)
SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(NO_HORIZONTAL_TRACEBACK)
SEQAN_SIMD_TRACE_OUT_OF_CLASS_INIT_(NO_VERTICAL_TRACEBACK)

// Type alias to choose between scalar and simd version of trace value.
template <typename TValue = __uint8>
using TraceBitMap_ = TraceValue_<TValue, typename Is<SimdVectorConcept<TValue> >::Type >;

// ----------------------------------------------------------------------------
// Tag GapsLeft
// ----------------------------------------------------------------------------

struct GapsLeft_;
typedef Tag<GapsLeft_> GapsLeft;

// ----------------------------------------------------------------------------
// Tag GapsRight
// ----------------------------------------------------------------------------

struct GapsRight_;
typedef Tag<GapsRight_> GapsRight;


// ----------------------------------------------------------------------------
// Tag SingleTrace
// ----------------------------------------------------------------------------

struct SingleTrace_;
typedef Tag<SingleTrace_> SingleTrace;

// ----------------------------------------------------------------------------
// Tag CompleteTrace
// ----------------------------------------------------------------------------

struct CompleteTrace_;
typedef Tag<CompleteTrace_> CompleteTrace;

// ----------------------------------------------------------------------------
// Tag TracebackConfig_
// ----------------------------------------------------------------------------

template <typename TTracesSpec, typename TGapsPlacement>
struct TracebackConfig_ {};

// ----------------------------------------------------------------------------
// Tag TracebackOn
// ----------------------------------------------------------------------------

template <typename TSpec = TracebackConfig_<CompleteTrace, GapsLeft> >
struct TracebackOn {};

// ----------------------------------------------------------------------------
// Tag TracebackOff
// ----------------------------------------------------------------------------

struct TracebackOff_ {};
typedef Tag<TracebackOff_> TracebackOff;

// ----------------------------------------------------------------------------
// Tag LinearGaps
// ----------------------------------------------------------------------------

/*!
 * @tag AlignmentAlgorithmTags#LinearGaps
 * @headerfile <seqan/align.h>
 * @brief Tag for selecting linear gap cost model. This tag can be used for all standard DP algorithms.
 *
 * @signature struct LinearGaps_;
 * @signature typedef Tag<LinearGaps_> LinearGaps;
 */
struct LinearGaps_;
typedef Tag<LinearGaps_> LinearGaps;

// ----------------------------------------------------------------------------
// Tag AffineGaps
// ----------------------------------------------------------------------------

/*!
 * @tag AlignmentAlgorithmTags#AffineGaps
 * @headerfile <seqan/align.h>
 * @brief Tag for selecting affine gap cost model. This tag can be used for all standard DP algorithms.
 *
 * @signature struct AffineGaps_;
 * @signature typedef Tag<AffineGaps_> AffineGaps;
 */
struct AffineGaps_;
typedef Tag<AffineGaps_> AffineGaps;

// ----------------------------------------------------------------------------
// Tag DynamicGaps
// ----------------------------------------------------------------------------

/*!
 * @tag AlignmentAlgorithmTags#DynamicGaps
 * @headerfile <seqan/align.h>
 * @brief Tag for selecting dynamic gap cost model. This tag can be used for all standard DP algorithms.
 *
 * @signature struct DynamicGaps_;
 * @signature typedef Tag<DynamicGaps_> DynamicGaps;
 */

struct DynamicGaps_;
typedef Tag<DynamicGaps_> DynamicGaps;

// ----------------------------------------------------------------------------
// Class DPProfile
// ----------------------------------------------------------------------------

// This meta-object takes three types to be specialized.
//
// TAlignment: The type to select the pairwise alignment algorithm.
// TGapCosts:  The gap cost function (LinearGaps or AffineGaps).
// TTraceback: The traceback switch (TracebackOn or TracebackOff).
template <typename TAlignment, typename TGapCosts, typename TTraceback>
struct DPProfile_ {};


// ----------------------------------------------------------------------------
// Tag DPFirstRow
// ----------------------------------------------------------------------------

// These tags are used to specify the four locations of a dp-matrix where
// free gaps can occur.
struct DPFirstRow_;
typedef Tag<DPFirstRow_> DPFirstRow;

// ----------------------------------------------------------------------------
// Tag DPFirstColumn
// ----------------------------------------------------------------------------

struct DPFirstColumn_;
typedef Tag<DPFirstColumn_> DPFirstColumn;

// ----------------------------------------------------------------------------
// Tag DPLastRow
// ----------------------------------------------------------------------------

struct DPLastRow_;
typedef Tag<DPLastRow_> DPLastRow;

// ----------------------------------------------------------------------------
// Tag DPLastColumn
// ----------------------------------------------------------------------------

struct DPLastColumn_;
typedef Tag<DPLastColumn_> DPLastColumn;

template <typename TDPType, typename TBand, typename TFreeEndGaps = FreeEndGaps_<False, False, False, False>,
          typename TTraceConfig = TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > >
class AlignConfig2
{
public:
    TBand _band;

    AlignConfig2() : _band()
    {}

    template <typename TPosition>
    AlignConfig2(TPosition const & lDiag, TPosition const & uDiag) : _band(lDiag, uDiag)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction IsGlobalAlignment
// ----------------------------------------------------------------------------

// Checks if the dp profile is a global alignment.
template <typename T>
struct IsGlobalAlignment_ :
    False {};

template <typename TSpec>
struct IsGlobalAlignment_<GlobalAlignment_<TSpec> >:
    True {};

template <typename TSpec>
struct IsGlobalAlignment_<GlobalAlignment_<TSpec> const>:
    True {};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag>
struct IsGlobalAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag> >:
    IsGlobalAlignment_<TAlgoSpec>{};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag>
struct IsGlobalAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag> const>:
    IsGlobalAlignment_<TAlgoSpec>{};

// ----------------------------------------------------------------------------
// Metafunction TraceTail_
// ----------------------------------------------------------------------------

// define whether to include the 'tail' of an alignment in the trace
template <typename TSpec>
struct TraceTail_ :
    IsGlobalAlignment_<TSpec>{};

// ----------------------------------------------------------------------------
// Metafunction TraceHead_
// ----------------------------------------------------------------------------

// define whether to include the 'head' of an alignment in the trace
template <typename TSpec>
struct TraceHead_ :
    IsGlobalAlignment_<TSpec>{};

// ----------------------------------------------------------------------------
// Metafunction HasTerminationCriterium_
// ----------------------------------------------------------------------------

// check whether an algorithm has an early termination criterium
// if an algorithm has this, it will get a DPscout that can be terminated
// see dp_scout.h for more info
template <typename TSpec>
struct HasTerminationCriterium_ :
    False {};

// ----------------------------------------------------------------------------
// Metafunction IsLocalAlignment_
// ----------------------------------------------------------------------------

// Checks if the dp profile is a local alignment.
template <typename T>
struct IsLocalAlignment_ :
    False {};

template <typename TSpec>
struct IsLocalAlignment_<LocalAlignment_<TSpec> >:
    True {};

template <typename TSpec>
struct IsLocalAlignment_<LocalAlignment_<TSpec> const>:
    True {};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag>
struct IsLocalAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag> >:
    IsLocalAlignment_<TAlgoSpec>{};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag>
struct IsLocalAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag> const>:
    IsLocalAlignment_<TAlgoSpec>{};

// ----------------------------------------------------------------------------
// Metafunction IsTracebackEnabled_
// ----------------------------------------------------------------------------

// Checks if the trace-back for the current dp profile is enabled.
template <typename T>
struct IsTracebackEnabled_ :
    False {};

template <typename TTracebackConfig>
struct IsTracebackEnabled_<TracebackOn<TTracebackConfig> >:
    True {};

template <typename TTracebackConfig>
struct IsTracebackEnabled_<TracebackOn<TTracebackConfig> const>:
    True {};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag>
struct IsTracebackEnabled_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag> >:
    IsTracebackEnabled_<TTraceFlag>{};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag>
struct IsTracebackEnabled_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag> const>:
    IsTracebackEnabled_<TTraceFlag>{};

// ----------------------------------------------------------------------------
// Metafunction IsGapsLeft_
// ----------------------------------------------------------------------------

template <typename TTraceConfig>
struct IsGapsLeft_ : False{};

template <typename TTraceSpec>
struct IsGapsLeft_<TracebackOn<TracebackConfig_<TTraceSpec, GapsLeft > > >
        : True{};

template <typename TAlgorithm, typename TGapSpec, typename TTraceConfig>
struct IsGapsLeft_<DPProfile_<TAlgorithm, TGapSpec, TTraceConfig> >
        : IsGapsLeft_<TTraceConfig>{};

// ----------------------------------------------------------------------------
// Metafunction IsSingleTrace_
// ----------------------------------------------------------------------------

template <typename TTraceConfig>
struct IsSingleTrace_ : False{};

template <typename TGapsPlacement>
struct IsSingleTrace_<TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > >
: True{};

template <typename TAlgorithm, typename TGapSpec, typename TTraceConfig>
struct IsSingleTrace_<DPProfile_<TAlgorithm, TGapSpec, TTraceConfig> >
: IsSingleTrace_<TTraceConfig>{};

// ----------------------------------------------------------------------------
// Metafunction IsFreeEndGap_
// ----------------------------------------------------------------------------

// Checks if for the current dp profile and a given gap location the algorithm uses free gaps.
template <typename TAlignmentSpec, typename TDPSide>
struct IsFreeEndGap_ :
    False {};

template <typename TAlignmentSpec, typename TGapSpec, typename TTracebackSpec, typename TDPSide>
struct IsFreeEndGap_<DPProfile_<TAlignmentSpec, TGapSpec, TTracebackSpec> const, TDPSide>:
    IsFreeEndGap_<TAlignmentSpec, TDPSide>{};

template <typename TAlignmentSpec, typename TGapSpec, typename TTracebackSpec, typename TDPSide>
struct IsFreeEndGap_<DPProfile_<TAlignmentSpec, TGapSpec, TTracebackSpec>, TDPSide>:
    IsFreeEndGap_<TAlignmentSpec, TDPSide>{};

template <typename TLocalSpec, typename TDPSide>
struct IsFreeEndGap_<LocalAlignment_<TLocalSpec> const, TDPSide>:
    True
{};

template <typename TLocalSpec, typename TDPSide>
struct IsFreeEndGap_<LocalAlignment_<TLocalSpec>, TDPSide>:
    True
{};

template <typename TFreeEndGaps, typename TDPSide>
struct IsFreeEndGap_<GlobalAlignment_<TFreeEndGaps> const, TDPSide>:
    IsFreeEndGap_<TFreeEndGaps, TDPSide>
{};

template <typename TFreeEndGaps, typename TDPSide>
struct IsFreeEndGap_<GlobalAlignment_<TFreeEndGaps>, TDPSide>:
    IsFreeEndGap_<TFreeEndGaps const, TDPSide>
{};

template <typename TFirstColumn, typename TLastRow, typename TLastColumn>
struct IsFreeEndGap_<FreeEndGaps_<True, TFirstColumn, TLastRow, TLastColumn>, DPFirstRow>:
    True
{};

template <typename TFirstColumn, typename TLastRow, typename TLastColumn>
struct IsFreeEndGap_<FreeEndGaps_<True, TFirstColumn, TLastRow, TLastColumn> const, DPFirstRow>:
    True
{};

template <typename TFirstRow, typename TLastRow, typename TLastColumn>
struct IsFreeEndGap_<FreeEndGaps_<TFirstRow, True, TLastRow, TLastColumn>, DPFirstColumn>:
    True
{};

template <typename TFirstRow, typename TLastRow, typename TLastColumn>
struct IsFreeEndGap_<FreeEndGaps_<TFirstRow, True, TLastRow, TLastColumn> const, DPFirstColumn>:
    True
{};

template <typename TFirstRow, typename TFirstColumn, typename TLastColumn>
struct IsFreeEndGap_<FreeEndGaps_<TFirstRow, TFirstColumn, True, TLastColumn>, DPLastRow>:
    True
{};

template <typename TFirstRow, typename TFirstColumn, typename TLastColumn>
struct IsFreeEndGap_<FreeEndGaps_<TFirstRow, TFirstColumn, True, TLastColumn> const, DPLastRow>:
    True
{};

template <typename TFirstRow, typename TFirstColumn, typename TLastRow>
struct IsFreeEndGap_<FreeEndGaps_<TFirstRow, TFirstColumn, TLastRow, True>, DPLastColumn>:
    True
{};

template <typename TFirstRow, typename TFirstColumn, typename TLastRow>
struct IsFreeEndGap_<FreeEndGaps_<TFirstRow, TFirstColumn, TLastRow, True> const, DPLastColumn>:
    True
{};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_PROFILE_H_
