// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
struct GlobalAlignment_{};

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
struct LocalAlignment_{};

typedef LocalAlignment_<> DPLocal;
typedef LocalAlignment_<SuboptimalAlignment> DPLocalEnumerate;

// ----------------------------------------------------------------------------
// Class TraceBitMap_
// ----------------------------------------------------------------------------

// Defines static const/constexpr tables for the traceback directions.
// We use TraceValue_ as a helper to distinguish between vector and
// scalar version.
// For the vector version we use some compile time hackery to initialize the
// const vector values with the corresponding values generically for
// seqan simd vector types and UME::SIMD vector types.
template <typename TValue, typename TIsSimdVector>
struct TraceValue_
{
    typedef uint8_t Type;
    static constexpr Type NONE = 0u;                         //0000000
    static constexpr Type DIAGONAL = 1u;                     //0000001
    static constexpr Type HORIZONTAL = 2u;                   //0000010
    static constexpr Type VERTICAL = 4u;                     //0000100
    static constexpr Type HORIZONTAL_OPEN = 8u;              //0001000
    static constexpr Type VERTICAL_OPEN = 16u;               //0010000
    static constexpr Type MAX_FROM_HORIZONTAL_MATRIX = 32u;  //0100000
    static constexpr Type MAX_FROM_VERTICAL_MATRIX = 64u;    //1000000
    static constexpr Type NO_VERTICAL_TRACEBACK = ~(VERTICAL | VERTICAL_OPEN);
    static constexpr Type NO_HORIZONTAL_TRACEBACK = ~(HORIZONTAL | HORIZONTAL_OPEN);
};

template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::NONE;
template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::DIAGONAL;
template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::HORIZONTAL;
template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::VERTICAL;
template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::HORIZONTAL_OPEN;
template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::VERTICAL_OPEN;
template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::MAX_FROM_HORIZONTAL_MATRIX;
template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::MAX_FROM_VERTICAL_MATRIX;
template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::NO_VERTICAL_TRACEBACK;
template <typename TValue, typename TIsSimdVector>
constexpr typename TraceValue_<TValue, TIsSimdVector>::Type TraceValue_<TValue, TIsSimdVector>::NO_HORIZONTAL_TRACEBACK;

// Recursion anchor to return the generated tuple with all fill values.
template <typename ...TValues>
constexpr auto _fillTraceValueVector(std::index_sequence<0> const &, std::tuple<TValues...> const & t)
{
    return t;
}

// Helper function to fill the vector with the correct number of values.
// We use the std::tuple (which supports constexpr functions) to make it evaluate at compile time.
template <size_t ...I, typename ...TValues>
constexpr auto _fillTraceValueVector(std::index_sequence<I...> const &, std::tuple<TValues...> const & t)
{
    // Expand the tuple by one and return the next tuple while reducing the number of elements to add.
    return _fillTraceValueVector(std::make_index_sequence<sizeof...(I) - 1>{},
                                 std::tuple_cat(std::make_tuple(std::get<0>(t)), t));
}

// Helper class to used to expand the elements from the returned tuple with the fill values.
// NOTE(rrahn): Might be easier to solve with fold expressions in SeqAn3.
template <typename TVector, typename TIndexSequence>
struct TraceValueVectorBase_;

template <typename TVector, size_t ...I>
struct TraceValueVectorBase_<TVector, std::index_sequence<I...>>
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

template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::NONE =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::NONE)))...};
template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::DIAGONAL =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::DIAGONAL)))...};
template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::HORIZONTAL =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::HORIZONTAL)))...};
template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::VERTICAL =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::VERTICAL)))...};
template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::HORIZONTAL_OPEN =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::HORIZONTAL_OPEN)))...};
template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::VERTICAL_OPEN =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::VERTICAL_OPEN)))...};
template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::MAX_FROM_HORIZONTAL_MATRIX =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::MAX_FROM_HORIZONTAL_MATRIX)))...};
template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::MAX_FROM_VERTICAL_MATRIX =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::MAX_FROM_VERTICAL_MATRIX)))...};
template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::NO_VERTICAL_TRACEBACK =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::NO_HORIZONTAL_TRACEBACK)))...};
template <typename TVector, size_t ...I>
const typename TraceValueVectorBase_<TVector, std::index_sequence<I...>>::Type TraceValueVectorBase_<TVector, std::index_sequence<I...>>::NO_HORIZONTAL_TRACEBACK =
    {std::get<I>(_fillTraceValueVector(std::make_index_sequence<LENGTH<TVector>::VALUE>{}, std::make_tuple(TraceValue_<uint8_t, False>::NO_VERTICAL_TRACEBACK)))...};

// SIMD Vector version.
// Simply delegates to the base class by passing the index_sequence with the corresponding length of the vector.
template <typename TVector>
struct TraceValue_<TVector, True> : public TraceValueVectorBase_<TVector, std::make_index_sequence<LENGTH<TVector>::VALUE>>
{};

// Type alias to choose between scalar and simd version of trace value.
template <typename TValue = uint8_t>
using TraceBitMap_ = TraceValue_<TValue, typename Is<SimdVectorConcept<TValue> >::Type>;

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
template <typename TAlignment, typename TGapCosts, typename TTraceback, typename TExecPolicy = Serial>
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

enum class DPProfileTypeId : uint8_t
{
    ALGORITHM = 0,
    GAP_MODEL = 1,
    TRACE_CONFIG = 2,
    EXEC_POLICY = 3
};

// ----------------------------------------------------------------------------
// Metafunction DPContextSpec
// ----------------------------------------------------------------------------

template <typename TDPProfile, DPProfileTypeId ID>
struct DPProfileType;

template <typename TAlignment, typename TGapCosts, typename TTraceback, typename TExecPolicy>
struct DPProfileType<DPProfile_<TAlignment, TGapCosts, TTraceback, TExecPolicy>, DPProfileTypeId::ALGORITHM>
{
    using Type = TAlignment;
};

template <typename TAlignment, typename TGapCosts, typename TTraceback, typename TExecPolicy>
struct DPProfileType<DPProfile_<TAlignment, TGapCosts, TTraceback, TExecPolicy>, DPProfileTypeId::GAP_MODEL>
{
    using Type = TGapCosts;
};

template <typename TAlignment, typename TGapCosts, typename TTraceback, typename TExecPolicy>
struct DPProfileType<DPProfile_<TAlignment, TGapCosts, TTraceback, TExecPolicy>, DPProfileTypeId::TRACE_CONFIG>
{
    using Type = TTraceback;
};

template <typename TAlignment, typename TGapCosts, typename TTraceback, typename TExecPolicy>
struct DPProfileType<DPProfile_<TAlignment, TGapCosts, TTraceback, TExecPolicy>, DPProfileTypeId::EXEC_POLICY>
{
    using Type = TExecPolicy;
};

// ----------------------------------------------------------------------------
// Metafunction GapTraits
// ----------------------------------------------------------------------------

template <typename T>
struct GapTraits;

template <typename T>
struct GapTraits<T const> :
    GapTraits<T>
{};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag, typename TExecPolicy>
struct GapTraits<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag, TExecPolicy> >
{
    typedef TGapCosts Type;
};

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

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag, typename TExecPolicy>
struct IsGlobalAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag, TExecPolicy> >:
    IsGlobalAlignment_<TAlgoSpec>{};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag, typename TExecPolicy>
struct IsGlobalAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag, TExecPolicy> const>:
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

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag, typename TExecPolicy>
struct IsLocalAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag, TExecPolicy> >:
    IsLocalAlignment_<TAlgoSpec>{};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag, typename TExecPolicy>
struct IsLocalAlignment_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag, TExecPolicy> const>:
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

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag, typename TExecPolicy>
struct IsTracebackEnabled_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag, TExecPolicy> >:
    IsTracebackEnabled_<TTraceFlag>{};

template <typename TAlgoSpec, typename TGapCosts, typename TTraceFlag, typename TExecPolicy>
struct IsTracebackEnabled_<DPProfile_<TAlgoSpec, TGapCosts, TTraceFlag, TExecPolicy> const>:
    IsTracebackEnabled_<TTraceFlag>{};

// ----------------------------------------------------------------------------
// Metafunction IsGapsLeft_
// ----------------------------------------------------------------------------

template <typename TTraceConfig>
struct IsGapsLeft_ : False{};

template <typename TTraceSpec>
struct IsGapsLeft_<TracebackOn<TracebackConfig_<TTraceSpec, GapsLeft > > >
        : True{};

template <typename TAlgorithm, typename TGapSpec, typename TTraceConfig, typename TExecPolicy>
struct IsGapsLeft_<DPProfile_<TAlgorithm, TGapSpec, TTraceConfig, TExecPolicy> >
        : IsGapsLeft_<TTraceConfig>{};

// ----------------------------------------------------------------------------
// Metafunction IsSingleTrace_
// ----------------------------------------------------------------------------

template <typename TTraceConfig>
struct IsSingleTrace_ : False
{};

template <typename TGapsPlacement>
struct IsSingleTrace_<TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > > : True
{};

template <typename TAlgorithm, typename TGapSpec, typename TTraceConfig, typename TExecPolicy>
struct IsSingleTrace_<DPProfile_<TAlgorithm, TGapSpec, TTraceConfig, TExecPolicy> > : IsSingleTrace_<TTraceConfig>
{};

// ----------------------------------------------------------------------------
// Metafunction IsFreeEndGap_
// ----------------------------------------------------------------------------

// Checks if for the current dp profile and a given gap location the algorithm uses free gaps.
template <typename TAlignmentSpec, typename TDPSide>
struct IsFreeEndGap_ :
    False {};

template <typename TAlignmentSpec, typename TGapSpec, typename TTracebackSpec, typename TExecPolicy, typename TDPSide>
struct IsFreeEndGap_<DPProfile_<TAlignmentSpec, TGapSpec, TTracebackSpec, TExecPolicy> const, TDPSide>:
    IsFreeEndGap_<TAlignmentSpec, TDPSide>{};

template <typename TAlignmentSpec, typename TGapSpec, typename TTracebackSpec, typename TExecPolicy, typename TDPSide>
struct IsFreeEndGap_<DPProfile_<TAlignmentSpec, TGapSpec, TTracebackSpec, TExecPolicy>, TDPSide>:
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

namespace impl
{
    template <typename TStream>
    inline TStream &
    printTraceValue(TStream & stream, char traceValue)
    {
        if (traceValue & (TraceBitMap_<char>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<char>::VERTICAL))
        {
            stream << "|";
        }
        if (traceValue & (TraceBitMap_<char>::MAX_FROM_HORIZONTAL_MATRIX | TraceBitMap_<char>::HORIZONTAL))
        {
            stream << "-";
        }
        if (traceValue & TraceBitMap_<char>::DIAGONAL)
        {
            stream << "\\";
        }
        return stream;
    }
}  // namespace impl
}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_PROFILE_H_
