// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Stefan Budach <budach@molgen.mpg.de>
// ==========================================================================
// DPScout_ specialization for the SIMD alignment implementation.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_SIMD_DP_SCOUT_SIMD_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_SIMD_DP_SCOUT_SIMD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct SimdAlignEqualLength_;
typedef Tag<SimdAlignEqualLength_> SimdAlignEqualLength;

template <typename TTraits>
struct SimdAlignVariableLength
{};

template <typename TSpec = SimdAlignEqualLength>
struct SimdAlignmentScout {};

// ----------------------------------------------------------------------------
// Class DPScoutState_
// ----------------------------------------------------------------------------

template <>
class DPScoutState_<SimdAlignEqualLength>
{};

template <typename TTraits>
class DPScoutState_<SimdAlignVariableLength<TTraits> >
{
public:
    using TSizeH = typename Size<typename TTraits::TSeqH>::Type;
    using TSizeV = typename Size<typename TTraits::TSeqV>::Type;

    using TIterator = typename Iterator<String<size_t>, Rooted>::Type;

    typename TTraits::TSimdVector endPosVecH;
    typename TTraits::TSimdVector endPosVecV;

    size_t posH;
    size_t posV;

    String<size_t> lengthsH;
    String<size_t> lengthsV;
    String<size_t> endsH;
    String<size_t> endsV;
    TModString     sortedEndsH;
    TModString     sortedEndsV;

    TIterator nextEndsH;
    TIterator nextEndsV;
};

// ----------------------------------------------------------------------------
// Class DPScout_
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
class DPScout_<TDPCell, SimdAlignmentScout<TSpec> > :
    public DPScout_<TDPCell, Default>
{
public:
    using TBase       = DPScout_<TDPCell, Default>;
    using TScoutState = DPScoutState_<TSpec>;

    //used in the SIMD version to keep track of all host positions
    //SIMD register size divided by 16bit is the amount of alignments
    //so we need two vectors of type 32bit to save the host for all alignments

    // TODO(rrahn): Abstract into a struct, so we can model different configurations.
    SimdVector<int32_t>::Type _maxHostLow; //first half of alignments
    SimdVector<int32_t>::Type _maxHostHigh; //other half
    TScoutState * state = nullptr;
    unsigned _simdLane  = 0;

    DPScout_(TScoutState & pState) : TBase(), state(&pState)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForSimdAlignment_
// ----------------------------------------------------------------------------

template<typename TAlignmentAlgorithm>
struct ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm, DPScoutState_<SimdAlignEqualLength> >
{
    typedef SimdAlignmentScout<SimdAlignEqualLength> Type;
};

template<typename TAlignmentAlgorithm, typename TTraits>
struct ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm, DPScoutState_<SimdAlignVariableLength<TTraits> > >
{
    typedef SimdAlignmentScout<SimdAlignVariableLength<TTraits> > Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _copySimdCell()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TScoreValue>
inline void
_copySimdCell(DPScout_<DPCell_<TValue, LinearGaps>, SimdAlignmentScout<TSpec> > & dpScout,
              DPCell_<TValue, LinearGaps> const & activeCell,
              TScoreValue const & cmp)
{
    dpScout._maxScore._score = blend(dpScout._maxScore._score, activeCell._score, cmp);
}

template <typename TValue, typename TSpec, typename TScoreValue>
inline void
_copySimdCell(DPScout_<DPCell_<TValue, AffineGaps>, SimdAlignmentScout<TSpec> > & dpScout,
              DPCell_<TValue, AffineGaps> const & activeCell,
              TScoreValue const & cmp)
{
    dpScout._maxScore._score = blend(dpScout._maxScore._score, activeCell._score, cmp);
    dpScout._maxScore._horizontalScore = blend(dpScout._maxScore._horizontalScore, activeCell._horizontalScore, cmp);
    dpScout._maxScore._verticalScore = blend(dpScout._maxScore._verticalScore, activeCell._verticalScore, cmp);
}

template <typename TValue, typename TSpec, typename TScoreValue>
inline void
_copySimdCell(DPScout_<DPCell_<TValue, DynamicGaps>, SimdAlignmentScout<TSpec> > & dpScout,
              DPCell_<TValue, DynamicGaps> const & activeCell,
              TScoreValue const & cmp)
{
    dpScout._maxScore._score = blend(dpScout._maxScore._score, activeCell._score, cmp);
    dpScout._maxScore._flagMask = blend(dpScout._maxScore._flagMask, activeCell._flagMask, cmp);
}

// ----------------------------------------------------------------------------
// Function _updateHostPositions()
// ----------------------------------------------------------------------------

template<typename TDPCell, typename TScoutSpec, typename TSimdVec>
inline void
_updateHostPositions(DPScout_<TDPCell, TScoutSpec> & dpScout,
                     TSimdVec & cmp,
                     SimdVector<int32_t>::Type positionNavigator)
{
// TODO(rrahn): Refactor!
#if SEQAN_UMESIMD_ENABLED
    using TSimdHalfVec = typename UME::SIMD::SIMDTraits<TSimdVec>::HALF_LEN_VEC_T;
    TSimdHalfVec cmpLow, cmpHigh;
    cmp.unpack(cmpLow, cmpHigh);

    dpScout._maxHostLow = blend(dpScout._maxHostLow, positionNavigator,
                                static_cast<SimdVector<int32_t>::Type>(cmpLow));
    dpScout._maxHostHigh = blend(dpScout._maxHostHigh, positionNavigator,
                                 static_cast<SimdVector<int32_t>::Type>(cmpHigh));
#elif defined(__AVX2__)
    dpScout._maxHostLow = blend(dpScout._maxHostLow, positionNavigator,
                                _mm256_cvtepi16_epi32(_mm256_castsi256_si128(reinterpret_cast<__m256i&>(cmp))));
    dpScout._maxHostHigh = blend(dpScout._maxHostHigh, positionNavigator,
                                 _mm256_cvtepi16_epi32(_mm256_extractf128_si256(reinterpret_cast<__m256i&>(cmp),1)));
#elif defined(__SSE3__)
    dpScout._maxHostLow = blend(dpScout._maxHostLow, positionNavigator,
                                _mm_unpacklo_epi16(reinterpret_cast<__m128i&>(cmp), reinterpret_cast<__m128i&>(cmp)));
    dpScout._maxHostHigh = blend(dpScout._maxHostHigh, positionNavigator,
                                 _mm_unpackhi_epi16(reinterpret_cast<__m128i&>(cmp), reinterpret_cast<__m128i&>(cmp)));
#endif
}

// ----------------------------------------------------------------------------
// Function _scoutBestScore()                            [SimdAlignEqualLength]
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraceMatrixNavigator, typename TIsLastColumn, typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignEqualLength> > & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/,
                TIsLastRow const & /**/)
{
    auto cmp = cmpGt(_scoreOfCell(activeCell), _scoreOfCell(dpScout._maxScore));
    _copySimdCell(dpScout, activeCell, cmp);
    _updateHostPositions(dpScout, cmp, createVector<SimdVector<int32_t>::Type>(position(navigator)));
}

// ----------------------------------------------------------------------------
// Function _getCompareMask()
// ----------------------------------------------------------------------------

// Helper functions to resolve the correct tracking of cells for different
// alignment modes.

// Standard global alignment.
template <typename TDPCell, typename TTraits,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & dpScout,
                            True const & /*lastCol*/,
                            True const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, False, False>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return ((dpScout.state->endPosVecH) == createVector<TSimdVec>(dpScout.state->posH)) &
           ((dpScout.state->endPosVecV) == createVector<TSimdVec>(dpScout.state->posV));
}

template <typename TDPCell, typename TTraits,
          typename TIsLastColumn,
          typename TIsLastRow,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & /*dpScout*/,
                            TIsLastColumn const & /*lastCol*/,
                            TIsLastRow const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, False, False>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return createVector<TSimdVec>(0);
}

// Tracking the last row is enabled
template <typename TDPCell, typename TTraits,
          typename TIsLastColumn,
          typename TIsLastRow,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & dpScout,
                            TIsLastColumn const & /*lastCol*/,
                            TIsLastRow const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, True, False>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return (createVector<TSimdVec>(dpScout.state->posH) <= (dpScout.state->endPosVecH)) &
           (createVector<TSimdVec>(dpScout.state->posV) == (dpScout.state->endPosVecV));
}

template <typename TDPCell, typename TTraits,
          typename TIsLastColumn,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & /*dpScout*/,
                            TIsLastColumn const & /*lastCol*/,
                            False const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, True, False>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return createVector<TSimdVec>(0);
}

// Tracking if the last column is enabled
template <typename TDPCell, typename TTraits,
          typename TIsLastColumn,
          typename TIsLastRow,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & dpScout,
                            TIsLastColumn const & /*lastCol*/,
                            TIsLastRow const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, False, True>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return (createVector<TSimdVec>(dpScout.state->posH) == (dpScout.state->endPosVecH)) &
           (createVector<TSimdVec>(dpScout.state->posV) <= (dpScout.state->endPosVecV));
}

template <typename TDPCell, typename TTraits,
          typename TIsLastRow,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & /*dpScout*/,
                            False const & /*lastCol*/,
                            TIsLastRow const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, False, True>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return createVector<TSimdVec>(0);
}

// Tracking if the last column and last row is enabled
template <typename TDPCell, typename TTraits,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & dpScout,
                            True const & /*lastCol*/,
                            True const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, True, True>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return ((createVector<TSimdVec>(dpScout.state->posH) == (dpScout.state->endPosVecH)) &
            (createVector<TSimdVec>(dpScout.state->posV) <= (dpScout.state->endPosVecV))) |
           ((createVector<TSimdVec>(dpScout.state->posH) <= (dpScout.state->endPosVecH)) &
            (createVector<TSimdVec>(dpScout.state->posV) == (dpScout.state->endPosVecV)));
}

template <typename TDPCell, typename TTraits,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & dpScout,
                            False const & /*lastCol*/,
                            True const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, True, True>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return (createVector<TSimdVec>(dpScout.state->posH) <= (dpScout.state->endPosVecH)) &
           (createVector<TSimdVec>(dpScout.state->posV) == (dpScout.state->endPosVecV));
}

template <typename TDPCell, typename TTraits,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & dpScout,
                            True const & /*lastCol*/,
                            False const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, True, True>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return (createVector<TSimdVec>(dpScout.state->posH) == (dpScout.state->endPosVecH)) &
           (createVector<TSimdVec>(dpScout.state->posV) <= (dpScout.state->endPosVecV));
}

template <typename TDPCell, typename TTraits,
          typename TTop, typename TLeft, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & /*dpScout*/,
                            False const & /*lastCol*/,
                            False const & /*lastRow*/,
                            DPProfile_<GlobalAlignment_<FreeEndGaps_<TTop, TLeft, True, True>>,
                                       TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return createVector<TSimdVec>(0);
}

// If local alignment.
template <typename TDPCell, typename TTraits,
          typename TIsLastColumn,
          typename TIsLastRow,
          typename TAlgoSpec, typename TGapModel, typename TTraceConfig>
inline auto _getCompareMask(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & dpScout,
                            TIsLastColumn const & /*lastCol*/,
                            TIsLastRow const & /*lastRow*/,
                            DPProfile_<LocalAlignment_<TAlgoSpec>, TGapModel, TTraceConfig> const &)
{
    using TSimdVec = typename TTraits::TSimdVector;
    return (createVector<TSimdVec>(dpScout.state->posH) <= (dpScout.state->endPosVecH)) &
           (createVector<TSimdVec>(dpScout.state->posV) <= (dpScout.state->endPosVecV));
}

// ----------------------------------------------------------------------------
// Function _scoutBestScore()                         [SimdAlignVariableLength]
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraits,
          typename TTraceMatrixNavigator,
          typename TIsLastColumn,
          typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/,
                TIsLastRow const & /**/)
{
    auto mask = cmpGt(_scoreOfCell(activeCell), _scoreOfCell(dpScout._maxScore)) &
                _getCompareMask(dpScout, TIsLastColumn{}, TIsLastRow{}, typename TTraits::TDPProfile{});
    _copySimdCell(dpScout, activeCell, mask);
    _updateHostPositions(dpScout, mask, createVector<SimdVector<int32_t>::Type>(position(navigator)));
}

// ----------------------------------------------------------------------------
// Function maxHostPosition()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TScoutSpec>
inline unsigned int
maxHostPosition(DPScout_<TDPCell, SimdAlignmentScout<TScoutSpec> > const & dpScout)
{
    if(dpScout._simdLane < LENGTH<SimdVector<int32_t>::Type>::VALUE)
        return value(dpScout._maxHostLow, dpScout._simdLane);
    else
        return value(dpScout._maxHostHigh, dpScout._simdLane - LENGTH<SimdVector<int32_t>::Type>::VALUE);
}

// ----------------------------------------------------------------------------
// Function _setSimdLane()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TScoutSpec, typename TPosition>
inline void
_setSimdLane(DPScout_<TDPCell, TScoutSpec> & dpScout, TPosition const pos)
{
    dpScout._simdLane = pos;
}

// ----------------------------------------------------------------------------
// Function _preInitScoutHorizontal()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraits>
inline void
_preInitScoutHorizontal(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & scout)
{
    goBegin(scout.state->nextEndsH);
    scout.state->posH = 0;
}

// ----------------------------------------------------------------------------
// Function _preInitScoutVertical()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraits>
inline void
_preInitScoutVertical(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & scout)
{
    // scout.state->updateMasks();
    goBegin(scout.state->nextEndsV);
    scout.state->posV = 0;
}

// ----------------------------------------------------------------------------
// Function _reachedHorizontalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraits, typename TIter>
inline bool
_reachedHorizontalEndPoint(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & scout,
                           TIter const & hIt)
{
    return *(scout.state->nextEndsH) == position(hIt) + 1;
}

// ----------------------------------------------------------------------------
// Function _reachedVerticalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraits, typename TIter>
inline bool
_reachedVerticalEndPoint(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & scout,
                         TIter const & vIt)
{
    return *(scout.state->nextEndsV) == position(vIt) + 1;
}

// ----------------------------------------------------------------------------
// Function _nextHorizontalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraits>
inline void
_nextHorizontalEndPos(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & scout)
{
    ++scout.state->nextEndsH;
}

// ----------------------------------------------------------------------------
// Function _nextVerticalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraits>
inline void
_nextVerticalEndPos(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & scout)
{
    ++scout.state->nextEndsV;
}

// ----------------------------------------------------------------------------
// Function _incHorizontalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraits>
inline void
_incHorizontalPos(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & scout)
{
    ++scout.state->posH;
}

// ----------------------------------------------------------------------------
// Function _incVerticalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TTraits>
inline void
_incVerticalPos(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > & scout)
{
    ++scout.state->posV;
}

// ----------------------------------------------------------------------------
// Function _hostLengthH()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSeqH>
inline auto
_hostLengthH(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignEqualLength> > const & /*scout*/,
             TSeqH const & seqH)
{
    return length(seqH);
}

template <typename TDPCell, typename TTraits, typename TSeqH>
inline auto
_hostLengthH(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & scout,
             TSeqH const & /*seqH*/)
{
    return (scout.state->endPosVecH)[scout._simdLane];
}

// ----------------------------------------------------------------------------
// Function _hostLengthV()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSeqV>
inline auto
_hostLengthV(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignEqualLength> > const & /*scout*/,
             TSeqV const & seqV)
{
    return length(seqV);
}

template <typename TDPCell, typename TTraits, typename TSeqV>
inline auto
_hostLengthV(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TTraits> > > const & scout,
             TSeqV const & /*seqV*/)
{
    return (scout.state->endPosVecV)[scout._simdLane];
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_SIMD_DP_SCOUT_SIMD_H_
