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

template <typename TSimdVec>
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

template <typename TSimdVector>
class DPScoutState_<SimdAlignVariableLength<TSimdVector> >
{
public:
    String<TSimdVector, Alloc<OverAligned> > masksH;
    String<TSimdVector, Alloc<OverAligned> > masksV;
    String<TSimdVector, Alloc<OverAligned> > masks;

    std::vector<size_t> endsH;
    std::vector<size_t> endsV;

    decltype(endsH.begin()) nextEndsH;
    decltype(endsH.begin()) nextEndsV;

    size_t dimV;
    size_t posH;
    size_t posV;
    bool right;
    bool bottom;
    bool isLocalAlignment;

    // ----------------------------------------------------------------------------
    // Function DPScout_#updateMasksRight()
    // ----------------------------------------------------------------------------

    inline void updateMasksRight()
    {
        for(size_t pos = dimV - 2; pos != MaxValue<size_t>::VALUE; --pos)
            masks[pos] |= masks[pos + 1];
    }

    // ----------------------------------------------------------------------------
    // Function DPScout_#updateMasksBottom()
    // ----------------------------------------------------------------------------

    inline void updateMasksBottom()
    {
        for(auto pos : endsV)
            for(auto it = nextEndsH; it != endsH.end(); ++it)
                masks[pos] |= (masksH[*it] & masksV[pos]);
    }

    // ----------------------------------------------------------------------------
    // Function DPScout_#updateMasks()
    // ----------------------------------------------------------------------------

    inline void updateMasks()
    {
        for(size_t pos = 0; pos < dimV; ++pos)
            masks[pos] = masksH[posH] & masksV[pos];
        //for local alignments the BOTTOM parameter must be checked first
        if(isLocalAlignment)
        {
            updateMasksBottom();
            updateMasksRight();
        }
        else
        {
            if(right && posH == *nextEndsH)
                updateMasksRight();
            if(bottom)
                updateMasksBottom();
        }
    }
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

template<typename TAlignmentAlgorithm, typename TSpec>
struct ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm, DPScoutState_<SimdAlignVariableLength<TSpec> > >
{
    typedef SimdAlignmentScout<SimdAlignVariableLength<TSpec> > Type;
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
#if defined(__AVX2__)
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
// Function _scoutBestScore()
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

template <typename TDPCell, typename TSimdVec,
          typename TTraceMatrixNavigator,
          typename TIsLastColumn,
          typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TSimdVec> > > & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/,
                TIsLastRow const & /**/)
{
    auto cmp = cmpGt(_scoreOfCell(activeCell), _scoreOfCell(dpScout._maxScore));
    cmp &= dpScout.state->masks[dpScout.state->posV];
    _copySimdCell(dpScout, activeCell, cmp);
    _updateHostPositions(dpScout, cmp, createVector<SimdVector<int32_t>::Type>(position(navigator)));
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

template <typename TDPCell, typename TSimdVec>
inline void
_preInitScoutHorizontal(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TSimdVec> > > & scout)
{
    scout.state->nextEndsH = scout.state->endsH.begin();
    scout.state->posH = 0;
}

// ----------------------------------------------------------------------------
// Function _preInitScoutVertical()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSimdVec>
inline void
_preInitScoutVertical(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TSimdVec> > > & scout)
{
    scout.state->updateMasks();
    scout.state->nextEndsV = scout.state->endsV.begin();
    scout.state->posV = 0;
}

// ----------------------------------------------------------------------------
// Function _reachedHorizontalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSimdVec, typename TIter>
inline bool
_reachedHorizontalEndPoint(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TSimdVec> > > & scout,
                           TIter const & hIt)
{
    return *(scout.state->nextEndsH) == position(hIt);
}

// ----------------------------------------------------------------------------
// Function _reachedVerticalEndPoint()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSimdVec, typename TIter>
inline bool
_reachedVerticalEndPoint(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TSimdVec> > > & scout,
                         TIter const & vIt)
{
    return *(scout.state->nextEndsV) == position(vIt);
}

// ----------------------------------------------------------------------------
// Function _nextHorizontalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSimdVec>
inline void
_nextHorizontalEndPos(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TSimdVec> > > & scout)
{
    ++scout.state->nextEndsH;
}

// ----------------------------------------------------------------------------
// Function _nextVerticalEndPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSimdVec>
inline void
_nextVerticalEndPos(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TSimdVec> > > & scout)
{
    ++scout.state->nextEndsV;
}

// ----------------------------------------------------------------------------
// Function _incHorizontalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSimdVec>
inline void
_incHorizontalPos(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TSimdVec> > > & scout)
{
    ++scout.state->posH;
}

// ----------------------------------------------------------------------------
// Function _incVerticalPos()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSimdVec>
inline void
_incVerticalPos(DPScout_<TDPCell, SimdAlignmentScout<SimdAlignVariableLength<TSimdVec> > > & scout)
{
    ++scout.state->posV;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_SIMD_DP_SCOUT_SIMD_H_
