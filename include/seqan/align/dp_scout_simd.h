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

template <typename TSpec = Default>
struct SimdAlignmentDefault_ {};

template <typename TSpec = Default>
struct SimdAlignmentVariable_ {};

struct SimdAlignmentScoutDefault_;
typedef Tag<SimdAlignmentScoutDefault_> SimdAlignmentScoutDefault;

struct SimdAlignmentScoutVariable_;
typedef Tag<SimdAlignmentScoutVariable_> SimdAlignmentScoutVariable;

// ----------------------------------------------------------------------------
// Class DPScoutState_
// ----------------------------------------------------------------------------

template <>
class DPScoutState_<SimdAlignmentScoutDefault> : public Nothing
{};

template <>
class DPScoutState_<SimdAlignmentScoutVariable>
{
public:
    String<TSimdAlign> masksH, masksV, masks;
    std::vector<size_t> endsH, endsV;
    decltype(endsH.begin()) nextEndsH, nextEndsV;

    size_t dimV, posH, posV;
    bool RIGHT, BOTTOM, isLocalAlignment;

    DPScoutState_() {}

    inline void updateMasksRight()
    {
        for(size_t pos = dimV-2; pos != MaxValue<size_t>::VALUE; --pos)
            masks[pos] |= masks[pos+1];
    }

    inline void updateMasksBottom()
    {
        for(auto pos: endsV)
            for(auto it = nextEndsH; it != endsH.end(); ++it)
                masks[pos] |= (masksH[*it] & masksV[pos]);
    }

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
            if(RIGHT && posH == *nextEndsH)
                updateMasksRight();
            if(BOTTOM)
                updateMasksBottom();
        }

    }
};

// ----------------------------------------------------------------------------
// Class DPScout_
// ----------------------------------------------------------------------------

template <typename TDPCell>
class DPScout_<TDPCell, SimdAlignmentScoutDefault> : public DPScout_<TDPCell, Default>
{
public:
    //used in the SIMD version to keep track of all host positions
    //SIMD register size divided by 16bit is the amount of alignments
    //so we need two vectors of type 32bit to save the host for all alignments
    SimdVector<int32_t>::Type _maxHostLow; //first half of alignments
    SimdVector<int32_t>::Type _maxHostHigh; //other half

    DPScout_(DPScoutState_<SimdAlignmentScoutDefault> const & /*state*/) {}
};

template <typename TDPCell>
class DPScout_<TDPCell, SimdAlignmentScoutVariable> : public DPScout_<TDPCell, Default>
{
public:
    //as above in DPScout_<TDPCell, SimdAlignmentScoutDefault
    SimdVector<int32_t>::Type _maxHostLow;
    SimdVector<int32_t>::Type _maxHostHigh;

    DPScoutState_<SimdAlignmentScoutVariable> * state;

    DPScout_() : state(nullptr) {}
    DPScout_(DPScoutState_<SimdAlignmentScoutVariable> & state) : state(&state) {}
};

// ----------------------------------------------------------------------------
// Concepts
// ----------------------------------------------------------------------------

SEQAN_CONCEPT(SimdVariableConcept, (T)) {};
template <typename TDPCell>
SEQAN_CONCEPT_IMPL((DPScout_<TDPCell, SimdAlignmentScoutVariable>), (SimdVariableConcept));

SEQAN_CONCEPT(SimdDefaultConcept, (T)) {};
template <typename TDPCell>
SEQAN_CONCEPT_IMPL((DPScout_<TDPCell, SimdAlignmentScoutDefault>), (SimdDefaultConcept));

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForSimdAlignment_
// ----------------------------------------------------------------------------

template<typename TSpec>
struct ScoutSpecForSimdAlignment_ {};

template<>
struct ScoutSpecForSimdAlignment_<DPScoutState_<SimdAlignmentScoutDefault> >
{
    typedef SimdAlignmentScoutDefault Type;
};

template<>
struct ScoutSpecForSimdAlignment_<DPScoutState_<SimdAlignmentScoutDefault> const>
{
    typedef SimdAlignmentScoutDefault Type;
};

template<>
struct ScoutSpecForSimdAlignment_<DPScoutState_<SimdAlignmentScoutVariable> >
{
    typedef SimdAlignmentScoutVariable Type;
};

template<>
struct ScoutSpecForSimdAlignment_<DPScoutState_<SimdAlignmentScoutVariable_> const>
{
    typedef SimdAlignmentScoutVariable Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _copySimdCell()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec, typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<LinearGapCosts<TDPCell> >, void)
_copySimdCell(DPScout_<TDPCell, TSpec> & dpScout,
              TDPCell const & activeCell,
              TScoreValue const & cmp)
{
    dpScout._maxScore._score = blend(dpScout._maxScore._score, activeCell._score, cmp);
}

template <typename TDPCell, typename TSpec, typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<AffineGapCosts<TDPCell> >, void)
_copySimdCell(DPScout_<TDPCell, TSpec> & dpScout,
              TDPCell const & activeCell,
              TScoreValue const & cmp)
{
    dpScout._maxScore._score = blend(dpScout._maxScore._score, activeCell._score, cmp);
    dpScout._maxScore._horizontalScore = blend(dpScout._maxScore._horizontalScore, activeCell._horizontalScore, cmp);
    dpScout._maxScore._verticalScore = blend(dpScout._maxScore._verticalScore, activeCell._verticalScore, cmp);
}

template <typename TDPCell, typename TSpec, typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<DynamicGapCosts<TDPCell> >, void)
_copySimdCell(DPScout_<TDPCell, TSpec> & dpScout,
              TDPCell const & activeCell,
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
_scoutBestScore(DPScout_<TDPCell, SimdAlignmentScoutDefault> & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/,
                TIsLastRow const & /**/)
{
    TSimdAlign cmp = cmpGt(_scoreOfCell(activeCell), _scoreOfCell(dpScout._maxScore));
    _copySimdCell(dpScout, activeCell, cmp);
    _updateHostPositions(dpScout, cmp, createVector<SimdVector<int32_t>::Type>(position(navigator)));
}

template <typename TDPCell, typename TTraceMatrixNavigator, typename TIsLastColumn, typename TIsLastRow>
inline void
_scoutBestScore(DPScout_<TDPCell, SimdAlignmentScoutVariable> & dpScout,
                TDPCell const & activeCell,
                TTraceMatrixNavigator const & navigator,
                TIsLastColumn const & /**/,
                TIsLastRow const & /**/)
{
    TSimdAlign cmp = cmpGt(_scoreOfCell(activeCell), _scoreOfCell(dpScout._maxScore));
    cmp &= dpScout.state->masks[dpScout.state->posV];
    _copySimdCell(dpScout, activeCell, cmp);
    _updateHostPositions(dpScout, cmp, createVector<SimdVector<int32_t>::Type>(position(navigator)));
}

// ----------------------------------------------------------------------------
// Function maxHostPosition()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TScoutSpec>
inline unsigned int
maxHostPosition(DPScout_<TDPCell, TScoutSpec> const & dpScout, size_t pos)
{
    if(pos < LENGTH<SimdVector<int32_t>::Type>::VALUE)
        return value(dpScout._maxHostLow, pos);
    else
        return value(dpScout._maxHostHigh, pos - LENGTH<SimdVector<int32_t>::Type>::VALUE);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_SIMD_DP_SCOUT_SIMD_H_
