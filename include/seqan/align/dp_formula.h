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
// Defines the recursion formula for the dp-alignment algorithms.
// ==========================================================================

// TODO(holtgrew): Documentation in this header necessary or internal only?

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag RecursionDirectionDiagonal
// ----------------------------------------------------------------------------

struct RecursionDirectionDiagonal_;
typedef Tag<RecursionDirectionDiagonal_> RecursionDirectionDiagonal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionHorizontal
// ----------------------------------------------------------------------------

struct RecursionDirectionHorizontal_;
typedef Tag<RecursionDirectionHorizontal_> RecursionDirectionHorizontal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionVertical
// ----------------------------------------------------------------------------

struct RecursionDirectionVertical_;
typedef Tag<RecursionDirectionVertical_> RecursionDirectionVertical;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionAll
// ----------------------------------------------------------------------------

struct RecursionDirectionAll_;
typedef Tag<RecursionDirectionAll_> RecursionDirectionAll;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionUpperDiagonal
// ----------------------------------------------------------------------------

struct RecursionDirectionUpperDiagonal_;
typedef Tag<RecursionDirectionUpperDiagonal_> RecursionDirectionUpperDiagonal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionLowerDiagonal
// ----------------------------------------------------------------------------

struct RecursionDirectionLowerDiagonal_;
typedef Tag<RecursionDirectionLowerDiagonal_> RecursionDirectionLowerDiagonal;

// ----------------------------------------------------------------------------
// Tag RecursionDirectionZero
// ----------------------------------------------------------------------------

struct RecursionDirectionZero_;
typedef Tag<RecursionDirectionZero_> RecursionDirectionZero;

// ============================================================================
// Metafunctions
// ============================================================================

// Helper typedef to get the correct score value type from the score-matrix navigator.
template <typename TCellTuple>
using ExtractedScoreValueType_ = std::decay_t<decltype(_scoreOfCell(std::get<0>(std::declval<TCellTuple>())))>;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _computeScore
// ----------------------------------------------------------------------------

template <typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TRecursionDirection,
          typename TDPProfile>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<ExtractedScoreValueType_<TRecursionCellTuple>>>>,
                            typename TraceBitMap_<ExtractedScoreValueType_<TRecursionCellTuple>>::Type)
_computeScore(TRecursionCellTuple && recursionCells,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              TRecursionDirection const & recDir,
              TDPProfile const & dpProfile)
{
    auto traceDir = _doComputeScore(std::get<0>(recursionCells),
                                    std::get<1>(recursionCells),
                                    std::get<2>(recursionCells),
                                    std::get<3>(recursionCells),
                                    seqHVal, seqVVal, scoringScheme, recDir, dpProfile);
    if (IsLocalAlignment_<TDPProfile>::VALUE)
    {
        if (_scoreOfCell(std::get<0>(recursionCells)) <= 0)
        {
            _setScoreOfCell(std::get<0>(recursionCells), static_cast<ExtractedScoreValueType_<TRecursionCellTuple>>(0));
            return TraceBitMap_<ExtractedScoreValueType_<TRecursionCellTuple>>::NONE;
        }
    }
    return traceDir;
}

template <typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TRecursionDirection,
          typename TDPProfile>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<ExtractedScoreValueType_<TRecursionCellTuple>>>,
                            typename TraceBitMap_<ExtractedScoreValueType_<TRecursionCellTuple>>::Type)
_computeScore(TRecursionCellTuple && recursionCells,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              TRecursionDirection const & recDir,
              TDPProfile const & dpProfile)
{
    using TScoreValue = ExtractedScoreValueType_<TRecursionCellTuple>;

    auto traceDir = _doComputeScore(std::get<0>(recursionCells),
                                    std::get<1>(recursionCells),
                                    std::get<2>(recursionCells),
                                    std::get<3>(recursionCells),
                                    seqHVal, seqVVal, scoringScheme, recDir, dpProfile);
    if (IsLocalAlignment_<TDPProfile>::VALUE)
    {
        TScoreValue cmp = cmpGt(createVector<TScoreValue>(1), _scoreOfCell(std::get<0>(recursionCells)));
        _setScoreOfCell(std::get<0>(recursionCells), TraceBitMap_<TScoreValue>::NONE, cmp);
        return blend(traceDir, TraceBitMap_<TScoreValue>::NONE, cmp);
    }
    return traceDir;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                        [RecursionDirectionDiagonal]
// ----------------------------------------------------------------------------
// Independent of gap cost model.
template <typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TDPProfile>
inline auto
_doComputeScore(DPCell_<TScoreValue, TGapCosts> & activeCell,
                DPCell_<TScoreValue, TGapCosts> const & previousDiagonal,
                DPCell_<TScoreValue, TGapCosts> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, TGapCosts> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionDiagonal const &,
                TDPProfile const &)
{
    _scoreOfCell(activeCell) = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, False(), False(), createVector<TScoreValue>(-1));

    if (!IsTracebackEnabled_<TDPProfile>::VALUE)
        return TraceBitMap_<TScoreValue>::NONE;
    return TraceBitMap_<TScoreValue>::DIAGONAL;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                            [RecursionDirectionZero]
// ----------------------------------------------------------------------------
// Independent of gap cost model.
template <typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgoTag, typename TTraceFlag>
inline auto
_doComputeScore(DPCell_<TScoreValue, TGapCosts> & activeCell,
                DPCell_<TScoreValue, TGapCosts> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, TGapCosts> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, TGapCosts> const & /*previousVertical*/,
                TSequenceHValue const & /*seqHVal*/,
                TSequenceVValue const & /*seqVVal*/,
                TScoringScheme const & /*scoringScheme*/,
                RecursionDirectionZero const &,
                DPProfile_<TAlgoTag, TGapCosts, TTraceFlag> const &)
{
    _scoreOfCell(activeCell) = createVector<TScoreValue>(0);
    return TraceBitMap_<TScoreValue>::NONE;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_
