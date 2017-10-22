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
// Implements the affine gap cost functions.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_

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

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _internalComputeScore      [RecursionDirectionDiagonal, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueGap,
                      TracebackOff const &,
                      RecursionDirectionDiagonal const &)
{
    using std::max;
    _scoreOfCell(activeCell) = max(_scoreOfCell(activeCell), rightCompare);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueGap,
                      TracebackOff const &,
                      RecursionDirectionDiagonal const &)
{
    _scoreOfCell(activeCell) = max(_scoreOfCell(activeCell), rightCompare);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueGap gapTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &,
                      RecursionDirectionDiagonal const &)
{
    return (activeCell._score <= rightCompare)
        ? (activeCell._score = rightCompare, TraceBitMap_<TScoreValue>::DIAGONAL | leftTrace)
        : (leftTrace | gapTrace);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL const & leftTrace,
                      TTraceValueGap const & gapTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &,
                      RecursionDirectionDiagonal const &)
{
    auto cmp = cmpGt(activeCell._score, rightCompare);
    activeCell._score = blend(rightCompare, activeCell._score, cmp);
    return blend(TraceBitMap_<TScoreValue>::DIAGONAL | leftTrace,
                 leftTrace | gapTrace,
                 cmp);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueGap gapTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionDiagonal const &)
{
    return (activeCell._score <= rightCompare)
        ? ((activeCell._score == rightCompare)
            ? (leftTrace | TraceBitMap_<TScoreValue>::DIAGONAL | gapTrace)
            : (activeCell._score = rightCompare, TraceBitMap_<TScoreValue>::DIAGONAL | leftTrace))
        : (leftTrace | gapTrace);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL const & leftTrace,
                      TTraceValueGap const & gapTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionDiagonal const &)
{
    auto cmpG = cmpGt(rightCompare, activeCell._score);
    auto cmpE = cmpEq(rightCompare, activeCell._score);
    TScoreValue result = leftTrace | gapTrace;
    activeCell._score = blend(activeCell._score, rightCompare, cmpG);
    result = blend(result, TraceBitMap_<TScoreValue>::DIAGONAL | leftTrace, cmpG);
    return blend(result, leftTrace | TraceBitMap_<TScoreValue>::DIAGONAL | gapTrace, cmpE);
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore    [RecursionDirectionHorizontal, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueR,
                      TracebackOff const &,
                      RecursionDirectionHorizontal const &)
{
    using std::max;
    _horizontalScoreOfCell(activeCell) = max(_horizontalScoreOfCell(activeCell), rightCompare);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueR,
                      TracebackOff const &,
                      RecursionDirectionHorizontal const &)
{
    _horizontalScoreOfCell(activeCell) = max(_horizontalScoreOfCell(activeCell), rightCompare);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR  rightTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &,
                      RecursionDirectionHorizontal const &)
{
    return (activeCell._horizontalScore < rightCompare)
        ? (activeCell._horizontalScore = rightCompare, rightTrace)
        : (leftTrace);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR  rightTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &,
                      RecursionDirectionHorizontal const &)
{
    auto cmp = cmpGt(rightCompare, activeCell._horizontalScore);
    activeCell._horizontalScore = blend(activeCell._horizontalScore, rightCompare, cmp);
    return blend(leftTrace, rightTrace, cmp);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionHorizontal const &)
{
    return (activeCell._horizontalScore <= rightCompare)
        ? ((activeCell._horizontalScore == rightCompare)
            ? (leftTrace | rightTrace)
            : (activeCell._horizontalScore = rightCompare, rightTrace))
        : (leftTrace);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL const & leftTrace,
                      TTraceValueR const & rightTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionHorizontal const &)
{
    auto cmpG = cmpGt(rightCompare, activeCell._horizontalScore);
    auto cmpE = cmpEq(rightCompare, activeCell._horizontalScore);
    activeCell._horizontalScore = blend(activeCell._horizontalScore, rightCompare, cmpG);

    auto result = blend(leftTrace, rightTrace, cmpG);
    return blend(result, leftTrace | rightTrace, cmpE);
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore      [RecursionDirectionVertical, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueR,
                      TracebackOff const &,
                      RecursionDirectionVertical const &)
{
    using std::max;
    _verticalScoreOfCell(activeCell) = max(_verticalScoreOfCell(activeCell), rightCompare);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueR,
                      TracebackOff const &,
                      RecursionDirectionVertical const &)
{
    _verticalScoreOfCell(activeCell) = max(_verticalScoreOfCell(activeCell), rightCompare);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &,
                      RecursionDirectionVertical const &)
{
    return (activeCell._verticalScore < rightCompare)
        ? (activeCell._verticalScore = rightCompare, rightTrace)
        : (leftTrace);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL const & leftTrace,
                      TTraceValueR const & rightTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &,
                      RecursionDirectionVertical const &)
{
    auto cmp = cmpGt(rightCompare, activeCell._verticalScore);
    activeCell._verticalScore = blend(activeCell._verticalScore, rightCompare, cmp);
    return blend(leftTrace, rightTrace, cmp);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionVertical const &)
{
    return (activeCell._verticalScore <= rightCompare)
        ? ((activeCell._verticalScore == rightCompare)
            ? (leftTrace | rightTrace)
            : (activeCell._verticalScore = rightCompare, rightTrace))
        : (leftTrace);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionVertical const &)
{
    auto cmpG = cmpGt(rightCompare, activeCell._verticalScore);
    auto cmpE = cmpEq(rightCompare, activeCell._verticalScore);
    activeCell._verticalScore = blend(activeCell._verticalScore, rightCompare, cmpG);

    TScoreValue result = leftTrace;
    result = blend(result, rightTrace, cmpG);
    return blend(result, leftTrace | rightTrace, cmpE);
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore          [Vertical vs Horizontal, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue>>>,
                            typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOff const &)
{
    using std::max;
    _scoreOfCell(activeCell) = max(_horizontalScoreOfCell(activeCell), _verticalScoreOfCell(activeCell));
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue>>,
                            typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOff const &)
{
    _scoreOfCell(activeCell) = max(_horizontalScoreOfCell(activeCell), _verticalScoreOfCell(activeCell));
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &)
{
    return (activeCell._verticalScore < activeCell._horizontalScore)
        ? (activeCell._score = activeCell._horizontalScore, TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX)
        : (activeCell._score = activeCell._verticalScore, TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX);
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &)
{
    auto cmp = cmpGt(activeCell._horizontalScore, activeCell._verticalScore);
    activeCell._score = blend(activeCell._verticalScore, activeCell._horizontalScore, cmp);

    return blend(TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                 TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                 cmp);
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &)
{
    return (activeCell._horizontalScore <= activeCell._verticalScore)
        ? ((activeCell._horizontalScore == activeCell._verticalScore)
            ? (activeCell._score = activeCell._horizontalScore,
               TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX |
               TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX)
            : (activeCell._score = activeCell._verticalScore, TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX))
        : (activeCell._score = activeCell._horizontalScore, TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX);
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &)
{
    auto cmpG = cmpGt(activeCell._horizontalScore, activeCell._verticalScore);
    auto cmpE = cmpEq(activeCell._horizontalScore, activeCell._verticalScore);
    activeCell._score = blend(activeCell._verticalScore, activeCell._horizontalScore, cmpG);

    return blend(blend(TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                       TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                       cmpG),
                 TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                 cmpE);
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                 [RecursionAllDirection, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & previousDiagonal,
                DPCell_<TScoreValue, AffineGaps> const & previousHorizontal,
                DPCell_<TScoreValue, AffineGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionAll const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig, TExecPolicy> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;

    // Compute horizontal direction.
    _horizontalScoreOfCell(activeCell) = _horizontalScoreOfCell(previousHorizontal) +
                                         scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    TTraceValue tvGap =
        _internalComputeScore(activeCell,
                              static_cast<TScoreValue>(_scoreOfCell(previousHorizontal) +
                                                       scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal)),
                              TraceBitMap_<TScoreValue>::HORIZONTAL,
                              TraceBitMap_<TScoreValue>::HORIZONTAL_OPEN,
                              TTracebackConfig(),
                              RecursionDirectionHorizontal());

    // Compute vertical direction.
    _verticalScoreOfCell(activeCell) = _verticalScoreOfCell(previousVertical) +
                                       scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    tvGap |=
        _internalComputeScore(activeCell,
                              static_cast<TScoreValue>(_scoreOfCell(previousVertical) +
                                                       scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal)),
                              TraceBitMap_<TScoreValue>::VERTICAL,
                              TraceBitMap_<TScoreValue>::VERTICAL_OPEN,
                              TTracebackConfig(),
                              RecursionDirectionVertical());

    // Get max from horiztonal and/or vertical direction and compare with diagonal direction.
    TTraceValue tvMax = _internalComputeScore(activeCell, TTracebackConfig());
    return _internalComputeScore(activeCell,
                                 static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                                          score(scoringScheme, seqHVal, seqVVal)),
                                 tvGap, tvMax, TTracebackConfig(), RecursionDirectionDiagonal());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore       [RecursionUpperDiagonalDirection, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & previousDiagonal,
                DPCell_<TScoreValue, AffineGaps> const & previousHorizontal,
                DPCell_<TScoreValue, AffineGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionUpperDiagonal const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig, TExecPolicy> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;

    // Compute horiztonal direction.
    _horizontalScoreOfCell(activeCell) = _horizontalScoreOfCell(previousHorizontal) +
                                         scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    TTraceValue tv = _internalComputeScore(activeCell,
                                           static_cast<TScoreValue>(_scoreOfCell(previousHorizontal) +
                                                                    scoreGapOpenHorizontal(scoringScheme,
                                                                                           seqHVal,
                                                                                           seqVVal)),
                                           TraceBitMap_<TScoreValue>::HORIZONTAL,
                                           TraceBitMap_<TScoreValue>::HORIZONTAL_OPEN,
                                           TTracebackConfig(),
                                           RecursionDirectionHorizontal());
    // Ignore vertical direction in upper diagonal.
    _verticalScoreOfCell(activeCell) = DPCellDefaultInfinity<DPCell_<TScoreValue, AffineGaps> >::VALUE;
    // Compute diagonal direction and compare with horizontal.
    _scoreOfCell(activeCell) = _horizontalScoreOfCell(activeCell);
    return _internalComputeScore(activeCell,
                                 static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                                          score(scoringScheme, seqHVal, seqVVal)),
                                 tv,
                                 TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                                 TTracebackConfig(),
                                 RecursionDirectionDiagonal());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore       [RecursionDirectionLowerDiagonal, AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & previousDiagonal,
                DPCell_<TScoreValue, AffineGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, AffineGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionLowerDiagonal const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig, TExecPolicy> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;

    // Compute vertical direction.
    _verticalScoreOfCell(activeCell) = _verticalScoreOfCell(previousVertical) +
                                       scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    TTraceValue tv =
        _internalComputeScore(activeCell,
                              static_cast<TScoreValue>(_scoreOfCell(previousVertical) +
                                                       scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal)),
                              TraceBitMap_<TScoreValue>::VERTICAL,
                              TraceBitMap_<TScoreValue>::VERTICAL_OPEN,
                              TTracebackConfig(),
                              RecursionDirectionVertical());
    // Ignore horizontal direction in lower diagonal.
    _horizontalScoreOfCell(activeCell) = DPCellDefaultInfinity<DPCell_<TScoreValue, AffineGaps> >::VALUE;
    // Compute diagonal direction and compare with vertical.
    _scoreOfCell(activeCell) = _verticalScoreOfCell(activeCell);
    return _internalComputeScore(activeCell,
                                 static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                                          score(scoringScheme, seqHVal, seqVVal)),
                                 tv,
                                 TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                                 TTracebackConfig(),
                                 RecursionDirectionDiagonal());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                      [RecursionHorizontalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, AffineGaps> const & previousHorizontal,
                DPCell_<TScoreValue, AffineGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionHorizontal const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig, TExecPolicy> const &)
{
    // Compute horizontal direction.
    _horizontalScoreOfCell(activeCell) = _horizontalScoreOfCell(previousHorizontal) +
                                         scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    auto traceDir = _internalComputeScore(activeCell,
                                static_cast<TScoreValue>(_scoreOfCell(previousHorizontal) +
                                                         scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal)),
                                TraceBitMap_<TScoreValue>::HORIZONTAL,
                                TraceBitMap_<TScoreValue>::HORIZONTAL_OPEN,
                                TTracebackConfig(),
                                RecursionDirectionHorizontal()) | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX;
    // Ignore vertical direction.
    _verticalScoreOfCell(activeCell) = DPCellDefaultInfinity<DPCell_<TScoreValue, AffineGaps> >::VALUE;
    _scoreOfCell(activeCell) = _horizontalScoreOfCell(activeCell);
    return traceDir;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                        [RecursionVerticalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, AffineGaps> & activeCell,
                DPCell_<TScoreValue, AffineGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, AffineGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, AffineGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionVertical const &,
                DPProfile_<TAlgorithm, AffineGaps, TTracebackConfig, TExecPolicy> const &)
{
    // Compute vertical direction.
    _verticalScoreOfCell(activeCell) = _verticalScoreOfCell(previousVertical) +
                                       scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    auto traceDir = _internalComputeScore(activeCell,
                                 static_cast<TScoreValue>(_scoreOfCell(previousVertical) +
                                                          scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal)),
                                 TraceBitMap_<TScoreValue>::VERTICAL,
                                 TraceBitMap_<TScoreValue>::VERTICAL_OPEN,
                                 TTracebackConfig(),
                                 RecursionDirectionVertical()) | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX;
    // Ignore horizontal direction.
    _horizontalScoreOfCell(activeCell) = DPCellDefaultInfinity<DPCell_<TScoreValue, AffineGaps> >::VALUE;
    _scoreOfCell(activeCell) = _verticalScoreOfCell(activeCell);
    return traceDir;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_
