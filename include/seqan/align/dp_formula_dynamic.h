// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
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
// Implements the score function for dynamic gap costs published in
// "Dynamic Gaps Selector: A Smith Waterman Sequence Alignment Algorithm with
// Affine Gap Model Optimization" by Gianvito Urgese et al.
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_DP_FORMULA_DYNAMIC_H_
#define INCLUDE_SEQAN_ALIGN_DP_FORMULA_DYNAMIC_H_

namespace seqan2 {

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
// Function _internalComputeScore      [RecursionDirectionDiagonal, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & diagCompare,
                      TTraceValueL,
                      TTraceValueGap,
                      TracebackOff const &,
                      RecursionDirectionDiagonal const &)
{
    if(_scoreOfCell(activeCell) < diagCompare)
    {
        activeCell._score = diagCompare;
        setGapExtension(activeCell, False(), False());
        return TraceBitMap_<TScoreValue>::NONE;
    }
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & diagCompare,
                      TTraceValueL,
                      TTraceValueGap,
                      TracebackOff const &,
                      RecursionDirectionDiagonal const &)
{
    auto cmp = cmpGt(diagCompare, _scoreOfCell(activeCell));
    activeCell._score = blend(activeCell._score, diagCompare, cmp);
    setGapExtension(activeCell, False(), False(), cmp);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & diagCompare,
                      TTraceValueL leftTrace,
                      TTraceValueGap gapTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &,
                      RecursionDirectionDiagonal const &)
{
    if(_scoreOfCell(activeCell) <= diagCompare)
    {
        activeCell._score = diagCompare;
        setGapExtension(activeCell, False(), False());
        return TraceBitMap_<TScoreValue>::DIAGONAL | leftTrace;
    }
    return leftTrace | gapTrace;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & diagCompare,
                      TTraceValueL const & leftTrace,
                      TTraceValueGap const & gapTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &,
                      RecursionDirectionDiagonal const &)
{
    auto cmp = cmpGt(_scoreOfCell(activeCell), diagCompare);
    activeCell._score = blend(diagCompare, activeCell._score, cmp);
    setGapExtension(activeCell, False(), False(), cmp);
    return blend(TraceBitMap_<TScoreValue>::DIAGONAL | leftTrace,
                 leftTrace | gapTrace,
                 cmp);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & diagCompare,
                      TTraceValueL leftTrace,
                      TTraceValueGap gapTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionDiagonal const &)
{
    if(_scoreOfCell(activeCell) < diagCompare)
    {
        activeCell._score = diagCompare;  // Maximal score comes from diagonal.
        setGapExtension(activeCell, False(), False());
        return TraceBitMap_<TScoreValue>::DIAGONAL | leftTrace;  // Return trace for Diagonal.
    }
    if (_scoreOfCell(activeCell) == diagCompare)  // Maximal score comes from previous computed directions and diagonal.
        return leftTrace | TraceBitMap_<TScoreValue>::DIAGONAL | gapTrace;  // Return all directions inclusively the flag indicating max from gap.

    return leftTrace | gapTrace;  // Maximum comes from gap. Return gap value inclusively the flag indicating max from gap.
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueGap, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & diagCompare,
                      TTraceValueL const & leftTrace,
                      TTraceValueGap const & gapTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &,
                      RecursionDirectionDiagonal const &)
{
    auto cmp = cmpGt(diagCompare, _scoreOfCell(activeCell));
    activeCell._score = blend(activeCell._score, diagCompare, cmp);
    setGapExtension(activeCell, False(), False(), cmp);
    return blend(blend(leftTrace | gapTrace,
                       TraceBitMap_<TScoreValue>::DIAGONAL | leftTrace,
                       cmp),
                 leftTrace | TraceBitMap_<TScoreValue>::DIAGONAL | gapTrace,
                 cmpEq(_scoreOfCell(activeCell), diagCompare));
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore    [RecursionDirectionHorizontal, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOff const &,
                      RecursionDirectionHorizontal const &)
{
    if(!isGapExtension(prevCell, DynamicGapExtensionHorizontal()))
        activeCell._score = _scoreOfCell(prevCell) + scoreGapOpenHorizontal(score, valH, valV);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOff const &,
                      RecursionDirectionHorizontal const &)
{
    activeCell._score = blend(_scoreOfCell(prevCell) + scoreGapOpenHorizontal(score, valH, valV), activeCell._score,
                              isGapExtension(prevCell, DynamicGapExtensionHorizontal()));
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore, typename TTraceConfig>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOn<TTraceConfig> const &,
                      RecursionDirectionHorizontal const &)
{
    if (!isGapExtension(prevCell, DynamicGapExtensionHorizontal()))
    {
        activeCell._score = _scoreOfCell(prevCell) + scoreGapOpenHorizontal(score, valH, valV);
        return TraceBitMap_<TScoreValue>::HORIZONTAL_OPEN;
    }
    return TraceBitMap_<TScoreValue>::HORIZONTAL;
}

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore, typename TTraceConfig>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOn<TTraceConfig> const &,
                      RecursionDirectionHorizontal const &)
{
    activeCell._score = blend(_scoreOfCell(prevCell) + scoreGapOpenHorizontal(score, valH, valV), activeCell._score,
                              isGapExtension(prevCell, DynamicGapExtensionHorizontal()));
    return TraceBitMap_<TScoreValue>::HORIZONTAL;
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore      [RecursionDirectionVertical, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOff const &,
                      RecursionDirectionVertical const &)
{
    if(!isGapExtension(prevCell, DynamicGapExtensionVertical()))
        activeCell._score = _scoreOfCell(prevCell) + scoreGapOpenVertical(score, valH, valV);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOff const &,
                      RecursionDirectionVertical const &)
{
    activeCell._score = blend(_scoreOfCell(prevCell) + scoreGapOpenVertical(score, valH, valV), activeCell._score,
                              isGapExtension(prevCell, DynamicGapExtensionVertical()));
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore, typename TTraceConfig>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOn<TTraceConfig> const &,
                      RecursionDirectionVertical const &)
{
    if (!isGapExtension(prevCell, DynamicGapExtensionVertical()))
    {
        activeCell._score = _scoreOfCell(prevCell) + scoreGapOpenVertical(score, valH, valV);
        return TraceBitMap_<TScoreValue>::VERTICAL_OPEN;
    }
    return TraceBitMap_<TScoreValue>::VERTICAL;
}

template <typename TScoreValue, typename TValueH, typename TValueV, typename TScore, typename TTraceConfig>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      DPCell_<TScoreValue, DynamicGaps> const & prevCell,
                      TValueH const & valH,
                      TValueV const & valV,
                      TScore const & score,
                      TracebackOn<TTraceConfig> const &,
                      RecursionDirectionVertical const &)
{
    auto cmp = isGapExtension(prevCell, DynamicGapExtensionVertical());
    activeCell._score = blend(_scoreOfCell(prevCell) + scoreGapOpenVertical(score, valH, valV), activeCell._score, cmp);
    return blend(TraceBitMap_<TScoreValue>::VERTICAL_OPEN,
                 TraceBitMap_<TScoreValue>::VERTICAL,
                 cmp);
}

// ----------------------------------------------------------------------------
// Function _internalComputeScore          [Vertical vs Horizontal, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & horizontalComp,
                      TracebackOff const &)
{
    if(_scoreOfCell(activeCell) < horizontalComp)
    {
        activeCell._score = horizontalComp;
        setGapExtension(activeCell, False(), True());
        return TraceBitMap_<TScoreValue>::NONE;
    }
    setGapExtension(activeCell, True(), False());
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & horizontalComp,
                      TracebackOff const &)
{
    auto cmp = cmpGt(horizontalComp, _scoreOfCell(activeCell));
    activeCell._score = blend(activeCell._score, horizontalComp, cmp);
    setGapExtension(activeCell, False(), True(), cmp);
    setGapExtension(activeCell, True(), False(), cmpEq(cmp, TraceBitMap_<TScoreValue>::NONE));
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & horizontalComp,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &)
{
    if(_scoreOfCell(activeCell) < horizontalComp)
    {
        activeCell._score = horizontalComp;
        setGapExtension(activeCell, False(), True());
        return TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX;
    }
    setGapExtension(activeCell, True(), False());
    return TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX;
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & horizontalComp,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> >  const &)
{
    auto cmp = cmpGt(horizontalComp, _scoreOfCell(activeCell));
    activeCell._score = blend(activeCell._score, horizontalComp, cmp);
    setGapExtension(activeCell, False(), True(), cmp);
    setGapExtension(activeCell, True(), False(), cmpEq(cmp, TraceBitMap_<TScoreValue>::NONE));
    return blend(TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                 TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                 cmp);
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & horizontalComp,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &)
{
    if(_scoreOfCell(activeCell) < horizontalComp)
    {
        setGapExtension(activeCell, False(), True());
        activeCell._score = horizontalComp;
        return TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX;
    }
    if (_scoreOfCell(activeCell) == horizontalComp)
    {
        setGapExtension(activeCell, True(), True());
        return TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX;
    }
    setGapExtension(activeCell, True(), False());
    return TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX;
}

template <typename TScoreValue, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
                      TScoreValue const & horizontalComp,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> >  const &)
{
    auto cmpG = cmpGt(horizontalComp, _scoreOfCell(activeCell));
    auto cmpE = cmpEq(horizontalComp, _scoreOfCell(activeCell));
    setGapExtension(activeCell, True(), False(), createVector<TScoreValue>(-1));
    setGapExtension(activeCell, False(), True(), cmpG);
    setGapExtension(activeCell, True(), True(), cmpE);

    TScoreValue result = TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX;
    result = blend(result,
                   TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                   cmpG);
    return blend(result,
                 TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX
                    | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                 cmpE);
}

// ----------------------------------------------------------------------------
// Function _computeScore                 [RecursionAllDirection, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
              DPCell_<TScoreValue, DynamicGaps> & previousDiagonal,
              DPCell_<TScoreValue, DynamicGaps> const & previousHorizontal,
              DPCell_<TScoreValue, DynamicGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionAll const &,
              DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig, TExecPolicy> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;
    typedef DPCell_<TScoreValue, DynamicGaps> TCell;

    // Compute intermediate diagonal result.
    TScoreValue intermediate = static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                                        score(scoringScheme, seqHVal, seqVVal));
    // Cache previous Diagonal
    _scoreOfCell(previousDiagonal) = _scoreOfCell(previousHorizontal);
    // Compute best alignment from either horizontal open or extension.
    TCell tmpScore(_scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal));
    TTraceValue tvGap = _internalComputeScore(tmpScore, previousHorizontal, seqHVal, seqVVal, scoringScheme,
                                              TTracebackConfig(), RecursionDirectionHorizontal());

    // Compute best alignment between vertical and vertical open gap.
    activeCell._score = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    tvGap |= _internalComputeScore(activeCell, previousVertical, seqHVal, seqVVal, scoringScheme,
                                   TTracebackConfig(), RecursionDirectionVertical());

    // Finds the maximum between the vertical and the horizontal matrix. Stores the flag for coming from a potential direction.
    TTraceValue tvMax = _internalComputeScore(activeCell, tmpScore._score, TTracebackConfig());  // Stores from where the maximal score comes.
    tmpScore._score = intermediate;
    tvMax = _internalComputeScore(activeCell, tmpScore._score, tvGap, tvMax, TTracebackConfig(), RecursionDirectionDiagonal());

    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        using TMaxScorePolicy = typename TraceConfigForLocalAlignment_<TTracebackConfig>::Type;
        tvMax = _maxScore(_scoreOfCell(activeCell),
                          TraceBitMap_<TScoreValue>::NONE,
                          _scoreOfCell(activeCell),
                          TraceBitMap_<TScoreValue>::NONE,
                          tvMax,
                          TMaxScorePolicy{});
    }
    previousVertical = activeCell;
    return tvMax;
}

// ----------------------------------------------------------------------------
// Function _computeScore       [RecursionUpperDiagonalDirection, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
              DPCell_<TScoreValue, DynamicGaps> & previousDiagonal,
              DPCell_<TScoreValue, DynamicGaps> const & previousHorizontal,
              DPCell_<TScoreValue, DynamicGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionUpperDiagonal const &,
              DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig, TExecPolicy> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;
    // Compute intermediate diagonal result.
    TScoreValue intermediate = static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                                        score(scoringScheme, seqHVal, seqVVal));
    // Cache previous Diagonal
    _scoreOfCell(previousDiagonal) = _scoreOfCell(previousHorizontal);

    // This computes the difference between the horizontal extend and horizontal open.
    activeCell._score = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    TTraceValue tv = _internalComputeScore(activeCell, previousHorizontal, seqHVal, seqVVal, scoringScheme,
                                           TTracebackConfig(), RecursionDirectionHorizontal());

    setGapExtension(activeCell, False(), True());
    tv = _internalComputeScore(activeCell, intermediate, tv, TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                               TTracebackConfig(),RecursionDirectionDiagonal());

    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        using TMaxScorePolicy = typename TraceConfigForLocalAlignment_<TTracebackConfig>::Type;
        tv = _maxScore(_scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TMaxScorePolicy{});
    }
    previousVertical = activeCell;
    return tv;
}

// ----------------------------------------------------------------------------
// Function _computeScore      [RecursionDirectionLowerDiagonal, DynamicGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
              DPCell_<TScoreValue, DynamicGaps> const & previousDiagonal,
              DPCell_<TScoreValue, DynamicGaps> const & /*previousHorizontal*/,
              DPCell_<TScoreValue, DynamicGaps> const & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionLowerDiagonal const &,
              DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig, TExecPolicy> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;

    // This computes the difference between the vertical extend and vertical open.
    activeCell._score = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    TTraceValue tv = _internalComputeScore(activeCell, previousVertical, seqHVal, seqVVal, scoringScheme,
                                           TTracebackConfig(), RecursionDirectionVertical());
    setGapExtension(activeCell, True(), False());
    TScoreValue tmpScore = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    tv = _internalComputeScore(activeCell, tmpScore, tv, TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                                 TTracebackConfig(), RecursionDirectionDiagonal());
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        using TMaxScorePolicy = typename TraceConfigForLocalAlignment_<TTracebackConfig>::Type;
        tv = _maxScore(_scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TMaxScorePolicy{});
    }
    return tv;
}

// ----------------------------------------------------------------------------
// Function _computeScore                      [RecursionHorizontalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
              DPCell_<TScoreValue, DynamicGaps> & previousDiagonal,
              DPCell_<TScoreValue, DynamicGaps> const previousHorizontal,  // NOTE(rrahn): We want the copy here. Don't change!!!
              DPCell_<TScoreValue, DynamicGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionHorizontal const & tag,
              DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig, TExecPolicy> const &)
{
    // Cache previous diagonal value.
    _scoreOfCell(previousDiagonal) = _scoreOfCell(previousHorizontal);
    activeCell._score = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, False(), True());
    auto tv = _internalComputeScore(activeCell, previousHorizontal, seqHVal, seqVVal, scoringScheme,
                                    TTracebackConfig(), tag) | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX;
    previousVertical = activeCell;
    return tv;
}

// NOTE(rrahn): Here we copy the previousCellHorizontal as it might refer to the same value as acticeCell.
// Since we cannot assure here to read the value before we write it, we have to take a copy from it.
template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_computeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
              DPCell_<TScoreValue, DynamicGaps> & previousDiagonal,
              DPCell_<TScoreValue, DynamicGaps> const previousHorizontal, // NOTE(rrahn): Don't change!!!
              DPCell_<TScoreValue, DynamicGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionHorizontal const & tag,
              DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig, TExecPolicy> const &)
{
    // Cache previous diagonal value.
    _scoreOfCell(previousDiagonal) = _scoreOfCell(previousHorizontal);
    activeCell._score = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, False(), True(), createVector<TScoreValue>(-1));
    auto tv = _internalComputeScore(activeCell, previousHorizontal, seqHVal, seqVVal, scoringScheme,
                                 TTracebackConfig(), tag) | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX;
    previousVertical = activeCell;
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        using TMaxScorePolicy = typename TraceConfigForLocalAlignment_<TTracebackConfig>::Type;
        tv = _maxScore(_scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TMaxScorePolicy{});
    }
    return tv;
}

// ----------------------------------------------------------------------------
// Function _computeScore                        [RecursionVerticalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, DynamicGaps> & activeCell,
              DPCell_<TScoreValue, DynamicGaps> & /*previousDiagonal*/,
              DPCell_<TScoreValue, DynamicGaps> const & /*previousHorizontal*/,
              DPCell_<TScoreValue, DynamicGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionVertical const & tag,
              DPProfile_<TAlgorithm, DynamicGaps, TTracebackConfig, TExecPolicy> const &)
{
    activeCell._score = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    setGapExtension(activeCell, True(), False());
    auto tv = _internalComputeScore(activeCell, previousVertical, seqHVal, seqVVal, scoringScheme,
                                 TTracebackConfig(), tag) | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX;
    previousVertical = activeCell;
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        using TMaxScorePolicy = typename TraceConfigForLocalAlignment_<TTracebackConfig>::Type;
        tv = _maxScore(_scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TMaxScorePolicy{});
    }
    return tv;
}

// Independent of gap cost model.
template <typename TScoreValue,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgoTag, typename TTraceFlag, typename TExecPolicy>
inline auto
_computeScore(DPCell_<TScoreValue, DynamicGaps> & current,
              DPCell_<TScoreValue, DynamicGaps> & cacheDiagonal,
              DPCell_<TScoreValue, DynamicGaps> const & /*cacheHorizontal*/,
              DPCell_<TScoreValue, DynamicGaps> & cacheVertical,
              TSequenceHValue const & /*seqHVal*/,
              TSequenceVValue const & /*seqVVal*/,
              TScoringScheme const & /*scoringScheme*/,
              RecursionDirectionZero const &,
              DPProfile_<TAlgoTag, DynamicGaps, TTraceFlag, TExecPolicy> const &)
{
    _scoreOfCell(current) = TraceBitMap_<TScoreValue>::NONE;
    setGapExtension(current, False(), False());
    _scoreOfCell(cacheDiagonal) = _scoreOfCell(current);
    cacheVertical = current;
    return TraceBitMap_<TScoreValue>::NONE;
}

}  // namespace seqan2

#endif  // INCLUDE_SEQAN_ALIGN_DP_FORMULA_DYNAMIC_H_
