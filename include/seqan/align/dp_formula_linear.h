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
// Defines the methods to compute the score when using linear gap costs.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_

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
// Function _internalComputeScore()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueR,
                      TracebackOff const &)
{
    if (activeCell._score < rightCompare)
        activeCell._score = rightCompare;
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueR,
                      TracebackOff const &)
{
    TScoreValue cmp = cmpGt(rightCompare, activeCell._score);
    activeCell._score = blend(activeCell._score, rightCompare, cmp);
    return TraceBitMap_<TScoreValue>::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &)
{
    if (activeCell._score < rightCompare)
    {
        activeCell._score = rightCompare;
        return rightTrace;
    }
    return leftTrace;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL const & leftTrace,
                      TTraceValueR const & rightTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &)
{
    TScoreValue cmp = cmpGt(rightCompare, activeCell._score);
    activeCell._score = blend(activeCell._score, rightCompare, cmp);
    return blend(leftTrace, rightTrace, cmp);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
    inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> > const &)
{
    if (activeCell._score < rightCompare)
    {
        activeCell._score = rightCompare;
        return rightTrace;
    }
    return (activeCell._score == rightCompare) ? (rightTrace | leftTrace) : (leftTrace);
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >, typename TraceBitMap_<TScoreValue>::Type)
_internalComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL const & leftTrace,
                      TTraceValueR const & rightTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> > const &)
{
    // Check for greater values.
    TScoreValue cmp = cmpGt(activeCell._score, rightCompare);  // cmp greater
    activeCell._score = blend(rightCompare, activeCell._score, cmp);  // activeCell._score

    // Check for equality.
    return blend(blend(rightTrace, leftTrace, cmp), leftTrace | rightTrace, cmpEq(rightCompare, activeCell._score));
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                 [RecursionDirectionAll, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & previousDiagonal,
                DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
                DPCell_<TScoreValue, LinearGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionAll const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;
    activeCell._score = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);

    TScoreValue tmpScore = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);

    TTraceValue tv = _internalComputeScore(activeCell,
                                           tmpScore,
                                           TraceBitMap_<TScoreValue>::DIAGONAL,
                                           TraceBitMap_<TScoreValue>::VERTICAL | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                                           TTracebackConfig());
    tmpScore = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    return _internalComputeScore(activeCell,
                                 tmpScore,
                                 tv,
                                 TraceBitMap_<TScoreValue>::HORIZONTAL | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                                 TTracebackConfig());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore       [RecursionDirectionUpperDiagonal, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & previousDiagonal,
                DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
                DPCell_<TScoreValue, LinearGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionUpperDiagonal const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig> const &)
{
    activeCell._score = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    TScoreValue tmpScore = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    return _internalComputeScore(activeCell,
                                 tmpScore,
                                 TraceBitMap_<TScoreValue>::DIAGONAL,
                                 TraceBitMap_<TScoreValue>::HORIZONTAL | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                                 TTracebackConfig());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore       [RecursionDirectionLowerDiagonal, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & previousDiagonal,
                DPCell_<TScoreValue, LinearGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, LinearGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionLowerDiagonal const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig> const &)
{
    activeCell._score = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    TScoreValue tmpScore = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);

    return _internalComputeScore(activeCell,
                                 tmpScore,
                                 TraceBitMap_<TScoreValue>::DIAGONAL,
                                 TraceBitMap_<TScoreValue>::VERTICAL | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                                 TTracebackConfig());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                      [RecursionDirectionHorizontal]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
                DPCell_<TScoreValue, LinearGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionHorizontal const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig> const &)
{
    activeCell._score = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);

    if (!IsTracebackEnabled_<TTracebackConfig>::VALUE)
        return TraceBitMap_<TScoreValue>::NONE;

    return TraceBitMap_<TScoreValue>::HORIZONTAL
                | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                        [RecursionDirectionVertical]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig>
inline typename TraceBitMap_<TScoreValue>::Type
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, LinearGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, LinearGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionVertical const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig> const &)
{
    activeCell._score = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);

    if (!IsTracebackEnabled_<TTracebackConfig>::VALUE)
        return TraceBitMap_<TScoreValue>::NONE;

    return TraceBitMap_<TScoreValue>::VERTICAL
                | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_
