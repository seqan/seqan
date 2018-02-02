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
// Function _computeScore                 [RecursionDirectionAll, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, LinearGaps> & current,
              DPCell_<TScoreValue, LinearGaps> & previousDiagonal,
              DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
              DPCell_<TScoreValue, LinearGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionAll const &,
              DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    // Cache next diagonal.
    auto intermediate = static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal));
    previousDiagonal = _scoreOfCell(previousHorizontal);

    auto tv = _maxScore(_scoreOfCell(current),
                        _scoreOfCell(previousVertical)+
                            scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal),
                        _scoreOfCell(previousHorizontal) +
                            scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal),
                        TraceBitMap_<TScoreValue>::VERTICAL | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                        TraceBitMap_<TScoreValue>::HORIZONTAL | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                        TTracebackConfig{});
    // Compute the intermediate value.
    tv = _maxScore(_scoreOfCell(current),
                   intermediate,
                   _scoreOfCell(current),
                   TraceBitMap_<TScoreValue>::DIAGONAL,
                   tv,
                   TTracebackConfig());
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        tv = _maxScore(_scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TTracebackConfig{});
    }
    previousVertical = current;
    return tv;
}

// ----------------------------------------------------------------------------
// Function _computeScore       [RecursionDirectionUpperDiagonal, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, LinearGaps> & current,
              DPCell_<TScoreValue, LinearGaps> & previousDiagonal,
              DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
              DPCell_<TScoreValue, LinearGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionUpperDiagonal const &,
              DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    // Precalculate diagonal direction.
    auto intermediate = static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal));
    // Cache next diagonal value.
    previousDiagonal = _scoreOfCell(previousHorizontal);

    auto tv = _maxScore(_scoreOfCell(current),
                        intermediate,
                        _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal),
                        TraceBitMap_<TScoreValue>::DIAGONAL,
                        TraceBitMap_<TScoreValue>::HORIZONTAL | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                        TTracebackConfig());
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        tv = _maxScore(_scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TTracebackConfig{});
    }
    previousVertical = current;
    return tv;
}

// ----------------------------------------------------------------------------
// Function _computeScore       [RecursionDirectionLowerDiagonal, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, LinearGaps> & current,
              DPCell_<TScoreValue, LinearGaps> const & previousDiagonal,
              DPCell_<TScoreValue, LinearGaps> const & /*previousHorizontal*/,
              DPCell_<TScoreValue, LinearGaps> const & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionLowerDiagonal const &,
              DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    auto tv  = _maxScore(_scoreOfCell(current),
                         static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                              score(scoringScheme, seqHVal, seqVVal)),
                         static_cast<TScoreValue>(_scoreOfCell(previousVertical) +
                                              scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal)),
                         TraceBitMap_<TScoreValue>::DIAGONAL,
                         TraceBitMap_<TScoreValue>::VERTICAL | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                         TTracebackConfig());
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        tv = _maxScore(_scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TTracebackConfig{});
    }
    return tv;
}

// ----------------------------------------------------------------------------
// Function _computeScore                      [RecursionDirectionHorizontal]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
              DPCell_<TScoreValue, LinearGaps> & previousDiagonal,
              DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
              DPCell_<TScoreValue, LinearGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionHorizontal const &,
              DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    // Cache previous diagonal.
    previousDiagonal = previousHorizontal;
    // Compute current value.
    _scoreOfCell(activeCell) = _scoreOfCell(previousHorizontal) +
                               scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    auto tv = TraceBitMap_<TScoreValue>::HORIZONTAL | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX;
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        tv = _maxScore(_scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(activeCell),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TTracebackConfig{});
    }
    // Cache next vertical.
    previousVertical = activeCell;
    return tv;
}

// ----------------------------------------------------------------------------
// Function _computeScore                        [RecursionDirectionVertical]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, LinearGaps> & current,
              DPCell_<TScoreValue, LinearGaps> const & /*previousDiagonal*/,
              DPCell_<TScoreValue, LinearGaps> const & /*previousHorizontal*/,
              DPCell_<TScoreValue, LinearGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionVertical const &,
              DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    _scoreOfCell(current) = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);
    auto tv = TraceBitMap_<TScoreValue>::VERTICAL | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX;
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        tv = _maxScore(_scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TTracebackConfig{});
    }
    // Cache previous vertical.
    previousVertical = current;
    return tv;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_
