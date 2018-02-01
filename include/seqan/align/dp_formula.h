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
// Defines the recursion formula for the dp-alignment algorithms.
// ==========================================================================

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

template <typename TTarget,
          typename TSourceLeft,
          typename TSourceRight,
          typename TTraceLeft,
          typename TTraceRight>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TTarget>>>, typename TraceBitMap_<TTarget>::Type)
_maxScore(TTarget & target,
          TSourceLeft const & srcLeft,
          TSourceRight const & srcRight,
          TTraceRight const /**/,
          TTraceLeft const  /**/,
          TracebackOff const & /*tag*/)
{
    using std::max;
    target = max(static_cast<TTarget>(srcLeft), static_cast<TTarget>(srcRight));
    return TraceBitMap_<TTarget>::NONE;
}

template <typename TTarget,
          typename TSourceLeft,
          typename TSourceRight,
          typename TTraceLeft,
          typename TTraceRight>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TTarget>>, typename TraceBitMap_<TTarget>::Type)
_maxScore(TTarget & target,
          TSourceLeft const & srcLeft,
          TSourceRight const & srcRight,
          TTraceRight const /**/,
          TTraceLeft const  /**/,
          TracebackOff const & /*tag*/)
{
    target = max(srcLeft, srcRight);
    return TraceBitMap_<TTarget>::NONE;
}

template <typename TTarget,
          typename TSourceLeft,
          typename TSourceRight,
          typename TTraceLeft,
          typename TTraceRight,
          typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TTarget>>>, typename TraceBitMap_<TTarget>::Type)
_maxScore(TTarget & target,
          TSourceLeft const & srcLeft,
          TSourceRight const & srcRight,
          TTraceLeft const  traceLeft,
          TTraceRight const traceRight,
          TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const & /*tag*/)
{
    return (srcLeft < srcRight)
        ? (target = srcRight, traceRight)
        : (target = srcLeft, traceLeft);
}

template <typename TTarget,
          typename TSourceLeft,
          typename TSourceRight,
          typename TTraceLeft,
          typename TTraceRight,
          typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TTarget>>, typename TraceBitMap_<TTarget>::Type)
_maxScore(TTarget & target,
          TSourceLeft const & srcLeft,
          TSourceRight const & srcRight,
          TTraceLeft const  traceLeft,
          TTraceRight const traceRight,
          TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const & /*tag*/)
{
    auto cmp = cmpGt(srcRight, srcLeft);
    target = blend(srcLeft, srcRight, cmp);
    return blend(traceLeft, traceRight, cmp);
}

template <typename TTarget,
          typename TSourceLeft,
          typename TSourceRight,
          typename TTraceLeft,
          typename TTraceRight,
          typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TTarget>>>, typename TraceBitMap_<TTarget>::Type)
_maxScore(TTarget & target,
          TSourceLeft const & srcLeft,
          TSourceRight const & srcRight,
          TTraceLeft const  traceLeft,
          TTraceRight const traceRight,
          TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> > const & /*tag*/)
{
    return (srcLeft == srcRight)
        ? (target = srcLeft, traceLeft | traceRight)
        : (srcLeft < srcRight)
            ? (target = srcRight, traceRight)
            : (target = srcLeft, traceLeft);
}

template <typename TTarget,
          typename TSourceLeft,
          typename TSourceRight,
          typename TTraceLeft,
          typename TTraceRight,
          typename TGapsPlacement>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TTarget>>, typename TraceBitMap_<TTarget>::Type)
_maxScore(TTarget & target,
          TSourceLeft const & srcLeft,
          TSourceRight const & srcRight,
          TTraceLeft const  traceLeft,
          TTraceRight const traceRight,
          TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> > const & /*tag*/)
{
    auto cmpG = cmpGt(srcRight, srcLeft);
    auto cmpE = cmpEq(srcRight, srcLeft);
    target = blend(srcLeft, srcRight, cmpG);
    auto result = blend(traceLeft, traceRight, cmpG);
    return blend(result, traceLeft | traceRight, cmpE);
}

// ----------------------------------------------------------------------------
// Function _computeScore                          [RecursionDirectionDiagonal]
// ----------------------------------------------------------------------------
// Independent of gap cost model.
template <typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline auto
_computeScore(DPCell_<TScoreValue, TGapCosts> & current,
              DPCell_<TScoreValue, TGapCosts> & diagonal,
              DPCell_<TScoreValue, TGapCosts> const & /*horizontal*/,
              DPCell_<TScoreValue, TGapCosts> const & /*vertical*/,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionDiagonal const &,
              DPProfile_<TAlgorithm, TGapCosts, TTracebackConfig, TExecPolicy> const &)
{
    _scoreOfCell(current) = _scoreOfCell(diagonal) + score(scoringScheme, seqHVal, seqVVal);
    setGapExtension(current, False(), False(), createVector<TScoreValue>(-1));

    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        return _maxScore(_scoreOfCell(current),
                         TraceBitMap_<TScoreValue>::NONE,
                         _scoreOfCell(current),
                         TraceBitMap_<TScoreValue>::NONE,
                         TraceBitMap_<TScoreValue>::DIAGONAL,
                         TTracebackConfig{});
    }
    else
    {
        return TraceBitMap_<TScoreValue>::DIAGONAL;
    }
}

// ----------------------------------------------------------------------------
// Function _computeScore                              [RecursionDirectionZero]
// ----------------------------------------------------------------------------

// Independent of gap cost model.
template <typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgoTag, typename TTraceFlag, typename TExecPolicy>
inline auto
_computeScore(DPCell_<TScoreValue, TGapCosts> & current,
              DPCell_<TScoreValue, TGapCosts> & diagonal,
              DPCell_<TScoreValue, TGapCosts> const & /*horizontal*/,
              DPCell_<TScoreValue, TGapCosts> & vertical,
              TSequenceHValue const & /*seqHVal*/,
              TSequenceVValue const & /*seqVVal*/,
              TScoringScheme const & /*scoringScheme*/,
              RecursionDirectionZero const &,
              DPProfile_<TAlgoTag, TGapCosts, TTraceFlag, TExecPolicy> const &)
{
    _scoreOfCell(current) = TraceBitMap_<TScoreValue>::NONE;
    _scoreOfCell(diagonal) = _scoreOfCell(current);
    _scoreOfCell(vertical) = _scoreOfCell(current);
    return TraceBitMap_<TScoreValue>::NONE;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_H_
