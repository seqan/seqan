// ==========================================================================
//                        test_alignment_dp_formula.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_FORMULA_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_FORMULA_H_


#include <seqan/basic.h>

#include <seqan/score.h>
#include <seqan/align.h>

template <typename TBand>
void testDPFormulaNoTraceLocalLinearDiagonalDirection(TBand const &)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOff> TDPProfile;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionDiagonal(), TDPProfile());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
    SEQAN_ASSERT_EQ(activeCell._score, 4);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 10);

    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                               scoringScheme, RecursionDirectionDiagonal(), TDPProfile());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
    SEQAN_ASSERT_EQ(activeCell._score, 0);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 10);
}


SEQAN_DEFINE_TEST(test_dp_formula_trace_global_linear_diagonal_direction)
{
    using namespace seqan;

    typedef TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = -10;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);
    TTraceValue traceValue;
    {
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }

    {
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                   scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                   scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_linear_horizontal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = -10;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionHorizontal(), TDPProfileSingleTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
    SEQAN_ASSERT_EQ(activeCell._score, 6);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 6);

    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionHorizontal(), TDPProfileCompleteTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
    SEQAN_ASSERT_EQ(activeCell._score, 6);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 6);
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_linear_vertical_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = -10;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionVertical(), TDPProfileSingleTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
    SEQAN_ASSERT_EQ(activeCell._score, 6);

    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
    SEQAN_ASSERT_EQ(activeCell._score, 2);
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_linear_upper_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = -10;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);

        prevDiagonal._score = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);
    }

    {
        prevDiagonal._score = -10;
        prevHorizontal._score = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, -8);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, -8);
    }

    {
        prevHorizontal._score = -4;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -4);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -4);
        SEQAN_ASSERT_EQ(prevVertical._score, -8);

        prevDiagonal._score = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -4);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -4);
        SEQAN_ASSERT_EQ(prevVertical._score, -8);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_linear_lower_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = -10;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }

    {
        prevVertical._score = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, -10);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, -10);
    }

    {
        prevVertical._score = -8;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, -8);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, -8);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_linear_all_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = -10;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::HORIZONTAL )|
                        TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX);
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);
    }

    {
        prevDiagonal._score = -10;
        prevVertical._score = 5;
        prevHorizontal._score = -5;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 1);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -5);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -5);
        SEQAN_ASSERT_EQ(prevVertical._score, 1);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -5);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -5);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);
    }

    {
        prevDiagonal._score = -10;
        prevVertical._score = -10;
        prevHorizontal._score = 5;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 1);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 5);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 5);
        SEQAN_ASSERT_EQ(prevVertical._score, 1);

        prevDiagonal._score = -10;
        prevVertical._score = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 1);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 5);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 5);
        SEQAN_ASSERT_EQ(prevVertical._score, 1);
    }

    {
        prevDiagonal._score = 4;
        prevVertical._score = 2;
        prevHorizontal._score = -4;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -4);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -4);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);

        prevDiagonal._score = 4;
        prevHorizontal._score = 10;
        prevVertical._score = -4;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);
    }

    {
        prevDiagonal._score = 10;
        prevVertical._score = 10;
        prevHorizontal._score = -12;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -12);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -12);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -12);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -12);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = 10;
        prevVertical._score = 8;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -12);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -12);
        SEQAN_ASSERT_EQ(prevVertical._score, 8);

        prevDiagonal._score = 10;
        prevVertical._score = 12;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -12);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -12);
        SEQAN_ASSERT_EQ(prevVertical._score, 8);
    }

    {
        prevDiagonal._score = 10;
        prevVertical._score = 12;
        prevHorizontal._score = 12;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 12);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 12);
        SEQAN_ASSERT_EQ(prevVertical._score, 8);

        prevDiagonal._score = 10;
        prevVertical._score = 12;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::HORIZONTAL |
                                     TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 12);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 12);
        SEQAN_ASSERT_EQ(prevVertical._score, 8);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_affine_diagonal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = -10;
    prevDiagonal._verticalScore = -10;
    prevDiagonal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._verticalScore = -10;
    prevHorizontal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = -10;
    prevVertical._horizontalScore = -10;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, -12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_affine_horizontal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = -10;
    prevDiagonal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._horizontalScore = -10;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }

    {
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevHorizontal._horizontalScore = 10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);

        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_affine_vertical_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = -10;
    prevDiagonal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = -10;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }

    {
        prevVertical._score = 10;
        prevVertical._verticalScore = 10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);

        prevVertical._score = 10;
        prevVertical._verticalScore = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_affine_upper_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = -10;
    prevDiagonal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = -10;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        prevHorizontal._horizontalScore = 16;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        prevHorizontal._horizontalScore = 16;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL |
                                     TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        prevHorizontal._score = 18;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 18);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 18);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX |
                                     TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 18);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 18);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        prevHorizontal._horizontalScore = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 18);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 18);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN |
                                     TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 18);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 18);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_affine_lower_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = -10;
    prevDiagonal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = -10;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = 16;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = 16;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = 16;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = 16;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = 16;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 18);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = 16;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX |
                                     TraceBitMap_<>::VERTICAL | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 18);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 18);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN |
                                     TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 18);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_affine_all_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    // We need to distinguish between all traces and single trace.
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = -10;
    prevDiagonal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = -10;
    prevHorizontal._horizontalScore = -10;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = -10;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = 16;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -14);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = 16;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -14);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = 16;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -14);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 10;
        prevVertical._verticalScore = 16;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL |
                                     TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -14);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = 16;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -14);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = 16;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::VERTICAL_OPEN |
                                     TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -14);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN |
                                     TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -14);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);

        prevDiagonal._score = 10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN |
                                     TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 12);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -14);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 12);
    }

    {
        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = -10;
        prevHorizontal._horizontalScore = 20;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL_OPEN |
                                     TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);

        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL_OPEN |
                                     TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, -10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, -10);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);
    }

    {
        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = -10;
        prevHorizontal._score = 22;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX | TraceBitMap_<>::HORIZONTAL |
                                     TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 22);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 22);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);

        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 18;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX | TraceBitMap_<>::HORIZONTAL |
                                     TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 22);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 22);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);
    }

    {
        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 22;
        prevVertical._verticalScore = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX |
                                     TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 22);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 22);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);

        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 22;
        prevVertical._verticalScore = -10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX |
                                     TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 22);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 22);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);
    }

    {
        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 22;
        prevVertical._verticalScore = 20;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL |
                                     TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 22);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 22);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);

        prevDiagonal._score = -10;
        prevDiagonal._horizontalScore = -10;
        prevVertical._score = 22;
        prevVertical._verticalScore = 20;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL |
                                     TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL_OPEN |
                                     TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 22);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 22);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);
    }

    {
        prevDiagonal._score = 14;
        prevVertical._score = 22;
        prevVertical._verticalScore = 20;
        prevHorizontal._score = 22;
        prevHorizontal._horizontalScore = 20;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 22);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 22);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);

        prevDiagonal._score = 14;
        prevVertical._score = 22;
        prevVertical._verticalScore = 20;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL |
                                     TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::VERTICAL_OPEN |
                                     TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX |
                                     TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 16);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 16);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 22);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 22);
        SEQAN_ASSERT_EQ(prevVertical._score, 16);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_dynamic_diagonal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, DynamicGaps> activeCell;
    DPCell_<int, DynamicGaps> prevDiagonal;
    prevDiagonal._score = -10;
    DPCell_<int, DynamicGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, DynamicGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);
    TTraceValue traceValue;
    {
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(_scoreOfCell(activeCell), -8);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevDiagonal), -10);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevHorizontal), 10);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevVertical), 10);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(_scoreOfCell(activeCell), -8);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevDiagonal), -10);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevHorizontal), 10);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevVertical), 10);
    }

    {
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                   scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(_scoreOfCell(activeCell), -12);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevDiagonal), -10);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevHorizontal), 10);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevVertical), 10);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                   scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(_scoreOfCell(activeCell), -12);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevDiagonal), -10);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevHorizontal), 10);
        SEQAN_ASSERT_EQ(_scoreOfCell(prevVertical), 10);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_dynamic_horizontal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;
    typedef DPCell_<int, DynamicGaps> TDPCell;

    TDPCell activeCell;
    TDPCell prevDiagonal;
    prevDiagonal._score = -10;
    TDPCell prevHorizontal;
    prevHorizontal._score = 10;
    TDPCell prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionHorizontal(), TDPProfileSingleTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
    SEQAN_ASSERT_EQ(_scoreOfCell(activeCell), 6);
    SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
    SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                           scoringScheme, RecursionDirectionHorizontal(), TDPProfileCompleteTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
    SEQAN_ASSERT_EQ(_scoreOfCell(activeCell), 6);
    SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
    SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

    _setBit(prevHorizontal, True(), DynamicGapExtensionHorizontal());
    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                           scoringScheme, RecursionDirectionHorizontal(), TDPProfileCompleteTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
    SEQAN_ASSERT_EQ(_scoreOfCell(activeCell), 8);
    SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
    SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_dynamic_vertical_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;
    typedef DPCell_<int, DynamicGaps> TCell;

    TCell activeCell;
    TCell prevDiagonal;
    prevDiagonal._score = -10;
    TCell prevHorizontal;
    prevHorizontal._score = 10;
    TCell prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionVertical(), TDPProfileSingleTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
    SEQAN_ASSERT_EQ(activeCell._score, 6);

    prevDiagonal._score = -10;
    prevVertical._score = 10;
    _setBit(prevVertical, False(), DynamicGapExtensionVertical());
    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
    SEQAN_ASSERT_EQ(activeCell._score, 6);

    prevDiagonal._score = -10;
    prevVertical._score = 10;
    _setBit(prevVertical, True(), DynamicGapExtensionVertical());
    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
    SEQAN_ASSERT_EQ(activeCell._score, 8);
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_dynamic_upper_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;
    typedef DPCell_<int, DynamicGaps> TDPCell;

    TDPCell activeCell;
    TDPCell prevDiagonal;
    prevDiagonal._score = -10;
    TDPCell prevHorizontal;
    prevHorizontal._score = 10;
    TDPCell prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        _setBit(prevHorizontal, True(), DynamicGapExtensionHorizontal());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 8);
    }

    {
        prevDiagonal._score = -10;
        prevVertical._score = 10;
        prevHorizontal._score = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
    }

    {
        prevDiagonal._score = -10;
        prevVertical._score = 10;
        prevHorizontal._score = -4;
        _setBit(prevHorizontal, False(), DynamicGapExtensionHorizontal());
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        _setBit(prevHorizontal, False(), DynamicGapExtensionVertical());
        prevHorizontal._score = -6;
        _setBit(prevHorizontal, True(), DynamicGapExtensionHorizontal());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                                       scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        _setBit(prevHorizontal, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        prevDiagonal._score = -10;
        prevVertical._score = 10;
        _setBit(prevHorizontal, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_dynamic_lower_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;
    typedef DPCell_<int, DynamicGaps> TDPCell;

    TDPCell activeCell;
    TDPCell prevDiagonal;
    prevDiagonal._score = -10;
    TDPCell prevHorizontal;
    prevHorizontal._score = 10;
    TDPCell prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);

        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 8);
    }

    {
        prevVertical._score = -10;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
    }

    {
        prevVertical._score = -4;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        prevVertical._score = -6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                                       scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_global_dynamic_all_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<GlobalAlignment_<>, DynamicGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;
    typedef DPCell_<int, DynamicGaps> TDPCell;

    TDPCell activeCell;
    TDPCell prevDiagonal;
    prevDiagonal._score = -10;
    TDPCell prevHorizontal;
    prevHorizontal._score = -10;
    TDPCell prevVertical;
    prevVertical._score = -10;

    Score<int, Simple> scoringScheme(2, -2, -2, -4);

    {  // From diagonal only.
        // Single trace.
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, -8);

        // Complete trace.
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, -8);
    }

    {  // From horizontal only.
        // Single Trace + horizontal open.
        prevHorizontal._score = 0;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, -4);

        // Complete Trace + horizontal open.
        prevDiagonal._score = -10;
        prevVertical._score = -10;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, -4);

        // Single Trace + horizontal extend.
        prevDiagonal._score = -10;
        prevVertical._score = -10;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        _setBit(prevHorizontal, True(), DynamicGapExtensionHorizontal());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, -2);

        // Complete Trace + horizontal extend.
        prevDiagonal._score = -10;
        prevVertical._score = -10;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, -2);
    }

    {  // From vertical only.
        // Single Trace + vertical open.
        prevDiagonal._score = -10;
        prevVertical._score = 4;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 0);

        prevDiagonal._score = -10;
        prevVertical._score = 4;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 0);

        // Single Trace + horizontal extend.
        prevDiagonal._score = -10;
        prevVertical._score = 4;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 2);

        // Complete Trace + horizontal extend.
        prevDiagonal._score = -10;
        prevVertical._score = 4;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 2);
    }

    {  // From horizontal + diagonal.
        // Single trace + horizontal extend.
        prevHorizontal._score = 6;
        prevDiagonal._score = 2;
        prevVertical._score = 4;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

        // Complete trace + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 4;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

        // Single Trace + horizontal open
        prevDiagonal._score = 2;
        prevVertical._score = 4;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        prevHorizontal._score = 8;
        _setBit(prevHorizontal, False(), DynamicGapExtensionHorizontal());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

        // Complete trace + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 4;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);
    }

    {  // From horziontal + vertical.
        // Single trace + vertical extend + horizontal open.
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Complete trace + vertical extend + horizontal open.
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Single trace + vertical open + horizontal open.
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Complete trace + vertical open + horizontal open.
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Single trace + vertical open + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        _setBit(prevHorizontal, True(), DynamicGapExtensionHorizontal());
        prevHorizontal._score = 6;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Complete trace + vertical open + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Single trace + vertical extend + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Complete trace + vertical extend + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'C',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);
    }

    {  // From vertical + diagonal.
        // Single trace + vertical extend.
        prevHorizontal._score = 0;
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

        // Complete trace + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Single Trace + vertical open
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

        // Complete trace + vertical extend.
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);
    }

    {  // From horizontal + vertical + diagonal.
        // Single trace + vertical open + horizontal extend.
        prevHorizontal._score = 6;
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

        // Complete trace + vertical open + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Single trace + vertical open + horizontal open.
        _setBit(prevHorizontal, False(), DynamicGapExtensionHorizontal());
        prevHorizontal._score = 8;
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

        // Complete trace + vertical open + horizontal open.
        prevDiagonal._score = 2;
        prevVertical._score = 8;
        _setBit(prevVertical, False(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Single trace + vertical extend + horizontal open.
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

        // Complete trace + vertical open + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);

        // Single trace + vertical extend + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        _setBit(prevHorizontal, True(), DynamicGapExtensionHorizontal());
        prevHorizontal._score = 6;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), false);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), false);

        // Complete trace + vertical extend + horizontal extend.
        prevDiagonal._score = 2;
        prevVertical._score = 6;
        _setBit(prevVertical, True(), DynamicGapExtensionVertical());
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionHorizontal()), true);
        SEQAN_ASSERT_EQ(isGapExtension(activeCell, DynamicGapExtensionVertical()), true);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_linear_diagonal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
    SEQAN_ASSERT_EQ(activeCell._score, 4);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 10);

    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
    SEQAN_ASSERT_EQ(activeCell._score, 4);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 10);

    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
    SEQAN_ASSERT_EQ(activeCell._score, 0);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 10);

    prevVertical._score = 10;
    traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
    SEQAN_ASSERT_EQ(activeCell._score, 0);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 10);
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_linear_horizontal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);
    }
    {
        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevHorizontal._score = 2;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);
    }
}
SEQAN_DEFINE_TEST(test_dp_formula_trace_local_linear_vertical_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 6);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 6);
    }
    {
        prevDiagonal._score = 2;
        prevVertical._score = 2;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);

        prevDiagonal._score = 2;
        prevVertical._score = 2;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);
    }
}
SEQAN_DEFINE_TEST(test_dp_formula_trace_local_linear_upper_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 8;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 2;

    Score<int, Simple> scoringScheme(2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 8);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 8);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevDiagonal._score = 2;
        prevVertical._score = 2;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 8);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 8);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }
    {
        prevHorizontal._score = 2;
        prevDiagonal._score = 2;
        prevVertical._score = 2;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);

        prevDiagonal._score = 2;
        prevVertical._score = 2;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_linear_lower_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 2;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 8;

    Score<int, Simple> scoringScheme(2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 8);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 8);
    }
    {
        prevDiagonal._score = 2;
        prevVertical._score = 2;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 2);

        prevDiagonal._score = 2;
        prevVertical._score = 2;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 2);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_linear_all_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 8;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 8;

    Score<int, Simple> scoringScheme(2, -2, -4);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 8);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 8);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevDiagonal._score = 2;
        prevVertical._score = 8;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL |
                                     TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 8);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 8);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }
    {
        prevDiagonal._score = 0;
        prevVertical._score = 2;
        prevHorizontal._score = 2;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);

        prevDiagonal._score = 0;
        prevVertical._score = 2;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 2);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_affine_diagonal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = 4;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = 4;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);
    int inf = DPCellDefaultInfinity<DPCell_<int, AffineGaps> >::VALUE;

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._verticalScore, inf);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._verticalScore, inf);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }


    {
        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 4;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 4;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_affine_horizontal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = 8;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = 4;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 4;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionHorizontal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }

    {
        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 4;
        prevHorizontal._score = 1;
        prevHorizontal._horizontalScore = 1;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -3);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 1);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 1);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 4;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionHorizontal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -3);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 1);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 1);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_affine_vertical_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = 8;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = 8;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);
    int inf = DPCellDefaultInfinity<DPCell_<int, AffineGaps> >::VALUE;

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 8;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL | TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }

    {
        prevDiagonal._score = 2;
        prevVertical._score = 1;
        prevVertical._verticalScore = 1;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);

        prevDiagonal._score = 2;
        prevVertical._score = 1;
        prevVertical._verticalScore = 1;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionVertical(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_affine_upper_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = 8;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = 8;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 8;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX |
                                     TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }

    {
        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 8;
        prevHorizontal._score = 1;
        prevHorizontal._horizontalScore = 1;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -3);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 1);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 1);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 8;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                   scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, -3);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 1);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 1);
        SEQAN_ASSERT_EQ(prevVertical._score, 0);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_affine_lower_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = 2;
    prevDiagonal._horizontalScore = 4;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = 8;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = 8;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);
    int inf = DPCellDefaultInfinity<DPCell_<int, AffineGaps> >::VALUE;

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL | TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 8;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX |
                                     TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 10);
    }

    {
        prevDiagonal._score = 2;
        prevVertical._score = 1;
        prevVertical._verticalScore = 1;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 1);

        prevDiagonal._score = 2;
        prevVertical._score = 1;
        prevVertical._verticalScore = 1;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                   scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, inf);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 1);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_trace_local_affine_all_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPProfileSingleTrace;
    typedef DPProfile_<LocalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsLeft> > > TDPProfileCompleteTrace;

    DPCell_<int, AffineGaps> activeCell;
    DPCell_<int, AffineGaps> prevDiagonal;
    prevDiagonal._score = 2;
    prevDiagonal._horizontalScore = 4;
    DPCell_<int, AffineGaps> prevHorizontal;
    prevHorizontal._score = 10;
    prevHorizontal._horizontalScore = 8;
    DPCell_<int, AffineGaps> prevVertical;
    prevVertical._score = 10;
    prevVertical._verticalScore = 8;

    Score<int, Simple> scoringScheme(2, -2, -4, -6);

    {
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                               scoringScheme, RecursionDirectionAll(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);

        prevDiagonal._score = 2;
        prevVertical._score = 10;
        prevVertical._verticalScore = 8;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                   scoringScheme, RecursionDirectionAll(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(traceValue, (TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX |
                                     TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX |
                                     TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN |
                                     TraceBitMap_<>::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(activeCell._score, 4);
        SEQAN_ASSERT_EQ(activeCell._horizontalScore, 4);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
        SEQAN_ASSERT_EQ(prevVertical._score, 4);
    }

    {
        prevDiagonal._score = 2;
        prevHorizontal._score = 1;
        prevHorizontal._horizontalScore = 1;
        prevVertical._score = 1;
        prevVertical._verticalScore = 1;
        TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileSingleTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 1);
        SEQAN_ASSERT_EQ(prevVertical._score, 1);

        prevDiagonal._score = 2;
        prevVertical._score = 1;
        prevVertical._verticalScore = 1;
        traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'C', 'A',
                                               scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfileCompleteTrace());

        SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(activeCell._score, 0);
        SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
        SEQAN_ASSERT_EQ(prevHorizontal._score, 1);
        SEQAN_ASSERT_EQ(prevVertical._score, 1);
    }
}

SEQAN_DEFINE_TEST(test_dp_formula_notrace_diagonal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOff> TDPProfile;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionDiagonal(), TDPProfile());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
    SEQAN_ASSERT_EQ(activeCell._score, 4);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 10);
}

SEQAN_DEFINE_TEST(test_dp_formula_notrace_horizontal_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOff> TDPProfile;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionHorizontal(), TDPProfile());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
    SEQAN_ASSERT_EQ(activeCell._score, 6);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 6);
}

SEQAN_DEFINE_TEST(test_dp_formula_notrace_vertical_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOff> TDPProfile;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 2;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionVertical(), TDPProfile());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
    SEQAN_ASSERT_EQ(activeCell._score, 6);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 2);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 6);
}

SEQAN_DEFINE_TEST(test_dp_formula_notrace_upper_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOff> TDPProfile;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 4;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionUpperDiagonal(), TDPProfile());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
    SEQAN_ASSERT_EQ(activeCell._score, 6);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 6);
}

SEQAN_DEFINE_TEST(test_dp_formula_notrace_lower_band_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOff> TDPProfile;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 4;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionLowerDiagonal(), TDPProfile());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
    SEQAN_ASSERT_EQ(activeCell._score, 6);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 4);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 10);
}

SEQAN_DEFINE_TEST(test_dp_formula_notrace_all_direction)
{
    using namespace seqan;

    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPProfile_<LocalAlignment_<>, LinearGaps, TracebackOff> TDPProfile;

    DPCell_<int, LinearGaps> activeCell;
    DPCell_<int, LinearGaps> prevDiagonal;
    prevDiagonal._score = 4;
    DPCell_<int, LinearGaps> prevHorizontal;
    prevHorizontal._score = 10;
    DPCell_<int, LinearGaps> prevVertical;
    prevVertical._score = 10;

    Score<int, Simple> scoringScheme(2, -2, -4);

    TTraceValue traceValue = _computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, 'A', 'A',
                                           scoringScheme, RecursionDirectionAll(), TDPProfile());

    SEQAN_ASSERT_EQ(static_cast<int>(traceValue), static_cast<int>(TraceBitMap_<>::NONE));
    SEQAN_ASSERT_EQ(activeCell._score, 6);
    SEQAN_ASSERT_EQ(prevDiagonal._score, 10);
    SEQAN_ASSERT_EQ(prevHorizontal._score, 10);
    SEQAN_ASSERT_EQ(prevVertical._score, 6);
}


#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_FORMULA_H_
