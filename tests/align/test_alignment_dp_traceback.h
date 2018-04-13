// ==========================================================================
//                       test_alignment_dp_traceback.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_TRACEBACK_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_TRACEBACK_H_

#include <seqan/basic.h>
#include <seqan/align.h>

SEQAN_DEFINE_TEST(test_align2_traceback_affine)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<> > TDPProfile;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;
    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TTraceMatrix;

    typedef DPMatrixNavigator_<TTraceMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> TDPTraceNavigator;

    String<TTraceSegment> target;
    TTraceMatrix traceMatrix;

    setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 4);

    resize(traceMatrix);

    value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 0) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 2, 0) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 3, 0) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;

    value(traceMatrix, 0, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 1, 1) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 2, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::VERTICAL_OPEN;
    value(traceMatrix, 3, 1) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;

    value(traceMatrix, 0, 2) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 1, 2) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 2, 2) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 3, 2) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;

    value(traceMatrix, 0, 3) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 1, 3) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 2, 3) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 3, 3) = TraceBitMap_<>::DIAGONAL;

    TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOff>()};

    DnaString str0 = "ACG";
    DnaString str1 = "ACG";

    DPScout_<int, Default> dpScout;
    dpScout._maxHostPosition = 15;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());
    SEQAN_ASSERT_EQ(length(target), 1u);
    SEQAN_ASSERT_EQ(target[0], TTraceSegment(0, 0, 3, TraceBitMap_<>::DIAGONAL));

    clear(target);
    dpScout._maxHostPosition = 14;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());
    SEQAN_ASSERT_EQ(length(target), 3u);
    SEQAN_ASSERT_EQ(target[0], TTraceSegment(3, 2, 1, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(0, 2, 3, TraceBitMap_<>::HORIZONTAL));
    SEQAN_ASSERT_EQ(target[2], TTraceSegment(0, 0, 2, TraceBitMap_<>::VERTICAL));

    clear(target);
    dpScout._maxHostPosition = 13;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());
    SEQAN_ASSERT_EQ(length(target), 4u);
    SEQAN_ASSERT_EQ(target[0], TTraceSegment(3, 1, 2, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(2, 1, 1, TraceBitMap_<>::HORIZONTAL));
    SEQAN_ASSERT_EQ(target[2], TTraceSegment(1, 0, 1, TraceBitMap_<>::DIAGONAL));
    SEQAN_ASSERT_EQ(target[3], TTraceSegment(0, 0, 1, TraceBitMap_<>::HORIZONTAL));

    clear(target);
    dpScout._maxHostPosition = 12;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());
    SEQAN_ASSERT_EQ(length(target), 2u);
    SEQAN_ASSERT_EQ(target[0], TTraceSegment(3, 0, 3, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(0, 0, 3, TraceBitMap_<>::HORIZONTAL));

    clear(target);
    dpScout._maxHostPosition = 11;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());
    SEQAN_ASSERT_EQ(length(target), 4u);
    SEQAN_ASSERT_EQ(target[0], TTraceSegment(2, 3, 1, TraceBitMap_<>::HORIZONTAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(2, 1, 2, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[2], TTraceSegment(1, 0, 1, TraceBitMap_<>::DIAGONAL));
    SEQAN_ASSERT_EQ(target[3], TTraceSegment(0, 0, 1, TraceBitMap_<>::HORIZONTAL));

    clear(target);
    dpScout._maxHostPosition = 7;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());
    SEQAN_ASSERT_EQ(length(target), 3u);
    SEQAN_ASSERT_EQ(target[0], TTraceSegment(1, 3, 2, TraceBitMap_<>::HORIZONTAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(1, 1, 2, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[2], TTraceSegment(0, 0, 1, TraceBitMap_<>::DIAGONAL));

    clear(target);
    dpScout._maxHostPosition = 3;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());
    SEQAN_ASSERT_EQ(length(target), 2u);
    SEQAN_ASSERT_EQ(target[0], TTraceSegment(0, 3, 3, TraceBitMap_<>::HORIZONTAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(0, 0, 3, TraceBitMap_<>::VERTICAL));
}

SEQAN_DEFINE_TEST(test_align2_traceback_linear_unbanded_alignment)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;
    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TTraceMatrix;

    typedef DPMatrixNavigator_<TTraceMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> TDPTraceNavigator;

    String<TTraceSegment> target;
    TTraceMatrix traceMatrix;

    setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 4);

    resize(traceMatrix);

    value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 0) = TraceBitMap_<>::VERTICAL;
    value(traceMatrix, 2, 0) = TraceBitMap_<>::VERTICAL;
    value(traceMatrix, 3, 0) = TraceBitMap_<>::VERTICAL;

    value(traceMatrix, 0, 1) = TraceBitMap_<>::HORIZONTAL;
    value(traceMatrix, 1, 1) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 2, 1) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 3, 1) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 2) = TraceBitMap_<>::HORIZONTAL;
    value(traceMatrix, 1, 2) = TraceBitMap_<>::NONE;
    value(traceMatrix, 2, 2) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 3, 2) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 3) = TraceBitMap_<>::HORIZONTAL;
    value(traceMatrix, 1, 3) = TraceBitMap_<>::NONE;
    value(traceMatrix, 2, 3) = TraceBitMap_<>::NONE;
    value(traceMatrix, 3, 3) = TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL  | TraceBitMap_<>::HORIZONTAL;

    TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOff>{}};

    DnaString str0 = "ACG";
    DnaString str1 = "ACG";

    DPScout_<int, Default> dpScout;
    dpScout._maxHostPosition = 15;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());

    SEQAN_ASSERT_EQ(target[0], TTraceSegment(2, 2, 1, TraceBitMap_<>::DIAGONAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(1, 2, 1, TraceBitMap_<>::HORIZONTAL));
    SEQAN_ASSERT_EQ(target[2], TTraceSegment(1, 1, 1, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[3], TTraceSegment(0, 0, 1, TraceBitMap_<>::DIAGONAL));
}

SEQAN_DEFINE_TEST(test_align2_traceback_linear_normal_banded_alignment)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;
    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TTraceMatrix;

    typedef DPMatrixNavigator_<TTraceMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWiseBanded> TDPTraceNavigator;

    String<TTraceSegment> target;
    TTraceMatrix traceMatrix;

    setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 5);
    setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 3);

    resize(traceMatrix);

    value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 0) = TraceBitMap_<>::VERTICAL;
    value(traceMatrix, 2, 0) = TraceBitMap_<>::VERTICAL;

    value(traceMatrix, 0, 1) = TraceBitMap_<>::HORIZONTAL;
    value(traceMatrix, 1, 1) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 2, 1) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;

    value(traceMatrix, 0, 2) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 2) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 2, 2) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 3) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 1, 3) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 2, 3) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 4) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 4) = TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::VERTICAL  | TraceBitMap_<>::HORIZONTAL;
    value(traceMatrix, 2, 4) = TraceBitMap_<>::NONE;


    TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOn>{-1, 1}};

    DnaString str0 = "ACGT";
    DnaString str1 = "ACGT";

    DPScout_<int, Default> dpScout;
    dpScout._maxHostPosition = 13;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOn>(-1, 1), TDPProfile());

    SEQAN_ASSERT_EQ(target[0], TTraceSegment(3, 3, 1, TraceBitMap_<>::DIAGONAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(3, 2, 1, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[2], TTraceSegment(1, 2, 2, TraceBitMap_<>::HORIZONTAL));
    SEQAN_ASSERT_EQ(target[3], TTraceSegment(1, 1, 1, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[4], TTraceSegment(0, 0, 1, TraceBitMap_<>::DIAGONAL));
}

SEQAN_DEFINE_TEST(test_align2_traceback_linear_wide_banded_alignment)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;
    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TTraceMatrix;

    typedef DPMatrixNavigator_<TTraceMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWiseBanded> TDPTraceNavigator;

    String<TTraceSegment> target;
    TTraceMatrix traceMatrix;

    setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 7);
    setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 7);

    resize(traceMatrix);

    value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 0) = TraceBitMap_<>::NONE;
    value(traceMatrix, 2, 0) = TraceBitMap_<>::NONE;
    value(traceMatrix, 3, 0) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 4, 0) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 5, 0) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 6, 0) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;

    value(traceMatrix, 0, 1) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 1) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 2, 1) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 3, 1) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 4, 1) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 5, 1) = TraceBitMap_<>::NONE;
    value(traceMatrix, 6, 1) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 2) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 1, 2) = TraceBitMap_<>::NONE;
    value(traceMatrix, 2, 2) = TraceBitMap_<>::NONE;
    value(traceMatrix, 3, 2) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 4, 2) = TraceBitMap_<>::NONE;
    value(traceMatrix, 5, 2) = TraceBitMap_<>::NONE;
    value(traceMatrix, 6, 2) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 3) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 3) = TraceBitMap_<>::NONE;
    value(traceMatrix, 2, 3) = TraceBitMap_<>::NONE;
    value(traceMatrix, 3, 3) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 4, 3) = TraceBitMap_<>::NONE;
    value(traceMatrix, 5, 3) = TraceBitMap_<>::NONE;
    value(traceMatrix, 6, 3) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 4) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 4) = TraceBitMap_<>::NONE;
    value(traceMatrix, 2, 4) = TraceBitMap_<>::NONE;
    value(traceMatrix, 3, 4) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 4, 4) = TraceBitMap_<>::NONE;
    value(traceMatrix, 5, 4) = TraceBitMap_<>::NONE;
    value(traceMatrix, 6, 4) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 5) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 5) = TraceBitMap_<>::NONE;
    value(traceMatrix, 2, 5) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 3, 5) = TraceBitMap_<>::NONE;
    value(traceMatrix, 4, 5) = TraceBitMap_<>::NONE;
    value(traceMatrix, 5, 5) = TraceBitMap_<>::NONE;
    value(traceMatrix, 6, 5) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 6) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 6) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 2, 6) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 3, 6) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 4, 6) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 5, 6) = TraceBitMap_<>::NONE;
    value(traceMatrix, 6, 6) = TraceBitMap_<>::NONE;

    TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOn>{-4, 4}};

    DnaString str0 = "ACGTAC";
    DnaString str1 = "ACGTAC";

    DPScout_<int, Default> dpScout;
    dpScout._maxHostPosition = 46;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOn>(-4, 4), TDPProfile());

    SEQAN_ASSERT_EQ(target[0], TTraceSegment(6, 3, 3, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(1, 3, 5, TraceBitMap_<>::HORIZONTAL));
    SEQAN_ASSERT_EQ(target[2], TTraceSegment(1, 1, 2, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[3], TTraceSegment(0, 0, 1, TraceBitMap_<>::DIAGONAL));
}

SEQAN_DEFINE_TEST(test_align2_traceback_linear_small_banded_alignment)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;
    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TTraceMatrix;

    typedef DPMatrixNavigator_<TTraceMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWiseBanded> TDPTraceNavigator;

    String<TTraceSegment> target;
    TTraceMatrix traceMatrix;

    setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 1);

    resize(traceMatrix);

    value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;

    value(traceMatrix, 0, 1) = TraceBitMap_<>::DIAGONAL;

    value(traceMatrix, 0, 2) = TraceBitMap_<>::DIAGONAL;

    value(traceMatrix, 0, 3) = TraceBitMap_<>::DIAGONAL;

    TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOn>{0, 0}};

    DnaString str0 = "ACG";
    DnaString str1 = "ACG";

    DPScout_<int, Default> dpScout;
    dpScout._maxHostPosition = 3;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOn>(0, 0), TDPProfile());

    SEQAN_ASSERT_EQ(target[0], TTraceSegment(0, 0, 3, TraceBitMap_<>::DIAGONAL));
}

SEQAN_DEFINE_TEST(test_align2_traceback_gaps_left_linear_gaps)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;
    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TTraceMatrix;

    typedef DPMatrixNavigator_<TTraceMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> TDPTraceNavigator;

    String<TTraceSegment> target;
    TTraceMatrix traceMatrix;

    setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 3);
    setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 4);

    resize(traceMatrix);

    value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 0) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 2, 0) = TraceBitMap_<>::VERTICAL;
    value(traceMatrix, 3, 0) = TraceBitMap_<>::VERTICAL;

    value(traceMatrix, 0, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 1, 1) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 2, 1) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 3, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::HORIZONTAL;

    value(traceMatrix, 0, 2) = TraceBitMap_<>::HORIZONTAL;
    value(traceMatrix, 1, 2) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 2, 2) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 3, 2) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::DIAGONAL;


    TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOff>()};

    DnaString str0 = "AC";
    DnaString str1 = "CCC";

    DPScout_<int, Default> dpScout;
    dpScout._maxHostPosition = 11;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());

    SEQAN_ASSERT_EQ(target[0], TTraceSegment(0, 1, 2, TraceBitMap_<>::DIAGONAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(0, 0, 1, TraceBitMap_<>::VERTICAL));
}

SEQAN_DEFINE_TEST(test_align2_traceback_gaps_right_linear_gaps)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsRight> > > TDPProfile;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;
    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TTraceMatrix;

    typedef DPMatrixNavigator_<TTraceMatrix, DPTraceMatrix<TracebackOn<TracebackConfig_<CompleteTrace, GapsRight> > >, NavigateColumnWise> TDPTraceNavigator;

    String<TTraceSegment> target;
    TTraceMatrix traceMatrix;

    setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 3);
    setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 4);

    resize(traceMatrix);

    value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
    value(traceMatrix, 1, 0) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 2, 0) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 3, 0) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;

    value(traceMatrix, 0, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 1, 1) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 2, 1) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 3, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;

    value(traceMatrix, 0, 2) = TraceBitMap_<>::HORIZONTAL | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
    value(traceMatrix, 1, 2) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
    value(traceMatrix, 2, 2) = TraceBitMap_<>::DIAGONAL;
    value(traceMatrix, 3, 2) = TraceBitMap_<>::VERTICAL | TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;


    TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOff>()};

    DnaString str0 = "AC";
    DnaString str1 = "CCC";

    DPScout_<int, Default> dpScout;
    dpScout._maxHostPosition = 11;
    _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());

    SEQAN_ASSERT_EQ(target[0], TTraceSegment(2, 2, 1, TraceBitMap_<>::VERTICAL));
    SEQAN_ASSERT_EQ(target[1], TTraceSegment(0, 0, 2, TraceBitMap_<>::DIAGONAL));
}

SEQAN_DEFINE_TEST(test_align2_traceback_gaps_left_affine_gaps)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<> > TDPProfile;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;
    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TTraceMatrix;

    typedef DPMatrixNavigator_<TTraceMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> TDPTraceNavigator;

    String<TTraceSegment> target;
    TTraceMatrix traceMatrix;


    {   // Tests gaps at end
        setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 3);
        setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 4);

        resize(traceMatrix);

        value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
        value(traceMatrix, 1, 0) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
        value(traceMatrix, 2, 0) = TraceBitMap_<>::VERTICAL;
        value(traceMatrix, 3, 0) = TraceBitMap_<>::VERTICAL;

        value(traceMatrix, 0, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
        value(traceMatrix, 1, 1) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 2, 1) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 3, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::HORIZONTAL;

        value(traceMatrix, 0, 2) = TraceBitMap_<>::HORIZONTAL;
        value(traceMatrix, 1, 2) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
        value(traceMatrix, 2, 2) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 3, 2) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::DIAGONAL;


        TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOff>()};

        DnaString str0 = "AC";
        DnaString str1 = "CCC";

        DPScout_<int, Default> dpScout;
        dpScout._maxHostPosition = 11;
        _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());

        // TODO(rmaerker): This is disabled by default for the affine gap costs.
        SEQAN_ASSERT_EQ(target[0], TTraceSegment(2, 2, 1, TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(target[1], TTraceSegment(0, 0, 2, TraceBitMap_<>::DIAGONAL));
    }

    {   // Tests inner gaps.
        clear(target);
        setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 5);
        setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 4);

        resize(traceMatrix);

        value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
        value(traceMatrix, 1, 0) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
        value(traceMatrix, 2, 0) = TraceBitMap_<>::VERTICAL;
        value(traceMatrix, 3, 0) = TraceBitMap_<>::VERTICAL;

        value(traceMatrix, 0, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
        value(traceMatrix, 1, 1) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 2, 1) = TraceBitMap_<>::NONE;
        value(traceMatrix, 3, 1) = TraceBitMap_<>::NONE;

        value(traceMatrix, 0, 2) = TraceBitMap_<>::HORIZONTAL;
        value(traceMatrix, 1, 2) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
        value(traceMatrix, 2, 2) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 3, 2) = TraceBitMap_<>::NONE;

        value(traceMatrix, 0, 3) = TraceBitMap_<>::HORIZONTAL;
        value(traceMatrix, 1, 3) = TraceBitMap_<>::NONE;
        value(traceMatrix, 2, 3) = TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
        value(traceMatrix, 3, 3) = TraceBitMap_<>::NONE;

        value(traceMatrix, 0, 4) = TraceBitMap_<>::HORIZONTAL;
        value(traceMatrix, 1, 4) = TraceBitMap_<>::NONE;
        value(traceMatrix, 2, 4) = TraceBitMap_<>::NONE;
        value(traceMatrix, 3, 4) = TraceBitMap_<>::DIAGONAL;

        TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOff>()};

        DnaString str0 = "ACCA";
        DnaString str1 = "ACA";

        DPScout_<int, Default> dpScout;
        dpScout._maxHostPosition = 19;
        _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());

        // TODO(rmaerker): This is disabled by default for the affine gap costs.
        SEQAN_ASSERT_EQ(target[0], TTraceSegment(2, 1, 2, TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(target[1], TTraceSegment(1, 1, 1, TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(target[2], TTraceSegment(0, 0, 1, TraceBitMap_<>::DIAGONAL));
    }

}

SEQAN_DEFINE_TEST(test_align2_traceback_gaps_right_affine_gaps)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, AffineGaps, TracebackOn<TracebackConfig_<CompleteTrace, GapsRight> > > TDPProfile;
    typedef TraceSegment_<unsigned, unsigned> TTraceSegment;
    typedef typename TraceBitMap_<>::Type TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TTraceMatrix;

    typedef DPMatrixNavigator_<TTraceMatrix, DPTraceMatrix<TracebackOn<TracebackConfig_<CompleteTrace, GapsRight> > >, NavigateColumnWise> TDPTraceNavigator;

    String<TTraceSegment> target;
    TTraceMatrix traceMatrix;

    {   // Tests gaps at end
        setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 3);
        setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 4);

        resize(traceMatrix);

        value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
        value(traceMatrix, 1, 0) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
        value(traceMatrix, 2, 0) = TraceBitMap_<>::VERTICAL;
        value(traceMatrix, 3, 0) = TraceBitMap_<>::VERTICAL;

        value(traceMatrix, 0, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
        value(traceMatrix, 1, 1) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 2, 1) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 3, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::HORIZONTAL;

        value(traceMatrix, 0, 2) = TraceBitMap_<>::HORIZONTAL;
        value(traceMatrix, 1, 2) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
        value(traceMatrix, 2, 2) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 3, 2) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX | TraceBitMap_<>::DIAGONAL;


        TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOff>()};

        DnaString str0 = "AC";
        DnaString str1 = "CCC";

        DPScout_<int, Default> dpScout;
        dpScout._maxHostPosition = 11;
        _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());

        SEQAN_ASSERT_EQ(target[0], TTraceSegment(2, 2, 1, TraceBitMap_<>::VERTICAL));
        SEQAN_ASSERT_EQ(target[1], TTraceSegment(0, 0, 2, TraceBitMap_<>::DIAGONAL));
    }

    {   // Test inner gaps.
        clear(target);
        setLength(traceMatrix, DPMatrixDimension_::HORIZONTAL, 5);
        setLength(traceMatrix, DPMatrixDimension_::VERTICAL, 4);

        resize(traceMatrix);

        value(traceMatrix, 0, 0) = TraceBitMap_<>::NONE;
        value(traceMatrix, 1, 0) = TraceBitMap_<>::VERTICAL_OPEN | TraceBitMap_<>::MAX_FROM_VERTICAL_MATRIX;
        value(traceMatrix, 2, 0) = TraceBitMap_<>::VERTICAL;
        value(traceMatrix, 3, 0) = TraceBitMap_<>::VERTICAL;

        value(traceMatrix, 0, 1) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
        value(traceMatrix, 1, 1) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 2, 1) = TraceBitMap_<>::NONE;
        value(traceMatrix, 3, 1) = TraceBitMap_<>::NONE;

        value(traceMatrix, 0, 2) = TraceBitMap_<>::HORIZONTAL;
        value(traceMatrix, 1, 2) = TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
        value(traceMatrix, 2, 2) = TraceBitMap_<>::DIAGONAL;
        value(traceMatrix, 3, 2) = TraceBitMap_<>::NONE;

        value(traceMatrix, 0, 3) = TraceBitMap_<>::HORIZONTAL;
        value(traceMatrix, 1, 3) = TraceBitMap_<>::NONE;
        value(traceMatrix, 2, 3) = TraceBitMap_<>::DIAGONAL | TraceBitMap_<>::HORIZONTAL_OPEN | TraceBitMap_<>::MAX_FROM_HORIZONTAL_MATRIX;
        value(traceMatrix, 3, 3) = TraceBitMap_<>::NONE;

        value(traceMatrix, 0, 4) = TraceBitMap_<>::HORIZONTAL;
        value(traceMatrix, 1, 4) = TraceBitMap_<>::NONE;
        value(traceMatrix, 2, 4) = TraceBitMap_<>::NONE;
        value(traceMatrix, 3, 4) = TraceBitMap_<>::DIAGONAL;

        TDPTraceNavigator navigator{traceMatrix, DPBandConfig<BandOff>()};

        DnaString str0 = "ACCA";
        DnaString str1 = "ACA";

        DPScout_<int, Default> dpScout;
        dpScout._maxHostPosition = 19;
        _computeTraceback(target, navigator, dpScout, str0, str1, DPBandConfig<BandOff>(), TDPProfile());

        // TODO(rmaerker): This is disabled by default for the affine gap costs.
        SEQAN_ASSERT_EQ(target[0], TTraceSegment(3, 2, 1, TraceBitMap_<>::DIAGONAL));
        SEQAN_ASSERT_EQ(target[1], TTraceSegment(2, 2, 1, TraceBitMap_<>::HORIZONTAL));
        SEQAN_ASSERT_EQ(target[2], TTraceSegment(0, 0, 2, TraceBitMap_<>::DIAGONAL));
    }
}


#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_TRACEBACK_H_
