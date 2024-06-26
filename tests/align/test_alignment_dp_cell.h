// ==========================================================================
//                         test_alignment_dp_cell.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_CELL_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_CELL_H_

#include <seqan/basic.h>

#include <seqan/align.h>

template <typename TGapCosts>
void testDPCellValue(TGapCosts const &)
{
    using namespace seqan2;

    typedef DPCell_<int, TGapCosts> TDPCell;
    typedef DPCell_<int, TGapCosts> const TDPConstCell;

    typedef typename Value<TDPCell>::Type TDPCellResult;
    typedef typename Value<TDPConstCell>::Type TDPConstCellResult;
    bool result1 = IsSameType<TDPCellResult, int>::VALUE;
    bool result2 = IsSameType<TDPConstCellResult, const int>::VALUE;
    SEQAN_ASSERT_EQ(result1, true);
    SEQAN_ASSERT_EQ(result2, true);
}

template <typename TGapCosts>
void testDPCellReference(TGapCosts const &)
{
    using namespace seqan2;

    typedef DPCell_<int, TGapCosts> TDPCell;
    typedef DPCell_<int, TGapCosts> const TDPConstCell;

    typedef typename Reference<TDPCell>::Type TDPCellResult;
    typedef typename Reference<TDPConstCell>::Type TDPConstCellResult;
    bool result1 = IsSameType<TDPCellResult, int &>::VALUE;
    bool result2 = IsSameType<TDPConstCellResult, const int &>::VALUE;

    SEQAN_ASSERT_EQ(result1, true);
    SEQAN_ASSERT_EQ(result2, true);
}

template <typename TGapCosts>
void testDPCellDefaultInfinity(TGapCosts const &)
{
    using namespace seqan2;

    int    result_int       = DPCellDefaultInfinity<DPCell_<int, TGapCosts>>::VALUE;
    int    result_int_const = DPCellDefaultInfinity<DPCell_<int, TGapCosts> const>::VALUE;
    float  result_flt       = DPCellDefaultInfinity<DPCell_<float, TGapCosts>>::VALUE;
    float  result_flt_const = DPCellDefaultInfinity<DPCell_<float, TGapCosts> const>::VALUE;
    double result_dbl       = DPCellDefaultInfinity<DPCell_<double, TGapCosts>>::VALUE;
    double result_dbl_const = DPCellDefaultInfinity<DPCell_<double, TGapCosts> const>::VALUE;

    SEQAN_ASSERT_EQ      (result_int,       std::numeric_limits<int>::lowest() / 2);
    SEQAN_ASSERT_EQ      (result_int_const, std::numeric_limits<int>::lowest() / 2);
    SEQAN_ASSERT_IN_DELTA(result_flt,       std::numeric_limits<float>::lowest() / 2.0f, 1e-3);
    SEQAN_ASSERT_IN_DELTA(result_flt_const, std::numeric_limits<float>::lowest() / 2.0f, 1e-3);
    SEQAN_ASSERT_IN_DELTA(result_dbl,       std::numeric_limits<double>::lowest() / 2.0, 1e-3);
    SEQAN_ASSERT_IN_DELTA(result_dbl_const, std::numeric_limits<double>::lowest() / 2.0, 1e-3);
}

template <typename TGaps>
void testDPCellConstructor(TGaps const & /*tag*/)
{
    using namespace seqan2;

    DPCell_<int, TGaps> dpValue;

    int score = DPCellDefaultInfinity<DPCell_<int, TGaps> >::VALUE;
    SEQAN_ASSERT_EQ(dpValue._score, score);

}

template <typename TGaps>
void testDPCellCopyConstructor(TGaps const & /*tag*/)
{
    using namespace seqan2;

    DPCell_<int, TGaps> dpValue;

    DPCell_<int, TGaps> dpValue2 = dpValue;
    int score = DPCellDefaultInfinity<DPCell_<int, TGaps> >::VALUE;
    SEQAN_ASSERT_EQ(dpValue2._score, score);

    dpValue._score = 10;
    DPCell_<int, TGaps> dpValue3 = dpValue;
    SEQAN_ASSERT_EQ(dpValue3._score, 10);
}

template <typename TGaps>
void testDPCellAssignment(TGaps const & /*tag*/)
{
    using namespace seqan2;

    DPCell_<int, TGaps> dpValue;

    DPCell_<int, TGaps> dpValue2;
    int score = DPCellDefaultInfinity<DPCell_<int, TGaps> >::VALUE;
    SEQAN_ASSERT_EQ(dpValue2._score, score);


    dpValue._score = 10;
    dpValue2 = dpValue;
    SEQAN_ASSERT_EQ(dpValue2._score, 10);
}

void testDPCellAffineConstructor()
{
    using namespace seqan2;

    DPCell_<int, AffineGaps> dpValue;

    int score = DPCellDefaultInfinity<DPCell_<int, AffineGaps> >::VALUE;
    SEQAN_ASSERT_EQ(dpValue._score, score);
    SEQAN_ASSERT_EQ(dpValue._horizontalScore, score);
    SEQAN_ASSERT_EQ(dpValue._verticalScore, score);
}

void testDPCellAffineCopyConstructor()
{
    using namespace seqan2;

    DPCell_<int, AffineGaps> dpValue;

    DPCell_<int, AffineGaps> dpValue2 = dpValue;
    int score = DPCellDefaultInfinity<DPCell_<int, AffineGaps> >::VALUE;
    SEQAN_ASSERT_EQ(dpValue2._horizontalScore, score);
    SEQAN_ASSERT_EQ(dpValue2._verticalScore, score);
    SEQAN_ASSERT_EQ(dpValue2._score, score);

    dpValue._score = 10;
    DPCell_<int, AffineGaps> dpValue3 = dpValue;
    SEQAN_ASSERT_EQ(dpValue3._horizontalScore, score);
    SEQAN_ASSERT_EQ(dpValue3._verticalScore, score);
    SEQAN_ASSERT_EQ(dpValue3._score, 10);

    dpValue._horizontalScore = -10;
    dpValue._verticalScore = -20;

    DPCell_<int, AffineGaps> dpValue4 = dpValue;
    SEQAN_ASSERT_EQ(dpValue4._horizontalScore, -10);
    SEQAN_ASSERT_EQ(dpValue4._verticalScore, -20);
    SEQAN_ASSERT_EQ(dpValue4._score, 10);
}

void testDPCellAffineAssignment()
{
    using namespace seqan2;

    DPCell_<int, AffineGaps> dpValue;

    DPCell_<int, AffineGaps> dpValue2;
    int score = DPCellDefaultInfinity<DPCell_<int, AffineGaps> >::VALUE;
    SEQAN_ASSERT_EQ(dpValue2._score, score);
    SEQAN_ASSERT_EQ(dpValue2._horizontalScore, score);
    SEQAN_ASSERT_EQ(dpValue2._verticalScore, score);


    dpValue._score = 10;
    dpValue2 = dpValue;
    SEQAN_ASSERT_EQ(dpValue2._score, 10);
    SEQAN_ASSERT_EQ(dpValue2._horizontalScore, score);
    SEQAN_ASSERT_EQ(dpValue2._verticalScore, score);

    dpValue._horizontalScore = -10;
    dpValue._verticalScore = -20;

    dpValue2 = dpValue;
    SEQAN_ASSERT_EQ(dpValue2._score, 10);
    SEQAN_ASSERT_EQ(dpValue2._horizontalScore, -10);
    SEQAN_ASSERT_EQ(dpValue2._verticalScore, -20);
}

template <typename TGapCosts>
void testDPCellScore(TGapCosts const &)
{
    using namespace seqan2;

    DPCell_<int, TGapCosts> dpValue;

    int testScore = DPCellDefaultInfinity<DPCell_<int, LinearGaps> >::VALUE;
    SEQAN_ASSERT_EQ(_scoreOfCell(dpValue), testScore);

    dpValue._score = 10;
    SEQAN_ASSERT_EQ(_scoreOfCell(dpValue), 10);
}

void testDPCellVerticalScore()
{
    using namespace seqan2;

    DPCell_<int, AffineGaps> dpValue;

    int testScore = DPCellDefaultInfinity<DPCell_<int, AffineGaps> >::VALUE;
    SEQAN_ASSERT_EQ(_verticalScoreOfCell(dpValue), testScore);

    _verticalScoreOfCell(dpValue) = 10;
    SEQAN_ASSERT_EQ(_verticalScoreOfCell(dpValue), 10);

    DPCell_<int, AffineGaps> const dpValue2(dpValue);
    SEQAN_ASSERT_EQ(_verticalScoreOfCell(dpValue2), 10);
}


void testDPCellHorizontalScore()
{
    using namespace seqan2;

    DPCell_<int, AffineGaps> dpValue;

    int testScore = DPCellDefaultInfinity<DPCell_<int, AffineGaps> >::VALUE;
    SEQAN_ASSERT_EQ(_horizontalScoreOfCell(dpValue), testScore);

    _horizontalScoreOfCell(dpValue) = 10;
    SEQAN_ASSERT_EQ(_horizontalScoreOfCell(dpValue), 10);

    DPCell_<int, AffineGaps> const dpValue2(dpValue);
    SEQAN_ASSERT_EQ(_horizontalScoreOfCell(dpValue2), 10);
}

SEQAN_DEFINE_TEST(test_dp_cell_value)
{
    testDPCellValue(seqan2::LinearGaps());
    testDPCellValue(seqan2::AffineGaps());
    testDPCellValue(seqan2::DynamicGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_reference)
{
    testDPCellReference(seqan2::LinearGaps());
    testDPCellReference(seqan2::AffineGaps());
    testDPCellReference(seqan2::DynamicGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_default_infinity)
{
    testDPCellDefaultInfinity(seqan2::LinearGaps());
    testDPCellDefaultInfinity(seqan2::AffineGaps());
    testDPCellDefaultInfinity(seqan2::DynamicGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_linear_constructor)
{
    testDPCellConstructor(seqan2::LinearGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_linear_copy_constructor)
{
    testDPCellCopyConstructor(seqan2::LinearGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_linear_assignment)
{
    testDPCellAssignment(seqan2::LinearGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_linear_score)
{
    testDPCellScore(seqan2::LinearGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_affine_constructor)
{
    testDPCellAffineConstructor();
}

SEQAN_DEFINE_TEST(test_dp_cell_affine_copy_constructor)
{
    testDPCellAffineCopyConstructor();
}

SEQAN_DEFINE_TEST(test_dp_cell_affine_assignment)
{
    testDPCellAffineAssignment();
}

SEQAN_DEFINE_TEST(test_dp_cell_affine_score)
{
    testDPCellScore(seqan2::AffineGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_affine_vertical_score)
{
    testDPCellVerticalScore();
}

SEQAN_DEFINE_TEST(test_dp_cell_affine_horizontal_score)
{
    testDPCellHorizontalScore();
}

SEQAN_DEFINE_TEST(test_dp_cell_dynamic_constructor)
{
    testDPCellConstructor(seqan2::DynamicGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_dynamic_copy_constructor)
{
    testDPCellCopyConstructor(seqan2::DynamicGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_dynamic_assignment)
{
    testDPCellAssignment(seqan2::DynamicGaps());
}

SEQAN_DEFINE_TEST(test_dp_cell_dynamic_score)
{
    testDPCellScore(seqan2::DynamicGaps());
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_CELL_H_
