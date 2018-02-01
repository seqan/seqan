// ==========================================================================
//                    test_alignment_dp_matrix_navigator.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_MATRIX_NAVIGATOR_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_MATRIX_NAVIGATOR_H_

#include <algorithm>

#include <seqan/basic.h>
#include <seqan/align.h>


template <typename TSpec>
void testAlignmentDPMatrixTraceNavigatorConstructorDefault(TSpec const &)
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, FullDPMatrix> TDPMatrix;
    typedef typename Iterator<TDPMatrix, Standard>::Type TIterator;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TSpec>, NavigateColumnWise> dpTraceMatrixNavigator;

    // Test if default constructor sets NULL pointer
    bool resultPointer = true;
    if (dpTraceMatrixNavigator._ptrDataContainer)
        resultPointer = false;

    bool resultActiveColIter = true;
    if (dpTraceMatrixNavigator._activeColIterator != TIterator())
        resultActiveColIter = false;

    SEQAN_ASSERT_EQ(resultPointer, true);
    SEQAN_ASSERT_EQ(resultActiveColIter, true);
    SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 0);
}

void testAlignmentDPMatrixNavigatorScoreMarixInitUnbanded()
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, FullDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, TDPCell{});

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        dpScoreMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -10);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
}

void testAlignmentDPMatrixNavigatorScoreMarixSparseInitUnbanded()
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, SparseDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;
    typedef DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> TDPMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, TDPCell{});

    TDPMatrixNavigator dpScoreMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
    SEQAN_ASSERT_EQ(value(dpScoreMatrixNavigator._activeColIterator), TDPCell{});
    SEQAN_ASSERT_EQ(value(dpScoreMatrixNavigator._prevColIterator), TDPCell{});
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -9);
}

void testAlignmentDPMatrixNavigatorScoreMarixInitBanded()
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, FullDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;
    typedef DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWiseBanded> TDPMatrixNavigator;


    { // Case1: Band intersects with poit of origin.
        TDPMatrix dpMatrix;
        setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix, TDPCell{});

        TDPMatrixNavigator dpScoreMatrixNavigator{dpMatrix, DPBandConfig<BandOn>{-4, 3}};

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 4);
    }

    {
        TDPMatrix dpMatrix2;
        setLength(dpMatrix2, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix2, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix2);

        value(host(dpMatrix2), 0) = 10;

        TDPMatrixNavigator dpScoreMatrixNavigator{dpMatrix2, DPBandConfig<BandOn>{0, 7}};

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix2, Standard()), 7);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix2, Standard()), -1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 8);
    }

    {
        TDPMatrix dpMatrix3;
        setLength(dpMatrix3, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix3, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix3);

        TDPMatrixNavigator dpScoreMatrixNavigator{dpMatrix3, DPBandConfig<BandOn>{-7, 0}};

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix3, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix3, Standard()), -8);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }
}

void testAlignmentDPMatrixNavigatorScoreMarixSparseInitBanded()
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, SparseDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;
    typedef DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWiseBanded> TDPMatrixNavigator;

    {
        TDPMatrix dpMatrix;
        setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix, TDPCell{});

        TDPMatrixNavigator dpScoreMatrixNavigator{dpMatrix, DPBandConfig<BandOn>{-4, 3}};

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -4);
    }

    {
        TDPMatrix dpMatrix2;
        setLength(dpMatrix2, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix2, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix2);

        TDPMatrixNavigator dpScoreMatrixNavigator{dpMatrix2, DPBandConfig<BandOn>{0, 7}};

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix2, Standard()), 7);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix2, Standard()), 7);  // Behind the last cell.
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 0);
    }

    {
        TDPMatrix dpMatrix3;
        setLength(dpMatrix3, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix3, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix3);

        TDPMatrixNavigator dpScoreMatrixNavigator{dpMatrix3, DPBandConfig<BandOn>{-7, 0}};

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix3, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix3, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -7);
    }
}

void testAlignmentDPMatrixNavigatorScoreMarixGoNextCell()
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, FullDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 4);
    resize(dpMatrix, TDPCell{});

    int score = 0;
    auto generator = [&score]()
    {
        TDPCell cell;
        cell = score++;
        return cell;
    };

    std::generate(begin(host(dpMatrix), Standard()), end(host(dpMatrix), Standard()), generator);

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        navi{dpMatrix, DPBandConfig<BandOff>{}};

    // Initialize with 0.
    _setScoreOfCell(value(navi), 0);
    // Move along initial column.
    _goNextCell(navi, MetaColumnDescriptor<DPInitialColumn, FullColumn>{}, FirstCell{});
    SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 0);
    for (int v_pos = 1; v_pos < 3; ++v_pos)
    {
        _goNextCell(navi, MetaColumnDescriptor<DPInitialColumn, FullColumn>{}, InnerCell{});
        SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), v_pos);
    }
    _goNextCell(navi, MetaColumnDescriptor<DPInitialColumn, FullColumn>{}, LastCell{});
    SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 3);

    // Move along inner columns.
    for (int h_pos = 1; h_pos < 3; ++h_pos)
    {
        _goNextCell(navi, MetaColumnDescriptor<DPInnerColumn, FullColumn>{}, FirstCell{});
        SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), h_pos * 4);
        SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), (h_pos - 1) * 4);
        for (int v_pos = 1; v_pos < 3; ++v_pos)
        {
            _goNextCell(navi, MetaColumnDescriptor<DPInnerColumn, FullColumn>{}, InnerCell{});
            SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), h_pos * 4 + v_pos);
            SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), (h_pos - 1) * 4 + v_pos);
        }
        _goNextCell(navi, MetaColumnDescriptor<DPInnerColumn, FullColumn>{}, LastCell{});
        SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), h_pos * 4 + 3);
        SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), (h_pos - 1) * 4 + 3);
    }

    // Move along last column.
    _goNextCell(navi, MetaColumnDescriptor<DPFinalColumn, FullColumn>{}, FirstCell{});
    SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 12);
    SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), 8);
    for (int v_pos = 1; v_pos < 3; ++v_pos)
    {
        _goNextCell(navi, MetaColumnDescriptor<DPFinalColumn, FullColumn>{}, InnerCell{});
        SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 12 + v_pos);
        SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), 8 + v_pos);
    }
    _goNextCell(navi, MetaColumnDescriptor<DPFinalColumn, FullColumn>{}, LastCell{});
    SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 15);
    SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), 11);
}

void testAlignmentDPMatrixNavigatorScoreMarixSparseGoNextCell()
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, SparseDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 4);
    resize(dpMatrix, TDPCell{});

    int score = 0;
    auto generator = [&score]()
    {
        TDPCell cell;
        cell = score++;
        return cell;
    };

    std::generate(begin(host(dpMatrix), Standard()), end(host(dpMatrix), Standard()), generator);

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        navi{dpMatrix, DPBandConfig<BandOff>()};

    // Initialize with 0.
    _setScoreOfCell(value(navi), 0);
    // Move along initial column.
    _goNextCell(navi, MetaColumnDescriptor<DPInitialColumn, FullColumn>{}, FirstCell{});
    SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 0);
    for (int v_pos = 1; v_pos < 3; ++v_pos)
    {
        _goNextCell(navi, MetaColumnDescriptor<DPInitialColumn, FullColumn>{}, InnerCell{});
        SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), v_pos);
    }
    _goNextCell(navi, MetaColumnDescriptor<DPInitialColumn, FullColumn>{}, LastCell{});
    SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 3);

    // Move along inner columns.
    for (int h_pos = 1; h_pos < 3; ++h_pos)
    {
        _goNextCell(navi, MetaColumnDescriptor<DPInnerColumn, FullColumn>{}, FirstCell{});
        SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 0);
        SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), 0);
        for (int v_pos = 1; v_pos < 3; ++v_pos)
        {
            _goNextCell(navi, MetaColumnDescriptor<DPInnerColumn, FullColumn>{}, InnerCell{});
            SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), v_pos);
            SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), v_pos);
        }
        _goNextCell(navi, MetaColumnDescriptor<DPInnerColumn, FullColumn>{}, LastCell{});
        SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 3);
        SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), 3);
    }

    // Move along last column.
    _goNextCell(navi, MetaColumnDescriptor<DPFinalColumn, FullColumn>{}, FirstCell{});
    SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 0);
    SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), 0);
    for (int v_pos = 1; v_pos < 3; ++v_pos)
    {
        _goNextCell(navi, MetaColumnDescriptor<DPFinalColumn, FullColumn>{}, InnerCell{});
        SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), v_pos);
        SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), v_pos);
    }
    _goNextCell(navi, MetaColumnDescriptor<DPFinalColumn, FullColumn>{}, LastCell{});
    SEQAN_ASSERT_EQ(_scoreOfCell(value(navi)), 3);
    SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(navi)), 3);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorAssignValue(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, TDPMatrixSpec> TDPMatrix;
    typedef typename Value<TDPMatrix>::Type TDPCell;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, TDPCell{});

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        dpScoreMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};
    TDPCell tmp;
    tmp = 20;

    assignValue(dpScoreMatrixNavigator, tmp);
    SEQAN_ASSERT_EQ(_scoreOfCell(value(dpScoreMatrixNavigator._activeColIterator)), 20);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorValue(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, TDPMatrixSpec> TDPMatrix;
    typedef typename Value<TDPMatrix>::Type TDPCell;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, TDPCell{});

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        dpScoreMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    TDPCell tmp;
    tmp = 20;

    assignValue(dpScoreMatrixNavigator, tmp);
    SEQAN_ASSERT_EQ(_scoreOfCell(value(dpScoreMatrixNavigator)), 20);

    const DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        dpScoreMatrixNavigatorConst(dpScoreMatrixNavigator);
    SEQAN_ASSERT_EQ(_scoreOfCell(value(dpScoreMatrixNavigatorConst)), 20);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorPreviousCellHorizontal(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, TDPMatrixSpec> TDPMatrix;
    typedef typename Value<TDPMatrix>::Type TDPCell;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, TDPCell{});

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        dpScoreMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    // A hacky way to test the previous cell function.
    dpScoreMatrixNavigator._prevColIterator = dpScoreMatrixNavigator._activeColIterator;
    *dpScoreMatrixNavigator._activeColIterator = 20;
    SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(dpScoreMatrixNavigator)), 20);

    const DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        dpScoreMatrixNavigatorConst(dpScoreMatrixNavigator);
    SEQAN_ASSERT_EQ(_scoreOfCell(previousCellHorizontal(dpScoreMatrixNavigatorConst)), 20);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorCoordinate(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, TDPMatrixSpec> TDPMatrix;
    typedef typename Value<TDPMatrix>::Type TDPCell;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, TDPCell{});

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        dpScoreMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    dpScoreMatrixNavigator._activeColIterator += 7;
    dpScoreMatrixNavigator._prevColIterator += 7;

    SEQAN_ASSERT_EQ(coordinate(dpScoreMatrixNavigator, +DPMatrixDimension_::HORIZONTAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpScoreMatrixNavigator, +DPMatrixDimension_::VERTICAL), 7u);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorContainer(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<DPCell_<int, LinearGaps>, TDPMatrixSpec> TDPMatrix;
    typedef typename Value<TDPMatrix>::Type TDPCell;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, TDPCell{});

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        dpScoreMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    SEQAN_ASSERT_EQ(&container(dpScoreMatrixNavigator), &dpMatrix);

    const DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise>
        dpScoreMatrixNavigatorConst(dpScoreMatrixNavigator);
    SEQAN_ASSERT_EQ(&container(dpScoreMatrixNavigatorConst), &dpMatrix);
}


// ----------------------------------------------------------------------------
// Test constructor                                  [DPScoreMatrix, FullDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_init_unbanded)
{
    testAlignmentDPMatrixNavigatorScoreMarixInitUnbanded();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_init_banded)
{
    testAlignmentDPMatrixNavigatorScoreMarixInitBanded();
}

// ----------------------------------------------------------------------------
// Test functions                                  [DPScoreMatrix, FullDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_value)
{
    testAlignmentDPScoreMatrixNavigatorValue(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_previous_cell_horizontal)
{
    testAlignmentDPScoreMatrixNavigatorPreviousCellHorizontal(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_go_next_cell)
{
    testAlignmentDPMatrixNavigatorScoreMarixGoNextCell();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_assign_value)
{
    testAlignmentDPScoreMatrixNavigatorAssignValue(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_coordinate)
{
    testAlignmentDPScoreMatrixNavigatorCoordinate(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_container)
{
    testAlignmentDPScoreMatrixNavigatorContainer(seqan::FullDPMatrix());
}

// ----------------------------------------------------------------------------
// Test constructor                                [DPScoreMatrix, SparseDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_init_unbanded)
{
    testAlignmentDPMatrixNavigatorScoreMarixSparseInitUnbanded();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_init_banded)
{
    testAlignmentDPMatrixNavigatorScoreMarixSparseInitBanded();
}

// ----------------------------------------------------------------------------
// Test functions                                [DPScoreMatrix, SparseDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_value)
{
    testAlignmentDPScoreMatrixNavigatorValue(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_previous_cell_horizontal)
{
    testAlignmentDPScoreMatrixNavigatorPreviousCellHorizontal(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_go_next)
{
    testAlignmentDPMatrixNavigatorScoreMarixSparseGoNextCell();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_assign_value)
{
    testAlignmentDPScoreMatrixNavigatorAssignValue(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_coordinate)
{
    testAlignmentDPScoreMatrixNavigatorCoordinate(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_container)
{
    testAlignmentDPScoreMatrixNavigatorContainer(seqan::SparseDPMatrix());
}

// ----------------------------------------------------------------------------
// Test functions                                  [DPTraceMatrix, FullDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_init_unbanded)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix);

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise>
        dpTraceMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix);
    SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
    SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_disabled_init_unbanded)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix);

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOff>, NavigateColumnWise>
        dpTraceMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    SEQAN_ASSERT_NEQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix);
    SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 0);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_init_banded)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    { // Case1: Band intersects with poit of origin.
        TDPMatrix dpMatrix;
        setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix);

        DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWiseBanded>
            dpTraceMatrixNavigator{dpMatrix, DPBandConfig<BandOn>{-4, 3}};

        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 4);
    }

    {
        TDPMatrix dpMatrix2;
        setLength(dpMatrix2, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix2, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix2);

        value(host(dpMatrix2), 0) = 10;
        DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWiseBanded>
            dpTraceMatrixNavigator{dpMatrix2, DPBandConfig<BandOn>{0, 7}};

        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix2);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix2, Standard()), 7);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 8);
    }

    {
        TDPMatrix dpMatrix3;
        setLength(dpMatrix3, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix3, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix3);


        DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWiseBanded>
            dpTraceMatrixNavigator{dpMatrix3, DPBandConfig<BandOn>{-7, 0}};

        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix3, Standard()), 0);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_disabled_init_banded)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    { // Case1: Band intersects with poit of origin.
        TDPMatrix dpMatrix;
        setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix);

        DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOff>, NavigateColumnWiseBanded>
            dpTraceMatrixNavigator{dpMatrix, DPBandConfig<BandOn>{-4, 3}};

        SEQAN_ASSERT_NEQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 0);
    }

    {
        TDPMatrix dpMatrix2;
        setLength(dpMatrix2, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix2, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix2);

        value(host(dpMatrix2), 0) = 10;
        DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOff>, NavigateColumnWiseBanded>
            dpTraceMatrixNavigator{dpMatrix2, DPBandConfig<BandOn>{0, 7}};

        SEQAN_ASSERT_NEQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix2);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 0);
    }

    {
        TDPMatrix dpMatrix3;
        setLength(dpMatrix3, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix3, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix3);


        DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOff>, NavigateColumnWiseBanded>
            dpTraceMatrixNavigator{dpMatrix3, DPBandConfig<BandOn>{-7, 0}};

        SEQAN_ASSERT_NEQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 0);
    }
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_go_next)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 4);
    resize(dpMatrix);

    std::iota(begin(host(dpMatrix), Standard()), end(host(dpMatrix), Standard()), 0);

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise>
        navi{dpMatrix, DPBandConfig<BandOff>()};

    // Move along initial column.
    _goNextCell(navi, MetaColumnDescriptor<DPInitialColumn, FullColumn>{}, FirstCell{});
    SEQAN_ASSERT_EQ(value(navi), 0);
    for (int v_pos = 1; v_pos < 3; ++v_pos)
    {
        _goNextCell(navi, MetaColumnDescriptor<DPInitialColumn, FullColumn>{}, InnerCell{});
        SEQAN_ASSERT_EQ(value(navi), v_pos);
    }
    _goNextCell(navi, MetaColumnDescriptor<DPInitialColumn, FullColumn>{}, LastCell{});
    SEQAN_ASSERT_EQ(value(navi), 3);

    // Move along inner columns.
    for (int h_pos = 1; h_pos < 3; ++h_pos)
    {
        _goNextCell(navi, MetaColumnDescriptor<DPInnerColumn, FullColumn>{}, FirstCell{});
        SEQAN_ASSERT_EQ(value(navi), h_pos * 4);
        for (int v_pos = 1; v_pos < 3; ++v_pos)
        {
            _goNextCell(navi, MetaColumnDescriptor<DPInnerColumn, FullColumn>{}, InnerCell{});
            SEQAN_ASSERT_EQ(value(navi), h_pos * 4 + v_pos);
        }
        _goNextCell(navi, MetaColumnDescriptor<DPInnerColumn, FullColumn>{}, LastCell{});
        SEQAN_ASSERT_EQ(value(navi), h_pos * 4 + 3);
    }

    // Move along last column.
    _goNextCell(navi, MetaColumnDescriptor<DPFinalColumn, FullColumn>{}, FirstCell{});
    SEQAN_ASSERT_EQ(value(navi), 12);
    for (int v_pos = 1; v_pos < 3; ++v_pos)
    {
        _goNextCell(navi, MetaColumnDescriptor<DPFinalColumn, FullColumn>{}, InnerCell{});
        SEQAN_ASSERT_EQ(value(navi), 12 + v_pos);
    }
    _goNextCell(navi, MetaColumnDescriptor<DPFinalColumn, FullColumn>{}, LastCell{});
    SEQAN_ASSERT_EQ(value(navi), 15);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_assign_value)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise>
        dpTraceMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    assignValue(dpTraceMatrixNavigator, 20);
    SEQAN_ASSERT_EQ(value(dpTraceMatrixNavigator._activeColIterator), 20);
}


SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_value)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise>
        dpTraceMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    assignValue(dpTraceMatrixNavigator, 20);
    SEQAN_ASSERT_EQ(value(dpTraceMatrixNavigator), 20);
    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> const
        dpTraceMatrixNavigatorConst(dpTraceMatrixNavigator);
    SEQAN_ASSERT_EQ(value(dpTraceMatrixNavigatorConst), 20);
}


SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_coordinate)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise>
        dpTraceMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    dpTraceMatrixNavigator._activeColIterator += 7;

    SEQAN_ASSERT_EQ(coordinate(dpTraceMatrixNavigator, +DPMatrixDimension_::HORIZONTAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpTraceMatrixNavigator, +DPMatrixDimension_::VERTICAL), 7u);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_container)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise>
        dpTraceMatrixNavigator{dpMatrix, DPBandConfig<BandOff>{}};

    SEQAN_ASSERT_EQ(&container(dpTraceMatrixNavigator), &dpMatrix);

    const DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpTraceMatrixNavigatorConst(dpTraceMatrixNavigator);
    SEQAN_ASSERT_EQ(&container(dpTraceMatrixNavigatorConst), &dpMatrix);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_to_global_position)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix);

    DPMatrixNavigator_<TDPMatrix,
                       DPTraceMatrix<TracebackOn<> >,
                       NavigateColumnWise> navi{dpMatrix, DPBandConfig<BandOff>{}};

    for (unsigned i = 0; i < length(dpMatrix, +DPMatrixDimension_::HORIZONTAL); ++i)
    {
        for (unsigned j = 0; j < length(dpMatrix, +DPMatrixDimension_::VERTICAL); ++j)
        {
            _setToPosition(navi, toGlobalPosition(navi, i, j));
            SEQAN_ASSERT_EQ(i, coordinate(navi, +DPMatrixDimension_::HORIZONTAL));
            SEQAN_ASSERT_EQ(j, coordinate(navi, +DPMatrixDimension_::VERTICAL));
        }
    }
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_MATRIX_NAVIGATOR_H_
