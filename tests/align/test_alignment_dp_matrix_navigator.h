// ==========================================================================
//                    test_alignment_dp_matrix_navigator.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_MATRIX_NAVIGATOR_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_MATRIX_NAVIGATOR_H_

#include <seqan/basic.h>

#include <seqan/align.h>


template <typename TSpec>
void testAlignmentDPMatrixScoreNavigatorConstructorDefault(TSpec const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, TSpec> TDPMatrix;
    typedef typename Value<TDPMatrix>::Type TDPCell;
    typedef typename Iterator<TDPMatrix, Standard>::Type TIterator;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    // Test if default constructor sets NULL pointer
    bool resultPointer = true;
    if (dpScoreMatrixNavigator._ptrDataContainer)
        resultPointer = false;

    bool resultPrevColIter = true;
    if (dpScoreMatrixNavigator._prevColIterator != TIterator())
        resultPrevColIter = false;

    bool resultActiveColIter = true;
    if (dpScoreMatrixNavigator._activeColIterator != TIterator())
        resultActiveColIter = false;

    SEQAN_ASSERT_EQ(resultPointer, true);
    SEQAN_ASSERT_EQ(resultPrevColIter, true);
    SEQAN_ASSERT_EQ(resultActiveColIter, true);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 0);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
}

template <typename TSpec>
void testAlignmentDPMatrixTraceNavigatorConstructorDefault(TSpec const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;
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

void testAlignmentDPMatrixNavigatorScoreMarixSparseContructor()
{
    using namespace seqan;

    typedef DPMatrix_<int, SparseDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;
    typedef Iterator<TDPMatrix, Standard>::Type TIterator;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    // Test if default constructor sets NULL pointer
    bool resultPointer = true;
    if (dpScoreMatrixNavigator._ptrDataContainer)
        resultPointer = false;

    bool resultActiveColIter = true;
    if (dpScoreMatrixNavigator._activeColIterator != TIterator())
        resultActiveColIter = false;

    bool resultPrevColIter = true;
    if (dpScoreMatrixNavigator._prevColIterator != TIterator())
        resultPrevColIter = false;

    SEQAN_ASSERT_EQ(resultPointer, true);
    SEQAN_ASSERT_EQ(resultActiveColIter, true);
    SEQAN_ASSERT_EQ(resultPrevColIter, true);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 0);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
}


void testAlignmentDPMatrixNavigatorScoreMarixInitUnbanded()
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 0);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -10);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
}

void testAlignmentDPMatrixNavigatorScoreMarixSparseInitUnbanded()
{
    using namespace seqan;

    typedef DPMatrix_<int, SparseDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;
    typedef DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> TDPMatrixNavigator;

    TDPMatrixNavigator dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 0);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
    SEQAN_ASSERT_EQ(value(dpScoreMatrixNavigator._activeColIterator), 0);
    SEQAN_ASSERT_EQ(value(dpScoreMatrixNavigator._prevColIterator), 0);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -9);
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
}

void testAlignmentDPMatrixNavigatorScoreMarixInitBanded()
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;
    typedef DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> TDPMatrixNavigator;

    TDPMatrixNavigator dpScoreMatrixNavigator;

    { // Case1: Band intersects with poit of origin.
        TDPMatrix dpMatrix;
        setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix, 0);

        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOn>(-4, 3));

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
    }

    {
        TDPMatrix dpMatrix2;
        setLength(dpMatrix2, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix2, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix2);

        value(host(dpMatrix2), 0) = 10;
        _init(dpScoreMatrixNavigator, dpMatrix2, DPBandConfig<BandOn>(0, 7));

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix2, Standard()), 7);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix2, Standard()), -1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 8);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
    }

    {
        TDPMatrix dpMatrix3;
        setLength(dpMatrix3, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix3, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix3);


        _init(dpScoreMatrixNavigator, dpMatrix3, DPBandConfig<BandOn>(-7, 0));

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix3, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix3, Standard()), -8);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
    }

}

void testAlignmentDPMatrixNavigatorScoreMarixSparseInitBanded()
{
    using namespace seqan;

    typedef DPMatrix_<int, SparseDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;
    typedef DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> TDPMatrixNavigator;

    TDPMatrixNavigator dpScoreMatrixNavigator;

    {
        TDPMatrix dpMatrix;
        setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix, 0);

        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOn>(-4, 3));

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
    }

    {
        TDPMatrix dpMatrix2;
        setLength(dpMatrix2, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix2, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix2);

        _init(dpScoreMatrixNavigator, dpMatrix2, DPBandConfig<BandOn>(0, 7));

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix2, Standard()), 7);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix2, Standard()), 7);  // Behind the last cell.
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
    }

    {
        TDPMatrix dpMatrix3;
        setLength(dpMatrix3, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix3, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix3);

        _init(dpScoreMatrixNavigator, dpMatrix3, DPBandConfig<BandOn>(-7, 0));

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix3, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix3, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -7);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
    }
}

void testAlignmentDPMatrixNavigatorScoreMarixGoNextCell()
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);
    resize(dpMatrix, 0);
    host(dpMatrix)[0] = 0;
    host(dpMatrix)[1] = 1;
    host(dpMatrix)[2] = 2;
    host(dpMatrix)[3] = 3;
    host(dpMatrix)[4] = 4;
    host(dpMatrix)[5] = 5;

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), -1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpScoreMatrixNavigator._activeColIterator += 2;
        dpScoreMatrixNavigator._prevColIterator += 2;
        ++dpScoreMatrixNavigator._laneLeap;  // For the test we simulate as if we were in a band.

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpScoreMatrixNavigator._activeColIterator += 2;
        dpScoreMatrixNavigator._prevColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpScoreMatrixNavigator._activeColIterator += 2;
        dpScoreMatrixNavigator._prevColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 2);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpScoreMatrixNavigator._activeColIterator += 2;
        dpScoreMatrixNavigator._prevColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpScoreMatrixNavigator._activeColIterator += 2;
        dpScoreMatrixNavigator._prevColIterator += 2;
        ++dpScoreMatrixNavigator._laneLeap;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpScoreMatrixNavigator._activeColIterator += 2;
        dpScoreMatrixNavigator._prevColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpScoreMatrixNavigator._activeColIterator += 2;
        dpScoreMatrixNavigator._prevColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 2);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpScoreMatrixNavigator._activeColIterator += 2;
        dpScoreMatrixNavigator._prevColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), FirstCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), InnerCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 4);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, 1);
    }
}

void testAlignmentDPMatrixNavigatorScoreMarixSparseGoNext()
{
    using namespace seqan;

    typedef DPMatrix_<int, SparseDPMatrix> TDPMatrix;
    typedef Value<TDPMatrix>::Type TDPCell;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 3);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);
    resize(dpMatrix, 3);

    host(dpMatrix)[1] = 1;
    host(dpMatrix)[2] = 2;

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);     // Was never set to 0.

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {

        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        // Need to update iterator just for the test.
        dpScoreMatrixNavigator._activeColIterator += 3;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        // Need to update Iterator just for the test.
        dpScoreMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        // Need to update Iterator just for the test.
        dpScoreMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        // Need to set iterator to correct position just for the test.
        --dpScoreMatrixNavigator._prevColIterator;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), LastCell());
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        dpScoreMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        dpScoreMatrixNavigator._activeColIterator += 3;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -3);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        dpScoreMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        dpScoreMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        --dpScoreMatrixNavigator._prevColIterator;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }

    {
        _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        assignValue(dpScoreMatrixNavigator._activeColIterator, 0);

        dpScoreMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), FirstCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), InnerCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 0);

        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), LastCell());

        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._laneLeap, -2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, 1);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, 2);
        SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, 1);
    }
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorAssignValue(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<int, TDPMatrixSpec> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    assignValue(dpScoreMatrixNavigator, 20);
    SEQAN_ASSERT_EQ(value(dpScoreMatrixNavigator._activeColIterator), 20);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorValue(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<int, TDPMatrixSpec> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    assignValue(dpScoreMatrixNavigator, 20);
    SEQAN_ASSERT_EQ(value(dpScoreMatrixNavigator), 20);

    const DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigatorConst(dpScoreMatrixNavigator);
    SEQAN_ASSERT_EQ(value(dpScoreMatrixNavigatorConst), 20);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorPreviousCellDiagonal(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<int, TDPMatrixSpec> TDPMatrix;
    typedef typename Value<TDPMatrix>::Type TDPCell;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellDiagonal, TDPCell());

    dpScoreMatrixNavigator._prevCellDiagonal = 20;
    SEQAN_ASSERT_EQ(previousCellDiagonal(dpScoreMatrixNavigator), 20);

    const DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigatorConst(dpScoreMatrixNavigator);
    SEQAN_ASSERT_EQ(previousCellDiagonal(dpScoreMatrixNavigatorConst), 20);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorPreviousCellHorizontal(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<int, TDPMatrixSpec> TDPMatrix;
    typedef typename Value<TDPMatrix>::Type TDPCell;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellHorizontal, TDPCell());

    dpScoreMatrixNavigator._prevCellHorizontal = 20;
    SEQAN_ASSERT_EQ(previousCellHorizontal(dpScoreMatrixNavigator), 20);

    const DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigatorConst(dpScoreMatrixNavigator);
    SEQAN_ASSERT_EQ(previousCellHorizontal(dpScoreMatrixNavigatorConst), 20);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorPreviousCellVertical(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<int, TDPMatrixSpec> TDPMatrix;
    typedef typename Value<TDPMatrix>::Type TDPCell;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
    SEQAN_ASSERT_EQ(dpScoreMatrixNavigator._prevCellVertical, TDPCell());

    dpScoreMatrixNavigator._prevCellVertical = 20;
    SEQAN_ASSERT_EQ(previousCellVertical(dpScoreMatrixNavigator), 20);

    const DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigatorConst(dpScoreMatrixNavigator);
    SEQAN_ASSERT_EQ(previousCellVertical(dpScoreMatrixNavigatorConst), 20);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorCoordinate(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<int, TDPMatrixSpec> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    dpScoreMatrixNavigator._activeColIterator += 7;
    dpScoreMatrixNavigator._prevColIterator += 7;

    SEQAN_ASSERT_EQ(coordinate(dpScoreMatrixNavigator, +DPMatrixDimension_::HORIZONTAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpScoreMatrixNavigator, +DPMatrixDimension_::VERTICAL), 7u);
}

template <typename TDPMatrixSpec>
void testAlignmentDPScoreMatrixNavigatorContainer(TDPMatrixSpec const)
{
    using namespace seqan;

    typedef DPMatrix_<int, TDPMatrixSpec> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    SEQAN_ASSERT_EQ(&container(dpScoreMatrixNavigator), &dpMatrix);

    const DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, NavigateColumnWise> dpScoreMatrixNavigatorConst(dpScoreMatrixNavigator);
    SEQAN_ASSERT_EQ(&container(dpScoreMatrixNavigatorConst), &dpMatrix);
}


// ----------------------------------------------------------------------------
// Test constructor                                  [DPScoreMatrix, FullDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_constructor)
{
    testAlignmentDPMatrixScoreNavigatorConstructorDefault(seqan::FullDPMatrix());
}

// ----------------------------------------------------------------------------
// Test functions                                  [DPScoreMatrix, FullDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_init_unbanded)
{
    testAlignmentDPMatrixNavigatorScoreMarixInitUnbanded();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_init_banded)
{
    testAlignmentDPMatrixNavigatorScoreMarixInitBanded();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_go_next_cell)
{
    testAlignmentDPMatrixNavigatorScoreMarixGoNextCell();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_assign_value)
{
    testAlignmentDPScoreMatrixNavigatorAssignValue(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_value)
{
    testAlignmentDPScoreMatrixNavigatorValue(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_previous_cell_diagonal)
{
    testAlignmentDPScoreMatrixNavigatorPreviousCellDiagonal(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_previous_cell_horizontal)
{
    testAlignmentDPScoreMatrixNavigatorPreviousCellHorizontal(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_full_previous_cell_vertical)
{
    testAlignmentDPScoreMatrixNavigatorPreviousCellVertical(seqan::FullDPMatrix());
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

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_constructor)
{
    testAlignmentDPMatrixScoreNavigatorConstructorDefault(seqan::SparseDPMatrix());
}

// ----------------------------------------------------------------------------
// Test functions                                [DPScoreMatrix, SparseDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_init_unbanded)
{
    testAlignmentDPMatrixNavigatorScoreMarixSparseInitUnbanded();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_init_banded)
{
    testAlignmentDPMatrixNavigatorScoreMarixSparseInitBanded();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_go_next)
{
    testAlignmentDPMatrixNavigatorScoreMarixSparseGoNext();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_assign_value)
{
    testAlignmentDPScoreMatrixNavigatorAssignValue(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_value)
{
    testAlignmentDPScoreMatrixNavigatorValue(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_previous_cell_diagonal)
{
    testAlignmentDPScoreMatrixNavigatorPreviousCellDiagonal(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_previous_cell_horizontal)
{
    testAlignmentDPScoreMatrixNavigatorPreviousCellHorizontal(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_score_matrix_sparse_previous_cell_vertical)
{
    testAlignmentDPScoreMatrixNavigatorPreviousCellVertical(seqan::SparseDPMatrix());
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
// Test constructor                                  [DPTraceMatrix, FullDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_constructor)
{
    testAlignmentDPMatrixTraceNavigatorConstructorDefault(seqan::TracebackOn<seqan::GapsLeft>());
    testAlignmentDPMatrixTraceNavigatorConstructorDefault(seqan::TracebackOff());
}

// ----------------------------------------------------------------------------
// Test functions                                  [DPTraceMatrix, FullDPMatrix]
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_init_unbanded)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpTraceMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix);

    _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix);
    SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
    SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_disabled_init_unbanded)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOff>, NavigateColumnWise> dpTraceMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix);

    _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    SEQAN_ASSERT_NEQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix);
    SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 0);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_init_banded)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpTraceMatrixNavigator;

    { // Case1: Band intersects with poit of origin.
        TDPMatrix dpMatrix;
        setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix);

        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOn>(-4, 3));

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
        _init(dpTraceMatrixNavigator, dpMatrix2, DPBandConfig<BandOn>(0, 7));

        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix2);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix2, Standard()), 7);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 8);
    }

    {
        TDPMatrix dpMatrix3;
        setLength(dpMatrix3, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix3, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix3);


        _init(dpTraceMatrixNavigator, dpMatrix3, DPBandConfig<BandOn>(-7, 0));

        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix3, Standard()), 0);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_disabled_init_banded)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOff>, NavigateColumnWise> dpTraceMatrixNavigator;

    { // Case1: Band intersects with poit of origin.
        TDPMatrix dpMatrix;
        setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix);

        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOn>(-4, 3));

        SEQAN_ASSERT_NEQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 0);
    }

    {
        TDPMatrix dpMatrix2;
        setLength(dpMatrix2, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix2, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix2);

        value(host(dpMatrix2), 0) = 10;
        _init(dpTraceMatrixNavigator, dpMatrix2, DPBandConfig<BandOn>(0, 7));

        SEQAN_ASSERT_NEQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix2);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 0);
    }

    {
        TDPMatrix dpMatrix3;
        setLength(dpMatrix3, DPMatrixDimension_::HORIZONTAL, 10);
        setLength(dpMatrix3, DPMatrixDimension_::VERTICAL, 8);
        resize(dpMatrix3);


        _init(dpTraceMatrixNavigator, dpMatrix3, DPBandConfig<BandOn>(-7, 0));

        SEQAN_ASSERT_NEQ(dpTraceMatrixNavigator._ptrDataContainer, &dpMatrix3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 0);
    }
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_go_next)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpTraceMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);
    resize(dpMatrix, 0);
    host(dpMatrix)[0] = 0;
    host(dpMatrix)[1] = 1;
    host(dpMatrix)[2] = 2;
    host(dpMatrix)[3] = 3;
    host(dpMatrix)[4] = 4;
    host(dpMatrix)[5] = 5;

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 0);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 1);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 2);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpTraceMatrixNavigator._activeColIterator += 2;
        ++dpTraceMatrixNavigator._laneLeap;  // For the test we simulate as if we were in a band.

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpTraceMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpTraceMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 2);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpTraceMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpTraceMatrixNavigator._activeColIterator += 2;
        ++dpTraceMatrixNavigator._laneLeap;

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpTraceMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpTraceMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 2);
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpTraceMatrixNavigator._activeColIterator += 2;

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), FirstCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 3);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), InnerCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 4);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);

        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), LastCell());
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._activeColIterator - begin(dpMatrix, Standard()), 5);
        SEQAN_ASSERT_EQ(dpTraceMatrixNavigator._laneLeap, 1);
    }
}

template <typename TNavigator, typename TMetaDescriptor, typename TCell>
void _testGoNextCell(TNavigator & navi, TMetaDescriptor const &, TCell const &)
{
    using namespace seqan;

    _goNextCell(navi, TMetaDescriptor(), TCell());
    SEQAN_ASSERT_EQ(navi._laneLeap, 0);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_disabled_go_next)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOff>, NavigateColumnWise> dpTraceMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);
    resize(dpMatrix, 0);
    host(dpMatrix)[0] = 0;
    host(dpMatrix)[1] = 1;
    host(dpMatrix)[2] = 2;
    host(dpMatrix)[3] = 3;
    host(dpMatrix)[4] = 4;
    host(dpMatrix)[5] = 5;

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, FullColumn>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        dpTraceMatrixNavigator._activeColIterator += 2;

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), LastCell());
    }

    {
        _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), FirstCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), InnerCell());
        _testGoNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, FullColumn>(), LastCell());
    }
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_assign_value)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    assignValue(dpScoreMatrixNavigator, 20);
    SEQAN_ASSERT_EQ(value(dpScoreMatrixNavigator._activeColIterator), 20);
}


SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_value)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpTraceMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    assignValue(dpTraceMatrixNavigator, 20);
    SEQAN_ASSERT_EQ(value(dpTraceMatrixNavigator), 20);

    const DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpTraceMatrixNavigatorConst(dpTraceMatrixNavigator);
    SEQAN_ASSERT_EQ(value(dpTraceMatrixNavigatorConst), 20);
}


SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_coordinate)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpScoreMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpScoreMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    dpScoreMatrixNavigator._activeColIterator += 7;

    SEQAN_ASSERT_EQ(coordinate(dpScoreMatrixNavigator, +DPMatrixDimension_::HORIZONTAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpScoreMatrixNavigator, +DPMatrixDimension_::VERTICAL), 7u);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_navigator_trace_matrix_enabled_container)
{
    using namespace seqan;

    typedef DPMatrix_<int, FullDPMatrix> TDPMatrix;

    DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpTraceMatrixNavigator;

    TDPMatrix dpMatrix;
    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 10);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 10);
    resize(dpMatrix, 3);

    _init(dpTraceMatrixNavigator, dpMatrix, DPBandConfig<BandOff>());

    SEQAN_ASSERT_EQ(&container(dpTraceMatrixNavigator), &dpMatrix);

    const DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TracebackOn<> >, NavigateColumnWise> dpTraceMatrixNavigatorConst(dpTraceMatrixNavigator);
    SEQAN_ASSERT_EQ(&container(dpTraceMatrixNavigatorConst), &dpMatrix);
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_MATRIX_NAVIGATOR_H_
