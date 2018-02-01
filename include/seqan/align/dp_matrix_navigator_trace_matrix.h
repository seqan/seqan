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
// This specialization is used to navigate through the traceback matrix
// of any standard dp-alignment algorithm. The DPTraceMatrix gets the
// traceback flag TracebackOn or TracebackOff. A traceback is only computed
// if the traceback is switched on. If this is not the case, the void
// functions will be compiled as no-op functions, while in functions that try
// to access a value of the underlying matrix via the navigator an assertion
// is thrown.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_TRACE_MATRIX_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_TRACE_MATRIX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPMatrixNavigator                        [FullDPMatrix, DPTraceMatrix]
// ----------------------------------------------------------------------------

// The matrix navigator for the trace-back matrix.
//
// It takes three types to be specialized. The first type defines the underlying
// dp-matrix it is working on. This has to be a FullDPMatrix. The second type,
// specifies that this is a trace-matrix navigator while the TTraceFlag can either
// be TracebackOn to enable the navigator or TracebackOff to disable it.
// The last parameter specifies the kind of navigation.
template <typename TValue, typename THost, typename TTraceFlag, typename TNavigationSpec>
class DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>, DPTraceMatrix<TTraceFlag>, TNavigationSpec>
{
public:

    typedef  DPMatrix_<TValue, FullDPMatrix, THost> TDPMatrix_;
    typedef typename Pointer_<TDPMatrix_>::Type TDPMatrixPointer_;
    typedef typename Iterator<TDPMatrix_, Standard>::Type TDPMatrixIterator;

    template <typename TBandSpec,
              std::enable_if_t<std::is_same<TBandSpec, BandOff>::value, int> = 0>
    DPMatrixNavigator_(TDPMatrix_ & matrix,
                       DPBandConfig<TBandSpec> const & /*band*/)
    {
        if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
            return;  // Leave navigator uninitialized because it is never used.

        _ptrDataContainer = &matrix;
        _activeColIterator = begin(matrix, Standard());
        _laneLeap = 1;
        *_activeColIterator = TValue();
    }

    template <typename TBandSpec,
              std::enable_if_t<!std::is_same<TBandSpec, BandOff>::value, int> = 0>
    DPMatrixNavigator_(TDPMatrix_ & matrix,
                       DPBandConfig<TBandSpec> const & band)
    {
        using TMatrixSize = typename Size<TDPMatrix_>::Type;
        using TSignedSize = std::make_signed_t<TMatrixSize>;

        if (std::is_same<TTraceFlag, TracebackOff>::value)
            return;  // Leave navigator as is because it should never be used.

        _ptrDataContainer = &matrix;

        // Band begins within the first row.
        if (lowerDiagonal(band) >= 0)
        {
            // The first cell of the first column starts at the last cell in the matrix of the current column.
            _laneLeap = _min(length(matrix, DPMatrixDimension_::VERTICAL), bandSize(band));
            _activeColIterator = begin(matrix, Standard()) + _dataLengths(matrix)[DPMatrixDimension_::VERTICAL] - 1;
        }
        else if (upperDiagonal(band) <= 0)  // Band begins within the first column.
        {
            // The first cell starts at the beginning of the current column.
            _laneLeap = 1;
            _activeColIterator = begin(matrix, Standard());
        }
        else  // Band intersects with the point of origin.
        {
            // First cell starts at position i, such that i + abs(lowerDiagonal) = length(seqV).
            TMatrixSize lengthVertical = length(matrix, DPMatrixDimension_::VERTICAL);
            int lastPos = _max(-static_cast<TSignedSize>(lengthVertical - 1), lowerDiagonal(band));
            _laneLeap = lengthVertical + lastPos;
            _activeColIterator = begin(matrix, Standard()) + _laneLeap - 1;
        }
        *_activeColIterator = TValue();
    }

    TDPMatrixPointer_   _ptrDataContainer{nullptr}; // The pointer to the underlying Matrix.
    int                 _laneLeap{0};               // Keeps track of the jump size from one column to another.
    unsigned            _simdLane{0};               // Used for tracing the correct cell in case of simd vectors.
    TDPMatrixIterator   _activeColIterator{};       // The current column iterator.
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _goNextCell()                                   [banded, FirstCell]
// ----------------------------------------------------------------------------

// In the initial column we don't need to do anything because, the navigagtor is already initialized.
template <typename TValue, typename THost, typename TTraceFlag,
          typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWiseBanded> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            FirstCell const &)
{
    // no-op
}

// overload needed to avoid ambigious overload.
template <typename TValue, typename THost, typename TTraceFlag>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWiseBanded> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, PartialColumnTop> const &,
            FirstCell const &)
{
    // no-op
}

// We are in the banded case, where the band crosses the first row.
// The left cell of the active cell is not valid, beacause we only can come from horizontal direction.
// The lower left cell of the active cell is the horizontal direction.
template <typename TValue, typename THost, typename TTraceFlag,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, PartialColumnTop> const &,
            FirstCell const &)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    --dpNavigator._laneLeap;
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
}

// We are in the banded case.
// The left cell of the active cell represents diagonal direction. The lower left diagonal represents the horizontal direction.
template <typename TValue, typename THost, typename TTraceFlag,
          typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            FirstCell const &)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    // go to begin of next column.
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
}

// ----------------------------------------------------------------------------
// Function _goNextCell()                                 [unbanded, FirstCell]
// ----------------------------------------------------------------------------

//
template <typename TValue, typename THost, typename TTraceFlag>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWise> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            FirstCell const &)
{
    // no-op
}

template <typename TValue, typename THost, typename TTraceFlag,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            FirstCell const &)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    dpNavigator._activeColIterator += dpNavigator._laneLeap;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                     [banded, InnerCell]
// ----------------------------------------------------------------------------

// For any other column type and location we can use the same navigation procedure.
template <typename TValue, typename THost, typename TTraceFlag,
          typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            InnerCell const &)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                   [unbanded, InnerCell]
// ----------------------------------------------------------------------------

// For any other column type and location we can use the same navigation procedure.
template <typename TValue, typename THost, typename TTraceFlag,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            InnerCell const &)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                      [banded, LastCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename THost, typename TTraceFlag>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom> const &,
            LastCell const &)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    ++dpNavigator._activeColIterator;
}

// If we are in banded case and the band crosses the last row, we have to update
// the additional leap for the current track.
template <typename TValue, typename THost, typename TTraceFlag,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, PartialColumnBottom> const &,
            LastCell const &)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    ++dpNavigator._activeColIterator;
    ++dpNavigator._laneLeap;
}

// If we are in the banded case the left cell of the active represents the diagonal direction.
template <typename TValue, typename THost, typename TTraceFlag,
          typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWiseBanded> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            LastCell const &)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                    [unbanded, LastCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename THost, typename TTraceFlag,
          typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            LastCell const &)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _traceHorizontal()
// ----------------------------------------------------------------------------

template <typename TValue, typename THost, typename TTraceFlag, typename TNavigationSpec>
inline void
_traceHorizontal(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  TNavigationSpec> & dpNavigator,
                 bool isBandShift)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    if (isBandShift)
        dpNavigator._activeColIterator -= _dataFactors(*dpNavigator._ptrDataContainer)[DPMatrixDimension_::HORIZONTAL] - 1;
    else
        dpNavigator._activeColIterator -= _dataFactors(*dpNavigator._ptrDataContainer)[DPMatrixDimension_::HORIZONTAL];

}

// ----------------------------------------------------------------------------
// Function _traceDiagonal()
// ----------------------------------------------------------------------------

template <typename TValue, typename THost, typename TTraceFlag, typename TNavigationSpec>
inline void
_traceDiagonal(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  TNavigationSpec> & dpNavigator,
               bool isBandShift)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    if (isBandShift)
        dpNavigator._activeColIterator -= _dataFactors(*dpNavigator._ptrDataContainer)[DPMatrixDimension_::HORIZONTAL];
    else
        dpNavigator._activeColIterator -= _dataFactors(*dpNavigator._ptrDataContainer)[DPMatrixDimension_::HORIZONTAL] + 1;

}

// ----------------------------------------------------------------------------
// Function _traceVertical()
// ----------------------------------------------------------------------------

template <typename TValue, typename THost, typename TTraceFlag, typename TNavigationSpec>
inline void
_traceVertical(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  TNavigationSpec> & dpNavigator,
               bool /*isBandShift*/)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    dpNavigator._activeColIterator -= _dataFactors(*dpNavigator._ptrDataContainer)[DPMatrixDimension_::VERTICAL];
}

// ----------------------------------------------------------------------------
// Function setToPosition()
// ----------------------------------------------------------------------------

template <typename TValue, typename THost, typename TTraceFlag, typename TNavigationSpec,
          typename TPosition>
inline void
_setToPosition(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  TNavigationSpec> & dpNavigator,
              TPosition const & hostPosition)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;

    SEQAN_ASSERT_LT(hostPosition, static_cast<TPosition>(length(container(dpNavigator))));

    dpNavigator._activeColIterator = begin(*dpNavigator._ptrDataContainer, Standard()) + hostPosition;
}


// Sets the host position based on the given horizontal and vertical position. Note that the horizontal and
// vertical positions must correspond to the correct size of the underlying matrix.
// For banded matrices the vertical dimension might not equal the length of the vertical sequence.
template <typename TValue, typename THost, typename TTraceFlag, typename TNavigationSpec,
          typename TPositionH, typename TPositionV>
inline void
_setToPosition(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix, THost>,
                                  DPTraceMatrix<TTraceFlag>,
                                  TNavigationSpec> & dpNavigator,
              TPositionH const & horizontalPosition,
              TPositionV const & verticalPosition)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;
    SEQAN_ASSERT_LT(horizontalPosition, static_cast<TPositionH>(length(container(dpNavigator), +DPMatrixDimension_::HORIZONTAL)));
    SEQAN_ASSERT_LT(verticalPosition, static_cast<TPositionV>(length(container(dpNavigator), +DPMatrixDimension_::VERTICAL)));

    TPositionH  hostPosition = horizontalPosition * _dataFactors(container(dpNavigator))[+DPMatrixDimension_::HORIZONTAL] + verticalPosition;
    dpNavigator._activeColIterator = begin(*dpNavigator._ptrDataContainer, Standard()) + hostPosition;
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec,
          typename TValue>
inline void
assignValue(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> & dpNavigator,
            TValue const & element)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    *dpNavigator._activeColIterator = element;
}

// ----------------------------------------------------------------------------
// Function scalarValue(); Helper to switch between simd and scalar.
// ----------------------------------------------------------------------------

// SIMD Version. Returns always a copy and never a reference.
template <typename TValue, typename TPos>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TValue> >, typename TraceBitMap_<>::Type)
_scalarValue(TValue const & vec,
             TPos const pos)
{
    return value(vec, pos);
}

// Non-simd variant. Identity version.
template <typename TValue, typename TPos>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TValue> > >, TValue)
_scalarValue(TValue & val,
             TPos const /*pos*/)
{
    return val;
}

// Wrapper to get the scalar value.
template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec>
inline auto
scalarValue(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> const & dpNavigator)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");

    return _scalarValue(*dpNavigator._activeColIterator, dpNavigator._simdLane);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

// Current position.
template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> >::Type
value(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> & dpNavigator)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");

    return *dpNavigator._activeColIterator;
}

template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> const >::Type
value(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> const & dpNavigator)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");

    return *dpNavigator._activeColIterator;
}

// At specified position.
template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec, typename TPosition>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> >::Type
value(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> & dpNavigator,
      TPosition const & position)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");

    return *(begin(*dpNavigator._ptrDataContainer) + position);
}

template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec,
          typename TPosition>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> const >::Type
value(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> const & dpNavigator,
      TPosition const & position)
{
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");

    return *(begin(*dpNavigator._ptrDataContainer) + position);
}

// ----------------------------------------------------------------------------
// Function coordinate()
// ----------------------------------------------------------------------------

// Returns the coordinate of the given dimension for the current position of the
// navigator within the matrix.
template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec>
inline size_t
coordinate(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> const & dpNavigator,
           typename DPMatrixDimension_::TValue const & dimension)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return 0;  // Returns default 0, when traceback is set off.
    SEQAN_ASSERT_EQ(_checkCorrectDimension(dimension), true);

    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return _dataLengths(*dpNavigator._ptrDataContainer)[dimension];  // Return lengths of given dimension.

    return coordinate(value(dpNavigator._ptrDataContainer), position(dpNavigator), dimension); // Simply delegate to coordinate of underlying matrix.
}

// ----------------------------------------------------------------------------
// Function toGlobalPosition()
// ----------------------------------------------------------------------------

// Returns the current position of the navigator within the matrix.
template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec,
          typename TPosH,
          typename TPosV>
inline typename Position<DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> >::Type
toGlobalPosition(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> const & dpNavigator,
                 TPosH const horizontalCoordinate,
                 TPosV const verticalCoordinate)
{
    // Return 0 when traceback is not enabled. This is necessary to still track the score even
    // the traceback is not enabled.
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return 0;

    return  toGlobalPosition(*dpNavigator._ptrDataContainer, horizontalCoordinate, verticalCoordinate);
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

// Returns the current position of the navigator within the matrix.
template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec>
inline typename Position<DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> >::Type
position(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> const & dpNavigator)
{
    // Return 0 when traceback is not enabled. This is necessary to still track the score even
    // the traceback is not enabled.
    SEQAN_IF_CONSTEXPR (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return 0;

    return position(dpNavigator._activeColIterator, *dpNavigator._ptrDataContainer);
}

// ----------------------------------------------------------------------------
// Function _setSimdLane()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TTraceFlag, typename TNavigationSpec,
          typename TPos>
inline void
_setSimdLane(DPMatrixNavigator_<TDPMatrix, DPTraceMatrix<TTraceFlag>, TNavigationSpec> & dpNavigator,
             TPos const pos)
{
    SEQAN_ASSERT_LT(pos, static_cast<TPos>(LENGTH<typename Value<TDPMatrix>::Type>::VALUE));

    dpNavigator._simdLane = pos;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_TRACE_MATRIX_H_
