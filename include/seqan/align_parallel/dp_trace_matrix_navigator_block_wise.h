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
// This specialization is used to navigate through the traceback matrix
// of any standard dp-alignment algorithm. The DPTraceMatrix gets the
// traceback flag TracebackOn or TracebackOff. A traceback is only computed
// if the traceback is switched on. If this is not the case, the void
// functions will be compiled as no-op functions, while in functions that try
// to access a value of the underlying matrix via the navigator an assertion
// is thrown.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TRACE_MATRIX_NAVIGATOR_BLOCK_WISE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TRACE_MATRIX_NAVIGATOR_BLOCK_WISE_H_

namespace seqan
{

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
template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
class BlockTraceNavigator
{
public:
    using TBlock       = typename Value<TMatrix>::Type;
    using TTraceValue  = typename Value<TBlock>::Type;
    using TBlockMatrix = DPMatrix_<TTraceValue, FullDPMatrix, TBlock>;
    using TNavigator   = DPMatrixNavigator_<TBlockMatrix, DPTraceMatrix<TTraceFlag>, NavigateColumnWise>;
    using TMatrixPtr_  = typename Pointer_<TMatrix>::Type;

    TNavigator                    mBlockNavigator  = TNavigator{};
    TBlockMatrix                  mBlockMatrix     = TBlockMatrix{};
    TMatrixPtr_                   mHostPtr         = nullptr;
    TSeqH&                        mSeqH;
    TSeqV&                        mSeqV;
    typename Size<TMatrix>::Type  mBlockH          = 0;
    typename Size<TMatrix>::Type  mBlockV          = 0;
    typename Size<TMatrix>::Type  mBlockSizeH      = 0;
    typename Size<TMatrix>::Type  mBlockSizeV      = 0;

    template <typename TSize>
    BlockTraceNavigator(TMatrix &matrix,
                        TSeqH &seqH,
                        TSeqV &seqV,
                        TSize const blockH,
                        TSize const blockV,
                        TSize const blockPos) :
        mHostPtr(&matrix),
        mSeqH(seqH),
        mSeqV(seqV),
        mBlockH(blockH),
        mBlockV(blockV),
        mBlockSizeH(length(front(mSeqH))),
        mBlockSizeV(length(front(mSeqV)))
    {
        _initNextBlock(*this, blockPos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
struct Position<BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> >
{
    using Type = typename Position<TMatrix>::Type;
};

namespace impl
{

template <typename T>
struct LocalPosition : Position<T>
{};

template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
struct LocalPosition<BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> >
{
    using Type = Pair<typename Position<BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> >::Type>;
};

}  // namespace impl

// ============================================================================
// Functions
// ============================================================================

namespace impl
{
template <typename TBlockId, typename TBlockPos,
          typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto toGlobalPosition(Pair<TBlockId, TBlockPos> const & localPos,
                             BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> const & me,
                             typename DPMatrixDimension_::TValue const & dim)
{
    if(dim == DPMatrixDimension_::HORIZONTAL)
        return localPos.i2 + localPos.i1 * me.mBlockSizeH;
    else
        return localPos.i2 + localPos.i1 * me.mBlockSizeV;
}
}  // namespace impl

// ----------------------------------------------------------------------------
// Function _goNextCell
// ----------------------------------------------------------------------------

template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TBlockPos>
inline void
_initNextBlock(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> & me,
               TBlockPos const blockPos)
{
    setLength(me.mBlockMatrix, +DPMatrixDimension_::HORIZONTAL, length(me.mSeqH[me.mBlockH]) + 1);
    setLength(me.mBlockMatrix, +DPMatrixDimension_::VERTICAL, length(me.mSeqV[me.mBlockV]) + 1);
    setHost(me.mBlockMatrix, (*me.mHostPtr)[me.mBlockH * length(me.mSeqV) + me.mBlockV]);
    resize(me.mBlockMatrix);  // Needed to actually update the internal variables.

    me.mBlockNavigator._ptrDataContainer = &me.mBlockMatrix;
    me.mBlockNavigator._activeColIterator = begin(me.mBlockMatrix, Standard());
    _setToPosition(me.mBlockNavigator, blockPos);  // Set the navigator to the correct position within the block.
}

// ----------------------------------------------------------------------------
// Function _traceHorizontal()
// ----------------------------------------------------------------------------

template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline void
_traceHorizontal(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> & me,
                 bool isBandShift)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    // Update the block.
    if (coordinate(me.mBlockNavigator, +DPMatrixDimension_::HORIZONTAL) == 0)
    {
        --me.mBlockH;
        auto blockPos = (length(me.mSeqH[me.mBlockH]) - 1) * (length(me.mSeqV[me.mBlockV]) + 1) +
                        coordinate(me.mBlockNavigator, +DPMatrixDimension_::VERTICAL);
        _initNextBlock(me, blockPos);
        // Block pos is last column times length of the column + coordinate in vertical dimension.
        return;
    }
    _traceHorizontal(me.mBlockNavigator, isBandShift);
}

// ----------------------------------------------------------------------------
// Function _traceDiagonal()
// ----------------------------------------------------------------------------

template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline void
_traceDiagonal(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> & me,
               bool isBandShift)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    if (coordinate(me.mBlockNavigator, +DPMatrixDimension_::VERTICAL) == 0 &&
        coordinate(me.mBlockNavigator, +DPMatrixDimension_::HORIZONTAL) == 0)
    {  // Continue in previous diagonal block.
        --me.mBlockV;
        --me.mBlockH;
        auto blockPos = (length(me.mSeqH[me.mBlockH]) - 1) * (length(me.mSeqV[me.mBlockV]) + 1) +
                        length(me.mSeqV[me.mBlockV]) - 1;
        _initNextBlock(me, blockPos);
    }
    else if(coordinate(me.mBlockNavigator, +DPMatrixDimension_::VERTICAL) == 0)
    {  // Continue in previous vertical block.
        --me.mBlockV;
        auto blockPos = (coordinate(me.mBlockNavigator, +DPMatrixDimension_::HORIZONTAL) - 1) *
                        (length(me.mSeqV[me.mBlockV]) + 1) + length(me.mSeqV[me.mBlockV]) - 1;
        _initNextBlock(me, blockPos);
    }
    else if (coordinate(me.mBlockNavigator, +DPMatrixDimension_::HORIZONTAL) == 0)
    {  // Continue in previous horizontal block.
        --me.mBlockH;
        auto blockPos = (length(me.mSeqH[me.mBlockH]) - 1) * (length(me.mSeqV[me.mBlockV]) + 1) +
                        coordinate(me.mBlockNavigator, +DPMatrixDimension_::VERTICAL) - 1;
        _initNextBlock(me, blockPos);
    }
    else
    {
        _traceDiagonal(me.mBlockNavigator, isBandShift);
    }
}

// ----------------------------------------------------------------------------
// Function _traceVertical()
// ----------------------------------------------------------------------------

template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline void
_traceVertical(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> & me,
               bool isBandShift)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    if (coordinate(me.mBlockNavigator, +DPMatrixDimension_::VERTICAL) == 0)
    {
        --me.mBlockV;
        auto blockPos = coordinate(me.mBlockNavigator, +DPMatrixDimension_::HORIZONTAL) *
                        (length(me.mSeqV[me.mBlockV]) + 1) + length(me.mSeqV[me.mBlockV]) - 1;
        _initNextBlock(me, blockPos);
    }
    else
    {
        _traceVertical(me.mBlockNavigator, isBandShift);
    }
}

// ----------------------------------------------------------------------------
// Function scalarValue(); Helper to switch between simd and scalar.
// ----------------------------------------------------------------------------

// Wrapper to get the scalar value.
template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto
scalarValue(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> const & dpNavigator)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");

    return scalarValue(dpNavigator.mBlockNavigator);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

// Current position.
template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto&
value(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> & dpNavigator)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");

    return value(dpNavigator.mBlockNavigator);
}

template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto&
value(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> const & dpNavigator)
{
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");

    return value(dpNavigator.mBlockNavigator);
}

// ----------------------------------------------------------------------------
// Function coordinate()
// ----------------------------------------------------------------------------

// Returns the coordinate of the given dimension for the current position of the
// navigator within the matrix.
template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto
coordinate(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> const & dpNavigator,
           typename DPMatrixDimension_::TValue const & dimension)
{
    using TPos = typename impl::LocalPosition<BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> >::Type;

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");
    SEQAN_ASSERT_EQ(_checkCorrectDimension(dimension), true);

    if (dimension == DPMatrixDimension_::HORIZONTAL)
        return impl::toGlobalPosition(TPos(dpNavigator.mBlockH, coordinate(dpNavigator.mBlockNavigator, dimension)),
                                      dpNavigator, +DPMatrixDimension_::HORIZONTAL);
    else
        return impl::toGlobalPosition(TPos(dpNavigator.mBlockV, coordinate(dpNavigator.mBlockNavigator, dimension)),
                                      dpNavigator, +DPMatrixDimension_::VERTICAL);
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

// Returns the current position of the navigator within the matrix.
template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto
position(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> const & dpNavigator)
{
    using TPos = typename Position<BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> >::Type;
    // Return 0 when traceback is not enabled. This is necessary to still track the score even
    // the traceback is not enabled.
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return TPos(0, 0);

    return TPos(dpNavigator.mBlockH * length(dpNavigator.mSeqV) + dpNavigator.mBlockV,
                          position(dpNavigator.mBlockNavigator));
}

// ----------------------------------------------------------------------------
// Function _setSimdLane()
// ----------------------------------------------------------------------------

template <typename TMatrix, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TPos>
inline auto&
_setSimdLane(BlockTraceNavigator<TMatrix, TSeqH, TSeqV, TTraceFlag> & dpNavigator,
             TPos const pos)
{
    _setSimdLane(dpNavigator.mBlockNavigator, pos);
}
    
}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TRACE_MATRIX_NAVIGATOR_BLOCK_WISE_H_
