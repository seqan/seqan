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
template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
class BlockTraceNavigator
{
public:
    using TLocalTraceStore = typename TTraceProxy::TLocalTraceStore;
    using TScalarTrace     = typename TLocalTraceStore::TScalarTraceValue;
    using TSimdTrace       = typename TLocalTraceStore::TSimdTraceValue;
    using TScalarTraceMat  = typename TLocalTraceStore::TScalarTraceMatrix;
    using TSimdTraceMat    = typename TLocalTraceStore::TSimdTraceMatrix;
    using TScalarDPMatrix  = DPMatrix_<TScalarTrace, FullDPMatrix, TScalarTraceMat>;
    using TSimdDPMatrix    = DPMatrix_<TSimdTrace, FullDPMatrix, TSimdTraceMat>;
    using TScalarMatNavi   = DPMatrixNavigator_<TScalarDPMatrix, DPTraceMatrix<TTraceFlag>, NavigateColumnWise>;
    using TSimdMatNavi     = DPMatrixNavigator_<TSimdDPMatrix, DPTraceMatrix<TTraceFlag>, NavigateColumnWise>;

    enum class NavigatorSwitch : uint8_t
    {
        SCALAR,
        SIMD
    };

    TScalarMatNavi  mScalarNavigator = TScalarMatNavi{};
    TSimdMatNavi    mSimdNavigator   = TSimdMatNavi{};
    TScalarDPMatrix mScalarMatrix    = TScalarDPMatrix{};
    TSimdDPMatrix   mSimdMatrix      = TSimdDPMatrix{};
    TTraceProxy *   mHostPtr         = nullptr;
    TSeqH&          mSeqH;
    TSeqV&          mSeqV;
    size_t          mBlockH          = 0;
    size_t          mBlockV          = 0;
    size_t          mBlockSizeH      = 0;
    size_t          mBlockSizeV      = 0;
    NavigatorSwitch mSwitch          = NavigatorSwitch::SCALAR;

    template <typename TSize>
    BlockTraceNavigator(TTraceProxy &matrix,
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

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
struct Position<BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> >
{
    using Type = typename Position<TTraceProxy>::Type;
};

namespace impl
{

template <typename T>
struct LocalPosition : Position<T>
{};

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
struct LocalPosition<BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> >
{
    using Type = Pair<typename Position<BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> >::Type>;
};

}  // namespace impl

// ============================================================================
// Functions
// ============================================================================

namespace impl
{
template <typename TBlockId, typename TBlockPos,
          typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto toGlobalPosition(Pair<TBlockId, TBlockPos> const & localPos,
                             BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> const & me,
                             typename DPMatrixDimension_::TValue const & dim)
{
    if(dim == DPMatrixDimension_::HORIZONTAL)
        return localPos.i2 + localPos.i1 * me.mBlockSizeH;
    else
        return localPos.i2 + localPos.i1 * me.mBlockSizeV;
}

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TNavigator,
          typename TDimension>
inline auto
coordinate(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> const & dpNavigator,
           TNavigator const & navi,
           TDimension const dimension)
{
    using TPos = typename impl::LocalPosition<BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> >::Type;

    if (dimension == DPMatrixDimension_::HORIZONTAL)
        return impl::toGlobalPosition(TPos(dpNavigator.mBlockH, coordinate(navi, dimension)), dpNavigator,
                                      +DPMatrixDimension_::HORIZONTAL);
    else
        return impl::toGlobalPosition(TPos(dpNavigator.mBlockV, coordinate(navi, dimension)), dpNavigator,
                                      +DPMatrixDimension_::VERTICAL);
}

}  // namespace impl

// ----------------------------------------------------------------------------
// Function _initNextBlock
// ----------------------------------------------------------------------------

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TBlockId,
          typename TBlockPos>
inline void
_initNextScalarBlock(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & me,
                     TBlockId const & blockId,
                     TBlockPos const blockPos)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    me.mSwitch = TNavi::NavigatorSwitch::SCALAR;

    setLength(me.mScalarMatrix, +DPMatrixDimension_::HORIZONTAL, length(me.mSeqH[me.mBlockH]) + 1);
    setLength(me.mScalarMatrix, +DPMatrixDimension_::VERTICAL, length(me.mSeqV[me.mBlockV]) + 1);
    setHost(me.mScalarMatrix, std::get<0>(blockId)->mScalarTraceVec[std::get<1>(blockId).mBlockPos]);
    resize(me.mScalarMatrix);  // Needed to actually update the internal variables.

    me.mScalarNavigator._ptrDataContainer = &me.mScalarMatrix;
    me.mScalarNavigator._activeColIterator = begin(me.mScalarMatrix, Standard());
    _setToPosition(me.mScalarNavigator, blockPos);  // Set the navigator to the correct position within the block.
}

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TBlockId,
          typename TBlockPos>
inline void
_initNextSimdBlock(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & me,
                   TBlockId const & blockId,
                   TBlockPos const blockPos)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    me.mSwitch = TNavi::NavigatorSwitch::SIMD;

    setLength(me.mSimdMatrix, +DPMatrixDimension_::HORIZONTAL, length(me.mSeqH[me.mBlockH]) + 1);
    setLength(me.mSimdMatrix, +DPMatrixDimension_::VERTICAL, length(me.mSeqV[me.mBlockV]) + 1);
    setHost(me.mSimdMatrix, std::get<0>(blockId)->mSimdTraceVec[std::get<1>(blockId).mBlockPos]);
    resize(me.mSimdMatrix);  // Needed to actually update the internal variables.

    me.mSimdNavigator._ptrDataContainer = &me.mSimdMatrix;
    me.mSimdNavigator._activeColIterator = begin(me.mSimdMatrix, Standard());
    _setSimdLane(me.mSimdNavigator, std::get<2>(blockId));
    _setToPosition(me.mSimdNavigator, blockPos);  // Set the navigator to the correct position within the block.
}

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TBlockPos>
inline void
_initNextBlock(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & me,
               TBlockPos const blockPos)
{
    auto blockId = me.mHostPtr->mTraceBlockMap[me.mBlockH * length(me.mSeqV) + me.mBlockV];

    if (std::get<1>(blockId).mBlockId == 0)
        _initNextScalarBlock(me, blockId, blockPos);
    else
        _initNextSimdBlock(me, blockId, blockPos);
}

// ----------------------------------------------------------------------------
// Function _traceHorizontal()
// ----------------------------------------------------------------------------

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TInternalNavigator>
inline void
_traceHorizontal(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & me,
                 TInternalNavigator & navi,
                 bool const isBandShift)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    // Update the block.
    if (coordinate(navi, +DPMatrixDimension_::HORIZONTAL) == 0)
    {
        --me.mBlockH;
        auto blockPos = length(me.mSeqH[me.mBlockH]) * (length(me.mSeqV[me.mBlockV]) + 1) +
                        coordinate(navi, +DPMatrixDimension_::VERTICAL);
        _initNextBlock(me, blockPos);
    }

    if (me.mSwitch == TNavi::NavigatorSwitch::SCALAR)
        _traceHorizontal(me.mScalarNavigator, isBandShift);
    else
        _traceHorizontal(me.mSimdNavigator, isBandShift);
}

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline void
_traceHorizontal(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & me,
                 bool const isBandShift)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    if (me.mSwitch == TNavi::NavigatorSwitch::SCALAR)
        _traceHorizontal(me, me.mScalarNavigator, isBandShift);
    else
        _traceHorizontal(me, me.mSimdNavigator, isBandShift);
}

// ----------------------------------------------------------------------------
// Function _traceDiagonal()
// ----------------------------------------------------------------------------

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TInternalNavigator>
inline void
_traceDiagonal(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & me,
               TInternalNavigator & navi,
               bool const isBandShift)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    if (coordinate(navi, +DPMatrixDimension_::VERTICAL) == 0 &&
        coordinate(navi, +DPMatrixDimension_::HORIZONTAL) == 0)
    {  // Continue in previous diagonal block.
        --me.mBlockV;
        --me.mBlockH;
        auto blockPos = length(me.mSeqH[me.mBlockH]) * (length(me.mSeqV[me.mBlockV]) + 1) +
                        length(me.mSeqV[me.mBlockV]);
        _initNextBlock(me, blockPos);
    }
    else if(coordinate(navi, +DPMatrixDimension_::VERTICAL) == 0)
    {  // Continue in previous vertical block.
        --me.mBlockV;
        auto blockPos = coordinate(navi, +DPMatrixDimension_::HORIZONTAL) *
                        (length(me.mSeqV[me.mBlockV]) + 1) + length(me.mSeqV[me.mBlockV]);
        _initNextBlock(me, blockPos);
    }
    else if (coordinate(navi, +DPMatrixDimension_::HORIZONTAL) == 0)
    {  // Continue in previous horizontal block.
        --me.mBlockH;
        auto blockPos = length(me.mSeqH[me.mBlockH]) * (length(me.mSeqV[me.mBlockV]) + 1) +
                        coordinate(navi, +DPMatrixDimension_::VERTICAL);
        _initNextBlock(me, blockPos);
    }

    if (me.mSwitch == TNavi::NavigatorSwitch::SCALAR)
        _traceDiagonal(me.mScalarNavigator, isBandShift);
    else
        _traceDiagonal(me.mSimdNavigator, isBandShift);
}

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline void
_traceDiagonal(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & me,
               bool const isBandShift)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    if (me.mSwitch == TNavi::NavigatorSwitch::SCALAR)
        _traceDiagonal(me, me.mScalarNavigator, isBandShift);
    else
        _traceDiagonal(me, me.mSimdNavigator, isBandShift);
}

// ----------------------------------------------------------------------------
// Function _traceVertical()
// ----------------------------------------------------------------------------

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TInternalNavigator>
inline void
_traceVertical(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & me,
               TInternalNavigator & navi,
               bool const isBandShift)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    if (coordinate(navi, +DPMatrixDimension_::VERTICAL) == 0)
    {
        --me.mBlockV;
        auto blockPos = coordinate(navi, +DPMatrixDimension_::HORIZONTAL) *
                        (length(me.mSeqV[me.mBlockV]) + 1) + length(me.mSeqV[me.mBlockV]);
        _initNextBlock(me, blockPos);
    }

    if (me.mSwitch == TNavi::NavigatorSwitch::SCALAR)
        _traceVertical(me.mScalarNavigator, isBandShift);
    else
        _traceVertical(me.mSimdNavigator, isBandShift);
}

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline void
_traceVertical(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & me,
               bool const isBandShift)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return;  // Do nothing since no trace back is computed.

    if (me.mSwitch == TNavi::NavigatorSwitch::SCALAR)
        _traceVertical(me, me.mScalarNavigator, isBandShift);
    else
        _traceVertical(me, me.mSimdNavigator, isBandShift);
}

// ----------------------------------------------------------------------------
// Function scalarValue(); Helper to switch between simd and scalar.
// ----------------------------------------------------------------------------

// Wrapper to get the scalar value.
template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto
scalarValue(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> const & dpNavigator)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");

    if (dpNavigator.mSwitch == TNavi::NavigatorSwitch::SCALAR)
        return scalarValue(dpNavigator.mScalarNavigator);

    return scalarValue(dpNavigator.mSimdNavigator);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

// Current position.
// NOTE(rrahn): This cannot work, since the return value depends on the internal state of the dpNavigator.
//template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
//inline auto&
//value(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & dpNavigator)
//{
//    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;
//
//    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
//        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");
//
//    if (dpNavigator.mSwitch == TNavi::NavigatorSwitch::SCALAR)
//        return value(dpNavigator.mScalarNavigator);
//
//    return value(dpNavigator.mSimdNavigator);
//}
//
//template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
//inline auto&
//value(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> const & dpNavigator)
//{
//    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;
//
//    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
//        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");
//
//    if (dpNavigator.mSwitch == TNavi::NavigatorSwitch::SCALAR)
//        return value(dpNavigator.mScalarNavigator);
//
//    return value(dpNavigator.mSimdNavigator);
//}

// ----------------------------------------------------------------------------
// Function coordinate()
// ----------------------------------------------------------------------------

// Returns the coordinate of the given dimension for the current position of the
// navigator within the matrix.
template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto
coordinate(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> const & dpNavigator,
           typename DPMatrixDimension_::TValue const & dimension)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        SEQAN_ASSERT_FAIL("Try to access uninitialized object!");
    SEQAN_ASSERT_EQ(_checkCorrectDimension(dimension), true);

    if (dpNavigator.mSwitch == TNavi::NavigatorSwitch::SCALAR)
        return impl::coordinate(dpNavigator, dpNavigator.mScalarNavigator, dimension);

    return impl::coordinate(dpNavigator, dpNavigator.mSimdNavigator, dimension);
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

// Returns the current position of the navigator within the matrix.
template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag>
inline auto
position(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> const & dpNavigator)
{
    using TPos = typename Position<BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> >::Type;
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;
    // Return 0 when traceback is not enabled. This is necessary to still track the score even
    // the traceback is not enabled.
    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return TPos(0, 0);

    if (dpNavigator.mSwitch == TNavi::NavigatorSwitch::SCALAR)
        return TPos(dpNavigator.mBlockH * length(dpNavigator.mSeqV) + dpNavigator.mBlockV,
                    position(dpNavigator.mScalarNavigator));

    return TPos(dpNavigator.mBlockH * length(dpNavigator.mSeqV) + dpNavigator.mBlockV,
                position(dpNavigator.mSimdNavigator));
}

// ----------------------------------------------------------------------------
// Function _setSimdLane()
// ----------------------------------------------------------------------------

template <typename TTraceProxy, typename TSeqH, typename TSeqV, typename TTraceFlag,
          typename TPos>
inline auto&
_setSimdLane(BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag> & dpNavigator,
             TPos const pos)
{
    using TNavi = BlockTraceNavigator<TTraceProxy, TSeqH, TSeqV, TTraceFlag>;

    if (dpNavigator.mSwitch == TNavi::NavigatorSwitch::SIMD)
        _setSimdLane(dpNavigator.mSimdNavigator, pos);
}
    
}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_PARALLEL_DP_TRACE_MATRIX_NAVIGATOR_BLOCK_WISE_H_
