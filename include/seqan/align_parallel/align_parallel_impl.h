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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_PARALLEL_IMPL_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_PARALLEL_IMPL_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TCond1, typename TCond2>
struct CorrectLastColumn_ : False
{};

template <>
struct CorrectLastColumn_<True, True> : True
{};

template <typename TCond1, typename TCond2>
struct CorrectLastRow_ : False
{};

template <>
struct CorrectLastRow_<True, True> : True
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _computeCell(); InitialCol;
// ----------------------------------------------------------------------------

// The _computeCell function is the basic interface that is called to comute
// the score for each cell and to store the corresponding traceback.
// The MetaColumnDescriptor and the CellDescriptor describe which cell in the dp matrix
// is computed. We use this information to overload the functions in order
// to initialize from the passed buffer and to store the last row/column in the buffer.

// Vertical initialization values are copied from buffer.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TCellDescriptor,
          typename TAlgo, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const & /*previousDiagonal*/,
             DPCell_<TScoreValue, TGapCosts> const & /*previousHorizontal*/,
             DPCell_<TScoreValue, TGapCosts> const & /*previousVertical*/,
             TSequenceHValue const & /*seqHVal*/,
             TSequenceVValue const & /*seqVVal*/,
             TScoringScheme const & /*scoringScheme*/,
             MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
             TCellDescriptor const &,   // One of FirstCell, InnerCell or LastCell.
             DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel> const &)
{
//    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
//    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInitialColumn, FullColumn> >   TMetaColumn;

//    auto pos = coordinate(traceMatrixNavigator, +DPMatrixDimension_::VERTICAL);
    activeCell = (*scout.state.ptrVerBuffer)[scout.mVerticalPos].i1;
    assignValue(traceMatrixNavigator, (*scout.state.ptrVerBuffer)[scout.mVerticalPos].i2);

//    if (TrackingEnabled_<TMetaColumn, TCellDescriptor>::VALUE)
//    {
//        typedef typename And<IsSameType<TCellDescriptor, LastCell>,
//                             Or<IsSameType<typename MetaColumnDescriptor<DPInitialColumn, FullColumn>::TLocation, PartialColumnBottom>,
//                                IsSameType<typename MetaColumnDescriptor<DPInitialColumn, FullColumn>::TLocation, FullColumn>
//                             >
//                         >::Type TIsLastRow;
//        _scoutBestScore(scout, activeCell, traceMatrixNavigator, False(), TIsLastRow());
//    }
}

// ----------------------------------------------------------------------------
// Function _computeCell(); InnerCol; FirstCell
// ----------------------------------------------------------------------------

// Horizontal initialization values are copied from buffer for all first cells.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const & /*previousDiagonal*/,
             DPCell_<TScoreValue, TGapCosts> const & /*previousHorizontal*/,
             DPCell_<TScoreValue, TGapCosts> const & /*previousVertical*/,
             TSequenceHValue const & /*seqHVal*/,
             TSequenceVValue const & /*seqVVal*/,
             TScoringScheme const & /*scoringScheme*/,
             MetaColumnDescriptor<DPInnerColumn, FullColumn> const &,
             FirstCell const &,   // One of FirstCell, InnerCell or LastCell.
             DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel> const &)
{
//    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
//    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, FullColumn> >     TMetaColumn;

//    auto pos = coordinate(traceMatrixNavigator, +DPMatrixDimension_::HORIZONTAL) - 1;
    activeCell = (*scout.state.ptrHorBuffer)[scout.mHorizontalPos - 1].i1;
    assignValue(traceMatrixNavigator, (*scout.state.ptrHorBuffer)[scout.mHorizontalPos - 1].i2);

//    if (TrackingEnabled_<TMetaColumn, FirstCell>::VALUE)
//    {
//        _scoutBestScore(scout, activeCell, traceMatrixNavigator, False(), False());
//    }
}

// ----------------------------------------------------------------------------
// Function _computeCell(); InnerCol; LastCell
// ----------------------------------------------------------------------------

// Values of last call are copied into the horizontal buffer for initializing next tile below.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const & previousDiagonal,
             DPCell_<TScoreValue, TGapCosts> const & previousHorizontal,
             DPCell_<TScoreValue, TGapCosts> const & previousVertical,
             TSequenceHValue const & seqHVal,
             TSequenceVValue const & seqVVal,
             TScoringScheme const & scoringScheme,
             MetaColumnDescriptor<DPInnerColumn, FullColumn> const &,
             LastCell const & /*cellDescriptor*/,
             DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel> const &)
{
    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, FullColumn> >     TMetaColumn;

    assignValue(traceMatrixNavigator,
                _computeScore(activeCell, previousDiagonal, previousHorizontal, previousVertical, seqHVal, seqVVal,
                              scoringScheme, typename RecursionDirection_<TMetaColumn, LastCell>::Type(),
                              TDPProfile()));
    // Copy values into horizontal buffer for the tile below this tile in vertical direction.
//    auto pos = coordinate(traceMatrixNavigator, +DPMatrixDimension_::HORIZONTAL) - 1;
    (*scout.state.ptrHorBuffer)[scout.mHorizontalPos - 1].i1 = activeCell;
    if (IsTracebackEnabled_<TTraceConfig>::VALUE)
    {
        (*scout.state.ptrHorBuffer)[scout.mHorizontalPos - 1].i2 = value(traceMatrixNavigator);
    }

    if (TrackingEnabled_<TMetaColumn, LastCell>::VALUE)
    {
        _scoutBestScore(scout, activeCell, traceMatrixNavigator, False(), True());
    }
}


// ----------------------------------------------------------------------------
// Function _computeCell(); FinalCol; FirstCell
// ----------------------------------------------------------------------------

// Horizontal initialization values are copied from buffer for all first cells.
// Vertical buffer is filled with value.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const & /*previousDiagonal*/,
             DPCell_<TScoreValue, TGapCosts> const & /*previousHorizontal*/,
             DPCell_<TScoreValue, TGapCosts> const & /*previousVertical*/,
             TSequenceHValue const & /*seqHVal*/,
             TSequenceVValue const & /*seqVVal*/,
             TScoringScheme const & /*scoringScheme*/,
             MetaColumnDescriptor<DPFinalColumn, FullColumn> const &,
             FirstCell const &,   // One of FirstCell, InnerCell or LastCell.
             DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel> const &)
{
    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, FullColumn> >     TMetaColumn;

//    auto pos = coordinate(traceMatrixNavigator, +DPMatrixDimension_::HORIZONTAL) - 1;
    activeCell = front(*scout.state.ptrVerBuffer).i1 = (*scout.state.ptrHorBuffer)[scout.mHorizontalPos - 1].i1;  // Copy horizontal buffer value in active cell and in
    assignValue(traceMatrixNavigator, (*scout.state.ptrHorBuffer)[scout.mHorizontalPos - 1].i2);
    if (IsTracebackEnabled_<TTraceConfig>::VALUE)
    {
        front(*scout.state.ptrVerBuffer).i2 = value(traceMatrixNavigator);   // Store trace value in vertical buffer.
    }

    if (TrackingEnabled_<TMetaColumn, FirstCell>::VALUE)
    {
        _scoutBestScore(scout, activeCell, traceMatrixNavigator, True(), False());
    }
}

// ----------------------------------------------------------------------------
// Function _computeCell(); FinalCol, InnerCell;
// ----------------------------------------------------------------------------

// Stores computed values in vertical buffer for initializing next tile right of the current.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const & previousDiagonal,
             DPCell_<TScoreValue, TGapCosts> const & previousHorizontal,
             DPCell_<TScoreValue, TGapCosts> const & previousVertical,
             TSequenceHValue const & seqHVal,
             TSequenceVValue const & seqVVal,
             TScoringScheme const & scoringScheme,
             MetaColumnDescriptor<DPFinalColumn, FullColumn> const &,
             InnerCell const &,
             DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel> const &)
{
    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, FullColumn> >     TMetaColumn;

    assignValue(traceMatrixNavigator,
                _computeScore(activeCell, previousDiagonal, previousHorizontal, previousVertical, seqHVal, seqVVal,
                              scoringScheme, typename RecursionDirection_<TMetaColumn, InnerCell>::Type(),
                              TDPProfile()));
    // Store values in vertical buffer.
//    auto pos = coordinate(traceMatrixNavigator, +DPMatrixDimension_::VERTICAL);
    (*scout.state.ptrVerBuffer)[scout.mVerticalPos].i1 = activeCell;
    if (IsTracebackEnabled_<TTraceConfig>::VALUE)
    {
        (*scout.state.ptrVerBuffer)[scout.mVerticalPos].i2 = value(traceMatrixNavigator);
    }

    if (TrackingEnabled_<TMetaColumn, InnerCell>::VALUE)
    {
        _scoutBestScore(scout, activeCell, traceMatrixNavigator, True(), False());
    }
}

// ----------------------------------------------------------------------------
// Function _computeCell(); FinalCol, LastCell;
// ----------------------------------------------------------------------------

// Stores computed values in vertical buffer for initializing next tile right of the current.
// Stores computed values in horizontal buffer for initializing next tile below.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const & previousDiagonal,
             DPCell_<TScoreValue, TGapCosts> const & previousHorizontal,
             DPCell_<TScoreValue, TGapCosts> const & previousVertical,
             TSequenceHValue const & seqHVal,
             TSequenceVValue const & seqVVal,
             TScoringScheme const & scoringScheme,
             MetaColumnDescriptor<DPFinalColumn, FullColumn> const &,
             LastCell const &,
             DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel> const &)
{
    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, FullColumn> >     TMetaColumn;

    assignValue(traceMatrixNavigator,
                _computeScore(activeCell, previousDiagonal, previousHorizontal, previousVertical, seqHVal, seqVVal,
                              scoringScheme, typename RecursionDirection_<TMetaColumn, LastCell>::Type(),
                              TDPProfile()));
    // Store values in vertical and horizontal buffer.
//    auto posH = coordinate(traceMatrixNavigator, +DPMatrixDimension_::HORIZONTAL) - 1;
//    auto posV = coordinate(traceMatrixNavigator, +DPMatrixDimension_::VERTICAL);
    (*scout.state.ptrHorBuffer)[scout.mHorizontalPos - 1].i1 = (*scout.state.ptrVerBuffer)[scout.mVerticalPos].i1 = activeCell;
    if (IsTracebackEnabled_<TTraceConfig>::VALUE)
    {
        (*scout.state.ptrHorBuffer)[scout.mHorizontalPos - 1].i2 =
            (*scout.state.ptrVerBuffer)[scout.mVerticalPos].i2 = value(traceMatrixNavigator);
    }

    if (TrackingEnabled_<TMetaColumn, LastCell>::VALUE)
    {
        _scoutBestScore(scout, activeCell, traceMatrixNavigator, True(), True());
    }
}

// ----------------------------------------------------------------------------
// Function _computeAlignment(); Parallel;
// ----------------------------------------------------------------------------

// Overloaded _computeAlignment interface for the parallel spec.
template <typename TDPScoreValue, typename TTraceValue, typename TScoreMatHost, typename TTraceMatHost,
          typename TTraceTarget,
          typename TBuffer, typename TThreadContext, typename TSpec,
          typename TSequenceH,
          typename TSequenceV,
          typename TScoreScheme,
          typename TBandSwitch,
          typename TAlignmentAlgorithm, typename TGapScheme, typename TTraceFlag>
inline void
_computeAlignment(DPContext<TDPScoreValue, TTraceValue, TScoreMatHost, TTraceMatHost> & dpContext,
                  TTraceTarget & /*traceSegments*/,
                  DPScoutState_<DPTiled<TBuffer, TThreadContext, TSpec> > & scoutState,
                  TSequenceH const & seqH,
                  TSequenceV const & seqV,
                  TScoreScheme const & scoreScheme,
                  DPBandConfig<TBandSwitch> const & band,
                  DPProfile_<TAlignmentAlgorithm, TGapScheme, TTraceFlag, Parallel> const & dpProfile,
                  bool const lastCol,
                  bool const lastRow)
{
    typedef typename DefaultScoreMatrixSpec_<TAlignmentAlgorithm>::Type TScoreMatrixSpec;

    typedef DPMatrix_<TDPScoreValue, TScoreMatrixSpec, TScoreMatHost>   TDPScoreMatrix;
    typedef DPMatrix_<TTraceValue, FullDPMatrix, TTraceMatHost>         TDPTraceMatrix;

    typedef DPMatrixNavigator_<TDPScoreMatrix, DPScoreMatrix, NavigateColumnWise> TDPScoreMatrixNavigator;
    typedef DPMatrixNavigator_<TDPTraceMatrix, DPTraceMatrix<TTraceFlag>, NavigateColumnWise> TDPTraceMatrixNavigator;

    typedef typename ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm, DPScoutState_<DPTiled<TBuffer, TThreadContext, TSpec> > >::Type TScoutSpec;
    typedef DPScout_<TDPScoreValue, TScoutSpec> TDPScout;

    // Check if current dp settings are valid. If not return infinity value for dp score value.
    if (!_isValidDPSettings(seqH, seqV, band, dpProfile))
        return; //createVector<TScoreValue>(MinValue<typename Value<TScoreValue>::Type>::VALUE);

    TDPScoreMatrix dpScoreMatrix;
    TDPTraceMatrix dpTraceMatrix;

    // TODO(rmaerker): Check whether the matrix allocation can be reduced if upperDiagonal < 0?
    setLength(dpScoreMatrix, +DPMatrixDimension_::HORIZONTAL, length(seqH) + 1 - std::max(0, lowerDiagonal(band)));
    setLength(dpTraceMatrix, +DPMatrixDimension_::HORIZONTAL, length(seqH) + 1 - std::max(0, lowerDiagonal(band)));

    if (IsSameType<TBandSwitch, BandOff>::VALUE)
    {
        setLength(dpScoreMatrix, +DPMatrixDimension_::VERTICAL, length(seqV) + 1);
        setLength(dpTraceMatrix, +DPMatrixDimension_::VERTICAL, length(seqV) + 1);
    }
    else
    {
        SEQAN_ASSERT_FAIL("Banded version not supported!");
    }

    // We set the host to the score matrix and the dp matrix.
    setHost(dpScoreMatrix, getDpScoreMatrix(dpContext));
    setHost(dpTraceMatrix, getDpTraceMatrix(dpContext));

    resize(dpScoreMatrix);
    // We do not need to allocate the memory for the trace matrix if the traceback is disabled.
    if (IsTracebackEnabled_<TTraceFlag>::VALUE)
    {
        resize(dpTraceMatrix);
    }

    TDPScoreMatrixNavigator dpScoreMatrixNavigator;
    TDPTraceMatrixNavigator dpTraceMatrixNavigator;

    _init(dpScoreMatrixNavigator, dpScoreMatrix, band);
    _init(dpTraceMatrixNavigator, dpTraceMatrix, band);

    TDPScout dpScout(scoutState, IsLocalAlignment_<TAlignmentAlgorithm>::VALUE);  // Now initalize the scout with the state from outside, which knows the buffer of the current tile.

    // Execute the alignment.
    if (!_isBandEnabled(band))
        _computeUnbandedAlignment(dpScout, dpScoreMatrixNavigator, dpTraceMatrixNavigator, seqH, seqV, scoreScheme, dpProfile);
    else
        SEQAN_ASSERT_FAIL("Banded version not supported!");

//    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
//        return;

    // Some check to get same behavior as old module.
    if (IsTracebackEnabled_<TTraceFlag>::VALUE &&
        ((IsFreeEndGap_<TAlignmentAlgorithm, DPLastColumn>::VALUE && lastCol) ||
        (IsFreeEndGap_<TAlignmentAlgorithm, DPLastRow>::VALUE && lastRow) ||
        (lastCol && lastRow)))
    {
        // Check if max was found at the bottom right corner of the matrix.
        // This is also true if in last row, and last column
        //        if ((maxHostPosition(dpScout) + 1) == (end(dpTraceMatrix) - begin(dpTraceMatrix)))

        //            maxHostPosition(dpScout); // We only have the trace value not the score value.
        _correctTraceValue(dpTraceMatrixNavigator, dpScout);
    }
    combineMaxScore(scoutState.mThreadContext, dpScout, dpTraceMatrix, IsTracebackEnabled_<TTraceFlag>::VALUE);
//    return dpScout; // Needed for local/semi-global alignment.
}

// ----------------------------------------------------------------------------
// Function implParallelTrace()
// ----------------------------------------------------------------------------

// The traceback wrapper to go backwards from block to block and call the internal traceback function.
template <typename TTarget,
          typename TTraceBlockIdH,
          typename TTraceBlockIdV,
          typename TTraceBlockPos,
          typename TTraceProxy,
          typename TSequenceH,
          typename TSequenceV,
          typename TAlgorithm, typename TGapCosts, typename TTracebackSpec>
void implParallelTrace(TTarget & target,
                       TTraceBlockIdH const idH,
                       TTraceBlockIdV const idV,
                       TTraceBlockPos const blockPos,
                       TTraceProxy const & traceProxy,
                       TSequenceH const & seqH,
                       TSequenceV const & seqV,
                       DPProfile_<TAlgorithm, TGapCosts, TTracebackSpec, Parallel> const & dpProfile)
{
    typedef DPProfile_<TAlgorithm, TGapCosts, TTracebackSpec, Parallel> DPProfile;

    using TNavigator = BlockTraceNavigator<TTraceProxy const, TSequenceH const, TSequenceV const, TTracebackSpec>;
    using TPos = typename impl::LocalPosition<TNavigator>::Type;

    // Initialized the blockNavi.
    TNavigator navi(traceProxy, seqH, seqV, idH, idV, blockPos);
    // extract firstValue?
    auto traceValue = scalarValue(navi);
    auto lastTraceValue = _retrieveInitialTraceDirection(traceValue, dpProfile);

    _computeTraceback(target, traceValue, lastTraceValue, navi,
                      impl::toGlobalPosition(TPos(length(seqH) - 1, length(back(seqH))), navi, +DPMatrixDimension_::HORIZONTAL),
                      impl::toGlobalPosition(TPos(length(seqV) - 1, length(back(seqV))), navi, +DPMatrixDimension_::VERTICAL),
                      DPBandConfig<BandOff>(), DPProfile(), True(), True());
}
// ----------------------------------------------------------------------------
// Function implParallelAlign()
// ----------------------------------------------------------------------------

// The block wise wrapper interface to arrange the dp matrix into blocks and process them via the minor diagonal.
template <typename TParSpec, typename TVecSpec,
          typename TTarget,
          typename TSeqH,
          typename TSeqV,
          typename TScoreValue, typename TScoreSpec>
inline TScoreValue
implParallelAlign(ExecutionPolicy<TParSpec, TVecSpec> const & execPolicy,
                  TTarget & target,
                  TSeqH const & seqH,
                  TSeqV const & seqV,
                  Score<TScoreValue, TScoreSpec> const & score)
{
    typedef typename TraceBitMap_<>::Type                   TTraceValue;
    typedef DPCell_<TScoreValue, AffineGaps>                TDPCell;        // Type of the DPCell; at moment this is fixed to AffineGaps.
    typedef Pair<TDPCell, TTraceValue>                      TBufferValue;   // Value type for the tile buffer: Used to store last row or column of a tile to intialize the neighboring tile.
    typedef String<TBufferValue>                            TBuffer;
    typedef DPTileBuffer<TBuffer>                           TTileBuffer;    // Represents the global buffer for all blocks.

    // ----------------------------------------------------------------------------
    // Standard configuration.

    // TODO(rrahn): Refine design. This design is based on the internal DP module and must be revised when changing the public API interface.
    typedef typename SubstituteAlignConfig_<AlignConfig<> >::Type TFreeEndGaps;
    typedef DPProfile_<GlobalAlignment_<TFreeEndGaps>, AffineGaps, TracebackOff/*TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> >*/, Parallel> TDPProfile;  // A type trait to configure the DP algorithm.

    // ----------------------------------------------------------------------------
    // Initialize TileBuffer.

    // The buffer is initialized, depending on the DPProfile. So we use the internal _doComputeScore function for this.
    TTileBuffer tileBuffer;
    resize(tileBuffer.horizontalBuffer, length(seqH), Exact());
    resize(tileBuffer.verticalBuffer, length(seqV), Exact());

    TBufferValue tmp;
    tmp.i2 = _doComputeScore(tmp.i1, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), score, RecursionDirectionZero(), TDPProfile());
    for (auto itH = begin(tileBuffer.horizontalBuffer, Standard()); itH != end(tileBuffer.horizontalBuffer, Standard()); ++itH)
    {
        resize(*itH, length(front(seqH)), Exact());
        for (auto it = begin(*itH, Standard()); it != end(*itH, Standard()); ++it)
        {
            it->i2 = _doComputeScore(it->i1, TDPCell(), tmp.i1, TDPCell(), Nothing(), Nothing(), score, RecursionDirectionHorizontal(), TDPProfile());
            tmp.i1 = it->i1;
        }
    }
    tmp.i2 = _doComputeScore(tmp.i1, TDPCell(), TDPCell(), TDPCell(), Nothing(), Nothing(), score, RecursionDirectionZero(), TDPProfile());
    for (auto itV = begin(tileBuffer.verticalBuffer, Standard()); itV != end(tileBuffer.verticalBuffer, Standard()); ++itV)
    {
        resize(*itV, length(front(seqV)) + 1, Exact());
        auto it = begin(*itV, Standard());
        it->i2 = tmp.i2;
        it->i1 = tmp.i1;
        ++it;
        for (; it != end(*itV, Standard()); ++it)
        {
            it->i2 = _doComputeScore(it->i1, TDPCell(), TDPCell(), tmp.i1, Nothing(), Nothing(), score, RecursionDirectionVertical(), TDPProfile());
            tmp.i1 = it->i1;
            tmp.i2 = it->i2;  // TODO(rrahn): Move out of loop.
        }
    }


    // DEBUG: DebugMatrix
//    impl::debug::DebugBuffer<TBuffer> debugMatrix;
//
//    resize(debugMatrix.matrix, length(seqH), Exact());
//    for (auto& column : debugMatrix.matrix)
//        resize(column, length(seqV), Exact());
    // DEBUG: End

    using TSimdVec = typename SimdVector<TScoreValue>::Type;
    using TLocalTraceStore = impl::dp::parallel::LocalTraceStore<TSimdVec>;
    impl::dp::parallel::TraceProxy<TLocalTraceStore> traceProxy(length(seqH), length(seqV));

    // Adaption of the dag code.
//    typedef DPContext<TDPCell, typename TraceBitMap_<>::Type, String<TDPCell>, TTraceTile> TDPContext;

    using TTaskContext = DPTaskContext<TSeqH const *, TSeqV const *,
                                       Score<TScoreValue, TScoreSpec> const *,
                                       DPBandConfig<BandOff> *,
                                       TTileBuffer *,
                                       decltype(traceProxy)*,
//                                       impl::debug::DebugBuffer<TBuffer>*,  // remove after debugging.
                                       TDPProfile,
                                       TBuffer>;

    DPBandConfig<BandOff> dpBand;
    TTaskContext taskContext(&seqH, &seqV, &score, &dpBand, &tileBuffer, &traceProxy); // , &debugMatrix);

    typename impl::dp::parallel::ThreadLocalStorage<impl::dp::parallel::DPLocalStorage<TScoreValue, TSimdVec>, TParSpec>::Type tls;

//    auto begin = sysTime();
    auto taskGraph = createGraph(taskContext, tls, execPolicy);
//    std::cout << "\nCreation: " << std::setw(15) << sysTime() - begin << "s\n";
    invoke(taskGraph, execPolicy);

    std::tuple<TScoreValue, size_t, size_t, size_t> dpMax{MinValue<TScoreValue>::VALUE, 0, 0, 0};
    combine(tls, [&dpMax](auto& localDPStore)
    {
        if (std::get<0>(dpMax) < localDPStore.mMaxScore)
        {
            std::get<0>(dpMax) = localDPStore.mMaxScore;
            std::get<1>(dpMax) = localDPStore.mMaxBlockHId;
            std::get<2>(dpMax) = localDPStore.mMaxBlockVId;
            std::get<3>(dpMax) = localDPStore.mMaxBlockPos;
        }
    });

    // TODO(rrahn): pair<score, tuple<blockH, blockV, blockPos>> combine(tls);

//    std::ofstream fStream("/Users/rmaerker/Documents/par_native.csv");
//    debugMatrix.write(fStream);
//    fStream.close();

//    std::cout << "Block: (" << std::get<1>(dpMax) << ", " << std::get<2>(dpMax) << ")\n";
//    std::cout << "BlockPos: " << std::get<3>(dpMax) << '\n';

    // TODO(rrahn): Fix to also support local and semi-global alignment!
    if (IsTracebackEnabled_<TDPProfile>::VALUE)
    {
        implParallelTrace(target, std::get<1>(dpMax), std::get<2>(dpMax), std::get<3>(dpMax), traceProxy,
                          seqH, seqV, TDPProfile());
    }

    return std::get<0>(dpMax);
}

// ----------------------------------------------------------------------------
// Function parallelAlign()
// ----------------------------------------------------------------------------

// Main function to be called by the user.
// The interface is minimalistic here, since we plan to revise the user interface
// in the future.
template <typename TExecPolicy,
          typename TTargetH,
          typename TTargetV,
          typename TScoreValue, typename TScoreSpec>
inline TScoreValue
parallelAlign(TExecPolicy const & policy,
              TTargetH & targetH,
              TTargetV & targetV,
              Score<TScoreValue, TScoreSpec> const & score,
              unsigned const blockSize = 100)
{
    typedef typename Source<TTargetH>::Type                                   TSeqH;
    typedef typename Source<TTargetV>::Type                                   TSeqV;
    typedef typename Iterator<typename Infix<TSeqH>::Type, Standard>::Type    TIterH;
    typedef typename Iterator<typename Infix<TSeqV>::Type, Standard>::Type    TIterV;


    typedef typename Position<TTargetH>::Type                                 TPosition;
    typedef typename Size<TTargetH>::Type                                     TSize;
    typedef TraceSegment_<TPosition, TSize>                                   TTraceSegment;

    // Variable used to set the length of a bock.
    // Blocks are quadratic. Only the blocks in the last block column or block row might differ in there size,
    // depending on the length of the original sequences.

    // Blocks are implemented as array of Ranges (storing begin and end iterator to corresponding block).
    String<Range<TIterH> >  seqHBlocks;
    String<Range<TIterV> >  seqVBlocks;

    // ----------------------------------------------------------------------------
    // Initialization of block structure.

    // Partition sequences into blocks
    //    std::cout << "Blocks: " << ((length(seqH) << 1) - 1) / BLOCK_SIZE << " == " << (length(seqH) + BLOCK_SIZE - 1) / BLOCK_SIZE << std::endl;
    resize(seqHBlocks, (length(source(targetH)) + blockSize - 1) / blockSize, Exact());
    resize(seqVBlocks, (length(source(targetV)) + blockSize - 1) / blockSize, Exact());

    for (unsigned id = 0; id < length(seqHBlocks); ++id)
        seqHBlocks[id] = toRange(infix(source(targetH), id * blockSize, _min(length(source(targetH)),(id + 1) * blockSize)));

    for (unsigned id = 0; id < length(seqVBlocks); ++id)
        seqVBlocks[id] = toRange(infix(source(targetV), id * blockSize, _min(length(source(targetV)),(id + 1) * blockSize)));

    String<TTraceSegment> trace;  // Stores the trace segments which are converted at the end into the final alignment.

    // The main blockwise wrapper for the dp algorithm.
    TScoreValue res = implParallelAlign(policy, trace, seqHBlocks, seqVBlocks, score);
    _adaptTraceSegmentsTo(targetH, targetV, trace);  // Convert to corresponding alignment representation.
    return res;
}
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_ALIGN_PARALLEL_IMPL_H_
