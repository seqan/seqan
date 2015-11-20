// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_ALIGN_DP_PARALLEL_H_
#define INCLUDE_SEQAN_ALIGN_DP_PARALLEL_H_

// #define __ORG_OMP_
#define __TBB_TASKS_

#ifdef __TBB_TASKS_
#include "tbb/task.h"
#include "tbb/enumerable_thread_specific.h"
#endif

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag DPTiled<TBuffer>
// ----------------------------------------------------------------------------

// Tag used to subclass DPScoutState and DPScout.
// T represents the buffer type.
template <typename T>
struct DPTiled;

// ----------------------------------------------------------------------------
// Class DPTileBuffer
// ----------------------------------------------------------------------------

// The structure owning the horizontal/vertical buffer.
template <typename TDPCellBuff, typename TBuffer = String<TDPCellBuff> >
struct DPTileBuffer
{
    TBuffer horizontalBuffer;
    TBuffer verticalBuffer;
};

// ----------------------------------------------------------------------------
// Class DPScoutState_; DPTiled
// ----------------------------------------------------------------------------

// The overloaded DPScoutState which simply stores the pointers to the corresponding buffer.
template <typename TBuffer>
class DPScoutState_<DPTiled<TBuffer> >
{
public:

    TBuffer* ptrHorBuffer;
    TBuffer* ptrVerBuffer;

    DPScoutState_() : ptrHorBuffer(nullptr), ptrVerBuffer(nullptr)
    {}

    DPScoutState_(TBuffer & horBuffer, TBuffer & verBuffer) :
        ptrHorBuffer(&horBuffer),
        ptrVerBuffer(&verBuffer)
    {}
};

// ----------------------------------------------------------------------------
// Class DPScout_; DPTiled
// ----------------------------------------------------------------------------

// Overloaded DPScout to store the corresponding buffer for the current dp tile.
template <typename TDPCell, typename TBuffer>
class DPScout_<TDPCell, DPTiled<TBuffer> > :
    public DPScout_<TDPCell, Default>
{
public:
    typedef DPScout_<TDPCell, Default>  TParent;

    DPScoutState_<DPTiled<TBuffer> >    state;

    DPScout_() :
        TParent(),
        state()
    {}

    DPScout_(DPScoutState_<DPTiled<TBuffer> > state) :
        TParent(),
        state(state)
    {}

    DPScout_(DPScout_ const & other) :
        TParent(other),
        state(other.state)
    {}

    DPScout_ & operator=(DPScout_ const & other)
    {
        if (this != &other)
        {
            TParent::operator=(other);
            state = other.state;
        }
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================


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
    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInitialColumn, FullColumn> >   TMetaColumn;

    auto pos = coordinate(traceMatrixNavigator, +DPMatrixDimension_::VERTICAL);
    activeCell = (*scout.state.ptrVerBuffer)[pos].i1;
    assignValue(traceMatrixNavigator, (*scout.state.ptrVerBuffer)[pos].i2);

    if (TrackingEnabled_<TMetaColumn, TCellDescriptor>::VALUE)
    {
        typedef typename And<IsSameType<TCellDescriptor, LastCell>,
                             Or<IsSameType<typename MetaColumnDescriptor<DPInitialColumn, FullColumn>::TLocation, PartialColumnBottom>,
                                IsSameType<typename MetaColumnDescriptor<DPInitialColumn, FullColumn>::TLocation, FullColumn>
                             >
                         >::Type TIsLastRow;
        _scoutBestScore(scout, activeCell, traceMatrixNavigator, False(), TIsLastRow());
    }
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
    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, FullColumn> >     TMetaColumn;

    auto pos = coordinate(traceMatrixNavigator, +DPMatrixDimension_::HORIZONTAL) - 1;
    activeCell = (*scout.state.ptrHorBuffer)[pos].i1;
    assignValue(traceMatrixNavigator, (*scout.state.ptrHorBuffer)[pos].i2);

    if (TrackingEnabled_<TMetaColumn, FirstCell>::VALUE)
    {
        _scoutBestScore(scout, activeCell, traceMatrixNavigator, False(), False());
    }
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
    auto pos = coordinate(traceMatrixNavigator, +DPMatrixDimension_::HORIZONTAL) - 1;
    (*scout.state.ptrHorBuffer)[pos].i1 = activeCell;
    (*scout.state.ptrHorBuffer)[pos].i2 = TraceBitMap_::NONE;

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

    activeCell = front(*scout.state.ptrVerBuffer).i1 = back(*scout.state.ptrHorBuffer).i1;  // Copy horizontal buffer value in active cell and in
    assignValue(traceMatrixNavigator, back(*scout.state.ptrHorBuffer).i2);
    front(*scout.state.ptrVerBuffer).i2 = TraceBitMap_::NONE;   // Store trace value in vertical buffer.

    if (TrackingEnabled_<TMetaColumn, FirstCell>::VALUE)
    {
        _scoutBestScore(scout, activeCell, traceMatrixNavigator, True(), False());
    }
}

// ----------------------------------------------------------------------------
// Function _computeCell(); FinalCol;
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
    auto pos = coordinate(traceMatrixNavigator, +DPMatrixDimension_::VERTICAL);
    (*scout.state.ptrVerBuffer)[pos].i1 = activeCell;
    (*scout.state.ptrVerBuffer)[pos].i2 = TraceBitMap_::NONE;

    if (TrackingEnabled_<TMetaColumn, InnerCell>::VALUE)
    {
        _scoutBestScore(scout, activeCell, traceMatrixNavigator, True(), False());
    }
}

// ----------------------------------------------------------------------------
// Function _computeCell(); FinalCol;
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
    back(*scout.state.ptrHorBuffer).i1 = back(*scout.state.ptrVerBuffer).i1 = activeCell;
    back(*scout.state.ptrHorBuffer).i2 = back(*scout.state.ptrVerBuffer).i2 = TraceBitMap_::NONE;

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
          typename TBuffer,
          typename TSequenceH,
          typename TSequenceV,
          typename TScoreScheme,
          typename TBandSwitch,
          typename TAlignmentAlgorithm, typename TGapScheme, typename TTraceFlag>
inline typename Value<TScoreScheme>::Type
_computeAlignment(DPContext<TDPScoreValue, TTraceValue, TScoreMatHost, TTraceMatHost> & dpContext,
                  TTraceTarget & /*traceSegments*/,
                  DPScoutState_<DPTiled<TBuffer> > & scoutState,
                  TSequenceH const & seqH,
                  TSequenceV const & seqV,
                  TScoreScheme const & scoreScheme,
                  DPBandConfig<TBandSwitch> const & band,
                  DPProfile_<TAlignmentAlgorithm, TGapScheme, TTraceFlag, Parallel> const & dpProfile)
{
    typedef typename DefaultScoreMatrixSpec_<TAlignmentAlgorithm>::Type TScoreMatrixSpec;

    typedef DPMatrix_<TDPScoreValue, TScoreMatrixSpec, TScoreMatHost>   TDPScoreMatrix;
    typedef DPMatrix_<TTraceValue, FullDPMatrix, TTraceMatHost>         TDPTraceMatrix;

    typedef DPMatrixNavigator_<TDPScoreMatrix, DPScoreMatrix, NavigateColumnWise> TDPScoreMatrixNavigator;
    typedef DPMatrixNavigator_<TDPTraceMatrix, DPTraceMatrix<TTraceFlag>, NavigateColumnWise> TDPTraceMatrixNavigator;

    typedef DPScout_<TDPScoreValue, DPTiled<TBuffer> > TDPScout;

    typedef typename Value<TScoreScheme>::Type TScoreValue;

    // Check if current dp settings are valid. If not return infinity value for dp score value.
    if (!_isValidDPSettings(seqH, seqV, band, dpProfile))
        return MinValue<TScoreValue>::VALUE;

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
        resize(dpTraceMatrix);  // Probably never needed here.

    TDPScoreMatrixNavigator dpScoreMatrixNavigator;
    TDPTraceMatrixNavigator dpTraceMatrixNavigator;

    _init(dpScoreMatrixNavigator, dpScoreMatrix, band);
    _init(dpTraceMatrixNavigator, dpTraceMatrix, band);

    TDPScout dpScout(scoutState);  // Now initalize the scout with the state from outside, which knows the buffer of the current tile.

    // Execute the alignment.
    if (!_isBandEnabled(band))
        _computeUnbandedAlignment(dpScout, dpScoreMatrixNavigator, dpTraceMatrixNavigator, seqH, seqV, scoreScheme, dpProfile);
    else
        SEQAN_ASSERT_FAIL("Banded version not supported!");

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return maxScore(dpScout);

    // Some check to get same behavior as old module.
    if (IsSingleTrace_<TTraceFlag>::VALUE)
    {
        // Check if max was found at the bottom right corner of the matrix.
        // This is also true if in last row, and last column
        //        if ((maxHostPosition(dpScout) + 1) == (end(dpTraceMatrix) - begin(dpTraceMatrix)))
        
        //            maxHostPosition(dpScout); // We only have the trace value not the score value.
        _correctTraceValue(dpTraceMatrixNavigator, dpScout);
    }
    return maxScore(dpScout); // Needed for local/semi-global alignment.
}

// ----------------------------------------------------------------------------
// Function implParallelTrace()
// ----------------------------------------------------------------------------

// The traceback wrapper to go backwards from block to block and call the internal traceback function.
template <typename TTarget,
          typename TPosition,
          typename TIdxH,
          typename TIdxV,
          typename TTracebackBlocks,
          typename TSequenceH,
          typename TSequenceV,
          typename TAlgorithm, typename TGapCosts, typename TTracebackSpec>
void implParallelTrace(TTarget & target,
                       TPosition blockPos,
                       TIdxH idxH,
                       TIdxV idxV,
                       TTracebackBlocks const & tracebackBlocks,
                       TSequenceH const & seqH,
                       TSequenceV const & seqV,
                       DPProfile_<TAlgorithm, TGapCosts, TTracebackSpec, Parallel> const & /*dpProfile*/,
                       typename TraceBitMap_::TTraceValue lastTraceValue = TraceBitMap_::NONE)
{
    typedef DPProfile_<TAlgorithm, TGapCosts, TTracebackSpec, Parallel> DPProfile;

    typedef typename Value<TTracebackBlocks>::Type                                                  TTraceMatHost;
    typedef typename Value<TTraceMatHost>::Type                                                     TTraceValue;
    typedef DPMatrix_<TTraceValue, FullDPMatrix, TTraceMatHost>                                     TDPTraceMatrix;
    typedef DPMatrixNavigator_<TDPTraceMatrix, DPTraceMatrix<TTracebackSpec>, NavigateColumnWise>   TDPTraceMatrixNavigator;


    // Emulate dp trace matrix and fill it with the corresponding trace block.
    // Thus, we can use the existing traceback engine to compute the traceback.
    TDPTraceMatrix dpTraceMatrix;
    setLength(dpTraceMatrix, +DPMatrixDimension_::HORIZONTAL, length(seqH[idxH]) + 1);
    setLength(dpTraceMatrix, +DPMatrixDimension_::VERTICAL, length(seqV[idxV]) + 1);
    setHost(dpTraceMatrix, tracebackBlocks[idxH * length(seqV) + idxV]);
    resize(dpTraceMatrix);  // Needed to actually update the internal variables.

    TDPTraceMatrixNavigator navi;
    _init(navi, dpTraceMatrix, DPBandConfig<BandOff>());  // Always assume Band off.
    _setToPosition(navi, blockPos);  // Set the navigator to the correct position within the block.


    // Delegation to the internal traceback engine.
    // TODO(rrahn): Fix this to work with local and semi-global as well.
    if ((idxH == length(seqH) - 1 && idxV == length(seqV) - 1) && IsGlobalAlignment_<TAlgorithm>::VALUE)
    {
        lastTraceValue = _retrieveInitialTraceDirection(value(navi), DPProfile());
        _computeTraceback(target, value(navi), lastTraceValue, navi, seqH[idxH], seqV[idxV], DPBandConfig<BandOff>(), DPProfile(), False(), True());
    }
    else if ((idxH == 0 && idxV == 0) && IsGlobalAlignment_<TAlgorithm>::VALUE)
        _computeTraceback(target, value(navi), lastTraceValue, navi, seqH[idxH], seqV[idxV], DPBandConfig<BandOff>(), DPProfile(), True(), False());
    else
        _computeTraceback(target, value(navi), lastTraceValue, navi, seqH[idxH], seqV[idxV], DPBandConfig<BandOff>(), DPProfile(), False(), False());

    // Check if we reached end of alignment.
    if (idxH == 0 && idxV == 0)  // TODO(rmaerker): Need to be updated for the local alignment and the semi-global alignment.
        return;

    // Now we need to determine the position where we start the tracing in the next block.
    // case 1: alignment ends in cell (0,0)
    if (position(navi) == 0)
    {
        --idxH; --idxV;
        blockPos = length(tracebackBlocks[idxH * length(seqV) + idxV]) - 1;
    }  // Case 2: alignment ends in cell (0,j).
    else if (coordinate(navi, +DPMatrixDimension_::HORIZONTAL) == 0)  // Continues in block left!
    {
        SEQAN_ASSERT_GT(coordinate(navi, +DPMatrixDimension_::VERTICAL), 0u);
        --idxH;
        // Block pos is last column times length of the column + coordinate in vertical dimension.
        blockPos = length(seqH[idxH]) * (length(seqV[idxV]) + 1) + coordinate(navi, +DPMatrixDimension_::VERTICAL);
    }
    else  // Case 3: alignment ends in cell (i,0).
    {
        SEQAN_ASSERT_EQ(coordinate(navi, +DPMatrixDimension_::VERTICAL), 0u);
        SEQAN_ASSERT_GT(coordinate(navi, +DPMatrixDimension_::HORIZONTAL), 0u);
        --idxV;
        blockPos = coordinate(navi, +DPMatrixDimension_::HORIZONTAL) * (length(seqV[idxV]) + 1) + length(seqV[idxV]);
    }

    // Recursive call to compute the traceback.
    implParallelTrace(target, blockPos, idxH, idxV, tracebackBlocks, seqH, seqV, DPProfile(), lastTraceValue);
}

#ifdef __ORG_OMP_

// ----------------------------------------------------------------------------
// Function implParallelAlign()
// ----------------------------------------------------------------------------

// The block wise wrapper interface to arrange the dp matrix into blocks and process them via the minor diagonal.
template <typename TTarget,
          typename TTraceHost,
          typename TSeqH,
          typename TSeqV,
          typename TScoreValue, typename TScoreSpec>
inline TScoreValue
implParallelAlign(TTarget & target,
                  TTraceHost & traceHost,
                  TSeqH const & seqH,
                  TSeqV const & seqV,
                  Score<TScoreValue, TScoreSpec> const & score)
{
    // Iterator over the blocks.
    typedef typename Iterator<TSeqH const, Standard>::Type  TBlockH;
    typedef typename Iterator<TSeqV const, Standard>::Type  TBlockV;

    // Iterator to the traceback host (the original memory for the trace matrix).
    typedef typename Iterator<TTraceHost, Standard>::Type   TTraceIter;
    typedef Range<TTraceIter>                               TTraceTile;  // Defines the tiles for the trace matrix.

    typedef typename Value<TTraceHost>::Type                TTraceValue;
    typedef DPCell_<TScoreValue, AffineGaps>                TDPCell;        // Type of the DPCell; at moment this is fixed to AffineGaps.
    typedef Pair<TDPCell, TTraceValue>                      TBufferValue;   // Value type for the tile buffer: Used to store last row or column of a tile to intialize the neighboring tile.
    typedef String<TBufferValue>                            TBuffer;
    typedef DPTileBuffer<TBuffer>                           TTileBuffer;    // Represents the global buffer for all blocks.

    // ----------------------------------------------------------------------------
    // Standard configuration.

    // TODO(rrahn): Refine design. This design is based on the internal DP module and must be revised when changing the public API interface.
    typedef typename SubstituteAlignConfig_<AlignConfig<> >::Type TFreeEndGaps;
    typedef DPProfile_<GlobalAlignment_<TFreeEndGaps>, AffineGaps, TracebackOn<>, Parallel> TDPProfile;  // A type trait to configure the DP algorithm.

    typedef DPScoutState_<DPTiled<TBuffer> > TDPScoutState;  // Through the dpScoutState we access the buffer per tile.


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
            tmp.i2 = it->i2;
        }
    }

    // ----------------------------------------------------------------------------
    // Initialize the block iterators.

    TBlockH blockHBeg = begin(seqH, Standard());
    TBlockV blockVBeg = begin(seqV, Standard());
    TBlockH blockHEnd = blockHBeg + 1;

    TBlockH blockHIt;
    TBlockV blockVIt;

    // ----------------------------------------------------------------------------
    // Initialize the trace blocks.

    TTraceIter itTraceBeg = begin(traceHost, Standard());
    TTraceIter itTraceEnd = itTraceBeg;

    // We use the tracebackBlock array as an index to access the correct part of the traceback host
    // for the current DP tile. Thus, we can arrange the traceback matrices in the order of the minor
    // diagonal, while accessing the correct block using the horziontal and vertical block index.
    // Maybe use an unordered map instead? Key would be the pair of h_block and v_block index.

    String<TTraceTile> tracebackBlock;
    resize(tracebackBlock, length(seqH) * length(seqV), Exact());

    // The number of minor diagonals.
    unsigned numDiag = length(seqH) + length(seqV);

    // Iterate over the blocks via the antidiagonals and store the current traceback block as range
    // in the corresponding field of the tracebackBlock array. This array will fit into the cache and
    // thus does not need to be arranged along the antidiagonals. This makes the access to the
    // traceback tile much easier.

    // TODO(rrahn): Check if the preprocessing is faster/slower than accessing the traceback blocks in the non-antidiagonal way.
    for (unsigned diag = 1; diag < numDiag; ++diag)
    {
        // Set iterator.
        blockVIt  = blockVBeg + _min(diag - 1, length(seqV) - 1);           // Vertical block starts at current minor diagonal or at the end of seqV.
        blockHIt  = blockHBeg + _max(0, static_cast<int>(diag) - static_cast<int>(length(seqV)));
        blockHEnd = blockHBeg + _min(diag, length(seqH));
        // Inner Loop:
        for (; blockHIt != blockHEnd; ++blockHIt, --blockVIt)
        {
            // Set the trace block.
            itTraceEnd += (length(*blockHIt) + 1) * (length(*blockVIt) + 1);
            tracebackBlock[((blockHIt - blockHBeg) * length(seqV)) + (blockVIt - blockVBeg)] = toRange(itTraceBeg, itTraceEnd);
            itTraceBeg = itTraceEnd;
        }
    }
    // Sanity check: The itTraceEnd must now be at the end of the traceHost.
    SEQAN_ASSERT_EQ(itTraceEnd, end(traceHost, Standard()));

    // ----------------------------------------------------------------------------
    // Iterate over the minor diagonal and compute the blocks.

    // TODO(rrahn): Check if dpContex should be copied to each thread or if we make an array of dpContext of the size of numThreads.
    DPContext<DPCell_<TScoreValue, AffineGaps>, typename TraceBitMap_::TTraceValue,
              String<DPCell_<TScoreValue, AffineGaps> >, TTraceTile> dpContext;

    // Outer Loop: iterate over all minor diagonals.
    for (unsigned diag = 1; diag < numDiag; ++diag)
    {
        // Update the block iterators for the minor diagonal.
        blockVIt  = blockVBeg + _min(diag - 1, length(seqV) - 1);
        blockHIt  = blockHBeg + _max(0, static_cast<int>(diag) - static_cast<int>(length(seqV)));
        blockHEnd = blockHBeg + _min(diag, length(seqH));

        // Start parallel section.
        SEQAN_OMP_PRAGMA(parallel)
        {
            int numJobs = _min(blockHEnd - blockHIt, omp_get_num_threads());   // Not more jobs than blocks per minor diagonal.
            SEQAN_OMP_PRAGMA(for firstprivate(dpContext))  // Only the dpContext needs to be copied. All other variables are arranged such that there are no data dependencies between the threads.
            for (int jobId = 0; jobId < numJobs; ++jobId)
            {
//                printf("Thread %d of %d!\n", jobId, omp_get_num_threads());
                unsigned diagChunk = (blockHEnd - blockHIt) / numJobs;  // Size of the diagonal chunk.
                TBlockH jobHIt  = blockHIt + diagChunk * jobId;         // begin of the horizontal block iterator per thread.
                TBlockH jobHEnd;                                        // end of the horizontal block iterator per thread.
                if (jobId == numJobs - 1)
                    jobHEnd = blockHEnd;
                else
                    jobHEnd = blockHIt + diagChunk * (jobId + 1);
                TBlockV jobVIt = blockVIt - diagChunk * jobId;          // begin of the vertical block iterator per thread.
                String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments. Only needed for the old interface. They are not filled.
                // Inner Loop: Iterate along the minor diagonal from lower left to upper right.
                for (; jobHIt != jobHEnd; ++jobHIt, --jobVIt)
                {
                    // Get the correct traceback block depending on the current horizontal and vertical block index.
                    getDpTraceMatrix(dpContext) = tracebackBlock[((jobHIt - blockHBeg) * length(seqV)) + (jobVIt - blockVBeg)];
                    // Define scout state and set up the correct buffer for the current tile.
                    TDPScoutState scoutState(tileBuffer.horizontalBuffer[jobHIt - blockHBeg], tileBuffer.verticalBuffer[jobVIt - blockVBeg]);
                    // Compute the alignment of the block.
                    _computeAlignment(dpContext, traceSegments, scoutState, *jobHIt, *jobVIt, score, DPBandConfig<BandOff>(), TDPProfile());
                }
            }
        }
    }

    // TODO(rrahn): Fix to also support local and semi-global alignment!
    implParallelTrace(target, length(back(tracebackBlock)) - 1, length(seqH) - 1, length(seqV) - 1, tracebackBlock,
                      seqH, seqV, TDPProfile());

    return _scoreOfCell(back(back(tileBuffer.horizontalBuffer)).i1);
}

#endif

#ifdef __TBB_TASKS_

// ----------------------------------------------------------------------------
// Function implParallelAlign()
// ----------------------------------------------------------------------------

// The block wise wrapper interface to arrange the dp matrix into blocks and process them via the minor diagonal.
template <typename TTarget,
          typename TTraceHost,
          typename TSeqH,
          typename TSeqV,
          typename TScoreValue, typename TScoreSpec>
inline TScoreValue
implParallelAlign(TTarget & target,
                  TTraceHost & traceHost,
                  TSeqH const & seqH,
                  TSeqV const & seqV,
                  Score<TScoreValue, TScoreSpec> const & score)
{
    // Iterator over the blocks.
    typedef typename Iterator<TSeqH const, Standard>::Type  TBlockH;
    typedef typename Iterator<TSeqV const, Standard>::Type  TBlockV;

    // Iterator to the traceback host (the original memory for the trace matrix).
    typedef typename Iterator<TTraceHost, Standard>::Type   TTraceIter;
    typedef Range<TTraceIter>                               TTraceTile;  // Defines the tiles for the trace matrix.

    typedef typename Value<TTraceHost>::Type                TTraceValue;
    typedef DPCell_<TScoreValue, AffineGaps>                TDPCell;        // Type of the DPCell; at moment this is fixed to AffineGaps.
    typedef Pair<TDPCell, TTraceValue>                      TBufferValue;   // Value type for the tile buffer: Used to store last row or column of a tile to intialize the neighboring tile.
    typedef String<TBufferValue>                            TBuffer;
    typedef DPTileBuffer<TBuffer>                           TTileBuffer;    // Represents the global buffer for all blocks.

    // ----------------------------------------------------------------------------
    // Standard configuration.

    // TODO(rrahn): Refine design. This design is based on the internal DP module and must be revised when changing the public API interface.
    typedef typename SubstituteAlignConfig_<AlignConfig<> >::Type TFreeEndGaps;
    typedef DPProfile_<GlobalAlignment_<TFreeEndGaps>, AffineGaps, TracebackOn<>, Parallel> TDPProfile;  // A type trait to configure the DP algorithm.

    typedef DPScoutState_<DPTiled<TBuffer> > TDPScoutState;  // Through the dpScoutState we access the buffer per tile.


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
            tmp.i2 = it->i2;
        }
    }

    // ----------------------------------------------------------------------------
    // Initialize the block iterators.

    TBlockH blockHBeg = begin(seqH, Standard());
    TBlockV blockVBeg = begin(seqV, Standard());
    TBlockH blockHEnd = blockHBeg + 1;

    TBlockH blockHIt;
    TBlockV blockVIt;

    // ----------------------------------------------------------------------------
    // Initialize the trace blocks.

    TTraceIter itTraceBeg = begin(traceHost, Standard());
    TTraceIter itTraceEnd = itTraceBeg;

    // We use the tracebackBlock array as an index to access the correct part of the traceback host
    // for the current DP tile. Thus, we can arrange the traceback matrices in the order of the minor
    // diagonal, while accessing the correct block using the horziontal and vertical block index.
    // Maybe use an unordered map instead? Key would be the pair of h_block and v_block index.

    String<TTraceTile> tracebackBlock;
    resize(tracebackBlock, length(seqH) * length(seqV), Exact());

    // The number of minor diagonals.
    unsigned numDiag = length(seqH) + length(seqV);

    // Iterate over the blocks via the antidiagonals and store the current traceback block as range
    // in the corresponding field of the tracebackBlock array. This array will fit into the cache and
    // thus does not need to be arranged along the antidiagonals. This makes the access to the
    // traceback tile much easier.

    // TODO(rrahn): Check if the preprocessing is faster/slower than accessing the traceback blocks in the non-antidiagonal way.
    for (unsigned diag = 1; diag < numDiag; ++diag)
    {
        // Set iterator.
        blockVIt  = blockVBeg + _min(diag - 1, length(seqV) - 1);           // Vertical block starts at current minor diagonal or at the end of seqV.
        blockHIt  = blockHBeg + _max(0, static_cast<int>(diag) - static_cast<int>(length(seqV)));
        blockHEnd = blockHBeg + _min(diag, length(seqH));
        // Inner Loop:
        for (; blockHIt != blockHEnd; ++blockHIt, --blockVIt)
        {
            // Set the trace block.
            itTraceEnd += (length(*blockHIt) + 1) * (length(*blockVIt) + 1);
            tracebackBlock[((blockHIt - blockHBeg) * length(seqV)) + (blockVIt - blockVBeg)] = toRange(itTraceBeg, itTraceEnd);
            itTraceBeg = itTraceEnd;
        }
    }
    // Sanity check: The itTraceEnd must now be at the end of the traceHost.
    SEQAN_ASSERT_EQ(itTraceEnd, end(traceHost, Standard()));

    // ----------------------------------------------------------------------------
    // Compute the block as a parallel wave front over the minor diagonal.
    // Parallelization modeled after https://software.intel.com/en-us/node/506110

    // TODO(rrahn): Check if dpContex should be copied to each thread or if we make an array of dpContext of the size of numThreads.
    // TODO(pcostanz): this tbb construct creates a copy per thread; check ETS_key_usage_type for potential additional performance improvements
    typedef tbb::enumerable_thread_specific<DPContext<DPCell_<TScoreValue, AffineGaps>, typename TraceBitMap_::TTraceValue,
      String<DPCell_<TScoreValue, AffineGaps> >, TTraceTile> > TDpContext;
    TDpContext dpContext;


    class DagTask: public tbb::task {
    public:
      DagTask* successor[2];

    private:
      const int i, j;
      TDpContext& dpContext;
      String<TTraceTile>& tracebackBlock;
      const TSeqH& seqH;
      const TSeqV& seqV;
      TTileBuffer& tileBuffer;
      const Score<TScoreValue, TScoreSpec>& score;

    public:
      DagTask(const int i, const int j,
	      TDpContext& dpContext,
	      String<TTraceTile>& tracebackBlock,
	      const TSeqH& seqH, const TSeqV& seqV,
	      TTileBuffer& tileBuffer,
	      const Score<TScoreValue, TScoreSpec>& score) :
        i(i), j(j),
        dpContext(dpContext),
	tracebackBlock(tracebackBlock),
	seqH(seqH), seqV(seqV),
	tileBuffer(tileBuffer),
	score(score)
      {}

      task* execute() {
	auto& localDpContext = dpContext.local();
	getDpTraceMatrix(localDpContext) = tracebackBlock[(i * length(seqV)) + j];
	TDPScoutState scoutState(tileBuffer.horizontalBuffer[i], tileBuffer.verticalBuffer[j]);
	String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments. Only needed for the old interface. They are not filled.
	_computeAlignment(localDpContext, traceSegments, scoutState, seqH[i], seqV[j], score, DPBandConfig<BandOff>(), TDPProfile());

	// spawning right and downward neighbors
	for (int k=0; k<2; ++k)
	  if (DagTask* t = successor[k])
	    if (t->decrement_ref_count() == 0) spawn(*t);

	return nullptr;
      }
    };

    DagTask* graph[length(seqH)][length(seqV)];
    for (int i=length(seqH); --i>=0;) {
      for (int j=length(seqV); --j>=0;) {
	graph[i][j] = new (tbb::task::allocate_root() ) DagTask(i, j, dpContext, tracebackBlock, seqH, seqV, tileBuffer, score);
	graph[i][j]->successor[0] = i+1<length(seqH) ? graph[i+1][j] : nullptr;
	graph[i][j]->successor[1] = j+1<length(seqV) ? graph[i][j+1] : nullptr;
	graph[i][j]->set_ref_count((i>0?1:0)+(j>0?1:0));
      }
    }
    auto lastTask = graph[length(seqH)-1][length(seqV)-1];
    lastTask->increment_ref_count();
    lastTask->spawn_and_wait_for_all(*graph[0][0]);
    lastTask->execute();
    tbb::task::destroy(*lastTask);

    // TODO(rrahn): Fix to also support local and semi-global alignment!
    implParallelTrace(target, length(back(tracebackBlock)) - 1, length(seqH) - 1, length(seqV) - 1, tracebackBlock,
                      seqH, seqV, TDPProfile());
    return _scoreOfCell(back(back(tileBuffer.horizontalBuffer)).i1);
}

#endif


// ----------------------------------------------------------------------------
// Function parallelAlign()
// ----------------------------------------------------------------------------

// Main function to be called by the user.
// The interface is minimalistic here, since we plan to revise the user interface
// in the future.
template <typename TTargetH,
          typename TTargetV,
          typename TScoreValue, typename TScoreSpec>
inline TScoreValue
parallelAlign(TTargetH & targetH,
              TTargetV & targetV,
              Score<TScoreValue, TScoreSpec> const & score)
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
    unsigned const BLOCK_SIZE = 80;

    // Blocks are implemented as array of Ranges (storing begin and end iterator to corresponding block).
    String<Range<TIterH> >  seqHBlocks;
    String<Range<TIterV> >  seqVBlocks;

    // ----------------------------------------------------------------------------
    // Initialization of block structure.

    // Partition sequences into blocks
//    std::cout << "Blocks: " << ((length(seqH) << 1) - 1) / BLOCK_SIZE << " == " << (length(seqH) + BLOCK_SIZE - 1) / BLOCK_SIZE << std::endl;
    resize(seqHBlocks, (length(source(targetH)) + BLOCK_SIZE - 1) / BLOCK_SIZE, Exact());
    resize(seqVBlocks, (length(source(targetV)) + BLOCK_SIZE - 1) / BLOCK_SIZE, Exact());

    for (unsigned id = 0; id < length(seqHBlocks); ++id)
        seqHBlocks[id] = toRange(infix(source(targetH), id * BLOCK_SIZE, _min(length(source(targetH)),(id + 1) * BLOCK_SIZE)));

    for (unsigned id = 0; id < length(seqVBlocks); ++id)
        seqVBlocks[id] = toRange(infix(source(targetV), id * BLOCK_SIZE, _min(length(source(targetV)),(id + 1) * BLOCK_SIZE)));

    // The original trace back martix is allocated here.
    // The future interface should pass a context from outside so this memory block can be reused for many calls.
    typedef String<typename TraceBitMap_::TTraceValue> TTraceMatrix;

    // Partition traceback into blocks
    TTraceMatrix traceMatrix;
    resize(traceMatrix, ((length(seqHBlocks) - 1) * (BLOCK_SIZE + 1) + length(back(seqHBlocks)) + 1) *
                        ((length(seqVBlocks) - 1) * (BLOCK_SIZE + 1) + length(back(seqVBlocks)) + 1), Exact());

    String<TTraceSegment> trace;  // Stores the trace segments which are converted at the end into the final alignment.

    // The main blockwise wrapper for the dp algorithm.
    TScoreValue res = implParallelAlign(trace, traceMatrix, seqHBlocks, seqVBlocks, score);
    _adaptTraceSegmentsTo(targetH, targetV, trace);  // Convert to corresponding alignment representation.
    return res;
}

}

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_DP_PARALLEL_H_
