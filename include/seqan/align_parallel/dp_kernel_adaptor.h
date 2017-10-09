// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_KERNEL_ADAPTOR_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_KERNEL_ADAPTOR_H_

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
          typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TCellDescriptor,
          typename TAlgo, typename TGapCosts, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             TRecursionCellTuple recursionCells,
             TSequenceHValue const & /*seqHVal*/,
             TSequenceVValue const & /*seqVVal*/,
             TScoringScheme const & /*scoringScheme*/,
             MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
             TCellDescriptor const &,   // One of FirstCell, InnerCell or LastCell.
             DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel> const &)
{
    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInitialColumn, FullColumn> >   TMetaColumn;

    std::get<0>(recursionCells) = (*scout.state.ptrVerBuffer)[scout.verticalPos].i1;
    assignValue(traceMatrixNavigator, (*scout.state.ptrVerBuffer)[scout.verticalPos].i2);

    if (TrackingEnabled_<TMetaColumn, TCellDescriptor>::VALUE)
    {
        _scoutBestScore(scout, std::get<0>(recursionCells), traceMatrixNavigator, False(), False());
    }
}

// ----------------------------------------------------------------------------
// Function _computeCell(); InnerCol; FirstCell
// ----------------------------------------------------------------------------

// Horizontal initialization values are copied from buffer for all first cells.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TGapCosts, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             TRecursionCellTuple recursionCells,
             TSequenceHValue const & /*seqHVal*/,
             TSequenceVValue const & /*seqVVal*/,
             TScoringScheme const & /*scoringScheme*/,
             MetaColumnDescriptor<DPInnerColumn, FullColumn> const &,
             FirstCell const &,   // One of FirstCell, InnerCell or LastCell.
             DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel> const &)
{
    std::get<0>(recursionCells) = (*scout.state.ptrHorBuffer)[scout.horizontalPos - 1].i1;
    assignValue(traceMatrixNavigator, (*scout.state.ptrHorBuffer)[scout.horizontalPos - 1].i2);
}

// ----------------------------------------------------------------------------
// Function _computeCell(); InnerCol; LastCell
// ----------------------------------------------------------------------------

// Values of last call are copied into the horizontal buffer for initializing next tile below.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TGapCosts, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             TRecursionCellTuple recursionCells,
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
                _computeScore(recursionCells, seqHVal, seqVVal,
                              scoringScheme, typename RecursionDirection_<TMetaColumn, LastCell>::Type(),
                              TDPProfile()));
    // Copy values into horizontal buffer for the tile below this tile in vertical direction.
    (*scout.state.ptrHorBuffer)[scout.horizontalPos - 1].i1 = std::get<0>(recursionCells);
    if (IsTracebackEnabled_<TTraceConfig>::VALUE)
    {
        (*scout.state.ptrHorBuffer)[scout.horizontalPos - 1].i2 = value(traceMatrixNavigator);
    }

    if (TrackingEnabled_<TMetaColumn, LastCell>::VALUE)
    {
        _scoutBestScore(scout, std::get<0>(recursionCells), traceMatrixNavigator, False(), True());
    }
}


// ----------------------------------------------------------------------------
// Function _computeCell(); FinalCol; FirstCell
// ----------------------------------------------------------------------------

// Horizontal initialization values are copied from buffer for all first cells.
// Vertical buffer is filled with value.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TGapCosts, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             TRecursionCellTuple recursionCells,
             TSequenceHValue const & /*seqHVal*/,
             TSequenceVValue const & /*seqVVal*/,
             TScoringScheme const & /*scoringScheme*/,
             MetaColumnDescriptor<DPFinalColumn, FullColumn> const &,
             FirstCell const &,   // One of FirstCell, InnerCell or LastCell.
             DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel> const &)
{
    typedef DPProfile_<TAlgo, TGapCosts, TTraceConfig, Parallel>                            TDPProfile;
    typedef DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, FullColumn> >     TMetaColumn;

    std::get<0>(recursionCells) =
        front(*scout.state.ptrVerBuffer).i1 = (*scout.state.ptrHorBuffer)[scout.horizontalPos - 1].i1;  // Copy horizontal buffer value in active cell and in
    assignValue(traceMatrixNavigator, (*scout.state.ptrHorBuffer)[scout.horizontalPos - 1].i2);
    if (IsTracebackEnabled_<TTraceConfig>::VALUE)
    {
        front(*scout.state.ptrVerBuffer).i2 = value(traceMatrixNavigator);   // Store trace value in vertical buffer.
    }

    if (TrackingEnabled_<TMetaColumn, FirstCell>::VALUE)
    {
        _scoutBestScore(scout, std::get<0>(recursionCells), traceMatrixNavigator, True(), False());
    }
}

// ----------------------------------------------------------------------------
// Function _computeCell(); FinalCol, InnerCell;
// ----------------------------------------------------------------------------

// Stores computed values in vertical buffer for initializing next tile right of the current.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TGapCosts, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             TRecursionCellTuple recursionCells,
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
                _computeScore(recursionCells, seqHVal, seqVVal,
                              scoringScheme, typename RecursionDirection_<TMetaColumn, InnerCell>::Type(),
                              TDPProfile()));
    // Store values in vertical buffer.
    (*scout.state.ptrVerBuffer)[scout.verticalPos].i1 = std::get<0>(recursionCells);
    if (IsTracebackEnabled_<TTraceConfig>::VALUE)
    {
        (*scout.state.ptrVerBuffer)[scout.verticalPos].i2 = value(traceMatrixNavigator);
    }

    if (TrackingEnabled_<TMetaColumn, InnerCell>::VALUE)
    {
        _scoutBestScore(scout, std::get<0>(recursionCells), traceMatrixNavigator, True(), False());
    }
}

// ----------------------------------------------------------------------------
// Function _computeCell(); FinalCol, LastCell;
// ----------------------------------------------------------------------------

// Stores computed values in vertical buffer for initializing next tile right of the current.
// Stores computed values in horizontal buffer for initializing next tile below.
template <typename TDPScout,
          typename TTraceMatrixNavigator,
          typename TRecursionCellTuple,
          typename TSequenceHValue,
          typename TSequenceVValue,
          typename TScoringScheme,
          typename TAlgo, typename TGapCosts, typename TTraceConfig>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             TRecursionCellTuple recursionCells,
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
                _computeScore(recursionCells, seqHVal, seqVVal,
                              scoringScheme, typename RecursionDirection_<TMetaColumn, LastCell>::Type(),
                              TDPProfile()));
    // Store values in vertical and horizontal buffer.
    (*scout.state.ptrHorBuffer)[scout.horizontalPos - 1].i1 = (*scout.state.ptrVerBuffer)[scout.verticalPos].i1 = std::get<0>(recursionCells);
    if (IsTracebackEnabled_<TTraceConfig>::VALUE)
    {
        (*scout.state.ptrHorBuffer)[scout.horizontalPos - 1].i2 =
            (*scout.state.ptrVerBuffer)[scout.verticalPos].i2 = value(traceMatrixNavigator);
    }

    if (TrackingEnabled_<TMetaColumn, LastCell>::VALUE)
    {
        _scoutBestScore(scout, std::get<0>(recursionCells), traceMatrixNavigator, True(), True());
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_KERNEL_ADAPTOR_H_
