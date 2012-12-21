// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Implements the traceback algorithm.
// ==========================================================================

#ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_ALIGN_DP_TRACEBACK_IMPL_H_
#define SANDBOX_RMAERKER_INCLUDE_SEQAN_ALIGN_DP_TRACEBACK_IMPL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Function _doTraceback()
// ----------------------------------------------------------------------------
template <typename TTarget, typename TDPTraceMatrixNavigator, typename TTraceValue, typename TSize, typename TPosition,
          typename TDPBand, typename TGapCosts>
inline void
_doTraceback(TTarget & target,
             TDPTraceMatrixNavigator & matrixNavigator,
             TTraceValue & traceValue,
             TTraceValue & lastTraceValue,
             TSize & fragmentLength,
             TPosition & seqHPos,
             TPosition & seqVPos,
             TDPBand const & band,
             TGapCosts const &)
{
    if (traceValue & TraceBitMap_::DIAGONAL)
    {
        if (!(lastTraceValue & TraceBitMap_::DIAGONAL)) // the old trace value was not diagonal
        {
            _recordSegment(target, seqHPos, seqVPos, fragmentLength, lastTraceValue);

            lastTraceValue = TraceBitMap_::DIAGONAL;
            fragmentLength = 0;
        }
        --seqHPos;
        --seqVPos;
        _traceDiagonal(matrixNavigator, band);
        ++fragmentLength;
        traceValue = value(matrixNavigator);
    }  // In case of Gotoh we prefer the longest possible way in this direction.
    else if (traceValue & TraceBitMap_::VERTICAL)
    {
        if (!(lastTraceValue & TraceBitMap_::VERTICAL)) // the old trace value was not diagonal
        {
            _recordSegment(target, seqHPos, seqVPos, fragmentLength, lastTraceValue);

            lastTraceValue = TraceBitMap_::VERTICAL;
            fragmentLength = 0;
        }
        // We are in a vertical gap. So continue after we reach the end of the vertical gap.
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            while (!(traceValue & TraceBitMap_::VERTICAL_OPEN) && (seqVPos != 1))
            {
                --seqVPos;
                _traceVertical(matrixNavigator, band);
                ++fragmentLength;
                traceValue = value(matrixNavigator);
            }
            // We have to ensure, that we do not continue in vertical direction if we reached a vertical_open sign.
            --seqVPos;
            _traceVertical(matrixNavigator, band);
            ++fragmentLength;
            // Forbid continuing in vertical direction.
            traceValue = value(matrixNavigator) & TraceBitMap_::NO_VERTICAL_TRACEBACK;
            SEQAN_ASSERT_NEQ(traceValue, +TraceBitMap_::NONE);
        }
        else
        {
            --seqVPos;
            _traceVertical(matrixNavigator, band);
            ++fragmentLength;
            traceValue = value(matrixNavigator);
        }
    }
    else if (traceValue & TraceBitMap_::MAX_FROM_VERTICAL_MATRIX)
    {
        if (!(lastTraceValue & TraceBitMap_::VERTICAL)) // the old trace value was not diagonal
        {
            _recordSegment(target, seqHPos, seqVPos, fragmentLength, lastTraceValue);
            lastTraceValue = TraceBitMap_::VERTICAL;
            fragmentLength = 0;
        }
        --seqVPos;
        _traceVertical(matrixNavigator, band);
        ++fragmentLength;
        // Forbid continuing in vertical direction.
        traceValue = value(matrixNavigator) & TraceBitMap_::NO_VERTICAL_TRACEBACK;
        SEQAN_ASSERT_NEQ(traceValue, +TraceBitMap_::NONE);
    }  // In case of Gotoh we prefer the longest possible way in this direction.
    else if (traceValue & TraceBitMap_::HORIZONTAL)
    {
        if (!(lastTraceValue & TraceBitMap_::HORIZONTAL)) // the old trace value was not diagonal
        {
            _recordSegment(target, seqHPos, seqVPos, fragmentLength, lastTraceValue);

            lastTraceValue = TraceBitMap_::HORIZONTAL;
            fragmentLength = 0;
        }
        if (IsSameType<TGapCosts, AffineGaps>::VALUE)
        {
            while (!(traceValue & TraceBitMap_::HORIZONTAL_OPEN) && (seqHPos != 1))
            {
                --seqHPos;
                _traceHorizontal(matrixNavigator, band);
                ++fragmentLength;
                traceValue = value(matrixNavigator);
            }
            --seqHPos;
            _traceHorizontal(matrixNavigator, band);
            ++fragmentLength;
            // Forbid continuing in horizontal direction.
            traceValue = value(matrixNavigator) & TraceBitMap_::NO_HORIZONTAL_TRACEBACK;
            SEQAN_ASSERT_NEQ(traceValue, +TraceBitMap_::NONE);
        }
        else
        {
            --seqHPos;
            _traceHorizontal(matrixNavigator, band);
            ++fragmentLength;
            traceValue = value(matrixNavigator);
        }
    }
    else if (traceValue & TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX)
    {
        if (!(lastTraceValue & TraceBitMap_::HORIZONTAL)) // the old trace value was not diagonal
        {
            _recordSegment(target, seqHPos, seqVPos, fragmentLength, lastTraceValue);
            lastTraceValue = TraceBitMap_::HORIZONTAL;
            fragmentLength = 0;
        }

        --seqHPos;
        _traceHorizontal(matrixNavigator, band);
        ++fragmentLength;
        // Forbid continuing in horizontal direction.
        traceValue = value(matrixNavigator) & TraceBitMap_::NO_HORIZONTAL_TRACEBACK;
        SEQAN_ASSERT_NEQ(traceValue, +TraceBitMap_::NONE);
    }
    else // the trace back is either NONE or something else
    {
        if (traceValue == TraceBitMap_::NONE)
        {
            return;
        }
        SEQAN_ASSERT_FAIL("Reached undefined traceback value!");
    }
}

// ----------------------------------------------------------------------------
// Function _computeTraceback()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TDPTraceMatrixNavigator, typename THostPosition, typename TSequenceH, typename TSequenceV,
          typename TBandFlag, typename TAlgorithm, typename TGapCosts, typename TTracebackSpec>
void _computeTraceback(TTarget & target,
                       TDPTraceMatrixNavigator & matrixNavigator,
                       THostPosition const & hostPosition,
                       TSequenceH const & seqH,
                       TSequenceV const & seqV,
                       DPBand_<TBandFlag> const & band,
                       DPProfile_<TAlgorithm, TGapCosts, TTracebackSpec> const &)
{
    typedef typename Container<TDPTraceMatrixNavigator>::Type TContainer;
    typedef typename Size<TContainer>::Type TSize;
    typedef typename Position<TContainer>::Type TPosition;
    typedef typename MakeSigned<TPosition>::Type TSignedPosition;
    typedef typename MakeSigned<TSize>::Type TSignedSize;
    typedef typename TraceBitMap_::TTraceValue TTraceValue;

    if (IsSameType<TTracebackSpec, TracebackOff>::VALUE)
        return;

    TSignedSize seqHSize = length(seqH);
    TSignedSize seqVSize = length(seqV);

    TSize fragmentLength = 0;
    // compute the sequence position:
    setToPosition(matrixNavigator, hostPosition);
    TSignedPosition currColumn = coordinate(matrixNavigator, +DPMatrixDimension_::HORIZONTAL);
    TSignedPosition currRow = coordinate(matrixNavigator, +DPMatrixDimension_::VERTICAL);

    SEQAN_ASSERT_LEQ(currColumn, seqHSize);
    SEQAN_ASSERT_LEQ(currRow, seqVSize);

    TTraceValue traceValue = value(matrixNavigator);
    TTraceValue lastTraceValue = TraceBitMap_::NONE;

    // we need to change the tb value here.
    if (IsSameType<TGapCosts, LinearGaps>::VALUE)
    {
        if (traceValue & TraceBitMap_::DIAGONAL)
            lastTraceValue = TraceBitMap_::DIAGONAL;
        else if (traceValue & (TraceBitMap_::VERTICAL | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX))
            lastTraceValue = TraceBitMap_::VERTICAL;
        else if (traceValue & (TraceBitMap_::HORIZONTAL | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX))
            lastTraceValue = TraceBitMap_::HORIZONTAL;
    }
    else
    {
        lastTraceValue = TraceBitMap_::DIAGONAL;
        if (traceValue & (TraceBitMap_::VERTICAL | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX))
        {
            traceValue &= (TraceBitMap_::VERTICAL | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX);
            lastTraceValue = TraceBitMap_::VERTICAL;
        }
        else if (traceValue & (TraceBitMap_::HORIZONTAL | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX))
        {
            traceValue &= (TraceBitMap_::HORIZONTAL | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX);
            lastTraceValue = TraceBitMap_::HORIZONTAL;
        }
    }

    // Correct the coordinates in banded case.
    TSignedPosition firstStop = 0;
    TSignedPosition secondStop = 0;

    if (IsSameType<TBandFlag, BandOn>::VALUE)
    {
        if (lowerDiagonal(band) >= 0)
            currColumn += lowerDiagonal(band);
        if (currColumn > upperDiagonal(band))
            currRow += currColumn - upperDiagonal(band);

        firstStop = _min(seqHSize, _max(0, upperDiagonal(band)));
        secondStop = _min(seqHSize, _max(0, static_cast<int>(seqVSize) + lowerDiagonal(band)));
    }

    if (IsGlobalAlignment_<TAlgorithm>::VALUE)
    {
        if (currRow != seqVSize)
            _recordSegment(target, seqHSize, currRow, seqVSize - currRow, +TraceBitMap_::VERTICAL);
        if (currColumn != seqHSize)
            _recordSegment(target, currColumn, currRow, seqHSize - currColumn, +TraceBitMap_::HORIZONTAL);  // if positions lie in middle of matrix, we go with manhatten distance to this point.
    }

    if (firstStop > secondStop)  // In this case we have to divide the traceback into three parts, which have different distances to the previous cells.
    {
        while (currColumn > firstStop && currRow > 0 && traceValue != TraceBitMap_::NONE)
            _doTraceback(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, currColumn, currRow, band, TGapCosts());

        while (currColumn > secondStop && currRow > 0 && traceValue != TraceBitMap_::NONE)
            _doTraceback(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, currColumn, currRow, DPBand_<BandOff>(), TGapCosts());

        while (currColumn > 0 && currRow > 0 && traceValue != TraceBitMap_::NONE)
            _doTraceback(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, currColumn, currRow, band, TGapCosts());

    }
    else  // This is the standard case for unbanded and small banded alignments.
        while (currColumn > 0 && currRow > 0 && traceValue != TraceBitMap_::NONE)
            _doTraceback(target, matrixNavigator, traceValue, lastTraceValue, fragmentLength, currColumn, currRow, band, TGapCosts());

    // Record last detected fragment.
    _recordSegment(target, currColumn, currRow, fragmentLength, lastTraceValue);
    if (IsGlobalAlignment_<TAlgorithm>::VALUE)
    {
        // Record trailing gaps if any.
        if (currRow != 0)
            _recordSegment(target, 0, 0, currRow, +TraceBitMap_::VERTICAL);
        if (currColumn != 0)
            _recordSegment(target, 0, 0, currColumn, +TraceBitMap_::HORIZONTAL);
    }

}

}  // namespace seqan

#endif  // #ifndef SANDBOX_RMAERKER_INCLUDE_SEQAN_ALIGN_DP_TRACEBACK_IMPL_H_
