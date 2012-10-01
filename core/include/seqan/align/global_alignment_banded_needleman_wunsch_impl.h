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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================
// Efficient banded implementations of the Needleman-Wunsch algorithm.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_BANDED_NEEDLEMAN_WUNSCH_IMPL_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_BANDED_NEEDLEMAN_WUNSCH_IMPL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _alignBandedNeedlemanWunschTrace()
// ----------------------------------------------------------------------------

template <typename TAlign, typename TSequenceH, typename TSequenceV, typename TId, typename TTrace, typename TValPair, typename TIndexPair, typename TDiagonal>
inline void
_alignBandedNeedlemanWunschTrace(TAlign& align,
                                 TSequenceH const & seqH,
                                 TSequenceV const & seqV,
                                 TId id1,
                                 TId id2,
                                 TTrace const& trace,
                                 TValPair const& overallMaxValue,
                                 TIndexPair const& overallMaxIndex,
                                 TDiagonal const diagL,
                                 TDiagonal const diagU)
{
	typedef typename Size<TTrace>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;

	// Traceback values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	// Initialization
	TSequenceH const& str1 = seqH;
	TSequenceV const& str2 = seqV;
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize lo_row = (diagU <= 0) ? -1 * diagU : 0;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	
	//// Debug stuff
	//TColumn originalMat;
	//resize(originalMat, len1 * len2);
	//TSize count = 0;
	//for(TSize i=0; i<len2; ++i) {
	//	for(TSize j=0; j<len1; ++j) {
	//		value(originalMat, i * len1 + j) = count;
	//		std::cout << count << ',';
	//		++count;
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	// Start the trace from the cell with the max value
	TSize row = (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxIndex[0] : overallMaxIndex[2];
	TSize col = (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxIndex[1] : overallMaxIndex[3];

	// Handle the skipped sequence parts
	TSize actualRow = row + lo_row;
	TSize actualCol = col + diagL + actualRow;
	if (actualCol + 1 < len1) _alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, (len1 - (actualCol + 1)),  Horizontal);
	if (actualRow + 1 < len2) _alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, (len2 - (actualRow + 1)),  Vertical);

	if ((actualRow != 0) && (actualCol != 0)) {
		// Find initial direction
		TTraceValue tv = trace[row * diagonalWidth + col];
		if (tv == Horizontal) --col;
		else if (tv == Vertical) {--row; ++col;} 
		else --row;
	
		// Walk until we hit a border
		TSize seqLen = 1;
		TTraceValue newTv = tv;
		while(true) {
			actualRow = row + lo_row;
			actualCol = col + diagL + actualRow;
			newTv = trace[row * diagonalWidth + col];

			// Check if we hit a border
			if ((actualRow == 0) || (actualCol == 0)) break;
			else {
				//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl; 
				if (tv == Diagonal) {
					if (newTv == Horizontal) {
						_alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, seqLen, tv);
						--col; seqLen = 1;
					} else if (newTv == Vertical) {
						_alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, seqLen, tv);
						--row; ++col; seqLen = 1;
					} else {
						--row; ++seqLen;
					}
				} else {
					if (tv == Horizontal) { 
						if (newTv == Diagonal) {
							_alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, seqLen, tv);
							--row; seqLen = 1;
						} else if (newTv == Vertical) {
							_alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, seqLen, tv);
							--row; ++col; seqLen = 1;
						} else {
							--col; ++seqLen;
						}
					} else { 
						if (newTv == Diagonal) {
							_alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, seqLen, tv);
							--row; seqLen = 1;
						} else if (newTv == Horizontal) {
							_alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, seqLen, tv);
							--col; seqLen = 1;
						} else {
							--row; ++col; ++seqLen;
						}
					}
				}
				tv = newTv;
			}
		}
	
		// Align left overs
		if (seqLen) _alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, seqLen, tv);
	}

	// Handle the remaining sequence
	if (actualCol != 0) _alignTracePrint(align, seqH, seqV, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol,  Horizontal);
	else if (actualRow != 0) _alignTracePrint(align, seqH, seqV, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow,  Vertical);

}

// ----------------------------------------------------------------------------
// Function _alignBandedNeedlemanWunsch()
// ----------------------------------------------------------------------------

template <typename TTrace, typename TSequenceH, typename TSequenceV, typename TScore, typename TValPair, typename TIndexPair, typename TDiagonal, typename TAlignConfig>
inline typename Value<TScore>::Type
_alignBandedNeedlemanWunsch(TTrace& trace,
                            TSequenceH const & seqH,
                            TSequenceV const & seqV,
                            TScore const & sc,
                            TValPair& overallMaxValue,
                            TIndexPair& overallMaxIndex,
                            TDiagonal diagL,
                            TDiagonal diagU,
                            TAlignConfig const)
{
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TTrace>::Type TSize;

	// Initialization
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	TSequenceH const& str1 = seqH;
	TSequenceV const& str2 = seqV;
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	TSize hi_diag = diagonalWidth;
	TSize lo_diag = 0;
	if (diagL > 0) lo_diag = 0;
	else lo_diag = (diagU < 0) ? hi_diag : (TSize) (1 - diagL); 
	TSize lo_row = (diagU <= 0) ? -diagU : 0;
	TSize hi_row = len2;
	if (len1 - diagL < hi_row) hi_row = len1 - diagL;
	TSize height = hi_row - lo_row;

	typedef String<TScoreValue> TRow;
	TRow mat;
	resize(mat, diagonalWidth);
	resize(trace, height * diagonalWidth);
	overallMaxValue[0] = MinValue<TScoreValue>::VALUE;
	overallMaxValue[1] = MinValue<TScoreValue>::VALUE;
	overallMaxIndex[0] = diagonalWidth;	overallMaxIndex[1] = height;
	overallMaxIndex[2] = diagonalWidth;	overallMaxIndex[3] = height;
	
	//// Debug stuff
	//String<TScoreValue> originalMat;
	//resize(originalMat, len1 * len2);
	//TSize count = 0;
	//for(TSize i=0; i<len2; ++i) {
	//	for(TSize j=0; j<len1; ++j) {
	//		value(originalMat, i * len1 + j) = count;
	//		std::cerr << count << ',';
	//		++count;
	//	}
	//	std::cerr << std::endl;
	//}
	//std::cerr << std::endl;

    // TODO(holtgrew): Surely, *linear* gap costs are meant in the comment below? The code below appears to be converted over from Gotoh, it might make sense to replace scoreGapExtend* with scoreGap to reflect that this is NW DP.
    
	// Classical DP with affine gap costs
	typedef typename Iterator<TRow, Standard>::Type TRowIter;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TSize actualCol = 0;
	TSize actualRow = 0;
	TScoreValue verti_val = 0;
	TScoreValue hori_val = 0;
	for(TSize row = 0; row < height; ++row) {
		actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if (row + lo_row >= len1 - diagU) --hi_diag;
		TTraceIter traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
		TRowIter matIt = begin(mat, Standard()) + lo_diag;
		hori_val = MinValue<TScoreValue>::VALUE;
		for(TSize col = lo_diag; col<hi_diag; ++col, ++matIt, ++traceIt) {
			actualCol = col + diagL + actualRow;
			//std::cerr << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl;

			if ((actualRow != 0) && (actualCol != 0)) {
				// Get the new maximum for mat
				*matIt += score(const_cast<TScore&>(sc), ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
				*traceIt = Diagonal;
				if ((verti_val = (col < diagonalWidth - 1) ? *(matIt+1) + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : MinValue<TScoreValue>::VALUE) > *matIt) {
					*matIt = verti_val;
					*traceIt = Vertical;
				}						
				if ((hori_val = (col > 0) ? hori_val + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : MinValue<TScoreValue>::VALUE) > *matIt) {
					*matIt = hori_val;
					*traceIt = Horizontal;
				}
				hori_val = *matIt;				
			} else {			
				// Usual initialization for first row and column
				if (actualRow == 0) _initFirstRow(TAlignConfig(), *matIt, (TScoreValue) actualCol * scoreGapExtendHorizontal(sc, ((int) actualCol - 1), -1, str1, str2));
				else {
					_initFirstColumn(TAlignConfig(), *matIt, (TScoreValue) actualRow * scoreGapExtendVertical(sc, -1, ((int) actualRow - 1), str1, str2));
					hori_val = *matIt;
				}
			}

			// Store the maximum
			if (actualCol == len1 - 1) _lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
			if (actualRow == len2 - 1) _lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
			//std::cerr << row << ',' << col << ':' << *matIt << std::endl;
		}
	}
	return (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxValue[0] : overallMaxValue[1];
}


// ----------------------------------------------------------------------------
// Function _globalAlignment()
// ----------------------------------------------------------------------------

// This function is called by the public interface functions.

template <typename TAlign, typename TSequenceH, typename TSequenceV,
          typename TSequenceId,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
inline TScoreValue
_globalAlignment(TAlign & align,
				 TSequenceH const & sequenceH,
                 TSequenceV const & sequenceV,
                 TSequenceId idH,
                 TSequenceId idV,
				 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
				 int lowerDiag,
				 int upperDiag,
				 NeedlemanWunsch const & /*algorithmTag*/)
{
	typedef typename Size<TSequenceH>::Type TSize;
	  
	// Maximum value
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[4];
	
	// Create the trace
	String<TraceBack> trace;
	TScoreValue maxScore = _alignBandedNeedlemanWunsch(trace, sequenceH, sequenceV, scoringScheme, overallMaxValue, overallMaxIndex, lowerDiag, upperDiag, alignConfig);
	
	// Follow the trace and create the graph
	_alignBandedNeedlemanWunschTrace(align, sequenceH, sequenceV, idH, idV, trace, overallMaxValue, overallMaxIndex, lowerDiag, upperDiag);

	return maxScore;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_BANDED_NEEDLEMAN_WUNSCH_IMPL_H_
