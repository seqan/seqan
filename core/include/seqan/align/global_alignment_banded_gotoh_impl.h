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
// Efficient banded implementations of Gotoh's algorithm.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_BANDED_GOTOH_IMPL_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_BANDED_GOTOH_IMPL_H_

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
// Function _alignBandedGotohTrace()
// ----------------------------------------------------------------------------

template <typename TAlign, typename TSequenceH, typename TSequenceV, typename TId, typename TTrace, typename TValPair, typename TIndexPair, typename TDiagonal>
inline void
_alignBandedGotohTrace(TAlign& align,
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
	typedef unsigned char TTraceValue;

	// Gotoh back-trace values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	// Initialization
	TSequenceH const& str1 = seqH;
	TSequenceV const& str2 = seqV;
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize lo_row = 0;
	if (diagU <= 0) lo_row = -1 * diagU;
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

	// Walk until we hit a border
	TTraceValue tv = (trace[row * diagonalWidth + col] & 3);
	TTraceValue oldTraceValue = tv;
	TSize seqLen = 0;
	bool hitBorder = false;
	do {
		actualRow = row + lo_row;
		actualCol = col + diagL + actualRow;

		// Direction changed, so make aligned segments
		if (oldTraceValue != tv) {
			_alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, seqLen, oldTraceValue);
			seqLen = 0;
		}

		// Check if we hit a border
		if ((actualRow == 0) || (actualCol == 0)) hitBorder = true;
		else {
			//std::cout << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl; 
			
			// Last value was diagonal
			if (tv == Diagonal) {
				oldTraceValue = Diagonal;
				tv = (trace[row * diagonalWidth + col] & 3);
				if (tv == Diagonal) {--row; ++seqLen;}
			} else if (tv == Horizontal) { // Last value was horizontal
				oldTraceValue = Horizontal;
				if ((col > 0) && ((trace[row * diagonalWidth + col] >> 2) & 1)) tv = Diagonal;
				--col; ++seqLen;
			} else { // Vertical
				oldTraceValue = Vertical;
				if ((col < diagonalWidth - 1) && ((trace[row * diagonalWidth + col] >> 3) & 1)) tv = Diagonal;
				--row; ++col; ++seqLen;
			}	
		}
	} while(!hitBorder);
	
	// Align left overs
	if (seqLen) _alignTracePrint(align, seqH, seqV, id1, actualCol, id2, actualRow, seqLen, tv);

	// Handle the remaining sequence
	if (actualCol != 0) _alignTracePrint(align, seqH, seqV, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol,  Horizontal);
	else if (actualRow != 0) _alignTracePrint(align, seqH, seqV, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow,  Vertical);
}

// ----------------------------------------------------------------------------
// Function _alignBandedGotoh()
// ----------------------------------------------------------------------------

template <typename TTrace, typename TSequenceH, typename TSequenceV, typename TScore, typename TValPair, typename TIndexPair, typename TDiagonal, typename TAlignConfig>
inline typename Value<TScore>::Type
_alignBandedGotoh(TTrace& trace,
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
	typedef String<TScoreValue> TColumn;
	TColumn mat;
	TColumn vertical;
	resize(mat, diagonalWidth);
	resize(vertical, diagonalWidth);
	resize(trace, height * diagonalWidth);
	overallMaxValue[0] = MinValue<TScoreValue>::VALUE;
	overallMaxValue[1] = MinValue<TScoreValue>::VALUE;
	overallMaxIndex[0] = diagonalWidth;
	overallMaxIndex[1] = height;
	overallMaxIndex[2] = diagonalWidth;
	overallMaxIndex[3] = height;
	
	//// Debug stuff
	//TColumn originalMat;
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

	// Classical DP with affine gap costs
	typedef typename Iterator<TColumn, Standard>::Type TColIter;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceValue tvMat = 0;
	TSize actualRow = 0;
	TSize actualCol = 0;
	TScoreValue a = 0;
	TScoreValue b = 0;
	TScoreValue hori_val = 0;
	for(TSize row = 0; row < height; ++row) {
		actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if (row + lo_row >= len1 - diagU) --hi_diag;
		TTraceIter traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
		TColIter vertIt = begin(vertical, Standard()) + lo_diag;
		TColIter matIt = begin(mat, Standard()) + lo_diag;
		hori_val = MinValue<TScoreValue>::VALUE;
		for(TSize col = lo_diag; col<hi_diag; ++col, ++vertIt, ++matIt, ++traceIt) {
			actualCol = col + diagL + actualRow;
			//std::cerr << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl;

			if ((actualRow != 0) && (actualCol != 0)) {
				// Get the new maximum for vertical
				*traceIt = 0;
				if (col < diagonalWidth - 1) {
					a = *(matIt + 1) + scoreGapOpenVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
					b = (*(vertIt + 1) != MinValue<TScoreValue>::VALUE) ? *(vertIt + 1) + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : MinValue<TScoreValue>::VALUE;
					if (a > b) {*vertIt = a; *traceIt = 1;}
					else *vertIt = b;
				} else *vertIt = MinValue<TScoreValue>::VALUE;

				// Get the new maximum for horizontal
				*traceIt <<= 1;
				if (col > 0) {
					a = *(matIt -1 ) + scoreGapOpenHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
					b = (hori_val != MinValue<TScoreValue>::VALUE) ? hori_val + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : MinValue<TScoreValue>::VALUE;
					if (a > b) {hori_val = a; *traceIt |= 1;}
					else hori_val = b;
				} else hori_val = MinValue<TScoreValue>::VALUE;

				*traceIt <<= 2;
				// Get the new maximum for mat
				*matIt += score(const_cast<TScore&>(sc), actualCol-1, actualRow-1, str1, str2);
				tvMat = Diagonal;
				if (*vertIt > *matIt) {
					*matIt = *vertIt;
					tvMat = Vertical;
				}
				if (hori_val > *matIt) {
					*matIt = hori_val;
					tvMat = Horizontal;
				}
				*traceIt |= tvMat;
			} else {
				// Usual initialization for first row and column
				if (actualRow == 0) {
					if (actualCol != 0) {
                        *traceIt = Horizontal;
						_initFirstRow(TAlignConfig(), *matIt, scoreGapOpenHorizontal(sc, 0, -1, str1, str2) + (actualCol - 1) * scoreGapExtendHorizontal(sc, ((int) actualCol - 1), -1, str1, str2));
						*vertIt = *matIt + scoreGapOpenVertical(sc, ((int) actualCol - 1), 0, str1, str2) - scoreGapExtendVertical(sc, ((int) actualCol - 1), 0, str1, str2);
						hori_val = MinValue<TScoreValue>::VALUE;
					} else {
						*matIt = 0;
						*vertIt = MinValue<TScoreValue>::VALUE;
						hori_val = MinValue<TScoreValue>::VALUE;
					}
				} else {
                    *traceIt = Vertical;
					_initFirstColumn(TAlignConfig(), *matIt, scoreGapOpenVertical(sc, -1, 0, str1, str2) + (actualRow - 1) * scoreGapExtendVertical(sc, -1, ((int) actualRow - 1), str1, str2));
					hori_val = *matIt + scoreGapOpenHorizontal(sc, 0, ((int) actualRow - 1), str1, str2) - scoreGapExtendHorizontal(sc, 0, ((int) actualRow - 1), str1, str2);
					*vertIt = MinValue<TScoreValue>::VALUE;
				}
			}

			// Store the maximum
			if (actualCol == len1 - 1) _lastColumn(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
			if (actualRow == len2 - 1) _lastRow(TAlignConfig(), overallMaxValue, overallMaxIndex, *matIt, row, col);
			//std::cerr << row << ',' << col << ':' << value(mat, row * diagonalWidth + col) << std::endl;
			//std::cerr << row << ',' << col << ':' << value(horizontal, row * diagonalWidth + col) << std::endl;
			//std::cerr << row << ',' << col << ':' << value(vertical, row * diagonalWidth + col) << std::endl;
		}
	}
	return (overallMaxValue[0] > overallMaxValue[1]) ? overallMaxValue[0] : overallMaxValue[1];
}

// ----------------------------------------------------------------------------
// Function _globalAlignment()
// ----------------------------------------------------------------------------

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
				 Gotoh const & /*algorithmTag*/)
{
	typedef typename Size<TSequenceH>::Type TSize;
	  
	// Maximum value
	TScoreValue overallMaxValue[2];
	TSize overallMaxIndex[4];

	// Create the trace
	String<unsigned char> trace;
	TScoreValue maxScore = _alignBandedGotoh(trace, sequenceH, sequenceV, scoringScheme, overallMaxValue, overallMaxIndex, lowerDiag, upperDiag, alignConfig);
	
	// Follow the trace and create the graph
	_alignBandedGotohTrace(align, sequenceH, sequenceV, idH, idV, trace, overallMaxValue, overallMaxIndex, lowerDiag, upperDiag);

	return maxScore;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_BANDED_GOTOH_IMPL_H_
