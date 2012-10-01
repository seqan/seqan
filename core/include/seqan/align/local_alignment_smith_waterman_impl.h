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
// Implementation of the Smith-Waterman algorithm.
//
// We carried this over from the old graph_align module because the Waterman-
// Eggert implementation from the old align module did not allow alignment
// graphs and no fragment strings.
//
// A future implementation of Waterman-Eggert should include this and this
// header could go away.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_SMITH_WATERMAN_IMPL_H_
#define CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_SMITH_WATERMAN_IMPL_H_

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
// Function _isClumping()
// ----------------------------------------------------------------------------

template<typename TSpec, typename TSize>
inline bool
_isClumping(String<bool, TSpec> const& forbidden,
			TSize row,
			TSize col,
			TSize len2)
{
	return forbidden[(col-1) * len2 + (row-1)];
}

template<typename TSize>
inline bool
_isClumping(Nothing&,
			TSize,
			TSize,
			TSize)
{
	return false;
}

// ----------------------------------------------------------------------------
// Function _alignSmithWatermanTrace()
// ----------------------------------------------------------------------------

template <typename TAlign, typename TSequenceH, typename TSequenceV, typename TId, typename TTrace, typename TVal, typename TIndexPair, typename TForbidden>
inline void
_alignSmithWatermanTrace(TAlign& align,
                         TSequenceH const & seqH,
                         TSequenceV const & seqV,
                         TId id1,
                         TId id2,
                         TTrace const& trace,
                         TVal const initialDir,
                         TIndexPair const& indexPair,
                         TForbidden& forbidden)
{
	typedef typename Size<TTrace>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;

	// TraceBack values for Gotoh
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

	TSize len1 = indexPair[1];
	TSize len2 = indexPair[0];
	if ((indexPair[0] == 0) || (indexPair[1] == 0)) return;
	TSize numCols = length(seqH);
	TSize numRowsOrig = length(seqV);
	if (len1 < numCols) _alignTracePrint(align, seqH, seqV, id1, len1, id2, len2, numCols - len1, Horizontal);
	if (len2 < numRowsOrig) _alignTracePrint(align, seqH, seqV, id1, len1, id2, len2, numRowsOrig - len2, Vertical);
	TSize numRows = (numRowsOrig >> 1) + (numRowsOrig & 1);
	
	

	// Initialize everything	
	TTraceValue nextTraceValue = (len2 & 1) ? trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)] >> 4 : trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)];
	TTraceValue tv = Diagonal;
	if (initialDir == Diagonal) tv = (nextTraceValue & 3);
	else if (initialDir == Horizontal) {
		if ((nextTraceValue >> 2) & 1) _alignTracePrint(align, seqH, seqV, id1, --len1, id2, len2, (TSize) 1, Horizontal);
		else tv = Horizontal;
	} else if (initialDir == Vertical) {
		if ((nextTraceValue >> 3) & 1) _alignTracePrint(align, seqH, seqV, id1, len1, id2, --len2, (TSize) 1, Vertical);
		else tv = Vertical;
	}
	TSize segLen = 0;
	TTraceValue tvOld = tv;

	// Now follow the trace
	do {
		nextTraceValue = (len2 & 1) ? trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)] >> 4 : trace[(len1 - 1)*numRows + ((len2 - 1) >> 1)];
		if ((nextTraceValue & 3) == Stop) break;
		_setForbiddenCell(forbidden, len1, len2, numRowsOrig);	
		if (tv == Diagonal) tv = (nextTraceValue & 3);
		else if (tv == Horizontal) {
			if ((nextTraceValue >> 2) & 1) tv = Diagonal; 
			else tv =  Horizontal;
		} else if (tv == Vertical) {
			if ((nextTraceValue >> 3) & 1) tv =  Diagonal; 
			else tv =  Vertical;
		}
		if (tv == Diagonal) {
			if (tv != tvOld) {
				if (tvOld == Vertical) --len2;
				else --len1;
				_alignTracePrint(align, seqH, seqV, id1, len1, id2, len2, ++segLen, tvOld);
				tvOld = tv; segLen = 0;
			} else {
				++segLen;
				--len1; --len2;
			}
		} else if (tv == Horizontal) {
			if (tv != tvOld) {
				_alignTracePrint(align, seqH, seqV, id1, len1, id2, len2, segLen, tvOld);
				if ((nextTraceValue >> 2) & 1) {
					_alignTracePrint(align, seqH, seqV, id1, --len1, id2, len2, (TSize) 1, Horizontal);
					tv = Diagonal; segLen = 0;
				} else {
					tvOld = tv; segLen = 1;
					--len1;
				}
			} else {
				++segLen;
				--len1;
			}
		} else if (tv == Vertical) {
			if (tv != tvOld) {
				_alignTracePrint(align, seqH, seqV, id1, len1, id2, len2, segLen, tvOld);
				if ((nextTraceValue >> 3) & 1) {
					_alignTracePrint(align, seqH, seqV, id1, len1, id2, --len2, (TSize) 1, Vertical);
					tv = Diagonal; segLen = 0;
				} else {
					tvOld = tv; segLen = 1;
					--len2;
				}
			} else {
				++segLen;
				--len2;
			}
		}
	} while ((len1 != 0) && (len2 !=0));
	// Process left-overs
	if (segLen) _alignTracePrint(align, seqH, seqV, id1, len1, id2, len2, segLen, tvOld);

	// Handle the remaining sequence
	if (len1 != 0) _alignTracePrint(align, seqH, seqV, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1, Horizontal);
	if (len2 != 0) _alignTracePrint(align, seqH, seqV, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2, Vertical);
}

// ----------------------------------------------------------------------------
// Function _alignSmithWaterman()
// ----------------------------------------------------------------------------

template <typename TTrace, typename TSequenceH, typename TSequenceV, typename TScore, typename TIndexPair, typename TForbidden>
inline typename Value<TScore>::Type
_alignSmithWaterman(TTrace& trace,
                    TSequenceH const & seqH,
                    TSequenceV const & seqV,
                    TScore const & sc,
                    typename Value<TTrace>::Type& initialDir,
                    TIndexPair& indexPair,
                    TForbidden& forbidden)
{
	typedef typename Size<TTrace>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;
	
	// TraceBack values for Smith Waterman
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

	// The DP Matrix for diagonal walks
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TScoreValue> TColumn;
	TColumn mat;
	// The DP Matrix for gaps from the left
	TColumn horizontal;
	// The DP Matrix for gaps from the top
	TScoreValue vert = 0;

	typedef typename Iterator<TColumn, Standard>::Type TMatIter;

	// Initialization
	TSequenceH const& str1 = seqH;
	TSequenceV const& str2 = seqV;
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	resize(mat, (len2+1));   // One column for the diagonal matrix
	resize(horizontal, (len2+1));   // One column for the horizontal matrix
	resize(trace, len1 * ((len2 >> 1) + (len2 & 1)), 0);
	TTraceValue tvMat= 0;

	// Record the max score
	TScoreValue score_max = 0;
	indexPair[0] = 0; indexPair[1] = 0;
	initialDir = Stop;
	
	// Classical DP
	TScoreValue max_val = 0;
	TScoreValue a = 0;
	TScoreValue b = 0;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TTraceIter it = begin(trace, Standard() );
	TMatIter matIt = begin(mat, Standard() );
	TMatIter horiIt = begin(horizontal, Standard() );
	*matIt = 0;
	for(TSize row = 1; row <= len2; ++row) {
		*(++matIt) = 0;
		*(++horiIt) = scoreGapOpenHorizontal(sc, 0, row-1, str1, str2) - scoreGapExtendHorizontal(sc, 0, row-1, str1, str2);
	}
	for(TSize col = 1; col <= len1; ++col) {
		matIt = begin(mat, Standard() );
		horiIt = begin(horizontal, Standard() );
		TScoreValue diagValMat = *matIt;
		*matIt = 0;
		vert = scoreGapOpenVertical(sc, col-1, 0, str1, str2) - scoreGapExtendVertical(sc, col-1, 0, str1, str2);
		TSize row = 1;
		while(row <= len2) {
			if (_isClumping(forbidden, row, col, len2)) {
				*it <<= 3;
				*it |= Stop;
				max_val = 0;
				vert = 0;
				*(++horiIt) = 0;
				++matIt;
			} else {
				// Get the new maximum for vertical
				a = *matIt + scoreGapOpenVertical(sc, col-1, row-1, str1, str2);
				b = vert + scoreGapExtendVertical(sc, col-1, row-1, str1, str2);
				if (a > b) { vert = a; *it |= 1;} 
				else vert = b;
	
				// Get the new maximum for horizontal
				*it <<= 1;
				a = *(++matIt) + scoreGapOpenHorizontal(sc, col-1, row-1, str1, str2);
				b = *(++horiIt) + scoreGapExtendHorizontal(sc, col-1, row-1, str1, str2);
				if (a > b) {*horiIt = a; *it |= 1; } 
				else *horiIt =  b;
	
				// Get the new maximum for mat
				*it <<= 2;
				max_val = diagValMat + score(const_cast<TScore&>(sc), col-1, row-1, str1, str2);
				tvMat =  Diagonal;
				if (vert > max_val) {
					max_val = vert;
					tvMat =  Vertical;
				}
				if (*horiIt > max_val) {
					max_val = *horiIt;
					tvMat =  Horizontal;
				}
				if (0 >= max_val) {
					max_val = 0;
					tvMat =  Stop;
				}
				*it |= tvMat;
			}

			// Assign the new diagonal values
			diagValMat = *matIt;
			*matIt = max_val;

			// Record the new best score
			if (max_val > score_max) {
				indexPair[0] = row; indexPair[1] = col;
				score_max = max_val;
				initialDir = tvMat;
			}

			if (row & 1) *it <<= 1; else ++it;
			++row;
		}
		if (!(row & 1)) {*it <<= 3; ++it; }
	}

	//// Debug code
	//for(TSize i= 0; i<len2;++i) {
	//	for(TSize j= 0; j<len1;++j) {
	//		std::cout << (TSize) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << "Max score: " << best_row << ',' << best_col << ':' << score_max << " (" << (TSize) initialDir << ")" << std::endl;

	return score_max;
}

// ----------------------------------------------------------------------------
// Function _localAlignment()
// ----------------------------------------------------------------------------

template<typename TAlign, typename TSequenceH, typename TId, typename TSequenceV, typename TScore>
inline typename Value<TScore>::Type
_localAlignment(TAlign& align,
                TSequenceH const & seqH,
                TSequenceV const & seqV,
                TId idH,
                TId idV,
				TScore const& sc,
				SmithWaterman)
{
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TSequenceH>::Type TSize;
	  
	TScoreValue maxScore;
	TSize indexPair[2];
	Nothing forbidden;

	// Trace
	String<unsigned char> trace;
	unsigned char initialDir;

	// Create the trace
	maxScore = _alignSmithWaterman(trace, seqH, seqV, sc, initialDir, indexPair, forbidden);	
	// Follow the trace and create the graph
	_alignSmithWatermanTrace(align, seqH, seqV, idH, idV, trace, initialDir, indexPair, forbidden);
	
	return maxScore;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_SMITH_WATERMAN_IMPL_H_
