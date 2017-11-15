// ==========================================================================
//                    STELLAR - SwifT Exact LocaL AligneR
//                   http://www.seqan.de/projects/stellar/
// ==========================================================================
// Copyright (C) 2010-2012 by Birte Kehr
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STELLAR_EXTENSION_H
#define SEQAN_HEADER_STELLAR_EXTENSION_H

#include <seqan/seeds.h>


///////////////////////////////////////////////////////////////////////////////
// Container for storing possible end positions in extension of eps-core
template<typename TPos_>
struct ExtensionEndPosition {
	typedef TPos_			TPosition;
	typedef Pair<TPosition>	TCoordinate;

	TPosition length;
	TCoordinate coord;

	ExtensionEndPosition():
		length(0), coord(TCoordinate(0,0)) {}

	ExtensionEndPosition(TPosition len, TPosition row, TPosition col):
		length(len), coord(TCoordinate(row, col)) {}
};

///////////////////////////////////////////////////////////////////////////////
// returns true if align has a match at pos, otherwise false
template<typename TSource, typename TSize>
inline bool
isMatch(Align<TSource> const & align, TSize pos) {

    if(isGap(row(align, 0), pos)) {
        return false;
    } else if(isGap(row(align, 1), pos)) {
        return false;
    } else if(row(align, 0)[pos] != row(align, 1)[pos]) {
        return false;
    } else {
        return true;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Computes possible end positions of an eps-match in a given alignment.
template<typename TSource, typename TPos>
void
_fillGapsString(Align<TSource> const & align,
                String<Triple<TPos, TPos, TPos> > & gaps) {
    typedef Triple<TPos, TPos, TPos> TGapInfo;
    TPos totalErrors = 0;
	typename Row<Align<TSource> >::Type row0 = row(align, 0);
    TPos i = 0;
	TPos endPos = length(row0);
    TPos gapBegin = i;

    // append gap starting at beginPosition (also if its length is 0!)
    while(i < endPos && !isMatch(align, i)) {
        ++i;
        ++totalErrors;
    }
    appendValue(gaps, TGapInfo(gapBegin, i, totalErrors));

    // iterate over alignment and append gaps
    while (i < endPos) {
        // skip matches
        while(i < endPos && isMatch(align, i)) {
            ++i;
        }
        gapBegin = i;
        // skip and count mismatches/indels
        while(i < endPos && !isMatch(align, i)) {
            ++i;
            ++totalErrors;
        }
        appendValue(gaps, TGapInfo(gapBegin, i, totalErrors));
    }
    /*for(unsigned l = 0; l < length(gaps); ++l) {
        std::cout << gaps[l].i1 << "  " << gaps[l].i2 << "  " << gaps[l].i3 << std::endl;
    }*/
}

///////////////////////////////////////////////////////////////////////////////
// Checks the error rate of the fragment between end of left and start of right.
template<typename TPos, typename TFloat>
inline bool
_isEpsMatch(Triple<TPos, TPos, TPos> const & left,
           Triple<TPos, TPos, TPos> const & right,
           TFloat eps) {
    // compute mismatches/indels and length
    TPos errors = right.i3 - left.i3 - (right.i2 - right.i1);
    TPos len = right.i1 - left.i2;

    // check error rate
	double const DELTA = 0.000001;  // Small delta against floating point precision problems.
    return errors/(TFloat)(len) <= eps + DELTA;
}

///////////////////////////////////////////////////////////////////////////////
// Identifies the longest epsilon match in align and sets the view positions of
// align to start and end position of the longest epsilon match
template<typename TSource, typename TSize, typename TFloat>
bool 
longestEpsMatch(Align<TSource> & align,
                TSize matchMinLength,
                TFloat epsilon) {
    // Preprocessing: compute and store gaps and lengths
    // A gap is a triple of gap begin position, gap end position, and total number of errors in sequence from begin
    //   to end position of this gap.
    typedef typename Position<Align<TSource> >::Type TPosition;
    typedef String<Triple<TPosition, TPosition, TPosition> > TGapsString;
    TGapsString gaps;
    _fillGapsString(align, gaps);

    // Identify longest eps match by iterating over combinations of left and right positions
    typename Iterator<TGapsString >::Type rightIt = end(gaps) - 1;
    typename Iterator<TGapsString >::Type leftIt = begin(gaps);

    TPosition beginPos = 0;
    TPosition endPos = 0;
    TSize minLength = matchMinLength - 1;
    
    while ((*leftIt).i2 + minLength < (*rightIt).i1) {
        while ((*leftIt).i2 + minLength < (*rightIt).i1) {
            if(_isEpsMatch(*leftIt, *rightIt, epsilon)) {
                beginPos = (*leftIt).i2;
                endPos = (*rightIt).i1;
                minLength = endPos - beginPos;
                break;
            }
            --rightIt;
        }
        rightIt = end(gaps) - 1;
        ++leftIt;
    }

    // Set view positions to the eps-match
    TPosition viewBeginRow0 = toSourcePosition(row(align, 0), 0);
    TPosition viewBeginRow1 = toSourcePosition(row(align, 1), 0);
	setClippedEndPosition(row(align, 0), viewBeginRow0 + endPos);
	setClippedEndPosition(row(align, 1), viewBeginRow1 + endPos);
    setClippedBeginPosition(row(align, 0), viewBeginRow0 + beginPos);
    setClippedBeginPosition(row(align, 1), viewBeginRow1 + beginPos);
	// setClippedBeginPosition(row(align, 0), toSourcePosition(row(align, 0), beginPos));
	// setClippedBeginPosition(row(align, 1), toSourcePosition(row(align, 1), beginPos));
	// setBeginPosition(row(align, 0), beginPos);
	// setBeginPosition(row(align, 1), beginPos);
	// setClippedEndPosition(row(align, 0), toSourcePosition(row(align, 0), endPos));
	// setClippedEndPosition(row(align, 1), toSourcePosition(row(align, 1), endPos));
    SEQAN_ASSERT_EQ(length(row(align, 0)), length(row(align, 1)));

	if (endPos == 0 && beginPos == 0) return 1;
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Identifies the longest epsilon match in align from possEndsLeft and possEndsRight and sets the view positions of
// align to start and end position of the longest epsilon match
template<typename TLength, typename TSize, typename TEps>
Pair<typename Iterator<String<ExtensionEndPosition<TLength> > const>::Type>
longestEpsMatch(String<ExtensionEndPosition<TLength> > const & possEndsLeft,
				String<ExtensionEndPosition<TLength> > const & possEndsRight,
				TLength const alignLen,
				TLength const alignErr,
				TSize const matchMinLength,
				TEps const epsilon) {
	typedef ExtensionEndPosition<TLength>               TEnd;
	typedef typename Iterator<String<TEnd> const>::Type TIterator;

    // Identify longest eps match by iterating over combinations of left and right positions
    TIterator rightIt = end(possEndsRight) - 1;
    TIterator leftIt = end(possEndsLeft) - 1;
	TIterator right = begin(possEndsRight);
    TIterator left = begin(possEndsLeft);

	/*for (int i = 0; i < length(possEndsRight); ++i) {
		std::cout << possEndsRight[i].length << "  " << possEndsRight[i].coord.i1 << "," << possEndsRight[i].coord.i2 << std::endl;
	}
	for (int i = 0; i < length(possEndsLeft); ++i) {
		std::cout << possEndsLeft[i].length << "  " << possEndsLeft[i].coord.i1 << "," << possEndsLeft[i].coord.i2 << std::endl;
	}*/

	TSize leftErr = length(possEndsLeft) - 1;

    TSize minLength = matchMinLength;
	bool found = false;
    // DELTA is used below against floating point rounding errors.
	double const DELTA = 0.000001;

	while (leftIt >= begin(possEndsLeft)) {
		TSize totalLen = (*leftIt).length + alignLen + (*rightIt).length;
		if (totalLen < minLength) break;
		TSize totalErr = leftErr + alignErr + length(possEndsRight) - 1;
		while (rightIt >= begin(possEndsRight)) {
			totalLen = (*leftIt).length + alignLen + (*rightIt).length;
			if (totalLen < minLength) break;
			if ((TEps)totalErr/(TEps)totalLen < epsilon + DELTA) {
				right = rightIt;
				left = leftIt;
				//std::cout << totalLen << std::endl;
				minLength = totalLen;
				found = true;
				break;
			}
			--rightIt;
			--totalErr;
		}
		rightIt = end(possEndsRight) - 1;
		--leftIt;
		--leftErr;
	}

	if (found)
        return Pair<TIterator>(left, right);
	else
        return Pair<TIterator>(0,0);
}

///////////////////////////////////////////////////////////////////////////////
// Computes the banded alignment matrix and additionally a string with the best
//   alignment end point for each alignment length.
template <typename TTrace, typename TEnd, typename TStringSet, typename TScore, typename TDiagonal>
inline void
_align_banded_nw_best_ends(TTrace& trace,
						   String<TEnd> & bestEnds,
						   TStringSet const & str,
						   TScore const & sc,
						   TDiagonal const diagL,
						   TDiagonal const diagU)
{
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TTrace>::Type TSize;

	SEQAN_ASSERT_GEQ(diagU, diagL);

	// Initialization
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
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
	TRow mat, len;
	resize(mat, diagonalWidth);
	resize(len, diagonalWidth);
	resize(trace, height * diagonalWidth);

	// Classical DP with affine gap costs
	typedef typename Iterator<TRow, Standard>::Type TRowIter;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TSize actualCol = 0;
	TSize actualRow = 0;
	TScoreValue verti_val = 0;
	TScoreValue hori_val = 0;
	TScoreValue hori_len = len1+len2+1;
	TSize errors;

	for(TSize row = 0; row < height; ++row) {
		actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if ((TDiagonal)actualRow >= (TDiagonal)len1 - diagU) --hi_diag;
		TTraceIter traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
		TRowIter matIt = begin(mat, Standard()) + lo_diag;
		TRowIter lenIt = begin(len, Standard()) + lo_diag;
		hori_val = std::numeric_limits<TScoreValue>::min();
		hori_len = len1+len2+1;
		for(TSize col = lo_diag; col<hi_diag; ++col, ++matIt, ++traceIt, ++lenIt) {
			actualCol = col + diagL + actualRow;
			if (actualCol >= len1) break;

			if ((actualRow != 0) && (actualCol != 0)) {
				// Get the new maximum for mat
				*matIt += score(const_cast<TScore&>(sc), sequenceEntryForScore(const_cast<TScore&>(sc), str1, ((int) actualCol - 1)),
				                sequenceEntryForScore(const_cast<TScore&>(sc), str2, ((int) actualRow - 1)));
				*traceIt = Diagonal;
				++(*lenIt);
				if ((verti_val = (col < diagonalWidth - 1) ? *(matIt+1) +
				    scoreGapExtendVertical(sc,sequenceEntryForScore(sc, str1, ((int) actualCol - 1)),
				                           sequenceEntryForScore(sc, str2, ((int) actualRow - 1))) : std::numeric_limits<TScoreValue>::min()) > *matIt)
				{
					*matIt = verti_val;
					*traceIt = Vertical;
					*lenIt = *(lenIt+1) + 1;
				}						
				if ((hori_val = (col > 0) ? hori_val +
				    scoreGapExtendHorizontal(sc, sequenceEntryForScore(sc, str1, ((int) actualCol - 1)),
				                             sequenceEntryForScore(sc, str2, ((int) actualRow - 1))) : std::numeric_limits<TScoreValue>::min()) > *matIt)
				{
					*matIt = hori_val;
					*traceIt = Horizontal;
					*lenIt = hori_len + 1;
				}
				hori_val = *matIt;
				hori_len = *lenIt;
			} else {			
				// Usual initialization for first row and column
				if (actualRow == 0) {
					*matIt = actualCol * scoreGapExtendHorizontal(sc, sequenceEntryForScore(sc, str1, ((int) actualCol - 1)),
					                                              sequenceEntryForScore(sc, str2, -1));
					*lenIt = actualCol;
				}
				else {
					*matIt = actualRow * scoreGapExtendVertical(sc, sequenceEntryForScore(sc, str1, -1),
					                                            sequenceEntryForScore(sc, str2, ((int) actualRow - 1)));
					*lenIt = actualRow;
					hori_val = *matIt;
					hori_len = actualRow;
				}
			}
			errors = (*matIt - (*lenIt * scoreMatch(const_cast<TScore&>(sc)))) /
						(scoreGap(const_cast<TScore&>(sc)) - scoreMatch(const_cast<TScore&>(sc)));
			SEQAN_ASSERT_LEQ(errors, length(bestEnds));
			if (errors == length(bestEnds)) {
					appendValue(bestEnds, TEnd(*lenIt, row, col));
			} else if (*lenIt > static_cast<TScoreValue>(value(bestEnds, errors).length))
				value(bestEnds, errors) = TEnd(*lenIt, row, col);
			//std::cerr << row << ',' << col << ':' << *matIt << std::endl;
		}
	}
	TSize newLength = length(bestEnds) - 1;
	while (newLength > 0 && bestEnds[newLength].length <= bestEnds[newLength-1].length) {
		--newLength;
	}
	resize(bestEnds, newLength + 1);
}

///////////////////////////////////////////////////////////////////////////////
// Reverses the infixes of the left extension in place in hosts of infH and infV.
template<typename TSequenceA, typename TSequenceB, typename TSeed>
void
_reverseLeftExtension(Segment<TSequenceA, InfixSegment> const & infH,
					  Segment<TSequenceB, InfixSegment> const & infV,
					  TSeed & seed,
					  TSeed & seedOld) {
	Segment<TSequenceA, InfixSegment> infixH(host(infH), beginPositionH(seed), beginPositionH(seedOld));
	Segment<TSequenceB, InfixSegment> infixV(host(infV), beginPositionV(seed), beginPositionV(seedOld));
	reverse(infixH);
	reverse(infixV);
}

///////////////////////////////////////////////////////////////////////////////
// Computes the banded alignment matrix for the left extension and 
//   returns a string with possible start positions of an eps-match.
// Caution: The infixes of the left extension is reversed in place in hosts of infH and infV!
template<typename TMatrix, typename TPossEnd, typename TSequence, typename TSeed, typename TScore>
void
_fillMatrixBestEndsLeft(TMatrix & matrixLeft,
						String<TPossEnd> & possibleEndsLeft,
						Segment<TSequence, InfixSegment> const & infH,
						Segment<TSequence, InfixSegment> const & infV,
						TSeed & seed,
						TSeed & seedOld,
						TScore const & scoreMatrix) {

	typedef Segment<TSequence, InfixSegment> TInfix;

	TInfix infixH(host(infH), beginPositionH(seed), beginPositionH(seedOld));
	TInfix infixV(host(infV), beginPositionV(seed), beginPositionV(seedOld));

	reverse(infixH);
	reverse(infixV);

	StringSet<TInfix> str;
	appendValue(str, infixH);
	appendValue(str, infixV);

	// _align_banded_nw_best_ends(matrixLeft, possibleEndsLeft, str, scoreMatrix,
	// 						   upperDiagonal(seedOld) - upperDiagonal(seed),
	// 						   upperDiagonal(seedOld) - lowerDiagonal(seed));

    // Compute diagonals for updated seeds module with infixH/first alignment row being in the horizontal direction.
    typedef typename Diagonal<TSeed>::Type TDiagonal;
    TDiagonal diagLower = lowerDiagonal(seedOld) - upperDiagonal(seed);
    TDiagonal diagUpper = lowerDiagonal(seedOld) - lowerDiagonal(seed);

    // std::cerr << "FILL MATRIX LEFT SEQS\n"
    //           << "0: " << infixH << "\n"
    //           << "1: " << infixV << "\n";

    // _align_banded_nw_best_ends(matrixLeft, possibleEndsLeft, str, scoreMatrix,
    //                            diagBegin - upperDiagonal(seed),
    //                            diagBegin - lowerDiagonal(seed));
                               // upperDiagonal(seedOld) - upperDiagonal(seed),
                               // upperDiagonal(seedOld) - lowerDiagonal(seed));

    // Use legacy adapted NW computation with infixH/first alignment row being in the vertical direction.
    // // TODO(holtgrew): When switching to DP from new alignment module, make sure to mirror diagonals.
    _align_banded_nw_best_ends(matrixLeft, possibleEndsLeft, str, scoreMatrix, -diagUpper, -diagLower);
}

///////////////////////////////////////////////////////////////////////////////
// Computes the banded alignment matrix for the right extension and 
//   returns a string with possible end positions of an eps-match.
template<typename TMatrix, typename TPossEnd, typename TSequence, typename TSeed, typename TScore>
void
_fillMatrixBestEndsRight(TMatrix & matrixRight,
						 String<TPossEnd> & possibleEndsRight,
						 Segment<TSequence, InfixSegment> const & infH,
						 Segment<TSequence, InfixSegment> const & infV,
						 TSeed & seed,
						 TSeed & seedOld,
						 TScore const & scoreMatrix) {
	typedef Segment<TSequence, InfixSegment> TInfix;

	TInfix infixH(host(infH), endPositionH(seedOld), endPositionH(seed));
	TInfix infixV(host(infV), endPositionV(seedOld), endPositionV(seed));

	StringSet<TInfix> str;
	appendValue(str, infixH);
	appendValue(str, infixV);

    // std::cerr << "FILL MATRIX RIGHT SEQS\n"
    //           << "0: " << infixH << "\n"
    //           << "1: " << infixV << "\n";

	// _align_banded_nw_best_ends(matrixRight, possibleEndsRight, str, scoreMatrix, 
	// 						   lowerDiagonal(seedOld) - upperDiagonal(seed),
	// 						   lowerDiagonal(seedOld) - lowerDiagonal(seed));

    // Compute diagonals for updated seeds module with infixH/first alignment row being in the horizontal direction.
    typedef typename Diagonal<TSeed>::Type TDiagonal;
    TDiagonal diagLower = upperDiagonal(seedOld) - upperDiagonal(seed);
    TDiagonal diagUpper = upperDiagonal(seedOld) - lowerDiagonal(seed);

    // _align_banded_nw_best_ends(matrixRight, possibleEndsRight, str, scoreMatrix,
    //                            diagEnd - upperDiagonal(seed),
    //                            diagEnd - lowerDiagonal(seed));
    //                            lowerDiagonal(seedOld) - upperDiagonal(seed),
    //                            lowerDiagonal(seedOld) - lowerDiagonal(seed));

    // Use legacy adapted NW computation with infixH/first alignment row being in the vertical direction.
    // TODO(holtgrew): When switching to DP from new alignment module, make sure to mirror diagonals.
    _align_banded_nw_best_ends(matrixRight, possibleEndsRight, str, scoreMatrix, -diagUpper, -diagLower);
}

///////////////////////////////////////////////////////////////////////////////
// Traceback from an arbitrary point (coordinate) in the banded alignment trace matrix (trace).
template <typename TAlign, typename TStringSet, typename TTrace, typename TCoord, typename TDiagonal>
inline void
_alignBandedNeedlemanWunschTrace(TAlign & align,
								 TStringSet const & str,
								 TTrace const & trace,
								 TCoord const & coordinate,
								 TDiagonal const diagL,
								 TDiagonal const diagU)
{
	//typedef typename Value<TStringSet>::Type TString;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TTrace>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;

	// Traceback values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;

	// Initialization
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize lo_row = (diagU <= 0) ? -1 * diagU : 0;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);

	// Start the trace from the cell with the max value
	TSize row = coordinate.i1;
	TSize col = coordinate.i2;

	// Handle the skipped sequence parts
	TSize actualRow = row + lo_row;
	TSize actualCol = col + diagL + actualRow;

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
						_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
						--col; seqLen = 1;
					} else if (newTv == Vertical) {
						_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
						--row; ++col; seqLen = 1;
					} else {
						--row; ++seqLen;
					}
				} else {
					if (tv == Horizontal) { 
						if (newTv == Diagonal) {
							_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
							--row; seqLen = 1;
						} else if (newTv == Vertical) {
							_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
							--row; ++col; seqLen = 1;
						} else {
							--col; ++seqLen;
						}
					} else { 
						if (newTv == Diagonal) {
							_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
							--row; seqLen = 1;
						} else if (newTv == Horizontal) {
							_alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
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
		if (seqLen) _alignTracePrint(align, str[0], str[1], id1, actualCol, id2, actualRow, seqLen, tv);
	}

	// Handle the remaining sequence
	if (actualCol != 0) _alignTracePrint(align, str[0], str[1], (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol,  Horizontal);
	else if (actualRow != 0) _alignTracePrint(align, str[0], str[1], (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow,  Vertical);

}

///////////////////////////////////////////////////////////////////////////////
// Conducts the traceback on the extension to the left from best start position
//   and writes the result into align.
template<typename TMatrix, typename TCoord, typename TSequence, typename TSeed, typename TPos, typename TAlign>
void
_tracebackLeft(TMatrix const & matrixLeft,
			   TCoord const & coordinate,
			   Segment<TSequence, InfixSegment> const & infH,
			   Segment<TSequence, InfixSegment> const & infV,
			   TSeed & seed,
			   TSeed & seedOld,
			   TPos const endLeftH,
			   TPos const endLeftV,
			   TAlign & align) {
	typedef Segment<TSequence, InfixSegment>			TInfix;

	StringSet<TInfix> str;
	TInfix infixH(host(infH), beginPositionH(seed), beginPositionH(seedOld));
	TInfix infixV(host(infV), beginPositionV(seed), beginPositionV(seedOld));
	appendValue(str, infixH);
	appendValue(str, infixV);

    typedef typename Diagonal<TSeed>::Type TDiagonal;
    TDiagonal diagLower = lowerDiagonal(seedOld) - upperDiagonal(seed);
    TDiagonal diagUpper = lowerDiagonal(seedOld) - lowerDiagonal(seed);

	AlignTraceback<TPos> traceBack;
	_alignBandedNeedlemanWunschTrace(traceBack, str, matrixLeft, coordinate,
                                     -diagUpper, -diagLower);
                                     // upperDiagonal(seedOld) - upperDiagonal(seed), upperDiagonal(seedOld) - lowerDiagonal(seed));
  //std::cerr << "TRACEBACK\n";
	//for (unsigned i = 0; i < length(traceBack.tvs); ++i)
  //  std::cerr << (int)traceBack.tvs[i] << "\t" << traceBack.sizes[i] << "\n";
  //std::cerr << "---------\n";

	reverse(traceBack.sizes);
	reverse(traceBack.tvs);

	Align<TInfix> infixAlign;
	resize(rows(infixAlign), 2);
	assignSource(row(infixAlign, 0), infix(str[0], length(str[0]) - endLeftH, length(str[0])));
	assignSource(row(infixAlign, 1), infix(str[1], length(str[1]) - endLeftV, length(str[1])));

    // std::cerr << "\nLEFT SEQS\n" << row(infixAlign, 0) << "\n" << row(infixAlign, 1) << "\n";
	_pumpTraceToGaps(row(infixAlign, 0), row(infixAlign, 1), traceBack);
    // std::cerr << "INFIX ALIGN AFTER LEFT TRACEBACK\n\n" << infixAlign << "\n";
    // std::cerr << "ALIGN BEFORE INTEGRATION WITH INFIX ALIGN\n\n" << align << "\n";
	integrateAlign(align, infixAlign);
    // std::cerr << "ALIGN AFTER INTEGRATION WITH INFIX ALIGN\n\n" << align << "\n";
}


///////////////////////////////////////////////////////////////////////////////
// Conducts the traceback on the extension to the right from best end position
//   and writes the result into align.
template<typename TMatrix, typename TCoord, typename TSequence, typename TSeed, typename TPos, typename TAlign>
void
_tracebackRight(TMatrix const & matrixRight,
			   TCoord const & coordinate,
			   Segment<TSequence, InfixSegment> const & infH,
			   Segment<TSequence, InfixSegment> const & infV,
			   TSeed & seed,
			   TSeed & seedOld,
			   TPos const endRightH,
			   TPos const endRightV,
			   TAlign & align) {
	typedef Segment<TSequence, InfixSegment>			TInfix;

	TInfix infixH(host(infH), endPositionH(seedOld), endPositionH(seed));
	TInfix infixV(host(infV), endPositionV(seedOld), endPositionV(seed));

	StringSet<TInfix> str;
	appendValue(str, infixH);
	appendValue(str, infixV);

    typedef typename Diagonal<TSeed>::Type TDiagonal;
    TDiagonal diagLower = upperDiagonal(seedOld) - upperDiagonal(seed);
    TDiagonal diagUpper = upperDiagonal(seedOld) - lowerDiagonal(seed);

	AlignTraceback<TPos> traceBack;
	_alignBandedNeedlemanWunschTrace(traceBack, str, matrixRight, coordinate,
                                     -diagUpper, -diagLower);
				// lowerDiagonal(seedOld) - upperDiagonal(seed), lowerDiagonal(seedOld) - lowerDiagonal(seed));
  //std::cerr << "TRACEBACK\n";
	//for (unsigned i = 0; i < length(traceBack.tvs); ++i)
  //  std::cerr << (int)traceBack.tvs[i] << "\t" << traceBack.sizes[i] << "\n";
  //std::cerr << "---------\n";

	Align<TInfix> infixAlign;
	resize(rows(infixAlign), 2);
	assignSource(row(infixAlign, 0), infix(str[0], 0, endRightH));
	assignSource(row(infixAlign, 1), infix(str[1], 0, endRightV));

    // std::cerr << "\nRIGHT SEQS\n" << row(infixAlign, 0) << "\n" << row(infixAlign, 1) << "\n";
	_pumpTraceToGaps(row(infixAlign, 0), row(infixAlign, 1), traceBack);
    // std::cerr << "INFIX ALIGN AFTER RIGHT TRACEBACK\n\n" << infixAlign << "\n";
    // std::cerr << "ALIGN BEFORE INTEGRATION WITH INFIX ALIGN\n\n" << align << "\n";
	integrateAlign(align, infixAlign);
    // std::cerr << "ALIGN AFTER INTEGRATION WITH INFIX ALIGN\n\n" << align << "\n";
}

///////////////////////////////////////////////////////////////////////////////
// Computes the banded alignment matrix and fills a string with possible start 
//   and end positions of an eps-match. Determines the optimal start and end
//   position for the longest eps-match and writes the trace into align.
template<typename TInfix, typename TSeed, typename TPos, typename TDir, typename TScore,
		 typename TSize, typename TEps, typename TAlign>
bool
_bestExtension(TInfix const & infH,
			   TInfix const & infV,
			   TSeed & seed,
			   TSeed & seedOld,
			   TPos const alignLen,
			   TPos const alignErr,
			   TScore const & scoreMatrix,
			   TDir const direction,
			   TSize const minLength,
			   TEps const eps,
			   TAlign & align)
{
	typedef String<TraceBack>							    TAlignmentMatrix;
	typedef ExtensionEndPosition<TPos>                      TEndInfo;
	typedef typename Iterator<String<TEndInfo> const>::Type TEndIterator;

	// variables for banded alignment and possible ends of match
	TAlignmentMatrix matrixRight, matrixLeft;
	String<TEndInfo> possibleEndsLeft, possibleEndsRight;

	// fill banded matrix and gaps string for ...
	if (direction == EXTEND_BOTH || direction == EXTEND_LEFT) { // ... extension to the left
		_fillMatrixBestEndsLeft(matrixLeft, possibleEndsLeft, infH, infV, seed, seedOld, scoreMatrix);
		// Caution: left extension infix is now reversed in host(infH and infV) !!!
        SEQAN_ASSERT_NOT(empty(possibleEndsLeft));
	} else appendValue(possibleEndsLeft, TEndInfo());
	if (direction == EXTEND_BOTH || direction == EXTEND_RIGHT) { // ... extension to the right
		_fillMatrixBestEndsRight(matrixRight, possibleEndsRight, infH, infV, seed, seedOld, scoreMatrix);
        SEQAN_ASSERT_NOT(empty(possibleEndsRight));
	} else appendValue(possibleEndsRight, TEndInfo());

	// longest eps match on poss ends string
	Pair<TEndIterator> endPair = longestEpsMatch(possibleEndsLeft, possibleEndsRight, alignLen, alignErr, minLength, eps);

	if (endPair == Pair<TEndIterator>(0, 0)) { // no eps-match found
		if (direction != 1) 
			_reverseLeftExtension(infH, infV, seed, seedOld); // back to original orientation
		return false;
	}

	// determine end positions of maximal eps-match in ...
	TPos endLeftH = 0, endLeftV = 0;
	TPos endRightH = 0, endRightV = 0;
	if((*endPair.i1).length != 0) { // ... extension to the left
		endLeftV = (*endPair.i1).coord.i1;
		// correction for banded coordinates to unbanded:
		if (upperDiagonal(seed) - lowerDiagonal(seedOld)  <= 0)
			endLeftV -= (TPos)(upperDiagonal(seed) - lowerDiagonal(seedOld));
		endLeftH = (TPos)((*endPair.i1).coord.i2 + endLeftV + lowerDiagonal(seed) - lowerDiagonal(seedOld));
	}
	if((*endPair.i2).length != 0) { // ... extension to the right
		endRightV = (*endPair.i2).coord.i1;
		// correction for banded coordinates to unbanded:
		if (upperDiagonal(seed) - upperDiagonal(seedOld) <= 0)
			endRightV -= (TPos)(upperDiagonal(seed) - upperDiagonal(seedOld));
		endRightH = (TPos)((*endPair.i2).coord.i2 + endRightV + lowerDiagonal(seed) - upperDiagonal(seedOld));
	}

	// set begin and end positions of align
	setBeginPosition(row(align, 0), beginPositionH(seedOld) - endLeftH);
	setBeginPosition(row(align, 1), beginPositionV(seedOld) - endLeftV);
	setEndPosition(row(align, 0), endPositionH(seedOld) + endRightH);
	setEndPosition(row(align, 1), endPositionV(seedOld) + endRightV);
	// setClippedBeginPosition(row(align, 0), beginPositionH(seedOld) - endLeftH);
	// setClippedBeginPosition(row(align, 1), beginPositionV(seedOld) - endLeftV);
	// setBeginPosition(row(align, 0), 0);
	// setBeginPosition(row(align, 1), 0);
	// setClippedEndPosition(row(align, 0), endPositionH(seedOld) + endRightH);
	// setClippedEndPosition(row(align, 1), endPositionV(seedOld) + endRightV);

	// traceback through matrix from begin/end pos on ...
	if((*endPair.i1).length != 0) { // ... extension to the left
		_tracebackLeft(matrixLeft, (*endPair.i1).coord, infH, infV, seed, seedOld, endLeftH, endLeftV, align);
	}
	if((*endPair.i2).length != 0) { // ... extension to the right
		_tracebackRight(matrixRight, (*endPair.i2).coord, infH, infV, seed, seedOld, endRightH, endRightV, align);
	}
    SEQAN_ASSERT_EQ(length(row(align, 0)), length(row(align, 1)));

	if (direction == EXTEND_BOTH || direction == EXTEND_LEFT)
		_reverseLeftExtension(infH, infV, seed, seedOld); // back to original orientation

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSpec1, typename TSpec2>
void
integrateAlign(Align<TSource, TSpec1> & align,
			   Align<Segment<Segment<TSource, InfixSegment>, InfixSegment>, TSpec2> const & infixAlign) {
	typedef typename Size<TSource>::Type TSize;
	typedef typename Position<typename Row<Align<TSource, TSpec1> >::Type>::Type TPos;

	String<TPos> viewPos;
	TPos pos;
	for (TSize i = 0; i < length(rows(infixAlign)); ++i) {
		pos = beginPosition(source(row(infixAlign, i))) + beginPosition(row(infixAlign, i));
		pos += beginPosition(host(source(row(infixAlign, i))));
		appendValue(viewPos, toViewPosition(row(align, i), pos));
	}

    // std::cerr << "HAHA infixAlign == \n" << row(infixAlign, 0) << "\n" << row(infixAlign, 1) << "\n";
	integrateAlign(align, infixAlign, viewPos);
    // std::cerr << "HAHA infixAlign == \n" << row(infixAlign, 0) << "\n" << row(infixAlign, 1) << "\n";
}

///////////////////////////////////////////////////////////////////////////////
// Conducts best X-drop extension and calls _bestExtension.
//  After the call align contains the longest eps-Match that spans the eps-core (localAlign).
template<typename TScoreValue, typename TScore, typename TSequence, typename TSize, typename TEps, typename TAlign>
bool
_extendAndExtract(Align<Segment<Segment<TSequence, InfixSegment>, InfixSegment> > const & localAlign,
				  TScoreValue scoreDropOff,
				  TScore const & scoreMatrix,
				  Segment<typename Infix<TSequence>::Type, InfixSegment> const & infH,
				  Segment<typename Infix<TSequence>::Type, InfixSegment>  const & infV,
				  ExtensionDirection direction,
				  TSize minLength,
				  TEps eps,
				  TAlign & align) {
    typedef typename Position<TSequence>::Type TPos;
    typedef Seed<Simple> TSeed;

    // std::cerr << "LOCAL ALIGN\n" << row(localAlign, 0) << "\n" << row(localAlign, 1) << "\n";
    // std::cerr << "ALIGN\n" << row(align, 0) << "\n" << row(align, 1) << "\n";
    integrateAlign(align, localAlign);
    //std::cerr << __LINE__ << "\tLOCAL ALIGN\n" << row(localAlign, 0) << "\n" << row(localAlign, 1) << "\n";
    //std::cerr << __LINE__ << "\tALIGN\n" << row(align, 0) << "\n" << row(align, 1) << "\n";

	// Get begin and end position of local alignment (seed) as source positions
	// in underlying sequences.
	TPos seedBeginH = beginPosition(row(localAlign, 0)) + beginPosition(infH);
	TPos seedBeginV = beginPosition(row(localAlign, 1)) + beginPosition(infV);
	TPos seedEndH = endPosition(row(localAlign, 0)) + beginPosition(infH);
	TPos seedEndV = endPosition(row(localAlign, 1)) + beginPosition(infV);

	if (direction == EXTEND_NONE) {
		// set begin and end positions of align
		setBeginPosition(row(align, 0), seedBeginH);
		setBeginPosition(row(align, 1), seedBeginV);
		setEndPosition(row(align, 0), seedEndH);
		setEndPosition(row(align, 1), seedEndV);

		if ((TSize)length(row(align, 0)) < minLength)
			return false;

		longestEpsMatch(align, minLength, eps);
	} else {
		// gapped X-drop extension of local alignment (seed)
		TSeed seed(seedBeginH, seedBeginV, seedEndH, seedEndV);
		TSeed seedOld(seed);
		extendSeed(seed, host(infH), host(infV), direction, scoreMatrix, scoreDropOff, GappedXDrop());

		if (static_cast<int64_t>(seedSize(seed)) < minLength - (int)floor(minLength*eps))
			return false;

		// determine length and number of error columns of local alignment (seed)
		TPos alignLen = _max(length(row(localAlign, 0)), length(row(localAlign, 1)));
		TPos alignErr = 0;
		for (TPos i = 0; i < alignLen; ++i) {
			if (!isMatch(localAlign, i)) ++alignErr;
		}

		// convert seeds from positions in host(seq) to positions in host(host(seq))
		setBeginPositionH(seedOld, beginPositionH(seedOld) + beginPosition(host(infH)));
		setEndPositionH(seedOld, endPositionH(seedOld) + beginPosition(host(infH)));
		setBeginPositionV(seedOld, beginPositionV(seedOld) + beginPosition(host(infV)));
		setEndPositionV(seedOld, endPositionV(seedOld) + beginPosition(host(infV)));
		setBeginPositionH(seed, beginPositionH(seed) + beginPosition(host(infH)));
		setEndPositionH(seed, endPositionH(seed) + beginPosition(host(infH)));
		setBeginPositionV(seed, beginPositionV(seed) + beginPosition(host(infV)));
		setEndPositionV(seed, endPositionV(seed) + beginPosition(host(infV)));

		// determine best extension lengths and write the trace into align
		typename Infix<TSequence>::Type infixH = infix(host(infH), beginPosition(infH), endPosition(infH));
		typename Infix<TSequence>::Type infixV = infix(host(infV), beginPosition(infV), endPosition(infV));
		if (!_bestExtension(infixH, infixV, seed, seedOld, alignLen, alignErr, scoreMatrix, direction, minLength, eps, align))
			return false;
		SEQAN_ASSERT_EQ(length(row(align, 0)), length(row(align, 1)));
	}
    SEQAN_ASSERT_EQ(length(row(align, 0)), length(row(align, 1)));
    //std::cerr << "extracted alignment\n-------------\n" << align << "----------------\n";
	return true;
}

#endif
