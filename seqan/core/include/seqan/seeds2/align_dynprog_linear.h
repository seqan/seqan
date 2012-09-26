// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Classic linear programming for sequence alignment.  Needleman Wunsch
// implementation, i.e. for linear gap costs only.
//
// The alignment code could use some optimization.  However, we cannot
// use the same optimization as in the graph alignment since we want to
// compute globally optimal trace through multiple connected alignment
// matrices.
// ==========================================================================

// TODO(holtgrew): Maybe adjust rausch's code so things can be copied in and use it instead since it is heavily tuned?

#ifndef SEQAN_SEEDS_ALIGN_DYNPROG_LINEAR_H_
#define SEQAN_SEEDS_ALIGN_DYNPROG_LINEAR_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TScoreValue, typename TSequence>
inline void
_alignResizeMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    setLength(matrix, 0, length(sequence0) + 1);
    setLength(matrix, 1, length(sequence1) + 1);
    resize(matrix, -420);
    // resize(matrix);
}


template <typename TScoreValue, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_alignInitGutter(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const & scoringScheme, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");

    // Init left gutter with zeroes if begin gaps are in dimension 0
    // free or with gap scores otherwise.
    if (BEGIN0_FREE) {
        TIterator it = begin(matrix);
        for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i) {
            *it = 0;
            goNext(it, 0);
        }
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        TIterator it = begin(matrix);
        for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i) {
            *it = x;
            x += gapScore;
            goNext(it, 0);
        }
    }

    // Init top gutter with zeroes if begin gaps are in dimension 1
    // free or with gap scores otherwise.
    if (BEGIN1_FREE) {
        TIterator it = begin(matrix);
        for (TPosition i = 0, iend = length(matrix, 1); i < iend; ++i) {
            *it = 0;
            goNext(it, 1);
        }
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        TIterator it = begin(matrix);
        for (TPosition i = 0, iend = length(matrix, 1); i < iend; ++i) {
            *it = x;
            x += gapScore;
            goNext(it, 1);
        }
    }
}

template <typename TScoreValue, typename TDiagonal, typename TOverlap>
inline void
_alignInitGutterFromBanded(Matrix<TScoreValue, 2> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Matrix<TScoreValue, 2> /*const*/ & otherMatrix, TOverlap overlap0, TOverlap overlap1, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;
    (void) scoringScheme; // Only used in assertions.

    // TODO(holtgrew): Really unnecessary? Remove along with all other unused parameters in all align_*.h files.
    (void) lowerDiagonal;
    (void) upperDiagonal;

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");

    // Copy over data for left gutter.
    // std::cout << "overlap0 = " << overlap0 << ", overlap1 = " << overlap1 << std::endl;
    TIterator it = begin(matrix);
    TIterator otherIt(otherMatrix);
    goTo(otherIt, length(otherMatrix, 0) - (overlap0 + 1), overlap0);
    TIterator srcIt = otherIt;
    for (TOverlap i = 0; i < overlap0 + 1; ++i) {
        *it = *srcIt;
        goNext(it, 0);
        goNext(srcIt, 0);
        goPrevious(srcIt, 1);
    }
    // Init the rest of the left gutter with infima.
    for (TPosition i = overlap0 + 1, iend = length(matrix, 0); i < iend; ++i, goNext(it, 0))
        *it = MinValue<TScoreValue>::VALUE / 2;

    // Copy over data for top gutter.
    it = begin(matrix);
    srcIt = otherIt;
    for (TOverlap i = 0; i < overlap1 + 1; ++i) {
        *it = *srcIt;
        goNext(it, 1);
        goNext(srcIt, 1);
    }
    // Init the rest of the top gutter with infima.
    for (TPosition i = overlap1 + 1, iend = length(matrix, 1); i < iend; ++i, goNext(it, 1))
        *it = MinValue<TScoreValue>::VALUE / 2;
}


template <typename TScoreValue, typename TSequence>
inline void
_alignFillMatrix(Matrix<TScoreValue, 2> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ_MSG(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme),
                        "Only linear gap costs allowed for Needleman-Wunsch.");

    // Fill matrix with standard NW DP programming.

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Iterator<TSequence, Standard>::Type TSequenceIterator;

    // We need three iterators in the alignment matrix to fill it.
    // itTop points to the cell in the top row of the current column.
    // itLeft points to the column to the top left of the current
    // cell.  itAbove points to the cell above the current cell.  We
    // can use itAbove for value assignment.
    TMatrixIterator itTop = begin(matrix);
    TMatrixIterator itAbove;
    TMatrixIterator itLeft;

    // Compute score values so they are available locally and it maybe
    // is easier to the compiler to save some indirect adressing.
    TScoreValue matchScore = scoreMatch(scoringScheme);
    TScoreValue mismatchScore = scoreMismatch(scoringScheme);
    TScoreValue gapScore = scoreGap(scoringScheme);
    
    // Perform the Needleman-Wunsch dynamic programming.
    // TODO(holtgrew): camelCase for itXend.
    for (TSequenceIterator it1 = begin(sequence1), it1end = end(sequence1); it1 != it1end; ++it1) {
        itLeft = itTop;
        goNext(itTop, 1);
        itAbove = itTop;
        for (TSequenceIterator it0 = begin(sequence0), it0end = end(sequence0); it0 != it0end; ++it0) {
            TScoreValue scoreMoveDiagonal = *itLeft + ((*it0 == *it1) ? matchScore : mismatchScore);
            goNext(itLeft, 0);
            TScoreValue scoreMoveRight = *itLeft + gapScore;
            TScoreValue scoreMoveDown = *itAbove + gapScore;
            goNext(itAbove, 0);
            *itAbove = _max(scoreMoveDiagonal, _max(scoreMoveRight, scoreMoveDown));
        }
    }

    // // TODO(holtgrew): Debug code, remove when working.
    // {
    //     for (int k = 0; k < 1; ++k) {
    //         std::cerr << ",-- *** filled unbanded alignment matrix " << k << std::endl;
    //         for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //             for (unsigned j = 0; j < length(matrix, 1); ++j) {
    //                 if (value(matrix, i, j, k) <= MinValue<int>::VALUE / 4)
    //                     std::cerr << "\tinf";
    //                 else
    //                     std::cerr << "\t" << value(matrix, i, j, k);
    //             }
    //             std::cerr << std::endl;
    //         }
    //         std::cerr << "`--" << std::endl;
    //     }
    // }
}


// Compute traceback in the given normal DP alignment matrix, starting
// at the lower right, shifted by the given overlap to the upper left.
// The matrix was filled with the Needleman-Wunschla lgorithm, end
// gaps are free as configured by the AlignConfig object.  Returns the
// best score.
template <typename TAlignmentIterator, typename TSequenceIterator, typename TPosition, typename TScoreValue, typename TScoringScheme, typename TOverlap, bool START0_FREE, bool START1_FREE, bool END0_FREE, bool END1_FREE>
TScoreValue
_alignTraceback(TAlignmentIterator & alignmentIt0, TAlignmentIterator & alignmentIt1, TSequenceIterator & sourceIt0, TSequenceIterator & sourceIt1, TPosition & finalPos0, TPosition & finalPos1, Matrix<TScoreValue, 2> /*const*/ & matrix, TScoringScheme const & scoringScheme, TOverlap lowerRightOverlap0, TOverlap lowerRightOverlap1, TOverlap upperLeftOverlap0, TOverlap upperLeftOverlap1, bool goToTopLeft, AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const &, NeedlemanWunsch const &)
{
    SEQAN_CHECKPOINT;

    // std::cout << "trace back unbanded" << std::endl;
    // // TODO(holtgrew): Debug output, remove when not needed any more.
    // {
    //     std::cout << ",-- filled unbanded alignment matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
    //     for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //         std::cout << "| ";
    //         for (unsigned j = 0; j < length(matrix, 1); ++j) {
    //             if (value(matrix, i, j) == MinValue<int>::VALUE / 2)
    //                 std::cout << "\tinf";
    //             else
    //                 std::cout << "\t" << value(matrix, i, j);
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "`--" << std::endl;
    // }
    // std::cout << "length(matrix, 0) = " << length(matrix, 0) << std::endl;
    // std::cout << "lowerRightOverlap0 = " << lowerRightOverlap0 << std::endl;
    // std::cout << "finalPos0 = " << finalPos0 << std::endl;
    // std::cout << "length(matrix, 1) = " << length(matrix, 1) << std::endl;
    // std::cout << "lowerRightOverlap1 = " << lowerRightOverlap1 << std::endl;
    // std::cout << "finalPos1 = " << finalPos1 << std::endl;
    // std::cout << "begin pos0 = " << (length(matrix, 0) - lowerRightOverlap0 + finalPos0) << std::endl;
    // std::cout << "begin pos1 = " << (length(matrix, 1) - lowerRightOverlap1 + finalPos1) << std::endl;

    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;

    // Initialization
    //
    // Precomputation of the score difference between mismatch and gap.
	TScoreValue scoreDifference = scoreMismatch(scoringScheme) - scoreGap(scoringScheme);
    TScoreValue matchScore = scoreMatch(scoringScheme);
    // Current position in the matrix.  Note that this is the position
    // in the matrix including the gutter.  When writing this out, we
    // want to have the position in the matrix, excluding the gutter,
    // i.e. everything is shifted to the upper left.
    TPosition pos0 = length(matrix, 0) - lowerRightOverlap0 + finalPos0 - 1;
    TPosition pos1 = length(matrix, 1) - lowerRightOverlap1 + finalPos1 - 1;
    TPosition origPos0 = pos0;
    (void) origPos0;  // In case we run without assertions.
    TPosition origPos1 = pos1;
    // Iterator to current entry in the matrix.
    TMatrixIterator matrixIt = begin(matrix);
    goTo(matrixIt, pos0, pos1);
    TMatrixIterator origMatrixIt = matrixIt;

    // Search for starting point of the trace.
    //
    // If end gaps are free in any direction, search for the position
    // with the highest score on the lower/right edge.  Prefer end
    // positions closer to the end and shifts in sequence 0 over
    // those in sequence 1.
    if (END0_FREE) {
        TMatrixIterator it = origMatrixIt;
        TPosition altPos1 = pos1;
        while (altPos1 > 1) {
            goPrevious(it, 1);
            altPos1 -= 1;
            if (*it > *origMatrixIt) {
                matrixIt = it;
                pos1 = altPos1;
            }
        }
    }
    if (END1_FREE) {
        TMatrixIterator it = origMatrixIt;
        TPosition altPos0 = pos0;
        while (altPos0 > 0) {
            goPrevious(it, 0);
            altPos0 -= 1;
            if (*it > *origMatrixIt) {
                matrixIt = it;
                pos1 = origPos1;
                pos0 = altPos0;
            }
        }
    }

    // Insert free end gaps if any.
    //
    // TODO(holtgrew): We could have 3 cases here but I guess it will not matter performance wise in the big picture.
    if (END0_FREE || END1_FREE) {
        SEQAN_ASSERT(pos0 == origPos0 || pos1 == origPos1);
        // std::cout << __FILE__ << ":" << __LINE__ << "-- Inserting " << (origPos0 - pos0) << " end gaps into sequence 1" << std::endl;
        insertGaps(alignmentIt1, origPos0 - pos0);
        // std::cout << __FILE__ << ":" << __LINE__ << "-- Inserting " << (origPos1 - pos1) << " end gaps into sequence 0" << std::endl;
        insertGaps(alignmentIt0, origPos1 - pos1);
    }

    // Now, perform the traceback.
    while (true) {
        // Break condition depends on whether this is the leading
        // rectangle, i.e. we want to go to the upper left corner.
        if (goToTopLeft) {
            if (pos0 == static_cast<TPosition>(0) || pos1 == static_cast<TPosition>(0))
                break;
        } else {
            // If we do not walk up to the upper left corner, we at
            // least have to move into the overlapping area.
            if ((pos0 <= static_cast<TPosition>(1) || pos1 <= static_cast<TPosition>(1)) || (pos0 < static_cast<TPosition>(upperLeftOverlap0) && pos1 < static_cast<TPosition>(upperLeftOverlap1)))
                break;
        }
        
        // Flags for "go horizontal" and "go vertical".
        bool gh = false;
        bool gv = false;

        // Determine whether to go vertical/horizontal.
        // std::cout << "Compare: " << *sourceIt0 << ", " << *sourceIt1 << std::endl;
        // std::cout << "  Score: " << *matrixIt << std::endl
        TMatrixIterator tmpIt = matrixIt;
        goPrevious(tmpIt, 0);
        goPrevious(tmpIt, 1);
        if ((*sourceIt0 == *sourceIt1) && (*tmpIt + matchScore == *matrixIt)) {
            gh = true;
            gv = true;
        } else {
			TMatrixIterator it = matrixIt;

			goPrevious(it, 0);
			TScoreValue v = *it;

			goPrevious(it, 1);
			TScoreValue d = *it;

			it = matrixIt;
			goPrevious(it, 1);
			TScoreValue h = *it;

			gv = (v >= h) || (d + scoreDifference >= h);
			gh = (h > v) || (d + scoreDifference >= v);
        }

        // if (gv && gh) {
        //     std::cout << "GO DIAGONAL (" << "*sourceIt0 == " << *sourceIt0 << ", *sourceIt1 == " << *sourceIt1 << ") *matrixIt == " << *matrixIt << ", pos0 == " << pos0 << ", pos1 == " << pos1 << std::endl;
        // } else if (gv) {
        //     std::cout << "GO VERTICAL (" << "*sourceIt0 == " << *sourceIt0 << ", *sourceIt1 == " << *sourceIt1 << ") *matrixIt == " << *matrixIt << ", pos0 == " << pos0 << ", pos1 == " << pos1 << std::endl;
        // } else if (gh) {
        //     std::cout << "GO HORIZONTAL (" << "*sourceIt0 == " << *sourceIt0 << ", *sourceIt1 == " << *sourceIt1 << ") *matrixIt == " << *matrixIt << ", pos0 == " << pos0 << ", pos1 == " << pos1 << std::endl;
        // }

        // Move iterators in source sequence, alignment rows, matrix
        // and possibly insert gaps.
        if (gv) {
            // std::cout << "  moving in dimension 0" << std::endl;
            goPrevious(sourceIt0);
            goPrevious(alignmentIt0);
            // std::cout << "  *alignmentIt0 == " << convert<char>(*alignmentIt0) << std::endl;
            goPrevious(matrixIt, 0);
            pos0 -= 1;
        } else {
            // std::cout << __FILE__ << ":" << __LINE__ << "-- gap in dimension 0" << std::endl;
            insertGap(alignmentIt0);
        }
        if (gh) {
            // std::cout << "  moving in dimension 1" << std::endl;
            goPrevious(sourceIt1);
            goPrevious(alignmentIt1);
            // std::cout << "  *alignmentIt1 == " << convert<char>(*alignmentIt1) << std::endl;
            goPrevious(matrixIt, 1);
            pos1 -= 1;
        } else {
            // std::cout << __FILE__ << ":" << __LINE__ << "-- gap in dimension 0" << std::endl;
            insertGap(alignmentIt1);
        }
    }
    // std::cout << "pos0 == " << pos0 << ", pos1 == " << pos1 << std::endl;

    // TODO(holtgrew): I guess we have to search for the maximum in the resulting row and column after all...
    
    // Go to the top left of the matrix if configured to do so.
    if (goToTopLeft) {
        if (pos0 > 0) {
            goPrevious(alignmentIt1);
            // std::cout << "Inserting " << pos0 << " gaps into alignment row 1" << std::endl;
            for (TPosition i = 0; i < pos0; ++i) {
                goPrevious(sourceIt0);
                goPrevious(alignmentIt0);
                // std::cout << "  *alignmentIt0 == " << convert<char>(*alignmentIt0) << std::endl;
                // std::cout << __FILE__ << ":" << __LINE__ << "-- gap in dimension 1" << std::endl;
                insertGap(alignmentIt1);
            }
            pos0 = 0;
        }
        if (pos1 > 0) {
            goPrevious(alignmentIt0);
            // std::cout << "Inserting " << pos1 << " gaps into alignment row 0" << std::endl;
            for (TPosition i = 0; i < pos1; ++i) {
                goPrevious(sourceIt1);
                goPrevious(alignmentIt1);
                // std::cout << "  *alignmentIt1 == " << convert<char>(*alignmentIt1) << std::endl;
                // std::cout << __FILE__ << ":" << __LINE__ << "-- gap in dimension 0" << std::endl;
                insertGap(alignmentIt0);
            }
            pos1 = 0;
        }
    }

    // Write out the final positions in the alignment matrix.
    // Position is in the current matrix, including gutter since we
    // can get there.
    finalPos0 = pos0;
    finalPos1 = pos1;
    // std::cout << "finalPos0 = " << finalPos0 << ", finalPos1 = " << finalPos1 << std::endl;

    return *origMatrixIt;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_LINEAR_H_
