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
// Classic Gotoh DP algorithm for affine gap costs.  The algorithm is
// split into matrix resizing, filling the gutter, filling the rest of
// the matrix and doing the traceback.  Otherwise, it is a text book
// implementation without any tricks.
// ==========================================================================
// We implement the three matrices as one 3-dimensional matrix.  The
// first two dimensions are for positions in sequences 0 and 1.  The
// third dimension differentiates between the three matrices.  They are
// -- in order -- M, I^a, I^b.
// ==========================================================================

// TODO(holtgrew): Maybe adjust rausch's code so things can be copied in and use it instead since it is heavily tuned?

#ifndef SEQAN_SEEDS_ALIGN_DYNPROG_AFFINE_H_
#define SEQAN_SEEDS_ALIGN_DYNPROG_AFFINE_H_

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
_alignResizeMatrix(Matrix<TScoreValue, 3> & matrix, TSequence const & sequence0, TSequence const & sequence1, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    setLength(matrix, 0, length(sequence0) + 1);
    setLength(matrix, 1, length(sequence1) + 1);
    setLength(matrix, 2, 3);
    resize(matrix);
}


template <typename TScoreValue, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_alignInitGutter(Matrix<TScoreValue, 3> & matrix, Score<TScoreValue, Simple> const scoringScheme, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ(length(matrix, 2), 3u);

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;

    // TODO(holtgrew): Support free begin gaps.
    SEQAN_ASSERT_NOT_MSG(BEGIN0_FREE, "Free begin gaps are not supported in seeds align module yet.");
    SEQAN_ASSERT_NOT_MSG(BEGIN1_FREE, "Free begin gaps are not supported in seeds align module yet.");

    // We do not take the real minimum here because of overflows.
    TScoreValue inf = MinValue<TScoreValue>::VALUE / 2;
    // Get shortcuts to gap related scores.
    TScoreValue gapOpenScore = scoreGapOpen(scoringScheme);
    TScoreValue gapExtendScore = scoreGapExtend(scoringScheme);

    // Get iterators to the top left corners of the matrices.
    TIterator itMBegin = begin(matrix);
    TIterator itIABegin = itMBegin;
    goNext(itIABegin, 2);
    TIterator itIBBegin = itIABegin;
    goNext(itIBBegin, 2);

    // Init top and left gutters.
    //
    // Left gutters...
    TIterator itM = itMBegin;
    *itM = 0;
    goNext(itM, 0);
    TIterator itIA = itIABegin;
    *itIA = inf;
    goNext(itIA, 0);
    TIterator itIB = itIBBegin;
    *itIB = inf;
    goNext(itIB, 0);
    for (TPosition pos = 1, posEnd = length(matrix, 0); pos < posEnd; ++pos) {
        *itM = gapOpenScore + (pos - 1) * gapExtendScore;
        *itIA = inf;
        *itIB = gapOpenScore + (pos - 1) * gapExtendScore;
        goNext(itM, 0);
        goNext(itIA, 0);
        goNext(itIB, 0);
    }
    // Top gutters...
    itM = itMBegin;
    itIA = itIABegin;
    itIB = itIBBegin;
    goNext(itM, 1);
    goNext(itIA, 1);
    goNext(itIB, 1);
    for (TPosition pos = 1, posEnd = length(matrix, 1); pos < posEnd; ++pos) {
        *itM = gapOpenScore + (pos - 1) * gapExtendScore;
        *itIA = gapOpenScore + (pos - 1) * gapExtendScore;
        *itIB = inf;
        goNext(itM, 1);
        goNext(itIA, 1);
        goNext(itIB, 1);
    }
}


template <typename TScoreValue, typename TDiagonal, typename TOverlap>
inline void
_alignInitGutterFromBanded(Matrix<TScoreValue, 3> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Matrix<TScoreValue, 3> /*const*/ & otherMatrix, TOverlap overlap0, TOverlap overlap1, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    // TODO(holtgrew): Really unnecessary? Remove along with all other unused parameters in all align_*.h files.
    (void) scoringScheme;
    (void) lowerDiagonal;
    (void) upperDiagonal;

    SEQAN_ASSERT_EQ(length(matrix, 2), 3u);

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;

    // Get iterators to the top left corners of the matrices.
    TIterator itMBegin = begin(matrix);
    TIterator itIABegin = itMBegin;
    goNext(itIABegin, 2);
    TIterator itIBBegin = itIABegin;
    goNext(itIBBegin, 2);
    TIterator otherItMBegin = begin(otherMatrix);
    goTo(otherItMBegin, length(otherMatrix, 0) - (overlap0 + 1), overlap0);
    TIterator otherItIABegin = otherItMBegin;
    goNext(otherItIABegin, 2);
    TIterator otherItIBBegin = otherItIABegin;
    goNext(otherItIBBegin, 2);

    // Copy over data for left gutter.
    {
        TIterator itM = itMBegin;
        TIterator otherItM = otherItMBegin;
        TIterator itIB = itIBBegin;
        TIterator otherItIB = otherItIBBegin;
        for (TOverlap i = 0; i < overlap0 + 1; ++i) {
            *itM = *otherItM;
            *itIB = *otherItIB;
            goNext(itM, 0);
            goNext(otherItM, 0);
            goPrevious(otherItM, 1);
            goNext(itIB, 0);
            goNext(otherItIB, 0);
            goPrevious(otherItIB, 1);
        }
        // Init the rest of the left gutter with infima.
        for (TPosition i = overlap0 + 1, iend = length(matrix, 0); i < iend; ++i, goNext(itM, 0), goNext(itIB, 0)) {
            *itM = MinValue<TScoreValue>::VALUE / 2;
            *itIB = MinValue<TScoreValue>::VALUE / 2;
        }
    }

    // Copy over data for top gutter.
    {
        TIterator itM = itMBegin;
        TIterator otherItM = otherItMBegin;
        TIterator itIA = itIABegin;
        TIterator otherItIA = otherItIABegin;
        for (TOverlap i = 0; i < overlap1 + 1; ++i) {
            *itM = *otherItM;
            *itIA = *otherItIA;
            goNext(itM, 1);
            goNext(otherItM, 1);
            goNext(itIA, 1);
            goNext(otherItIA, 1);
        }
        // Init the rest of the top gutter with infima.
        for (TPosition i = overlap1 + 1, iend = length(matrix, 1); i < iend; ++i, goNext(itM, 1), goNext(itIA, 1)) {
            *itM = MinValue<TScoreValue>::VALUE / 2;
            *itIA = MinValue<TScoreValue>::VALUE / 2;
        }
    }
}


template <typename TScoreValue, typename TSequence>
inline void
_alignFillMatrix(Matrix<TScoreValue, 3> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_EQ(length(matrix, 2), 3u);

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Position<TMatrix>::Type TPosition;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Iterator<TSequence, Standard>::Type TSequenceIterator;

    // We need three iterators in each alignment matrix to fill it.
    // it*Top point to the cell in the top row of the current column.
    // it*Left points to the column to the top left of the current
    // cell.  it*Above points to the cell above the current cell.  We
    // can use it*Above for value assignment.
    TMatrixIterator itMTop = begin(matrix);
    TMatrixIterator itMAbove;
    TMatrixIterator itMLeft;
    TMatrixIterator itIATop = itMTop;
    goNext(itIATop, 2);
    TMatrixIterator itIAAbove;
    TMatrixIterator itIALeft;
    TMatrixIterator itIBTop = itIATop;
    goNext(itIBTop, 2);
    TMatrixIterator itIBAbove;
    TMatrixIterator itIBLeft;

    // Compute score values so they are available locally and it maybe
    // is easier to the compiler to save some indirect adressing.
    TScoreValue matchScore = scoreMatch(scoringScheme);
    TScoreValue mismatchScore = scoreMismatch(scoringScheme);
    TScoreValue gapOpenScore = scoreGapOpen(scoringScheme);
    TScoreValue gapExtendScore = scoreGapExtend(scoringScheme);

    // Perform the matrix filling, column-wise for cache efficiency.
    for (TSequenceIterator it1 = begin(sequence1), it1End = end(sequence1); it1 != it1End; ++it1) {
        // std::cout << "iteration it1" << std::endl;
        itMLeft = itMTop;
        itIALeft = itIATop;
        itIBLeft = itIBTop;
        goNext(itMTop, 1);
        goNext(itIATop, 1);
        goNext(itIBTop, 1);
        itMAbove = itMTop;
        itIAAbove = itIATop;
        itIBAbove = itIBTop;
        for (TSequenceIterator it0 = begin(sequence0), it0End = end(sequence0); it0 != it0End; ++it0) {
            // std::cout << "  iteration it0" << std::endl;
            // Compute M_{i-1,j-1} + match/mismatch
            TScoreValue scoreMoveDiagonal = *itMLeft + ((*it0 == *it1) ? matchScore : mismatchScore);
            goNext(itMLeft, 0);
            // Compute I^a_{i,j}
            TScoreValue scoreIA = _max(*itMAbove + gapOpenScore, *itIAAbove + gapExtendScore);
            goNext(itIALeft, 0);  // TODO(holtgrew): Remove IALeft? Not necessary!
            goNext(itIAAbove, 0);
            *itIAAbove = scoreIA;
            // Compute I^b_{i,j}
            goNext(itIBLeft, 0);
            TScoreValue scoreIB = _max(*itMLeft + gapOpenScore, *itIBLeft + gapExtendScore);
            goNext(itIBAbove, 0);
            *itIBAbove = scoreIB;
            // Assign M_{i,j}
            goNext(itMAbove, 0);
            *itMAbove = _max(scoreMoveDiagonal, _max(*itIAAbove, *itIBAbove));
        }
    }

    // TODO(holtgrew): Debug code, remove when working.
    {
        for (int k = 0; k < 3; ++k) {
            std::cerr << ",-- *** filled alignment matrix " << k << std::endl;
            for (unsigned i = 0; i < length(matrix, 0); ++i) {
                for (unsigned j = 0; j < length(matrix, 1); ++j) {
                    if (value(matrix, i, j, k) <= MinValue<int>::VALUE / 4)
                        std::cerr << "\tinf";
                    else
                        std::cerr << "\t" << value(matrix, i, j, k);
                }
                std::cerr << std::endl;
            }
            std::cerr << "`--" << std::endl;
        }
    }
}

// Compute traceback in the given normal DP alignment matrix, starting
// at the lower right, shifted by the given overlap to the upper left.
// The matrix was filled with the Needleman-Wunschla lgorithm, end
// gaps are free as configured by the AlignConfig object.  Returns the
// best score.
template <typename TAlignmentIterator, typename TSequenceIterator, typename TPosition, typename TScoreValue, typename TScoringScheme, typename TOverlap, bool START0_FREE, bool START1_FREE, bool END0_FREE, bool END1_FREE>
TScoreValue
_alignTraceback(TAlignmentIterator & alignmentIt0, TAlignmentIterator & alignmentIt1, TSequenceIterator & sourceIt0, TSequenceIterator & sourceIt1, TPosition & finalPos0, TPosition & finalPos1, Matrix<TScoreValue, 3> /*const*/ & matrix, TScoringScheme const & scoringScheme, TOverlap lowerRightOverlap0, TOverlap lowerRightOverlap1, TOverlap upperLeftOverlap0, TOverlap upperLeftOverlap1, bool goToTopLeft, AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    std::cout << "unbanded traceback" << std::endl;

    // Suppress warnings about unused parameters.
    (void) scoringScheme;

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;

    // Initialization
    //
    // Current position in the matrix.  Note that this is the position
    // in the matrix including the gutter.  When writing this out, we
    // want to have the position in the matrix, excluding the gutter,
    // i.e. everything is shifted to the upper left.
    TPosition pos0 = length(matrix, 0) - lowerRightOverlap0 + (finalPos0 - 1);
    TPosition pos1 = length(matrix, 1) - lowerRightOverlap1 + (finalPos1 - 1);
    // Iterators to current entries in the matrices.
    TMatrixIterator diagonalIt = begin(matrix);
    goTo(diagonalIt, pos0, pos1);
    std::cout << "STARTING AT (" << pos0 << ", " << pos1 << ")" << std::endl;
    TMatrixIterator verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    TMatrixIterator horizontalIt = verticalIt;
    goNext(horizontalIt, 2);

    // TODO(holtgrew): At this point, we could search for the starting point of the backtrace if free end gaps were suppported.
    //
    SEQAN_ASSERT_NOT_MSG(END0_FREE, "Free end gaps are not supported yet.");
    SEQAN_ASSERT_NOT_MSG(END1_FREE, "Free end gaps are not supported yet.");

    // Flags for movement, also determines which matrix we are in.
    // TODO(holtgrew): Maybe better use an enum here?
    bool horizontal = false;
    bool vertical = false;
    bool diagonal = false;
    TScoreValue result = *diagonalIt;
    if (*diagonalIt > *horizontalIt) {
        if (*diagonalIt > *verticalIt)
            diagonal = true;
        else
            vertical = true;
    } else {
        if (*horizontalIt > *verticalIt)
            horizontal = true;
        else
            vertical = true;
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
            std::cerr << "pos0 == " << pos0 << " pos1 == " << pos1 << std::endl;
            if ((pos0 <= static_cast<TPosition>(1) || pos1 <= static_cast<TPosition>(1)) || (pos0 < static_cast<TPosition>(upperLeftOverlap0) && pos1 < static_cast<TPosition>(upperLeftOverlap1)))
                break;
        }

        if (diagonal) {
            std::cout << "DIAGONAL" << std::endl;
            // Move iterators in sequences, alignment rows and matrices.
            goPrevious(sourceIt0);
            goPrevious(sourceIt1);
            goPrevious(alignmentIt0);
            goPrevious(alignmentIt1);
            goPrevious(diagonalIt, 0);
            goPrevious(diagonalIt, 1);
            goPrevious(horizontalIt, 0);
            goPrevious(horizontalIt, 1);
            goPrevious(verticalIt, 0);
            goPrevious(verticalIt, 1);
            SEQAN_ASSERT_GEQ(pos0, 1u);
            SEQAN_ASSERT_GEQ(pos1, 1u);
            pos0 -= 1;
            pos1 -= 1;

            // Update the movement/matrix indicators.
            // TODO(holtgrew): This code will prefer gaps over matches, should be changed.
            if (*diagonalIt > *horizontalIt) {
                if (*diagonalIt <= *verticalIt) {
                    vertical = true;
                    diagonal = false;
                }
            } else {
                diagonal = false;
                if (*horizontalIt > *verticalIt)
                    horizontal = true;
                else
                    vertical = true;
            }
        } else if (vertical) {
            std::cout << "VERTICAL" << std::endl;
            // Insert gap.
            insertGap(alignmentIt1);

            // Move iterators in sequence, alignment rows and matrices.,
            goPrevious(sourceIt0);
            goPrevious(alignmentIt0);
            goPrevious(diagonalIt, 0);
            goPrevious(verticalIt, 0);
            goPrevious(horizontalIt, 0);
            pos0 -= 1;

            // std::cout << "  *verticalIt == " << *verticalIt << ", *horizontalIt == " << *horizontalIt << ", *diagonalIt == " << *diagonalIt << std::endl;

            // Update the movement/matrix indicators.
            if (*verticalIt >= *horizontalIt) {
                if (*diagonalIt > *verticalIt) {
                    diagonal = true;
                    vertical = false;
                }
            } else {
                vertical = false;
                if (*diagonalIt > *horizontalIt)
                    diagonal = true;
                else
                    horizontal = true;
            }
        } else {
            SEQAN_ASSERT(horizontal);
            std::cout << "HORIZONTAL" << std::endl;
            // Insert gap.
            insertGap(alignmentIt0);
                    
            // Move iterators in sequence, alignment rows and matrices.,
            goPrevious(sourceIt1);
            goPrevious(alignmentIt1);
            goPrevious(diagonalIt, 1);
            goPrevious(verticalIt, 1);
            goPrevious(horizontalIt, 1);
            pos1 -= 1;

            // std::cout << "  *verticalIt == " << *verticalIt << ", *horizontalIt == " << *horizontalIt << ", *diagonalIt == " << *diagonalIt << std::endl;
                    
            // Move iterators in sequence, alignment rows and matrices.,
            // Update the movement/matrix indicators.
            if (*horizontalIt >= *verticalIt) {
                if (*diagonalIt > *horizontalIt) {
                    diagonal = true;
                    horizontal = false;
                }
            } else {
                horizontal = false;
                if (*diagonalIt > *verticalIt)
                    diagonal = true;
                else
                    vertical = true;
            }
        }
    }
    // std::cout << "pos0 == " << pos0 << ", pos1 == " << pos1 << std::endl;

    // Go to the top left of the matrix if configured to do so.
    if (goToTopLeft) {
        if (pos0 > 0) {
            goPrevious(alignmentIt1);
            for (TPosition i = 0; i < pos0; ++i) {
                // std::cout << "Inserting " << pos0 << " gaps into alignment row 1" << std::endl;
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
            for (TPosition i = 0; i< pos1; ++i) {
                // std::cout << "Inserting " << pos0 << " gaps into alignment row 0" << std::endl;
                goPrevious(sourceIt1);
                goPrevious(alignmentIt1);
                // std::cout << "  *alignmentIt1 == " << convert<char>(*alignmentIt1) << std::endl;
                // std::cout << __FILE__ << ":" << __LINE__ << "-- gap in dimension 0" << std::endl;
                insertGap(alignmentIt0);
            }
            pos1 = 0;
        }
    }

    // Write out the final positions in the alignment matrix.  Convert
    // from position in the current alignment matrix with gutter to
    // positioin without gutter by shifting the position to the upper
    // left.
    finalPos0 = pos0;
    finalPos1 = pos1;
    std::cout << "finalPos0 = " << finalPos0 << ", finalPos1 = " << finalPos1 << std::endl;

    return result;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_AFFINE_H_

