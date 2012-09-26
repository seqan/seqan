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

// TODO(holtgrew): Maybe adjust rausch's code so things can be copied in and use it instead since it is heavily tuned?

#ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_AFFINE_H_
#define SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_AFFINE_H_

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

template <typename TScoreValue, typename TSequence, typename TDiagonal>
inline void
_alignBandedResizeMatrix(Matrix<TScoreValue, 3> & matrix, TSequence const & sequence0, TSequence const & sequence1, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    (void) sequence1;  // In case we do not run in debug mode.

    SEQAN_ASSERT_GEQ(upperDiagonal, 0);
    SEQAN_ASSERT_LEQ(lowerDiagonal, 0);
    SEQAN_ASSERT_GEQ_MSG(length(sequence1) - lowerDiagonal, length(sequence0), "Lower diagonal is not low enough.");
    SEQAN_ASSERT_GEQ_MSG(length(sequence0) + upperDiagonal, length(sequence1), "Upper diagonal is not high enough.");

    // We need space length(sequence0) x bandwidth for the alignment,
    // the top gutter and the left and right gutters.
    setLength(matrix, 0, length(sequence0) + 1);
    setLength(matrix, 1, upperDiagonal - lowerDiagonal + 3);
    setLength(matrix, 2, 3);
    resize(matrix);
}


template <typename TScoreValue, typename TDiagonal, bool BEGIN1_FREE, bool BEGIN0_FREE, bool END1_FREE, bool END0_FREE>
inline void
_alignBandedInitGutter(Matrix<TScoreValue, 3> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, AlignConfig<BEGIN1_FREE, BEGIN0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;
    
    // Initialize the diagonal below the lower one with infimas.
    TIterator diagonalIt = begin(matrix);
    TIterator verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    TIterator horizontalIt = verticalIt;
    goNext(horizontalIt, 2);
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(diagonalIt, 0), goNext(horizontalIt, 0)) {
        *diagonalIt = MinValue<TScoreValue>::VALUE / 2;
        *horizontalIt = MinValue<TScoreValue>::VALUE / 2;
    }
    // Initialize the left gutter according to the AlignConfig.
    setPosition(diagonalIt, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    horizontalIt = verticalIt;
    goNext(horizontalIt, 2);
    if (BEGIN0_FREE) {
        for (TPosition i = 0, iend = -lowerDiagonal; i <= iend; ++i, goNext(diagonalIt, 0), goPrevious(diagonalIt, 1)) {
            *diagonalIt = 0;
        }
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        for (TPosition i = 0, iend = -lowerDiagonal; i <= iend; ++i, goNext(diagonalIt, 0), goPrevious(diagonalIt, 1)) {
            *diagonalIt = x;
            x += gapScore;
        }
    }
    for (TPosition i = 0, iend = -lowerDiagonal; i <= iend; ++i, goNext(verticalIt, 0), goPrevious(verticalIt, 1), goNext(horizontalIt, 0), goPrevious(horizontalIt, 1)) {
        *verticalIt = MinValue<TScoreValue>::VALUE / 2;
        *horizontalIt = MinValue<TScoreValue>::VALUE / 2;
    }
    // Initialize the top gutter according to the AlignConfig.
    setPosition(diagonalIt, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    horizontalIt = verticalIt;
    goNext(horizontalIt, 2);
    if (BEGIN1_FREE) {
        for (TPosition i = 0, iend = upperDiagonal; i <= iend; ++i, goNext(diagonalIt, 1))
            *diagonalIt = 0;
    } else {
        TScoreValue gapScore = scoreGap(scoringScheme);
        TScoreValue x = 0;
        for (TPosition i = 0, iend = upperDiagonal; i <= iend; ++i, goNext(diagonalIt, 1)) {
            *diagonalIt = x;
            x += gapScore;
        }
    }
    for (TPosition i = 0, iend = upperDiagonal; i <= iend; ++i, goNext(verticalIt, 1), goNext(horizontalIt, 1)) {
        *verticalIt = MinValue<TScoreValue>::VALUE / 2;
        *horizontalIt = MinValue<TScoreValue>::VALUE / 2;
    }
    // Initialize the diagonal above the upper one with infimas.
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(diagonalIt, 0), goNext(verticalIt, 0)) {
        *diagonalIt = MinValue<TScoreValue>::VALUE / 2;
        *verticalIt = MinValue<TScoreValue>::VALUE / 2;
    }
}


template <typename TScoreValue, typename TDiagonal, typename TOverlap>
inline void
_alignBandedInitGutterFromUnbanded(Matrix<TScoreValue, 3> & matrix, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Matrix<TScoreValue, 3> /*const*/ & otherMatrix, TOverlap overlap0, TOverlap overlap1, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    // TODO(holtgrew): Really unnecessary? Remove along with all other unused parameters in all align_*.h files.
    (void) scoringScheme;
    (void) upperDiagonal;

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TIterator;
    typedef typename Position<TMatrix>::Type TPosition;

    // Copy over left gutter column.
    TIterator diagonalIt = begin(matrix);
    goTo(diagonalIt, 0, 1 - lowerDiagonal);
    TIterator verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    TIterator horizontalIt = verticalIt;
    goNext(horizontalIt, 2);
    TIterator srcDiagonalIt = begin(otherMatrix);
    goTo(srcDiagonalIt, length(otherMatrix, 0) - overlap0 - 1, length(otherMatrix, 1) - overlap1 - 1);
    TIterator srcVerticalIt = srcDiagonalIt;
    goNext(srcVerticalIt, 2);
    TIterator srcHorizontalIt = srcVerticalIt;
    goNext(srcHorizontalIt, 2);
    for (TOverlap i = 0; i < overlap0 + 1; ++i) {
        *diagonalIt = *srcDiagonalIt;
        // *verticalIt = *srcVerticalIt;
        *horizontalIt = *srcHorizontalIt;
        goNext(srcDiagonalIt, 0);
        // goNext(srcVerticalIt, 0);
        goNext(srcHorizontalIt, 0);
        goNext(diagonalIt, 0);
        goPrevious(diagonalIt, 1);
        // goNext(verticalIt, 0);
        // goPrevious(verticalIt, 1);
        goNext(horizontalIt, 0);
        goPrevious(horizontalIt, 1);
    }
    // Initialize the diagonal below the lower one with infimas.
    diagonalIt = begin(matrix);
    verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    horizontalIt = verticalIt;
    goNext(horizontalIt, 2);
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(diagonalIt, 0), goNext(horizontalIt, 0)) {
        *diagonalIt = MinValue<TScoreValue>::VALUE / 2;
        *horizontalIt = MinValue<TScoreValue>::VALUE / 2;
    }

    // Copy over a row from the other matrix into the top gutter.
    goTo(diagonalIt, 0, 1 - lowerDiagonal);
    verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    horizontalIt = verticalIt;
    goNext(horizontalIt, 2);
    srcDiagonalIt = begin(otherMatrix);
    goTo(srcDiagonalIt, length(otherMatrix, 0) - overlap0 - 1, length(otherMatrix, 1) - overlap1 - 1);
    srcVerticalIt = srcDiagonalIt;
    goNext(srcVerticalIt, 2);
    // srcHorizontalIt = srcVerticalIt;
    // goNext(srcHorizontalIt, 2);
    for (TOverlap i = 0; i < overlap1; ++i) {
        *diagonalIt = *srcDiagonalIt;
        *verticalIt = *srcVerticalIt;
        // *horizontalIt = *srcHorizontalIt;
        goNext(srcDiagonalIt, 1);
        goNext(srcVerticalIt, 1);
        // goNext(srcHorizontalIt, 1);
        goNext(diagonalIt, 1);
        goNext(verticalIt, 1);
        // goNext(horizontalIt, 1);
    }

    // Initialize the diagonal above the upper one with infimas.
    for (TPosition i = 0, iend = length(matrix, 0); i < iend; ++i, goNext(diagonalIt, 0), goNext(verticalIt, 0)) {
        *diagonalIt = MinValue<TScoreValue>::VALUE / 2;
        *verticalIt = MinValue<TScoreValue>::VALUE / 2;
    }
}


template <typename TScoreValue, typename TSequence, typename TDiagonal>
inline void
_alignBandedFillMatrix(Matrix<TScoreValue, 3> & matrix, TSequence const & sequence0, TSequence const & sequence1, Score<TScoreValue, Simple> const & scoringScheme, TDiagonal lowerDiagonal, TDiagonal upperDiagonal, Gotoh const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GEQ_MSG(upperDiagonal, 0, "Upper diagonal must not lie below main diagonal.");
    SEQAN_ASSERT_LEQ_MSG(lowerDiagonal, 0, "Upper diagonal must not lie above main diagonal.");

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;
    typedef typename Iterator<TSequence, Standard>::Type TSequenceIterator;
    typedef typename Size<TSequence>::Type TSize;
    typedef typename Position<TSequence>::Type TPosition;

    // We need three iterators in each alignment matrix to fill it.
    // it*Top points to the cell in the top row of the current column.
    // it*Left points to the column to the top left of the current
    // cell.  it*Above points to the cell above the current cell.  We
    // can use itAbove for value assignment.
    TMatrixIterator itMTop = begin(matrix);
    TMatrixIterator itMAbove;
    TMatrixIterator itMLeft;
    TMatrixIterator itIATop;
    TMatrixIterator itIAAbove;
    TMatrixIterator itIALeft;
    TMatrixIterator itIBTop;
    TMatrixIterator itIBAbove;
    TMatrixIterator itIBLeft;

    // Compute score values so they are available locally and it maybe
    // is easier to the compiler to save some indirect adressing.
    TScoreValue matchScore = scoreMatch(scoringScheme);
    TScoreValue mismatchScore = scoreMismatch(scoringScheme);
    TScoreValue gapOpenScore = scoreGapOpen(scoringScheme);
    TScoreValue gapExtendScore = scoreGapExtend(scoringScheme);

    // Fill matrix column wise for cache efficiency.  For this, we
    // need the length of the current column in the sheared alignment
    // matrix which complicates things.
    setPosition(itMTop, (1 - lowerDiagonal) * _dataFactors(matrix)[1]); // TODO(holtgrew): There should be a function that accepts two coordinates for the Matrix class.
    itIATop = itMTop;
    goNext(itIATop, 2);
    itIBTop = itIATop;
    goNext(itIBTop, 2);
    TPosition seq1Pos = 0;
    TSequenceIterator it0Begin = begin(sequence0);
    for (TSequenceIterator it1 = begin(sequence1), it1end = end(sequence1); it1 != it1end; ++it1, ++seq1Pos) {
        if (seq1Pos <= static_cast<TPosition>(upperDiagonal)) {
            itMLeft = itMTop;
            goNext(itMTop, 1);
            itIALeft = itIATop;
            goNext(itIATop, 1);
            itIBLeft = itIBTop;
            goNext(itIBTop, 1);
        } else {
            goNext(itMTop, 0);
            itMLeft = itMTop;
            goPrevious(itMLeft, 1);
            goNext(itIATop, 0);
            itIALeft = itIATop;
            goPrevious(itIALeft, 1);
            goNext(itIBTop, 0);
            itIBLeft = itIBTop;
            goPrevious(itIBLeft, 1);
        }
        itMAbove = itMTop;
        itIAAbove = itIATop;
        itIBAbove = itIBTop;
        TDiagonal from = _max(static_cast<TDiagonal>(seq1Pos) - upperDiagonal, TDiagonal(0));
        TDiagonal to = _min(seq1Pos - lowerDiagonal, length(sequence0) - 1);
        TDiagonal tmp = 1 + to - from;
        TSize runLength = _max(TDiagonal(0), tmp);
        // std::cout << "RUN LENGTH == " << runLength << std::endl;
        // std::cout << "from == " << from << std::endl << "to == " << to << std::endl;
        // std::cout << "seq1Pos = " << seq1Pos << std::endl;
        // std::cout << "lowerDiagonal = " << lowerDiagonal << std::endl;
        // std::cout << "length(sequence0) = " << length(sequence0) << std::endl;
        // std::cout << "seq1Pos = " << seq1Pos << std::endl;
        // std::cout << "upperDiagonal = " << upperDiagonal << std::endl;
        TSize i = 0;
        for (TSequenceIterator it0 = it0Begin; i < runLength; ++it0, ++i) {
            // Compute M_{i-1,j-1} + match/mismatch
            TScoreValue scoreMoveDiagonal = *itMLeft + ((*it0 == *it1) ? matchScore : mismatchScore);
            goNext(itMLeft, 0);
            goPrevious(itMLeft, 1);
            // Compute I^a_{i,j}
            TScoreValue scoreIA = _max(*itMAbove + gapOpenScore, *itIAAbove + gapExtendScore);
            goNext(itIALeft, 0);  // TODO(holtgrew): Remove IALeft? Not necessary!
            goPrevious(itIALeft, 1);
            goNext(itIAAbove, 0);
            goPrevious(itIAAbove, 1);
            *itIAAbove = scoreIA;
            // Compute I^b_{i,j}
            goNext(itIBLeft, 0);
            goPrevious(itIBLeft, 1);
            TScoreValue scoreIB = _max(*itMLeft + gapOpenScore, *itIBLeft + gapExtendScore);
            goNext(itIBAbove, 0);
            goPrevious(itIBAbove, 1);
            *itIBAbove = scoreIB;
            // Assign M_{i,j}
            goNext(itMAbove, 0);
            goPrevious(itMAbove, 1);
            *itMAbove = _max(scoreMoveDiagonal, _max(*itIAAbove, *itIBAbove));
        }
        // std::cout << "++it1" << std::endl;
        // We only need an offset for it0 when all diagonals up to the
        // upper one are aligned.
        if (seq1Pos >= static_cast<TPosition>(upperDiagonal))
            it0Begin += 1;
    }

    // // TODO(holtgrew): Debug code, remove when working.
    // {
    //     for (int k = 0; k < 3; ++k) {
    //         std::cerr << ",-- *** filled banded alignment matrix " << k << std::endl;
    //         for (unsigned i = 0; i < length(matrix, 0); ++i) {
    //             std::cerr << "| ";
    //             for (unsigned j = 0; j < i; ++j)
    //                 std::cerr << "\t";
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


// Compute traceback in the given banded alignment matrix, starting at
// the lower right, shifted by the given overlap to the upper left.
// The matrix was filled with the Needleman-Wunschla lgorithm, end
// gaps are free as configured by the AlignConfig object.  Returns the
// best score.
template <typename TAlignmentIterator, typename TSequenceIterator, typename TPosition, typename TScoreValue, typename TScoringScheme, typename TOverlap, bool START0_FREE, bool START1_FREE, bool END0_FREE, bool END1_FREE>
TScoreValue
_alignBandedTraceback(TAlignmentIterator & alignmentIt0, TAlignmentIterator & alignmentIt1, TSequenceIterator & sourceIt0, TSequenceIterator & sourceIt1, TPosition & finalPos0, TPosition & finalPos1, Matrix<TScoreValue, 3> /*const*/ & matrix, TScoringScheme const & scoringScheme, TOverlap overlap0, TOverlap overlap1, TOverlap upperTriangleEdgeLength, TOverlap lowerTriangleEdgeLength, bool goToTopLeft, AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const &, Gotoh const &)
{
    SEQAN_CHECKPOINT;
    // std::cout << "banded traceback" << std::endl;

    // Suppress unused parameters warning.
    (void) goToTopLeft;
    (void) scoringScheme;

    typedef Matrix<TScoreValue, 3> TMatrix;
    typedef typename Iterator<TMatrix>::Type TMatrixIterator;

    // Initialization
    //
    // Current position in the matrix.  Note that this is the position
    // in the matrix including the gutter.  When writing this out, we
    // want to have the position in the matrix, excluding the gutter,
    // i.e. everything is shifted to the upper left.
    TPosition pos0 = length(matrix, 0) - overlap0 + (finalPos0 - 1);
    TPosition pos1 = length(matrix, 1) - overlap1 - lowerTriangleEdgeLength - 1 + (length(matrix, 0) - pos0 - 1) + (finalPos1 - 1);
    // Iterators to current entries in the matrices.
    TMatrixIterator diagonalIt = begin(matrix);
    setPosition(diagonalIt, pos0 + pos1 * _dataFactors(matrix)[1]);  // TODO(holtgrew): Matrix class should have setPosition with coordinates.
    std::cout << "STARTING AT (" << pos0 << ", " << pos1 << ")" << std::endl;
    TMatrixIterator verticalIt = diagonalIt;
    goNext(verticalIt, 2);
    TMatrixIterator horizontalIt = verticalIt;
    goNext(horizontalIt, 2);

    // Search for starting point of the trace.
    //
    // Insert free end gaps if any.
    //
    // TODO(holtgrew): Currently not supported.
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
    // std::cerr << "before loop " << __LINE__ << std::endl;
    while (true) {
        // std::cerr << "in loop " << __LINE__ << std::endl;
        // std::cerr << "pos0 = " << pos0 << " pos1 = " << pos1 << std::endl;
        // std::cerr << "upperTriangleEdgeLength = " << upperTriangleEdgeLength << std::endl;
        if (pos0 <= static_cast<TPosition>(1))
          break;
        if (pos0 <= static_cast<TPosition>(upperTriangleEdgeLength) && pos1 <= static_cast<TPosition>(upperTriangleEdgeLength + 1 - pos0))
          break;

        SEQAN_ASSERT_GT(pos1, static_cast<TPosition>(0));

        // Determine
        if (diagonal) {
            // std::cout << "DIAGONAL" << std::endl;
            // Move iterators in sequences, alignment rows and matrices.
            goPrevious(sourceIt0);  // XXX
            goPrevious(sourceIt1);  // XXX
            goPrevious(alignmentIt0);  // XXX
            goPrevious(alignmentIt1);  // XXX
            goPrevious(diagonalIt, 0);
            // goPrevious(diagonalIt, 1);
            goPrevious(horizontalIt, 0);
            // goPrevious(horizontalIt, 1);
            goPrevious(verticalIt, 0);
            // goPrevious(verticalIt, 1);
            pos0 -= 1;
            // pos1 -= 1;

            // Update the movement/matrix indicators.
            // TODO(holtgrew): Do not access {diagonal,horizontal,vertical}It if invalid!
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
            insertGap(alignmentIt0);

            // Move iterators in sequence, alignment rows and matrices.,
            goPrevious(sourceIt0);
            goPrevious(alignmentIt0);
            goPrevious(diagonalIt, 0);
            goNext(diagonalIt, 1);  // XXX added
            goPrevious(verticalIt, 0);
            goNext(verticalIt, 1);  // XXX added
            goPrevious(horizontalIt, 0);
            goNext(horizontalIt, 1);  // XXX added
            pos0 -= 1;
            pos1 += 1;

            // Update the movement/matrix indicators.
            // TODO(holtgrew): Do not access {diagonal,horizontal,vertical}It if invalid!
            // TODO(holtgrew): This code will prefer gaps over matches, should be changed.
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
            insertGap(alignmentIt1);
                    
            // Move iterators in sequence, alignment rows and matrices.,
            goPrevious(sourceIt1);
            goPrevious(alignmentIt1);
            goPrevious(diagonalIt, 1);
            goPrevious(verticalIt, 1);
            goPrevious(horizontalIt, 1);
            pos1 -= 1;
                    
            // Move iterators in sequence, alignment rows and matrices.,
            // Update the movement/matrix indicators.
            // TODO(holtgrew): Do not access {diagonal,horizontal,vertical}It if invalid!
            // TODO(holtgrew): This code will prefer gaps over matches, should be changed.
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

    // Go to the top left of the matrix if configured to do so.
    //
    // TODO(holtgrew): Currently not supported.
    SEQAN_ASSERT_NOT_MSG(goToTopLeft, "goToTopLeft is currently not supported in banded alignment.");

    // Write out the final positions in the alignment matrix.  Convert
    // from position in the current alignment matrix with gutter to
    // positioin with gutter.  Adjust for the shear of the matrix.
    // std::cout << "@and: pos0 = " << pos0 << ", pos1 = " << pos1 << std::endl;
    finalPos0 = pos0;
    finalPos1 = pos1 - 1 - (upperTriangleEdgeLength - pos0);
    //SEQAN_ASSERT_EQ_MSG(finalPos0, 1u, "Must reach top of matrix in banded alignment traceback.");
#if SEQAN_DEBUGGING
    if (pos0 > static_cast<TPosition>(1) &&
        (pos0 > upperTriangleEdgeLength ||
         pos1 > static_cast<TPosition>(upperTriangleEdgeLength + 1 - pos0)))
        SEQAN_ASSERT_FAIL("Must reach upper left of matrix in banded alignment traceback.");
#endif  // #if SEQAN_DEBUGGING
    std::cout << "upper triangle length = " << upperTriangleEdgeLength << std::endl;
    std::cout << "finalPos0 = " << finalPos0 << ", finalPos1 = " << finalPos1 << std::endl;

    return result;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_DYNPROG_BANDED_AFFINE_H_

