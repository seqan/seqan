// FRAGMENT(header)
/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de
 ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
 ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ============================================================================
  Demonstration on how to initialize a scoring matrix programatically with:

   - one of the built-in matrices, here BLOSUM30
   - arbitrary values
   - a new, built-in matrix.
 ==========================================================================*/
// FRAGMENT(includes)
#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>   // For printing strings.
#include <seqan/score.h>  // The module score.

using namespace seqan;

// FRAGMENT(user-defined-matrix)
// Extend SeqAn by a user-define scoring matrix.
namespace seqan {

// We have to create a new specialization of the ScoringMatrix_ class
// for amino acids.  For this, we first create a new tag.
struct UserDefinedMatrix {};
// We also do this for the DNA alphabet.
struct AnotherUserDefinedMatrix {};

// Then, we specialize the class ScoringMatrix_.
template <>
struct ScoringMatrixData_<int, AminoAcid, UserDefinedMatrix> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        // The user defined data table.  In this case, we use the data from BLOSUM-30.
        static int const _data[TAB_SIZE] = {
             4, -1,  0,  0, -3,  1,  0,  0, -2,  0, -1,  0,  1, -2, -1,  1,  1, -5, -4,  1,  0,  0,  0, -7,
            -1,  8, -2, -1, -2,  3, -1, -2, -1, -3, -2,  1,  0, -1, -1, -1, -3,  0,  0, -1, -2,  0, -1, -7,
             0, -2,  8,  1, -1, -1, -1,  0, -1,  0, -2,  0,  0, -1, -3,  0,  1, -7, -4, -2,  4, -1,  0, -7,
             0, -1,  1,  9, -3, -1,  1, -1, -2, -4, -1,  0, -3, -5, -1,  0, -1, -4, -1, -2,  5,  0, -1, -7,
            -3, -2, -1, -3, 17, -2,  1, -4, -5, -2,  0, -3, -2, -3, -3, -2, -2, -2, -6, -2, -2,  0, -2, -7,
             1,  3, -1, -1, -2,  8,  2, -2,  0, -2, -2,  0, -1, -3,  0, -1,  0, -1, -1, -3, -1,  4,  0, -7,
             0, -1, -1,  1,  1,  2,  6, -2,  0, -3, -1,  2, -1, -4,  1,  0, -2, -1, -2, -3,  0,  5, -1, -7,
             0, -2,  0, -1, -4, -2, -2,  8, -3, -1, -2, -1, -2, -3, -1,  0, -2,  1, -3, -3,  0, -2, -1, -7,
            -2, -1, -1, -2, -5,  0,  0, -3, 14, -2, -1, -2,  2, -3,  1, -1, -2, -5,  0, -3, -2,  0, -1, -7,
             0, -3,  0, -4, -2, -2, -3, -1, -2,  6,  2, -2,  1,  0, -3, -1,  0, -3, -1,  4, -2, -3,  0, -7,
            -1, -2, -2, -1,  0, -2, -1, -2, -1,  2,  4, -2,  2,  2, -3, -2,  0, -2,  3,  1, -1, -1,  0, -7,
             0,  1,  0,  0, -3,  0,  2, -1, -2, -2, -2,  4,  2, -1,  1,  0, -1, -2, -1, -2,  0,  1,  0, -7,
             1,  0,  0, -3, -2, -1, -1, -2,  2,  1,  2,  2,  6, -2, -4, -2,  0, -3, -1,  0, -2, -1,  0, -7,
            -2, -1, -1, -5, -3, -3, -4, -3, -3,  0,  2, -1, -2, 10, -4, -1, -2,  1,  3,  1, -3, -4, -1, -7,
            -1, -1, -3, -1, -3,  0,  1, -1,  1, -3, -3,  1, -4, -4, 11, -1,  0, -3, -2, -4, -2,  0, -1, -7,
             1, -1,  0,  0, -2, -1,  0,  0, -1, -1, -2,  0, -2, -1, -1,  4,  2, -3, -2, -1,  0, -1,  0, -7,
             1, -3,  1, -1, -2,  0, -2, -2, -2,  0,  0, -1,  0, -2,  0,  2,  5, -5, -1,  1,  0, -1,  0, -7,
            -5,  0, -7, -4, -2, -1, -1,  1, -5, -3, -2, -2, -3,  1, -3, -3, -5, 20,  5, -3, -5, -1, -2, -7,
            -4,  0, -4, -1, -6, -1, -2, -3,  0, -1,  3, -1, -1,  3, -2, -2, -1,  5,  9,  1, -3, -2, -1, -7,
             1, -1, -2, -2, -2, -3, -3, -3, -3,  4,  1, -2,  0,  1, -4, -1,  1, -3,  1,  5, -2, -3,  0, -7,
             0, -2,  4,  5, -2, -1,  0,  0, -2, -2, -1,  0, -2, -3, -2,  0,  0, -5, -3, -2,  5,  0, -1, -7,
             0,  0, -1,  0,  0,  4,  5, -2,  0, -3, -1,  1, -1, -4,  0, -1, -1, -1, -2, -3,  0,  4,  0, -7,
             0, -1,  0, -1, -2,  0, -1, -1, -1,  0,  0,  0,  0, -1, -1,  0,  0, -2, -1,  0, -1,  0, -1, -7,
            -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7,  1,
        };
        return _data;
    }
};

// And we do this for the Dna5 alphabet.
template <>
struct ScoringMatrixData_<int, Dna5, AnotherUserDefinedMatrix> {
    enum {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        // The user defined data table.  In this case, we use the data from BLOSUM-30.
        static int const _data[TAB_SIZE] = {
          1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 0
        };
        return _data;
    }
};
}  // namespace seqan

// FRAGMENT(show-scoring-matrix)
// Print a scoring scheme matrix to stdout.
template <typename TScoreValue, typename TSequenceValue, typename TSpec>
void showScoringMatrix(Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & scoringScheme)
{
    // Print top row.
    for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
        std::cout << "\t" << TSequenceValue(i);
    std::cout << std::endl;
    // Print each row.
    for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i) {
        std::cout << TSequenceValue(i);
        for (unsigned j = 0; j < ValueSize<TSequenceValue>::VALUE; ++j) {
            std::cout << "\t" << score(scoringScheme, TSequenceValue(i), TSequenceValue(j));
        }
        std::cout << std::endl;
    }
}

// FRAGMENT(main)
int main()
{
    // 1. Define type and constants.
    //
    // Define types for the score value and the scoring scheme.
    typedef int TValue;
    typedef Score<TValue, ScoreMatrix<AminoAcid, Default> > TScoringScheme;
    // Define our gap scores in some constants.
    int const gapOpenScore = -1;
    int const gapExtendScore = -1;

    // 2. Construct scoring scheme with default/empty matrix.
    //
    // Construct new scoring scheme, alternatively only give one score
    // that is used for both opening and extension.
    TScoringScheme scoringScheme(gapOpenScore, gapExtendScore);

    // 3. Fill the now-existing ScoreMatrix
    //
    // The scoring scheme now already has a matrix of the size
    // ValueSize<AminoAcid>::VALUE x ValueSize<AminoAcid>::VALUE which
    // we can now fill.

    // 3.1 First, fill it with BLOSUM30.
    std::cout << "BLOSUM 30" << std::endl;
    setDefaultScoreMatrix(scoringScheme, Blosum30_());
    showScoringMatrix(scoringScheme);

    // 3.2 Now, we fill it with the product of the coordinates.
    std::cout << std::endl << "Coordinate Products" << std::endl;
    for (unsigned i = 0; i < ValueSize<AminoAcid>::VALUE; ++i) {
        for (unsigned j = 0; j < ValueSize<AminoAcid>::VALUE; ++j) {
            setScore(scoringScheme, AminoAcid(i), AminoAcid(j), i * j);
        }
    }
    showScoringMatrix(scoringScheme);

    // 3.3 Now, we fill it with the user defined matrix above.
    std::cout << "User defined matrix (also BLOSUM 30)..." << std::endl;
    setDefaultScoreMatrix(scoringScheme, UserDefinedMatrix());
    showScoringMatrix(scoringScheme);

    // 4. Create ScoreMatrix object with user-defined matrix.
    std::cout << "User scoring scheme..." << std::endl;
    Score<TValue, ScoreMatrix<AminoAcid, UserDefinedMatrix> > userScoringScheme;
    showScoringMatrix(userScoringScheme);

    // 5. Show our Dna5 scoring matrix.
    std::cout << "User DNA scoring scheme..." << std::endl;
    Score<TValue, ScoreMatrix<Dna5, AnotherUserDefinedMatrix> > userScoringSchemeDna;
    showScoringMatrix(userScoringSchemeDna);

	return 0;
}
