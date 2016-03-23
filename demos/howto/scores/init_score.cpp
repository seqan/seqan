//![header]
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
//![header]
//![includes]
#include <iostream>

#include <seqan/basic.h>
#include <seqan/stream.h>   // For printing strings.
#include <seqan/score.h>    // The module score.

using namespace seqan;
//![includes]

//![user-defined-matrix]
// Extend SeqAn by a user-define scoring matrix.
namespace seqan {

// We have to create a new specialization of the ScoringMatrix_ class
// for the DNA alphabet.  For this, we first create a new tag.
struct UserDefinedMatrix {};

// Then, we specialize the class ScoringMatrix_ for the Dna5 alphabet.
template <>
struct ScoringMatrixData_<int, Dna5, UserDefinedMatrix>
{
    enum
    {
        VALUE_SIZE = ValueSize<Dna5>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData()
    {
        // The user defined data table.  In this case, we use the data from BLOSUM-30.
        static int const _data[TAB_SIZE] =
        {
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
//![user-defined-matrix]

//![show-scoring-matrix]
// Print a scoring scheme matrix to stdout.
template <typename TScoreValue, typename TSequenceValue, typename TSpec>
void showScoringMatrix(Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & scoringScheme)
{
    // Print top row.
    for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
        std::cout << "\t" << TSequenceValue(i);
    std::cout << std::endl;
    // Print each row.
    for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
    {
        std::cout << TSequenceValue(i);
        for (unsigned j = 0; j < ValueSize<TSequenceValue>::VALUE; ++j)
        {
            std::cout << "\t" << score(scoringScheme, TSequenceValue(i), TSequenceValue(j));
        }
        std::cout << std::endl;
    }
}
//![show-scoring-matrix]

//![main]
int main()
{
    // 1. Define type and constants.
    //
    // Define types for the score value and the scoring scheme.
    typedef int TValue;
    typedef Score<TValue, ScoreMatrix<Dna5, Default> > TScoringScheme;
    // Define our gap scores in some constants.
    int const gapOpenScore = -1;
    int const gapExtendScore = -1;

    // 2. Construct scoring scheme with default/empty matrix.
    //
    // Construct new scoring scheme, alternatively only give one score
    // that is used for both opening and extension.
    TScoringScheme scoringScheme(gapExtendScore, gapOpenScore);

    // 3. Fill the now-existing ScoreMatrix
    //
    // The scoring scheme now already has a matrix of the size
    // ValueSize<Dna5>::VALUE x ValueSize<Dna5>::VALUE which
    // we can now fill.

    // 3.1 We fill the scoring scheme with the product of the coordinates.
    std::cout << std::endl << "Coordinate Products" << std::endl;
    for (unsigned i = 0; i < ValueSize<Dna5>::VALUE; ++i)
    {
        for (unsigned j = 0; j < ValueSize<Dna5>::VALUE; ++j)
        {
            setScore(scoringScheme, Dna5(i), Dna5(j), i * j);
        }
    }
    showScoringMatrix(scoringScheme);

    // 3.2 Now, we fill it with the user defined matrix above.
    std::cout << "User defined matrix (also Dna5 scoring matrix)..." << std::endl;
    setDefaultScoreMatrix(scoringScheme, UserDefinedMatrix());
    showScoringMatrix(scoringScheme);

    // 4. Show our user-defined Dna5 scoring matrix.
    std::cout << "User DNA scoring scheme..." << std::endl;
    Score<TValue, ScoreMatrix<Dna5, UserDefinedMatrix> > userScoringSchemeDna;
    showScoringMatrix(userScoringSchemeDna);

    return 0;
}
//![main]
