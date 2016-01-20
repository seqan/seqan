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
  Demonstration on how to load a score matrix from a file.
 ==========================================================================*/
//![header]
//![includes]
#include <iostream>

#include <seqan/basic.h>
#include <seqan/stream.h>   // For printing strings.
#include <seqan/score.h>  // The module score.

using namespace seqan;
//![includes]

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
int main(int argc, char ** argv)
{
    typedef int TScoreValue;

    Score<TScoreValue, ScoreMatrix<Dna5, Default> > scoreMatrix;
    loadScoreMatrix(scoreMatrix, toCString(getAbsolutePath("demos/howto/scores/dna_example.txt")));
    showScoringMatrix(scoreMatrix);

    return 0;
}
//![main]
