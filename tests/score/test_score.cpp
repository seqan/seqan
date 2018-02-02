// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrew@fu-berlin.de>
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Tests for the SeqAn module score.
// ==========================================================================

// TODO(holtgrew): There could be one test header for each module header.
// TODO(holtgrew): Write-out-read-in only necessary once.

#include <iostream>
#include <fstream>
#include <sstream>

#define SEQAN_DEBUG


#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/score.h>
// TODO(rmaerker): Remove when done with the tests.
#include <seqan/align.h>

using namespace std;
using namespace seqan;

// Helper function that compares two amino acid matrices for equality.
// TODO(holtgrew): If used somewhere else, put into some place to share.
template <typename TScore1, typename TScore2>
void assertAminoAcidMatricesAreEqual(TScore1 const & mat1, TScore2 const & mat2) {
    for (AminoAcid a = 'A'; a <= '*'; ++a) {
        for (AminoAcid b = 'A'; b <= '*'; ++b) {
            SEQAN_ASSERT_EQ_MSG(
                    score(mat1, a, b), score(mat2, a, b),
                    "a = %c, b = %c", static_cast<char>(a),
                    static_cast<char>(b));
        }
    }
}


// Performs SEQAN_ASSERT_EQ(0, x) for all cells x in matrix.
template<typename TScore>
void assertAminoAcidMatrixIsDefaultInitialized(const TScore &matrix) {
    for (AminoAcid a = 'A'; a <= '*'; ++a) {
        for (AminoAcid b = 'A'; b <= '*'; ++b) {
            SEQAN_ASSERT_EQ(0, score(matrix, a, b));
        }
    }
}


// Test the default implementation for scoreGap* from score_base.h.
SEQAN_DEFINE_TEST(test_score_gap_open) {
    // TODO(holtgrew): The following *crashes* with edit distance score!  References to temporaries returned.
    // We simply test this with an simple score, test default implementation.
    SimpleScore simpleScore;
    const CharString kSeq1 = "A";
    const CharString kSeq2 = "A";
    typedef Position<CharString>::Type TPos;
    const TPos kPos1 = 0;
    const TPos kPos2 = 0;

    SEQAN_ASSERT_EQ(scoreGapOpen(simpleScore),
                    scoreGapOpenHorizontal(simpleScore, sequenceEntryForScore(simpleScore, kSeq1, kPos1),
                                           sequenceEntryForScore(simpleScore, kSeq2, kPos2)));
    SEQAN_ASSERT_EQ(scoreGapOpen(simpleScore),
                    scoreGapOpenVertical(simpleScore, sequenceEntryForScore(simpleScore, kSeq1, kPos1),
                                         sequenceEntryForScore(simpleScore, kSeq2, kPos2)));
    SEQAN_ASSERT_EQ(scoreGapExtend(simpleScore),
                    scoreGapExtendHorizontal(simpleScore, sequenceEntryForScore(simpleScore, kSeq1, kPos1),
                                             sequenceEntryForScore(simpleScore, kSeq2, kPos2)));
    SEQAN_ASSERT_EQ(scoreGapExtend(simpleScore),
                    scoreGapExtendVertical(simpleScore, sequenceEntryForScore(simpleScore, kSeq1, kPos1),
                                           sequenceEntryForScore(simpleScore, kSeq2, kPos2)));
    SEQAN_ASSERT_EQ(scoreGap(simpleScore),
                    scoreGapHorizontal(simpleScore, sequenceEntryForScore(simpleScore, kSeq1, kPos1),
                                       sequenceEntryForScore(simpleScore, kSeq2, kPos2)));
    SEQAN_ASSERT_EQ(scoreGap(simpleScore),
                    scoreGapVertical(simpleScore, sequenceEntryForScore(simpleScore, kSeq1, kPos1),
                                     sequenceEntryForScore(simpleScore, kSeq2, kPos2)));
    SEQAN_ASSERT_EQ(score(simpleScore, kSeq1[kPos1], kSeq2[kPos2]),
                    score(simpleScore, sequenceEntryForScore(simpleScore, kSeq1, kPos1),
                          sequenceEntryForScore(simpleScore, kSeq2, kPos2)));
}


// Compare the built-in Blosum matrices against matrices from
// reference files.
SEQAN_DEFINE_TEST(test_score_matrix) {
    // We test the scoring matrix with amino acids as the underlying
    // sequence.  This should be the most common case.
    typedef int TValue;
    typedef Score<TValue, ScoreMatrix<AminoAcid, Default> > TScore;

    // Define path to BLOSUM62 matrix that we want to load.
    // TODO(holtgrew): It should be easier to construct these paths.
    String<char> pathToTestSrc = getAbsolutePath("/tests/score/");
    String<char> pathToBlosum62(pathToTestSrc);
    append(pathToBlosum62, "BLOSUM62");

    // Test with default constructor.
    {
        // Call appropriate constructor.
        TScore matrixScore;
        // Assert score state.
        SEQAN_ASSERT_EQ(-1, scoreGap(matrixScore));
        SEQAN_ASSERT_EQ(-1, scoreGapExtend(matrixScore));
        SEQAN_ASSERT_EQ(-1, scoreGapOpen(matrixScore));
        // Test function score().  The default score is TValue() == 0
        // for all matches.
        assertAminoAcidMatrixIsDefaultInitialized(matrixScore);
    }

    // Test the const accessor functions.
    {
        // Call appropriate constructor.
        const TScore matrixScore;
        // Assert score state.
        SEQAN_ASSERT_EQ(-1, scoreGap(matrixScore));
        SEQAN_ASSERT_EQ(-1, scoreGapExtend(matrixScore));
        SEQAN_ASSERT_EQ(-1, scoreGapOpen(matrixScore));
        // Test function score().  The default score is TValue() == 0
        // for all matches.
        assertAminoAcidMatrixIsDefaultInitialized(matrixScore);
    }

    // Test the setter functions.
    {
        // The test score value for the test.
        const TValue kScoreValue = 42;
        // Create a new score matrix.
        TScore matrixScore;
        // Test function score().  The default score is TValue() == 0
        // for all matches.
        assertAminoAcidMatrixIsDefaultInitialized(matrixScore);
        // Set the score for a substitution and check the result.
        setScore(matrixScore, AminoAcid('A'), AminoAcid('X'), kScoreValue);
        SEQAN_ASSERT_EQ(kScoreValue,
                        score(matrixScore, AminoAcid('A'), AminoAcid('X')));
    }

    // Test with gap extension constructor.
    {
        // Define constant test data and call appropriate constructor.
        const TValue kGapScore = 1;
        TScore matrixScore(kGapScore);
        // Assert score state.
        SEQAN_ASSERT_EQ(kGapScore, scoreGap(matrixScore));
        SEQAN_ASSERT_EQ(kGapScore, scoreGapExtend(matrixScore));
        SEQAN_ASSERT_EQ(kGapScore, scoreGapOpen(matrixScore));
        // Test function score().  The default score is TValue() == 0
        // for all matches.
        assertAminoAcidMatrixIsDefaultInitialized(matrixScore);
    }

    // Test with gap extension and gap open constructor.
    {
        // Define constant test data and call appropriate constructor.
        const TValue kGapExtensionScore = 1;
        const TValue kGapOpenScore = 2;
        TScore matrixScore(kGapExtensionScore, kGapOpenScore);
        // Assert score state.
        SEQAN_ASSERT_EQ(kGapExtensionScore, scoreGap(matrixScore));
        SEQAN_ASSERT_EQ(kGapExtensionScore, scoreGapExtend(matrixScore));
        SEQAN_ASSERT_EQ(kGapOpenScore, scoreGapOpen(matrixScore));
        // Test function score().  The default score is TValue() == 0
        // for all matches.
        assertAminoAcidMatrixIsDefaultInitialized(matrixScore);
    }

    // Test with path to file constructor.
    {
        // Call appropriate constructor.
        const TScore matrixScore(toCString(pathToBlosum62));
        // Assert score state.
        SEQAN_ASSERT_EQ(-1, scoreGap(matrixScore));
        SEQAN_ASSERT_EQ(-1, scoreGapExtend(matrixScore));
        SEQAN_ASSERT_EQ(-1, scoreGapOpen(matrixScore));
        // The resulting matrix from the file should be equal to the
        // built-in BLOSUM62 matrix.
        Blosum62 blosum62;
        assertAminoAcidMatricesAreEqual(blosum62, matrixScore);
    }

    // Test with path to file, gap extension constructor.
    {
        // Call appropriate constructor.
        const TValue kGapScore = 1;
        TScore matrixScore(toCString(pathToBlosum62), kGapScore);
        // Assert score state.
        SEQAN_ASSERT_EQ(kGapScore, scoreGap(matrixScore));
        SEQAN_ASSERT_EQ(kGapScore, scoreGapExtend(matrixScore));
        SEQAN_ASSERT_EQ(kGapScore, scoreGapOpen(matrixScore));
        // The resulting matrix from the file should be equal to the
        // built-in BLOSUM62 matrix.
        Blosum62 blosum62;
        assertAminoAcidMatricesAreEqual(blosum62, matrixScore);
    }

    // Test with path to file, gap extension and gap open constructor.
    {
        // Define constant test data and call appropriate constructor.
        const TValue kGapExtensionScore = 1;
        const TValue kGapOpenScore = 2;
        TScore matrixScore(toCString(pathToBlosum62), kGapExtensionScore, kGapOpenScore);
        // Assert score state.
        SEQAN_ASSERT_EQ(kGapExtensionScore, scoreGap(matrixScore));
        SEQAN_ASSERT_EQ(kGapExtensionScore, scoreGapExtend(matrixScore));
        SEQAN_ASSERT_EQ(kGapOpenScore, scoreGapOpen(matrixScore));
        // The resulting matrix from the file should be equal to the
        // built-in BLOSUM62 matrix.
        Blosum62 blosum62;
        assertAminoAcidMatricesAreEqual(blosum62, matrixScore);
    }
}


// Test the File I/O code for score matrices.
SEQAN_DEFINE_TEST(test_score_matrix_file) {
    // TODO(holtgrew): It should be easier to construct these paths.
    // The path to the directory with the test's sources and fixtures.
    String<char> pathToTestSrc = getAbsolutePath("/tests/score/");

    // Load fixture BLOSUM62 matrix.
    // TODO(holtgrew): Should be done in a function.
    Score<int, ScoreMatrix<> > sc;
    String<char> meta;
    String<char> pathToBlosum62(pathToTestSrc);
    append(pathToTestSrc, "BLOSUM62");
    // TODO(holtgrew): If the file does not exist, a bus error occurs, code should catch this case and print an error.
    loadScoreMatrix(meta, sc, toCString(pathToTestSrc));

    // Compare fixture BLOSUM62 matrix to built-in one.
    {
        Blosum62 blosum62;
        assertAminoAcidMatricesAreEqual(blosum62, sc);
    }

    // Perform assertions on scores of the built-in Blosum62 matrix.
    {
        Blosum62 blosum62;
        SEQAN_ASSERT_EQ(scoreGapExtend(blosum62), -1);
        SEQAN_ASSERT_EQ(scoreGapOpen(blosum62), scoreGapExtend(blosum62));
        SEQAN_ASSERT_EQ(scoreGap(blosum62), scoreGapExtend(blosum62));
    }

    // Store and load the fixture BLOSUM62 matrix again.
    {
        const char *temp_filename = SEQAN_TEMP_FILENAME();
        {
            std::ofstream fl(temp_filename);
            typename DirectionIterator<std::ofstream, Output>::Type it = directionIterator(fl, Output());
            write(it, sc, meta);
        }

        Score<int, ScoreMatrix<> > sc2;
        String<char> meta2;
        SEQAN_ASSERT(loadScoreMatrix(meta2, sc2, temp_filename));
        assertAminoAcidMatricesAreEqual(sc, sc2);
        SEQAN_ASSERT_EQ(meta, meta2);
    }

    // Store and load the built-in matrix again.
    {
        const char *temp_filename = SEQAN_TEMP_FILENAME();
        {
            std::ofstream fl(temp_filename);
            typename DirectionIterator<std::ofstream, Output>::Type it = directionIterator(fl, Output());
            write(it, Blosum62());
        }

        Score<int, ScoreMatrix<> > sc2;
        SEQAN_ASSERT(loadScoreMatrix(sc2, temp_filename));
        assertAminoAcidMatricesAreEqual(sc2, Blosum62());
    }

    // Test setScore()
    {
        setScore(sc, 'A', '*', 100);
        SEQAN_ASSERT_EQ(score(sc, 'A', '*'), 100);
    }
}


// Testing the edit distance score is simple, the functions simply
// return the edit distance values.
SEQAN_DEFINE_TEST(test_score_edit) {
    // TODO(holtgrew): Break out each block into a test of its own.
    // We will only test with int scores, the most common case.
    typedef int TValue;
    typedef Score<TValue, EditDistance> TScore;
    
    // Test the default constructor.
    {
        TScore editDistanceScore;
        SEQAN_ASSERT_EQ(0, scoreMatch(editDistanceScore));
        SEQAN_ASSERT_EQ(-1, scoreMismatch(editDistanceScore));
        // TODO(holtgrew): FIXME! Does not compile because of reference to temporary.
        // SEQAN_ASSERT_EQ(-1, scoreGap(editDistanceScore));
        SEQAN_ASSERT_EQ(-1, scoreGapExtend(editDistanceScore));
        SEQAN_ASSERT_EQ(-1, scoreGapOpen(editDistanceScore));
        // Test function score().
        SEQAN_ASSERT_EQ(0, score(editDistanceScore, 'A', 'A'));
        SEQAN_ASSERT_EQ(-1, score(editDistanceScore, 'A', 'C'));
    }

    // Test the const accessor functions.
    {
        const TScore editDistanceScore;
        SEQAN_ASSERT_EQ(0, scoreMatch(editDistanceScore));
        SEQAN_ASSERT_EQ(-1, scoreMismatch(editDistanceScore));
        // TODO(holtgrew): FIXME! Does not compile because of reference to temporary.
        // SEQAN_ASSERT_EQ(-1, scoreGap(editDistanceScore));
        SEQAN_ASSERT_EQ(-1, scoreGapExtend(editDistanceScore));
        SEQAN_ASSERT_EQ(-1, scoreGapOpen(editDistanceScore));
        // Test function score().
        SEQAN_ASSERT_EQ(0, score(editDistanceScore, 'A', 'A'));
        SEQAN_ASSERT_EQ(-1, score(editDistanceScore, 'A', 'C'));
    }

    // Test the shortcut.  Since the assignment operator is not
    // overloaded, the following should not compile if the shortcut is
    // not defined appropriately.
    {
        // TODO(holtgrew): Use a metaprogramming for type equality instead?, something like SEQAN_ASSERT(type-equals(x, y));?
        Score<int, EditDistance> scoreEditDistance;
        EditDistanceScore editDistanceScore;
        scoreEditDistance = editDistanceScore;
    }
}


// Score<TValue, Simple> is, you have guessed it, very simple.  Thus,
// we will simply use different constructors to construct Simple Score
// objects and call the scoring functions on them.
SEQAN_DEFINE_TEST(test_score_simple) {
    // TODO(holtgrew): Break out each block into a test of its own.
    // We will only test with int scores, the most common case.
    typedef int TValue;
    typedef Score<TValue, Simple> TScore;
    
    // Test the default constructor.
    {
        TScore simpleScore;
        SEQAN_ASSERT_EQ(0, scoreMatch(simpleScore));
        SEQAN_ASSERT_EQ(-1, scoreMismatch(simpleScore));
        SEQAN_ASSERT_EQ(-1, scoreGap(simpleScore));
        SEQAN_ASSERT_EQ(-1, scoreGapExtend(simpleScore));
        SEQAN_ASSERT_EQ(-1, scoreGapOpen(simpleScore));
        // Test function score().
        SEQAN_ASSERT_EQ(0, score(simpleScore, 'A', 'A'));
        SEQAN_ASSERT_EQ(-1, score(simpleScore, 'A', 'C'));
    }

    // Test the const member retrieval functions.
    // TODO(holtgrew): Should these functions not be called getFUNCNAME?
    {
        const TScore simpleScore;
        SEQAN_ASSERT_EQ(0, scoreMatch(simpleScore));
        SEQAN_ASSERT_EQ(-1, scoreMismatch(simpleScore));
        SEQAN_ASSERT_EQ(-1, scoreGap(simpleScore));
        SEQAN_ASSERT_EQ(-1, scoreGapExtend(simpleScore));
        SEQAN_ASSERT_EQ(-1, scoreGapOpen(simpleScore));
        // Test function score().
        SEQAN_ASSERT_EQ(0, score(simpleScore, 'A', 'A'));
        SEQAN_ASSERT_EQ(-1, score(simpleScore, 'A', 'C'));
    }

    // Test the non-const member retrieval functions with assignments.
    {
        const int kMatch = 1;
        const int kMismatch = 2;
        const int kGapExtension = 3;
        const int kGapOpen = 4;
        // Perform assignments.
        TScore simpleScore;
        setScoreMatch(simpleScore, kMatch);
        setScoreMismatch(simpleScore, kMismatch);
        setScoreGap(simpleScore, kGapExtension);
        setScoreGapExtend(simpleScore, kGapExtension);
        setScoreGapOpen(simpleScore, kGapOpen);
        // Check results.
        SEQAN_ASSERT_EQ(kMatch, scoreMatch(simpleScore));
        SEQAN_ASSERT_EQ(kMismatch, scoreMismatch(simpleScore));
        SEQAN_ASSERT_EQ(kGapExtension, scoreGap(simpleScore));
        SEQAN_ASSERT_EQ(kGapExtension, scoreGapExtend(simpleScore));
        SEQAN_ASSERT_EQ(kGapOpen, scoreGapOpen(simpleScore));
    }
 
    // Test the constructor with match, mismatch, gap arguments.
    {
        // Define constant test data.
        const int kMatch = 1;
        const int kMismatch = 2;
        const int kGap = 3;
        // Construct the score and make assertions about its state.
        TScore simpleScore(kMatch, kMismatch, kGap);
        SEQAN_ASSERT_EQ(kMatch, scoreMatch(simpleScore));
        SEQAN_ASSERT_EQ(kMismatch, scoreMismatch(simpleScore));
        SEQAN_ASSERT_EQ(kGap, scoreGap(simpleScore));
        SEQAN_ASSERT_EQ(kGap, scoreGapExtend(simpleScore));
        SEQAN_ASSERT_EQ(kGap, scoreGapOpen(simpleScore));
        // Test function score().
        SEQAN_ASSERT_EQ(kMatch, score(simpleScore, 'A', 'A'));
        SEQAN_ASSERT_EQ(kMismatch, score(simpleScore, 'A', 'C'));
    }

    // Test the constructor with match, mismatch, gap extension, gap
    // open arguments.
    {
        // Define constant test data.
        const int kMatch = 1;
        const int kMismatch = 2;
        const int kGapExtension = 3;
        const int kGapOpen = 4;
        // Construct the score and make assertions about its state.
        TScore simpleScore(kMatch, kMismatch, kGapExtension, kGapOpen);
        SEQAN_ASSERT_EQ(kMatch, scoreMatch(simpleScore));
        SEQAN_ASSERT_EQ(kMismatch, scoreMismatch(simpleScore));
        SEQAN_ASSERT_EQ(kGapExtension, scoreGap(simpleScore));
        SEQAN_ASSERT_EQ(kGapExtension, scoreGapExtend(simpleScore));
        SEQAN_ASSERT_EQ(kGapOpen, scoreGapOpen(simpleScore));
        // Test function score().
        SEQAN_ASSERT_EQ(kMatch, score(simpleScore, 'A', 'A'));
        SEQAN_ASSERT_EQ(kMismatch, score(simpleScore, 'A', 'C'));
    }

    // Test the shortcut.  Since the assignment operator is not
    // overloaded, the following should not compile if the shortcut is
    // not defined appropriately.
    {
        // TODO(holtgrew): Use a metaprogramming for type equality instead?, something like SEQAN_ASSERT(type-equals(x, y));?
        Score<int, Simple> scoreSimple;
        SimpleScore simpleScore;
        simpleScore = scoreSimple;
    }
}


// Test the built-in data matrices by comparing them to the matrices
// from test data files.
SEQAN_DEFINE_TEST(test_score_matrix_data) {
    // We test the scoring matrix with amino acids as the underlying
    // sequence.  This should be the most common case.
    typedef int TValue;
    typedef Score<TValue, ScoreMatrix<AminoAcid, ScoreMatrixFile> > TScore;

    // TODO(holtgrew): It should be easier to construct these paths.
    String<char> pathToTestSrc = getAbsolutePath("/tests/score/");

    // Test with BLOSUM30.
    {
        // The built-in BLOSUM30 matrix.
        Blosum30 blosum30;
        // Use a quick-and-dirty test that the matrix indeed is BLOSUM30.
        SEQAN_ASSERT_EQ(0, score(blosum30, AminoAcid('A'), AminoAcid('N')));
        SEQAN_ASSERT_EQ(1, score(blosum30, AminoAcid('V'), AminoAcid('A')));

        // Build path to BLOSUM30 matrix.
        String<char> pathToBlosum30(pathToTestSrc);
        append(pathToBlosum30, "BLOSUM30");
        // Load matrix.
        TScore loadedBlosum30(toCString(pathToBlosum30));
        // Compare loaded with built-in matrix.
        assertAminoAcidMatricesAreEqual(loadedBlosum30, blosum30);
    }

    // Test with BLOSUM45.
    {
        // The built-in BLOSUM45 matrix.
        Blosum45 blosum45;
        // Use a quick-and-dirty test that the matrix indeed is BLOSUM45.
        SEQAN_ASSERT_EQ(-1, score(blosum45, AminoAcid('A'), AminoAcid('N')));
        SEQAN_ASSERT_EQ(0, score(blosum45, AminoAcid('V'), AminoAcid('A')));

        // Build path to BLOSUM45 matrix.
        String<char> pathToBlosum45(pathToTestSrc);
        append(pathToBlosum45, "BLOSUM45");
        // Load matrix.
        TScore loadedBlosum45(toCString(pathToBlosum45));
        // Compare loaded with built-in matrix.
        assertAminoAcidMatricesAreEqual(loadedBlosum45, blosum45);
    }

    // Test with BLOSUM62.
    {
        // The built-in BLOSUM62 matrix.
        Blosum62 blosum62;
        // Use a quick-and-dirty test that the matrix indeed is BLOSUM62.
        SEQAN_ASSERT_EQ(-2, score(blosum62, AminoAcid('A'), AminoAcid('N')));
        SEQAN_ASSERT_EQ(0, score(blosum62, AminoAcid('V'), AminoAcid('A')));

        // Build path to BLOSUM62 matrix.
        String<char> pathToBlosum62(pathToTestSrc);
        append(pathToBlosum62, "BLOSUM62");
        // Load matrix.
        TScore loadedBlosum62(toCString(pathToBlosum62));
        // Compare loaded with built-in matrix.
        assertAminoAcidMatricesAreEqual(loadedBlosum62, blosum62);
    }

    // Test with BLOSUM80.
    {
        // The built-in BLOSUM80 matrix.
        Blosum80 blosum80;
        // Use a quick-and-dirty test that the matrix indeed is BLOSUM80.
        SEQAN_ASSERT_EQ(-3, score(blosum80, AminoAcid('A'), AminoAcid('N')));
        SEQAN_ASSERT_EQ(-1, score(blosum80, AminoAcid('V'), AminoAcid('A')));

        // Build path to BLOSUM80 matrix.
        String<char> pathToBlosum80(pathToTestSrc);
        append(pathToBlosum80, "BLOSUM80");
        // Load matrix.
        TScore loadedBlosum80(toCString(pathToBlosum80));
        // Compare loaded with built-in matrix.
        assertAminoAcidMatricesAreEqual(loadedBlosum80, blosum80);
    }

    // Test with PAM40.
    {
        // The built-in PAM40 matrix.
        Pam40 pam40;
        // Use a quick-and-dirty test that the matrix indeed is PAM40.
        SEQAN_ASSERT_EQ(-3, score(pam40, AminoAcid('A'), AminoAcid('N')));
        SEQAN_ASSERT_EQ(-2, score(pam40, AminoAcid('V'), AminoAcid('A')));

        // Build path to PAM40 matrix.
        String<char> pathToPam40(pathToTestSrc);
        append(pathToPam40, "PAM40");
        // Load matrix.
        TScore loadedPam40(toCString(pathToPam40));
        // Compare loaded with built-in matrix.
        assertAminoAcidMatricesAreEqual(loadedPam40, pam40);
    }

    // Test with PAM120.
    {
        // The built-in PAM120 matrix.
        Pam120 pam120;
        // Use a quick-and-dirty test that the matrix indeed is PAM120.
        SEQAN_ASSERT_EQ(-1, score(pam120, AminoAcid('A'), AminoAcid('N')));
        SEQAN_ASSERT_EQ(0, score(pam120, AminoAcid('V'), AminoAcid('A')));

        // Build path to PAM120 matrix.
        String<char> pathToPam120(pathToTestSrc);
        append(pathToPam120, "PAM120");
        // Load matrix.
        TScore loadedPam120(toCString(pathToPam120));
        // Compare loaded with built-in matrix.
        assertAminoAcidMatricesAreEqual(loadedPam120, pam120);
    }

    // Test with PAM200.
    {
        // The built-in PAM200 matrix.
        Pam200 pam200;
        // Use a quick-and-dirty test that the matrix indeed is PAM200.
        SEQAN_ASSERT_EQ(0, score(pam200, AminoAcid('A'), AminoAcid('N')));
        SEQAN_ASSERT_EQ(0, score(pam200, AminoAcid('V'), AminoAcid('A')));

        // Build path to PAM200 matrix.
        String<char> pathToPam200(pathToTestSrc);
        append(pathToPam200, "PAM200");
        // Load matrix.
        TScore loadedPam200(toCString(pathToPam200));
        // Compare loaded with built-in matrix.
        assertAminoAcidMatricesAreEqual(loadedPam200, pam200);
    }

    // Test with PAM250.
    {
        // The built-in PAM250 matrix.
        Pam250 pam250;
        // Use a quick-and-dirty test that the matrix indeed is PAM250.
        SEQAN_ASSERT_EQ(0, score(pam250, AminoAcid('A'), AminoAcid('N')));
        SEQAN_ASSERT_EQ(0, score(pam250, AminoAcid('V'), AminoAcid('A')));

        // Build path to PAM250 matrix.
        String<char> pathToPam250(pathToTestSrc);
        append(pathToPam250, "PAM250");
        // Load matrix.
        TScore loadedPam250(toCString(pathToPam250));
        // Compare loaded with built-in matrix.
        assertAminoAcidMatricesAreEqual(loadedPam250, pam250);
    }

    // Test with VTML200.
    {
        // The built-in VTML200 matrix.
        Vtml200 vtml200;
        // Use a quick-and-dirty test that the matrix indeed is VTML200.
        SEQAN_ASSERT_EQ(-1, score(vtml200, AminoAcid('A'), AminoAcid('N')));
        SEQAN_ASSERT_EQ(0, score(vtml200, AminoAcid('V'), AminoAcid('A')));

        // Build path to VTML200 matrix.
        String<char> pathToVTML200(pathToTestSrc);
        append(pathToVTML200, "VTML200I");
        // Load matrix.
        TScore loadedVTML200(toCString(pathToVTML200));
        // Compare loaded with built-in matrix.
        assertAminoAcidMatricesAreEqual(loadedVTML200, vtml200);
    }
}

template <typename TScoreSpec>
void testScoreSequenceEntryForScore()
{
    typedef Score<int, TScoreSpec> TScoringScheme;
    typedef typename SequenceEntryForScore<TScoringScheme, DnaString>::Type TEntry1;

    DnaString seq1 = "ACGTACG";
    Dna5String seq2 = "ACGTNACNA";

    TScoringScheme scoringScheme;

    TEntry1 val1 = sequenceEntryForScore(scoringScheme, seq1, 0);
    TEntry1 val2 = sequenceEntryForScore(scoringScheme, seq2, 0);

    SEQAN_ASSERT_EQ(val1, 'A');
    SEQAN_ASSERT_EQ(val2, 'A');

    val1 = sequenceEntryForScore(scoringScheme, seq1, 4);
    val2 = sequenceEntryForScore(scoringScheme, seq2, 4);

    SEQAN_ASSERT_EQ(val1, 'A');
    SEQAN_ASSERT_EQ(val2, 'N');
}

SEQAN_DEFINE_TEST(test_score_sequence_entry_for_score)
{
    testScoreSequenceEntryForScore<Simple>();
    testScoreSequenceEntryForScore<ScoreMatrix<AminoAcid, Blosum62_> >();
}

SEQAN_DEFINE_TEST(test_score_dynamic_score_matrix)
{
    Blosum62 blo;
    SelectableAminoAcidMatrix sel;
    SEQAN_ASSERT(getScoreMatrixId(sel) != AminoAcidScoreMatrixID::BLOSUM62);

    setScoreMatrixById(sel, AminoAcidScoreMatrixID::BLOSUM62);
    SEQAN_ASSERT(getScoreMatrixId(sel) == AminoAcidScoreMatrixID::BLOSUM62);

    for (unsigned i = 0; i < ValueSize<AminoAcid>::VALUE; ++i)
        for (unsigned j = 0; j < ValueSize<AminoAcid>::VALUE; ++j)
            SEQAN_ASSERT_EQ(score(blo, static_cast<AminoAcid>(i), static_cast<AminoAcid>(j)),
                            score(sel, static_cast<AminoAcid>(i), static_cast<AminoAcid>(j)));
}

SEQAN_BEGIN_TESTSUITE(test_score) {
    // Call the tests for this module.
    SEQAN_CALL_TEST(test_score_gap_open);
    SEQAN_CALL_TEST(test_score_simple);
    SEQAN_CALL_TEST(test_score_edit);
    SEQAN_CALL_TEST(test_score_matrix);
    SEQAN_CALL_TEST(test_score_matrix_file);
    SEQAN_CALL_TEST(test_score_matrix_data);
    SEQAN_CALL_TEST(test_score_sequence_entry_for_score);
    SEQAN_CALL_TEST(test_score_dynamic_score_matrix);
}
SEQAN_END_TESTSUITE
