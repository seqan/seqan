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
// Tests for the classic DP programming algorithms in the seeds module.
// ==========================================================================

#ifndef TEST_SEEDS_TEST_ALIGN_DYNPROG_LINEAR_H_
#define TEST_SEEDS_TEST_ALIGN_DYNPROG_LINEAR_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

// Test matrix resizing.
SEQAN_DEFINE_TEST(test_align_dynprog_linear_resize_matrix)
{
    using namespace seqan;

    Matrix<int, 2> matrix;

    _alignResizeMatrix(matrix, CharString("length 9"), CharString("length    11"), NeedlemanWunsch());

    SEQAN_ASSERT_EQ(9u, length(matrix, 0));
    SEQAN_ASSERT_EQ(13u, length(matrix, 1));
}


// Test gutter initialization if gap costs are free.
SEQAN_DEFINE_TEST(test_align_dynprog_linear_init_gutter_free)
{
    using namespace seqan;

    Matrix<int, 2> matrix;
    setLength(matrix, 0, 2);
    setLength(matrix, 1, 3);
    resize(matrix);

    _alignInitGutter(matrix, Score<int, Simple>(1, -1, -2), AlignConfig<true, true, true, true>(), NeedlemanWunsch());

    SEQAN_ASSERT_EQ(0, value(matrix, 0, 0));
    SEQAN_ASSERT_EQ(0, value(matrix, 0, 1));
    SEQAN_ASSERT_EQ(0, value(matrix, 1, 0));
    SEQAN_ASSERT_EQ(0, value(matrix, 2, 0));
}


// Test gutter initialization if gap costs are not free.
SEQAN_DEFINE_TEST(test_align_dynprog_linear_init_gutter_not_free)
{
    using namespace seqan;

    Matrix<int, 2> matrix;
    setLength(matrix, 0, 2);
    setLength(matrix, 1, 3);
    resize(matrix);

    _alignInitGutter(matrix, Score<int, Simple>(1, -1, -2), AlignConfig<false, false, true, true>(), NeedlemanWunsch());

    SEQAN_ASSERT_EQ(0, value(matrix, 0, 0));
    SEQAN_ASSERT_EQ(-2, value(matrix, 0, 1));
    SEQAN_ASSERT_EQ(-4, value(matrix, 0, 2));
    SEQAN_ASSERT_EQ(-2, value(matrix, 1, 0));
}


// Test DP matrix filling
SEQAN_DEFINE_TEST(test_align_dynprog_linear_fill_matrix)
{
    using namespace seqan;

    Matrix<int, 2> matrix;
    DnaString const sequence0 = "CCA";
    DnaString const sequence1 = "CAA";
    Score<int, Simple> const scoringScheme(1, -1, -1);

    _alignResizeMatrix(matrix, sequence0, sequence1, NeedlemanWunsch());
    _alignInitGutter(matrix, scoringScheme, AlignConfig<false, false, false, false>(), NeedlemanWunsch());
    _alignFillMatrix(matrix, sequence0, sequence1, scoringScheme, NeedlemanWunsch());

    int const expected[16] = {
         0, -1, -2, -3,
        -1,  1,  0, -1,
        -2,  0,  0, -1,
        -3, -1,  1,  1
    };

    for (unsigned j = 0; j < 4; ++j) {
        for (unsigned i = 0; i < 4; ++i) {
            SEQAN_ASSERT_EQ_MSG(expected[i * 4 + j], value(matrix, i, j), "i = %d, j = %d", i, j);
        }
    }
}


SEQAN_DEFINE_TEST(test_align_dynprog_linear_traceback)
{
    using namespace seqan;

    typedef CharString TString;
    typedef Position<CharString>::Type TPosition;
    typedef Align<TString> TAlign;
    typedef Row<TAlign>::Type TAlignRow;
    typedef Iterator<TAlignRow, Standard>::Type TAlignRowIterator;
    typedef Iterator<TString, Standard>::Type TStringIterator;

    // Case: No free begin/end gaps.
    {
        // Fill the matrix with DP (tested by other tests).
        Matrix<int, 2> matrix;
        TString const sequence0 = "CCAAA";
        TString const sequence1 = "CAA";
        Score<int, Simple> const scoringScheme(1, -1, -1);
        
        _alignResizeMatrix(matrix, sequence0, sequence1, NeedlemanWunsch());
        _alignInitGutter(matrix, scoringScheme, AlignConfig<false, false, false, false>(), NeedlemanWunsch());
        _alignFillMatrix(matrix, sequence0, sequence1, scoringScheme, NeedlemanWunsch());

        // Perform the traceback.
        Align<TString> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequence0);
        assignSource(row(alignment, 1), sequence1);

        size_t finalPos0 = 0;
        size_t finalPos1 = 0;
        TStringIterator seq0It = end(sequence0) - 1;
        TStringIterator seq1It = end(sequence1) - 1;
        TAlignRowIterator align0It = end(row(alignment, 0));
        TAlignRowIterator align1It = end(row(alignment, 1));
        int score = _alignTraceback(align0It, align1It, seq0It, seq1It, finalPos0, finalPos1, matrix, scoringScheme, 0, 0, 0, 0, true, AlignConfig<false, false, false, false>(), NeedlemanWunsch());
        // std::cout << alignment;
        // // TODO(holtgrew): Debug code, remove when working.
        // {
        //     for (int k = 0; k < 1; ++k) {
        //         std::cout << ",-- *** filled alignment matrix " << k << std::endl;
        //         for (unsigned i = 0; i < length(matrix, 0); ++i) {
        //             std::cout << "| ";
        //             for (unsigned j = 0; j < length(matrix, 1); ++j) {
        //                 if (value(matrix, i, j, k) <= MinValue<int>::VALUE / 4)
        //                     std::cout << "\tinf";
        //                 else
        //                     std::cout << "\t" << value(matrix, i, j, k);
        //             }
        //             std::cout << std::endl;
        //         }
        //         std::cout << "`--" << std::endl;
        //     }
        // }

        SEQAN_ASSERT_EQ(score, 1);
        SEQAN_ASSERT(seq0It + 1 == begin(sequence0));
        SEQAN_ASSERT(seq1It + 1 == begin(sequence1));
        // TODO(holtgrew): Why does this not work?
        // SEQAN_ASSERT(align0It == begin(row(alignment, 0)));
        // SEQAN_ASSERT(align1It == begin(row(alignment, 1)));
        SEQAN_ASSERT_EQ(finalPos0, 0u);
        SEQAN_ASSERT_EQ(finalPos1, 0u);
        // Expected alignment:
        //
        //  CCAAA
        //  C-AA-
        SEQAN_ASSERT(row(alignment, 0) == "CCAAA");
        SEQAN_ASSERT(row(alignment, 1) == "C-AA"/*-*/);
        // Leading gaps are not shown, we test through the iterators.
        SEQAN_ASSERT(isGap(iter(row(alignment, 1), 0)));
        SEQAN_ASSERT_NOT(isGap(iter(row(alignment, 1), 1)));
        SEQAN_ASSERT(isGap(iter(row(alignment, 1), 2)));
        SEQAN_ASSERT_NOT(isGap(iter(row(alignment, 1), 3)));
        SEQAN_ASSERT_NOT(isGap(iter(row(alignment, 1), 4)));
    }
    // TODO(holtgrew): Case with free begin and end gaps.
}

#endif  // TEST_SEEDS_TEST_ALIGN_DYNPROG_LINEAR_H_
