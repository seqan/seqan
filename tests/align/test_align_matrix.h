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

#ifndef TESTS_ALIGN_TEST_ALIGN_MATRIX_H_
#define TESTS_ALIGN_TEST_ALIGN_MATRIX_H_

SEQAN_DEFINE_TEST(test_align_matrix)
{
    using namespace std;
    using namespace seqan;

    // Resize the matrix.
    Matrix<double,2> matrix1;

    setLength(matrix1, 0, 2);
    setLength(matrix1, 1, 2);
    resize(matrix1);
    value(matrix1,0,0) = 3.0;
    value(matrix1,0,1) = 3.5;
    value(matrix1,1,0) = 14.5;
    value(matrix1,1,1) = -6.0;

    Matrix<double,2> matrix2;

    setLength(matrix2, 0, 2);
    setLength(matrix2, 1, 2);
    resize(matrix2,1.0);

    // Basic 2D operations, A+B,A-B,A*a,A*B, A==B

    Matrix<double,2> matrix3;

    matrix3 = matrix1 + matrix2;
    SEQAN_ASSERT_EQ(value(matrix3,0,0), value(matrix1,0,0) + value(matrix2,0,0));
    SEQAN_ASSERT_EQ(value(matrix3,0,1), value(matrix1,0,1) + value(matrix2,0,1));
    SEQAN_ASSERT_EQ(value(matrix3,1,0), value(matrix1,1,0) + value(matrix2,1,0));
    SEQAN_ASSERT_EQ(value(matrix3,1,1), value(matrix1,1,1) + value(matrix2,1,1));
    SEQAN_ASSERT_EQ(matrix1 + matrix2, matrix2 + matrix1);
    matrix3 = matrix1 - matrix2;
    SEQAN_ASSERT_EQ(value(matrix3,0,0), value(matrix1,0,0) - value(matrix2,0,0));
    SEQAN_ASSERT_EQ(value(matrix3,0,1), value(matrix1,0,1) - value(matrix2,0,1));
    SEQAN_ASSERT_EQ(value(matrix3,1,0), value(matrix1,1,0) - value(matrix2,1,0));
    SEQAN_ASSERT_EQ(value(matrix3,1,1), value(matrix1,1,1) - value(matrix2,1,1));

    SEQAN_ASSERT_EQ(matrix1-matrix2, matrix1+(matrix2*(-1.0)));
    matrix3=matrix1*matrix2;
    matrix3=matrix1*5.0;
    SEQAN_ASSERT_EQ(matrix1 * 5.0, 5.0 * matrix1);

    // n-dimensional matrix
    Matrix<double> matrixN;
    setDimension(matrixN,2);
    setLength(matrix2, 0, 3);
    setLength(matrix2, 1, 2);
    resize(matrixN,1.0);
}

#endif  // TESTS_ALIGN_TEST_ALIGN_MATRIX_H_
