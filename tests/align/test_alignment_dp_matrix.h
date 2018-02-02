// ==========================================================================
//                         test_alignment_dp_matrix.h
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_MATRIX_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_MATRIX_H_

#include <seqan/basic.h>

#include <seqan/align.h>

void testAlignmentDPMatrixDataHostMF()
{
    using namespace seqan;

    typedef DPMatrix_<char, FullDPMatrix> TDPMatrix;
    typedef DPMatrix_<char, FullDPMatrix> const TDPMatrixConst;

    typedef Member<TDPMatrix, DPMatrixMember>::Type TDataHost;
    typedef Member<TDPMatrixConst, DPMatrixMember>::Type TDataHostConst;

    bool result1 = IsSameType<Matrix<char, 2>, TDataHost>::VALUE;
    bool result2 = IsSameType<Matrix<char, 2> const, TDataHostConst>::VALUE;

    SEQAN_ASSERT_EQ(result1, true);
    SEQAN_ASSERT_EQ(result2, true);
}

void testAlignmentDPMatrixSizeArrMF()
{
    using namespace seqan;

    typedef DPMatrix_<char, FullDPMatrix> TDPMatrix;
    typedef DPMatrix_<char, FullDPMatrix> const TDPMatrixConst;

    typedef SizeArr_<TDPMatrix>::Type TSizeArr;
    typedef SizeArr_<TDPMatrixConst>::Type TSizeArrConst;

    bool result1 = IsSameType<SizeArr_<Matrix<char, 2> >::Type, TSizeArr>::VALUE;
    bool result2 = IsSameType<SizeArr_<Matrix<char, 2> >::Type const, TSizeArrConst>::VALUE;

    SEQAN_ASSERT_EQ(result1, true);
    SEQAN_ASSERT_EQ(result2, true);
}

void testAlignmentDPMatrixDataHost()
{
    using namespace seqan;

    typedef DPMatrix_<char, FullDPMatrix> TDPMatrix;
    typedef DPMatrix_<char, FullDPMatrix> const TDPMatrixConst;

    typedef Member<TDPMatrix, DPMatrixMember>::Type TDataHost;
    typedef Member<TDPMatrixConst, DPMatrixMember>::Type TDataHostConst;

    TDPMatrix dpMatrix;

    value(dpMatrix.data_host).data_lengths[0] = 10;
    value(dpMatrix.data_host).data_lengths[1] = 20;

    TDPMatrixConst dpMatrixConst(dpMatrix);

    TDataHost dataHost = value(dpMatrix.data_host);
    TDataHostConst dataHostConst = value(dpMatrixConst.data_host);

    SEQAN_ASSERT_EQ(dataHost.data_lengths[0], 10u);
    SEQAN_ASSERT_EQ(dataHost.data_lengths[1], 20u);
    SEQAN_ASSERT_EQ(dataHostConst.data_lengths[0], 10u);
    SEQAN_ASSERT_EQ(dataHostConst.data_lengths[1], 20u);
}

void testAlignmentDPMatrixDataLengths()
{
    using namespace seqan;

    typedef DPMatrix_<char, FullDPMatrix> TDPMatrix;
    typedef DPMatrix_<char, FullDPMatrix> const TDPMatrixConst;

    TDPMatrix dpMatrix;

    value(dpMatrix.data_host).data_lengths[0] = 10;
    value(dpMatrix.data_host).data_lengths[1] = 20;

    TDPMatrixConst dpMatrixConst(dpMatrix);

    SEQAN_ASSERT_EQ(_dataLengths(dpMatrix)[0], 10u);
    SEQAN_ASSERT_EQ(_dataLengths(dpMatrix)[1], 20u);
    SEQAN_ASSERT_EQ(_dataLengths(dpMatrixConst)[0], 10u);
    SEQAN_ASSERT_EQ(_dataLengths(dpMatrixConst)[1], 20u);
}

void testAlignmentDPMatrixDataFactors()
{
    using namespace seqan;

    typedef DPMatrix_<char, FullDPMatrix> TDPMatrix;
    typedef DPMatrix_<char, FullDPMatrix> const TDPMatrixConst;

    TDPMatrix dpMatrix;

    value(dpMatrix.data_host).data_factors[0] = 10;
    value(dpMatrix.data_host).data_factors[1] = 20;

    TDPMatrixConst dpMatrixConst(dpMatrix);

    SEQAN_ASSERT_EQ(_dataFactors(dpMatrix)[0], 10u);
    SEQAN_ASSERT_EQ(_dataFactors(dpMatrix)[1], 20u);
    SEQAN_ASSERT_EQ(_dataFactors(dpMatrixConst)[0], 10u);
    SEQAN_ASSERT_EQ(_dataFactors(dpMatrixConst)[1], 20u);
}

void testAlignmentDPMatrixCheckDimension()
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(_checkCorrectDimension(0u), true);
    SEQAN_ASSERT_EQ(_checkCorrectDimension(1u), true);
    SEQAN_ASSERT_EQ(_checkCorrectDimension(2u), false);
    SEQAN_ASSERT_EQ(_checkCorrectDimension(10u), false);
}

template <typename TSpec>
void testAlignmentDPMatrixSetLength(TSpec const &)
{
    using namespace seqan;

    DPMatrix_<char, TSpec> dpMatrix;

    SEQAN_ASSERT_EQ(_dataLengths(dpMatrix)[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(_dataLengths(dpMatrix)[DPMatrixDimension_::VERTICAL], 0u);

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);

    SEQAN_ASSERT_EQ(_dataLengths(dpMatrix)[DPMatrixDimension_::HORIZONTAL], 4u);
    SEQAN_ASSERT_EQ(_dataLengths(dpMatrix)[DPMatrixDimension_::VERTICAL], 3u);
}

template <typename TSpec>
void testAlignmentDPMatrixLengthDimension(TSpec const &)
{
    using namespace seqan;

    DPMatrix_<char, TSpec> dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);

    SEQAN_ASSERT_EQ(length(dpMatrix, +DPMatrixDimension_::HORIZONTAL), 4u);
    SEQAN_ASSERT_EQ(length(dpMatrix, +DPMatrixDimension_::VERTICAL), 3u);
}

template <typename TSpec>
void testAlignmentDPMatrixEmpty(TSpec const &)
{
    using namespace seqan;

    DPMatrix_<char, TSpec> dpMatrix;

    SEQAN_ASSERT_EQ(empty(dpMatrix), true);

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);

    resize(dpMatrix, 'x');

    SEQAN_ASSERT_EQ(empty(dpMatrix), false);
}

void testAlignmentDPMatrixClear()
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;


    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);
    resize(dpMatrix);

    String<char> str = "Hello World!";
    setValue(value(dpMatrix.data_host).data_host, str);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrix.data_host).data_host), true);
    SEQAN_ASSERT_EQ(empty(dpMatrix), false);
    SEQAN_ASSERT_EQ(host(value(dpMatrix.data_host)), "Hello World!");
    clear(dpMatrix);

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrix.data_host).data_host), false);
    SEQAN_ASSERT_EQ(empty(dpMatrix), true);
}

void testAlignmentDPMatrixHost()
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);

    resize(dpMatrix, 'x');

    SEQAN_ASSERT_EQ(host(dpMatrix), "xxxxxxxxxxxx");
}

void testAlignmentDPMatrixSetHost()
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    String<char> str = "xxxxxxxxxxxx";
    setHost(dpMatrix, str);

    SEQAN_ASSERT_EQ(host(dpMatrix), "xxxxxxxxxxxx");
}

template <typename TIterSpec>
void testAlignmentDPMatrixBegin(TIterSpec const &)
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);
    resize(dpMatrix, 'x');

    value(dpMatrix, 0) = 'a';

    SEQAN_ASSERT_EQ(host(dpMatrix), "axxxxxxxxxxx");
    SEQAN_ASSERT_EQ(value(begin(dpMatrix, TIterSpec())), 'a');
}

template <typename TIterSpec>
void testAlignmentDPMatrixEnd(TIterSpec const &)
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);
    resize(dpMatrix, 'x');

    value(dpMatrix, length(dpMatrix) - 1) = 'a';

    SEQAN_ASSERT_EQ(host(dpMatrix), "xxxxxxxxxxxa");
    SEQAN_ASSERT_EQ(value(end(dpMatrix, TIterSpec()) - 1), 'a');
}

template <typename TSpec>
void testAlignmentDPMatrixGetValueMF(TSpec const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, TSpec> TDPMatrix;
    typedef DPMatrix_<int, TSpec> const TDPMatrixConst;

    bool result = IsSameType<typename GetValue<TDPMatrix>::Type, int const &>::VALUE;
    bool result2 = IsSameType<typename GetValue<TDPMatrixConst>::Type, int const &>::VALUE;

    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(result2, true);
}

template <typename TSpec>
void testAlignmentDPMatrixValueMF(TSpec const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, TSpec> TDPMatrix;
    typedef DPMatrix_<int, TSpec> const TDPMatrixConst;

    bool result = IsSameType<typename Value<TDPMatrix>::Type, int>::VALUE;
    bool result2 = IsSameType<typename Value<TDPMatrixConst>::Type, int const>::VALUE;

    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(result2, true);
}

template <typename TSpec>
void testAlignmentDPMatrixReferenceMF(TSpec const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, TSpec> TDPMatrix;
    typedef DPMatrix_<int, TSpec> const TDPMatrixConst;

    bool result = IsSameType<typename Reference<TDPMatrix>::Type, int &>::VALUE;
    bool result2 = IsSameType<typename Reference<TDPMatrixConst>::Type, int const &>::VALUE;

    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(result2, true);
}

template <typename TSpec>
void testAlignmentDPMatrixPositionMF(TSpec const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, TSpec> TDPMatrix;
    typedef DPMatrix_<int, TSpec> const TDPMatrixConst;

    typedef typename Position<Matrix<int, 2> >::Type TMatrixPos;

    bool result = IsSameType<typename Position<TDPMatrix>::Type, TMatrixPos>::VALUE;
    bool result2 = IsSameType<typename Position<TDPMatrixConst>::Type, TMatrixPos>::VALUE;

    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(result2, true);
}

template <typename TSpec>
void testAlignmentDPMatrixSizeMF(TSpec const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, TSpec> TDPMatrix;
    typedef DPMatrix_<int, TSpec> const TDPMatrixConst;

    typedef typename Size<Matrix<int, 2> >::Type TMatrixSize;

    bool result = IsSameType<typename Size<TDPMatrix>::Type, TMatrixSize>::VALUE;
    bool result2 = IsSameType<typename Size<TDPMatrixConst>::Type, TMatrixSize>::VALUE;

    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(result2, true);
}

template <typename TSpec>
void testAlignmentDPMatrixHostMF(TSpec const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, TSpec> TDPMatrix;
    typedef DPMatrix_<int, TSpec> const TDPMatrixConst;

    typedef String<int> THost;
    typedef String<int> const THostConst;

    bool result = IsSameType<typename Host<TDPMatrix>::Type, THost>::VALUE;
    bool result2 = IsSameType<typename Host<TDPMatrixConst>::Type, THostConst>::VALUE;

    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(result2, true);
}

template <typename TSpec>
void testAlignmentDPMatrixIteratorMF(TSpec const &, seqan::Standard const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, TSpec> TDPMatrix;
    typedef DPMatrix_<int, TSpec> const TDPMatrixConst;

    typedef String<int> THost;
    typedef String<int> const THostConst;

    typedef typename Iterator<THost>::Type TIterator;
    typedef typename Iterator<THostConst>::Type TIteratorConst;

    bool result = IsSameType<typename Iterator<TDPMatrix, Standard>::Type, TIterator>::VALUE;
    bool result2 = IsSameType<typename Iterator<TDPMatrixConst, Standard>::Type, TIteratorConst>::VALUE;
    bool result3 = IsSameType<typename Iterator<TDPMatrix>::Type, TIterator>::VALUE;
    bool result4 = IsSameType<typename Iterator<TDPMatrixConst>::Type, TIteratorConst>::VALUE;

    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(result2, true);
    SEQAN_ASSERT_EQ(result3, true);
    SEQAN_ASSERT_EQ(result4, true);
}

template <typename TSpec>
void testAlignmentDPMatrixIteratorMF(TSpec const &, seqan::Rooted const &)
{
    using namespace seqan;

    typedef DPMatrix_<int, TSpec> TDPMatrix;
    typedef DPMatrix_<int, TSpec> const TDPMatrixConst;

    typedef Matrix<int, 2> THost;
    typedef Matrix<int, 2> const THostConst;

    typedef typename Iterator<THost>::Type TIterator;
    typedef typename Iterator<THostConst>::Type TIteratorConst;

    bool result = IsSameType<typename Iterator<TDPMatrix, Rooted>::Type, TIterator>::VALUE;
    bool result2 = IsSameType<typename Iterator<TDPMatrixConst, Rooted>::Type, TIteratorConst>::VALUE;
    bool result3 = IsSameType<typename Iterator<TDPMatrix>::Type, TIterator>::VALUE;
    bool result4 = IsSameType<typename Iterator<TDPMatrixConst>::Type, TIteratorConst>::VALUE;

    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(result2, true);
    SEQAN_ASSERT_EQ(result3, false);
    SEQAN_ASSERT_EQ(result4, false);
}

// ----------------------------------------------------------------------------
// Test the constructors of the dp matrix.
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_constructor)
{
    using namespace seqan;

    DPMatrix_<int, FullDPMatrix> dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrix.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrix.data_host).data_host), false);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_copy_constructor)
{
    using namespace seqan;

    DPMatrix_<int, FullDPMatrix> dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrix.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrix.data_host).data_host), false);

    value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL] = 10;
    value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL] = 20;
    value(dpMatrix.data_host).data_factors[DPMatrixDimension_::HORIZONTAL] = 20;

    DPMatrix_<int, FullDPMatrix> dpMatrixCopy1 = dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 10u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 20u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 20u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrixCopy1.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrixCopy1.data_host).data_host), false);

    DPMatrix_<int, FullDPMatrix> const dpMatrixCopy2 = dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrixCopy2.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 10u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy2.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 20u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy2.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 20u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy2.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrixCopy2.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrixCopy2.data_host).data_host), false);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_assigment)
{
    using namespace seqan;

    DPMatrix_<int, FullDPMatrix> dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrix.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrix.data_host).data_host), false);

    value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL] = 10;
    value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL] = 20;
    value(dpMatrix.data_host).data_factors[DPMatrixDimension_::HORIZONTAL] = 20;

    DPMatrix_<int, FullDPMatrix> dpMatrixCopy1;
    dpMatrixCopy1 = dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 10u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 20u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 20u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrixCopy1.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrixCopy1.data_host).data_host), false);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_constructor)
{
    using namespace seqan;

    DPMatrix_<int, SparseDPMatrix> dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrix.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrix.data_host).data_host), false);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_copy_constructor)
{
    using namespace seqan;

    DPMatrix_<int, SparseDPMatrix> dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrix.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrix.data_host).data_host), false);

    value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL] = 10;
    value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL] = 20;

    DPMatrix_<int, SparseDPMatrix> dpMatrixCopy1 = dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 10u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 20u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrixCopy1.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrixCopy1.data_host).data_host), false);

    DPMatrix_<int, SparseDPMatrix> const dpMatrixCopy2 = dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrixCopy2.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 10u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy2.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 20u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy2.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy2.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrixCopy2.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrixCopy2.data_host).data_host), false);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_assigment)
{
    using namespace seqan;

    DPMatrix_<int, SparseDPMatrix> dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrix.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrix.data_host).data_host), false);

    value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL] = 10;
    value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL] = 20;

    DPMatrix_<int, SparseDPMatrix> dpMatrixCopy1;
    dpMatrixCopy1 = dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 10u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 20u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_factors[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrixCopy1.data_host).data_factors[DPMatrixDimension_::VERTICAL], 1u);
    SEQAN_ASSERT_EQ(empty(value(dpMatrixCopy1.data_host).data_host), false);
    SEQAN_ASSERT_EQ(dependent(value(dpMatrixCopy1.data_host).data_host), false);
}

// ----------------------------------------------------------------------------
// Test the meta-function of the dp matrix.
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_metafunction_data_host)
{
    testAlignmentDPMatrixDataHostMF();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_metafunction_size_arr)
{
    testAlignmentDPMatrixSizeArrMF();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_metafunction_value)
{
    testAlignmentDPMatrixValueMF(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_metafunction_value)
{
    testAlignmentDPMatrixValueMF(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_metafunction_reference)
{
    testAlignmentDPMatrixReferenceMF(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_metafunction_reference)
{
    testAlignmentDPMatrixReferenceMF(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_metafunction_getvalue)
{
    testAlignmentDPMatrixGetValueMF(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_metafunction_getvalue)
{
    testAlignmentDPMatrixGetValueMF(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_metafunction_position)
{
    testAlignmentDPMatrixPositionMF(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_metafunction_position)
{
    testAlignmentDPMatrixPositionMF(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_metafunction_size)
{
    testAlignmentDPMatrixSizeMF(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_metafunction_size)
{
    testAlignmentDPMatrixSizeMF(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_metafunction_host)
{
    testAlignmentDPMatrixHostMF(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_metafunction_host)
{
    testAlignmentDPMatrixHostMF(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_metafunction_iterator_standard)
{
    testAlignmentDPMatrixIteratorMF(seqan::FullDPMatrix(), seqan::Standard());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_metafunction_iterator_rooted)
{
    testAlignmentDPMatrixIteratorMF(seqan::FullDPMatrix(), seqan::Rooted());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_metafunction_iterator_standard)
{
    testAlignmentDPMatrixIteratorMF(seqan::SparseDPMatrix(), seqan::Standard());
}


SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_metafunction_iterator_rooted)
{
    testAlignmentDPMatrixIteratorMF(seqan::SparseDPMatrix(), seqan::Rooted());
}

// ----------------------------------------------------------------------------
// Test the global interfaces of the dp matrix.
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_data_host)
{
    testAlignmentDPMatrixDataHost();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_data_lengths)
{
    testAlignmentDPMatrixDataLengths();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_data_factors)
{
    testAlignmentDPMatrixDataFactors();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_resize)
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    setLength(dpMatrix, +DPMatrixDimension_::HORIZONTAL, 3);
    setLength(dpMatrix, +DPMatrixDimension_::VERTICAL, 4);

    resize(dpMatrix);

    SEQAN_ASSERT_EQ(length(host(dpMatrix)), 12u);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_resize_with_value)
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    setLength(dpMatrix, +DPMatrixDimension_::HORIZONTAL, 3);
    setLength(dpMatrix, +DPMatrixDimension_::VERTICAL, 4);

    resize(dpMatrix, 'x');

    SEQAN_ASSERT_EQ(length(host(dpMatrix)), 12u);
    SEQAN_ASSERT_EQ(host(dpMatrix), "xxxxxxxxxxxx");
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_resize)
{
    using namespace seqan;

    DPMatrix_<char, SparseDPMatrix> dpMatrix;

    setLength(dpMatrix, +DPMatrixDimension_::HORIZONTAL, 3);
    setLength(dpMatrix, +DPMatrixDimension_::VERTICAL, 4);

    resize(dpMatrix);

    SEQAN_ASSERT_EQ(length(host(dpMatrix)), 4u);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_resize_with_value)
{
    using namespace seqan;

    DPMatrix_<char, SparseDPMatrix> dpMatrix;

    setLength(dpMatrix, +DPMatrixDimension_::HORIZONTAL, 3);
    setLength(dpMatrix, +DPMatrixDimension_::VERTICAL, 4);

    resize(dpMatrix, 'x');

    SEQAN_ASSERT_EQ(length(host(dpMatrix)), 4u);
    SEQAN_ASSERT_EQ(host(dpMatrix), "xxxx");
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_value)
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    setLength(dpMatrix, +DPMatrixDimension_::HORIZONTAL, 3);
    setLength(dpMatrix, +DPMatrixDimension_::VERTICAL, 4);

    resize(dpMatrix, 'x');

    SEQAN_ASSERT_EQ(value(dpMatrix, 0), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 4), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 8), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 11), 'x');

    value(dpMatrix, 1) = 'a';
    SEQAN_ASSERT_EQ(value(dpMatrix, 1), 'a');
    SEQAN_ASSERT_EQ(value(value(dpMatrix.data_host).data_host), "xaxxxxxxxxxx");

    DPMatrix_<char, FullDPMatrix> const dpMatrixConst(dpMatrix);

    SEQAN_ASSERT_EQ(value(dpMatrixConst, 1), 'a');
    SEQAN_ASSERT_EQ(value(value(dpMatrixConst.data_host).data_host), "xaxxxxxxxxxx");
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_value)
{
    using namespace seqan;

    DPMatrix_<char, SparseDPMatrix> dpMatrix;

    setLength(dpMatrix, +DPMatrixDimension_::HORIZONTAL, 3);
    setLength(dpMatrix, +DPMatrixDimension_::VERTICAL, 4);

    resize(dpMatrix, 'x');

    SEQAN_ASSERT_EQ(value(dpMatrix, 0), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 1), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 2), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 3), 'x');

    value(dpMatrix, 1) = 'a';
    SEQAN_ASSERT_EQ(value(dpMatrix, 1), 'a');
    SEQAN_ASSERT_EQ(value(value(dpMatrix.data_host).data_host), "xaxx");

    DPMatrix_<char, SparseDPMatrix> const dpMatrixConst(dpMatrix);

    SEQAN_ASSERT_EQ(value(dpMatrixConst, 1), 'a');
    SEQAN_ASSERT_EQ(value(value(dpMatrixConst.data_host).data_host), "xaxx");
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_value_with_coordinates)
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    setLength(dpMatrix, +DPMatrixDimension_::HORIZONTAL, 3);
    setLength(dpMatrix, +DPMatrixDimension_::VERTICAL, 4);

    resize(dpMatrix, 'x');

    SEQAN_ASSERT_EQ(value(dpMatrix, 0, 0), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 3, 2), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 1, 1), 'x');

    value(dpMatrix, 1, 0) = 'a';
    value(dpMatrix, 3, 2) = 'o';
    SEQAN_ASSERT_EQ(value(dpMatrix, 1, 0), 'a');
    SEQAN_ASSERT_EQ(value(dpMatrix, 3, 2), 'o');
    SEQAN_ASSERT_EQ(value(value(dpMatrix.data_host).data_host), "xaxxxxxxxxxo");

    DPMatrix_<char, FullDPMatrix> const dpMatrixConst(dpMatrix);

    SEQAN_ASSERT_EQ(value(dpMatrixConst, 1, 0), 'a');
    SEQAN_ASSERT_EQ(value(dpMatrixConst, 3, 2), 'o');
    SEQAN_ASSERT_EQ(value(value(dpMatrixConst.data_host).data_host), "xaxxxxxxxxxo");
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_value_with_coordinates)
{
    using namespace seqan;

    DPMatrix_<char, SparseDPMatrix> dpMatrix;

    setLength(dpMatrix, +DPMatrixDimension_::HORIZONTAL, 3);
    setLength(dpMatrix, +DPMatrixDimension_::VERTICAL, 4);

    resize(dpMatrix, 'x');

    SEQAN_ASSERT_EQ(value(dpMatrix, 0, 0), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 3, 2), 'x');
    SEQAN_ASSERT_EQ(value(dpMatrix, 1, 3), 'x');

    value(dpMatrix, 1, 0) = 'a';
    value(dpMatrix, 3, 2) = 'o';
    SEQAN_ASSERT_EQ(value(dpMatrix, 1, 0), 'a');
    SEQAN_ASSERT_EQ(value(dpMatrix, 3, 2), 'o');
    SEQAN_ASSERT_EQ(value(value(dpMatrix.data_host).data_host), "xaxo");

    DPMatrix_<char, SparseDPMatrix> const dpMatrixConst(dpMatrix);

    SEQAN_ASSERT_EQ(value(dpMatrixConst, 1, 0), 'a');
    SEQAN_ASSERT_EQ(value(dpMatrixConst, 3, 2), 'o');
    SEQAN_ASSERT_EQ(value(value(dpMatrixConst.data_host).data_host), "xaxo");
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_check_dimension)
{
    testAlignmentDPMatrixCheckDimension();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_set_length)
{
    testAlignmentDPMatrixSetLength(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_set_length)
{
    testAlignmentDPMatrixSetLength(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_length_dimension)
{
    testAlignmentDPMatrixLengthDimension(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_length_dimension)
{
    testAlignmentDPMatrixLengthDimension(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_length)
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);

    SEQAN_ASSERT_EQ(length(dpMatrix), 0u);
    resize(dpMatrix);
    SEQAN_ASSERT_EQ(length(dpMatrix), 12u);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_length)
{
    using namespace seqan;

    DPMatrix_<char, SparseDPMatrix> dpMatrix;

    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::HORIZONTAL], 0u);
    SEQAN_ASSERT_EQ(value(dpMatrix.data_host).data_lengths[DPMatrixDimension_::VERTICAL], 0u);

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);

    SEQAN_ASSERT_EQ(length(dpMatrix), 0u);
    resize(dpMatrix);
    SEQAN_ASSERT_EQ(length(dpMatrix), 3u);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_empty)
{
    testAlignmentDPMatrixEmpty(seqan::FullDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_empty)
{
    testAlignmentDPMatrixEmpty(seqan::SparseDPMatrix());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_clear)
{
    testAlignmentDPMatrixClear();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_host)
{
    testAlignmentDPMatrixHost();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_set_host)
{
    testAlignmentDPMatrixSetHost();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_full_coordinate)
{
    using namespace seqan;

    DPMatrix_<char, FullDPMatrix> dpMatrix;

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);
    resize(dpMatrix);

    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 0, DPMatrixDimension_::HORIZONTAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 0, DPMatrixDimension_::VERTICAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 4, DPMatrixDimension_::HORIZONTAL), 1u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 4, DPMatrixDimension_::VERTICAL), 1u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 8, DPMatrixDimension_::HORIZONTAL), 2u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 8, DPMatrixDimension_::VERTICAL), 2u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 11, DPMatrixDimension_::HORIZONTAL), 3u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 11, DPMatrixDimension_::VERTICAL), 2u);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_sparse_coordinate)
{
    using namespace seqan;

    DPMatrix_<char, SparseDPMatrix> dpMatrix;

    setLength(dpMatrix, DPMatrixDimension_::HORIZONTAL, 4);
    setLength(dpMatrix, DPMatrixDimension_::VERTICAL, 3);
    resize(dpMatrix);

    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 0, DPMatrixDimension_::HORIZONTAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 0, DPMatrixDimension_::VERTICAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 1, DPMatrixDimension_::HORIZONTAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 1, DPMatrixDimension_::VERTICAL), 1u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 2, DPMatrixDimension_::HORIZONTAL), 0u);
    SEQAN_ASSERT_EQ(coordinate(dpMatrix, 2, DPMatrixDimension_::VERTICAL), 2u);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_begin_standard)
{
    testAlignmentDPMatrixBegin(seqan::Standard());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_begin_rooted)
{
    testAlignmentDPMatrixBegin(seqan::Rooted());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_end_standard)
{
    testAlignmentDPMatrixEnd(seqan::Standard());
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_end_rooted)
{
    testAlignmentDPMatrixEnd(seqan::Rooted());
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_MATRIX_H_
