// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Tests for SIMD vectors.
// ==========================================================================

#ifndef SEQAN_CORE_TESTS_BASIC_TEST_BASIC_SIMD_VECTOR_H_
#define SEQAN_CORE_TESTS_BASIC_TEST_BASIC_SIMD_VECTOR_H_

#include <random>

#include <seqan/sequence.h>
#include <seqan/misc/bit_twiddling.h>
#include <seqan/basic/basic_simd_vector.h>

#if defined(SEQAN_SIMD_ENABLED)
namespace seqan {

template <int ROWS, typename TVector>
inline void test_matrix_transpose()
{
    typedef typename Value<TVector>::Type TValue;
    typedef TVector TMatrix[LENGTH<TVector>::VALUE];
    const int COLS = LENGTH<TVector>::VALUE;
    
    String<TValue> random;
    resize(random, ROWS * COLS);

    std::mt19937 rng;
    std::uniform_int_distribution<TValue> pdf(0, MaxValue<TValue>::VALUE);
    for (unsigned i = 0; i < length(random); ++i)
        random[i] = pdf(rng);

    TMatrix tmp;
    for (int i = 0; i < ROWS; ++i)
        for (int j = 0; j < COLS; ++j)
            tmp[i][j] = random[i * COLS + j];

//    for(int i=0;i<ROWS;++i)
//        print(std::cout, tmp[i]) << std::endl;

    transpose<ROWS>(tmp);

//    std::cout << std::endl;
//    std::cout << std::endl;
//    for(int i=0;i<DIM;++i)
//        print(std::cout, tmp[i]) << std::endl;
    
    for (int i = 0; i < ROWS; ++i)
        for (int j = 0; j < COLS; ++j)
            SEQAN_ASSERT_EQ(tmp[i][j], random[j * ROWS + i]);
}

}

#ifdef SEQAN_SSE4

SEQAN_DEFINE_TEST(test_basic_simd_shuffle)
{
    seqan::SimdVector<unsigned short, 8>::Type vec;
    seqan::SimdVector<unsigned char, 8>::Type  indices;

    for (int i = 0; i < 8; ++i)
        vec[i] = i * 259 + 3;

    for (int i = 0; i < 8; ++i)
        indices[i] = 7 - i;

    vec = seqan::shuffleVector(vec, indices);

    for (int i = 0; i < 8; ++i)
        SEQAN_ASSERT_EQ(vec[i], (7 - i) * 259 + 3);
}

SEQAN_DEFINE_TEST(test_basic_simd_transpose_8x8)
{
    seqan::test_matrix_transpose<8, seqan::SimdVector<unsigned char, 8>::Type>();
}

SEQAN_DEFINE_TEST(test_basic_simd_transpose_16x16)
{
    seqan::test_matrix_transpose<16, seqan::SimdVector<unsigned char, 16>::Type>();
}

#ifdef __AVX2__

SEQAN_DEFINE_TEST(test_basic_simd_shuffle_avx)
{
    seqan::SimdVector<unsigned short, 16>::Type vec;
    seqan::SimdVector<unsigned char, 16>::Type  indices;
    
    const int perm[] = {1,4,2,6,3,5,0,7};

    for (int i = 0; i < 8; ++i)
    {
        vec[i]   = i * 259 + 3;
        vec[i+8] = i * 432 + 9;
    }

    for (int i = 0; i < 8; ++i)
    {
        indices[i]   = 7 - i;
        indices[i+8] = perm[i];
    }

    vec = seqan::shuffleVector(vec, indices);

    for (int i = 0; i < 8; ++i)
    {
        SEQAN_ASSERT_EQ(vec[i],   (7 - i) * 259 + 3);
        SEQAN_ASSERT_EQ(vec[i+8], perm[i] * 432 + 9);
    }
}


SEQAN_DEFINE_TEST(test_basic_simd_transpose_32x32)
{
    seqan::test_matrix_transpose<32, seqan::SimdVector<unsigned char, 32>::Type >();
}

#endif  // #ifdef __AVX2__
#endif  // #ifdef SEQAN_SSE4
#endif  // SEQAN_SIMD_ENABLED

#endif  // #ifndef SEQAN_CORE_TESTS_BASIC_TEST_BASIC_SIMD_VECTOR_H_
