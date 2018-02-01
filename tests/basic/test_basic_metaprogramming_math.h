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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the math metaprogramming code.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_MATH_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_MATH_H_

// In the following, we use the unary plus to create an lvalue such
// that the comparison functionwith "T const &" typed parameter does
// not try to access the address of ::VALUE.

SEQAN_DEFINE_TEST(test_basic_metaprogramming_math_log2)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(+Log2<1>::VALUE, 0u);
    SEQAN_ASSERT_EQ(+Log2<2>::VALUE, 1u);
    SEQAN_ASSERT_EQ(+Log2<3>::VALUE, 2u);
    SEQAN_ASSERT_EQ(+Log2<4>::VALUE, 2u);
    SEQAN_ASSERT_EQ(+Log2<5>::VALUE, 3u);
    SEQAN_ASSERT_EQ(+Log2<7>::VALUE, 3u);
    SEQAN_ASSERT_EQ(+Log2<8>::VALUE, 3u);
    SEQAN_ASSERT_EQ(+Log2<9>::VALUE, 4u);

    SEQAN_ASSERT_EQ(+Log2<127>::VALUE, 7u);
    SEQAN_ASSERT_EQ(+Log2<128>::VALUE, 7u);
    SEQAN_ASSERT_EQ(+Log2<129>::VALUE, 8u);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_math_log2_floor)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(+Log2Floor<1>::VALUE, 0u);
    SEQAN_ASSERT_EQ(+Log2Floor<2>::VALUE, 1u);
    SEQAN_ASSERT_EQ(+Log2Floor<3>::VALUE, 1u);
    SEQAN_ASSERT_EQ(+Log2Floor<4>::VALUE, 2u);
    SEQAN_ASSERT_EQ(+Log2Floor<5>::VALUE, 2u);
    SEQAN_ASSERT_EQ(+Log2Floor<7>::VALUE, 2u);
    SEQAN_ASSERT_EQ(+Log2Floor<8>::VALUE, 3u);
    SEQAN_ASSERT_EQ(+Log2Floor<9>::VALUE, 3u);

    SEQAN_ASSERT_EQ(+Log2Floor<127>::VALUE, 6u);
    SEQAN_ASSERT_EQ(+Log2Floor<128>::VALUE, 7u);
    SEQAN_ASSERT_EQ(+Log2Floor<129>::VALUE, 7u);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_math_log2_power)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ((+Power<0, 0>::VALUE), 1u);
    SEQAN_ASSERT_EQ((+Power<0, 1>::VALUE), 0u);
    SEQAN_ASSERT_EQ((+Power<0, 42>::VALUE), 0u);

    SEQAN_ASSERT_EQ((+Power<1, 0>::VALUE), 1u);
    SEQAN_ASSERT_EQ((+Power<1, 1>::VALUE), 1u);
    SEQAN_ASSERT_EQ((+Power<1, 42>::VALUE), 1u);

    SEQAN_ASSERT_EQ((+Power<2, 0>::VALUE), 1u);
    SEQAN_ASSERT_EQ((+Power<2, 1>::VALUE), 2u);
    SEQAN_ASSERT_EQ((+Power<2, 10>::VALUE), 1024u);
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_MATH_H_
