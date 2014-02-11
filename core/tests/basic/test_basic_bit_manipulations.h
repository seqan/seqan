// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Implements test cases for basic bit manipulations.
// ==========================================================================

#ifndef CORE_TESTS_BASIC_TEST_BASIC_BIT_MANIPULATIONS_H_
#define CORE_TESTS_BASIC_TEST_BASIC_BIT_MANIPULATIONS_H_

#include <seqan/basic/basic_debug.h>
#include <seqan/basic/basic_alphabet.h>
#include <seqan/basic/alphabet_storage.h>

using namespace seqan;

template <typename TValue>
TValue testBasicBitManipulationsBitwiseAnd(TValue const & lhs, TValue const & rhs)
{
    TValue res;
    bitwiseAnd(res, lhs, rhs);
    return res;
}

template <typename TValue>
TValue testBasicBitManipulationsBitwiseOr(TValue const & lhs, TValue const & rhs)
{
    TValue res;
    bitwiseOr(res, lhs, rhs);
    return res;
}

template <typename TValue>
TValue testBasicBitManipulationsBitwiseAndNot(TValue const & lhs, TValue const & rhs)
{
    TValue res;
    bitwiseAndNot(res, lhs, rhs);
    return res;
}

template <typename TValue>
TValue testBasicBitManipulationsBitwiseNot(TValue const & val)
{
    TValue res;
    bitwiseNot(res, val);
    return res;
}

template <typename TValue>
bool testBasicBitManipulationsTestAllZeros(TValue const & val)
{
    return testAllZeros(val);
}

template <typename TValue>
TValue testBasicBitManipulationsSetAllZeros(TValue const & val)
{
    TValue res = val;
    setAllZeros(res);
    return res;
}

template <typename TValue>
bool testBasicBitManipulationsTestAllOnes(TValue const & val)
{
    return testAllOnes(val);
}

template <typename TValue>
TValue testBasicBitManipulationsSetAllOnes(TValue const & val)
{
    TValue res = val;
    setAllOnes(res);
    return res;
}

template <typename TValue>
int testBasicBitManipulationsBitScanForward(TValue const & val)
{
    int index;
    bitScanForward(index, val);
    return index;
}

template <typename TValue>
int testBasicBitManipulationsBitScanReverse(TValue const & val)
{
    int index;
    bitScanReverse(index, val);
    return index;
}


SEQAN_DEFINE_TEST(test_basic_bit_manipulations_bitwise_and)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseAnd(3, 8), 0);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseAnd(3, 9), 1);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseAnd(3, 3), 3);
}

SEQAN_DEFINE_TEST(test_basic_bit_manipulations_bitwise_or)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseOr(3, 8), 11);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseOr(3, 9), 11);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseOr(3, 3), 3);
}

SEQAN_DEFINE_TEST(test_basic_bit_manipulations_bitwise_and_not)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseAndNot(3, 8), 3);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseAndNot(3, 9), 2);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseAndNot(3, 3), 0);
}

SEQAN_DEFINE_TEST(test_basic_bit_manipulations_bitwise_not)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseNot(0u), MaxValue<unsigned>::VALUE);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseNot(3u), MaxValue<unsigned>::VALUE << 2);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitwiseNot(MaxValue<unsigned>::VALUE), 0u);
}

SEQAN_DEFINE_TEST(test_basic_bit_manipulations_test_all_zeros)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsTestAllZeros(0), true);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsTestAllZeros(3), false);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsTestAllZeros(MaxValue<unsigned>::VALUE), false);
}

SEQAN_DEFINE_TEST(test_basic_bit_manipulations_set_all_zeros)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsSetAllZeros(0), 0);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsSetAllZeros(3), 0);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsSetAllZeros(MaxValue<unsigned>::VALUE), 0u);
}

SEQAN_DEFINE_TEST(test_basic_bit_manipulations_test_all_ones)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsTestAllOnes(0), false);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsTestAllOnes(3), false);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsTestAllOnes(MaxValue<unsigned>::VALUE), true);
}

SEQAN_DEFINE_TEST(test_basic_bit_manipulations_set_all_ones)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsSetAllOnes(0u), MaxValue<unsigned>::VALUE);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsSetAllOnes(3u), MaxValue<unsigned>::VALUE);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsSetAllOnes(MaxValue<unsigned>::VALUE), MaxValue<unsigned>::VALUE);
}

SEQAN_DEFINE_TEST(test_basic_bit_manipulations_bit_scan_forward)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitScanForward(1), 0);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitScanForward(MaxValue<unsigned>::VALUE), 0);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitScanForward(10), 1);
}

SEQAN_DEFINE_TEST(test_basic_bit_manipulations_bit_scan_reverse)
{
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitScanReverse(1), 0);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitScanReverse(MaxValue<unsigned>::VALUE), (int)+BitsPerValue<int>::VALUE -1);
    SEQAN_ASSERT_EQ(testBasicBitManipulationsBitScanReverse(10), 3);
}

#endif // CORE_TESTS_BASIC_TEST_BASIC_BIT_MANIPULATIONS_H_
