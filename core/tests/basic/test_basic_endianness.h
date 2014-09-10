// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Unit tests for basic endianness.
// ==========================================================================

#ifndef EXTRAS_TESTS_BASIC_TEST_BASIC_ENDIANNESS_H_
#define EXTRAS_TESTS_BASIC_TEST_BASIC_ENDIANNESS_H_

#include <seqan/basic.h>

template <typename TVal>
TVal _testToHostByteOrder(TVal val)
{
    return endianSwap(val, seqan::HostByteOrder(), seqan::HostByteOrder());
}

template <typename TVal>
TVal _testToOtherByteOrder(TVal val)
{
    return endianSwap(val, seqan::LittleEndian(), seqan::BigEndian());
}

SEQAN_DEFINE_TEST(test_basic_endianness_endian_swap_8)
{
    __uint8 x = 0x12;
    SEQAN_ASSERT_EQ(_testToHostByteOrder(x), (__uint8)0x12);
    SEQAN_ASSERT_EQ(_testToOtherByteOrder(x), (__uint8)0x12);
}

SEQAN_DEFINE_TEST(test_basic_endianness_endian_swap_16)
{
    __uint16 x = 0x1234;
    SEQAN_ASSERT_EQ(_testToHostByteOrder(x), (__uint16)0x1234);
    SEQAN_ASSERT_EQ(_testToOtherByteOrder(x), (__uint16)0x3412);
}

SEQAN_DEFINE_TEST(test_basic_endianness_endian_swap_32)
{
    __uint32 x = 0x12345678;
    SEQAN_ASSERT_EQ(_testToHostByteOrder(x), (__uint32)0x12345678);
    SEQAN_ASSERT_EQ(_testToOtherByteOrder(x), (__uint32)0x78563412);
}

SEQAN_DEFINE_TEST(test_basic_endianness_endian_swap_64)
{
    __uint64 x = 0x123456789ABCDEF0;
    SEQAN_ASSERT_EQ(_testToHostByteOrder(x), (__uint64)0x123456789ABCDEF0);
    SEQAN_ASSERT_EQ(_testToOtherByteOrder(x), (__uint64)0xF0DEBC9A78563412);
}

#endif  // EXTRAS_TESTS_BASIC_TEST_BASIC_ENDIANNESS_H_
