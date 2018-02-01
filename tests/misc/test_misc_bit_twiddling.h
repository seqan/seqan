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

#ifndef SEQAN_TESTS_MISC_TEST_MISC_BIT_TWIDDLING_H_
#define SEQAN_TESTS_MISC_TEST_MISC_BIT_TWIDDLING_H_

#include <sstream>
#include <string>

#include <seqan/basic.h>
#include <seqan/misc/bit_twiddling.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_char)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<char>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<char>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<char>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_signed_char)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed char>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed char>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed char>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_unsigned_char)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned char>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned char>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned char>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_short)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<short>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<short>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<short>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_signed_short)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed short>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed short>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed short>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_unsigned_short)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned short>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned short>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned short>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_int)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_signed_int)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed int>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed int>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed int>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_unsigned_int)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned int>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned int>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned int>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_long)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<long>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<long>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<long>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_signed_long)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed long>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed long>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed long>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_unsigned_long)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned long>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned long>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned long>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_long_long)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<long long>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<long long>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<long long>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_signed_long_long)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed long long>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed long long>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<signed long long>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_unsigned_long_long)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned long long>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned long long>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<unsigned long long>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_int8)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int8_t>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int8_t>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int8_t>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_uint8)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint8_t>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint8_t>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint8_t>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_int16)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int16_t>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int16_t>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int16_t>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_uint16)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint16_t>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint16_t>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint16_t>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_int32)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int32_t>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int32_t>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int32_t>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_uint32)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint32_t>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint32_t>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint32_t>(19)), 3u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_int64)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int64_t>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int64_t>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int64_t>(19)), 3u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int64_t>(0x00000000FFFFFFFFll)), 32u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int64_t>(0xFFFFFFFF00000000ll)), 32u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<int64_t>(0xFFFFFFFFFFFFFFFFll)), 64u);
}

SEQAN_DEFINE_TEST(test_misc_bit_twiddling_pop_count_uint64)
{
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint64_t>(0)), 0u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint64_t>(1)), 1u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint64_t>(19)), 3u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint64_t>(0x00000000FFFFFFFFll)), 32u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint64_t>(0xFFFFFFFF00000000ll)), 32u);
    SEQAN_ASSERT_EQ(seqan::popCount(static_cast<uint64_t>(0xFFFFFFFFFFFFFFFFll)), 64u);
}

#endif  // SEQAN_TESTS_MISC_TEST_MISC_BIT_TWIDDLING_H_
