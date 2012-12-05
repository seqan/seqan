// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Tests for adaption implementations of builtin types to alphabet concepts.
//
// The fulfillment of concepts is already tested in
// test_basic_alphabet_concepts.h.
// ==========================================================================

#ifndef SEQAN_CORE_TESTS_BASIC_TEST_BASIC_ALPHABET_ADAPT_BUILTINS_H_
#define SEQAN_CORE_TESTS_BASIC_TEST_BASIC_ALPHABET_ADAPT_BUILTINS_H_

// Test metafunction IsCharType<>.
SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_metafunction_is_char_type)
{
    using namespace seqan;

    SEQAN_ASSERT_NOT(+(IsCharType<bool>::VALUE));
    SEQAN_ASSERT_NOT(+(IsCharType<bool const>::VALUE));
    SEQAN_ASSERT(+(IsCharType<char>::VALUE));
    SEQAN_ASSERT(+(IsCharType<char const>::VALUE));
    SEQAN_ASSERT(+(IsCharType<wchar_t>::VALUE));
    SEQAN_ASSERT(+(IsCharType<wchar_t const>::VALUE));
}

// ---------------------------------------------------------------------------
// Test Alphabet Concept Implementation
// ---------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_bool)
{
    using namespace seqan;

    // Alphabet Concept
    {
        bool b = false;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<bool>::VALUE), 1);
        assign(b, true);
        SEQAN_ASSERT_EQ(b, true);
    }

    // Ordered Alphabet Concept
    {
        bool b = false, c = true;
        
        SEQAN_ASSERT_EQ(minValue(bool()), false);
        SEQAN_ASSERT_EQ(minValue<bool>(), false);
        SEQAN_ASSERT_EQ(+(MinValue<bool>::VALUE), 0/*false*/);
        SEQAN_ASSERT_EQ(maxValue(bool()), true);
        SEQAN_ASSERT_EQ(maxValue<bool>(), true);
        SEQAN_ASSERT_EQ(+(MaxValue<bool>::VALUE), 1/*true*/);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        bool b = false;

        SEQAN_ASSERT_EQ(ordValue(false), 0u);
        SEQAN_ASSERT_EQ(ordValue(true), 1u);
        SEQAN_ASSERT_EQ(+ValueSize<bool>::VALUE, 2u);
        SEQAN_ASSERT_EQ(valueSize<bool>(), 2u);

        b = 1;
        SEQAN_ASSERT_EQ(b, true);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_char)
{
    using namespace seqan;

    // Alphabet Concept
    {
        bool b = false;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<char>::VALUE), 8u);
        assign(b, true);
        SEQAN_ASSERT_EQ(b, true);
    }

    // Ordered Alphabet Concept
    {
        char b = false, c = true;
#if !defined(_MSC_VER) || defined(_CHAR_UNSIGNED)
        SEQAN_ASSERT_EQ(minValue(char()), '\0');
        SEQAN_ASSERT_EQ(minValue<char>(), '\0');
        SEQAN_ASSERT_EQ(+(MinValue<char>::VALUE), '\0');
        // TODO(holtgrew): Is the following correct?
        SEQAN_ASSERT_EQ(maxValue(char()), char(-1));
        SEQAN_ASSERT_EQ(maxValue<char>(), char(-1));
        SEQAN_ASSERT_EQ(+(MaxValue<char>::VALUE), char(-1));
        SEQAN_ASSERT(b < c);
#else  // #if !defined(_MSC_VER) || defined(_CHAR_UNSIGNED)
        SEQAN_ASSERT_EQ(minValue(char()), -128);
        SEQAN_ASSERT_EQ(minValue<char>(), -128);
        SEQAN_ASSERT_EQ(+(MinValue<char>::VALUE), -128);
        // TODO(holtgrew): Is the following correct?
        SEQAN_ASSERT_EQ(maxValue(char()), 127);
        SEQAN_ASSERT_EQ(maxValue<char>(), 127);
        SEQAN_ASSERT_EQ(+(MaxValue<char>::VALUE), 127);
        SEQAN_ASSERT(b < c);
#endif  // #if !defined(_MSC_VER) || defined(_CHAR_UNSIGNED)
    }

    // Finite Ordered Alphabet
    {
        char b = 0;

        SEQAN_ASSERT_EQ(ordValue('-'), 45u);
        SEQAN_ASSERT_EQ(ordValue('Z'), 90u);
        SEQAN_ASSERT_EQ(+ValueSize<char>::VALUE, 256u);
        SEQAN_ASSERT_EQ(valueSize<char>(), 256u);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Alphabet With Gaps
    {
        SEQAN_ASSERT_EQ(gapValueImpl((char *)0), '-');
        SEQAN_ASSERT_EQ(gapValue<char>(), '-');
    }

    // Alphabet With Unknown Value
    {
        SEQAN_ASSERT_EQ(unknownValueImpl((char *)0), 'N');
        SEQAN_ASSERT_EQ(unknownValue<char>(), 'N');
    }

    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_short)
{
    using namespace seqan;

    // Alphabet Concept
    {
        short b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<short>::VALUE), 16u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        short b = 0, c = 42;

        SEQAN_ASSERT_LEQ(minValue(short()), -32768);
        SEQAN_ASSERT_LEQ(minValue<short>(), -32768);
        SEQAN_ASSERT_LEQ(+(MinValue<short>::VALUE), -32768);
        SEQAN_ASSERT_GEQ(maxValue(short()), 32767);
        SEQAN_ASSERT_GEQ(maxValue<short>(), 32767);
        SEQAN_ASSERT_GEQ(+(MaxValue<short>::VALUE), 32767);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        short b = 0;

        // TODO(holtgrew): This is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(short(-1)), 65535u);
        SEQAN_ASSERT_EQ(ordValue(short(1)), 1u);
        SEQAN_ASSERT_EQ(+ValueSize<short>::VALUE, 65536u);
        SEQAN_ASSERT_EQ(valueSize<short>(), 65536u);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_int)
{
    using namespace seqan;

    // Alphabet Concept
    {
        int b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<int>::VALUE), 32u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        int b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(int()), -(int)2147483648u);
        SEQAN_ASSERT_EQ(minValue<int>(), -(int)2147483648u);
        SEQAN_ASSERT_EQ(+(MinValue<int>::VALUE), -(int)2147483648u);
        SEQAN_ASSERT_EQ(maxValue(int()), 2147483647);
        SEQAN_ASSERT_EQ(maxValue<int>(), 2147483647);
        SEQAN_ASSERT_EQ(+(MaxValue<int>::VALUE), 2147483647);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        int b = 0;

        // TODO(holtgrew): This is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(int(-1)), 4294967295u);
        SEQAN_ASSERT_EQ(ordValue(int(1)), 1u);
        SEQAN_ASSERT_EQ(+ValueSize<int>::VALUE, 4294967296ull);
        SEQAN_ASSERT_EQ(valueSize<int>(), 4294967296ull);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_long)
{
    using namespace seqan;

    // Alphabet Concept
    {
        long b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<long>::VALUE), sizeof(long) * 8);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        long b = 0, c = 42;

        SEQAN_ASSERT_LEQ(minValue(long()), +(MinValue<int>::VALUE));
        SEQAN_ASSERT_LEQ(minValue<long>(), +(MinValue<int>::VALUE));
        SEQAN_ASSERT_LEQ(+(MinValue<long>::VALUE), +(MinValue<int>::VALUE));
        SEQAN_ASSERT_GEQ(maxValue(long()), +(MaxValue<int>::VALUE));
        SEQAN_ASSERT_GEQ(maxValue<long>(), +(MaxValue<int>::VALUE));
        SEQAN_ASSERT_GEQ(+(MaxValue<long>::VALUE), +(MaxValue<int>::VALUE));
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        long b = 0;

        // TODO(holtgrew): This is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(long(-1)), 4294967295u);
        SEQAN_ASSERT_EQ(ordValue(long(1)), 1u);
#if SEQAN_IS_32_BIT
        SEQAN_ASSERT_EQ(+ValueSize<long>::VALUE, 4294967296ull);
        SEQAN_ASSERT_EQ(valueSize<long>(), 4294967296ull);
#else  // #if SEQAN_IS_32_BIT
        if (sizeof(long) == 8u)  // long has 64 bit
        {
            SEQAN_ASSERT_EQ(+ValueSize<long>::VALUE, 0u);
            SEQAN_ASSERT_EQ(valueSize<long>(), 0u);
        }
        else
        {
            SEQAN_ASSERT_EQ(+ValueSize<long>::VALUE, 4294967296ull);
            SEQAN_ASSERT_EQ(valueSize<long>(), 4294967296ull);
        }
#endif  // #if SEQAN_IS_32_BIT

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_int8)
{
    using namespace seqan;

    // Alphabet Concept
    {
        __int8 b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<__int8>::VALUE), 8u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        __int8 b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(__int8()), -128);
        SEQAN_ASSERT_EQ(minValue<__int8>(), -128);
        SEQAN_ASSERT_EQ(+(MinValue<__int8>::VALUE), -128);
        SEQAN_ASSERT_GEQ(maxValue(__int8()), 127);
        SEQAN_ASSERT_GEQ(maxValue<__int8>(), 127);
        SEQAN_ASSERT_GEQ(+(MaxValue<__int8>::VALUE), 127);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        __int8 b = 0;

        // TODO(holtgrew): This is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(__int8(-1)), 255u);
        SEQAN_ASSERT_EQ(ordValue(__int8(1)), 1u);
        SEQAN_ASSERT_EQ(+ValueSize<__int8>::VALUE, 256u);
        SEQAN_ASSERT_EQ(valueSize<__int8>(), 256u);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_uint8)
{
    using namespace seqan;

    // Alphabet Concept
    {
        __uint8 b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<__uint8>::VALUE), 8u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        __uint8 b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(__uint8()), 0u);
        SEQAN_ASSERT_EQ(minValue<__uint8>(), 0u);
        SEQAN_ASSERT_EQ(+(MinValue<__uint8>::VALUE), 0);
        SEQAN_ASSERT_EQ(maxValue(__uint8()), 255u);
        SEQAN_ASSERT_EQ(maxValue<__uint8>(), 255u);
        SEQAN_ASSERT_EQ(+(MaxValue<__uint8>::VALUE), 255);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        __uint8 b = 0;

        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(__uint8(23)), 23u);
        SEQAN_ASSERT_EQ(ordValue(__uint8(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<__uint8>::VALUE, 256u);
        SEQAN_ASSERT_EQ(valueSize<__uint8>(), 256u);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_int16)
{
    using namespace seqan;

    // Alphabet Concept
    {
        __int16 b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<__int16>::VALUE), 16u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        __int16 b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(__int16()), -32768);
        SEQAN_ASSERT_EQ(minValue<__int16>(), -32768);
        SEQAN_ASSERT_EQ(+(MinValue<__int16>::VALUE), -32768);
        SEQAN_ASSERT_EQ(maxValue(__int16()), 32767);
        SEQAN_ASSERT_EQ(maxValue<__int16>(), 32767);
        SEQAN_ASSERT_EQ(+(MaxValue<__int16>::VALUE), 32767);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        __int16 b = 0;

        // The following is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(__int16(-23)), 65513u);
        SEQAN_ASSERT_EQ(ordValue(__int16(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<__int16>::VALUE, 65536u);
        SEQAN_ASSERT_EQ(valueSize<__int16>(), 65536u);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_uint16)
{
    using namespace seqan;

    // Alphabet Concept
    {
        __uint16 b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<__uint16>::VALUE), 16u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        __uint16 b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(__uint16()), 0u);
        SEQAN_ASSERT_EQ(minValue<__uint16>(), 0u);
        SEQAN_ASSERT_EQ(+(MinValue<__uint16>::VALUE), 0);
        SEQAN_ASSERT_EQ(maxValue(__uint16()), 65535u);
        SEQAN_ASSERT_EQ(maxValue<__uint16>(), 65535u);
        SEQAN_ASSERT_EQ(+(MaxValue<__uint16>::VALUE), 65535);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        __uint16 b = 0;

        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(__uint16(23)), 23u);
        SEQAN_ASSERT_EQ(ordValue(__uint16(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<__uint16>::VALUE, 65536u);
        SEQAN_ASSERT_EQ(valueSize<__uint16>(), 65536u);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_int32)
{
    using namespace seqan;

    // Alphabet Concept
    {
        __int32 b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<__int32>::VALUE), 32u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        __int32 b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(__int32()), -(int)2147483648u);
        SEQAN_ASSERT_EQ(minValue<__int32>(), -(int)2147483648u);
        SEQAN_ASSERT_EQ(+(MinValue<__int32>::VALUE), -(int)2147483648u);
        SEQAN_ASSERT_EQ(maxValue(__int32()), 2147483647);
        SEQAN_ASSERT_EQ(maxValue<__int32>(), 2147483647);
        SEQAN_ASSERT_EQ(+(MaxValue<__int32>::VALUE), 2147483647);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        __int32 b = 0;

        // The following is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(__int32(-23)), 4294967273u);
        SEQAN_ASSERT_EQ(ordValue(__int32(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<__int32>::VALUE, 4294967296ull);
        SEQAN_ASSERT_EQ(valueSize<__int32>(), 4294967296ull);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_uint32)
{
    using namespace seqan;

    // Alphabet Concept
    {
        __uint32 b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<__uint32>::VALUE), 32u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1u);
    }

    // Ordered Alphabet Concept
    {
        __uint32 b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(__uint32()), 0u);
        SEQAN_ASSERT_EQ(minValue<__uint32>(), 0u);
        SEQAN_ASSERT_EQ(+(MinValue<__uint32>::VALUE), 0u);
        SEQAN_ASSERT_EQ(maxValue(__uint32()), 4294967295u);
        SEQAN_ASSERT_EQ(maxValue<__uint32>(), 4294967295u);
        SEQAN_ASSERT_EQ(+(MaxValue<__uint32>::VALUE), 4294967295u);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        __uint32 b = 0;

        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(__uint32(23)), 23u);
        SEQAN_ASSERT_EQ(ordValue(__uint32(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<__uint32>::VALUE, 4294967296ull);
        SEQAN_ASSERT_EQ(valueSize<__uint32>(), 4294967296ull);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1u);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_int64)
{
    using namespace seqan;

    // Alphabet Concept
    {
        __int64 b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<__int64>::VALUE), 64u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        __int64 b = 0, c = 42;

        SEQAN_ASSERT_LT(minValue((__int64)(0)), minValue((__int32)(0)));
        SEQAN_ASSERT_LT(minValue<__int64>(), minValue<__int32>());
        SEQAN_ASSERT_LT(+(MinValue<__int64>::VALUE), +(MinValue<__int32>::VALUE));
        SEQAN_ASSERT_GT(maxValue((__int64)(0)), maxValue((__int32)(0)));
        SEQAN_ASSERT_GT(maxValue<__int64>(), maxValue<__int32>());
        SEQAN_ASSERT_GT(+(MaxValue<__int64>::VALUE), +(MaxValue<__int32>::VALUE));
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        __int64 b = 0;

        // The following is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        // SEQAN_ASSERT_EQ(ordValue((__int64)(-23)), 65513u);
        SEQAN_ASSERT_EQ(ordValue((__int64)(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<__int64>::VALUE, 0u);
        SEQAN_ASSERT_EQ(valueSize<__int64>(), 0u);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_uint64)
{
    using namespace seqan;

    // Alphabet Concept
    {
        __uint64 b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<__uint64>::VALUE), 64u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1u);
    }

    // Ordered Alphabet Concept
    {
        __uint64 b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(__uint64()), 0u);
        SEQAN_ASSERT_EQ(minValue<__uint64>(), 0u);
        SEQAN_ASSERT_EQ(+(MinValue<__uint64>::VALUE), 0u);
        SEQAN_ASSERT_GT(maxValue(__uint64()), maxValue(__uint32()));
        SEQAN_ASSERT_GT(maxValue<__uint64>(), maxValue<__uint32>());
        SEQAN_ASSERT_GT(+(MaxValue<__uint64>::VALUE), +(MaxValue<__uint32>::VALUE));
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        __uint64 b = 0;

        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(__uint64(23)), 23u);
        SEQAN_ASSERT_EQ(ordValue(__uint64(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<__uint64>::VALUE, 0u);
        SEQAN_ASSERT_EQ(valueSize<__uint64>(), 0u);

        b = 1;
        SEQAN_ASSERT_EQ(b, 1u);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}


SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_float)
{
    using namespace seqan;

    // Alphabet Concept
    {
        float b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<float>::VALUE), 32u);
        assign(b, 1.0);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        float b = 0, c = 42;

        SEQAN_ASSERT_LT(minValue(float()), 0);
        SEQAN_ASSERT_LT(minValue<float>(), 0);
        SEQAN_ASSERT_LT(+(MinValue<float>::VALUE), 0);
        SEQAN_ASSERT_GT(maxValue(float()), 0);
        SEQAN_ASSERT_GT(maxValue<float>(), 0);
        SEQAN_ASSERT_GT(+(MaxValue<float>::VALUE), 0);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet - Not Applicable
    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

SEQAN_DEFINE_TEST(test_basic_alphabet_adapt_builtins_concepts_double)
{
    using namespace seqan;

    // Alphabet Concept
    {
        double b = 0;
        
        SEQAN_ASSERT_EQ(+(BitsPerValue<double>::VALUE), 64u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        double b = 0, c = 42;

        SEQAN_ASSERT_LT(minValue(double()), 0);
        SEQAN_ASSERT_LT(minValue<double>(), 0);
        SEQAN_ASSERT_LT(+(MinValue<double>::VALUE), 0);
        SEQAN_ASSERT_GT(maxValue(double()), 0);
        SEQAN_ASSERT_GT(maxValue<double>(), 0);
        SEQAN_ASSERT_GT(+(MaxValue<double>::VALUE), 0);
        SEQAN_ASSERT(b < c);
    }

    // Alphabet With Gaps - Not Applicable
    // Alphabet With Unknown Value - Not Applicable
    // Alphabet With Qualities - Not Applicable
}

#endif  // SEQAN_CORE_TESTS_BASIC_TEST_BASIC_ALPHABET_ADAPT_BUILTINS_H_
