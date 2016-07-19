// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_ADAPT_BUILTINS_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_ADAPT_BUILTINS_H_

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
        SEQAN_ASSERT_EQ(minValue(char()), '\0');
        SEQAN_ASSERT_EQ(minValue<char>(), '\0');
        SEQAN_ASSERT_EQ(+(MinValue<char>::VALUE), '\0');
        // TODO(holtgrew): Is the following correct?
        SEQAN_ASSERT_EQ(maxValue(char()), char(-1));
        SEQAN_ASSERT_EQ(maxValue<char>(), char(-1));
        SEQAN_ASSERT_EQ(+(MaxValue<char>::VALUE), char(-1));
        SEQAN_ASSERT(b < c);
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

        // note(marehr): -2147483648 would produce the compiler warning C4146:
        // unary minus operator applied to unsigned type, result still unsigned
        SEQAN_ASSERT_EQ(minValue(int()), -2147483647 - 1);
        SEQAN_ASSERT_EQ(minValue<int>(), -2147483647 - 1);
        SEQAN_ASSERT_EQ(+(MinValue<int>::VALUE), -2147483647 - 1);
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
        int8_t b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<int8_t>::VALUE), 8u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        int8_t b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(int8_t()), -128);
        SEQAN_ASSERT_EQ(minValue<int8_t>(), -128);
        SEQAN_ASSERT_EQ(+(MinValue<int8_t>::VALUE), -128);
        SEQAN_ASSERT_GEQ(maxValue(int8_t()), 127);
        SEQAN_ASSERT_GEQ(maxValue<int8_t>(), 127);
        SEQAN_ASSERT_GEQ(+(MaxValue<int8_t>::VALUE), 127);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        int8_t b = 0;

        // TODO(holtgrew): This is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(int8_t(-1)), 255u);
        SEQAN_ASSERT_EQ(ordValue(int8_t(1)), 1u);
        SEQAN_ASSERT_EQ(+ValueSize<int8_t>::VALUE, 256u);
        SEQAN_ASSERT_EQ(valueSize<int8_t>(), 256u);

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
        uint8_t b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<uint8_t>::VALUE), 8u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        uint8_t b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(uint8_t()), 0u);
        SEQAN_ASSERT_EQ(minValue<uint8_t>(), 0u);
        SEQAN_ASSERT_EQ(+(MinValue<uint8_t>::VALUE), 0);
        SEQAN_ASSERT_EQ(maxValue(uint8_t()), 255u);
        SEQAN_ASSERT_EQ(maxValue<uint8_t>(), 255u);
        SEQAN_ASSERT_EQ(+(MaxValue<uint8_t>::VALUE), 255);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        uint8_t b = 0;

        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(uint8_t(23)), 23u);
        SEQAN_ASSERT_EQ(ordValue(uint8_t(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<uint8_t>::VALUE, 256u);
        SEQAN_ASSERT_EQ(valueSize<uint8_t>(), 256u);

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
        int16_t b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<int16_t>::VALUE), 16u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        int16_t b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(int16_t()), -32768);
        SEQAN_ASSERT_EQ(minValue<int16_t>(), -32768);
        SEQAN_ASSERT_EQ(+(MinValue<int16_t>::VALUE), -32768);
        SEQAN_ASSERT_EQ(maxValue(int16_t()), 32767);
        SEQAN_ASSERT_EQ(maxValue<int16_t>(), 32767);
        SEQAN_ASSERT_EQ(+(MaxValue<int16_t>::VALUE), 32767);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        int16_t b = 0;

        // The following is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(int16_t(-23)), 65513u);
        SEQAN_ASSERT_EQ(ordValue(int16_t(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<int16_t>::VALUE, 65536u);
        SEQAN_ASSERT_EQ(valueSize<int16_t>(), 65536u);

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
        uint16_t b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<uint16_t>::VALUE), 16u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        uint16_t b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(uint16_t()), 0u);
        SEQAN_ASSERT_EQ(minValue<uint16_t>(), 0u);
        SEQAN_ASSERT_EQ(+(MinValue<uint16_t>::VALUE), 0);
        SEQAN_ASSERT_EQ(maxValue(uint16_t()), 65535u);
        SEQAN_ASSERT_EQ(maxValue<uint16_t>(), 65535u);
        SEQAN_ASSERT_EQ(+(MaxValue<uint16_t>::VALUE), 65535);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        uint16_t b = 0;

        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(uint16_t(23)), 23u);
        SEQAN_ASSERT_EQ(ordValue(uint16_t(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<uint16_t>::VALUE, 65536u);
        SEQAN_ASSERT_EQ(valueSize<uint16_t>(), 65536u);

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
        int32_t b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<int32_t>::VALUE), 32u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        int32_t b = 0, c = 42;

        // note(marehr): -2147483648 would produce the compiler warning C4146:
        // unary minus operator applied to unsigned type, result still unsigned
        SEQAN_ASSERT_EQ(minValue(int32_t()), -2147483647 - 1);
        SEQAN_ASSERT_EQ(minValue<int32_t>(), -2147483647 - 1);
        SEQAN_ASSERT_EQ(+(MinValue<int32_t>::VALUE), -2147483647 - 1);
        SEQAN_ASSERT_EQ(maxValue(int32_t()), 2147483647);
        SEQAN_ASSERT_EQ(maxValue<int32_t>(), 2147483647);
        SEQAN_ASSERT_EQ(+(MaxValue<int32_t>::VALUE), 2147483647);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        int32_t b = 0;

        // The following is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(int32_t(-23)), 4294967273u);
        SEQAN_ASSERT_EQ(ordValue(int32_t(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<int32_t>::VALUE, 4294967296ull);
        SEQAN_ASSERT_EQ(valueSize<int32_t>(), 4294967296ull);

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
        uint32_t b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<uint32_t>::VALUE), 32u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1u);
    }

    // Ordered Alphabet Concept
    {
        uint32_t b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(uint32_t()), 0u);
        SEQAN_ASSERT_EQ(minValue<uint32_t>(), 0u);
        SEQAN_ASSERT_EQ(+(MinValue<uint32_t>::VALUE), 0u);
        SEQAN_ASSERT_EQ(maxValue(uint32_t()), 4294967295u);
        SEQAN_ASSERT_EQ(maxValue<uint32_t>(), 4294967295u);
        SEQAN_ASSERT_EQ(+(MaxValue<uint32_t>::VALUE), 4294967295u);
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        uint32_t b = 0;

        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(uint32_t(23)), 23u);
        SEQAN_ASSERT_EQ(ordValue(uint32_t(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<uint32_t>::VALUE, 4294967296ull);
        SEQAN_ASSERT_EQ(valueSize<uint32_t>(), 4294967296ull);

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
        int64_t b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<int64_t>::VALUE), 64u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1);
    }

    // Ordered Alphabet Concept
    {
        int64_t b = 0, c = 42;

        SEQAN_ASSERT_LT(minValue((int64_t)(0)), minValue((int32_t)(0)));
        SEQAN_ASSERT_LT(minValue<int64_t>(), minValue<int32_t>());
        SEQAN_ASSERT_LT(+(MinValue<int64_t>::VALUE), +(MinValue<int32_t>::VALUE));
        SEQAN_ASSERT_GT(maxValue((int64_t)(0)), maxValue((int32_t)(0)));
        SEQAN_ASSERT_GT(maxValue<int64_t>(), maxValue<int32_t>());
        SEQAN_ASSERT_GT(+(MaxValue<int64_t>::VALUE), +(MaxValue<int32_t>::VALUE));
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        int64_t b = 0;

        // The following is most probably wrong.
        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        // SEQAN_ASSERT_EQ(ordValue((int64_t)(-23)), 65513u);
        SEQAN_ASSERT_EQ(ordValue((int64_t)(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<int64_t>::VALUE, 0u);
        SEQAN_ASSERT_EQ(valueSize<int64_t>(), 0u);

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
        uint64_t b = 0;

        SEQAN_ASSERT_EQ(+(BitsPerValue<uint64_t>::VALUE), 64u);
        assign(b, 1);
        SEQAN_ASSERT_EQ(b, 1u);
    }

    // Ordered Alphabet Concept
    {
        uint64_t b = 0, c = 42;

        SEQAN_ASSERT_EQ(minValue(uint64_t()), 0u);
        SEQAN_ASSERT_EQ(minValue<uint64_t>(), 0u);
        SEQAN_ASSERT_EQ(+(MinValue<uint64_t>::VALUE), 0u);
        SEQAN_ASSERT_GT(maxValue(uint64_t()), maxValue(uint32_t()));
        SEQAN_ASSERT_GT(maxValue<uint64_t>(), maxValue<uint32_t>());
        SEQAN_ASSERT_GT(+(MaxValue<uint64_t>::VALUE), +(MaxValue<uint32_t>::VALUE));
        SEQAN_ASSERT(b < c);
    }

    // Finite Ordered Alphabet
    {
        uint64_t b = 0;

        SEQAN_ASSERT_EQ(ordValue(0), 0u);
        SEQAN_ASSERT_EQ(ordValue(uint64_t(23)), 23u);
        SEQAN_ASSERT_EQ(ordValue(uint64_t(42)), 42u);
        SEQAN_ASSERT_EQ(+ValueSize<uint64_t>::VALUE, 0u);
        SEQAN_ASSERT_EQ(valueSize<uint64_t>(), 0u);

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

#endif  // SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_ADAPT_BUILTINS_H_
