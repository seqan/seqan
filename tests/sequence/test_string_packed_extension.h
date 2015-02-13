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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements additional tests for the packed string.
// ==========================================================================

#ifndef TESTS_SEQUENCE_TEST_STRING_PACKED_EXTENSION_H_
#define TESTS_SEQUENCE_TEST_STRING_PACKED_EXTENSION_H_

#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

void testStringPackedExtensionBitwiseAnd()
{
    typedef String<bool, Packed<> > TBitString;

    {  // Test standard case.
        TBitString str1;
        TBitString str2;
        TBitString res;

        resize(str1, 120, false, Exact());
        resize(str2, 120, false, Exact());

        str1[0] = true;
        str1[10] = true;
        str1[64] = true;
        str1[119] = true;

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr1;
        originStr1 << str1;
        std::stringstream originStr2;
        originStr2 << str2;

        res = str1 & str2;

        SEQAN_ASSERT_EQ(res[0], false);
        SEQAN_ASSERT_EQ(res[10], true);
        SEQAN_ASSERT_EQ(res[63], false);
        SEQAN_ASSERT_EQ(res[64], false);
        SEQAN_ASSERT_EQ(res[119], true);

        std::stringstream testStr1;
        testStr1 << str1;
        std::stringstream testStr2;
        testStr2 << str2;

        SEQAN_ASSERT_EQ(originStr1.str(), testStr1.str());
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());

    }

    {  // Test empty strings.
        TBitString str1;
        TBitString str2;
        TBitString res;

        resize(str2, 120, false, Exact());

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr2;
        originStr2 << str2;

        res = str1 & str2;
        SEQAN_ASSERT_EQ(empty(res), true);
        SEQAN_ASSERT_EQ(empty(host(res)), true);

        std::stringstream testStr2;
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());

        res = str2 & str1;
        SEQAN_ASSERT_EQ(empty(res), true);
        SEQAN_ASSERT_EQ(empty(host(res)), true);

        testStr2.str("");
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }
}

void testStringPackedExtensionBitwiseAndAssign()
{
    typedef String<bool, Packed<> > TBitString;

    {  // Test standard case.
        TBitString str1;
        TBitString str2;

        resize(str1, 120, false, Exact());
        resize(str2, 120, false, Exact());

        str1[0] = true;
        str1[10] = true;
        str1[64] = true;
        str1[119] = true;

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr1;
        originStr1 << str1;
        std::stringstream originStr2;
        originStr2 << str2;

        str1 &= str2;

        SEQAN_ASSERT_EQ(str1[0], false);
        SEQAN_ASSERT_EQ(str1[10], true);
        SEQAN_ASSERT_EQ(str1[63], false);
        SEQAN_ASSERT_EQ(str1[64], false);
        SEQAN_ASSERT_EQ(str1[119], true);

        std::stringstream testStr1;
        testStr1 << str1;
        std::stringstream testStr2;
        testStr2 << str2;

        SEQAN_ASSERT_NEQ(originStr1.str(), testStr1.str());
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }

    {  // Test empty strings.
        TBitString str1;
        TBitString str2;

        resize(str2, 120, false, Exact());

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr2;
        originStr2 << str2;

        str1 &= str2;
        SEQAN_ASSERT_EQ(empty(str1), true);
        SEQAN_ASSERT_EQ(empty(host(str1)), true);

        std::stringstream testStr2;
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());

        str2 &= str1;
        testStr2.str("");
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }
}

void testStringPackedExtensionBitwiseOr()
{
    typedef String<bool, Packed<> > TBitString;

    {
        TBitString str1;
        TBitString str2;
        TBitString res;

        resize(str1, 120, false, Exact());
        resize(str2, 120, false, Exact());

        str1[0] = true;
        str1[10] = true;
        str1[64] = true;
        str1[119] = true;

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr1;
        originStr1 << str1;
        std::stringstream originStr2;
        originStr2 << str2;

        res = str1 | str2;

        SEQAN_ASSERT_EQ(res[0], true);
        SEQAN_ASSERT_EQ(res[10], true);
        SEQAN_ASSERT_EQ(res[63], true);
        SEQAN_ASSERT_EQ(res[64], true);
        SEQAN_ASSERT_EQ(res[119], true);
        std::stringstream testStr1;
        testStr1 << str1;
        std::stringstream testStr2;
        testStr2 << str2;

        SEQAN_ASSERT_EQ(originStr1.str(), testStr1.str());
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }

    {  // Test empty strings.
        TBitString str1;
        TBitString str2;
        TBitString res;

        resize(str2, 120, false, Exact());

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr2;
        originStr2 << str2;

        res = str1 | str2;
        SEQAN_ASSERT_EQ(empty(res), true);
        SEQAN_ASSERT_EQ(empty(host(res)), true);

        std::stringstream testStr2;
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());

        res = str2 | str1;
        SEQAN_ASSERT_EQ(empty(res), true);
        SEQAN_ASSERT_EQ(empty(host(res)), true);

        testStr2.str("");
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }
}

void testStringPackedExtensionBitwiseOrAssign()
{
    typedef String<bool, Packed<> > TBitString;

    {
        TBitString str1;
        TBitString str2;

        resize(str1, 120, false, Exact());
        resize(str2, 120, false, Exact());

        str1[0] = true;
        str1[10] = true;
        str1[64] = true;
        str1[119] = true;

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr1;
        originStr1 << str1;
        std::stringstream originStr2;
        originStr2 << str2;

        str1 |= str2;

        SEQAN_ASSERT_EQ(str1[0], true);
        SEQAN_ASSERT_EQ(str1[10], true);
        SEQAN_ASSERT_EQ(str1[63], true);
        SEQAN_ASSERT_EQ(str1[64], true);
        SEQAN_ASSERT_EQ(str1[119], true);
        std::stringstream testStr1;
        testStr1 << str1;
        std::stringstream testStr2;
        testStr2 << str2;

        SEQAN_ASSERT_NEQ(originStr1.str(), testStr1.str());
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }

    {  // Test empty strings.
        TBitString str1;
        TBitString str2;

        resize(str2, 120, false, Exact());

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr2;
        originStr2 << str2;

        str1 |= str2;
        SEQAN_ASSERT_EQ(empty(str1), true);
        SEQAN_ASSERT_EQ(empty(host(str1)), true);

        std::stringstream testStr2;
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());

        str2 |= str1;
        SEQAN_ASSERT_EQ(empty(str1), true);
        SEQAN_ASSERT_EQ(empty(host(str1)), true);

        testStr2.str("");
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }
}

void testStringPackedExtensionBitwiseXor()
{
    typedef String<bool, Packed<> > TBitString;

    {
        TBitString str1;
        TBitString str2;
        TBitString res;

        resize(str1, 120, false, Exact());
        resize(str2, 120, false, Exact());

        str1[0] = true;
        str1[10] = true;
        str1[64] = true;
        str1[119] = true;

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr1;
        originStr1 << str1;
        std::stringstream originStr2;
        originStr2 << str2;

        res = str1 ^ str2;

        SEQAN_ASSERT_EQ(res[0], true);
        SEQAN_ASSERT_EQ(res[10], false);
        SEQAN_ASSERT_EQ(res[20], false);
        SEQAN_ASSERT_EQ(res[63], true);
        SEQAN_ASSERT_EQ(res[64], true);
        SEQAN_ASSERT_EQ(res[80], false);
        SEQAN_ASSERT_EQ(res[119], false);

        std::stringstream testStr1;
        testStr1 << str1;
        std::stringstream testStr2;
        testStr2 << str2;

        SEQAN_ASSERT_EQ(originStr1.str(), testStr1.str());
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }

    {  // Test empty strings.
        TBitString str1;
        TBitString str2;
        TBitString res;

        resize(str2, 120, false, Exact());

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr2;
        originStr2 << str2;

        res = str1 ^ str2;
        SEQAN_ASSERT_EQ(empty(res), true);
        SEQAN_ASSERT_EQ(empty(host(res)), true);

        std::stringstream testStr2;
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());

        res = str2 ^ str1;
        SEQAN_ASSERT_EQ(empty(res), true);
        SEQAN_ASSERT_EQ(empty(host(res)), true);

        testStr2.str("");
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }
}

void testStringPackedExtensionBitwiseXorAssign()
{
    typedef String<bool, Packed<> > TBitString;

    {
        TBitString str1;
        TBitString str2;

        resize(str1, 120, false, Exact());
        resize(str2, 120, false, Exact());

        str1[0] = true;
        str1[10] = true;
        str1[64] = true;
        str1[119] = true;

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr1;
        originStr1 << str1;
        std::stringstream originStr2;
        originStr2 << str2;

        str1 ^= str2;

        SEQAN_ASSERT_EQ(str1[0], true);
        SEQAN_ASSERT_EQ(str1[10], false);
        SEQAN_ASSERT_EQ(str1[20], false);
        SEQAN_ASSERT_EQ(str1[63], true);
        SEQAN_ASSERT_EQ(str1[64], true);
        SEQAN_ASSERT_EQ(str1[80], false);
        SEQAN_ASSERT_EQ(str1[119], false);
        std::stringstream testStr1;
        testStr1 << str1;
        std::stringstream testStr2;
        testStr2 << str2;

        SEQAN_ASSERT_NEQ(originStr1.str(), testStr1.str());
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }

    {  // Test empty strings.
        TBitString str1;
        TBitString str2;

        resize(str2, 120, false, Exact());

        str2[1] = true;
        str2[10] = true;
        str2[63] = true;
        str2[119] = true;

        std::stringstream originStr2;
        originStr2 << str2;

        str1 ^= str2;
        SEQAN_ASSERT_EQ(empty(str1), true);
        SEQAN_ASSERT_EQ(empty(host(str1)), true);

        std::stringstream testStr2;
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());

        str2 ^= str1;
        SEQAN_ASSERT_EQ(empty(str1), true);
        SEQAN_ASSERT_EQ(empty(host(str1)), true);

        testStr2.str("");
        testStr2 << str2;
        SEQAN_ASSERT_EQ(originStr2.str(), testStr2.str());
    }
}

void testStringPackedExtensionBitwiseNot()
{
    typedef String<bool, Packed<> > TBitString;

    {
        TBitString str;
        TBitString res;

        resize(str, 120, true, Exact());

        str[0] = false;
        str[10] = false;
        str[64] = false;
        str[119] = false;

        res = ~str;

        SEQAN_ASSERT_EQ(res[0], true);
        SEQAN_ASSERT_EQ(res[1], false);
        SEQAN_ASSERT_EQ(res[9], false);
        SEQAN_ASSERT_EQ(res[10], true);
        SEQAN_ASSERT_EQ(res[11], false);
        SEQAN_ASSERT_EQ(res[63], false);
        SEQAN_ASSERT_EQ(res[64], true);
        SEQAN_ASSERT_EQ(res[65], false);
        SEQAN_ASSERT_EQ(res[118], false);
        SEQAN_ASSERT_EQ(res[119], true);
    }

    {  // Empty string
        TBitString str;
        TBitString res;

        res = ~str;

        SEQAN_ASSERT_EQ(empty(host(res)), true);
        SEQAN_ASSERT_EQ(empty(res), true);
    }
}

void testStringPackedExtensionTestAllZeros()
{
    typedef String<bool, Packed<> > TBitString;

    {
        TBitString str;
        resize(str, 120, false, Exact());

        SEQAN_ASSERT_EQ(testAllZeros(str), true);
        str[0] = true;
        SEQAN_ASSERT_EQ(testAllZeros(str), false);
        str[0] = false;
        str[63] = true;
        SEQAN_ASSERT_EQ(testAllZeros(str), false);
        str[63] = false;
        str[64] = true;
        SEQAN_ASSERT_EQ(testAllZeros(str), false);
        str[64] = false;
        str[119] = true;
        SEQAN_ASSERT_EQ(testAllZeros(str), false);
        str[119] = false;
        SEQAN_ASSERT_EQ(testAllZeros(str), true);
    }

    {
        TBitString str;
        SEQAN_ASSERT_EQ(testAllZeros(str), false);
    }
}

void testStringPackedExtensionTestAllOnes()
{
    typedef String<bool, Packed<> > TBitString;
    {
        TBitString str;
        resize(str, 120, true, Exact());

        SEQAN_ASSERT_EQ(testAllOnes(str), true);
        str[0] = false;
        SEQAN_ASSERT_EQ(testAllOnes(str), false);
        str[0] = true;
        str[63] = false;
        SEQAN_ASSERT_EQ(testAllOnes(str), false);
        str[63] = true;
        str[64] = false;
        SEQAN_ASSERT_EQ(testAllOnes(str), false);
        str[64] = true;
        str[119] = false;
        SEQAN_ASSERT_EQ(testAllOnes(str), false);
        str[119] = true;
        SEQAN_ASSERT_EQ(testAllOnes(str), true);
    }

    {
        TBitString str;
        SEQAN_ASSERT_EQ(testAllOnes(str), false);
    }
}

void testStringPackedExtensionTestEqual()
{
    typedef String<bool, Packed<> > TPackedString;
    typedef Value<TPackedString>::Type TValue;
    typedef PackedTraits_<TPackedString> TTraits;
    typedef TTraits::THostValue THostValue;
    typedef THostValue::TBitVector TBitVector;
    typedef Iterator<TPackedString, Standard>::Type TIter;

    {  // Empty host.
        TPackedString strL, strR;

        resize(strL, 1243, Exact());
        SEQAN_ASSERT_NOT(isEqual(strL, strR));
        SEQAN_ASSERT_NOT(strL == strR);
        SEQAN_ASSERT(strL != strR);
    }

    {  // Different Size.
        TPackedString strL, strR;
        resize(strL, +TTraits::VALUES_PER_HOST_VALUE, TValue());
        resize(strR, +TTraits::VALUES_PER_HOST_VALUE * 2, TValue());

        SEQAN_ASSERT_NOT(isEqual(strL, strR));
        SEQAN_ASSERT_NOT(strL == strR);
        SEQAN_ASSERT(strL != strR);
    }

    {  // Test Equal.

        // Random string with size: TTraits::VALUES_PER_HOST_VALUE - 3
        CharString testString = "1110000011010011110011110110011001011111101011010111000000010";
        Iterator<CharString, Standard>::Type testStrIt = begin(testString, Standard());

        TPackedString strL, strR;
        resize(strL, length(testString), Exact());
        for (TIter it = begin(strL, Standard()); it != end(strL, Standard()); ++it, ++testStrIt)
            assignValue(it, (*testStrIt == '1') ? true : false);

        strR = strL;
        // Manipulate wasted bits.
        *(begin(host(strR), Standard()) + 1) |= ~(~static_cast<TBitVector>(0) >> TTraits::WASTED_BITS);
        // Manipulate inactive bits at end.
        *(end(host(strL), Standard()) - 1) |= ((1 << (3 * BitsPerValue<TValue>::VALUE)) - 1);

        SEQAN_ASSERT(isEqual(strL, strR));
        SEQAN_ASSERT(strL == strR);
        SEQAN_ASSERT_NOT(strL != strR);

        // Random string with size: 8 * TTraits::VALUES_PER_HOST_VALUE - 3
        testString = "0101110010111101011101111101001100001001110101100111110010001001"
                     "0100001001110110101000110001101101010100000000010000111110111011"
                     "0001110011001000111111110010111101011010100001111110110100011000"
                     "1111110111011001111000100000001010101111111101010000001001110001"
                     "1001111001110111000001100101000111010001111101111111111010100001"
                     "1100101001010000100101101100101100101000001000101111000110101010"
                     "0011100011100101100111101110010110110111100000000010101101101110"
                     "0011011011111000101111000110011101111001011000010010010000111";

        testStrIt = begin(testString, Standard());
        resize(strL, length(testString), TValue(), Exact());
        for (TIter it = begin(strL, Standard()); it != end(strL, Standard()); ++it)
            assignValue(it, (*testStrIt == '1') ? true : false);

        strR = strL;
        // Manipulate wasted bits.
        *(begin(host(strR), Standard()) + 4) |= ~(~static_cast<TBitVector>(0) >> TTraits::WASTED_BITS);
        // Manipulate inactive bits at end.
        *(end(host(strL), Standard()) - 1) |= ((1 << (3 * BitsPerValue<TValue>::VALUE)) - 1);

        SEQAN_ASSERT(isEqual(strL, strR));
        SEQAN_ASSERT(strL == strR);
        SEQAN_ASSERT_NOT(strL != strR);

        unsigned pos = 3 * TTraits::VALUES_PER_HOST_VALUE + 2;
        TValue tmp = strR[pos];
        strR[pos] = convert<TValue>(~ordValue(tmp));  // Invert value.

        SEQAN_ASSERT_NOT(isEqual(strL, strR));
        SEQAN_ASSERT_NOT(strL == strR);
        SEQAN_ASSERT(strL != strR);

        strR[pos] = tmp;
        *(end(strL, Standard()) - 1) = convert<TValue>(~ordValue(getValue(end(strL, Standard()) - 1)));  // Invert last value.

        SEQAN_ASSERT_NOT(isEqual(strL, strR));
        SEQAN_ASSERT_NOT(strL == strR);
        SEQAN_ASSERT(strL != strR);
    }
}

void testStringPackedExtensionBitScanForward()
{
    typedef String<bool, Packed<> > TBitString;
    typedef Position<TBitString>::Type TPosition;

    {
        TBitString str;
        resize(str, 120, false, Exact());

        str[0] = true;
        str[10] = true;
        str[63] = true;
        str[64] = true;
        str[119] = true;

        SEQAN_ASSERT_EQ(bitScanForward(str), 0u);
        str[0] = false;
        SEQAN_ASSERT_EQ(bitScanForward(str), 10u);
        str[10] = false;
        SEQAN_ASSERT_EQ(bitScanForward(str), 63u);
        str[63] = false;
        SEQAN_ASSERT_EQ(bitScanForward(str), 64u);
        str[64] = false;
        SEQAN_ASSERT_EQ(bitScanForward(str), 119u);
        str[119] = false;
        SEQAN_ASSERT_EQ(bitScanForward(str), 120u);
    }

    {
        TBitString str;
        SEQAN_ASSERT_EQ(bitScanForward(str), MaxValue<TPosition>::VALUE);
    }
}

void testStringPackedExtensionBitScanReverse()
{
    typedef String<bool, Packed<> > TBitString;
    typedef Position<TBitString>::Type TPosition;

    {
        TBitString str;
        resize(str, 120, false, Exact());

        str[0] = true;
        str[10] = true;
        str[63] = true;
        str[64] = true;
        str[119] = true;

        SEQAN_ASSERT_EQ(bitScanReverse(str), 119u);
        str[119] = false;
        SEQAN_ASSERT_EQ(bitScanReverse(str), 64u);
        str[64] = false;
        SEQAN_ASSERT_EQ(bitScanReverse(str), 63u);
        str[63] = false;
        SEQAN_ASSERT_EQ(bitScanReverse(str), 10u);
        str[10] = false;
        SEQAN_ASSERT_EQ(bitScanReverse(str), 0u);
        str[0] = false;
        SEQAN_ASSERT_EQ(bitScanReverse(str), 120u);
    }

    {
        TBitString str;
        SEQAN_ASSERT_EQ(bitScanForward(str), MaxValue<TPosition>::VALUE);
    }
}


SEQAN_DEFINE_TEST(String_Packed_Extension)
{
    testStringPackedExtensionBitwiseAnd();
    testStringPackedExtensionBitwiseAndAssign();
    testStringPackedExtensionBitwiseOr();
    testStringPackedExtensionBitwiseOrAssign();
    testStringPackedExtensionBitwiseXor();
    testStringPackedExtensionBitwiseXorAssign();
    testStringPackedExtensionBitwiseNot();
    testStringPackedExtensionTestAllZeros();
    testStringPackedExtensionTestAllOnes();
    testStringPackedExtensionTestEqual();
    testStringPackedExtensionBitScanForward();
    testStringPackedExtensionBitScanReverse();
}

#endif // TESTS_SEQUENCE_TEST_STRING_PACKED_EXTENSION_H_
