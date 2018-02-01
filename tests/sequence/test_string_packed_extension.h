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
        SEQAN_ASSERT_EQ(bitScanForward(str), std::numeric_limits<TPosition>::max());
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
        SEQAN_ASSERT_EQ(bitScanForward(str), std::numeric_limits<TPosition>::max());
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
    testStringPackedExtensionBitScanForward();
    testStringPackedExtensionBitScanReverse();
}

#endif // TESTS_SEQUENCE_TEST_STRING_PACKED_EXTENSION_H_
