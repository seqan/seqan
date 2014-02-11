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
// Implements additional tests for the packed string.
// ==========================================================================

#ifndef CORE_TESTS_SEQUENCE_TEST_STRING_PACKED_EXTENSION_H_
#define CORE_TESTS_SEQUENCE_TEST_STRING_PACKED_EXTENSION_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

using namespace seqan;

void testStringPackedExtensionBitwiseAnd()
{
    typedef String<bool, Packed<> > TBitString;
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

    bitwiseAnd(res, str1, str2);
    SEQAN_ASSERT_EQ(res[0], false);
    SEQAN_ASSERT_EQ(res[10], true);
    SEQAN_ASSERT_EQ(res[63], false);
    SEQAN_ASSERT_EQ(res[64], false);
    SEQAN_ASSERT_EQ(res[119], true);
}

void testStringPackedExtensionBitwiseOr()
{
    typedef String<bool, Packed<> > TBitString;
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

    bitwiseOr(res, str1, str2);
    SEQAN_ASSERT_EQ(res[0], true);
    SEQAN_ASSERT_EQ(res[10], true);
    SEQAN_ASSERT_EQ(res[63], true);
    SEQAN_ASSERT_EQ(res[64], true);
    SEQAN_ASSERT_EQ(res[119], true);
}

void testStringPackedExtensionBitwiseAndNot()
{
    typedef String<bool, Packed<> > TBitString;
    TBitString str1;
    TBitString str2;
    TBitString res;

    resize(str1, 120, false, Exact());
    resize(str2, 120, true, Exact());

    str1[0] = true;
    str1[10] = true;
    str1[64] = true;
    str1[119] = true;

    str2[1] = false;
    str2[10] = false;
    str2[63] = false;
    str2[119] = false;

    bitwiseAndNot(res, str1, str2);
    SEQAN_ASSERT_EQ(res[0], false);
    SEQAN_ASSERT_EQ(res[10], true);
    SEQAN_ASSERT_EQ(res[63], false);
    SEQAN_ASSERT_EQ(res[64], false);
    SEQAN_ASSERT_EQ(res[119], true);
}

void testStringPackedExtensionBitwiseNot()
{
    typedef String<bool, Packed<> > TBitString;
    TBitString str;
    TBitString res;

    resize(str, 120, true, Exact());

    str[0] = false;
    str[10] = false;
    str[64] = false;
    str[119] = false;

    bitwiseNot(res, str);
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

void testStringPackedExtensionTestAllZeros()
{
    typedef String<bool, Packed<> > TBitString;
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

void testStringPackedExtensionSetAllZeros()
{
    typedef String<bool, Packed<> > TBitString;
    TBitString str;

    resize(str, 120, true, Exact());
    SEQAN_ASSERT_EQ(testAllZeros(str), false);
    setAllZeros(str);
    SEQAN_ASSERT_EQ(testAllZeros(str), true);
}

void testStringPackedExtensionTestAllOnes()
{
    typedef String<bool, Packed<> > TBitString;
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

void testStringPackedExtensionSetAllOnes()
{
    typedef String<bool, Packed<> > TBitString;
    TBitString str;

    resize(str, 120, false, Exact());
    SEQAN_ASSERT_EQ(testAllOnes(str), false);
    setAllOnes(str);
    SEQAN_ASSERT_EQ(testAllOnes(str), true);
}

void testStringPackedExtensionBitScanForward()
{
    typedef String<bool, Packed<> > TBitString;
    TBitString str;

    resize(str, 120, false, Exact());

    str[0] = true;
    str[10] = true;
    str[63] = true;
    str[64] = true;
    str[119] = true;

    unsigned index;
    bitScanForward(index, str);
    SEQAN_ASSERT_EQ(index, 0u);
    str[0] = false;
    bitScanForward(index, str);
    SEQAN_ASSERT_EQ(index, 10u);
    str[10] = false;
    bitScanForward(index, str);
    SEQAN_ASSERT_EQ(index, 63u);
    str[63] = false;
    bitScanForward(index, str);
    SEQAN_ASSERT_EQ(index, 64u);
    str[64] = false;
    bitScanForward(index, str);
    SEQAN_ASSERT_EQ(index, 119u);
}

void testStringPackedExtensionBitScanReverse()
{
    typedef String<bool, Packed<> > TBitString;
    TBitString str;

    resize(str, 120, false, Exact());

    str[0] = true;
    str[10] = true;
    str[63] = true;
    str[64] = true;
    str[119] = true;

    unsigned index;
    bitScanReverse(index, str);
    SEQAN_ASSERT_EQ(index, 119u);
    str[119] = false;
    bitScanReverse(index, str);
    SEQAN_ASSERT_EQ(index, 64u);
    str[64] = false;
    bitScanReverse(index, str);
    SEQAN_ASSERT_EQ(index, 63u);
    str[63] = false;
    bitScanReverse(index, str);
    SEQAN_ASSERT_EQ(index, 10u);
    str[10] = false;
    bitScanReverse(index, str);
    SEQAN_ASSERT_EQ(index, 0u);
}


SEQAN_DEFINE_TEST(String_Packed_Extension)
{
    testStringPackedExtensionBitwiseAnd();
    testStringPackedExtensionBitwiseAndNot();
    testStringPackedExtensionBitwiseOr();
    testStringPackedExtensionTestAllZeros();
    testStringPackedExtensionSetAllZeros();
    testStringPackedExtensionTestAllOnes();
    testStringPackedExtensionSetAllOnes();
    testStringPackedExtensionBitScanForward();
    testStringPackedExtensionBitScanReverse();
}

#endif // CORE_TESTS_SEQUENCE_TEST_STRING_PACKED_EXTENSION_H_
