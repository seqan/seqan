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
// Tests for math functions and metafunctions.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_MATH_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_MATH_H_

// Test math-related metafunctions.  We only need to check that the
// metafunctions exist.  The actual implementations are given for each type,
// or adaptions to the builtins.

SEQAN_DEFINE_TEST(test_basic_alphabet_math_metafunctions)
{
    typedef seqan::MaxValue<int> TMaxValue SEQAN_UNUSED_TYPEDEF;
    typedef seqan::MinValue<int> TMinValue SEQAN_UNUSED_TYPEDEF;
}

// Now, test that the forwards for minValue() and maxValue() work correctly.

struct MyNumber_
{
    int value;

    MyNumber_() : value(0)
    {}

    explicit
    MyNumber_(int v) : value(v)
    {}
};

namespace seqan
{
template <>
struct MinValue<MyNumber_>
{
    static const MyNumber_ VALUE;
};

const MyNumber_ MinValue<MyNumber_>::VALUE = MyNumber_(-1);

template <>
struct MaxValue<MyNumber_>
{
    static const MyNumber_ VALUE;
};

const MyNumber_ MaxValue<MyNumber_>::VALUE = MyNumber_(1);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_math_min_value)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(MinValue<MyNumber_>::VALUE.value, -1);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_math_max_value)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(MaxValue<MyNumber_>::VALUE.value, 1);
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_MATH_H_
