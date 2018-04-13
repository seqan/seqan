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
// Tests for alphabet storage related code.
//
// Test default implementations and existance of functions only, the
// individual implementations are tested when checking for concept
// fulfillment.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_STORAGE_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_STORAGE_H_

struct LongStruct_
{
    int64_t a, b;
};

namespace seqan {

template <>
struct ValueSize<LongStruct_>
{
    typedef unsigned Type;
    static const unsigned VALUE = 42;
};

}

SEQAN_DEFINE_TEST(test_basic_alphabet_storage_bits_per_value_metafunction)
{
    using namespace seqan;

    typedef BitsPerValue<int> TBitsPerValue SEQAN_UNUSED_TYPEDEF;  // Check existance.

    SEQAN_ASSERT_EQ(+BitsPerValue<LongStruct_>::VALUE, 128u);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_storage_value_size_metafunction)
{
    using namespace seqan;

    typedef ValueSize<int> TValueSize SEQAN_UNUSED_TYPEDEF;  // Check existance.

    SEQAN_ASSERT_EQ(+ValueSize<LongStruct_>::VALUE, 42u);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_storage_value_size_function)
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(valueSize<LongStruct_>(), 42u);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_storage_integral_for_value_metafunction)
{
}

SEQAN_DEFINE_TEST(test_basic_alphabet_storage_bytes_per_value_metafunction)
{
    using namespace seqan;

    typedef BytesPerValue<int> TBitsPerValue SEQAN_UNUSED_TYPEDEF;  // Check existance.

    SEQAN_ASSERT_EQ(+BytesPerValue<LongStruct_>::VALUE, 16);
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_STORAGE_H_
