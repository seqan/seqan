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
// Tests for quality related functions and metafunctions.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_QUALITIES_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_QUALITIES_H_

// ---------------------------------------------------------------------------
// Test Metafunctions.  Test the default implementations where possible.
// ---------------------------------------------------------------------------

namespace seqan {

struct MyType_;

// For testing the const -> non-const implementation.
template <>
struct QualityValueSize<MyType_>
{
    enum { VALUE = 3 };
};

}  // namespace seqan

SEQAN_DEFINE_TEST(test_basic_alphabet_qualities_quality_value_size_metafunction)
{
    using namespace seqan;

    // Make sure the symbol exist.
    typedef QualityValueSize<int> TQualityValueSize SEQAN_UNUSED_TYPEDEF;

    SEQAN_ASSERT_EQ(+(QualityValueSize<char>::VALUE), 256);  // TODO(holtgrew): Possibly remove default implementation.
    SEQAN_ASSERT_EQ(+(QualityValueSize<MyType_ const>::VALUE), 3);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_qualities_quality_has_qualities_metafunction)
{
    using namespace seqan;

    SEQAN_ASSERT_NOT(+HasQualities<int>::VALUE);
    SEQAN_ASSERT(+HasQualities<Dna5Q>::VALUE);
    SEQAN_ASSERT(+HasQualities<DnaQ>::VALUE);
}

// ---------------------------------------------------------------------------
// Test Functions.
// ---------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_alphabet_qualities_convert_quality)
{
    using namespace seqan;

    char c1 = '\0';
    convertQuality(c1, 33);
    SEQAN_ASSERT_EQ(c1, 'B');
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_QUALITIES_H_
