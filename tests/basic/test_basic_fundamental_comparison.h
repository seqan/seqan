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
// Tests for fundamental comparison code.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_COMPARISON_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_COMPARISON_H_

// Test for the CompareType default implementation.

namespace seqan {

// Extend CompareTypeImpl for int and unsigned.
template <>
struct CompareTypeImpl<int, unsigned>
{
    typedef long Type;
};

template <>
struct CompareTypeImpl<unsigned, int>
{
    typedef long Type;
};

}  // namespace seqan

SEQAN_DEFINE_TEST(test_basic_fundamental_comparison_compare_type)
{
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<int, int>::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<int, int const>::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<int const, int>::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<int const, int const>::Type, int>::VALUE));

    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<unsigned, int>::Type, long>::VALUE));
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<unsigned, int const>::Type, long>::VALUE));
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<unsigned const, int>::Type, long>::VALUE));
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<unsigned const, int const>::Type, long>::VALUE));

    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<int, unsigned>::Type, long>::VALUE));
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<int, unsigned const>::Type, long>::VALUE));
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<int const, unsigned>::Type, long>::VALUE));
    SEQAN_ASSERT((+SameType_<typename seqan::CompareType<int const, unsigned const>::Type, long>::VALUE));
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_COMPARISON_H_
