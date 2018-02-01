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
// Tests for the conditional enabling part of the metaprogramming library.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_ENABLE_IF_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_ENABLE_IF_H_

// TODO(holtgrew): Test the macros.

template <typename T>
bool testForEnableIf(
    T const & /*x*/,
    typename seqan::EnableIf2<seqan::IsSameType<T, int>::VALUE>::Type* /*dummy*/ = 0)
{
    return true;
}

template <typename T>
bool testForEnableIf(
    T const & /*x*/,
    typename seqan::DisableIf2<seqan::IsSameType<T, int>::VALUE>::Type* /*dummy*/ = 0)
{
    return false;
}

template <typename T>
bool testForEnableIf2(
    T const & /*x*/,
    typename seqan::EnableIf<typename seqan::IsSameType<T, int>::Type>::Type* /*dummy*/ = 0)
{
    return true;
}

template <typename T>
bool testForEnableIf2(
    T const & /*x*/,
    typename seqan::DisableIf<typename seqan::IsSameType<T, int>::Type>::Type* /*dummy*/ = 0)
{
    return false;
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_enable_if_disable_if)
{
    SEQAN_ASSERT_EQ(true, testForEnableIf(1));
    SEQAN_ASSERT_EQ(false, testForEnableIf(1u));
    SEQAN_ASSERT_EQ(false, testForEnableIf(1.0));
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_enable_if2_disable_if2)
{
    SEQAN_ASSERT_EQ(true, testForEnableIf2(1));
    SEQAN_ASSERT_EQ(false, testForEnableIf2(1u));
    SEQAN_ASSERT_EQ(false, testForEnableIf2(1.0));
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_ENABLE_IF_H_
