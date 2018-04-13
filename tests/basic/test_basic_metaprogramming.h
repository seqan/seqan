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
// Tests for the SeqAn metaprogramming header.
// ==========================================================================

#ifndef TEST_BASIC_TEST_BASIC_METAPROGRAMMING_H_
#define TEST_BASIC_TEST_BASIC_METAPROGRAMMING_H_

#include <seqan/basic.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_basic_metaprogramming_true)
{
    bool b = True::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_false)
{
    bool b = False::VALUE;
    SEQAN_ASSERT_EQ(b, false);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_eval)
{
    bool b = Eval<true>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = Eval<false>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_or)
{
    bool b = Or<True, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = Or<True, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = Or<False, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = Or<False, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_and)
{
    bool b = And<True, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
    b = And<True, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
    b = And<False, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
    b = And<False, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_is_same_type)
{
    bool b = IsSameType<True, False>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, false);
    b = IsSameType<True, True>::Type::VALUE;
    SEQAN_ASSERT_EQ(b, true);
}

// Helper functions and struct for the test for metaprogrammign Switch.

SEQAN_DEFINE_TEST(test_basic_metaprogramming_log2)
{
    uint64_t x = Log2<1>::VALUE;
    SEQAN_ASSERT_EQ(x, 0u);
    x = Log2<2>::VALUE;
    SEQAN_ASSERT_EQ(x, 1u);
    x = Log2<3>::VALUE;
    SEQAN_ASSERT_EQ(x, 2u);
    x = Log2<4>::VALUE;
    SEQAN_ASSERT_EQ(x, 2u);
    x = Log2<5>::VALUE;
    SEQAN_ASSERT_EQ(x, 3u);
    x = Log2<15>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
    x = Log2<16>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
    x = Log2<17>::VALUE;
    SEQAN_ASSERT_EQ(x, 5u);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_log2_floor)
{
    uint64_t x = Log2Floor<1>::VALUE;
    SEQAN_ASSERT_EQ(x, 0u);
    x = Log2Floor<2>::VALUE;
    SEQAN_ASSERT_EQ(x, 1u);
    x = Log2Floor<3>::VALUE;
    SEQAN_ASSERT_EQ(x, 1u);
    x = Log2Floor<4>::VALUE;
    SEQAN_ASSERT_EQ(x, 2u);
    x = Log2Floor<5>::VALUE;
    SEQAN_ASSERT_EQ(x, 2u);
    x = Log2Floor<15>::VALUE;
    SEQAN_ASSERT_EQ(x, 3u);
    x = Log2Floor<16>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
    x = Log2Floor<17>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_power)
{
    uint64_t x = Power<2, 2>::VALUE;
    SEQAN_ASSERT_EQ(x, 4u);
    x = Power<3, 2>::VALUE;
    SEQAN_ASSERT_EQ(x, 9u);
    x = Power<10, 0>::VALUE;
    SEQAN_ASSERT_EQ(x, 1u);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_unsigned)
{
    bool b = IsSameType<unsigned int, typename MakeUnsigned_<int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<uint64_t, typename MakeUnsigned_<int64_t>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_signed)
{
    bool b = IsSameType<int, typename MakeSigned_<unsigned>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int64_t, typename MakeSigned_<uint64_t>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_remove_const)
{
    bool b = IsSameType<int, const int>::Type::VALUE;
    SEQAN_ASSERT_NOT(b);
    b = IsSameType<int, RemoveConst_<const int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int, RemoveConst_<int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_copy_const)
{
    bool b = IsSameType<int, CopyConst_<int, int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int const, CopyConst_<int, int const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int const, CopyConst_<int const, int const>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
    b = IsSameType<int const, CopyConst_<int const, int>::Type>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_make_is_const)
{
    bool b = IsConst_<int>::Type::VALUE;
    SEQAN_ASSERT_NOT(b);
    b = IsConst_<int const>::Type::VALUE;
    SEQAN_ASSERT(b);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_class_identifier)
{
    void * id1 = ClassIdentifier_<int>::getID();
    void * id2 = ClassIdentifier_<signed int>::getID();
    void * id3 = ClassIdentifier_<double>::getID();
    SEQAN_ASSERT_EQ(id1, id2);
    SEQAN_ASSERT_NEQ(id1, id3);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_memset)
{
    int i[] = { 1, 2, 3 };

    memset<3 * sizeof(int), '\0'>(i);

    SEQAN_ASSERT_EQ(i[0], 0);
    SEQAN_ASSERT_EQ(i[1], 0);
    SEQAN_ASSERT_EQ(i[2], 0);
}

// Functions to demonstrate EnableIf and related functions.
template <typename T>
bool testForEnableIf(
    T const & /*x*/,
    typename seqan::EnableIf<typename IsSameType<T, int>::Type>::Type* /*dummy*/ = 0)
{
    return true;
}

template <typename T>
bool testForEnableIf(
    T const & /*x*/,
    typename seqan::DisableIf<typename IsSameType<T, int>::Type>::Type* /*dummy*/ = 0)
{
    return false;
}

template <typename T>
bool testForEnableIf2(
    T const & /*x*/,
    typename seqan::EnableIf2<IsSameType<T, int>::VALUE>::Type* /*dummy*/ = 0)
{
    return true;
}

template <typename T>
bool testForEnableIf2(
    T const & /*x*/,
    typename seqan::DisableIf2<IsSameType<T, int>::VALUE>::Type* /*dummy*/ = 0)
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

#endif  // TEST_BASIC_TEST_BASIC_METAPROGRAMMING_H_
