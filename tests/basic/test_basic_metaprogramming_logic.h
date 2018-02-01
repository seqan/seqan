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
// Tests for the logic part of the metaprogramming library.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_LOGIC_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_LOGIC_H_

// --------------------------------------------------------------------------
// Helper metafunction to check for type equalness.
// --------------------------------------------------------------------------

template <typename T1, typename T2>
struct TestTypeEq
{
    typedef bool TValue;
    static const bool VALUE;
};

template <typename T>
struct TestTypeEq<T, T>
{
    typedef bool TValue;
    static const bool VALUE;
};

template <typename T1, typename T2>
const bool TestTypeEq<T1, T2>::VALUE = false;

template <typename T>
const bool TestTypeEq<T, T>::VALUE = true;

// --------------------------------------------------------------------------
// Actual Tests.
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_metaprogramming_logic_bool_type)
{
    using namespace seqan;

    // Test for the values of the VALUE members.
    SEQAN_ASSERT_EQ(+False::VALUE, 0);
    SEQAN_ASSERT_EQ(+True::VALUE, 1);

    // Test for type of the Type members.
    SEQAN_ASSERT_EQ((TestTypeEq<typename False::Type, False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<typename True::Type,  True>::VALUE), true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_logic_eval)
{
    using namespace seqan;

    // Test for the VALUE members.
    SEQAN_ASSERT_EQ(+Eval<false>::VALUE, +false);
    SEQAN_ASSERT_EQ(+Eval<true>::VALUE, +true);

    // Test for the Type members.
    SEQAN_ASSERT_EQ((+TestTypeEq<Eval<false>::Type, False>::VALUE), +true);
    SEQAN_ASSERT_EQ((+TestTypeEq<Eval<true>::Type,  True>::VALUE),  +true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_logic_or)
{
    using namespace seqan;

    // Test for the VALUE members.
    SEQAN_ASSERT_EQ((+Or<False, False>::VALUE), +false);
    SEQAN_ASSERT_EQ((+Or<False, True>::VALUE),  +true);
    SEQAN_ASSERT_EQ((+Or<True, False>::VALUE),  +true);
    SEQAN_ASSERT_EQ((+Or<True, True>::VALUE),   +true);

    // Test for the Type members.
    SEQAN_ASSERT_EQ((TestTypeEq<Or<False, False>::Type, False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<Or<False, True>::Type,  True>::VALUE),  true);
    SEQAN_ASSERT_EQ((TestTypeEq<Or<True, False>::Type,  True>::VALUE),  true);
    SEQAN_ASSERT_EQ((TestTypeEq<Or<True, True>::Type,   True>::VALUE),  true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_logic_or_c)
{
    using namespace seqan;

    // Test for the VALUE members.
    SEQAN_ASSERT_EQ((+OrC<false, false>::VALUE), +false);
    SEQAN_ASSERT_EQ((+OrC<false, true>::VALUE),  +true);
    SEQAN_ASSERT_EQ((+OrC<true, false>::VALUE),  +true);
    SEQAN_ASSERT_EQ((+OrC<true, true>::VALUE),   +true);

    // Test for the Type members.
    SEQAN_ASSERT_EQ((TestTypeEq<OrC<false, false>::Type, False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<OrC<false, true>::Type,  True>::VALUE),  true);
    SEQAN_ASSERT_EQ((TestTypeEq<OrC<true, false>::Type,  True>::VALUE),  true);
    SEQAN_ASSERT_EQ((TestTypeEq<OrC<true, true>::Type,   True>::VALUE),  true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_logic_and)
{
    using namespace seqan;

    // Test for the VALUE members.
    SEQAN_ASSERT_EQ((+And<False, False>::VALUE), +false);
    SEQAN_ASSERT_EQ((+And<False, True>::VALUE),  +false);
    SEQAN_ASSERT_EQ((+And<True, False>::VALUE),  +false);
    SEQAN_ASSERT_EQ((+And<True, True>::VALUE),   +true);

    // Test for the Type members.
    SEQAN_ASSERT_EQ((TestTypeEq<And<False, False>::Type, False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<And<False, True>::Type,  False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<And<True, False>::Type,  False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<And<True, True>::Type,   True>::VALUE),  true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_logic_and_c)
{
    using namespace seqan;

    // Test for the VALUE members.
    SEQAN_ASSERT_EQ((+AndC<false, false>::VALUE), +false);
    SEQAN_ASSERT_EQ((+AndC<false, true>::VALUE),  +false);
    SEQAN_ASSERT_EQ((+AndC<true, false>::VALUE),  +false);
    SEQAN_ASSERT_EQ((+AndC<true, true>::VALUE),   +true);

    // Test for the Type members.
    SEQAN_ASSERT_EQ((TestTypeEq<AndC<false, false>::Type, False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<AndC<false, true>::Type,  False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<AndC<true, false>::Type,  False>::VALUE), true);
    SEQAN_ASSERT_EQ((TestTypeEq<AndC<true, true>::Type,   True>::VALUE),  true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_logic_if)
{
    using namespace seqan;

    // Test for the Type members.
    SEQAN_ASSERT_EQ((TestTypeEq<If<True,  True, False>::Type, True>::VALUE),  true);
    SEQAN_ASSERT_EQ((TestTypeEq<If<False, True, False>::Type, False>::VALUE), true);
}

SEQAN_DEFINE_TEST(test_basic_metaprogramming_logic_if_c)
{
    using namespace seqan;

    // Test for the Type members.
    SEQAN_ASSERT_EQ((TestTypeEq<If<True,  True, False>::Type, True>::VALUE),  true);
    SEQAN_ASSERT_EQ((TestTypeEq<If<False, True, False>::Type, False>::VALUE), true);
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_METAPROGRAMMING_LOGIC_H_
