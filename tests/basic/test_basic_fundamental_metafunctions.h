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
// Tests for fundamental metafunctions.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_METAFUNCTIONS_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_METAFUNCTIONS_H_

// Test availability of the symbol Value<T, I>.
SEQAN_DEFINE_TEST(test_basic_fundamental_metafunctions_value)
{
    typedef seqan::Value<int> TIntValue     SEQAN_UNUSED_TYPEDEF;
    typedef seqan::Value<int, 0> TIntValue0 SEQAN_UNUSED_TYPEDEF;
    typedef seqan::Value<int, 1> TIntValue1 SEQAN_UNUSED_TYPEDEF;
    typedef seqan::Value<int, 2> TIntValue2 SEQAN_UNUSED_TYPEDEF;
}

// Test availability of the symbol GetValue<T>.
SEQAN_DEFINE_TEST(test_basic_fundamental_metafunctions_get_value)
{
    typedef seqan::GetValue<int> TIntGetValue SEQAN_UNUSED_TYPEDEF;
}

// Test availability of the symbol Reference<T>.
SEQAN_DEFINE_TEST(test_basic_fundamental_metafunctions_reference)
{
    typedef seqan::Reference<int> TIntReference SEQAN_UNUSED_TYPEDEF;
}

// Test availability of the symbol Size<T>.
SEQAN_DEFINE_TEST(test_basic_fundamental_metafunctions_size)
{
    typedef seqan::Size<int> TIntSize SEQAN_UNUSED_TYPEDEF;
}

// Test availability of the symbol Difference<T>.
SEQAN_DEFINE_TEST(test_basic_fundamental_metafunctions_difference)
{
    typedef seqan::Difference<int> TIntDifference SEQAN_UNUSED_TYPEDEF;
}

// Test availability of the symbol Position<T>.
SEQAN_DEFINE_TEST(test_basic_fundamental_metafunctions_position)
{
    typedef seqan::Position<int> TIntPosition SEQAN_UNUSED_TYPEDEF;
}

// Test the Spec<> metafunction.

// Helper class.
template <typename TSpec>
struct ClassName;

SEQAN_DEFINE_TEST(test_basic_fundamental_metafunctions_spec)
{
    using namespace seqan;

    SEQAN_ASSERT((+SameType_<typename Spec<int>::Type, void>::VALUE));
    SEQAN_ASSERT((+SameType_<typename Spec<int const>::Type, void>::VALUE));
    SEQAN_ASSERT((+SameType_<typename Spec<ClassName<int> >::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename Spec<ClassName<int> const>::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename Spec<ClassName<int const> >::Type, int const>::VALUE));
    SEQAN_ASSERT((+SameType_<typename Spec<ClassName<int const> const>::Type, int const>::VALUE));
}

// Test the DeepestSpec<> metafunction.

// Helper classes.

class Nothing_ {};

SEQAN_DEFINE_TEST(test_basic_fundamental_metafunctions_deepest_spec)
{
    using namespace seqan;

    // No nesting.
    // TODO(holtgrew): Should the result for int not be void?
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<int>::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<int const>::Type, int>::VALUE));

    // On level of nesting.
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<ClassName<int> >::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<ClassName<int> const>::Type, int>::VALUE));
    // TODO(holtgrew): Should the const really be removed?
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<ClassName<int const> >::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<ClassName<int const> const>::Type, int>::VALUE));

    // Two levels of nesting.
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<ClassName<ClassName<int> > >::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<ClassName<ClassName<int> > >::Type, int>::VALUE));

    // Three levels of nesting.
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<ClassName<ClassName<ClassName<int> > > >::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename DeepestSpec<ClassName<ClassName<ClassName<int> > > >::Type, int>::VALUE));
}


#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_METAFUNCTIONS_H_
