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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for fundamental concepts.
//
// All checks are compile time checks, as long as this file compiles, the
// checks pass.
// ==========================================================================

#ifndef SEQAN_CORE_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_CONCEPTS_H_
#define SEQAN_CORE_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_CONCEPTS_H_

// ============================================================================
// Test Fundamental Concepts
// ============================================================================

template <typename T>
inline void testInteger()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((IntegerConcept<T>));
    SEQAN_STATIC_ASSERT_MSG(Is< IntegerConcept<T> >::VALUE, "Type is not marked to be an integer");
}

template <typename T>
inline void testSignedInteger()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<T>));
    SEQAN_STATIC_ASSERT_MSG(Is< SignedIntegerConcept<T> >::VALUE, "Type is not marked to be a signed integer");
    testInteger<T>();
}

template <typename T>
inline void testUnsignedInteger()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((UnsignedIntegerConcept<T>));
    SEQAN_STATIC_ASSERT_MSG(Is< UnsignedIntegerConcept<T> >::VALUE, "Type is not marked to be an unsigned integer");
    testInteger<T>();
}

inline void testIntegers()
{
    using namespace seqan;

    testInteger<char>();
    testSignedInteger<signed char>();
    testSignedInteger<signed short>();
    testSignedInteger<signed int>();
    testSignedInteger<signed long>();
    testSignedInteger<__int64>();
    testUnsignedInteger<unsigned char>();
    testUnsignedInteger<unsigned short>();
    testUnsignedInteger<unsigned int>();
    testUnsignedInteger<unsigned long>();
    testUnsignedInteger<__uint64>();
//  testSignedInteger<unsigned long>();   // <== this should fail
}

namespace seqan
{

    SEQAN_CONCEPT(TestConceptA_, (T)){};

    SEQAN_CONCEPT(TestConceptB_, (T)){};

    struct TestConceptModelA_ {};

    template <>
    SEQAN_CONCEPT_IMPL((TestConceptA_)(TestConceptB_), TestConceptModelA_);

    template <typename T1, typename T2>
    struct TestConceptModelBSpec_{};

    template <typename T, typename TSpec>
    struct TestConceptModelB_{};

    template <typename T, typename TSpec1, typename TSpec2>
    SEQAN_CONCEPT_IMPL((TestConceptA_)(TestConceptB_), TestConceptModelB_<T, TestConceptModelBSpec_<TSpec1, TSpec2> >);
}

SEQAN_DEFINE_TEST(test_basic_concepts_concept_impl)
{
    using namespace seqan;

    typedef TestConceptModelA_ TModelA;
    typedef TestConceptModelB_<int, TestConceptModelBSpec_<void, void> > TModelB;

    SEQAN_CONCEPT_ASSERT((TestConceptA_<TModelA>));
    SEQAN_CONCEPT_ASSERT((TestConceptB_<TModelA>));

    SEQAN_STATIC_ASSERT_MSG(Is< TestConceptA_<TModelA> >::VALUE, "Type is not marked to be a TestConceptA_");
    SEQAN_STATIC_ASSERT_MSG(Is< TestConceptB_<TModelA> >::VALUE, "Type is not marked to be a TestConceptB_");

    SEQAN_CONCEPT_ASSERT((TestConceptA_<TModelB>));
    SEQAN_CONCEPT_ASSERT((TestConceptB_<TModelB>));

    SEQAN_STATIC_ASSERT_MSG(Is< TestConceptA_<TModelB> >::VALUE, "Type is not marked to be a TestConceptA_");
    SEQAN_STATIC_ASSERT_MSG(Is< TestConceptB_<TModelB> >::VALUE, "Type is not marked to be a TestConceptB_");

}

#endif  // #ifndef SEQAN_CORE_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_CONCEPTS_H_
