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

#ifndef TESTS_BASIC_TEST_BASIC_CONCEPTS_H_
#define TESTS_BASIC_TEST_BASIC_CONCEPTS_H_


//#include <boost/concept_check.hpp>
//#include <boost/concept/assert.hpp>
//#include <boost/concept/requires.hpp>

struct ConceptsInStruct
{
    SEQAN_CONCEPT_ASSERT((seqan::IntegerConcept<char>));
    // SEQAN_CONCEPT_ASSERT((seqan::IntegerConcept<double>));

    // BOOST_CONCEPT_ASSERT((boost::Integer<char>));
    // BOOST_CONCEPT_ASSERT((boost::Integer<double>));
};

template <typename T>
SEQAN_CONCEPT_REQUIRES(
        ((seqan::IntegerConcept<T>))
        ((seqan::SignedIntegerConcept<T>)),
        (void)) // return type
foo(T const &)
{}

// template <typename T>
// BOOST_CONCEPT_REQUIRES(
//         ((boost::Integer<T>))
//         ((boost::SignedInteger<T>)),
//         (void)) // return type
// foo2(T const &)
// {}

// A test for strings.
SEQAN_DEFINE_TEST(test_basic_concepts_integer_concept)
{
    using namespace seqan;

    // foo(10);
    // foo(10u);
    // foo(3.0);

    // SEQAN_ASSERT((IntegerConcept<char>));
    // SEQAN_ASSERT_EQ((seqan::IntegerConcept<double>::VALUE), 0);

    // foo2(10);
    // foo2(10u);

    SEQAN_CONCEPT_ASSERT((IntegerConcept<char>));
    // SEQAN_CONCEPT_ASSERT((seqan::IntegerConcept<double>));

    // BOOST_CONCEPT_ASSERT((boost::Integer<char>));
    // BOOST_CONCEPT_ASSERT((boost::Integer<double>));
}

SEQAN_DEFINE_TEST(test_basic_concepts_move_construtible_concept)
{
    using namespace seqan;

    struct TestMoveable
    {
        TestMoveable() = default;

        TestMoveable(TestMoveable && other)
        {
            ignoreUnusedVariableWarning(other);
        }
    };

    struct TestNotMoveable
    {
        TestNotMoveable() = default;

        TestNotMoveable(TestNotMoveable const & other)
        {
            ignoreUnusedVariableWarning(other);
        }
    };

    struct TestNotCopyableAndMoveable
    {
        TestNotCopyableAndMoveable() = default;

        TestNotCopyableAndMoveable(TestNotCopyableAndMoveable const &) = delete;
        TestNotCopyableAndMoveable(TestNotCopyableAndMoveable &&) = delete;
    };

    SEQAN_CONCEPT_ASSERT((MoveConstructible<char>));
    SEQAN_CONCEPT_ASSERT((MoveConstructible<TestMoveable>));
    SEQAN_CONCEPT_ASSERT((MoveConstructible<TestNotMoveable>));
    //NOTE: Fails compiling because of deleted copy and move c'tor, which is the expected behavior.
    //SEQAN_CONCEPT_ASSERT((MoveConstructible<TestNotCopyableAndMoveable>));
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_CONCEPTS_H_
