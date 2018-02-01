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

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_TEST_BASIC_ALPHABET_PROFILE_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_TEST_BASIC_ALPHABET_PROFILE_H_

// --------------------------------------------------------------------------
// Check Concept Conformance
// --------------------------------------------------------------------------

// Not called, checked at compile time.

void testProfileCharConcepts()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((AlphabetConcept<ProfileChar<Dna> >));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<ProfileChar<Dna5> >));
}

// --------------------------------------------------------------------------
// Check Metafunction and Function Implementation
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_alphabet_profile_metafunctions)
{
    using namespace seqan;

    typedef ProfileChar<bool> TProfileChar;

    // Metafunction BitsPerValue

    SEQAN_ASSERT(+(SameType_<typename BitsPerValue<TProfileChar>::Type, unsigned>::VALUE));
    SEQAN_ASSERT_EQ(+(BitsPerValue<TProfileChar>::VALUE), 96u);

    // Metafunction ValueSize

    SEQAN_ASSERT(+(SameType_<typename ValueSize<TProfileChar>::Type, unsigned>::VALUE));
    SEQAN_ASSERT_EQ(+(ValueSize<TProfileChar>::VALUE), 3);

    // Metafunction SourceValue

    SEQAN_ASSERT(+(SameType_<typename SourceValue<TProfileChar>::Type, bool>::VALUE));
    SEQAN_ASSERT(+(SameType_<typename SourceValue<TProfileChar const>::Type, bool>::VALUE));
}

SEQAN_DEFINE_TEST(test_basic_alphabet_profile_constructors)
{
    using namespace seqan;

    typedef ProfileChar<Dna> TProfileChar;

    // Test default constructor.
    {
        TProfileChar c;

        SEQAN_ASSERT_EQ(c.count[0], 0u);
        SEQAN_ASSERT_EQ(c.count[1], 0u);
        SEQAN_ASSERT_EQ(c.count[2], 0u);
        SEQAN_ASSERT_EQ(c.count[3], 0u);
    }

    // Test copy constructor.
    {
        TProfileChar c;
        c.count[2] = 32;

        TProfileChar d(c);

        SEQAN_ASSERT_EQ(d.count[0], 0u);
        SEQAN_ASSERT_EQ(d.count[1], 0u);
        SEQAN_ASSERT_EQ(d.count[2], 32u);
        SEQAN_ASSERT_EQ(d.count[3], 0u);
    }

    // Test construction from source value.
    {
        Dna d = 'C';
        TProfileChar c(d);

        SEQAN_ASSERT_EQ(c.count[0], 0u);
        SEQAN_ASSERT_EQ(c.count[1], 1u);
        SEQAN_ASSERT_EQ(c.count[2], 0u);
        SEQAN_ASSERT_EQ(c.count[3], 0u);
    }
}

SEQAN_DEFINE_TEST(test_basic_alphabet_profile_relations)
{
    using namespace seqan;

    ProfileChar<Dna> c, d;

    SEQAN_ASSERT(c == d);
    SEQAN_ASSERT_NOT(c != d);

    c.count[1] = 3;
    SEQAN_ASSERT(c != d);
    SEQAN_ASSERT_NOT(c == d);

    d.count[1] = 3;
    SEQAN_ASSERT(c == d);
    SEQAN_ASSERT_NOT(c != d);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_profile_empty)
{
    using namespace seqan;

    ProfileChar<Dna> c;

    SEQAN_ASSERT(empty(c));

    c.count[0] = 3;
    SEQAN_ASSERT_NOT(empty(c));
}

// TODO(holtgrew): Test assign(), convertImpl(), operator<<, conversions.

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_PROFILE_H_
