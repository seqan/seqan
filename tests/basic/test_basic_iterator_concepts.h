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
// Tests for alphabet concepts.
//
// All checks are compile time checks, as long as this file compiles, the
// checks pass.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_CONCEPTS_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_CONCEPTS_H_

// ============================================================================
// Test basic concepts
// ============================================================================

// Test the conformance to the Alphabet concept for (1) built-in types, (2) aggregates, (3) simple types.

inline void testAlphabetConcepts()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((AlphabetConcept<bool>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<char>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<unsigned>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<int>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<double>));

    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Pair<int, double> >));

    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Dna>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Dna5>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<DnaQ>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Dna5Q>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Rna>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Rna5>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Iupac>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<AminoAcid>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<char>));
    // SEQAN_CONCEPT_ASSERT((AlphabetConcept<Finite<10> >));
}

// Test the conformance to the OrderedAlphabet concept for (1) built-in types, (2) aggregates, (3) simple types.

inline void testOrderedAlphabetConcepts()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<bool>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<char>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<unsigned>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<int>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<double>));

    // TODO(holtgrew): With lexicographic ordering, pairs are also ordered. We would need complete implementations of the functions then, though.
    // SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<Pair<int, double> >));

    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<Dna>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<Dna5>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<DnaQ>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<Dna5Q>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<Rna>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<Rna5>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<Iupac>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<AminoAcid>));
    SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<char>));
    // SEQAN_CONCEPT_ASSERT((OrderedAlphabetConcept<Finite<10> >));
}

// Test the conformance to the FiniteOrderedAlphabet concept for (1) built-in types, (2) aggregates, (3) simple types.

inline void testFiniteOrderedAlphabetConcepts()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<bool>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<char>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<unsigned>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<int>));

    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<Dna>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<Dna5>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<DnaQ>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<Dna5Q>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<Rna>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<Rna5>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<Iupac>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<AminoAcid>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<char>));
    // SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<Finite<10> >));
}

// Test the conformance to the AlphabetWithGaps concept for char and modified alphabets.

inline void testAlphabetWithGapsConcept()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((AlphabetWithGapsConcept<char>));
    // // TODO(holtgrew): <seqan/modifier.h> is required for this.
    // SEQAN_CONCEPT_ASSERT((AlphabetWithGapsConcept<ModifiedAlphabet<Dna, ModExpand<'-'> > >));
}

// Test the conformance to the AlphabetWithGaps concept for char and the *5 alphabets.

inline void testAlphabetWithUnknownValueConcept()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((AlphabetWithUnknownValueConcept<char>));
    SEQAN_CONCEPT_ASSERT((AlphabetWithUnknownValueConcept<Dna5>));
    SEQAN_CONCEPT_ASSERT((AlphabetWithUnknownValueConcept<Rna5>));
    SEQAN_CONCEPT_ASSERT((AlphabetWithUnknownValueConcept<Dna5Q>));
}

// Test the conformance to the AlphabetWithQualities concept.

inline void testAlphabetWithQualitiesConcept()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((AlphabetWithQualitiesConcept<DnaQ>));
    SEQAN_CONCEPT_ASSERT((AlphabetWithQualitiesConcept<Dna5Q>));
}

// Test conformance of pointer with TrivialIteratorConcept.
/*
inline void testTrivialIteratorConcept()
{
    SEQAN_CONCEPT_ASSERT((TrivialIteratorConcept<int *>));
}

// Test conformance of pointer with InputIteratorConcept.

inline void testInputIteratorConcept()
{
    SEQAN_CONCEPT_ASSERT((InputIteratorConcept<int *>));
}

// Test conformance of pointer with OutputIteratorConcept.

inline void testOutputIteratorConcept()
{
    SEQAN_CONCEPT_ASSERT((OutputIteratorConcept<int *>));
}

// Test conformance of pointer with ForwardIteratorConcept.

inline void testForwardIteratorConcept()
{
    SEQAN_CONCEPT_ASSERT((ForwardIteratorConcept<int *>));
}

// Test conformance of pointer with BidirectionalIteratorConcept.

inline void testBidirectionalIteratorConcept()
{
    SEQAN_CONCEPT_ASSERT((BidirectionalIteratorConcept<int *>));
}

// Test conformance of pointer with RandomAccessIteratorConcept.

inline void testRandomAccessIteratorConcept()
{
    SEQAN_CONCEPT_ASSERT((RandomAccessIteratorConcept<int *>));
}
*/

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
    testSignedInteger<int64_t>();
    testUnsignedInteger<unsigned char>();
    testUnsignedInteger<unsigned short>();
    testUnsignedInteger<unsigned int>();
    testUnsignedInteger<unsigned long>();
    testUnsignedInteger<uint64_t>();
//  testSignedInteger<unsigned long>();   // <== this should fail
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_CONCEPTS_H_
