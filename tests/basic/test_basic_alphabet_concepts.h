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

// Test the conformance to the Alphabet concept for (1) built-in types, (2) aggregates, (3) simple types.

inline void testAlphabetConcepts()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((AlphabetConcept<bool>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<char>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<signed short>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<short>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<unsigned short>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<unsigned>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<signed int>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<int>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<unsigned int>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<signed long>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<long>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<unsigned long>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<int8_t>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<uint8_t >));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<int16_t>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<uint16_t>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<int32_t>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<uint32_t>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<int64_t>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<uint64_t>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<float>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<double>));

    // TODO(holtgrew): This has to be checked in aggregate module
    // SEQAN_CONCEPT_ASSERT((AlphabetConcept<Pair<int, double> >));

    SEQAN_CONCEPT_ASSERT((AlphabetConcept<char>));
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
}

// Test the conformance to the FiniteOrderedAlphabet concept for (1) built-in types, (2) aggregates, (3) simple types.

inline void testFiniteOrderedAlphabetConcepts()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<bool>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<char>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<unsigned>));
    SEQAN_CONCEPT_ASSERT((FiniteOrderedAlphabetConcept<int>));
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
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_CONCEPTS_H_
