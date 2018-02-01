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

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_RESIDUE_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_ALPHABET_RESIDUE_H_

// TODO(holtgrew): Tests for Finite.

#include <sstream>

// --------------------------------------------------------------------------
// Check Concept Conformance
// --------------------------------------------------------------------------

inline void testResidueConcepts()
{
    using namespace seqan;

    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Dna>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Dna5>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<DnaQ>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Dna5Q>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Rna>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Rna5>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Iupac>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<AminoAcid>));
    // SEQAN_CONCEPT_ASSERT((AlphabetConcept<Finite<10> >));

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

    SEQAN_CONCEPT_ASSERT((AlphabetWithUnknownValueConcept<Dna5>));
    SEQAN_CONCEPT_ASSERT((AlphabetWithUnknownValueConcept<Rna5>));
    SEQAN_CONCEPT_ASSERT((AlphabetWithUnknownValueConcept<Dna5Q>));

    SEQAN_CONCEPT_ASSERT((AlphabetWithQualitiesConcept<DnaQ>));
    SEQAN_CONCEPT_ASSERT((AlphabetWithQualitiesConcept<Dna5Q>));
}

// --------------------------------------------------------------------------
// Check Metafunction and Function Implementation
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_metafunctions_dna)
{
    using namespace seqan;

    // Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+(SameType_<typename BitsPerValue<Dna>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+(BitsPerValue<Dna>::VALUE), 2);

    // Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(std::numeric_limits<Dna>::max(), Dna('T'));
    SEQAN_ASSERT_EQ(std::numeric_limits<Dna>::min(), Dna('A'));

    // Finited Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(ordValue(Dna('A')), 0);
    SEQAN_ASSERT_EQ(ordValue(Dna('C')), 1);
    SEQAN_ASSERT_EQ(ordValue(Dna('G')), 2);
    SEQAN_ASSERT_EQ(ordValue(Dna('T')), 3);

    SEQAN_ASSERT(+(SameType_<typename ValueSize<Dna>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+ValueSize<Dna>::VALUE, 4);
    SEQAN_ASSERT_EQ(valueSize<Dna>(), 4u);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_metafunctions_dna5)
{
    using namespace seqan;

    // Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+(SameType_<typename BitsPerValue<Dna5>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+(BitsPerValue<Dna5>::VALUE), 3);

    // Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(std::numeric_limits<Dna5>::max(), Dna5('N'));
    SEQAN_ASSERT_EQ(std::numeric_limits<Dna5>::min(), Dna5('A'));

    // Finited Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(ordValue(Dna5('A')), 0);
    SEQAN_ASSERT_EQ(ordValue(Dna5('C')), 1);
    SEQAN_ASSERT_EQ(ordValue(Dna5('G')), 2);
    SEQAN_ASSERT_EQ(ordValue(Dna5('T')), 3);
    SEQAN_ASSERT_EQ(ordValue(Dna5('N')), 4);

    SEQAN_ASSERT(+(SameType_<typename ValueSize<Dna5>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+ValueSize<Dna5>::VALUE, 5);
    SEQAN_ASSERT_EQ(valueSize<Dna5>(), 5u);

    // Alphabet With Unknown Value Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(unknownValue<Dna5>(), Dna5('N'));
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_metafunctions_dna_q)
{
    using namespace seqan;

    // Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+(SameType_<typename BitsPerValue<DnaQ>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+(BitsPerValue<DnaQ>::VALUE), 8);

    // Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(std::numeric_limits<DnaQ>::max(), DnaQ('T'));
    SEQAN_ASSERT_EQ(std::numeric_limits<DnaQ>::min(), DnaQ('A'));

    // Finited Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(ordValue(DnaQ('A')), 0);
    SEQAN_ASSERT_EQ(ordValue(DnaQ('C')), 1);
    SEQAN_ASSERT_EQ(ordValue(DnaQ('G')), 2);
    SEQAN_ASSERT_EQ(ordValue(DnaQ('T')), 3);

    SEQAN_ASSERT(+(SameType_<typename ValueSize<DnaQ>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+ValueSize<DnaQ>::VALUE, 4);
    SEQAN_ASSERT_EQ(valueSize<DnaQ>(), 4u);

    // Alphabet With Qualities Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+HasQualities<DnaQ>::VALUE);
    SEQAN_ASSERT_EQ(+QualityValueSize<DnaQ>::VALUE, 63);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_metafunctions_dna5_q)
{
    using namespace seqan;

    // Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+(SameType_<typename BitsPerValue<Dna5Q>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+(BitsPerValue<Dna5Q>::VALUE), 8);

    // Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(std::numeric_limits<Dna5Q>::max(), Dna5Q('N'));
    SEQAN_ASSERT_EQ(std::numeric_limits<Dna5Q>::min(), Dna5Q('A'));

    // Finited Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(ordValue(Dna5Q('A')), 0);
    SEQAN_ASSERT_EQ(ordValue(Dna5Q('C')), 1);
    SEQAN_ASSERT_EQ(ordValue(Dna5Q('G')), 2);
    SEQAN_ASSERT_EQ(ordValue(Dna5Q('T')), 3);
    SEQAN_ASSERT_EQ(ordValue(Dna5Q('N')), 4);

    SEQAN_ASSERT(+(SameType_<typename ValueSize<Dna5Q>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+ValueSize<Dna5Q>::VALUE, 5);
    SEQAN_ASSERT_EQ(valueSize<Dna5Q>(), 5u);

    // Alphabet With Unknown Value Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(unknownValue<Dna5Q>(), Dna5Q('N'));

    // Alphabet With Qualities Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+HasQualities<Dna5Q>::VALUE);
    SEQAN_ASSERT_EQ(+QualityValueSize<Dna5Q>::VALUE, 63);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_metafunctions_rna)
{
    using namespace seqan;

    // Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+(SameType_<typename BitsPerValue<Rna>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+(BitsPerValue<Rna>::VALUE), 2);

    // Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(std::numeric_limits<Rna>::max(), Rna('U'));
    SEQAN_ASSERT_EQ(std::numeric_limits<Rna>::min(), Rna('A'));

    // Finited Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(ordValue(Rna('A')), 0);
    SEQAN_ASSERT_EQ(ordValue(Rna('C')), 1);
    SEQAN_ASSERT_EQ(ordValue(Rna('G')), 2);
    SEQAN_ASSERT_EQ(ordValue(Rna('U')), 3);

    SEQAN_ASSERT(+(SameType_<typename ValueSize<Rna>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+ValueSize<Rna>::VALUE, 4);
    SEQAN_ASSERT_EQ(valueSize<Rna>(), 4u);
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_metafunctions_rna5)
{
    using namespace seqan;

    // Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+(SameType_<typename BitsPerValue<Rna5>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+(BitsPerValue<Rna5>::VALUE), 3);

    // Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(std::numeric_limits<Rna5>::max(), Rna5('N'));
    SEQAN_ASSERT_EQ(std::numeric_limits<Rna5>::min(), Rna5('A'));

    // Finited Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(ordValue(Rna5('A')), 0);
    SEQAN_ASSERT_EQ(ordValue(Rna5('C')), 1);
    SEQAN_ASSERT_EQ(ordValue(Rna5('G')), 2);
    SEQAN_ASSERT_EQ(ordValue(Rna5('U')), 3);
    SEQAN_ASSERT_EQ(ordValue(Rna5('N')), 4);

    SEQAN_ASSERT(+(SameType_<typename ValueSize<Rna5>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+ValueSize<Rna5>::VALUE, 5);
    SEQAN_ASSERT_EQ(valueSize<Rna5>(), 5u);

    // Alphabet With Unknown Value Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(unknownValue<Rna5>(), Rna5('N'));
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_metafunctions_iupac)
{
    using namespace seqan;

    // Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+(SameType_<typename BitsPerValue<Iupac>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+(BitsPerValue<Iupac>::VALUE), 4);

    // Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(std::numeric_limits<Iupac>::max(), Iupac('N'));
    SEQAN_ASSERT_EQ(std::numeric_limits<Iupac>::min(), Iupac('U'));

    // Finited Ordered Alphabet Concept Metafunctions / Type Queries

    for (unsigned i = 0; i < 5; ++i)
    {
        char x = Dna5(i);
        int o = (i < 4)? 1 << i : 15;
        SEQAN_ASSERT_EQ(ordValue(Iupac(x)), o);
        SEQAN_ASSERT_EQ(ordValue(Iupac(Dna5(i))), o);
        SEQAN_ASSERT_EQ((char)(Iupac(o)), (char)x);
        SEQAN_ASSERT_EQ((Dna5)(Iupac(o)), x);
    }

    for (unsigned i = 0; i < 16; ++i)
        SEQAN_ASSERT_EQ(ordValue((Iupac)(char)(Iupac(i))), i);

    SEQAN_ASSERT(+(SameType_<typename ValueSize<Iupac>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+ValueSize<Iupac>::VALUE, 16);
    SEQAN_ASSERT_EQ(valueSize<Iupac>(), 16u);

    // Alphabet With Unknown Value Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(unknownValue<Iupac>(), Iupac('N'));
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_metafunctions_amino_acid)
{
    using namespace seqan;

    // Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT(+(SameType_<typename BitsPerValue<AminoAcid>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+(BitsPerValue<AminoAcid>::VALUE), 5);

    // Ordered Alphabet Concept Metafunctions / Type Queries
    SEQAN_ASSERT_EQ(std::numeric_limits<AminoAcid>::max(), 26u);
    SEQAN_ASSERT_EQ(std::numeric_limits<AminoAcid>::max(), AminoAcid('*'));
    SEQAN_ASSERT_EQ(std::numeric_limits<AminoAcid>::min(), AminoAcid('A'));

    // Finited Ordered Alphabet Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(ordValue(AminoAcid('A')), 0);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('B')), 1);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('C')), 2);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('D')), 3);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('E')), 4);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('F')), 5);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('G')), 6);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('H')), 7);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('I')), 8);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('J')), 9);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('K')), 10);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('L')), 11);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('M')), 12);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('N')), 13);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('O')), 14);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('P')), 15);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('Q')), 16);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('R')), 17);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('S')), 18);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('T')), 19);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('U')), 20);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('V')), 21);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('W')), 22);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('Y')), 23);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('Z')), 24);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('X')), 25);
    SEQAN_ASSERT_EQ(ordValue(AminoAcid('*')), 26);

    SEQAN_ASSERT(+(SameType_<typename ValueSize<AminoAcid>::Type, uint8_t>::VALUE));
    SEQAN_ASSERT_EQ(+ValueSize<AminoAcid>::VALUE, 27);
    SEQAN_ASSERT_EQ(valueSize<AminoAcid>(), 27u);

    // Alphabet With Unknown Value Concept Metafunctions / Type Queries

    SEQAN_ASSERT_EQ(unknownValue<AminoAcid>(), AminoAcid('X'));
}

// --------------------------------------------------------------------------
// Check Constructors
// --------------------------------------------------------------------------

// constructors: default, copy same type, copy through assign different type
// assignment: same type, different type
// conversion operators: to builtin types.

// --------------------------------------------------------------------------
// Check Assignment
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Check Conversion
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Check Quality Related Functions
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Check Remaining Functionality
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_usage_dna)
{
    using namespace seqan;

    // Test Increment / Decrement Operators
    {
        Dna a = 'A', t = 'T';

        SEQAN_ASSERT_EQ(++a, Dna('C'));
        SEQAN_ASSERT_EQ(a, Dna('C'));
        SEQAN_ASSERT_EQ(a++, Dna('C'));
        SEQAN_ASSERT_EQ(a, Dna('G'));

        SEQAN_ASSERT_EQ(--t, Dna('G'));
        SEQAN_ASSERT_EQ(t, Dna('G'));
        SEQAN_ASSERT_EQ(t--, Dna('G'));
        SEQAN_ASSERT_EQ(t, Dna('C'));
    }

    // Test Relations, Same Type
    {
        Dna c = 'C', c2 = 'C', t = 'T';

        SEQAN_ASSERT(c == c2);
        SEQAN_ASSERT_NOT(c == t);

        SEQAN_ASSERT(c != t);
        SEQAN_ASSERT_NOT(c != c2);

        SEQAN_ASSERT_NOT(c < c2);
        SEQAN_ASSERT(c < t);

        SEQAN_ASSERT(c <= c2);
        SEQAN_ASSERT(c <= t);
        SEQAN_ASSERT_NOT(t <= c);

        SEQAN_ASSERT_NOT(c > c2);
        SEQAN_ASSERT(t > c);

        SEQAN_ASSERT(c >= c2);
        SEQAN_ASSERT(t >= c);
        SEQAN_ASSERT_NOT(c >= t);
    }

    // Test Stream Operator.
    {
        std::stringstream ss;
        ss << Dna('A') << Dna('C') << Dna('G') << Dna('T');
        SEQAN_ASSERT_EQ(ss.str(), "ACGT");
    }
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_usage_dna_q)
{
    using namespace seqan;

    // Test Increment / Decrement Operators
    {
        DnaQ a = 'A', t = 'T';

        SEQAN_ASSERT_EQ(++a, DnaQ('C'));
        SEQAN_ASSERT_EQ(a, DnaQ('C'));
        SEQAN_ASSERT_EQ(a++, DnaQ('C'));
        SEQAN_ASSERT_EQ(a, DnaQ('G'));

        SEQAN_ASSERT_EQ(--t, DnaQ('G'));
        SEQAN_ASSERT_EQ(t, DnaQ('G'));
        SEQAN_ASSERT_EQ(t--, DnaQ('G'));
        SEQAN_ASSERT_EQ(t, DnaQ('C'));
    }

    // Test Relations, Same Type
    {
        DnaQ c = 'C', c2 = 'C', t = 'T';

        SEQAN_ASSERT(c == c2);
        SEQAN_ASSERT_NOT(c == t);

        SEQAN_ASSERT(c != t);
        SEQAN_ASSERT_NOT(c != c2);

        SEQAN_ASSERT_NOT(c < c2);
        SEQAN_ASSERT(c < t);

        SEQAN_ASSERT(c <= c2);
        SEQAN_ASSERT(c <= t);
        SEQAN_ASSERT_NOT(t <= c);

        SEQAN_ASSERT_NOT(c > c2);
        SEQAN_ASSERT(t > c);

        SEQAN_ASSERT(c >= c2);
        SEQAN_ASSERT(t >= c);
        SEQAN_ASSERT_NOT(c >= t);
    }

    // Test Stream Operator.
    {
        std::stringstream ss;
        ss << DnaQ('A') << DnaQ('C') << DnaQ('G') << DnaQ('T');
        SEQAN_ASSERT_EQ(ss.str(), "ACGT");
    }
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_usage_dna5)
{
    using namespace seqan;

    // Test Increment / Decrement Operators
    {
        Dna5 a = 'A', n = 'N';

        SEQAN_ASSERT_EQ(++a, Dna5('C'));
        SEQAN_ASSERT_EQ(a, Dna5('C'));
        SEQAN_ASSERT_EQ(a++, Dna5('C'));
        SEQAN_ASSERT_EQ(a, Dna5('G'));

        SEQAN_ASSERT_EQ(--n, Dna5('T'));
        SEQAN_ASSERT_EQ(n, Dna5('T'));
        SEQAN_ASSERT_EQ(n--, Dna5('T'));
        SEQAN_ASSERT_EQ(n, Dna5('G'));
    }

    // Test Relations, Same Type
    {
        Dna5 c = 'C', c2 = 'C', t = 'T';

        SEQAN_ASSERT(c == c2);
        SEQAN_ASSERT_NOT(c == t);

        SEQAN_ASSERT(c != t);
        SEQAN_ASSERT_NOT(c != c2);

        SEQAN_ASSERT_NOT(c < c2);
        SEQAN_ASSERT(c < t);

        SEQAN_ASSERT(c <= c2);
        SEQAN_ASSERT(c <= t);
        SEQAN_ASSERT_NOT(t <= c);

        SEQAN_ASSERT_NOT(c > c2);
        SEQAN_ASSERT(t > c);

        SEQAN_ASSERT(c >= c2);
        SEQAN_ASSERT(t >= c);
        SEQAN_ASSERT_NOT(c >= t);
    }

    // Test Stream Operator.
    {
        std::stringstream ss;
        ss << Dna5('A') << Dna5('C') << Dna5('G') << Dna5('T') << Dna5('N');
        SEQAN_ASSERT_EQ(ss.str(), "ACGTN");
    }
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_usage_dna5_q)
{
    using namespace seqan;

    // Test Increment / Decrement Operators
    {
        Dna5Q a = 'A', n = 'N';

        SEQAN_ASSERT_EQ(++a, Dna5Q('C'));
        SEQAN_ASSERT_EQ(a, Dna5Q('C'));
        SEQAN_ASSERT_EQ(a++, Dna5Q('C'));
        SEQAN_ASSERT_EQ(a, Dna5Q('G'));

        SEQAN_ASSERT_EQ(--n, Dna5Q('T'));
        SEQAN_ASSERT_EQ(n, Dna5Q('T'));
        SEQAN_ASSERT_EQ(n--, Dna5Q('T'));
        SEQAN_ASSERT_EQ(n, Dna5Q('G'));
    }

    // Test Relations, Same Type
    {
        Dna5Q c = 'C', c2 = 'C', t = 'T';

        SEQAN_ASSERT(c == c2);
        SEQAN_ASSERT_NOT(c == t);

        SEQAN_ASSERT(c != t);
        SEQAN_ASSERT_NOT(c != c2);

        SEQAN_ASSERT_NOT(c < c2);
        SEQAN_ASSERT(c < t);

        SEQAN_ASSERT(c <= c2);
        SEQAN_ASSERT(c <= t);
        SEQAN_ASSERT_NOT(t <= c);

        SEQAN_ASSERT_NOT(c > c2);
        SEQAN_ASSERT(t > c);

        SEQAN_ASSERT(c >= c2);
        SEQAN_ASSERT(t >= c);
        SEQAN_ASSERT_NOT(c >= t);
    }

    // Test Stream Operator.
    {
        std::stringstream ss;
        ss << Dna5Q('A') << Dna5Q('C') << Dna5Q('G') << Dna5Q('T') << Dna5Q('N');
        SEQAN_ASSERT_EQ(ss.str(), "ACGTN");
    }
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_usage_rna)
{
    using namespace seqan;

    // Test Increment / Decrement Operators
    {
        Rna a = 'A', u = 'U';

        SEQAN_ASSERT_EQ(++a, Rna('C'));
        SEQAN_ASSERT_EQ(a, Rna('C'));
        SEQAN_ASSERT_EQ(a++, Rna('C'));
        SEQAN_ASSERT_EQ(a, Rna('G'));

        SEQAN_ASSERT_EQ(--u, Rna('G'));
        SEQAN_ASSERT_EQ(u, Rna('G'));
        SEQAN_ASSERT_EQ(u--, Rna('G'));
        SEQAN_ASSERT_EQ(u, Rna('C'));
    }

    // Test Relations, Same Type
    {
        Rna c = 'C', c2 = 'C', u = 'u';

        SEQAN_ASSERT(c == c2);
        SEQAN_ASSERT_NOT(c == u);

        SEQAN_ASSERT(c != u);
        SEQAN_ASSERT_NOT(c != c2);

        SEQAN_ASSERT_NOT(c < c2);
        SEQAN_ASSERT(c < u);

        SEQAN_ASSERT(c <= c2);
        SEQAN_ASSERT(c <= u);
        SEQAN_ASSERT_NOT(u <= c);

        SEQAN_ASSERT_NOT(c > c2);
        SEQAN_ASSERT(u > c);

        SEQAN_ASSERT(c >= c2);
        SEQAN_ASSERT(u >= c);
        SEQAN_ASSERT_NOT(c >= u);
    }

    // Test Stream Operator.
    {
        std::stringstream ss;
        ss << Rna('A') << Rna('C') << Rna('G') << Rna('U');
        SEQAN_ASSERT_EQ(ss.str(), "ACGU");
    }
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_usage_rna5)
{
    using namespace seqan;

    // Test Increment / Decrement Operators
    {
        Rna5 a = 'A', n = 'N';

        SEQAN_ASSERT_EQ(++a, Rna5('C'));
        SEQAN_ASSERT_EQ(a, Rna5('C'));
        SEQAN_ASSERT_EQ(a++, Rna5('C'));
        SEQAN_ASSERT_EQ(a, Rna5('G'));

        SEQAN_ASSERT_EQ(--n, Rna5('U'));
        SEQAN_ASSERT_EQ(n, Rna5('U'));
        SEQAN_ASSERT_EQ(n--, Rna5('U'));
        SEQAN_ASSERT_EQ(n, Rna5('G'));
    }

    // Test Relations, Same Type
    {
        Rna5 c = 'C', c2 = 'C', u = 'U';

        SEQAN_ASSERT(c == c2);
        SEQAN_ASSERT_NOT(c == u);

        SEQAN_ASSERT(c != u);
        SEQAN_ASSERT_NOT(c != c2);

        SEQAN_ASSERT_NOT(c < c2);
        SEQAN_ASSERT(c < u);

        SEQAN_ASSERT(c <= c2);
        SEQAN_ASSERT(c <= u);
        SEQAN_ASSERT_NOT(u <= c);

        SEQAN_ASSERT_NOT(c > c2);
        SEQAN_ASSERT(u > c);

        SEQAN_ASSERT(c >= c2);
        SEQAN_ASSERT(u >= c);
        SEQAN_ASSERT_NOT(c >= u);
    }

    // Test Stream Operator.
    {
        std::stringstream ss;
        ss << Rna5('A') << Rna5('C') << Rna5('G') << Rna5('U') << Rna5('N');
        SEQAN_ASSERT_EQ(ss.str(), "ACGUN");
    }
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_usage_iupac)
{
    using namespace seqan;

    // Test Increment / Decrement Operators
    {
        Iupac a = 'A', v = 'V';

        SEQAN_ASSERT_EQ(++a, Iupac('C'));
        SEQAN_ASSERT_EQ(a, Iupac('C'));
        SEQAN_ASSERT_EQ(a++, Iupac('C'));
        SEQAN_ASSERT_EQ(a, Iupac('M'));

        SEQAN_ASSERT_EQ(--v, Iupac('S'));
        SEQAN_ASSERT_EQ(v, Iupac('S'));
        SEQAN_ASSERT_EQ(v--, Iupac('S'));
        SEQAN_ASSERT_EQ(v, Iupac('R'));
    }

    // Test Relations, Same Type
    {
        Iupac c = 'C', c2 = 'C', v = 'V';

        SEQAN_ASSERT(c == c2);
        SEQAN_ASSERT_NOT(c == v);

        SEQAN_ASSERT(c != v);
        SEQAN_ASSERT_NOT(c != c2);

        SEQAN_ASSERT_NOT(c < c2);
        SEQAN_ASSERT(c < v);

        SEQAN_ASSERT(c <= c2);
        SEQAN_ASSERT(c <= v);
        SEQAN_ASSERT_NOT(v <= c);

        SEQAN_ASSERT_NOT(c > c2);
        SEQAN_ASSERT(v > c);

        SEQAN_ASSERT(c >= c2);
        SEQAN_ASSERT(v >= c);
        SEQAN_ASSERT_NOT(c >= v);
    }

    // Test Stream Operator.
    {
        std::stringstream ss;
        ss << Iupac('A') << Iupac('C') << Iupac('G') << Iupac('U') << Iupac('N');
        SEQAN_ASSERT_EQ(ss.str(), "ACGUN");
    }
    {
        std::stringstream ss;
        ss << Iupac('A') << Iupac('C') << Iupac('G') << Iupac('=') << Iupac('N');
        SEQAN_ASSERT_EQ(ss.str(), "ACGUN");
    }
}

SEQAN_DEFINE_TEST(test_basic_alphabet_residue_usage_amino_acid)
{
    using namespace seqan;

    // Test Increment / Decrement Operators
    {
        AminoAcid a = 'A', n = 'N';

        SEQAN_ASSERT_EQ(++a, AminoAcid('B'));
        SEQAN_ASSERT_EQ(a, AminoAcid('B'));
        SEQAN_ASSERT_EQ(a++, AminoAcid('B'));
        SEQAN_ASSERT_EQ(a, AminoAcid('C'));

        SEQAN_ASSERT_EQ(--n, AminoAcid('M'));
        SEQAN_ASSERT_EQ(n, AminoAcid('M'));
        SEQAN_ASSERT_EQ(n--, AminoAcid('M'));
        SEQAN_ASSERT_EQ(n, AminoAcid('L'));
    }

    // Test Relations, Same Type
    {
        AminoAcid c = 'C', c2 = 'C', g = 'g';

        SEQAN_ASSERT(c == c2);
        SEQAN_ASSERT_NOT(c == g);

        SEQAN_ASSERT(c != g);
        SEQAN_ASSERT_NOT(c != c2);

        SEQAN_ASSERT_NOT(c < c2);
        SEQAN_ASSERT(c < g);

        SEQAN_ASSERT(c <= c2);
        SEQAN_ASSERT(c <= g);
        SEQAN_ASSERT_NOT(g <= c);

        SEQAN_ASSERT_NOT(c > c2);
        SEQAN_ASSERT(g > c);

        SEQAN_ASSERT(c >= c2);
        SEQAN_ASSERT(g >= c);
        SEQAN_ASSERT_NOT(c >= g);
    }

    // Test Stream Operator.
    {
        std::stringstream ss;
        ss << AminoAcid('A') << AminoAcid('C') << AminoAcid('G') << AminoAcid('X') << AminoAcid('N');
        SEQAN_ASSERT_EQ(ss.str(), "ACGXN");
    }
}

#endif  // #ifndef SEQAN_TESTS_BASIC_ALPHABET_RESIDUE_TYPE_H_
