// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2021, Knut Reinert, FU Berlin
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
// Tests for modifier/modifier_shortcuts.h
// ==========================================================================

#ifndef SEQAN_TESTS_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_
#define SEQAN_TESTS_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_reverse)
{
    seqan2::DnaString str = "CGAT";
    seqan2::DnaString const EXPECTED_STRING = "ATCG";

    seqan2::DnaStringReverseComplement modifiedString(str);
    seqan2::DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_reverse)
{
    seqan2::Dna5String str = "CGATN";
    seqan2::Dna5String const EXPECTED_STRING = "NATCG";

    seqan2::Dna5StringReverseComplement modifiedString(str);
    seqan2::Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_reverse)
{
    seqan2::RnaString str = "CGAU";
    seqan2::RnaString const EXPECTED_STRING = "AUCG";

    seqan2::RnaStringReverseComplement modifiedString(str);
    seqan2::RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_reverse)
{
    seqan2::Rna5String str = "CGAUN";
    seqan2::Rna5String const EXPECTED_STRING = "NAUCG";

    seqan2::Rna5StringReverseComplement modifiedString(str);
    seqan2::Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_complement)
{
    seqan2::DnaString str = "CGAT";
    seqan2::DnaString const EXPECTED_STRING = "ATCG";

    seqan2::DnaStringReverseComplement modifiedString(str);
    seqan2::DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_complement)
{
    seqan2::Dna5String str = "CGATN";
    seqan2::Dna5String const EXPECTED_STRING = "NATCG";

    seqan2::Dna5StringReverseComplement modifiedString(str);
    seqan2::Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_complement)
{
    seqan2::RnaString str = "CGAU";
    seqan2::RnaString const EXPECTED_STRING = "AUCG";

    seqan2::RnaStringReverseComplement modifiedString(str);
    seqan2::RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_complement)
{
    seqan2::Rna5String str = "CGAUN";
    seqan2::Rna5String const EXPECTED_STRING = "NAUCG";

    seqan2::Rna5StringReverseComplement modifiedString(str);
    seqan2::Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_reverse_complement)
{
    seqan2::DnaString str = "CGAT";
    seqan2::DnaString const EXPECTED_STRING = "ATCG";

    seqan2::DnaStringReverseComplement modifiedString(str);
    seqan2::DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_reverse_complement)
{
    seqan2::Dna5String str = "CGATN";
    seqan2::Dna5String const EXPECTED_STRING = "NATCG";

    seqan2::Dna5StringReverseComplement modifiedString(str);
    seqan2::Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_reverse_complement)
{
    seqan2::RnaString str = "CGAU";
    seqan2::RnaString const EXPECTED_STRING = "AUCG";

    seqan2::RnaStringReverseComplement modifiedString(str);
    seqan2::RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_reverse_complement)
{
    seqan2::Rna5String str = "CGAUN";
    seqan2::Rna5String const EXPECTED_STRING = "NAUCG";

    seqan2::Rna5StringReverseComplement modifiedString(str);
    seqan2::Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_complement_in_place_string)
{
    seqan2::Dna5String str = "CGATN";
    seqan2::Dna5String const STR = str;
    seqan2::Dna5String const EXPECTED_RESULT = "GCTAN";

    // Test non-const version.
    complement(str);
    SEQAN_ASSERT_EQ(EXPECTED_RESULT, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_complement_in_place_string_set)
{
    seqan2::Dna5String str1 = "CCGGTTAANN";
    seqan2::Dna5String str2 = "CGTANCGTAN";
    seqan2::Dna5String const EXPECTED_STRING1 = "GGCCAATTNN";
    seqan2::Dna5String const EXPECTED_STRING2 = "GCATNGCATN";

    seqan2::StringSet<seqan2::Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    // Test non-const version.
    {
        seqan2::StringSet<seqan2::Dna5String> strSetCopy(strSet);

        complement(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_complement_in_place_string)
{
    seqan2::Dna5String str = "CGATN";
    seqan2::Dna5String const kStr = str;
    seqan2::Dna5String const kExpectedResult = "NATCG";

    // Test non-const version.
    reverseComplement(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_complement_in_place_string_set)
{
    seqan2::Dna5String str1 = "CCGGTTAANN";
    seqan2::Dna5String str2 = "CGTANCGTAN";
    seqan2::Dna5String const EXPECTED_STRING1 = "NNTTAACCGG";
    seqan2::Dna5String const EXPECTED_STRING2 = "NTACGNTACG";

    seqan2::StringSet<seqan2::Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    // Test non-const version.
    {
        seqan2::StringSet<seqan2::Dna5String> strSetCopy = strSet;

        reverseComplement(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_in_place_string)
{
    seqan2::Dna5String str = "CGATN";
    seqan2::Dna5String const kStr = str;
    seqan2::Dna5String const kExpectedResult = "NTAGC";

    // Test non-const version.
    reverse(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_in_place_string_set)
{
    seqan2::Dna5String str1 = "CCGGTTAANN";
    seqan2::Dna5String str2 = "CGTANCGTAN";
    seqan2::Dna5String const EXPECTED_STRING1 = "NNAATTGGCC";
    seqan2::Dna5String const EXPECTED_STRING2 = "NATGCNATGC";

    seqan2::StringSet<seqan2::Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    // Test non-const version.
    {
        seqan2::StringSet<seqan2::Dna5String> strSetCopy(strSet);

        reverse(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_lower_in_place_string)
{
    seqan2::CharString str = "This is a test!";
    seqan2::CharString const kStr = str;
    seqan2::CharString const kExpectedResult = "this is a test!";

    // Test non-const version.
    toLower(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_lower_in_place_string_set)
{
    seqan2::CharString str1 = "This is a test!";
    seqan2::CharString str2 = "This is also a test!";
    seqan2::CharString const EXPECTED_STRING1 = "this is a test!";
    seqan2::CharString const EXPECTED_STRING2 = "this is also a test!";

    // Test non-const version.
    seqan2::StringSet<seqan2::CharString> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    toLower(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSet[0]);
    SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSet[1]);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_upper_in_place_string)
{
    seqan2::CharString str = "This is a test!";
    seqan2::CharString const kStr = str;
    seqan2::CharString const kExpectedResult = "THIS IS A TEST!";

    // Test non-const version.
    toUpper(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_upper_in_place_string_set)
{
    seqan2::CharString str1 = "This is a test!";
    seqan2::CharString str2 = "This is also a test!";
    seqan2::CharString const EXPECTED_STRING1 = "THIS IS A TEST!";
    seqan2::CharString const EXPECTED_STRING2 = "THIS IS ALSO A TEST!";

    // Test non-const version.
    seqan2::StringSet<seqan2::CharString> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    toUpper(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSet[0]);
    SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSet[1]);
}

// https://github.com/seqan/seqan/pull/2433
SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_packed_concat_direct_stringset)
{
    // It might require a few repetitions to encounter the error.
    size_t const test_repetitions{4};

    for (size_t repetition = 0; repetition < test_repetitions; ++repetition)
    {
        size_t const string_repetitions{4};
        seqan2::DnaString str1("ACGTAGCTCGTACGATCGATCGTAGCATCGATCG");
        seqan2::DnaString str2("ACTACGATGCTAGCTGACTGAC");
        seqan2::DnaString const EXPECTED_STRING1 = "GCTAGCTACGATGCTAGCTAGCATGCTCGATGCA";
        seqan2::DnaString const EXPECTED_STRING2 = "CAGTCAGTCGATCGTAGCATCA";

        // Test non-const version.
        seqan2::StringSet<seqan2::String<seqan2::Dna, seqan2::Packed<> >, seqan2::Owner<seqan2::ConcatDirect<> > > strSet;
        for (size_t i = 0; i < string_repetitions; ++i)
        {
            appendValue(strSet, str1);
            appendValue(strSet, str2);
        }

        seqan2::reverse(strSet);
        SEQAN_ASSERT_EQ(string_repetitions * 2u, length(strSet));

        for (size_t i = 0; i < string_repetitions; ++i)
        {
            SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSet[i]);
            SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSet[++i]);
        }
    }
}

// https://github.com/seqan/seqan/pull/2433
SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_complement_packed_concat_direct_stringset)
{
    // It might require a few repetitions to encounter the error.
    size_t const test_repetitions{4};

    for (size_t repetition = 0; repetition < test_repetitions; ++repetition)
    {
        size_t const string_repetitions{4};
        seqan2::DnaString str1("ACGTAGCTCGTACGATCGATCGTAGCATCGATCG");
        seqan2::DnaString str2("ACTACGATGCTAGCTGACTGAC");
        seqan2::DnaString const EXPECTED_STRING1 = "CGATCGATGCTACGATCGATCGTACGAGCTACGT";
        seqan2::DnaString const EXPECTED_STRING2 = "GTCAGTCAGCTAGCATCGTAGT";

        // Test non-const version.
        seqan2::StringSet<seqan2::String<seqan2::Dna, seqan2::Packed<> >, seqan2::Owner<seqan2::ConcatDirect<> > > strSet;
        for (size_t i = 0; i < string_repetitions; ++i)
        {
            appendValue(strSet, str1);
            appendValue(strSet, str2);
        }

        seqan2::reverseComplement(strSet, seqan2::Parallel());
        SEQAN_ASSERT_EQ(string_repetitions * 2u, length(strSet));

        for (size_t i = 0; i < string_repetitions; ++i)
        {
            SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSet[i]);
            SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSet[++i]);
        }
    }
}

#endif  // SEQAN_TESTS_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_
