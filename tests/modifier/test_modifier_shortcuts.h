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
// Tests for modifier/modifier_shortcuts.h
// ==========================================================================

#ifndef SEQAN_TESTS_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_
#define SEQAN_TESTS_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_reverse)
{
    seqan::DnaString str = "CGAT";
    seqan::DnaString const EXPECTED_STRING = "ATCG";

    seqan::DnaStringReverseComplement modifiedString(str);
    seqan::DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_reverse)
{
    seqan::Dna5String str = "CGATN";
    seqan::Dna5String const EXPECTED_STRING = "NATCG";

    seqan::Dna5StringReverseComplement modifiedString(str);
    seqan::Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_reverse)
{
    seqan::RnaString str = "CGAU";
    seqan::RnaString const EXPECTED_STRING = "AUCG";

    seqan::RnaStringReverseComplement modifiedString(str);
    seqan::RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_reverse)
{
    seqan::Rna5String str = "CGAUN";
    seqan::Rna5String const EXPECTED_STRING = "NAUCG";

    seqan::Rna5StringReverseComplement modifiedString(str);
    seqan::Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_complement)
{
    seqan::DnaString str = "CGAT";
    seqan::DnaString const EXPECTED_STRING = "ATCG";

    seqan::DnaStringReverseComplement modifiedString(str);
    seqan::DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_complement)
{
    seqan::Dna5String str = "CGATN";
    seqan::Dna5String const EXPECTED_STRING = "NATCG";

    seqan::Dna5StringReverseComplement modifiedString(str);
    seqan::Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_complement)
{
    seqan::RnaString str = "CGAU";
    seqan::RnaString const EXPECTED_STRING = "AUCG";

    seqan::RnaStringReverseComplement modifiedString(str);
    seqan::RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_complement)
{
    seqan::Rna5String str = "CGAUN";
    seqan::Rna5String const EXPECTED_STRING = "NAUCG";

    seqan::Rna5StringReverseComplement modifiedString(str);
    seqan::Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_reverse_complement)
{
    seqan::DnaString str = "CGAT";
    seqan::DnaString const EXPECTED_STRING = "ATCG";

    seqan::DnaStringReverseComplement modifiedString(str);
    seqan::DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_reverse_complement)
{
    seqan::Dna5String str = "CGATN";
    seqan::Dna5String const EXPECTED_STRING = "NATCG";

    seqan::Dna5StringReverseComplement modifiedString(str);
    seqan::Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_reverse_complement)
{
    seqan::RnaString str = "CGAU";
    seqan::RnaString const EXPECTED_STRING = "AUCG";

    seqan::RnaStringReverseComplement modifiedString(str);
    seqan::RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_reverse_complement)
{
    seqan::Rna5String str = "CGAUN";
    seqan::Rna5String const EXPECTED_STRING = "NAUCG";

    seqan::Rna5StringReverseComplement modifiedString(str);
    seqan::Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_complement_in_place_string)
{
    seqan::Dna5String str = "CGATN";
    seqan::Dna5String const STR = str;
    seqan::Dna5String const EXPECTED_RESULT = "GCTAN";

    // Test non-const version.
    complement(str);
    SEQAN_ASSERT_EQ(EXPECTED_RESULT, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_complement_in_place_string_set)
{
    seqan::Dna5String str1 = "CCGGTTAANN";
    seqan::Dna5String str2 = "CGTANCGTAN";
    seqan::Dna5String const EXPECTED_STRING1 = "GGCCAATTNN";
    seqan::Dna5String const EXPECTED_STRING2 = "GCATNGCATN";

    seqan::StringSet<seqan::Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    // Test non-const version.
    {
        seqan::StringSet<seqan::Dna5String> strSetCopy(strSet);

        complement(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_complement_in_place_string)
{
    seqan::Dna5String str = "CGATN";
    seqan::Dna5String const kStr = str;
    seqan::Dna5String const kExpectedResult = "NATCG";

    // Test non-const version.
    reverseComplement(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_complement_in_place_string_set)
{
    seqan::Dna5String str1 = "CCGGTTAANN";
    seqan::Dna5String str2 = "CGTANCGTAN";
    seqan::Dna5String const EXPECTED_STRING1 = "NNTTAACCGG";
    seqan::Dna5String const EXPECTED_STRING2 = "NTACGNTACG";

    seqan::StringSet<seqan::Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    // Test non-const version.
    {
        seqan::StringSet<seqan::Dna5String> strSetCopy = strSet;

        reverseComplement(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_in_place_string)
{
    seqan::Dna5String str = "CGATN";
    seqan::Dna5String const kStr = str;
    seqan::Dna5String const kExpectedResult = "NTAGC";

    // Test non-const version.
    reverse(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_in_place_string_set)
{
    seqan::Dna5String str1 = "CCGGTTAANN";
    seqan::Dna5String str2 = "CGTANCGTAN";
    seqan::Dna5String const EXPECTED_STRING1 = "NNAATTGGCC";
    seqan::Dna5String const EXPECTED_STRING2 = "NATGCNATGC";

    seqan::StringSet<seqan::Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    // Test non-const version.
    {
        seqan::StringSet<seqan::Dna5String> strSetCopy(strSet);

        reverse(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_lower_in_place_string)
{
    seqan::CharString str = "This is a test!";
    seqan::CharString const kStr = str;
    seqan::CharString const kExpectedResult = "this is a test!";

    // Test non-const version.
    toLower(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_lower_in_place_string_set)
{
    seqan::CharString str1 = "This is a test!";
    seqan::CharString str2 = "This is also a test!";
    seqan::CharString const EXPECTED_STRING1 = "this is a test!";
    seqan::CharString const EXPECTED_STRING2 = "this is also a test!";

    // Test non-const version.
    seqan::StringSet<seqan::CharString> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    toLower(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSet[0]);
    SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSet[1]);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_upper_in_place_string)
{
    seqan::CharString str = "This is a test!";
    seqan::CharString const kStr = str;
    seqan::CharString const kExpectedResult = "THIS IS A TEST!";

    // Test non-const version.
    toUpper(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);
}

SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_upper_in_place_string_set)
{
    seqan::CharString str1 = "This is a test!";
    seqan::CharString str2 = "This is also a test!";
    seqan::CharString const EXPECTED_STRING1 = "THIS IS A TEST!";
    seqan::CharString const EXPECTED_STRING2 = "THIS IS ALSO A TEST!";

    // Test non-const version.
    seqan::StringSet<seqan::CharString> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    toUpper(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSet[0]);
    SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSet[1]);
}

#endif  // SEQAN_TESTS_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_
