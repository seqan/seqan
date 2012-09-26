// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef TEST_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_
#define TEST_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>

using namespace seqan;

// TODO(holtgrew): The usage of the prefix "k" to indicate const-ness is obsolete and should be replaced by ALL_UPPERCASE.

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_reverse) {
    DnaString str = "CGAT";
    DnaString const EXPECTED_STRING = "ATCG";

    DnaStringReverseComplement modifiedString(str);
    DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_reverse) {
    Dna5String str = "CGATN";
    Dna5String const EXPECTED_STRING = "NATCG";

    Dna5StringReverseComplement modifiedString(str);
    Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_reverse) {
    RnaString str = "CGAU";
    RnaString const EXPECTED_STRING = "AUCG";

    RnaStringReverseComplement modifiedString(str);
    RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_reverse) {
    Rna5String str = "CGAUN";
    Rna5String const EXPECTED_STRING = "NAUCG";

    Rna5StringReverseComplement modifiedString(str);
    Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_complement) {
    DnaString str = "CGAT";
    DnaString const EXPECTED_STRING = "ATCG";

    DnaStringReverseComplement modifiedString(str);
    DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_complement) {
    Dna5String str = "CGATN";
    Dna5String const EXPECTED_STRING = "NATCG";

    Dna5StringReverseComplement modifiedString(str);
    Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_complement) {
    RnaString str = "CGAU";
    RnaString const EXPECTED_STRING = "AUCG";

    RnaStringReverseComplement modifiedString(str);
    RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_complement) {
    Rna5String str = "CGAUN";
    Rna5String const EXPECTED_STRING = "NAUCG";

    Rna5StringReverseComplement modifiedString(str);
    Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_reverse_complement) {
    DnaString str = "CGAT";
    DnaString const EXPECTED_STRING = "ATCG";

    DnaStringReverseComplement modifiedString(str);
    DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_reverse_complement) {
    Dna5String str = "CGATN";
    Dna5String const EXPECTED_STRING = "NATCG";

    Dna5StringReverseComplement modifiedString(str);
    Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_reverse_complement) {
    RnaString str = "CGAU";
    RnaString const EXPECTED_STRING = "AUCG";

    RnaStringReverseComplement modifiedString(str);
    RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_reverse_complement) {
    Rna5String str = "CGAUN";
    Rna5String const EXPECTED_STRING = "NAUCG";

    Rna5StringReverseComplement modifiedString(str);
    Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedStringCopy);
}


// TODO(holtgrew): The following could be made non-redundant with a helper function.


SEQAN_DEFINE_TEST(test_modifer_shortcuts_complement_in_place_string) {
    Dna5String str = "CGATN";
    Dna5String const kStr = str;
    Dna5String const kExpectedResult = "GCTAN";

    // Test non-const version.
    complement(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    complement(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_complement_in_place_string_set) {
    Dna5String str1 = "CCGGTTAANN";
    Dna5String str2 = "CGTANCGTAN";
    Dna5String const EXPECTED_STRING1 = "GGCCAATTNN";
    Dna5String const EXPECTED_STRING2 = "GCATNGCATN";

    StringSet<Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    // Test non-const version.
    {
        StringSet<Dna5String> strSetCopy(strSet);

        complement(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }

    // Test const version.
    {
        StringSet<Dna5String> const strSetCopy(strSet);
        
        complement(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_complement_in_place_string) {
    Dna5String str = "CGATN";
    Dna5String const kStr = str;
    Dna5String const kExpectedResult = "NATCG";

    // Test non-const version.
    reverseComplement(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    reverseComplement(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_complement_in_place_string_set) {
    Dna5String str1 = "CCGGTTAANN";
    Dna5String str2 = "CGTANCGTAN";
    Dna5String const EXPECTED_STRING1 = "NNTTAACCGG";
    Dna5String const EXPECTED_STRING2 = "NTACGNTACG";

    StringSet<Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    // Test non-const version.
    {
        StringSet<Dna5String> strSetCopy = strSet;

        reverseComplement(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }

    // Test const version.
    {
        StringSet<Dna5String> const strSetCopy = strSet;
        
        reverseComplement(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_in_place_string) {
    Dna5String str = "CGATN";
    Dna5String const kStr = str;
    Dna5String const kExpectedResult = "NTAGC";

    // Test non-const version.
    reverse(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    reverse(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_in_place_string_set) {
    Dna5String str1 = "CCGGTTAANN";
    Dna5String str2 = "CGTANCGTAN";
    Dna5String const EXPECTED_STRING1 = "NNAATTGGCC";
    Dna5String const EXPECTED_STRING2 = "NATGCNATGC";

    StringSet<Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    // Test non-const version.
    {
        StringSet<Dna5String> strSetCopy(strSet);

        reverse(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }

    // Test const version.
    {
        StringSet<Dna5String> strSetCopy(strSet);
        
        reverse(strSetCopy);
        SEQAN_ASSERT_EQ(2u, length(strSetCopy));
        SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSetCopy[0]);
        SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSetCopy[1]);
    }
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_lower_in_place_string) {
    CharString str = "This is a test!";
    CharString const kStr = str;
    CharString const kExpectedResult = "this is a test!";

    // Test non-const version.
    toLower(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    toLower(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_lower_in_place_string_set) {
    CharString str1 = "This is a test!";
    CharString str2 = "This is also a test!";
    CharString const EXPECTED_STRING1 = "this is a test!";
    CharString const EXPECTED_STRING2 = "this is also a test!";

    // Test non-const version.
    StringSet<CharString> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    toLower(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSet[0]);
    SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSet[1]);

    // Test const version.
    StringSet<CharString> const kStrSet = strSet;

    toLower(kStrSet);
    SEQAN_ASSERT_EQ(2u, length(kStrSet));
    SEQAN_ASSERT_EQ(EXPECTED_STRING1, kStrSet[0]);
    SEQAN_ASSERT_EQ(EXPECTED_STRING2, kStrSet[1]);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_upper_in_place_string) {
    CharString str = "This is a test!";
    CharString const kStr = str;
    CharString const kExpectedResult = "THIS IS A TEST!";

    // Test non-const version.
    toUpper(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    toUpper(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_upper_in_place_string_set) {
    CharString str1 = "This is a test!";
    CharString str2 = "This is also a test!";
    CharString const EXPECTED_STRING1 = "THIS IS A TEST!";
    CharString const EXPECTED_STRING2 = "THIS IS ALSO A TEST!";

    // Test non-const version.
    StringSet<CharString> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    toUpper(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(EXPECTED_STRING1, strSet[0]);
    SEQAN_ASSERT_EQ(EXPECTED_STRING2, strSet[1]);

    // Test const version.
    StringSet<CharString> const kStrSet = strSet;

    toUpper(kStrSet);
    SEQAN_ASSERT_EQ(2u, length(kStrSet));
    SEQAN_ASSERT_EQ(EXPECTED_STRING1, kStrSet[0]);
    SEQAN_ASSERT_EQ(EXPECTED_STRING2, kStrSet[1]);
}



#endif  // TEST_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_
