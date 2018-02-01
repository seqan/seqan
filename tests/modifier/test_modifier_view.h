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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_TESTS_MODIFIER_TEST_MODIFIER_VIEW_H_
#define SEQAN_TESTS_MODIFIER_TEST_MODIFIER_VIEW_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>

#include "helpers.h"


// Test the modifier view metafunctions.
SEQAN_DEFINE_TEST(test_modifier_view_iterator_metafunctions)
{
    using namespace seqan;

    typedef CaesarChiffre<char> TFunctor;
    typedef ModifiedIterator<CharString, ModView<TFunctor> > TModifiedIterator;

    // TODO(holtgrew): We really want a static assertion here.
    {
        typedef ModViewCargo<TFunctor> TExpected;
        typedef Cargo<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
    }
    {
        typedef char TExpected;
        typedef Value<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
    }
    {
        typedef char TExpected;
        typedef GetValue<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
    }
    {
        typedef char TExpected;
        typedef Reference<TModifiedIterator>::Type TResult;
        bool res = IsSameType<TExpected, TResult>::VALUE;
        SEQAN_ASSERT(res);
    }
}

// Test the modifier view iterator.
SEQAN_DEFINE_TEST(test_modifier_view_iterator)
{
    using namespace seqan;

    typedef CaesarChiffre<char> TFunctor;
    typedef Iterator<CharString, Standard>::Type    TIterator;
    typedef Iterator<CharString, Rooted>::Type      TRootedIterator;
    typedef ModifiedIterator<TIterator, ModView<TFunctor> >         TModifiedIterator;
    typedef ModifiedIterator<TRootedIterator, ModView<TFunctor> >   TModifiedRootedIterator;

    // The string and functor we will work with.
    CharString myString = "This is a nice string!";
    TFunctor myFunctor(1);

    // Manually shift characters in string by one.
    CharString rotatedString = "This is a nice string!";
    for (size_t i = 0; i < length(rotatedString); ++i) {
        if (rotatedString[i] == ' ' || rotatedString[i] == '!')
            continue;
        rotatedString[i] += 1;
    }

    // Test the various ways to construct the iterator with both a
    // container and a functor.
//    {
//        TModifiedIterator it(myFunctor);  // (weese:) removed this one, as it makes no sense to have this constructor
//    }
    {
        TModifiedIterator it(begin(myString, Standard()), myFunctor);
        TModifiedIterator it2(begin(myString, Rooted()), myFunctor);
        TModifiedRootedIterator itR(begin(myString, Rooted()), myFunctor);
        TModifiedIterator it3(itR);
    }
    {
        TModifiedIterator it;
        TModifiedIterator it2(it);
    }
//    {
//        TModifiedIterator const IT(myFunctor);  // (weese:) removed this one, as it makes no sense to have this constructor
//        TModifiedIterator it(IT);
//    }

    // Test value() and getValue().
    {
        TModifiedIterator it(begin(myString, Rooted()), myFunctor);
        TModifiedIterator const IT = it;

        // TODO(holtgrew): The following does not compile.
        SEQAN_ASSERT_EQ('U', value(it));
        SEQAN_ASSERT_EQ('U', value(IT));

        SEQAN_ASSERT_EQ('U', getValue(it));
        SEQAN_ASSERT_EQ('U', getValue(IT));
    }
}

SEQAN_DEFINE_TEST(test_modifier_view_const_iterator)
{
    SEQAN_ASSERT_FAIL("Implement me!");
}

// Test the modified string class with caesar chiffre.
SEQAN_DEFINE_TEST(test_modifier_view_string_caesar_chiffre)
{
    using namespace seqan;

    typedef CaesarChiffre<char> TFunctor;
    TFunctor myFunctor(1);

    CharString originalStr = "This is a test!";
    CharString const EXPECTED_RESULT = "Uijt jt b uftu!";

    // Test the various ways to initialize a ModifiedString.
    // TODO(holtgrew): Should modified strings not be const to the outside?  Lots of non-const functions are superfluous, right?
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(myFunctor);
        setHost(modifiedStr, originalStr);

        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStr);
        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStrCopy);
    }
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStr);
        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStrCopy);
    }
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStr);
        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStrCopy);
    }
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStr);
        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStrCopy);

        ModifiedString<CharString, ModView<TFunctor> > modifiedStr2(modifiedStr);

        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStr2);
        CharString modifiedStrCopy2 = modifiedStr2;
        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStrCopy2);
    }

    // Test operator[].
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);
        SEQAN_ASSERT_EQ('U', modifiedStr[0]);
        SEQAN_ASSERT_EQ('i', modifiedStr[1]);
        SEQAN_ASSERT_EQ('j', modifiedStr[2]);
        SEQAN_ASSERT_EQ('t', modifiedStr[3]);
        SEQAN_ASSERT_EQ(' ', modifiedStr[4]);
    }

    // Test value() and getValue().
    {
        ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);
        SEQAN_ASSERT_EQ('U', value(modifiedStr, 0));
        SEQAN_ASSERT_EQ('U', getValue(modifiedStr, 0));
    }
}

// Test the modified string class with upper case functor.
SEQAN_DEFINE_TEST(test_modifier_view_string_upper_case)
{
    using namespace seqan;

    typedef FunctorUpcase<char> TFunctor;
    TFunctor myFunctor;

    CharString originalStr = "This is a test!";
    CharString const EXPECTED_RESULT = "THIS IS A TEST!";

    ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

    SEQAN_ASSERT_EQ(length(EXPECTED_RESULT), length(modifiedStr));
    for (size_t i = 0; i < length(originalStr); ++i)
        SEQAN_ASSERT_EQ_MSG(EXPECTED_RESULT[i], modifiedStr[i], "i = %lu", i);

    SEQAN_ASSERT_EQ(modifiedStr, EXPECTED_RESULT);
    CharString modifiedStrCopy = modifiedStr;
    SEQAN_ASSERT_EQ(modifiedStrCopy, EXPECTED_RESULT);

    // We do not test the whole interface as in _caesar_chiffre, so this is
    // more a test of FunctorUpcase.
}

// Test the modified string class with low case functor.
SEQAN_DEFINE_TEST(test_modifier_view_string_low_case)
{
    using namespace seqan;

    // TODO(holtgrew): Would it make more sense to name this FunctorDowncase?
    typedef FunctorLowcase<char> TFunctor;
    TFunctor myFunctor;

    CharString originalStr = "This is a test!";
    CharString const EXPECTED_RESULT = "this is a test!";

    ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

    SEQAN_ASSERT_EQ(length(EXPECTED_RESULT), length(modifiedStr));
    for (size_t i = 0; i < length(originalStr); ++i)
        SEQAN_ASSERT_EQ_MSG(EXPECTED_RESULT[i], modifiedStr[i], "i = %lu", i);

    SEQAN_ASSERT_EQ(modifiedStr, EXPECTED_RESULT);
    CharString modifiedStrCopy = modifiedStr;
    SEQAN_ASSERT_EQ(modifiedStrCopy, EXPECTED_RESULT);

    // We do not test the whole interface as in _caesar_chiffre, so this is
    // more a test of FunctorLowcase
}

// Test the modified string class with alphabet conversion.
SEQAN_DEFINE_TEST(test_modifier_view_string_alphabet_conversion)
{
    using namespace seqan;

    typedef FunctorConvert<char, Dna5> TFunctor;
    TFunctor myFunctor;

    CharString originalStr = "acgtnACGTN";
    Dna5String const EXPECTED_RESULT = "ACGTNACGTN";

    ModifiedString<CharString, ModView<TFunctor> > modifiedStr(originalStr, myFunctor);

    SEQAN_ASSERT_EQ(length(EXPECTED_RESULT), length(modifiedStr));
    for (size_t i = 0; i < length(originalStr); ++i)
        SEQAN_ASSERT_EQ_MSG(EXPECTED_RESULT[i], modifiedStr[i], "i = %lu", i);

    // TODO(holtgrew): This does not compile.
//     SEQAN_ASSERT_EQ(modifiedStr, EXPECTED_RESULT);
    CharString modifiedStrCopy = modifiedStr;
    SEQAN_ASSERT_EQ(modifiedStrCopy, EXPECTED_RESULT);

    // We do not test the whole interface as in _caesar_chiffre, so this is
    // more a test of FunctorConvert.
}

// Test the modified string class with nested modifier.
SEQAN_DEFINE_TEST(test_modifier_view_string_nested_modifier)
{
    using namespace seqan;

    {
        typedef CaesarChiffre<char> TFunctor;
        typedef ModifiedString<ModifiedString<CharString, ModView<TFunctor> >, ModView<TFunctor> > TModifiedString;

        CharString originalStr = "This is a test!";
        CharString const EXPECTED_RESULT = "Wklv lv d whvw!";

        TFunctor func1(1), func2(2);
        TModifiedString modifiedStr(originalStr, func1);
        assignModViewFunctor(modifiedStr, func1);
        assignModViewFunctor(host(modifiedStr), func2);

        SEQAN_ASSERT_EQ(modifiedStr, EXPECTED_RESULT);

        CharString modifiedStrCopy = modifiedStr;
        SEQAN_ASSERT_EQ(EXPECTED_RESULT, modifiedStrCopy);
    }
    // test nested reverse/view modifiers and independence of order
    {
        Dna5String str = "CGATN";
        Dna5String const EXPECTED_STRING = "NATCG";

        ModifiedString<
            ModifiedString<    Dna5String, ModView< FunctorComplement<Dna5> > >,
            ModReverse
        > modifiedString1(str);

        ModifiedString<
            ModifiedString<    Dna5String, ModReverse >,
            ModView< FunctorComplement<Dna5> >
        > modifiedString2(str);

        SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedString1);
        SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedString2);
    }
    // test const in different nesting levels
    {
        Dna5String str = "CGATN";
        Dna5String const EXPECTED_STRING = "NATCG";

        // test independence of nesting order
        ModifiedString<
            ModifiedString<    Dna5String const, ModView< FunctorComplement<Dna5> > >,
            ModReverse
        > modifiedString1(str);

        ModifiedString<
            ModifiedString<    Dna5String const, ModReverse >,
            ModView< FunctorComplement<Dna5> >
        > modifiedString2(str);

        ModifiedString<
            ModifiedString<    Dna5String const, ModView< FunctorComplement<Dna5> > > const,
            ModReverse
        > modifiedString3(str);

        ModifiedString<
            ModifiedString<    Dna5String const, ModReverse > const,
            ModView< FunctorComplement<Dna5> >
        > modifiedString4(str);

        ModifiedString<
            ModifiedString<    Dna5String, ModReverse > const,
            ModView< FunctorComplement<Dna5> >
        > modifiedString5(str);

        ModifiedString<
            ModifiedString<    Dna5String, ModReverse > const,
            ModView< FunctorComplement<Dna5> >
        > modifiedString6(str);

        SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedString1);
        SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedString2);
        SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedString3);
        SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedString4);
        SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedString5);
        SEQAN_ASSERT_EQ(EXPECTED_STRING, modifiedString6);
    }
}

// Test the convert() function.
SEQAN_DEFINE_TEST(test_modifier_convert_in_place)
{
    using namespace seqan;

    CharString const originalStr = "This is a test!";
    CharString const expectedResult = "Uijt jt b uftu!";

    // Non-const variant on string.
    {
        CharString strCopy = originalStr;
        convert(strCopy, CaesarChiffre<char>(1));
        SEQAN_ASSERT_EQ(expectedResult, strCopy);
    }

    // Non-const variant on segment.
    {
        CharString strCopy = originalStr;
        Segment<CharString, InfixSegment> stringInfix(strCopy);
        convert(strCopy, CaesarChiffre<char>(1));
        SEQAN_ASSERT_EQ(expectedResult, strCopy);
    }

    // Const variant on segment.
    {
        CharString strCopy = originalStr;
        Segment<CharString, InfixSegment> const stringInfix(strCopy);
        convert(strCopy, CaesarChiffre<char>(1));
        SEQAN_ASSERT_EQ(expectedResult, strCopy);
    }
}

#endif  // SEQAN_TESTS_MODIFIER_TEST_MODIFIER_VIEW_H_
