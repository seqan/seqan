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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
// This file contains functions to test the functionality of the sequence
// module.
// ==========================================================================

#ifndef CORE_TESTS_SEQUENCE_TEST_STRINGSET_H_
#define CORE_TESTS_SEQUENCE_TEST_STRINGSET_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include "test_sequence.h"

// TODO(singer): The Value metafunction of a concat direct string set returns an infix. 
// This should be changed!
// The following metafunction is just a workaround.

template <typename TStringSet>
struct TestStringSetValue_
{
    typedef typename Value<TStringSet>::Type Type;
};

template <typename TString>
struct TestStringSetValue_<StringSet<TString, Owner<ConcatDirect<> > > >
{
    typedef TString Type;
};

template <typename TAlphabetSpecPair>
class StringSetTest : public seqan::Test
{
public:
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 1>::Type TAlphabet;
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 2>::Type TStringSpec;
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 3>::Type TSetSpec;
    typedef seqan::StringSet<seqan::String<TAlphabet, TStringSpec>, TSetSpec> TStringSet;
};

// ((a (b (c)))
//  ((d (e (f)))
//   ((g (h (i)))
// )))

typedef seqan::TagList<
//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
//             seqan::TagList<CountingChar, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
//             seqan::TagList<CountingChar, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
//             seqan::TagList<CountingChar, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
//             seqan::TagList<CountingChar, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<

            seqan::TagList<seqan::Dna, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::External<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<

            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<

            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<

            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<

// TODO(Singer): 7 errors and about 400 warnings (deprecated ...)
//             seqan::TagList<char, seqan::TagList<seqan::CStyle, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::CStyle, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::CStyle, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::CStyle, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<

            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Owner<> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Owner<ConcatDirect<> > > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Dependent<Tight> > > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Dependent<Generous> > > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::Dependent<Generous> > > >//, seqan::TagList<
//         > > > > > > > > > > > > > > > >
        > > > > > > > > > > > > > > > >
        > > > > > > > > > > > >
        > > > > > > > > > > > > > > > >
        > > > > > > > > > > > > > > > >
//         > > > >
        > > > > > > > > > > > > > > > >
        StringSetTestTypes;

template <typename T>
class StringSetTestCommon : public StringSetTest<T>
{};

SEQAN_TYPED_TEST_CASE(StringSetTestCommon, StringSetTestTypes);

// --------------------------------------------------------------------------
// Generic String Tests
// --------------------------------------------------------------------------

template <typename TReturnStringSet, typename TStringSet>
TStringSet createStringSet(TStringSet & stringSet)
{
    return TReturnStringSet(stringSet);
}

// Test whether PrefixSegments are copy constructible.
template <typename TStringSet>
void testStringSetCopyConstructible(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TString str1("AAAA");
    TString str2("CC");
    TString str3("GGG");

    TStringSet stringSet1;
    appendValue(stringSet1, str1);
    appendValue(stringSet1, str2);
    appendValue(stringSet1, str3);

    {
        TStringSet stringSet2(stringSet1);

        SEQAN_ASSERT_EQ(getValue(stringSet2, 0), "AAAA");
        SEQAN_ASSERT_EQ(getValue(stringSet2, 1), "CC");
        SEQAN_ASSERT_EQ(getValue(stringSet2, 2), "GGG");
    }
    {
        TStringSet const stringSet2(stringSet1);

        SEQAN_ASSERT_EQ(getValue(stringSet2, 0), "AAAA");
        SEQAN_ASSERT_EQ(getValue(stringSet2, 1), "CC");
        SEQAN_ASSERT_EQ(getValue(stringSet2, 2), "GGG");
    }
}

// TODO(singer): AppendValue is not available for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetCopyConstructible(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetCopyConstructible(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetCopyConstructible(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}

// Test whether sequences are copy constructible.
SEQAN_TYPED_TEST(StringSetTestCommon, CopyConstructible)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetCopyConstructible(strSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test whether sequences are default constructible.
template <typename TStringSet>
void testStringSetDefaultConstructible(TStringSet & /*Tag*/)
{
    using namespace seqan;

    TStringSet stringSet;
    SEQAN_ASSERT(begin(stringSet) == end(stringSet));
}

// Test whether sequences are copy constructible.
SEQAN_TYPED_TEST(StringSetTestCommon, DefaultConstructible)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetDefaultConstructible(strSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}


// Test operator<().
template <typename TStringSet>
void testStringSetLessGreaterEqual(TStringSet & /*Tag*/)
{
    using namespace seqan;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    // nothing - nohing   
    {
        TNonConstStringSet nonConstStringSet1;
        TNonConstStringSet nonConstStringSet2;
        TStringSet stringSet1, stringSet2;

        std::cerr << (stringSet1 < stringSet2) << std::endl;
        SEQAN_ASSERT_EQ((stringSet1 < stringSet2), false);
        SEQAN_ASSERT((stringSet1 <= stringSet2) == true);
        SEQAN_ASSERT((stringSet1 > stringSet2) == false);
        SEQAN_ASSERT((stringSet1 >= stringSet2) == true);
        SEQAN_ASSERT((stringSet1 == stringSet2) == true);
        SEQAN_ASSERT((stringSet1 != stringSet2) == false);
    }

    // something - nothing.
    {
        TNonConstStringSet nonConstStringSet1;
        TNonConstStringSet nonConstStringSet2;
        resize(nonConstStringSet1, 3);
        TStringSet stringSet1(nonConstStringSet1);
        TStringSet stringSet2(nonConstStringSet2);
        SEQAN_ASSERT((stringSet1 < stringSet2) == false);
        SEQAN_ASSERT((stringSet1 <= stringSet2) == false);
        SEQAN_ASSERT((stringSet1 > stringSet2) == true);
        SEQAN_ASSERT((stringSet1 >= stringSet2) == true);
        SEQAN_ASSERT((stringSet1 == stringSet2) == false);
        SEQAN_ASSERT((stringSet1 != stringSet2) == true);
    }

    // something (initialized) - nothing.
    {
        TNonConstStringSet nonConstStringSet1;
        TNonConstStringSet nonConstStringSet2;
        resize(nonConstStringSet1, 3);
        nonConstStringSet1[0] = "AAAA";
        nonConstStringSet1[1] = "CCCC";
        nonConstStringSet1[2] = "GGGG";
        TStringSet stringSet1(nonConstStringSet1);
        TStringSet stringSet2(nonConstStringSet2);
        SEQAN_ASSERT((stringSet1 < stringSet2) == false);
        SEQAN_ASSERT((stringSet1 <= stringSet2) == false);
        SEQAN_ASSERT((stringSet1 > stringSet2) == true);
        SEQAN_ASSERT((stringSet1 >= stringSet2) == true);
        SEQAN_ASSERT((stringSet1 == stringSet2) == false);
        SEQAN_ASSERT((stringSet1 != stringSet2) == true);
    }

    // nothing - something
    {
        TNonConstStringSet nonConstStringSet1;
        TNonConstStringSet nonConstStringSet2;
        resize(nonConstStringSet2, 3);
        TStringSet stringSet1(nonConstStringSet1);
        TStringSet stringSet2(nonConstStringSet2);
        SEQAN_ASSERT((stringSet1 < stringSet2) == false);
        SEQAN_ASSERT((stringSet1 <= stringSet2) == false);
        SEQAN_ASSERT((stringSet1 > stringSet2) == true);
        SEQAN_ASSERT((stringSet1 >= stringSet2) == true);
        SEQAN_ASSERT((stringSet1 == stringSet2) == false);
        SEQAN_ASSERT((stringSet1 != stringSet2) == true);
    }
    // nothing - something
    {
        TNonConstStringSet nonConstStringSet1;
        TNonConstStringSet nonConstStringSet2;
        resize(nonConstStringSet2, 3);
        nonConstStringSet2[0] = "AAAA";
        nonConstStringSet2[1] = "CCCC";
        nonConstStringSet2[2] = "GGGG";
        TStringSet stringSet1(nonConstStringSet1);
        TStringSet stringSet2(nonConstStringSet2);
        SEQAN_ASSERT((stringSet1 < stringSet2) == true);
        SEQAN_ASSERT((stringSet1 <= stringSet2) == true);
        SEQAN_ASSERT((stringSet1 > stringSet2) == false);
        SEQAN_ASSERT((stringSet1 >= stringSet2) == false);
        SEQAN_ASSERT((stringSet1 == stringSet2) == false);
        SEQAN_ASSERT((stringSet1 != stringSet2) == true);
    }
    // something - something
    {
        TNonConstStringSet nonConstStringSet1;
        TNonConstStringSet nonConstStringSet2;
        resize(nonConstStringSet1, 3);
        nonConstStringSet1[0] = "AAAA";
        nonConstStringSet1[1] = "CCCC";
        nonConstStringSet1[2] = "GGGG";
        resize(nonConstStringSet2, 3);
        nonConstStringSet2[0] = "AAAA";
        nonConstStringSet2[1] = "CCCC";
        nonConstStringSet2[2] = "GGGT";
        TStringSet stringSet1(nonConstStringSet1);
        TStringSet stringSet2(nonConstStringSet2);
        SEQAN_ASSERT((stringSet1 < stringSet2) == true);
        SEQAN_ASSERT((stringSet1 <= stringSet2) == true);
        SEQAN_ASSERT((stringSet1 > stringSet2) == false);
        SEQAN_ASSERT((stringSet1 >= stringSet2) == false);
        SEQAN_ASSERT((stringSet1 == stringSet2) == false);
        SEQAN_ASSERT((stringSet1 != stringSet2) == true);
    }
}

// TODO (singer): comparison operators are not implemented for string sets.
// Test whether sequences are copy constructible.
SEQAN_TYPED_TEST(StringSetTestCommon, Comparison)
{
//     CountingChar::clear();
// 
//     typename TestFixture::TStringSet strSet;
//     testStringSetLessGreaterEqual(strSet);
// 
//     typename TestFixture::TStringSet const constStrSet;
//     testStringSetLessGreaterEqual(constStrSet);
// 
//     testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of append().
template <typename TStringSet>
void testStringSetAppend(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TStringSet stringSet1;
    TStringSet stringSet2;

    // Test the append function
    append(stringSet1, stringSet2);
    SEQAN_ASSERT_EQ(length(stringSet1), 0u);

    // Test the appendValue function
    resize(stringSet2, 3);
    append(stringSet1, stringSet2);
    SEQAN_ASSERT_EQ(length(stringSet1), 3u);
}

// TODO(singer): append not implemented for string sets
SEQAN_TYPED_TEST(StringSetTestCommon, Append)
{
//     CountingChar::clear();
// 
//     typename TestFixture::TStringSet strSet;
//     testStringSetAppend(strSet);
// 
//     testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of appendValue().
template <typename TStringSet>
void testStringSetAppendValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TStringSet stringSet;
    TStringSet stringSet2;

    // Test the appendValue function
    TString str("ACGT");
    appendValue(stringSet, str);
    SEQAN_ASSERT_EQ(length(stringSet), 1u);
    SEQAN_ASSERT_EQ(stringSet[0], "ACGT");

    // Test the appendValue function
    str = TString("CGTA");
    appendValue(stringSet, str);
    SEQAN_ASSERT_EQ(length(stringSet), 2u);
    SEQAN_ASSERT_EQ(stringSet[0], "ACGT");
    SEQAN_ASSERT_EQ(stringSet[1], "CGTA");
}

// TODO(singer): AppendValue is not available for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetAppendValue(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAppendValue(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAppendValue(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, AppendValue)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetAppendValue(strSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of assign().
template <typename TStringSet>
void testStringSetAssign(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    {
        TStringSet stringSet1;
        TStringSet stringSet2;
        
        assign(stringSet1, stringSet2);
        SEQAN_ASSERT_EQ(length(stringSet1), 0u);
        SEQAN_ASSERT(begin(stringSet1) == end(stringSet1));
    }
    {
        TStringSet stringSet1;
        TStringSet stringSet2;

        resize(stringSet2, 3u);
        assign(stringSet1, stringSet2);
        for(unsigned i = 0; i < length(stringSet1); ++i)
            SEQAN_ASSERT(stringSet1[1] == stringSet2[i]);
    }
}

// TODO(singer): Seg fault
template <typename TValue, typename TStringSetSpec>
void testStringSetAssign(StringSet<String<TValue, External<> >, Owner<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSpec, typename TStringSetSpec>
void testStringSetAssign(StringSet<String<TValue, TStringSpec>, Dependent<TStringSetSpec> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Assign)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetAssign(strSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of assignValue().
template <typename TStringSet>
void testStringSetAssignValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TStringSet stringSet;
    TString string;

    // Test the assignValue function
    resize(stringSet, 3u);
    assignValue(stringSet, 1u, string);
    SEQAN_ASSERT_EQ(length(stringSet), 3u);
    SEQAN_ASSERT_EQ(stringSet[1], string);
}

template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValue(StringSet<String<TValue, Block<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValue(StringSet<String<TValue, Block<> >, Owner<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValue(StringSet<String<TValue, Block<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValue(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValue(StringSet<String<TValue, Packed<> >, Owner<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValue(StringSet<String<TValue, Packed<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}

// TODO(singer): Seg fault
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValue(StringSet<String<TValue, External<> >, Owner<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSpec, typename TStringSetSpec>
void testStringSetAssignValue(StringSet<String<TValue, TStringSpec>, Dependent<TStringSetSpec> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, AssignValue)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetAssignValue(strSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of assignValueById().
template <typename TStringSet>
void testStringSetAssignValueById(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TStringSet stringSet1;

    // Assigning a string.
    TString string("ACGT");
    unsigned id = assignValueById(stringSet1, string);
    SEQAN_ASSERT_EQ(length(stringSet1), 1u);
    SEQAN_ASSERT_EQ(stringSet1[0], string);
    SEQAN_ASSERT_EQ(id, 0u);

    // TODO (singer): this is not documented, such that it can be easily understood!
    // Assigning a string from a string set.
    TString str2("AAAA");
    TString str3("TTTT");
    
    TStringSet stringSet2;
    appendValue(stringSet2, str2);
    appendValue(stringSet2, str3);
    stringSet2[1] = str3;
    id = assignValueById(stringSet1, stringSet2, 1u);

    SEQAN_ASSERT_EQ(length(stringSet1), 2u);
    SEQAN_ASSERT_EQ(stringSet1[0], string); 
    SEQAN_ASSERT_EQ(stringSet1[1], str3);
    SEQAN_ASSERT_EQ(id, 1u);
}

// TODO(singer): Seg fault
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValueById(StringSet<String<TValue, Block<> >, Owner<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValueById(StringSet<String<TValue, Array<100> >, Owner<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValueById(StringSet<String<TValue, Packed<> >, Owner<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetAssignValueById(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue>
void testStringSetAssignValueById(StringSet<String<TValue, External<> >, Owner<ConcatDirect<> > > & /*Tag*/) {}
template <typename TValue>
void testStringSetAssignValueById(StringSet<String<TValue, Alloc<> >, Owner<ConcatDirect<> > > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, AssignValueById)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetAssignValueById(strSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of back() for non const strings.
// TODO (singer): not in docu.
template <typename TStringSet>
void testStringSetBack(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    
    TString str("CCCC");
    TString str2("ACGT");
    TStringSet stringSet;
    appendValue(stringSet, str);
    appendValue(stringSet, str2);
    
    // val is a reference in contrast to the const version of back()
    TString & val = back(stringSet);
    val = "TTTT";
    SEQAN_ASSERT_EQ(val, stringSet[1]);
}

// Test of back() for const strings.
// TODO (singer): not in docu.
template <typename TStringSet>
void testStringSetBack(TStringSet const & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TString str("CCCC");
    TString str2("ACGT");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    TStringSet stringSet(nonConstStringSet);

    // val is a reference in contrast to the const version of back()
    TString val = back(stringSet);
    SEQAN_ASSERT_EQ(val, stringSet[1]);
    val = "TTTT";
    SEQAN_ASSERT_EQ(stringSet[1], str2);
}

// TODO(singer)
template <typename TValue, typename TStringSpec>
void testStringSetBack(StringSet<String<TValue, TStringSpec>, Owner<ConcatDirect<> > > & /*Tag*/) {}
template <typename TValue, typename TStringSpec>
void testStringSetBack(StringSet<String<TValue, TStringSpec>, Owner<ConcatDirect<> > > const & /*Tag*/) {}
template <typename TValue> 
void testStringSetBack(StringSet<String<TValue, External<> >, Owner<> > & /*Tag*/) {}
template <typename TValue> 
void testStringSetBack(StringSet<String<TValue, External<> >, Owner<> > const & /*Tag*/) {}
template <typename TValue, typename TStringSpec, typename TStringSetSpec>
void testStringSetBack(StringSet<String<TValue, TStringSpec>, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSpec, typename TStringSetSpec>
void testStringSetBack(StringSet<String<TValue, TStringSpec>, Dependent<TStringSetSpec> > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Back)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetBack(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetBack(constStrSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of begin().
template <typename TStringSet>
void testStringSetBegin(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TString str("ACGT");
    TString str2("ATTT");

    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    TStringSet stringSet(nonConstStringSet);
    SEQAN_ASSERT_EQ(*begin(stringSet), str);
    SEQAN_ASSERT_EQ(*begin(stringSet, Standard()), str);
    SEQAN_ASSERT_EQ(*begin(stringSet, Rooted()), str);
}
// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetBegin(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetBegin(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetBegin(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetBegin(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}


template <typename TValue, typename TStringSetSpec>
void testStringSetBegin(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetBegin(StringSet<String<TValue, MMap<> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Begin)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetBegin(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetBegin(constStrSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of beginPosition().
template <typename TStringSet>
void testStringSetBeginPosition(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    // Test on an empty string.
    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3u);
    {
        TStringSet stringSet(nonConstStringSet);
        SEQAN_ASSERT_EQ(beginPosition(stringSet), 0u);
    }

    // Test on a non empty string.
    {
        resize(nonConstStringSet, 0u);
        TString str("ACGT");
        appendValue(nonConstStringSet, str); 
        appendValue(nonConstStringSet, str); 
        TStringSet stringSet(nonConstStringSet);
        SEQAN_ASSERT_EQ(beginPosition(stringSet), 0u);
    }
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetBeginPosition(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetBeginPosition(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetBeginPosition(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetBeginPosition(StringSet<String<TValue, MMap<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetBeginPosition(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetBeginPosition(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, BeginPosition)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetBeginPosition(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetBeginPosition(constStrSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of clear().
template <typename TStringSet>
void testStringSetClear(TStringSet & /*Tag*/)
{
    using namespace seqan;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    // Test on an empty string set.
    {
        TStringSet stringSet;
        clear(stringSet);
        SEQAN_ASSERT(begin(stringSet) == end(stringSet));
    }

    // Test on a non empty string set.
    {
        TString str("ACGT");
        TStringSet stringSet;
        appendValue(stringSet, str);
        clear(stringSet);
        SEQAN_ASSERT(begin(stringSet) == end(stringSet));
    }
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetClear(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetClear(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetClear(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetClear(StringSet<String<TValue, MMap<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetClear(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetClear(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Clear)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetClear(strSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of concat().
template <typename TStringSet>
void testStringSetConcat(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename Concatenator<TStringSet>::Type TConcat;

    // TODO (singer): test fails for not initialized string sets.
    // Error message: Assertion failed : static_cast<TStringPos>(pos) < static_cast<TStringPos>(length(me)) was: 0 >= 0 (Trying to access an element behind the last one!).
    // Test on an empty string set.
    // {
    //     TStringSet stringSet;
    //     TConcat concatString = concat(stringSet);
    //     SEQAN_ASSERT(begin(concatString) == end(concatString));
    // }

    // Test on a non empty string set.
    {
        TString str1("AAAA");
        TString str2("CCCC");
        TString str3("GGGG");
        TString str4("TTTT");

    	TNonConstStringSet nonConstStringSet;
    	appendValue(nonConstStringSet, str1); 
    	appendValue(nonConstStringSet, str2); 
    	appendValue(nonConstStringSet, str3); 
    	appendValue(nonConstStringSet, str4); 
    	TString string("AAAACCCCGGGGTTTT");
        TStringSet stringSet(nonConstStringSet);
        TConcat concatString = concat(stringSet);
        for (unsigned i = 0; i < length(string); ++i)
        	SEQAN_ASSERT_EQ(string[i], concatString[i]);
    }
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetConcat(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetConcat(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetConcat(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetConcat(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Concat)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetConcat(strSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of end().
template <typename TStringSet>
void testStringSetEnd(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TNonConstStringSet nonConstStringSet;
    TString str("ACGT");
    TString str2("AAAA");
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    TStringSet stringSet(nonConstStringSet);

     typename Iterator<TStringSet>::Type iter = end(stringSet);
     typename Iterator<TStringSet, Standard>::Type standardIter = end(stringSet, Standard());
     typename Iterator<TStringSet, Rooted>::Type rootedIter = end(stringSet, Rooted());

    --iter;
    --standardIter;
    --rootedIter;

    SEQAN_ASSERT_EQ(*iter, str2);
    SEQAN_ASSERT_EQ(*standardIter, str2);
    SEQAN_ASSERT_EQ(*rootedIter, str2);
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetEnd(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetEnd(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetEnd(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetEnd(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, End)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetEnd(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetEnd(constStrSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of endPosition().
// TODO (singer): erase() is not in the docu.
template <typename TStringSet>
void testStringSetEndPosition(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    // Test on an empty string.
    TNonConstStringSet nonConstStringSet;
    resize(nonConstStringSet, 3u);
    {
        TStringSet stringSet(nonConstStringSet);
        SEQAN_ASSERT_EQ(endPosition(stringSet), 3u);
    }

    // Test on a non empty string.
    {
        TString str("ACGT");
        resize(nonConstStringSet, 0u);
        appendValue(nonConstStringSet, str);
        appendValue(nonConstStringSet, str);
        appendValue(nonConstStringSet, str);
        TStringSet stringSet(nonConstStringSet);
        SEQAN_ASSERT_EQ(endPosition(stringSet), 3u);
    }
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetEndPosition(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetEndPosition(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetEndPosition(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetEndPosition(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, EndPosition)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetEndPosition(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetEndPosition(constStrSet);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of erase().
// TODO (singer): erase() is not in the docu.
template <typename TStringSet>
void testStringSetErase(TStringSet & /*Tag*/)
{
    using namespace seqan;

    // Test on an empty string.

    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TString str("A");
    TString str2("ACGTACGTACGT");
    TString str3("CG");
    TStringSet stringSet;
    appendValue(stringSet, str);
    appendValue(stringSet, str2);
    appendValue(stringSet, str3);
    SEQAN_ASSERT_EQ(stringSet[1], str2);
    erase(stringSet, 1);
    SEQAN_ASSERT_EQ(length(stringSet), 2u);
    SEQAN_ASSERT_EQ(stringSet[0], TString());
    SEQAN_ASSERT_EQ(stringSet[1], TString());
}

// TODO(singer): Seg. fault
template <typename TValue>
void testStringSetErase(StringSet<String<TValue, MMap<> >, Owner<> > & /*Tag*/) {}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue>
void testStringSetErase(StringSet<String<TValue, Block<> >, Owner<> > & /*Tag*/) {}
// template <typename TValue, typename TStringSetSpec>
// void testStringSetErase(StringSet<String<TValue, Packed<> >, Owner<ConcatDirect<> > > & /*Tag*/) {}
// template <typename TValue, typename TStringSetSpec>
// void testStringSetErase(StringSet<String<TValue, Array<100> >, Owner<ConcatDirect<> > > & /*Tag*/) {}

// TODO(singer)
template <typename TValue, typename TStringSpec, typename TStringSetSpec>
void testStringSetErase(StringSet<String<TValue, TStringSpec>, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSpec>
void testStringSetErase(StringSet<String<TValue, TStringSpec>, Owner<ConcatDirect<> > > & /*Tag*/) {}
template <typename TValue>
void testStringSetErase(StringSet<String<TValue, External<> >, Owner<ConcatDirect<> > > & /*Tag*/) {}
template <typename TValue>
void testStringSetErase(StringSet<String<TValue, External<> >, Owner<ConcatDirect<> > > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Erase)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetErase(strSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of eraseBack().
// TODO (singer): eraseBack() is not in the docu.
template <typename TStringSet>
void testStringSetEraseBack(TStringSet & /*Tag*/)
{
    using namespace seqan;


    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TStringSet stringSet;
    // Test on an empty string.
    // TODO (singer): does not work on empty string
    // Assertion failed : length(me) > 0u was: 0 <= 0 (String must have more than 0 characters in eraseBack()!)
    //eraseBack(stringSet);

    // Test on a non empty string.
    TString str;
    TString str2("ACGTACGTACGT");
    appendValue(stringSet, str);
    appendValue(stringSet, str);
    appendValue(stringSet, str2);
    SEQAN_ASSERT_EQ(stringSet[2], str2);
    eraseBack(stringSet);
    SEQAN_ASSERT_EQ(length(stringSet), 2u);
    SEQAN_ASSERT_EQ(stringSet[0], TString());
    SEQAN_ASSERT_EQ(stringSet[1], TString());
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetEraseBack(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetEraseBack(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetEraseBack(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetEraseBack(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, EraseBack)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetEraseBack(strSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// TODO (singer): not in docu.
// Test of front().
template <typename TStringSet>
void testStringSetFront(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TString str("ACGT");
    TString str2("T");
    TStringSet stringSet;
    appendValue(stringSet, str);
    appendValue(stringSet, str2);

    // val is a reference in contrast to the const version of front()
    TString & val = front(stringSet);
    val = "TTTT";
    SEQAN_ASSERT_EQ(val, stringSet[0]);
}

// Test of front() for const strings.
// TODO (singer): not in docu.
template <typename TStringSet>
void testStringSetFront(TStringSet const & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TString str("ACGT");
    TString str2("T");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    TStringSet stringSet(nonConstStringSet);

    // val is a reference in contrast to the const version of front()
    TString val = front(stringSet);
    SEQAN_ASSERT_EQ(val, stringSet[0]);
    val = "TTTT";
    SEQAN_ASSERT_EQ(stringSet[0], str);
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetFront(StringSet<String<TValue, Packed<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetFront(StringSet<String<TValue, Packed<> >, Dependent<TStringSetSpec> > const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetFront(StringSet<String<TValue, Array<100> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetFront(StringSet<String<TValue, Array<100> >, Dependent<TStringSetSpec> > const & /*Tag*/) {}

// TODO(singer)
template <typename TValue>
void testStringSetFront(StringSet<String<TValue, MMap<> >, Owner<> > & /*Tag*/) {}
template <typename TValue>
void testStringSetFront(StringSet<String<TValue, MMap<> >, Owner<> > const & /*Tag*/) {}
template <typename TValue, typename TStringSpec>
void testStringSetFront(StringSet<String<TValue, TStringSpec>, Owner<ConcatDirect<> > > & /*Tag*/) {}
template <typename TValue, typename TStringSpec>
void testStringSetFront(StringSet<String<TValue, TStringSpec>, Owner<ConcatDirect<> > > const & /*Tag*/) {}
template <typename TValue>
void testStringSetFront(StringSet<String<TValue, External<> >, Owner<> > & /*Tag*/) {}
template <typename TValue>
void testStringSetFront(StringSet<String<TValue, External<> >, Owner<> > const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetFront(StringSet<String<TValue, External<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetFront(StringSet<String<TValue, External<> >, Dependent<TStringSetSpec> > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Front)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetFront(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetFront(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of getValue().
template <typename TStringSet>
void testStringSetGetValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

     // In contrast to value(), getValue() does not return a reference but a copy.
     // We test this using the variable value_.
    TString str("CG");
    TString str2("ACGT");
    TString str3("CGACGT");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    appendValue(nonConstStringSet, str3);
    TStringSet stringSet(nonConstStringSet);
    SEQAN_ASSERT_EQ(getValue(stringSet, 1), str2);
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetGetValue(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetGetValue(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetGetValue(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetGetValue(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, GetValue)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetGetValue(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetGetValue(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// TODO (singer): not defined for const string sets.
// Test of getValueById().
template <typename TStringSet>
void testStringSetGetValueById(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TString str("CG");
    TString str2("ACGT");
    TString str3("CGACGT");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    appendValue(nonConstStringSet, str3);
    TStringSet stringSet(nonConstStringSet);
    SEQAN_ASSERT_EQ(getValueById(stringSet, typename Id<TStringSet>::Type(1)), str2);
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetGetValueById(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetGetValueById(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetGetValueById(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetGetValueById(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, GetValueById)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetGetValueById(strSet);

// TODO (singer): not defined for const string sets.
//     typename TestFixture::TStringSet const constStrSet;
//     testStringSetGetValueById(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// TODO (singer): define behaviour and adjust test.
// Infix() compiles and does what it is supposed to do?!
// However, it is not very intuitive. For details see comments below.
// There is a need to improve the documentation of this!
// Test of infix()
template <typename TStringSet>
void testStringSetInfix(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TString str("AAAA");
    TString str2("CCC");
    TString str3("GGGG");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    appendValue(nonConstStringSet, str3);

    // Only non-const test for this scenario possible.
    TStringSet stringSet(nonConstStringSet);
    TString string = infix(stringSet, 0, 1); // Returns the first character not string!
    // std::cerr << string << std::endl; -> "A"

    // Only non-const test for this scenario possible.
    // TString str2("TT");
    // nonConstStringSet[0] = str2;

    // Since the infix (should) point to the fist element it should point to "T"
    // std::cerr << string << std::endl; -> "A"
    // Therefore there is a different behaviour to normal strings.
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetInfix(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfix(StringSet<String<TValue, MMap<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfix(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfix(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfix(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfix(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

// TODO(singer)
template <typename TValue, typename TStringSetSpec>
void testStringSetInfix(StringSet<String<TValue, External<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue,typename TStringSetSpec>
void testStringSetInfix(StringSet<String<TValue, External<> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Infix)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetInfix(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetInfix(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// TODO (singer): define behaviour and adjust test.
// Infix() compiles and does what it is supposed to do?!
// However, it is not very intuitive. For details see comments below.
// There is a need to improve the documentation of this!
// Test of infixWithLength()
template <typename TStringSet>
void testStringSetInfixWithLength(TStringSet & /*Tag*/)
{
using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;

    TString str("AAAA");
    TString str2("CCC");
    TString str3("GGGG");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    appendValue(nonConstStringSet, str3);

    // Only non-const test for this scenario possible.
    TStringSet stringSet(nonConstStringSet);
    TString string = infixWithLength(stringSet, 0, 1); // Returns the first character not string!
    // std::cerr << string << std::endl; -> "A"

    // Only non-const test for this scenario possible.
    // TString str2("TT");
    // nonConstStringSet[0] = str2;

    // Since the infix (should) point to the fist element it should point to "T"
    // std::cerr << string << std::endl; -> "A"
    // Therefore there is a different behaviour to normal strings.
}
// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetInfixWithLength(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfixWithLength(StringSet<String<TValue, MMap<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfixWithLength(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfixWithLength(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfixWithLength(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetInfixWithLength(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

// TODO(singer)
template <typename TValue, typename TStringSetSpec>
void testStringSetInfixWithLength(StringSet<String<TValue, External<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue,typename TStringSetSpec>
void testStringSetInfixWithLength(StringSet<String<TValue, External<> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, InfixWithLength)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetInfixWithLength(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetInfixWithLength(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of insert().
// TODO (singer): no insert function implemented.
// template <typename TStringSet>
// void testStringSetInsert(TStringSet & /*Tag*/)
// {
//     using namespace seqan;
// 
//     // Test of inserting an empty string.
//     TStringSet stringSet1;
//     resize(stringSet1, 1u);
//     TStringSet stringSet2;
//     insert(stringSet1, 0u, stringSet2);
//     SEQAN_ASSERT_EQ(length(stringSet1), 1u);
// 
//     resize(stringSet2, 3u);
//     stringSet2[0] = "ACGT";
//     insert(stringSet1, 0u, stringSet2);
//     SEQAN_ASSERT_EQ(length(stringSet1), 4u);
//     SEQAN_ASSERT_EQ(stringSet1[1], "ACGT");
// }
// 
// SEQAN_TYPED_TEST(StringSetTestCommon, Insert)
// {
//     CountingChar::clear();
// 
//     typename TestFixture::TStringSet strSet;
//     testStringSetInsert(strSet);
//     
//     typename TestFixture::TStringSet const constStrSet;
//     testStringSetInsert(constStrSet);
//     
//     testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
// }
// 
// Test of insertValue().
// TODO (singer): no insertValue function implemented.
// template <typename TStringSet>
// void testStringSetInsertValue(TStringSet & /*Tag*/)
// {
//     using namespace seqan;
//     
//     typedef typename TestStringSetValue_<TStringSet>::Type TString;
// 
//     // Test of inserting an empty string.
//     TStringSet stringSet;
//     resize(stringSet, 1u);
//     insertValue(stringSet, 0, "ACGT");
//     SEQAN_ASSERT_EQ(length(stringSet), 2u);
//     SEQAN_ASSERT_EQ(stringSet[0], "ACGT");
//     SEQAN_ASSERT_EQ(stringSet[1], TString());
// }
// 
// SEQAN_TYPED_TEST(StringSetTestCommon, InsertValue)
// {
//     CountingChar::clear();
// 
//     typename TestFixture::TStringSet strSet;
//     testStringSetInsertValue(strSet);
//     
//     typename TestFixture::TStringSet const constStrSet;
//     testStringSetInsertValue(constStrSet);
//     
//     testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
// }

// Test of iter().
template <typename TStringSet>
void testStringSetIter(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    typedef typename Iterator<TStringSet>::Type TIterator;
    typedef typename Iterator<TStringSet, Standard>::Type TStandardIterator;
    typedef typename Iterator<TStringSet, Rooted>::Type TRootedIterator;

    // Test on an empty string set.
    {
    	TStringSet stringSet;
        TIterator iterator = iter(stringSet, 0);
        TStandardIterator standardIterator = iter(stringSet, 0);
        TRootedIterator rootedIterator = iter(stringSet, 0);
        SEQAN_ASSERT(iterator == begin(stringSet));
        SEQAN_ASSERT(standardIterator == begin(stringSet, Standard()));
        SEQAN_ASSERT(rootedIterator == begin(stringSet, Rooted()));
    }

    // Test on a non empty stringSet.
    {
        TString str1("AAAA");
        TString str2("CCCC");
        TString str3("GGGG");
    	TNonConstStringSet nonConstStringSet;
    	appendValue(nonConstStringSet, str1);
    	appendValue(nonConstStringSet, str2);
    	appendValue(nonConstStringSet, str3);
    	TStringSet stringSet(nonConstStringSet);
        TIterator iterator = iter(stringSet, 0);
        TStandardIterator standardIterator = iter(stringSet, 0);
        TRootedIterator rootedIterator = iter(stringSet, 0);
        SEQAN_ASSERT_EQ(getValue(iterator), "AAAA");
        SEQAN_ASSERT(getValue(iterator) == getValue(stringSet, 0));
        SEQAN_ASSERT(getValue(standardIterator) == getValue(stringSet, 0));
        SEQAN_ASSERT(getValue(rootedIterator) == getValue(stringSet, 0));
    }

    // Test on a non empty stringSet.
    {
        TString str1("AAAA");
        TString str2("CCCC");
        TString str3("GGGG");
        TString str4("TTTT");
    	TNonConstStringSet nonConstStringSet;
    	appendValue(nonConstStringSet, str1);
    	appendValue(nonConstStringSet, str2);
    	appendValue(nonConstStringSet, str3);
    	appendValue(nonConstStringSet, str4);
    	TStringSet stringSet(nonConstStringSet);
        TIterator iterator = iter(stringSet, 3);
        TStandardIterator standardIterator = iter(stringSet, 3);
        TRootedIterator rootedIterator = iter(stringSet, 3);
        SEQAN_ASSERT_EQ(getValue(iterator), "TTTT");
        SEQAN_ASSERT(getValue(iterator) == getValue(stringSet, 3));
        SEQAN_ASSERT(getValue(standardIterator) == getValue(stringSet, 3));
        SEQAN_ASSERT(getValue(rootedIterator) == getValue(stringSet, 3));
    }
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetIter(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetIter(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetIter(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetIter(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}


SEQAN_TYPED_TEST(StringSetTestCommon, Iter)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetIter(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetIter(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of length().
template <typename TStringSet>
void testStringSetLength(TStringSet & /*Tag*/)
{
    using namespace seqan;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    // Test on an empty string.
    {
    	TStringSet stringSet;
    	SEQAN_ASSERT_EQ(length(stringSet), 0u);
    }

    // Test on a non empty string.
    {
    	TNonConstStringSet nonConstStringSet;
    	resize(nonConstStringSet, 10u);
    	TStringSet stringSet(nonConstStringSet);
    	SEQAN_ASSERT_EQ(length(stringSet), 10u);
    }
}

SEQAN_TYPED_TEST(StringSetTestCommon, Length)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetLength(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetLength(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// Test of moveValue().
template <typename TStringSet>
void testStringSetMoveValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    TStringSet stringSet;

    resize(stringSet, 2u);
    moveValue(stringSet, 1, "ACGT");
    SEQAN_ASSERT_EQ(stringSet[1], "ACGT");
}

// TODO(singer): No moveValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetMoveValue(StringSet<String<TValue, Alloc<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetMoveValue(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetMoveValue(StringSet<String<TValue, Block<> >, Owner<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetMoveValue(StringSet<String<TValue, Packed<> >, Owner<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetMoveValue(StringSet<String<TValue, Array<100> >, Owner<TStringSetSpec> > & /*Tag*/) {}

// TODO(singer): Seg fault
template <typename TValue, typename TStringSpec, typename TStringSetSpec>
void testStringSetMoveValue(StringSet<String<TValue, TStringSpec>, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetMoveValue(StringSet<String<TValue, External<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetMoveValue(StringSet<String<TValue, External<> >, Owner<TStringSetSpec> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, MoveValue)
{
    // TODO(singer); Simply not working at all!
//     CountingChar::clear();
// 
//     typename TestFixture::TStringSet strSet;
//     testStringSetMoveValue(strSet);
//     
//     testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// TODO (singer): see infix.
// Test of prefix().
template <typename TStringSet>
void testStringSetPrefix(TStringSet & /*Tag*/)
{
    using namespace seqan;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TString str("ACGTACGT");
    TString str2("GTACGT");
    TString str3("TACGT");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    appendValue(nonConstStringSet, str3);
    TStringSet stringSet(nonConstStringSet);
    TString pref = prefix(stringSet, 3);
    TString str4("ACG");
    SEQAN_ASSERT_EQ(pref, str4);
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetPrefix(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetPrefix(StringSet<String<TValue, MMap<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetPrefix(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetPrefix(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetPrefix(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetPrefix(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}

// TODO(singer)
template <typename TValue, typename TStringSetSpec>
void testStringSetPrefix(StringSet<String<TValue, External<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetPrefix(StringSet<String<TValue, External<> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Prefix)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetPrefix(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetPrefix(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// TODO (singer); replace is not defined for string sets.
//// Test of replace().
//template <typename TStringSet>
//void testStringSetReplace(TStringSet & /*Tag*/)
//{
//    using namespace seqan;
//
//    TStringSet stringSet1;
//    TStringSet stringSet2;
//
//    // This is problematic according to the documentation.
//    // 0 can be a position or an iterator causing compiler errors.
//    replace(stringSet1, 0, 0, stringSet2);
//    SEQAN_ASSERT_EQ(stringSet1, TStringSet());
//
//    appendValue(stringSet1, "AAAA");
//    appendValue(stringSet1, "CCCC");
//    appendValue(stringSet1, "GGGG");
//    appendValue(stringSet1, "TTTT");
//
//    replace(stringSet1, 1, 3, stringSet2);
//    SEQAN_ASSERT_EQ(length(stringSet1), 2u);
//    SEQAN_ASSERT_EQ(stringSet1[0], "AAAA");
//    SEQAN_ASSERT_EQ(stringSet1[1], "TTTT");
//
//    replace(stringSet2, 0, 0, stringSet1);
//    SEQAN_ASSERT_EQ(stringSet2[0], "AAAA");
//    SEQAN_ASSERT_EQ(stringSet2[1], "TTTT");
//
//    clear(stringSet1);
//    clear(stringSet2);
//    appendValue(stringSet1, "AAAA");
//    appendValue(stringSet1, "CCCC");
//    appendValue(stringSet2, "GGGG");
//    appendValue(stringSet2, "TTTT");
//
//    replace(stringSet1, 1, 2, stringSet2);
//    SEQAN_ASSERT_EQ(stringSet1[0], "AAAA");
//    SEQAN_ASSERT_EQ(stringSet2[1], "GGGG");
//    SEQAN_ASSERT_EQ(stringSet2[2], "TTTT");
//
//}
// 
// SEQAN_TYPED_TEST(StringSetTestCommon, Replace)
// {
//     CountingChar::clear();
// 
//     typename TestFixture::TStringSet strSet;
//     testStringSetReplace(strSet);
//         
//     testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
// }

// Test of resize().
template <typename TStringSet>
void testStringSetResize(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TStringSet stringSet;

    resize(stringSet, 0u);
    SEQAN_ASSERT_EQ(length(stringSet), 0u);

    resize(stringSet, 10u);
    SEQAN_ASSERT_EQ(length(stringSet), 10u);
    // TODO (singer): resize should initialize newly allocated memory,
    // which it does not at the moment!
    SEQAN_ASSERT_EQ(stringSet[0], TString());
    SEQAN_ASSERT_EQ(stringSet[9], TString());

    resize(stringSet, 0u);
    SEQAN_ASSERT_EQ(length(stringSet), 0u);

    // TODO (singer): resize with a provided value is not possible with string sets.
//    TString string = "ACGT";
//    resize(stringSet, 10, string);
//    SEQAN_ASSERT_EQ(stringSet[0], string);
//    SEQAN_ASSERT_EQ(stringSet[9], string);
}

template <typename TValue, typename TStringSetSpec>
void testStringSetResize(StringSet<String<TValue, External<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetResize(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetResize(StringSet<String<TValue, External<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetResize(StringSet<String<TValue, Packed<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSpec, typename TStringSetSpec>
void testStringSetResize(StringSet<String<TValue, TStringSpec>, Dependent<TStringSetSpec> > & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Resize)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetResize(strSet);
        
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}

// TODO (singer): see infix.
// Test of suffix().
template <typename TStringSet>
void testStringSetSuffix(TStringSet & /*Tag*/)
{
    using namespace seqan;
    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    TString str("ACGTACGT");
    TString str2("GTACGT");
    TString str3("TACGT");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str);
    appendValue(nonConstStringSet, str2);
    appendValue(nonConstStringSet, str3);
    TStringSet stringSet(nonConstStringSet);
    TString pref = suffix(stringSet, 5);
    TString str4("CGT");
    SEQAN_ASSERT_EQ(pref, str4);
}

// TODO(singer): No appendValue for string sets of packed strings
template <typename TValue, typename TStringSetSpec>
void testStringSetSuffix(StringSet<String<TValue, MMap<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetSuffix(StringSet<String<TValue, MMap<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetSuffix(StringSet<String<TValue, Packed<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetSuffix(StringSet<String<TValue, Packed<> >, TStringSetSpec> const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetSuffix(StringSet<String<TValue, Array<100> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetSuffix(StringSet<String<TValue, Array<100> >, TStringSetSpec> const & /*Tag*/) {}


// TODO(singer)
template <typename TValue, typename TStringSetSpec>
void testStringSetSuffix(StringSet<String<TValue, External<> >, TStringSetSpec> & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetSuffix(StringSet<String<TValue, External<> >, TStringSetSpec> const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Suffix)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetSuffix(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetSuffix(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}


// TODO (singer): swap is not working because a constructor string set (StringSet(StringSet, Move)) is missing.
// Test of swap().
// template <typename TStringSet>
// void testStringSetSwap(TStringSet & /*Tag*/)
// {
//    using namespace seqan;
//
//    TStringSet stringSet1;
//    TStringSet stringSet2;
//    TStringSet stringSet3 = stringSet1;
//    TStringSet stringSet4 = stringSet2;
//
//    swap(stringSet1, stringSet2);
//    SEQAN_ASSERT(stringSet1 == stringSet4);
//    SEQAN_ASSERT(stringSet2 == stringSet3);
//
//    appendValue(stringSet1, "ACGT");
//    stringSet3 = stringSet1;
//    stringSet4 = stringSet2;
//
//    swap(stringSet1, stringSet2);
//    SEQAN_ASSERT(stringSet1 == stringSet4);
//    SEQAN_ASSERT(stringSet2 == stringSet3);
//
//    clear(stringSet1);
//    appendValue(stringSet2, "ACGT");
//    stringSet3 = stringSet1;
//    stringSet4 = stringSet2;
//
//    swap(stringSet1, stringSet2);
//    SEQAN_ASSERT(stringSet1 == stringSet4);
//    SEQAN_ASSERT(stringSet2 == stringSet3);
//
//    stringSet1[0] = "ACAC";
//    stringSet2[0] = "GTGT";
//    stringSet3 = stringSet1;
//    stringSet4 = stringSet2;
//
//    swap(stringSet1, stringSet2);
//    SEQAN_ASSERT(stringSet1 == stringSet4);
//    SEQAN_ASSERT(stringSet2 == stringSet3);
// }

// SEQAN_TYPED_TEST(StringSetTestCommon, Swap)
// {
//     CountingChar::clear();
// 
//     typename TestFixture::TStringSet strSet;
//     testStringSetSwap(strSet);
//     
//     typename TestFixture::TStringSet const constStrSet;
//     testStringSetSwap(constStrSet);
//     
//     testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
// }

// Test of value().
template <typename TStringSet>
void testStringSetValue(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    
    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.

    TString str1("ACAC");
    TString str2("AAAA");
    TString str3("TTTT");
    TStringSet stringSet;
    appendValue(stringSet, str1);
    appendValue(stringSet, str2);
    appendValue(stringSet, str3);
    TString & value_ = value(stringSet, 0);
    SEQAN_ASSERT_EQ(value_, str1);

    value_ = "GGGG";
    TString str4("GGGG");
    SEQAN_ASSERT_EQ(value_, str4);
    SEQAN_ASSERT_EQ(stringSet[0], str4);
}

// Test of value().
template <typename TStringSet>
void testStringSetValue(TStringSet const & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;
    
    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.

    TString str1("ACAC");
    TString str2("AAAA");
    TString str3("TTTT");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str1);
    appendValue(nonConstStringSet, str2);
    appendValue(nonConstStringSet, str3);

    TStringSet stringSet(nonConstStringSet);
    TString value_ = value(stringSet, 0);
    SEQAN_ASSERT_EQ(value_, str1);

    value_ = "GGGG";
    TString str4("GGGG");
    SEQAN_ASSERT_EQ(value_, str4);
    SEQAN_ASSERT_EQ(stringSet[0], str1);
}
// TODO(singer): Seg. fault
template <typename TValue>
void testStringSetValue(StringSet<String<TValue, MMap<> >, Owner<> > & /*Tag*/) {}
template <typename TValue>
void testStringSetValue(StringSet<String<TValue, MMap<> >, Owner<> > const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetValue(StringSet<String<TValue, MMap<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetValue(StringSet<String<TValue, MMap<> >, Dependent<TStringSetSpec> > const & /*Tag*/) {}


// TODO(singer)
template <typename TValue>
void testStringSetValue(StringSet<String<TValue, External<> >, Owner<> > & /*Tag*/) {}
template <typename TValue>
void testStringSetValue(StringSet<String<TValue, External<> >, Owner<> > const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetValue(StringSet<String<TValue, External<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetValue(StringSet<String<TValue, External<> >, Dependent<TStringSetSpec> > const & /*Tag*/) {}

template <typename TValue, typename TStringSpec>
void testStringSetValue(StringSet<String<TValue, TStringSpec>, Owner<ConcatDirect<> > > & /*Tag*/) {}
template <typename TValue, typename TStringSpec>
void testStringSetValue(StringSet<String<TValue, TStringSpec>, Owner<ConcatDirect<> > > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, Value)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetValue(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetValue(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}


// Test of valueById().
template <typename TStringSet>
void testStringSetValueById(TStringSet & /*Tag*/)
{
    using namespace seqan;

    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
    TString str1("ACAC");
    TString str2("AAAA");
    TString str3("TTTT");
    TStringSet stringSet;
    appendValue(stringSet, str1);
    appendValue(stringSet, str2);
    appendValue(stringSet, str3);
    TString & value_ = valueById(stringSet, 0);
    SEQAN_ASSERT_EQ(value_, str1);

    value_ = "GGGG";
    TString str4("GGGG");
    SEQAN_ASSERT_EQ(value_, str4);
    SEQAN_ASSERT_EQ(stringSet[0], str4);
}

// Test of valueById().
template <typename TStringSet>
void testStringSetValueById(TStringSet const & /*Tag*/)
{
    using namespace seqan;

    typedef typename RemoveConst<TStringSet>::Type TNonConstStringSet;
    typedef typename TestStringSetValue_<TStringSet>::Type TString;

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
  
    TString str1("ACAC");
    TString str2("AAAA");
    TString str3("TTTT");
    TNonConstStringSet nonConstStringSet;
    appendValue(nonConstStringSet, str1);
    appendValue(nonConstStringSet, str2);
    appendValue(nonConstStringSet, str3);
    
    TStringSet stringSet(nonConstStringSet);
    TString value_ = valueById(stringSet, 0);
    SEQAN_ASSERT_EQ(value_, str1);

    value_ = "GGGG";
    TString str4("GGGG");
    SEQAN_ASSERT_EQ(value_, str4);
    SEQAN_ASSERT_EQ(stringSet[0], str4);
}

// TODO(singer): Seg. fault
template <typename TValue>
void testStringSetValueById(StringSet<String<TValue, MMap<> >, Owner<> > & /*Tag*/) {}
template <typename TValue>
void testStringSetValueById(StringSet<String<TValue, MMap<> >, Owner<> > const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetValueById(StringSet<String<TValue, MMap<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetValueById(StringSet<String<TValue, MMap<> >, Dependent<TStringSetSpec> > const & /*Tag*/) {}


// TODO(singer)
template <typename TValue>
void testStringSetValueById(StringSet<String<TValue, External<> >, Owner<> > & /*Tag*/) {}
template <typename TValue>
void testStringSetValueById(StringSet<String<TValue, External<> >, Owner<> > const & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetValueById(StringSet<String<TValue, External<> >, Dependent<TStringSetSpec> > & /*Tag*/) {}
template <typename TValue, typename TStringSetSpec>
void testStringSetValueById(StringSet<String<TValue, External<> >, Dependent<TStringSetSpec> > const & /*Tag*/) {}

template <typename TValue, typename TStringSpec>
void testStringSetValueById(StringSet<String<TValue, TStringSpec>, Owner<ConcatDirect<> > > & /*Tag*/) {}
template <typename TValue, typename TStringSpec>
void testStringSetValueById(StringSet<String<TValue, TStringSpec>, Owner<ConcatDirect<> > > const & /*Tag*/) {}

SEQAN_TYPED_TEST(StringSetTestCommon, ValueById)
{
    CountingChar::clear();

    typename TestFixture::TStringSet strSet;
    testStringSetValueById(strSet);
    
    typename TestFixture::TStringSet const constStrSet;
    testStringSetValueById(constStrSet);
    
    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TStringSet>::Type>::Type());
}
#endif // CORE_TESTS_SEQUENCE_TEST_STRINGSET_H_
