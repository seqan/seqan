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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
// This file contains functions to test the functionality of the sequence
// module.
// ==========================================================================

#ifndef TESTS_SEQUENCE_TEST_INFIX_H_
#define TESTS_SEQUENCE_TEST_INFIX_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include "test_sequence.h"

template <typename TAlphabetSpecPair>
class SegmentTest : public seqan::Test
{
public:
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 1>::Type TAlphabet;
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 2>::Type TStringSpec;
    typedef typename seqan::TagListValue<TAlphabetSpecPair, 3>::Type TSegmentSpec;
    typedef seqan::Segment<seqan::String<TAlphabet, TStringSpec>, TSegmentSpec> TSegment;
    typedef seqan::String<TAlphabet, TStringSpec> TString;
};

// TODO(singer): MMap causes lots of seg faults
typedef seqan::TagList<
//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
//             seqan::TagList<CountingChar, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
//             seqan::TagList<CountingChar, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
//             seqan::TagList<CountingChar, seqan::TagList<seqan::MMap<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<

            seqan::TagList<seqan::Dna, seqan::TagList<seqan::External<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::External<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::External<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::External<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::External<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::External<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::External<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::External<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::External<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::External<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::External<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::External<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<

//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
//             seqan::TagList<seqan::Dna, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
//             seqan::TagList<short, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::Packed<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<

            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Array<100>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<

            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Block<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<

// TODO(Singer): 7 errors and about 400 warnings (deprecated ...)
//             seqan::TagList<char, seqan::TagList<seqan::CStyle, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::CStyle, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
//             seqan::TagList<char, seqan::TagList<seqan::CStyle, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<

            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::InfixSegment> > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::PrefixSegment> > >, seqan::TagList<
            seqan::TagList<seqan::Dna, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<short, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<char, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::SuffixSegment> > >, seqan::TagList<
            seqan::TagList<CountingChar, seqan::TagList<seqan::Alloc<>, seqan::TagList<seqan::SuffixSegment> > >//, seqan::TagList<
//         > > > > > > > > > > > >
        > > > > > > > > > > > >
//         > > > > > > > > >
        > > > > > > > > > > > >
        > > > > > > > > > > > >
//         > > >
        > > > > > > > > > > > >
        SegmentTestTypes;

// TODO(singer): MMap causes lots of seg faults
typedef seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::External<> > >, seqan::TagList<
    seqan::TagList<char,        seqan::TagList<seqan::External<> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::External<> > >, seqan::TagList<
    seqan::TagList<CountingChar,seqan::TagList<seqan::External<> > >, seqan::TagList<
//     seqan::TagList<seqan::Dna5, seqan::TagList<seqan::MMap<> > >, seqan::TagList<
//     seqan::TagList<char,        seqan::TagList<seqan::MMap<> > >, seqan::TagList<
//     seqan::TagList<int,         seqan::TagList<seqan::MMap<> > >, seqan::TagList<
//     seqan::TagList<CountingChar,seqan::TagList<seqan::MMap<> > >, seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::Packed<> > >, seqan::TagList<
    seqan::TagList<char,        seqan::TagList<seqan::Packed<> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::Packed<> > >, seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::Array<100> > >, seqan::TagList<
    seqan::TagList<char,        seqan::TagList<seqan::Array<100> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::Array<100> > >, seqan::TagList<
    seqan::TagList<CountingChar,seqan::TagList<seqan::Array<100> > >, seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::Block<> > >, seqan::TagList<
    seqan::TagList<char,        seqan::TagList<seqan::Block<> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::Block<> > >, seqan::TagList<
    seqan::TagList<CountingChar,seqan::TagList<seqan::Block<> > >, seqan::TagList<
    seqan::TagList<seqan::Dna5, seqan::TagList<seqan::Alloc<> > >, seqan::TagList<
    seqan::TagList<char,        seqan::TagList<seqan::Alloc<> > >, seqan::TagList<
    seqan::TagList<int,         seqan::TagList<seqan::Alloc<> > >, seqan::TagList<
    seqan::TagList<CountingChar,seqan::TagList<seqan::Alloc<> > >, seqan::TagList<
    seqan::TagList<char,        seqan::TagList<seqan::CStyle> >
    > > > > > > > > > > > > > > > > > > > >// > > > >
    SegmentTestStringTypes;

template <typename T>
class SegmentTestCommon : public SegmentTest<T>
{};

template <typename T>
class SegmentTestString : public SegmentTest<T>
{};

SEQAN_TYPED_TEST_CASE(SegmentTestCommon, SegmentTestTypes);
SEQAN_TYPED_TEST_CASE(SegmentTestString, SegmentTestTypes);

// Test whether PrefixSegments are copy constructible.
template <typename TSegment>
void testPrefixConstructible(TSegment & /*Tag*/)
{
    using namespace seqan;
    typedef typename Host<TSegment>::Type TString;

    TString string("ACGTACGTACGT");

    // Test prefix() without start and end
    {
        TSegment prefixSeg(string);
        SEQAN_ASSERT_EQ(prefixSeg, string);
    }

    // Test prefix() with start and end
    {
        TSegment prefixSeg(string, 6);
        TString pre("ACGTAC");
        SEQAN_ASSERT_EQ(prefixSeg, pre);
    }

    // Test prefix() with an infix as host with start and end
    {
        TSegment prefixSeg1(string, 6);
        Segment<TSegment, PrefixSegment > prefixSeg2(prefixSeg1, 2);
        TString pre("AC");
        SEQAN_ASSERT_EQ(prefixSeg2, pre);
    }
}

// Test whether InfixSegments are copy constructible.
template <typename TSegment>
void testInfixConstructible(TSegment & /*Tag*/)
{
    using namespace seqan;
    typedef typename Host<TSegment>::Type TString;

    TString string("ACGTACGTACGT");

    // Test infix() without start and end
    {
        TSegment infixSeg(string);
        SEQAN_ASSERT(infixSeg == string);
    }

    // Test infix() with start and end
    {
        TSegment infixSeg(string, 4, 10);
        TString inf("ACGTAC");
        SEQAN_ASSERT(infixSeg == inf);
    }

    // Test infix() with an infix as host with start and end
    {
        TSegment infixSeg(string, 4, 10);
        Segment<TSegment> infixSeg2(infixSeg, 1, 3);
        TString inf("CG");
        SEQAN_ASSERT_EQ(infixSeg2, inf);
    }
}

// Test whether SuffixSegments are copy constructible.
template <typename TSegment>
void testSuffixConstructible(TSegment & /*Tag*/)
{
    using namespace seqan;
    typedef typename Host<TSegment>::Type TString;

    TString string("ACGTACGTACGT");

    // Test suffix() without start and end
    {
        TSegment suffixSeg(string);
        SEQAN_ASSERT_EQ(suffixSeg, string);
    }

    // Test suffix() with start and end
    {
        TSegment suffixSeg(string, 6);
        TString suf("GTACGT");
        SEQAN_ASSERT_EQ(suffixSeg, suf);
    }

    // Test suffix() with an infix as host with start and end
    {
        TSegment suffixSeg1(string, 6);
        Segment<TSegment, SuffixSegment > suffixSeg2(suffixSeg1, 4);
        TString suf("GT");
        SEQAN_ASSERT_EQ(suffixSeg2, suf);
    }
}

template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, MMap<> >, PrefixSegment> &) {}
template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, MMap<> >, PrefixSegment> const &) {}
template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, Packed<> >, PrefixSegment> &) {}
template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, Packed<> >, PrefixSegment> const &) {}

template <typename TString>
void testSegmentConstructible(Segment<TString, PrefixSegment> & seg)
{
    testPrefixConstructible(seg);
}

template <typename TString>
void testSegmentConstructible(Segment<TString, PrefixSegment> const & seg)
{
    testPrefixConstructible(seg);
}

template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, Packed<> >, InfixSegment> &) {}
template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, Packed<> >, InfixSegment> const &) {}

template <typename TString>
void testSegmentConstructible(Segment<TString, InfixSegment> & seg)
{
    testInfixConstructible(seg);
}

template <typename TString>
void testSegmentConstructible(Segment<TString, InfixSegment> const & seg)
{
    testInfixConstructible(seg);
}

template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, MMap<> >, SuffixSegment> &) {}
template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, MMap<> >, SuffixSegment> const &) {}
template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, Packed<> >, SuffixSegment> &) {}
template <typename TValue>
void testSegmentConstructible(Segment<String<TValue, Packed<> >, SuffixSegment> const &) {}


template <typename TString>
void testSegmentConstructible(Segment<TString, SuffixSegment> & seg)
{
    testSuffixConstructible(seg);
}

template <typename TString>
void testSegmentConstructible(Segment<TString, SuffixSegment> const & seg)
{
    testSuffixConstructible(seg);
}

// Test whether sequences are copy constructible.
SEQAN_TYPED_TEST(SegmentTestCommon, Constructible)
{
    CountingChar::clear();

    typename TestFixture::TSegment seg;
    testSegmentConstructible(seg);

    typename TestFixture::TSegment const constSeg;
    testSegmentConstructible(constSeg);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TSegment>::Type>::Type());
}


// TODO (singer): the constructor of segments has problems if 0 is one of the
// positions because 0 could be an iterator or a position.
// I would propose that segments are positional only.
// Test whether PrefixSegments are copy constructible.
template <typename TSegment>
void testPrefixCopyConstructible(TSegment & /*Tag*/)
{
    using namespace seqan;
    typedef typename Host<TSegment>::Type TString;

    TString string("ACGTACGTACGT");

    // Test copy constructible on empty template
    {
        TSegment prefixSeg(string);
        TSegment prefixSegCopy(prefixSeg);
        SEQAN_ASSERT_EQ(prefixSeg, prefixSegCopy);
    }

    // Test copy constructible on non empty template
    {
        TSegment prefixSeg(string, 6);
        TSegment prefixSegCopy(prefixSeg);
        TString pre("ACGTAC");
        SEQAN_ASSERT_EQ(prefixSeg, pre);
        SEQAN_ASSERT_EQ(prefixSeg, prefixSegCopy);
    }
}

// TODO (singer): the constructor of segments has problems if 0 is one of the
// positions because 0 could be an iterator or a position.
// I would propose that segments are positional only.
// Test whether InfixSegments are copy constructible.
template <typename TSegment>
void testInfixCopyConstructible(TSegment & /*Tag*/)
{
    using namespace seqan;
    typedef typename Host<TSegment>::Type TString;

    TString string("ACGTACGTACGT");

    // Test copy constructible on empty template
    {
        TSegment infixSeg(string);
        TSegment infixSegCopy(infixSeg);
        SEQAN_ASSERT_EQ(infixSeg, infixSegCopy);
    }

   // Test copy constructible on non empty template
    {
        TSegment infixSeg(string, 4, 10);
        TSegment infixSegCopy(infixSeg);
        TString inf("ACGTAC");
        SEQAN_ASSERT_EQ(infixSeg, inf);
        SEQAN_ASSERT_EQ(infixSeg, infixSegCopy);
    }
}

// TODO (singer): the constructor of segments has problems if 0 is one of the
// positions because 0 could be an iterator or a position.
// I would propose that segments are positional only.
// Test whether SuffixSegments are copy constructible.
template <typename TSegment>
void testSuffixCopyConstructible(TSegment & /*Tag*/)
{
    using namespace seqan;
    typedef typename Host<TSegment>::Type TString;

    TString string("ACGTACGTACGT");

    // Test copy constructible on empty template
    {
        TSegment suffixSeg(string);
        TSegment suffixSegCopy(suffixSeg);
        SEQAN_ASSERT_EQ(suffixSeg, suffixSegCopy);
    }

    // Test copy constructible on non empty template
    {
        TSegment suffixSeg(string, 6);
        TSegment suffixSegCopy(suffixSeg);
        TString suf("GTACGT");
        SEQAN_ASSERT_EQ(suffixSeg, suf);
        SEQAN_ASSERT_EQ(suffixSeg, suffixSegCopy);
    }
}

template <typename TString>
void testSegmentCopyConstructible(Segment<TString, PrefixSegment> & seg)
{
    testPrefixCopyConstructible(seg);
}

template <typename TString>
void testSegmentCopyConstructible(Segment<TString, PrefixSegment> const & seg)
{
    testPrefixConstructible(seg);
}

template <typename TValue>
void testSegmentCopyConstructible(Segment<String<TValue, MMap<> >, PrefixSegment> &) {}
template <typename TValue>
void testSegmentCopyConstructible(Segment<String<TValue, MMap<> >, PrefixSegment> const &) {}

template <typename TString>
void testSegmentCopyConstructible(Segment<TString, InfixSegment> & seg)
{
    testInfixCopyConstructible(seg);
}

template <typename TString>
void testSegmentCopyConstructible(Segment<TString, InfixSegment> const & seg)
{
    testInfixCopyConstructible(seg);
}

template <typename TValue>
void testSegmentCopyConstructible(Segment<String<TValue, MMap<> >, InfixSegment> &) {}
template <typename TValue>
void testSegmentCopyConstructible(Segment<String<TValue, MMap<> >, InfixSegment> const &) {}

template <typename TString>
void testSegmentCopyConstructible(Segment<TString, SuffixSegment> & seg)
{
    testSuffixCopyConstructible(seg);
}

template <typename TString>
void testSegmentCopyConstructible(Segment<TString, SuffixSegment> const & seg)
{
    testSuffixCopyConstructible(seg);
}

template <typename TValue>
void testSegmentCopyConstructible(Segment<String<TValue, MMap<> >, SuffixSegment> &) {}
template <typename TValue>
void testSegmentCopyConstructible(Segment<String<TValue, MMap<> >, SuffixSegment> const &) {}


// Test whether sequences are copy constructible.
SEQAN_TYPED_TEST(SegmentTestCommon, CopyConstructible)
{
    CountingChar::clear();

    typename TestFixture::TSegment seg;
    testSegmentCopyConstructible(seg);

    typename TestFixture::TSegment const constSeg;
    testSegmentCopyConstructible(constSeg);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TSegment>::Type>::Type());
}



// Test whether InfixSegments are default constructible.
template <typename TSegment>
void testSegmentDefaultConstructible(TSegment & /*Tag*/)
{
    using namespace seqan;
    typedef typename Host<TSegment>::Type TString;

    TSegment infixSeg;
    SEQAN_ASSERT(begin(infixSeg) == end(infixSeg));
}

template <typename TString>
void testSegmentDefaultConstructible(Segment<TString, PrefixSegment> & ) {}
template <typename TString>
void testSegmentDefaultConstructible(Segment<TString, PrefixSegment> const &) {}

template <typename TString>
void testSegmentDefaultConstructible(Segment<TString, InfixSegment> &) {}
template <typename TString>
void testSegmentDefaultConstructible(Segment<TString, InfixSegment> const &) {}

template <typename TString>
void testSegmentDefaultConstructible(Segment<TString, SuffixSegment> &) {}
template <typename TString>
void testSegmentDefaultConstructible(Segment<TString, SuffixSegment> const &) {}


// Test whether sequences are copy constructible.
SEQAN_TYPED_TEST(SegmentTestCommon, DefaultConstructible)
{
    CountingChar::clear();

    typename TestFixture::TSegment seg;
    testSegmentDefaultConstructible(seg);

    typename TestFixture::TSegment const constSeg;
    testSegmentDefaultConstructible(constSeg);

    testConstructDeconstruct(typename Value<typename Value<typename TestFixture::TSegment>::Type>::Type());
}

// Test operator<().
template <typename TString>
void testSegmentLess(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("ACGTACGTAAAA");

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 8);
        TInfix infixSeg(string, 5, 9);
        TSuffix suffixSeg(string, 3);

        SEQAN_ASSERT_EQ(prefixSeg1 < prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 < infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 < suffixSeg, true);
    }
    {
        TInfix infixSeg1(string, 4, 8);
        TInfix infixSeg2(string, 5, 9);
        TPrefix prefixSeg(string, 5);
        TSuffix suffixSeg(string, 3);

        SEQAN_ASSERT_EQ(infixSeg1 < infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 < prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 < suffixSeg, true);
    }
    {
        TSuffix suffixSeg1(string, 8);
        TSuffix suffixSeg2(string, 5);
        TPrefix prefixSeg(string, 5);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 < suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 < prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 < infixSeg, true);
    }
}

// TODO(singer): Seg fault
template <typename TValue>
void testSegmentLess(String<TValue, MMap<> > & /*Tag*/) {}
template <typename TValue>
void testSegmentLess(String<TValue, MMap<> > const & /*Tag*/) {}

SEQAN_TYPED_TEST(SegmentTestString, Less)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentLess(seg);

    typename TestFixture::TString const constSeg;
    testSegmentLess(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator<=().
template <typename TString>
void testSegmentLessEqual(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("ACGTACGTAAAA");

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 8);
        TInfix infixSeg(string, 4, 8);
        TSuffix suffixSeg(string, 4);

        SEQAN_ASSERT_EQ(prefixSeg1 <= prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 <= infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 <= suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 4, 8);
        TInfix infixSeg2(string, 0, 4);
        TPrefix prefixSeg(string, 5);
        TSuffix suffixSeg(string, 3);

        SEQAN_ASSERT_EQ(infixSeg1 <= infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 <= prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 <= suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 8);
        TSuffix suffixSeg2(string, 5);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 <= suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 <= prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 <= infixSeg, true);
    }
}

SEQAN_TYPED_TEST(SegmentTestString, LessEqual)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentLessEqual(seg);

    typename TestFixture::TString const constSeg;
    testSegmentLessEqual(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator>().
template <typename TString>
void testSegmentGreater(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("ACGTACGTAAAA");

    {
        TPrefix prefixSeg1(string, 5);
        TPrefix prefixSeg2(string, 4);
        TInfix infixSeg(string, 4, 8);
        TSuffix suffixSeg(string, 8);

        SEQAN_ASSERT_EQ(prefixSeg1 > prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 > infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 > suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 0, 12);
        TInfix infixSeg2(string, 0, 4);
        TPrefix prefixSeg(string, 5);
        TSuffix suffixSeg(string, 8);

        SEQAN_ASSERT_EQ(infixSeg1 > infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 > prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 > suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 4);
        TSuffix suffixSeg2(string, 9);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 > suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 > prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 > infixSeg, true);
    }
}

SEQAN_TYPED_TEST(SegmentTestString, Greater)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentGreater(seg);

    typename TestFixture::TString const constSeg;
    testSegmentGreater(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator>=().
template <typename TString>
void testSegmentGreaterEqual(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("ACGTACGTAAAA");

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 4);
        TInfix infixSeg(string, 4, 8);
        TSuffix suffixSeg(string, 8);

        SEQAN_ASSERT_EQ(prefixSeg1 >= prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 >= infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 >= suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 0, 12);
        TInfix infixSeg2(string, 0, 4);
        TPrefix prefixSeg(string, 5);
        TSuffix suffixSeg(string, 9);

        SEQAN_ASSERT_EQ(infixSeg1 >= infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 >= prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 >= suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 4);
        TSuffix suffixSeg2(string, 9);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 >= suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 >= prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 >= infixSeg, true);
    }
}

SEQAN_TYPED_TEST(SegmentTestString, GreaterEqual)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentGreaterEqual(seg);

    typename TestFixture::TString const constSeg;
    testSegmentGreaterEqual(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator==().
template <typename TString>
void testSegmentEqual(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("ACGTACGT");

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 4);
        TInfix infixSeg(string, 4, 8);
        TSuffix suffixSeg(string, 4);

        SEQAN_ASSERT_EQ(prefixSeg1 == prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 == infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 == suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 0, 4);
        TInfix infixSeg2(string, 4, 8);
        TPrefix prefixSeg(string, 4);
        TSuffix suffixSeg(string, 4);

        SEQAN_ASSERT_EQ(infixSeg1 == infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 == prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 == suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 4);
        TSuffix suffixSeg2(string, 4);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 == suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 == prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 == infixSeg, true);
    }
}

SEQAN_TYPED_TEST(SegmentTestString, Equal)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentEqual(seg);

    typename TestFixture::TString const constSeg;
    testSegmentEqual(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test operator!=().
template <typename TString>
void testSegmentUnequal(TString & /*Tag*/)
{
    using namespace seqan;

    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("ACGTACGT");

    {
        TPrefix prefixSeg1(string, 4);
        TPrefix prefixSeg2(string, 5);
        TInfix infixSeg(string, 4, 7);
        TSuffix suffixSeg(string, 3);

        SEQAN_ASSERT_EQ(prefixSeg1 != prefixSeg2, true);
        SEQAN_ASSERT_EQ(prefixSeg1 != infixSeg, true);
        SEQAN_ASSERT_EQ(prefixSeg1 != suffixSeg, true);
    }

    {
        TInfix infixSeg1(string, 1, 4);
        TInfix infixSeg2(string, 4, 8);
        TPrefix prefixSeg(string, 4);
        TSuffix suffixSeg(string, 4);

        SEQAN_ASSERT_EQ(infixSeg1 != infixSeg2, true);
        SEQAN_ASSERT_EQ(infixSeg1 != prefixSeg, true);
        SEQAN_ASSERT_EQ(infixSeg1 != suffixSeg, true);
    }

    {
        TSuffix suffixSeg1(string, 3);
        TSuffix suffixSeg2(string, 4);
        TPrefix prefixSeg(string, 4);
        TInfix infixSeg(string, 4, 8);

        SEQAN_ASSERT_EQ(suffixSeg1 != suffixSeg2, true);
        SEQAN_ASSERT_EQ(suffixSeg1 != prefixSeg, true);
        SEQAN_ASSERT_EQ(suffixSeg1 != infixSeg, true);
    }
}

SEQAN_TYPED_TEST(SegmentTestString, Unequal)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentUnequal(seg);

    typename TestFixture::TString const constSeg;
    testSegmentUnequal(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// // Test whether sequences are assignable.
// template <typename TString>
// void testSegmentAssignable(TString & /*Tag*/)
// {
//     using namespace seqan;
//
//     typedef Segment<TString, PrefixSegment> TPrefix;
//     typedef Segment<TString, InfixSegment> TInfix;
//     typedef Segment<TString, SuffixSegment> TSuffix;
//
//     TString string("ACGTACGT");
//
//     // Test assign() function.
//     {
//         // Test with default segment (the segment has no host).
//         {
//             TPrefix prefixSeg1;
//             TPrefix prefixSeg2(string, 5);
//             assign(prefixSeg1, prefixSeg2);
//             SEQAN_ASSERT_EQ(prefixSeg1, prefixSeg2);
//         }
//         {
//             TInfix infixSeg1;
//             TInfix infixSeg2(string, 0, 5);
//             assign(infixSeg1, infixSeg2);
//             SEQAN_ASSERT_EQ(infixSeg1, infixSeg2);
//
//             TInfix infixSeg3;
//             TPrefix prefixSeg(string, 5);
//             assign(infixSeg3, prefixSeg);
//             SEQAN_ASSERT_EQ(infixSeg3, prefixSeg);
//
//             TInfix infixSeg4;
//             TSuffix suffixSeg(string, 3);
//             assign(infixSeg4, suffixSeg);
//             SEQAN_ASSERT_EQ(infixSeg4, suffixSeg);
//         }
//         {
//             TSuffix suffixSeg1;
//             TSuffix suffixSeg2(string, 5);
//             assign(suffixSeg1, suffixSeg2);
//             SEQAN_ASSERT_EQ(suffixSeg1, suffixSeg2);
//         }
//     }
//     {
//         // Test with initialized segment (the segment has a host).
//         // In contrast to above the content of the host is changed.
//         // TODO (singer): the described behaviour is not desired.
//         {
//             TPrefix prefixSeg1(string, 4);
//             TPrefix prefixSeg2(string, 4);
//             assign(prefixSeg1, prefixSeg2);
//             SEQAN_ASSERT_EQ(prefixSeg1, prefixSeg2);
//         }
//         {
//             TInfix infixSeg1(string, 0, 4);
//             TInfix infixSeg2(string, 4, 8);
//             assign(infixSeg1, infixSeg2);
//             SEQAN_ASSERT_EQ(infixSeg1, infixSeg2);
//         }
//         {
//             TSuffix suffixSeg1(string, 4);
//             TSuffix suffixSeg2(string, 4);
//             assign(suffixSeg1, suffixSeg2);
//             SEQAN_ASSERT_EQ(suffixSeg1, suffixSeg2);
//         }
//     }
//
//     // Test operator=()
//     {
//         // Test with default segment (the segment has no host).
//         {
//             TPrefix prefixSeg1;
//             TPrefix prefixSeg2(string, 5);
//             prefixSeg1 = prefixSeg2;
//             SEQAN_ASSERT_EQ(prefixSeg1, prefixSeg2);
//         }
//         {
//             TInfix infixSeg1;
//             TInfix infixSeg2(string, 0, 5);
//             infixSeg1 = infixSeg2;
//             SEQAN_ASSERT_EQ(infixSeg1, infixSeg2);
//
//             TInfix infixSeg3;
//             TPrefix prefixSeg(string, 5);
//             infixSeg3 = prefixSeg;
//             SEQAN_ASSERT_EQ(infixSeg3, prefixSeg);
//
//             TInfix infixSeg4;
//             TSuffix suffixSeg(string, 3);
//             infixSeg4 = suffixSeg;
//             SEQAN_ASSERT_EQ(infixSeg4, suffixSeg);
//         }
//         {
//             TSuffix suffixSeg1;
//             TSuffix suffixSeg2(string, 5);
//             suffixSeg1 = suffixSeg2;
//             SEQAN_ASSERT_EQ(suffixSeg1, suffixSeg2);
//         }
//     }
//     {
//         // Test with initialized segment (the segment has a host).
//         // In contrast to above the content of the host is changed.
//         // TODO (singer): the described behaviour is not desired.
//         {
//             TPrefix prefixSeg1(string, 4);
//             TPrefix prefixSeg2(string, 4);
//             prefixSeg1 = prefixSeg2;
//             SEQAN_ASSERT_EQ(prefixSeg1, prefixSeg2);
//         }
//         {
//             TInfix infixSeg1(string, 0, 4);
//             TInfix infixSeg2(string, 4, 8);
//             infixSeg1 = infixSeg2;
//             SEQAN_ASSERT_EQ(infixSeg1, infixSeg2);
//         }
//         {
//             TSuffix suffixSeg1(string, 4);
//             TSuffix suffixSeg2(string, 4);
//             suffixSeg1 = suffixSeg2;
//             SEQAN_ASSERT_EQ(suffixSeg1, suffixSeg2);
//         }
//     }
// }
//
// // Test of assignValue().
// template <typename TString>
// void testSegmentAssignValue(TString & /*Tag*/)
// {
//     using namespace seqan;
//
//     typedef typename Value<TString>::Type TValue;
//     typedef Segment<TString, PrefixSegment> TPrefix;
//     typedef Segment<TString, InfixSegment> TInfix;
//     typedef Segment<TString, SuffixSegment> TSuffix;
//
//     TString string("ACGTACGTACGT");
//
//     TPrefix prefixSeg(string, 4);
//     TInfix infixSeg(string, 4, 8);
//     TSuffix suffixSeg(string, 8);
//
//     // TODO (singer): assignValue for segments not in docu.
//     assignValue(prefixSeg, 0, TValue('G'));
//     assignValue(infixSeg, 0, TValue('G'));
//     assignValue(suffixSeg, 0, TValue('G'));
//     SEQAN_ASSERT_EQ(string[0], TValue('G'));
//     SEQAN_ASSERT_EQ(string[4], TValue('G'));
//     SEQAN_ASSERT_EQ(string[8], TValue('G'));
// }
//
// // We need two back() tests, since back() returns a reference or a copy.
// // We check whether we can modify the reference.
// // Test of back() for non const strings.
// template <typename TString>
// void testSegmentBack(TString & /*Tag*/)
// {
//     using namespace seqan;
//
//     typedef typename Value<TString>::Type TValue;
//     typedef Segment<TString, PrefixSegment> TPrefix;
//     typedef Segment<TString, InfixSegment> TInfix;
//     typedef Segment<TString, SuffixSegment> TSuffix;
//
//     TString string("AAAAAAAAAA");
//
//     TPrefix prefixSeg(string, 3);
//     TInfix infixSeg(string, 3, 8);
//     TSuffix suffixSeg(string, 8);
//
//     // val is a reference in contrast to the const version of back().
//     TValue & preVal = back(prefixSeg);
//     preVal = 'C';
//     TValue & infVal = back(infixSeg);
//     infVal = 'C';
//     TValue & suffVal = back(suffixSeg);
//     suffVal = 'C';
//     SEQAN_ASSERT_EQ(string, "AACAAAACAC");
// }

// Test of back() for const strings.
template <typename TString>
void testSegmentBack(TString const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString const, PrefixSegment> TPrefix;
    typedef Segment<TString const, InfixSegment> TInfix;
    typedef Segment<TString const, SuffixSegment> TSuffix;

    TString string("AACAAAAGAT");

    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    // val is a reference in contrast to the const version of back().
    TValue preVal = back(prefixSeg);
    SEQAN_ASSERT_EQ(preVal, 'C');
    TValue infVal = back(infixSeg);
    SEQAN_ASSERT_EQ(infVal, 'G');
    TValue suffVal = back(suffixSeg);
    SEQAN_ASSERT_EQ(suffVal, 'T');
}

SEQAN_TYPED_TEST(SegmentTestString, Back)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentBack(seg);

    typename TestFixture::TString const constSeg;
    testSegmentBack(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}


// Test of begin().
template <typename TString>
void testSegmentBegin(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type const TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("CAAGAAAATA");
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    {
        typename Iterator<TPrefix>::Type it;

        it = begin(prefixSeg);
        SEQAN_ASSERT_EQ(getValue(it), TValue('C'));
        it = begin(prefixSeg, Standard());
        SEQAN_ASSERT_EQ(getValue(it), TValue('C'));
        it = begin(prefixSeg, Rooted());
        SEQAN_ASSERT_EQ(getValue(it), TValue('C'));
    }
    {
        typename Iterator<TInfix>::Type it;

        it = begin(infixSeg);
        SEQAN_ASSERT_EQ(getValue(it), TValue('G'));
        it = begin(infixSeg, Standard());
        SEQAN_ASSERT_EQ(getValue(it), TValue('G'));
        it = begin(infixSeg, Rooted());
        SEQAN_ASSERT_EQ(getValue(it), TValue('G'));
    }
    {
        typename Iterator<TSuffix>::Type it;

        it = begin(suffixSeg);
        SEQAN_ASSERT_EQ(getValue(it), TValue('T'));
        it = begin(suffixSeg, Standard());
        SEQAN_ASSERT_EQ(getValue(it), TValue('T'));
        it = begin(suffixSeg, Rooted());
        SEQAN_ASSERT_EQ(getValue(it), TValue('T'));
    }
}

SEQAN_TYPED_TEST(SegmentTestString, Begin)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentBegin(seg);

    typename TestFixture::TString const constSeg;
    testSegmentBegin(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of beginPosition().
template <typename TString>
void testSegmentBeginPosition(TString & /*Tag*/)
{
    using namespace seqan;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("CAAGAAAATA");
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    SEQAN_ASSERT_EQ(beginPosition(prefixSeg), 0u);
    SEQAN_ASSERT_EQ(beginPosition(infixSeg), 3u);
    SEQAN_ASSERT_EQ(beginPosition(suffixSeg), 8u);
}

SEQAN_TYPED_TEST(SegmentTestString, BeginPosition)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentBeginPosition(seg);

    typename TestFixture::TString const constSeg;
    testSegmentBeginPosition(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}
// TODO(singer): Clear uses assign, which we will replace with replace.
// Test of clear().
template <typename TString>
void testSegmentClear(TString & /*Tag*/)
{
//     using namespace seqan;
//     typedef Segment<TString, PrefixSegment> TPrefix;
//     typedef Segment<TString, InfixSegment> TInfix;
//     typedef Segment<TString, SuffixSegment> TSuffix;
//
//     // Since clear does not have to free anything it is only
//     // necessary that the code compiles.
//     // TODO (singer): at the moment clear erases parts of the host.
//     TString string("CAAGAAAATA");
//     TPrefix prefixSeg(string, 3);
//     clear(prefixSeg);
//     TInfix infixSeg(string, 3, 7);
//     clear(infixSeg);
//     TSuffix suffixSeg(string, 1);
//     clear(suffixSeg);
}
SEQAN_TYPED_TEST(SegmentTestString, Clear)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentClear(seg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of end().
template <typename TString>
void testSegmentEnd(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("AACAAAAGAT");
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    typename Iterator<TPrefix, Standard>::Type preIter = end(prefixSeg);
    typename Iterator<TPrefix, Standard>::Type preStandardIter = end(prefixSeg, Standard());
    typename Iterator<TPrefix, Rooted>::Type preRootedIter = end(prefixSeg, Rooted());
    --preIter;
    --preStandardIter;
    --preRootedIter;
    SEQAN_ASSERT_EQ(*preIter, TValue('C'));
    SEQAN_ASSERT_EQ(*preStandardIter, TValue('C'));
    SEQAN_ASSERT_EQ(*preRootedIter, TValue('C'));

    typename Iterator<TInfix, Standard>::Type infIter = end(infixSeg);
    typename Iterator<TInfix, Standard>::Type infStandardIter = end(infixSeg, Standard());
    typename Iterator<TInfix, Rooted>::Type infRootedIter = end(infixSeg, Rooted());
    --infIter;
    --infStandardIter;
    --infRootedIter;
    SEQAN_ASSERT_EQ(*infIter, TValue('G'));
    SEQAN_ASSERT_EQ(*infStandardIter, TValue('G'));
    SEQAN_ASSERT_EQ(*infRootedIter, TValue('G'));

    typename Iterator<TSuffix, Standard>::Type sufixIter = end(suffixSeg);
    typename Iterator<TSuffix, Standard>::Type suffStandardIter = end(suffixSeg, Standard());
    typename Iterator<TSuffix, Rooted>::Type suffRootedIter = end(suffixSeg, Rooted());
    --sufixIter;
    --suffStandardIter;
    --suffRootedIter;
    SEQAN_ASSERT_EQ(*sufixIter, TValue('T'));
    SEQAN_ASSERT_EQ(*suffStandardIter, TValue('T'));
    SEQAN_ASSERT_EQ(*suffRootedIter, TValue('T'));
}

SEQAN_TYPED_TEST(SegmentTestString, End)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentEnd(seg);

    typename TestFixture::TString const constSeg;
    testSegmentEnd(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}


// Test of endPosition().
template <typename TString>
void testSegmentEndPosition(TString & /*Tag*/)
{
    using namespace seqan;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("CAAGAAAATA");
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    SEQAN_ASSERT_EQ(endPosition(prefixSeg), 3u);
    SEQAN_ASSERT_EQ(endPosition(infixSeg), 8u);
    SEQAN_ASSERT_EQ(endPosition(suffixSeg), 10u);
}

SEQAN_TYPED_TEST(SegmentTestString, EndPosition)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentEndPosition(seg);

    typename TestFixture::TString const constSeg;
    testSegmentEndPosition(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// We need two front() tests, since back() returns a reference or a copy.
// We check whether we can modify the reference.
// Test of front() for non const strings.
template <typename TString>
void testSegmentFront(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("AAAAAAAAAA");

    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    // val is a reference in contrast to the const version of back().
    TValue & preVal = front(prefixSeg);
    preVal = 'C';
    TValue & infVal = front(infixSeg);
    infVal = 'C';
    TValue & suffVal = front(suffixSeg);
    suffVal = 'C';
    SEQAN_ASSERT_EQ(string, "CAACAAAACA");
}

// Test of back() for const strings.
template <typename TString>
void testSegmentFront(TString const & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString const, PrefixSegment> TPrefix;
    typedef Segment<TString const, InfixSegment> TInfix;
    typedef Segment<TString const, SuffixSegment> TSuffix;

    TString string("CAAGAAAATA");

    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    // val is a reference in contrast to the const version of back().
    TValue preVal = front(prefixSeg);
    SEQAN_ASSERT_EQ(preVal, 'C');
    TValue infVal = front(infixSeg);
    SEQAN_ASSERT_EQ(infVal, 'G');
    TValue suffVal = front(suffixSeg);
    SEQAN_ASSERT_EQ(suffVal, 'T');
}

// TODO(singer)
template <typename TValue>
void testSegmentFront(String<TValue, Packed<> > & /*Tag*/) {}
template <typename TValue>
void testSegmentFront(String<TValue, Packed<> > const & /*Tag*/) {}

SEQAN_TYPED_TEST(SegmentTestString, Front)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentFront(seg);

    typename TestFixture::TString const constSeg;
    testSegmentFront(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of getValue().
template <typename TString>
void testSegmentGetValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString const, PrefixSegment> TPrefix;
    typedef Segment<TString const, InfixSegment> TInfix;
    typedef Segment<TString const, SuffixSegment> TSuffix;

    TString string("CAAGAAAATA");

    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

     // In contrast to value(), getValue() does not return a reference but a copy.
    TValue preValue = getValue(prefixSeg, 0);
    TValue infValue = getValue(infixSeg, 0);
    TValue suffValue = getValue(suffixSeg, 0);

    SEQAN_ASSERT_EQ(preValue, TValue('C'));
    SEQAN_ASSERT_EQ(infValue, TValue('G'));
    SEQAN_ASSERT_EQ(suffValue, TValue('T'));
}

SEQAN_TYPED_TEST(SegmentTestString, GetValue)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentGetValue(seg);

    typename TestFixture::TString const constSeg;
    testSegmentGetValue(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of iter().
template <typename TString>
void testSegmentIter(TString & /*Tag*/)
{
    using namespace seqan;

    //typedef typename Value<TString>::Type TValue;

    TString string("ACGTACGTACGT");

    {
        typedef Segment<TString, PrefixSegment> TPrefix;
        typedef typename Iterator<TPrefix>::Type TIterator;
        typedef typename Iterator<TPrefix, Standard>::Type TStandardIterator;
        typedef typename Iterator<TPrefix, Rooted>::Type TRootedIterator;

        TPrefix prefixSeg(string, 4);

        TIterator iterator = iter(prefixSeg, 2);
        TStandardIterator standardIterator = iter(prefixSeg, 2);
        TRootedIterator rootedIterator = iter(prefixSeg, 2);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 2));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 2));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 2));
    }
    {
        typedef Segment<TString, InfixSegment> TInfix;
        typedef typename Iterator<TInfix>::Type TIterator;
        typedef typename Iterator<TInfix, Standard>::Type TStandardIterator;
        typedef typename Iterator<TInfix, Rooted>::Type TRootedIterator;

        TInfix infixSeg(string, 4, 8);

        TIterator iterator = iter(infixSeg, 2);
        TStandardIterator standardIterator = iter(infixSeg, 2);
        TRootedIterator rootedIterator = iter(infixSeg, 2);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 6));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 6));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 6));
    }
    {
        typedef Segment<TString, SuffixSegment> TSuffix;
        typedef typename Iterator<TSuffix>::Type TIterator;
        typedef typename Iterator<TSuffix, Standard>::Type TStandardIterator;
        typedef typename Iterator<TSuffix, Rooted>::Type TRootedIterator;

        TSuffix suffixSeg(string, 4);

        TIterator iterator = iter(suffixSeg, 2);
        TStandardIterator standardIterator = iter(suffixSeg, 2);
        TRootedIterator rootedIterator = iter(suffixSeg, 2);
        SEQAN_ASSERT_EQ(getValue(iterator), getValue(string, 6));
        SEQAN_ASSERT_EQ(getValue(standardIterator), getValue(string, 6));
        SEQAN_ASSERT_EQ(getValue(rootedIterator), getValue(string, 6));
    }
}

SEQAN_TYPED_TEST(SegmentTestString, Iter)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentIter(seg);

    typename TestFixture::TString const constSeg;
    testSegmentIter(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of length().
template <typename TString>
void testSegmentLength(TString & /*Tag*/)
{
    using namespace seqan;

    //typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("ACGTACGTACGT");
    TPrefix prefixSeg(string, 4);
    TInfix infixSeg(string, 4, 8);
    TSuffix suffixSeg(string, 8);
    SEQAN_ASSERT_EQ(length(infixSeg), 4u);
    SEQAN_ASSERT_EQ(length(infixSeg), 4u);
    SEQAN_ASSERT_EQ(length(infixSeg), 4u);
}

SEQAN_TYPED_TEST(SegmentTestString, Length)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentLength(seg);

    typename TestFixture::TString const constSeg;
    testSegmentLength(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// TODO(singer): Uses the assign which needs to be modified
// Test of MoveValue().
template <typename TString>
void testSegmentMoveValue(TString & /*Tag*/)
{
//     using namespace seqan;
//
//     typedef typename Value<TString>::Type TValue;
//     typedef Segment<TString, PrefixSegment> TPrefix;
//     typedef Segment<TString, InfixSegment> TInfix;
//     typedef Segment<TString, SuffixSegment> TSuffix;
//
//     TString string("AAAAAAAAAAAA");
//     TPrefix prefixSeg(string, 4);
//     TInfix infixSeg(string, 4, 8);
//     TSuffix suffixSeg(string, 8);
//
//     moveValue(prefixSeg, 3, 'C');
//     moveValue(infixSeg, 1, 'C');
//     moveValue(suffixSeg, 3, 'C');
//     SEQAN_ASSERT_EQ(string, "AAACACAAAAAC");
}

SEQAN_TYPED_TEST(SegmentTestString, Move)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentMoveValue(seg);

    typename TestFixture::TString const constSeg;
    testSegmentMoveValue(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}


// Test of replace().
template <typename TString>
void testSegmentReplace(TString & /*Tag*/)
{
    using namespace seqan;

    //typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    {
        // Segment do not overlap.
        TString string("AACCGGTT");
        TPrefix prefixSeg(string, 6);
        TInfix infixSeg(string, 4, 6);
        replace(prefixSeg, 0, 2, infixSeg);
        SEQAN_ASSERT_EQ(string, "GGCCGGTT");
    }
    {
        // Segment overlap.
        TString string("AACCGGTT");
        TPrefix prefixSeg(string, 6);
        TInfix infixSeg(string, 2, 6);
        replace(prefixSeg, 0, 4, infixSeg);
        SEQAN_ASSERT_EQ(string, "CCGGGGTT");
    }
    {
        // Segment do not overlap.
        TString string("AACCGGTT");
        TInfix infixSeg(string, 0, 2);
        TSuffix suffixSeg(string, 6);
        replace(infixSeg, 0, 2, suffixSeg);
        SEQAN_ASSERT_EQ(string, "TTCCGGTT");
    }
    {
        // Segment overlap.
        TString string("AACCGGTT");
        TInfix infixSeg(string, 0, 4);
        TSuffix suffixSeg(string, 6);
        replace(infixSeg, 0, 2, suffixSeg);
        SEQAN_ASSERT_EQ(string, "TTCCGGTT");
    }
    {
        // Segment do not overlap.
        TString string("AACCGGTT");
        TSuffix suffixSeg(string, 4);
        TPrefix prefixSeg(string, 2);
        replace(suffixSeg, 2, 4, prefixSeg);
        SEQAN_ASSERT_EQ(string, "AACCGGAA");
    }
    {
        // Segment overlap.
        TString string("AACCGGTT");
        TSuffix suffixSeg(string, 2);
        TPrefix prefixSeg(string, 4);
        replace(suffixSeg, 2, 6, prefixSeg);
        SEQAN_ASSERT_EQ(string, "AACCAACC");
    }
}
 // TODO(singer): error: no matching function for call to 'replace'
SEQAN_TYPED_TEST(SegmentTestString, Replace)
{
//     CountingChar::clear();
//
//     typename TestFixture::TString seg;
//     testSegmentReplace(seg);
//
//     typename TestFixture::TString const constSeg;
//     testSegmentReplace(constSeg);
//
//     testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of resize().
// TODO (singer): resize() is not defined for Segments (even though it
// appears in the docu)!
template <typename TString>
void testSegmentResize(TString & /*Tag*/)
{
    using namespace seqan;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("CAAGAAAATA");
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

//    resize(prefixSeg, 4);
//    resize(infixSeg, 4);
//    resize(suffixSeg, 4);
//    SEQAN_ASSERT_EQ(length(prefixSeg), 4u);
//    SEQAN_ASSERT_EQ(length(infixSeg), 4u);
//    SEQAN_ASSERT_EQ(length(suffixSeg), 4u);
}

SEQAN_TYPED_TEST(SegmentTestString, Resize)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentResize(seg);

    typename TestFixture::TString const constSeg;
    testSegmentResize(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of swap().
// TODO (singer): swap() is not defined for Segments!
// It appears in the docu, but it seams inappropriate
// for segments since it is not cleat how to deal with
// overlaps.
template <typename TString>
void testSegmentSwap(TString & /*Tag*/)
{
//    using namespace seqan;
//
//    typedef typename Value<TString>::Type TValue;
//    typedef Segment<TString, PrefixSegment> TPrefix;
//    typedef Segment<TString, InfixSegment> TInfix;
//    typedef Segment<TString, SuffixSegment> TSuffix;
//
//    {
//        TString string("AACCGGTT");
//        TPrefix prefixSeg(string, 4);
//        TInfix infixSeg(string, 4, 8);
//        swap(prefixSeg, infixSeg);
//        SEQAN_ASSERT_EQ(string, "GGTTAACC");
//    }
//    {
//        TString string("AACCGGTT");
//        TInfix infixSeg(string, 0, 4);
//        TSuffix suffixSeg(string, 4);
//        swap(infixSeg, suffixSeg);
//        SEQAN_ASSERT_EQ(string, "GGTTAACC");
//    }
//    {
//        TString string("AACCGGTT");
//        TSuffix suffixSeg(string, 4);
//        TPrefix prefixSeg(string, 4);
//        swap(suffixSeg, prefixSeg);
//        SEQAN_ASSERT_EQ(string, "GGTTAACC");
//    }
}

SEQAN_TYPED_TEST(SegmentTestString, Swap)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentSwap(seg);

    typename TestFixture::TString const constSeg;
    testSegmentSwap(constSeg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

// Test of value().
template <typename TString>
void testSegmentValue(TString & /*Tag*/)
{
    using namespace seqan;

    typedef typename Value<TString>::Type TValue;
    typedef Segment<TString, PrefixSegment> TPrefix;
    typedef Segment<TString, InfixSegment> TInfix;
    typedef Segment<TString, SuffixSegment> TSuffix;

    TString string("CAAGAAAATA");
    TPrefix prefixSeg(string, 3);
    TInfix infixSeg(string, 3, 8);
    TSuffix suffixSeg(string, 8);

    // In contrast to getValue(), value() does not return a copy but a reference.
    // We test this using the variable value_.
    TValue & preValue = value(prefixSeg, 0);
    SEQAN_ASSERT_EQ(preValue, 'C');
    TValue & infValue = value(infixSeg, 0);
    SEQAN_ASSERT_EQ(infValue, 'G');
    TValue & suffValue = value(suffixSeg, 0);
    SEQAN_ASSERT_EQ(suffValue, 'T');

    preValue = 'A';
    infValue = 'A';
    suffValue = 'A';
    SEQAN_ASSERT_EQ(string, "AAAAAAAAAA");
}

template <typename TValue>
void testSegmentValue(String<TValue, Packed<> > & /*Tag*/) {}

SEQAN_TYPED_TEST(SegmentTestString, Value)
{
    CountingChar::clear();

    typename TestFixture::TString seg;
    testSegmentValue(seg);

    testConstructDeconstruct(typename Value<typename TestFixture::TString>::Type());
}

//
// // --------------------------------------------------------------------------
// // Testing Alloc Strings With Simple Types
// // --------------------------------------------------------------------------
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_dna_constructible)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_char_constructible)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_short_constructible)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_counting_char_constructible)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_infix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_infix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_infix_short_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar, Array<100> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar,Array<100> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Block<> > string("A");
//         Segment<String<Dna, Block<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > string("A");
//         Segment<String<Dna, Block<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Block<> > string("A");
//         Segment<String<char, Block<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char,  Block<> > string("A");
//         Segment<String<char, Block<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Block<> > string("A");
//         Segment<String<short, Block<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<short,  Block<> > string("A");
//         Segment<String<short, Block<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc <> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, InfixSegment> tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, InfixSegment> const tag(string);
//         testInfixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_infix_char_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<char, CStyle> string("A");
// //         Segment<String<char, CStyle>, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, InfixSegment> tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> string("A");
// //         Segment<String<char, CStyle>, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, InfixSegment> const tag(string);
// //         testInfixConstructible(tag);
// //     }
// }
// //_______________________________________________________________________________
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_dna_constructible)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_char_constructible)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_short_constructible)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_counting_char_constructible)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_prefix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_prefix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_prefix_short_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar, Array<100> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar,Array<100> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Block<> > string("A");
//         Segment<String<Dna, Block<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > string("A");
//         Segment<String<Dna, Block<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Block<> > string("A");
//         Segment<String<char, Block<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char,  Block<> > string("A");
//         Segment<String<char, Block<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Block<> > string("A");
//         Segment<String<short, Block<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<short,  Block<> > string("A");
//         Segment<String<short, Block<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, PrefixSegment> tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, PrefixSegment> const tag(string);
//         testPrefixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_prefix_char_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<char, CStyle> string("A");
// //         Segment<String<char, CStyle>, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, PrefixSegment> tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> string("A");
// //         Segment<String<char, CStyle>, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, PrefixSegment> const tag(string);
// //         testPrefixConstructible(tag);
// //     }
// }
//
// //_______________________________________________________________________________
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_dna_constructible)
// {
// // TODO (singer): Segmenatation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_char_constructible)
// {
// // TODO (singer): Segmenatation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_short_constructible)
// {
// // TODO (singer): Segmenatation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_counting_char_constructible)
// {
// // TODO (singer): Segmenatation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_suffix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_suffix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_suffix_short_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar, Array<100> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar,Array<100> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Block<> > string("A");
//         Segment<String<Dna, Block<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > string("A");
//         Segment<String<Dna, Block<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Block<> > string("A");
//         Segment<String<char, Block<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char,  Block<> > string("A");
//         Segment<String<char, Block<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Block<> > string("A");
//         Segment<String<short, Block<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<short,  Block<> > string("A");
//         Segment<String<short, Block<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_dna_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_short_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc <> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_counting_char_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, SuffixSegment> tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, SuffixSegment> const tag(string);
//         testSuffixConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_suffix_char_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<char, CStyle> string("A");
// //         Segment<String<char, CStyle>, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, SuffixSegment> tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> string("A");
// //         Segment<String<char, CStyle>, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, SuffixSegment> const tag(string);
// //         testSuffixConstructible(tag);
// //     }
// }
//
// //_______________________________________________________________________________
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_dna_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_char_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_short_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_counting_char_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_infix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_infix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_infix_short_copy_constructible)
// {
// // TODO (singer): Shoult short work or not? Comparison with Proxy is problem.
// //     using namespace seqan;
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar, Array<100> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar,Array<100> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Block<> > string("A");
//         Segment<String<Dna, Block<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > string("A");
//         Segment<String<Dna, Block<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Block<> > string("A");
//         Segment<String<char, Block<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char,  Block<> > string("A");
//         Segment<String<char, Block<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Block<> > string("A");
//         Segment<String<short, Block<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<short,  Block<> > string("A");
//         Segment<String<short, Block<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc <> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, InfixSegment> tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, InfixSegment> const tag(string);
//         testInfixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_infix_char_copy_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<char, CStyle> string("A");
// //         Segment<String<char, CStyle>, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, InfixSegment> tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> string("A");
// //         Segment<String<char, CStyle>, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, InfixSegment> const tag(string);
// //         testInfixCopyConstructible(tag);
// //     }
// }
// //_______________________________________________________________________________
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_dna_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_char_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_short_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_counting_char_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_prefix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_prefix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_prefix_short_copy_constructible)
// {
// // TODO (singer): Should short work. Comparison with proxy is the problem.
// //     using namespace seqan;
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar, Array<100> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar,Array<100> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Block<> > string("A");
//         Segment<String<Dna, Block<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > string("A");
//         Segment<String<Dna, Block<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Block<> > string("A");
//         Segment<String<char, Block<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char,  Block<> > string("A");
//         Segment<String<char, Block<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Block<> > string("A");
//         Segment<String<short, Block<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<short,  Block<> > string("A");
//         Segment<String<short, Block<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, PrefixSegment> tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, PrefixSegment> const tag(string);
//         testPrefixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_prefix_char_copy_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<char, CStyle> string("A");
// //         Segment<String<char, CStyle>, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, PrefixSegment> tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> string("A");
// //         Segment<String<char, CStyle>, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, PrefixSegment> const tag(string);
// //         testPrefixCopyConstructible(tag);
// //     }
// }
//
// //_______________________________________________________________________________
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_dna_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_char_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_short_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_counting_char_copy_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_suffix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > string("A");
//         Segment<String<Dna, Packed<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Packed<> > const string("A");
//         Segment<String<Dna, Packed<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_suffix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char, Packed<> > string("A");
//         Segment<String<char, Packed<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char, Packed<> > const string("A");
//         Segment<String<char, Packed<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_suffix_short_copy_constructible)
// {
// // TODO (singer): Should short work. Comparison with proxy is the problem.
// //     using namespace seqan;
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar, Array<100> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar,Array<100> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Block<> > string("A");
//         Segment<String<Dna, Block<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > string("A");
//         Segment<String<Dna, Block<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Block<> > string("A");
//         Segment<String<char, Block<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char,  Block<> > string("A");
//         Segment<String<char, Block<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Block<> > string("A");
//         Segment<String<short, Block<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<short,  Block<> > string("A");
//         Segment<String<short, Block<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_dna_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_short_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Alloc<> > string("A");
//         Segment<String<Dna, Alloc<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > string("A");
//         Segment<String<Dna, Alloc <> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
//     {
//         String<Dna,  Alloc<> > const string("A");
//         Segment<String<Dna, Alloc<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_counting_char_copy_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, SuffixSegment> tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > string("A");
//         Segment<String<CountingChar, Alloc<> >, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Alloc<> > const string("A");
//         Segment<String<CountingChar, Alloc<> > const, SuffixSegment> const tag(string);
//         testSuffixCopyConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_suffix_char_copy_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<char, CStyle> string("A");
// //         Segment<String<char, CStyle>, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, SuffixSegment> tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> string("A");
// //         Segment<String<char, CStyle>, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, SuffixSegment> const tag(string);
// //         testSuffixCopyConstructible(tag);
// //     }
// }
//
// //_______________________________________________________________________________
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_dna_default_constructible)
// {
// // error: invalid operands to binary expression ('typename Iterator<Segment<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, InfixSegment>, typename DefaultGetIteratorSpec<Segment<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, InfixSegment> >::Type>::Type' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >') and 'typename Iterator<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, typename DefaultGetIteratorSpec<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > > >::Type>::Type' (aka 'ExtStringFwdIterator<String<seqan::SimpleType<unsigned char, seqan::Dna_>, External<seqan::ExternalConfig<seqan::File<seqan::Async<void> >, 4194304, 2> > > >'))
// //     using namespace seqan;
// //     {
// //         String<Dna, External<> > string("A");
// //         Segment<String<Dna, External<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, External<> > const string("A");
// //         Segment<String<Dna, External<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, External<> > string("A");
// //         Segment<String<Dna, External<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, External<> > const string("A");
// //         Segment<String<Dna, External<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_char_default_constructible)
// {
// // error: invalid operands to binary expression ('typename Iterator<Segment<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, InfixSegment>, typename DefaultGetIteratorSpec<Segment<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, InfixSegment> >::Type>::Type' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >') and 'typename Iterator<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, typename DefaultGetIteratorSpec<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > > >::Type>::Type' (aka 'ExtStringFwdIterator<String<seqan::SimpleType<unsigned char, seqan::Dna_>, External<seqan::ExternalConfig<seqan::File<seqan::Async<void> >, 4194304, 2> > > >'))
// //     using namespace seqan;
// //     {
// //         String<char, External<> > string("A");
// //         Segment<String<char, External<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, External<> > const string("A");
// //         Segment<String<char, External<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, External<> > string("A");
// //         Segment<String<char, External<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, External<> > const string("A");
// //         Segment<String<char, External<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_short_default_constructible)
// {
// // error: invalid operands to binary expression ('typename Iterator<Segment<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, InfixSegment>, typename DefaultGetIteratorSpec<Segment<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, InfixSegment> >::Type>::Type' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >') and 'typename Iterator<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, typename DefaultGetIteratorSpec<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > > >::Type>::Type' (aka 'ExtStringFwdIterator<String<seqan::SimpleType<unsigned char, seqan::Dna_>, External<seqan::ExternalConfig<seqan::File<seqan::Async<void> >, 4194304, 2> > > >'))
// //     using namespace seqan;
// //     {
// //         String<short, External<> > string("A");
// //         Segment<String<short, External<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, External<> > const string("A");
// //         Segment<String<short, External<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, External<> > string("A");
// //         Segment<String<short, External<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, External<> > const string("A");
// //         Segment<String<short, External<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_infix_counting_char_default_constructible)
// {
// // error: invalid operands to binary expression ('typename Iterator<Segment<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, InfixSegment>, typename DefaultGetIteratorSpec<Segment<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, InfixSegment> >::Type>::Type' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >') and 'typename Iterator<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > >, typename DefaultGetIteratorSpec<String<SimpleType<unsigned char, Dna_>, External<ExternalConfig<File<Async<void> >, 4194304, 2> > > >::Type>::Type' (aka 'ExtStringFwdIterator<String<seqan::SimpleType<unsigned char, seqan::Dna_>, External<seqan::ExternalConfig<seqan::File<seqan::Async<void> >, 4194304, 2> > > >'))
// //     using namespace seqan;
// //     {
// //         String<CountingChar, External<> > string("A");
// //         Segment<String<CountingChar, External<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, External<> > const string("A");
// //         Segment<String<CountingChar, External<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, External<> > string("A");
// //         Segment<String<CountingChar, External<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, External<> > const string("A");
// //         Segment<String<CountingChar, External<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_infix_counting_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_infix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Packed<> > string("A");
// //         Segment<String<Dna, Packed<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, Packed<> > const string("A");
// //         Segment<String<Dna, Packed<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no matching constructor for initialization of 'TIterator' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >')
// // //         String<Dna, Packed<> > string("A");
// // //         Segment<String<Dna, Packed<> >, InfixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, Packed<> > const string("A");
// //         Segment<String<Dna, Packed<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_infix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, Packed<> > string("A");
// //         Segment<String<char, Packed<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, Packed<> > const string("A");
// //         Segment<String<char, Packed<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no matching constructor for initialization of 'TIterator' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >')
// // //         String<char, Packed<> > string("A");
// // //         Segment<String<char, Packed<> >, InfixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, Packed<> > const string("A");
// //         Segment<String<char, Packed<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_infix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no matching constructor for initialization of 'TIterator' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >')
// // //         String<short, Packed<> > string("A");
// // //         Segment<String<short, Packed<> >, InfixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_dna_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_short_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_infix_counting_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar, Array<100> > const, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar,Array<100> > const, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_dna_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Block<> > string("A");
//         Segment<String<Dna, Block<> >, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         // TODO (singer):  error: no viable conversion from 'seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> > *')
// //         String<Dna,  Block<> > string("A");
// //         Segment<String<Dna, Block<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Block<> > string("A");
//         Segment<String<char, Block<> >, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         // TODO (singer): error: no viable conversion from 'seqan::String<char, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<char, seqan::Block<4096> > *')
// //         String<char,  Block<> > string("A");
// //         Segment<String<char, Block<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_short_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Block<> > string("A");
//         Segment<String<short, Block<> >, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         // TODO (singer): error: no viable conversion from 'seqan::String<short, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<short, seqan::Block<4096> > *')
// //         String<short,  Block<> > string("A");
// //         Segment<String<short, Block<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_infix_counting_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, InfixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         // TODO (singer): error: no viable conversion from 'seqan::String<CountingChar, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<CountingChar, seqan::Block<4096> > *')
// //         String<CountingChar,  Block<> > string("A");
// //         Segment<String<CountingChar, Block<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, InfixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > string("A");
// //         Segment<String<Dna, Alloc <> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_infix_counting_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, Alloc<> > string("A");
// //         Segment<String<CountingChar, Alloc<> >, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Alloc<> > const string("A");
// //         Segment<String<CountingChar, Alloc<> > const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Alloc<> > string("A");
// //         Segment<String<CountingChar, Alloc<> >, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Alloc<> > const string("A");
// //         Segment<String<CountingChar, Alloc<> > const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_infix_char_default_constructible)
// {
//     using namespace seqan;
// //     {
// //         String<char, CStyle> string("A");
// //         Segment<String<char, CStyle>, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, InfixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> string("A");
// //         Segment<String<char, CStyle>, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, InfixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
// //_______________________________________________________________________________
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_dna_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna, External<> > string("A");
//         Segment<String<Dna, External<> >, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna, External<> > const string("A");
//         Segment<String<Dna, External<> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char, External<> > string("A");
//         Segment<String<char, External<> >, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char, External<> > const string("A");
//         Segment<String<char, External<> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_short_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short, External<> > string("A");
//         Segment<String<short, External<> >, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short, External<> > const string("A");
//         Segment<String<short, External<> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_prefix_counting_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > string("A");
//         Segment<String<CountingChar, External<> >, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, External<> > const string("A");
//         Segment<String<CountingChar, External<> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_prefix_counting_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_prefix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Packed<> > string("A");
// //         Segment<String<Dna, Packed<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, Packed<> > const string("A");
// //         Segment<String<Dna, Packed<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no matching constructor for initialization of 'TIterator' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >')
// // //         String<Dna, Packed<> > string("A");
// // //         Segment<String<Dna, Packed<> >, PrefixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, Packed<> > const string("A");
// //         Segment<String<Dna, Packed<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_prefix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, Packed<> > string("A");
// //         Segment<String<char, Packed<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, Packed<> > const string("A");
// //         Segment<String<char, Packed<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no matching constructor for initialization of 'TIterator' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >')
// // //         String<char, Packed<> > string("A");
// // //         Segment<String<char, Packed<> >, PrefixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, Packed<> > const string("A");
// //         Segment<String<char, Packed<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_prefix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no matching constructor for initialization of 'TIterator' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >')
// // //         String<short, Packed<> > string("A");
// // //         Segment<String<short, Packed<> >, PrefixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_dna_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > string("A");
//         Segment<String<Dna, Array<100> >, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna, Array<100> > const string("A");
//         Segment<String<Dna, Array<100> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char, Array<100> > string("A");
//         Segment<String<char, Array<100> >, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char, Array<100> > const string("A");
//         Segment<String<char, Array<100> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_short_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short, Array<100> > string("A");
//         Segment<String<short, Array<100> >, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short, Array<100> > const string("A");
//         Segment<String<short, Array<100> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_prefix_counting_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar, Array<100> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > string("A");
//         Segment<String<CountingChar, Array<100> >, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar, Array<100> > const string("A");
//         Segment<String<CountingChar,Array<100> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_dna_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<Dna, Block<> > string("A");
//         Segment<String<Dna, Block<> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         // TODO (singer): error: no viable conversion from 'seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> > *')
// //         String<Dna,  Block<> > string("A");
// //         Segment<String<Dna, Block<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<Dna,  Block<> > const string("A");
//         Segment<String<Dna, Block<> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<char, Block<> > string("A");
//         Segment<String<char, Block<> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         // TODO (singer): error: no viable conversion from 'seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> > *')
// //         String<char,  Block<> > string("A");
// //         Segment<String<char, Block<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<char,  Block<> > const string("A");
//         Segment<String<char, Block<> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_short_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<short, Block<> > string("A");
//         Segment<String<short, Block<> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
//     {
//         // TODO (singer): error: no viable conversion from 'seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> > *')
// //         String<short,  Block<> > string("A");
// //         Segment<String<short, Block<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
//     }
//     {
//         String<short,  Block<> > const string("A");
//         Segment<String<short, Block<> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_prefix_counting_char_default_constructible)
// {
//     using namespace seqan;
//     {
//         String<CountingChar, Block<> > string("A");
//         Segment<String<CountingChar, Block<> >, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, PrefixSegment> tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         // TODO (singer): error: no viable conversion from 'seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> > *')
// //         String<CountingChar,  Block<> > string("A");
// //         Segment<String<CountingChar, Block<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
//     {
//         String<CountingChar,  Block<> > const string("A");
//         Segment<String<CountingChar, Block<> > const, PrefixSegment> const tag(string);
//         testSegmentDefaultConstructible(tag);
//         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_prefix_counting_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, Alloc<> > string("A");
// //         Segment<String<CountingChar, Alloc<> >, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Alloc<> > const string("A");
// //         Segment<String<CountingChar, Alloc<> > const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Alloc<> > string("A");
// //         Segment<String<CountingChar, Alloc<> >, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Alloc<> > const string("A");
// //         Segment<String<CountingChar, Alloc<> > const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_prefix_char_default_constructible)
// {
// //     using namespace seqan;
// //     {
// //         String<char, CStyle> string("A");
// //         Segment<String<char, CStyle>, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, PrefixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> string("A");
// //         Segment<String<char, CStyle>, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, PrefixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// //_______________________________________________________________________________
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, External<> > string("A");
// //         Segment<String<Dna, External<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, External<> > const string("A");
// //         Segment<String<Dna, External<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, External<> > string("A");
// //         Segment<String<Dna, External<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, External<> > const string("A");
// //         Segment<String<Dna, External<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, External<> > string("A");
// //         Segment<String<char, External<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, External<> > const string("A");
// //         Segment<String<char, External<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, External<> > string("A");
// //         Segment<String<char, External<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, External<> > const string("A");
// //         Segment<String<char, External<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, External<> > string("A");
// //         Segment<String<short, External<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, External<> > const string("A");
// //         Segment<String<short, External<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, External<> > string("A");
// //         Segment<String<short, External<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, External<> > const string("A");
// //         Segment<String<short, External<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_suffix_counting_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, External<> > string("A");
// //         Segment<String<CountingChar, External<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, External<> > const string("A");
// //         Segment<String<CountingChar, External<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, External<> > string("A");
// //         Segment<String<CountingChar, External<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, External<> > const string("A");
// //         Segment<String<CountingChar, External<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > string("A");
// //         Segment<String<Dna, MMap<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, MMap<> > const string("A");
// //         Segment<String<Dna, MMap<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > string("A");
// //         Segment<String<char, MMap<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, MMap<> > const string("A");
// //         Segment<String<char, MMap<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > string("A");
// //         Segment<String<short, MMap<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, MMap<> > const string("A");
// //         Segment<String<short, MMap<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_suffix_counting_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > string("A");
// //         Segment<String<CountingChar, MMap<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, MMap<> > const string("A");
// //         Segment<String<CountingChar, MMap<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_suffix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Packed<> > string("A");
// //         Segment<String<Dna, Packed<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, Packed<> > const string("A");
// //         Segment<String<Dna, Packed<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no matching constructor for initialization of 'TIterator' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >')
// // //         String<Dna, Packed<> > string("A");
// // //         Segment<String<Dna, Packed<> >, SuffixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, Packed<> > const string("A");
// //         Segment<String<Dna, Packed<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_suffix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, Packed<> > string("A");
// //         Segment<String<char, Packed<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, Packed<> > const string("A");
// //         Segment<String<char, Packed<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no matching constructor for initialization of 'TIterator' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >')
// // //         String<char, Packed<> > string("A");
// // //         Segment<String<char, Packed<> >, SuffixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, Packed<> > const string("A");
// //         Segment<String<char, Packed<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_suffix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, Packed<> > string("A");
// //         Segment<String<short, Packed<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no matching constructor for initialization of 'TIterator' (aka 'Iter<TSequence_, AdaptorIterator<TIterator_> >')
// // //         String<short, Packed<> > string("A");
// // //         Segment<String<short, Packed<> >, SuffixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, Packed<> > const string("A");
// //         Segment<String<short, Packed<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Array<100> > string("A");
// //         Segment<String<Dna, Array<100> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, Array<100> > const string("A");
// //         Segment<String<Dna, Array<100> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, Array<100> > string("A");
// //         Segment<String<Dna, Array<100> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna, Array<100> > const string("A");
// //         Segment<String<Dna, Array<100> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, Array<100> > string("A");
// //         Segment<String<char, Array<100> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, Array<100> > const string("A");
// //         Segment<String<char, Array<100> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, Array<100> > string("A");
// //         Segment<String<char, Array<100> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char, Array<100> > const string("A");
// //         Segment<String<char, Array<100> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, Array<100> > string("A");
// //         Segment<String<short, Array<100> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, Array<100> > const string("A");
// //         Segment<String<short, Array<100> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, Array<100> > string("A");
// //         Segment<String<short, Array<100> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short, Array<100> > const string("A");
// //         Segment<String<short, Array<100> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_suffix_counting_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, Array<100> > string("A");
// //         Segment<String<CountingChar, Array<100> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, Array<100> > const string("A");
// //         Segment<String<CountingChar, Array<100> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, Array<100> > string("A");
// //         Segment<String<CountingChar, Array<100> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar, Array<100> > const string("A");
// //         Segment<String<CountingChar,Array<100> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Block<> > string("A");
// //         Segment<String<Dna, Block<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Block<> > const string("A");
// //         Segment<String<Dna, Block<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no viable conversion from 'seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> > *')
// // //         String<Dna,  Block<> > string("A");
// // //         Segment<String<Dna, Block<> >, SuffixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Block<> > const string("A");
// //         Segment<String<Dna, Block<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<char, Block<> > string("A");
// //         Segment<String<char, Block<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char,  Block<> > const string("A");
// //         Segment<String<char, Block<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no viable conversion from 'seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> > *')
// // //         String<char,  Block<> > string("A");
// // //         Segment<String<char, Block<> >, SuffixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char,  Block<> > const string("A");
// //         Segment<String<char, Block<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<short, Block<> > string("A");
// //         Segment<String<short, Block<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short,  Block<> > const string("A");
// //         Segment<String<short, Block<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         // TODO (singer): error: no viable conversion from 'seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> > *')
// // //         String<short,  Block<> > string("A");
// // //         Segment<String<short, Block<> >, SuffixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<short,  Block<> > const string("A");
// //         Segment<String<short, Block<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_suffix_counting_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, Block<> > string("A");
// //         Segment<String<CountingChar, Block<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Block<> > const string("A");
// //         Segment<String<CountingChar, Block<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         // TODO (singer): error: no viable conversion from 'seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> >' to 'TContainerPointer_' (aka 'const seqan::String<seqan::SimpleType<unsigned char, seqan::Dna_>, seqan::Block<4096> > *')
// // //         String<CountingChar,  Block<> > string("A");
// // //         Segment<String<CountingChar, Block<> >, SuffixSegment> const tag(string);
// // //         testSegmentDefaultConstructible(tag);
// // //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Block<> > const string("A");
// //         Segment<String<CountingChar, Block<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_dna_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_short_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<Dna, Alloc<> > string("A");
// //         Segment<String<Dna, Alloc<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > string("A");
// //         Segment<String<Dna, Alloc <> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<Dna,  Alloc<> > const string("A");
// //         Segment<String<Dna, Alloc<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_suffix_counting_char_default_constructible)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     {
// //         String<CountingChar, Alloc<> > string("A");
// //         Segment<String<CountingChar, Alloc<> >, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Alloc<> > const string("A");
// //         Segment<String<CountingChar, Alloc<> > const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Alloc<> > string("A");
// //         Segment<String<CountingChar, Alloc<> >, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// //     {
// //         String<CountingChar,  Alloc<> > const string("A");
// //         Segment<String<CountingChar, Alloc<> > const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //         SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     }
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_suffix_char_default_constructible)
// {
//     using namespace seqan;
// //     {
// //         String<char, CStyle> string("A");
// //         Segment<String<char, CStyle>, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, SuffixSegment> tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char,  CStyle> string("A");
// //         Segment<String<char, CStyle>, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// //     {
// //         String<char,  CStyle> const string("A");
// //         Segment<String<char, CStyle> const, SuffixSegment> const tag(string);
// //         testSegmentDefaultConstructible(tag);
// //     }
// }
//
// //_______________________________________________________________________________
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_dna_less)
// {
//     using namespace seqan;
//     String<Dna, External<> > tag;
//     testSegmentLess(tag);
//
//     String<Dna, External<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_char_less)
// {
//     using namespace seqan;
//     String<char, External<> > tag;
//     testSegmentLess(tag);
//
//     String<char, External<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_short_less)
// {
//     using namespace seqan;
//     String<short, External<> > tag;
//     testSegmentLess(tag);
//
//     String<short, External<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_counting_char_less)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, External<> > tag;
//     testSegmentLess(tag);
//
//     String<CountingChar, External<> > const constTag;
//     testSegmentLess(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_dna_less)
// {
//     // TODO (singer): Segementation fault: 11
// //     using namespace seqan;
// //     String<Dna, MMap<> > tag;
// //     testSegmentLess(tag);
// //
// //     String<Dna, MMap<> > const constTag;
// //     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_char_less)
// {
//     // TODO (singer): Segementation fault: 11
// //     using namespace seqan;
// //     String<char, MMap<> > tag;
// //     testSegmentLess(tag);
// //
// //     String<char, MMap<> > const constTag;
// //     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_short_less)
// {
//     // TODO (singer): Segementation fault: 11
// //     using namespace seqan;
// //     String<short, MMap<> > tag;
// //     testSegmentLess(tag);
// //
// //     String<short, MMap<> > const constTag;
// //     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_counting_char_less)
// {
//     // TODO (singer): Segementation fault: 11
// //     using namespace seqan;
// //     CountingChar::clear();
// //
// //     String<CountingChar, MMap<> > tag;
// //     testSegmentLess(tag);
// //
// //     String<CountingChar, MMap<> > const constTag;
// //     testSegmentLess(constTag);
// //
// //     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_dna_less)
// {
//     using namespace seqan;
//     String<Dna, Packed<> > tag;
//     testSegmentLess(tag);
//
//     String<Dna, Packed<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_char_less)
// {
//     using namespace seqan;
//     String<char, Packed<> > tag;
//     testSegmentLess(tag);
//
//     String<char, Packed<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_short_less)
// {
//     using namespace seqan;
//     String<short, Packed<> > tag;
//     testSegmentLess(tag);
//
//     String<short, Packed<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_dna_less)
// {
//     using namespace seqan;
//     String<Dna, Array<100> > tag;
//     testSegmentLess(tag);
//
//     String<Dna, Array<100> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_char_less)
// {
//     using namespace seqan;
//     String<char, Array<100> > tag;
//     testSegmentLess(tag);
//
//     String<char, Array<100> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_short_less)
// {
//     using namespace seqan;
//     String<short, Array<100> > tag;
//     testSegmentLess(tag);
//
//     String<short, Array<100> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_counting_char_less)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Array<100> > tag;
//     testSegmentLess(tag);
//
//     String<CountingChar, Array<100> > const constTag;
//     testSegmentLess(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_dna_less)
// {
//     using namespace seqan;
//     String<Dna, Block<> > tag;
//     testSegmentLess(tag);
//
//     String<Dna, Block<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_char_less)
// {
//     using namespace seqan;
//     String<char, Block<> > tag;
//     testSegmentLess(tag);
//
//     String<char, Block<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_short_less)
// {
//     using namespace seqan;
//     String<short, Block<> > tag;
//     testSegmentLess(tag);
//
//     String<short, Block<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_counting_char_less)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Block<> > tag;
//     testSegmentLess(tag);
//
//     String<CountingChar, Block<> > const constTag;
//     testSegmentLess(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_less)
// {
//     using namespace seqan;
//     String<Dna, Alloc<> > tag;
//     testSegmentLess(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_char_less)
// {
//     using namespace seqan;
//     String<char, Alloc<> > tag;
//     testSegmentLess(tag);
//
//     String<char, Alloc<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_short_less)
// {
//     using namespace seqan;
//     String<short, Alloc<> > tag;
//     testSegmentLess(tag);
//
//     String<short, Alloc<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_counting_char_less)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Alloc<> > tag;
//     testSegmentLess(tag);
//
//     String<CountingChar, Alloc<> > const constTag;
//     testSegmentLess(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_char_less)
// {
// //     using namespace seqan;
// //     String<char, CStyle> tag;
// //     testSegmentLess(tag);
// //
// //     String<char, CStyle> const constTag;
// //     testSegmentLess(constTag);
// }
//
// // // Test operator<()
// // SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_less)
// // {
// //     using namespace seqan;
// //
// //     String<Dna, Alloc<> > tag;
// //     testSegmentLess(tag);
// //
// //     String<Dna, Alloc<> > const constTag;
// //     testSegmentLess(constTag);
// // }
//
// //################################################################################
// SEQAN_DEFINE_TEST(test_sequence_external_segment_dna_less_equal)
// {
//     using namespace seqan;
//     String<Dna, External<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<Dna, External<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_char_less_equal)
// {
//     using namespace seqan;
//     String<char, External<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<char, External<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_short_less_equal)
// {
//     using namespace seqan;
//     String<short, External<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<short, External<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_counting_char_less_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, External<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<CountingChar, External<> > const constTag;
//     testSegmentLessEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_dna_less_equal)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<Dna, MMap<> > tag;
// //     testSegmentLessEqual(tag);
// //
// //     String<Dna, MMap<> > const constTag;
// //     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_char_less_equal)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<char, MMap<> > tag;
// //     testSegmentLessEqual(tag);
// //
// //     String<char, MMap<> > const constTag;
// //     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_short_less_equal)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<short, MMap<> > tag;
// //     testSegmentLessEqual(tag);
// //
// //     String<short, MMap<> > const constTag;
// //     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_counting_char_less_equal)
// {
//     // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     CountingChar::clear();
// //
// //     String<CountingChar, MMap<> > tag;
// //     testSegmentLessEqual(tag);
// //
// //     String<CountingChar, MMap<> > const constTag;
// //     testSegmentLessEqual(constTag);
// //
// //     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_dna_less_equal)
// {
//     using namespace seqan;
//     String<Dna, Packed<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<Dna, Packed<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_char_less_equal)
// {
//     using namespace seqan;
//     String<char, Packed<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<char, Packed<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_short_less_equal)
// {
//     using namespace seqan;
//     String<short, Packed<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<short, Packed<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_dna_less_equal)
// {
//     using namespace seqan;
//     String<Dna, Array<100> > tag;
//     testSegmentLessEqual(tag);
//
//     String<Dna, Array<100> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_char_less_equal)
// {
//     using namespace seqan;
//     String<char, Array<100> > tag;
//     testSegmentLessEqual(tag);
//
//     String<char, Array<100> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_short_less_equal)
// {
//     using namespace seqan;
//     String<short, Array<100> > tag;
//     testSegmentLessEqual(tag);
//
//     String<short, Array<100> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_counting_char_less_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Array<100> > tag;
//     testSegmentLessEqual(tag);
//
//     String<CountingChar, Array<100> > const constTag;
//     testSegmentLessEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_dna_less_equal)
// {
//     using namespace seqan;
//     String<Dna, Block<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<Dna, Block<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_char_less_equal)
// {
//     using namespace seqan;
//     String<char, Block<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<char, Block<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_short_less_equal)
// {
//     using namespace seqan;
//     String<short, Block<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<short, Block<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_counting_char_less_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Block<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<CountingChar, Block<> > const constTag;
//     testSegmentLessEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_less_equal)
// {
//     using namespace seqan;
//     String<Dna, Alloc<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_char_less_equal)
// {
//     using namespace seqan;
//     String<char, Alloc<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<char, Alloc<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_short_less_equal)
// {
//     using namespace seqan;
//     String<short, Alloc<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<short, Alloc<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_counting_char_less_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Alloc<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<CountingChar, Alloc<> > const constTag;
//     testSegmentLessEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_char_less_equal)
// {
// //     using namespace seqan;
// //     String<char, CStyle> tag;
// //     testSegmentLessEqual(tag);
// //
// //     String<char, CStyle> const constTag;
// //     testSegmentLessEqual(constTag);
// }
//
// // // Test operator<=()
// // SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_less_equal)
// // {
// //     using namespace seqan;
// //
// //     String<Dna, Alloc<> > tag;
// //     testSegmentLessEqual(tag);
// //
// //     String<Dna, Alloc<> > const constTag;
// //     testSegmentLessEqual(constTag);
// // }
//
// // ########################################################
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_dna_greater)
// {
//     using namespace seqan;
//     String<Dna, External<> > tag;
//     testSegmentGreater(tag);
//
//     String<Dna, External<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_char_greater)
// {
//     using namespace seqan;
//     String<char, External<> > tag;
//     testSegmentGreater(tag);
//
//     String<char, External<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_short_greater)
// {
//     using namespace seqan;
//     String<short, External<> > tag;
//     testSegmentGreater(tag);
//
//     String<short, External<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_counting_char_greater)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, External<> > tag;
//     testSegmentGreater(tag);
//
//     String<CountingChar, External<> > const constTag;
//     testSegmentGreater(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_dna_greater)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<Dna, MMap<> > tag;
// //     testSegmentGreater(tag);
// //
// //     String<Dna, MMap<> > const constTag;
// //     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_char_greater)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<char, MMap<> > tag;
// //     testSegmentGreater(tag);
// //
// //     String<char, MMap<> > const constTag;
// //     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_short_greater)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<short, MMap<> > tag;
// //     testSegmentGreater(tag);
// //
// //     String<short, MMap<> > const constTag;
// //     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_counting_char_greater)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     CountingChar::clear();
// //
// //     String<CountingChar, MMap<> > tag;
// //     testSegmentGreater(tag);
// //
// //     String<CountingChar, MMap<> > const constTag;
// //     testSegmentGreater(constTag);
// //
// //     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_dna_greater)
// {
//     using namespace seqan;
//     String<Dna, Packed<> > tag;
//     testSegmentGreater(tag);
//
//     String<Dna, Packed<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_char_greater)
// {
//     using namespace seqan;
//     String<char, Packed<> > tag;
//     testSegmentGreater(tag);
//
//     String<char, Packed<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_short_greater)
// {
//     using namespace seqan;
//     String<short, Packed<> > tag;
//     testSegmentGreater(tag);
//
//     String<short, Packed<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_dna_greater)
// {
//     using namespace seqan;
//     String<Dna, Array<100> > tag;
//     testSegmentGreater(tag);
//
//     String<Dna, Array<100> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_char_greater)
// {
//     using namespace seqan;
//     String<char, Array<100> > tag;
//     testSegmentGreater(tag);
//
//     String<char, Array<100> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_short_greater)
// {
//     using namespace seqan;
//     String<short, Array<100> > tag;
//     testSegmentGreater(tag);
//
//     String<short, Array<100> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_counting_char_greater)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Array<100> > tag;
//     testSegmentGreater(tag);
//
//     String<CountingChar, Array<100> > const constTag;
//     testSegmentGreater(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_dna_greater)
// {
//     using namespace seqan;
//     String<Dna, Block<> > tag;
//     testSegmentGreater(tag);
//
//     String<Dna, Block<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_char_greater)
// {
//     using namespace seqan;
//     String<char, Block<> > tag;
//     testSegmentGreater(tag);
//
//     String<char, Block<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_short_greater)
// {
//     using namespace seqan;
//     String<short, Block<> > tag;
//     testSegmentGreater(tag);
//
//     String<short, Block<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_counting_char_greater)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Block<> > tag;
//     testSegmentGreater(tag);
//
//     String<CountingChar, Block<> > const constTag;
//     testSegmentGreater(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_greater)
// {
//     using namespace seqan;
//     String<Dna, Alloc<> > tag;
//     testSegmentGreater(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_char_greater)
// {
//     using namespace seqan;
//     String<char, Alloc<> > tag;
//     testSegmentGreater(tag);
//
//     String<char, Alloc<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_short_greater)
// {
//     using namespace seqan;
//     String<short, Alloc<> > tag;
//     testSegmentGreater(tag);
//
//     String<short, Alloc<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_counting_char_greater)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Alloc<> > tag;
//     testSegmentGreater(tag);
//
//     String<CountingChar, Alloc<> > const constTag;
//     testSegmentGreater(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_char_greater)
// {
// //     using namespace seqan;
// //     String<char, CStyle> tag;
// //     testSegmentGreater(tag);
// //
// //     String<char, CStyle> const constTag;
// //     testSegmentGreater(constTag);
// }
//
// // ########################################################
//
// // // Test operator>()
// // SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_greater)
// // {
// //     using namespace seqan;
// //
// //     String<Dna, Alloc<> > tag;
// //     testSegmentGreater(tag);
// //
// //     String<Dna, Alloc<> > const constTag;
// //     testSegmentGreater(constTag);
// // }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_dna_greater_equal)
// {
//     using namespace seqan;
//     String<Dna, External<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<Dna, External<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_char_greater_equal)
// {
//     using namespace seqan;
//     String<char, External<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<char, External<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_short_greater_equal)
// {
//     using namespace seqan;
//     String<short, External<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<short, External<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_counting_char_greater_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, External<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<CountingChar, External<> > const constTag;
//     testSegmentGreaterEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_dna_greater_equal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<Dna, MMap<> > tag;
// //     testSegmentGreaterEqual(tag);
// //
// //     String<Dna, MMap<> > const constTag;
// //     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_char_greater_equal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<char, MMap<> > tag;
// //     testSegmentGreaterEqual(tag);
// //
// //     String<char, MMap<> > const constTag;
// //     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_short_greater_equal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<short, MMap<> > tag;
// //     testSegmentGreaterEqual(tag);
// //
// //     String<short, MMap<> > const constTag;
// //     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_counting_char_greater_equal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     CountingChar::clear();
// //
// //     String<CountingChar, MMap<> > tag;
// //     testSegmentGreaterEqual(tag);
// //
// //     String<CountingChar, MMap<> > const constTag;
// //     testSegmentGreaterEqual(constTag);
// //
// //     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_dna_greater_equal)
// {
//     using namespace seqan;
//     String<Dna, Packed<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<Dna, Packed<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_char_greater_equal)
// {
//     using namespace seqan;
//     String<char, Packed<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<char, Packed<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_short_greater_equal)
// {
//     using namespace seqan;
//     String<short, Packed<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<short, Packed<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_dna_greater_equal)
// {
//     using namespace seqan;
//     String<Dna, Array<100> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<Dna, Array<100> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_char_greater_equal)
// {
//     using namespace seqan;
//     String<char, Array<100> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<char, Array<100> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_short_greater_equal)
// {
//     using namespace seqan;
//     String<short, Array<100> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<short, Array<100> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_counting_char_greater_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Array<100> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<CountingChar, Array<100> > const constTag;
//     testSegmentGreaterEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_dna_greater_equal)
// {
//     using namespace seqan;
//     String<Dna, Block<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<Dna, Block<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_char_greater_equal)
// {
//     using namespace seqan;
//     String<char, Block<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<char, Block<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_short_greater_equal)
// {
//     using namespace seqan;
//     String<short, Block<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<short, Block<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_counting_char_greater_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Block<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<CountingChar, Block<> > const constTag;
//     testSegmentGreaterEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_greater_equal)
// {
//     using namespace seqan;
//     String<Dna, Alloc<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_char_greater_equal)
// {
//     using namespace seqan;
//     String<char, Alloc<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<char, Alloc<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_short_greater_equal)
// {
//     using namespace seqan;
//     String<short, Alloc<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<short, Alloc<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_counting_char_greater_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Alloc<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<CountingChar, Alloc<> > const constTag;
//     testSegmentGreaterEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_char_greater_equal)
// {
// //     using namespace seqan;
// //     String<char, CStyle> tag;
// //     testSegmentGreaterEqual(tag);
// //
// //     String<char, CStyle> const constTag;
// //     testSegmentGreaterEqual(constTag);
// }
//
// // ########################################################
//
// // // Test operator>=()
// // SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_greater_equal)
// // {
// //     using namespace seqan;
// //
// //     String<Dna, Alloc<> > tag;
// //     testSegmentGreaterEqual(tag);
// //
// //     String<Dna, Alloc<> > const constTag;
// //     testSegmentGreaterEqual(constTag);
// // }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_dna_equal)
// {
//     using namespace seqan;
//     String<Dna, External<> > tag;
//     testSegmentEqual(tag);
//
//     String<Dna, External<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_char_equal)
// {
//     using namespace seqan;
//     String<char, External<> > tag;
//     testSegmentEqual(tag);
//
//     String<char, External<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_short_equal)
// {
//     using namespace seqan;
//     String<short, External<> > tag;
//     testSegmentEqual(tag);
//
//     String<short, External<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_counting_char_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, External<> > tag;
//     testSegmentEqual(tag);
//
//     String<CountingChar, External<> > const constTag;
//     testSegmentEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_dna_equal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<Dna, MMap<> > tag;
// //     testSegmentEqual(tag);
// //
// //     String<Dna, MMap<> > const constTag;
// //     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_char_equal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<char, MMap<> > tag;
// //     testSegmentEqual(tag);
// //
// //     String<char, MMap<> > const constTag;
// //     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_short_equal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<short, MMap<> > tag;
// //     testSegmentEqual(tag);
// //
// //     String<short, MMap<> > const constTag;
// //     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_counting_char_equal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     CountingChar::clear();
// //
// //     String<CountingChar, MMap<> > tag;
// //     testSegmentEqual(tag);
// //
// //     String<CountingChar, MMap<> > const constTag;
// //     testSegmentEqual(constTag);
// //
// //     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_dna_equal)
// {
//     using namespace seqan;
//     String<Dna, Packed<> > tag;
//     testSegmentEqual(tag);
//
//     String<Dna, Packed<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_char_equal)
// {
//     using namespace seqan;
//     String<char, Packed<> > tag;
//     testSegmentEqual(tag);
//
//     String<char, Packed<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_short_equal)
// {
//     using namespace seqan;
//     String<short, Packed<> > tag;
//     testSegmentEqual(tag);
//
//     String<short, Packed<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_dna_equal)
// {
//     using namespace seqan;
//     String<Dna, Array<100> > tag;
//     testSegmentEqual(tag);
//
//     String<Dna, Array<100> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_char_equal)
// {
//     using namespace seqan;
//     String<char, Array<100> > tag;
//     testSegmentEqual(tag);
//
//     String<char, Array<100> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_short_equal)
// {
//     using namespace seqan;
//     String<short, Array<100> > tag;
//     testSegmentEqual(tag);
//
//     String<short, Array<100> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_counting_char_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Array<100> > tag;
//     testSegmentEqual(tag);
//
//     String<CountingChar, Array<100> > const constTag;
//     testSegmentEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_dna_equal)
// {
//     using namespace seqan;
//     String<Dna, Block<> > tag;
//     testSegmentEqual(tag);
//
//     String<Dna, Block<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_char_equal)
// {
//     using namespace seqan;
//     String<char, Block<> > tag;
//     testSegmentEqual(tag);
//
//     String<char, Block<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_short_equal)
// {
//     using namespace seqan;
//     String<short, Block<> > tag;
//     testSegmentEqual(tag);
//
//     String<short, Block<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_counting_char_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Block<> > tag;
//     testSegmentEqual(tag);
//
//     String<CountingChar, Block<> > const constTag;
//     testSegmentEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_equal)
// {
//     using namespace seqan;
//     String<Dna, Alloc<> > tag;
//     testSegmentEqual(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_char_equal)
// {
//     using namespace seqan;
//     String<char, Alloc<> > tag;
//     testSegmentEqual(tag);
//
//     String<char, Alloc<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_short_equal)
// {
//     using namespace seqan;
//     String<short, Alloc<> > tag;
//     testSegmentEqual(tag);
//
//     String<short, Alloc<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_counting_char_equal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Alloc<> > tag;
//     testSegmentEqual(tag);
//
//     String<CountingChar, Alloc<> > const constTag;
//     testSegmentEqual(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_char_equal)
// {
// //     using namespace seqan;
// //     String<char, CStyle> tag;
// //     testSegmentEqual(tag);
// //
// //     String<char, CStyle> const constTag;
// //     testSegmentEqual(constTag);
// }
//
// // ########################################################
//
// // // Test operator==()
// // SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_equal)
// // {
// //     using namespace seqan;
// //
// //     String<Dna, Alloc<> > tag;
// //     testSegmentEqual(tag);
// //
// //     String<Dna, Alloc<> > const constTag;
// //     testSegmentEqual(constTag);
// // }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_dna_unequal)
// {
//     using namespace seqan;
//     String<Dna, External<> > tag;
//     testSegmentUnequal(tag);
//
//     String<Dna, External<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_char_unequal)
// {
//     using namespace seqan;
//     String<char, External<> > tag;
//     testSegmentUnequal(tag);
//
//     String<char, External<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_short_unequal)
// {
//     using namespace seqan;
//     String<short, External<> > tag;
//     testSegmentUnequal(tag);
//
//     String<short, External<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_external_segment_counting_char_unequal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, External<> > tag;
//     testSegmentUnequal(tag);
//
//     String<CountingChar, External<> > const constTag;
//     testSegmentUnequal(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_dna_unequal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<Dna, MMap<> > tag;
// //     testSegmentUnequal(tag);
// //
// //     String<Dna, MMap<> > const constTag;
// //     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_char_unequal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<char, MMap<> > tag;
// //     testSegmentUnequal(tag);
// //
// //     String<char, MMap<> > const constTag;
// //     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_short_unequal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     String<short, MMap<> > tag;
// //     testSegmentUnequal(tag);
// //
// //     String<short, MMap<> > const constTag;
// //     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_mmap_segment_counting_char_unequal)
// {
// // TODO (singer): Segmentation fault: 11
// //     using namespace seqan;
// //     CountingChar::clear();
// //
// //     String<CountingChar, MMap<> > tag;
// //     testSegmentUnequal(tag);
// //
// //     String<CountingChar, MMap<> > const constTag;
// //     testSegmentUnequal(constTag);
// //
// //     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
// //     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_dna_unequal)
// {
//     using namespace seqan;
//     String<Dna, Packed<> > tag;
//     testSegmentUnequal(tag);
//
//     String<Dna, Packed<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_char_unequal)
// {
//     using namespace seqan;
//     String<char, Packed<> > tag;
//     testSegmentUnequal(tag);
//
//     String<char, Packed<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_packed_segment_short_unequal)
// {
//     using namespace seqan;
//     String<short, Packed<> > tag;
//     testSegmentUnequal(tag);
//
//     String<short, Packed<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_dna_unequal)
// {
//     using namespace seqan;
//     String<Dna, Array<100> > tag;
//     testSegmentUnequal(tag);
//
//     String<Dna, Array<100> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_char_unequal)
// {
//     using namespace seqan;
//     String<char, Array<100> > tag;
//     testSegmentUnequal(tag);
//
//     String<char, Array<100> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_short_unequal)
// {
//     using namespace seqan;
//     String<short, Array<100> > tag;
//     testSegmentUnequal(tag);
//
//     String<short, Array<100> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_array_segment_counting_char_unequal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Array<100> > tag;
//     testSegmentUnequal(tag);
//
//     String<CountingChar, Array<100> > const constTag;
//     testSegmentUnequal(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_dna_unequal)
// {
//     using namespace seqan;
//     String<Dna, Block<> > tag;
//     testSegmentUnequal(tag);
//
//     String<Dna, Block<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_char_unequal)
// {
//     using namespace seqan;
//     String<char, Block<> > tag;
//     testSegmentUnequal(tag);
//
//     String<char, Block<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_short_unequal)
// {
//     using namespace seqan;
//     String<short, Block<> > tag;
//     testSegmentUnequal(tag);
//
//     String<short, Block<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_block_segment_counting_char_unequal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Block<> > tag;
//     testSegmentUnequal(tag);
//
//     String<CountingChar, Block<> > const constTag;
//     testSegmentUnequal(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_unequal)
// {
//     using namespace seqan;
//     String<Dna, Alloc<> > tag;
//     testSegmentUnequal(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_char_unequal)
// {
//     using namespace seqan;
//     String<char, Alloc<> > tag;
//     testSegmentUnequal(tag);
//
//     String<char, Alloc<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_short_unequal)
// {
//     using namespace seqan;
//     String<short, Alloc<> > tag;
//     testSegmentUnequal(tag);
//
//     String<short, Alloc<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_counting_char_unequal)
// {
//     using namespace seqan;
//     CountingChar::clear();
//
//     String<CountingChar, Alloc<> > tag;
//     testSegmentUnequal(tag);
//
//     String<CountingChar, Alloc<> > const constTag;
//     testSegmentUnequal(constTag);
//
//     SEQAN_ASSERT_EQ(CountingChar::numConstruct, CountingChar::numDeconstruct);
//     SEQAN_ASSERT_GT(CountingChar::numConstruct, 0u);
// }
//
// SEQAN_DEFINE_TEST(test_sequence_cstyle_segment_char_unequal)
// {
// //     using namespace seqan;
// //     String<char, CStyle> tag;
// //     testSegmentUnequal(tag);
// //
// //     String<char, CStyle> const constTag;
// //     testSegmentUnequal(constTag);
// }
//
// // ########################################################
//
// // // Test operator!=()
// // SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_unequal)
// // {
// //     using namespace seqan;
// //
// //     String<Dna, Alloc<> > tag;
// //     testSegmentUnequal(tag);
// //
// //     String<Dna, Alloc<> > const constTag;
// //     testSegmentUnequal(constTag);
// // }
//
// // Test whether sequences are assignable.
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_assignable)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentAssignable(tag);
// }
//
// // Test of assignValue().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_assign_value)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentAssignValue(tag);
// }
//
// // Test of back().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_back)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentBack(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentBack(constTag);
// }
//
// // Test of begin().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_begin)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentBegin(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentBegin(constTag);
// }
//
// // Test of beginPosition().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_begin_position)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentBeginPosition(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     //testSegmentBeginPosition(constTag);
// }
//
// // Test of clear().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_clear)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentClear(tag);
// }
//
// // Test of end().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_end)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentEnd(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentEnd(constTag);
// }
//
// // Test of endPosition().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_end_position)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentEndPosition(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentEndPosition(constTag);
// }
//
// // Test of front().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_front)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentFront(tag);
//
//     String<Dna, Alloc<> > const consTag;
//     testSegmentFront(consTag);
// }
//
// // Test of getValue().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_get_value)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentGetValue(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentGetValue(constTag);
// }
//
// // Test of iter().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_iter)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentIter(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentIter(constTag);
// }
//
// // Test of length().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_length)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentLength(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentLength(constTag);
// }
//
// // Test of moveValue().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_move_value)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentMoveValue(tag);
// }
//
// // Test of replace().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_replace)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentReplace(tag);
// }
//
// // Test of resize().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_resize)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentResize(tag);
// }
//
// // Test of swap().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_swap)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentSwap(tag);
// }
//
// // Test of value().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_value)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentValue(tag);
// }
//
// // TODO (singer): add tests to ensure that the underlying string of a segment is not modified.
//
// //        String<Dna, Packed<> > tag;
// //        testInfixConstructible(tag);
// //        testPrefixConstructible(tag);
// //        testSuffixConstructible(tag);
//     }
//
//     String<Dna, Alloc<> > const constTag;
//     testInfixConstructible(constTag);
//     testPrefixConstructible(constTag);
//     testSuffixConstructible(constTag);
// }
//
// // Test whether sequences are copy constructible.
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_copy_constructible)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testInfixCopyConstructible(tag);
//     testPrefixCopyConstructible(tag);
//     testSuffixCopyConstructible(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testInfixCopyConstructible(constTag);
//     testPrefixCopyConstructible(constTag);
//     testSuffixCopyConstructible(constTag);
// }
//
// // Test whether sequences are default constructible.
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_default_constructible)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testPrefixDefaultConstructible(tag);
//     testInfixDefaultConstructible(tag);
//     testSuffixDefaultConstructible(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testPrefixDefaultConstructible(constTag);
//     testInfixDefaultConstructible(constTag);
//     testSuffixDefaultConstructible(constTag);
// }
//
// // Test operator<()
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_less)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentLess(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentLess(constTag);
// }
//
// // Test operator<=()
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_less_equal)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentLessEqual(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentLessEqual(constTag);
// }
//
// // Test operator>()
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_greater)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentGreater(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentGreater(constTag);
// }
//
// // Test operator>=()
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_greater_equal)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentGreaterEqual(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentGreaterEqual(constTag);
// }
//
// // Test operator==()
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_equal)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentEqual(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentEqual(constTag);
// }
//
// // Test operator!=()
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_unequal)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentUnequal(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentUnequal(constTag);
// }
//
// // Test whether sequences are assignable.
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_assignable)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentAssignable(tag);
// }
//
// // Test of assignValue().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_assign_value)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentAssignValue(tag);
// }
//
// // Test of back().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_back)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentBack(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentBack(constTag);
// }
//
// // Test of begin().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_begin)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentBegin(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentBegin(constTag);
// }
//
// // Test of beginPosition().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_begin_position)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentBeginPosition(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     //testSegmentBeginPosition(constTag);
// }
//
// // // Test of clear().
// // SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_clear)
// // {
// //     using namespace seqan;
// //
// //     String<Dna, Alloc<> > tag;
// //     testSegmentClear(tag);
// // }
//
// // Test of end().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_end)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentEnd(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentEnd(constTag);
// }
//
// // Test of endPosition().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_end_position)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentEndPosition(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentEndPosition(constTag);
// }
//
// // Test of front().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_front)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentFront(tag);
//
//     String<Dna, Alloc<> > const consTag;
//     testSegmentFront(consTag);
// }
//
// // Test of getValue().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_get_value)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentGetValue(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentGetValue(constTag);
// }
//
// // Test of iter().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_iter)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentIter(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentIter(constTag);
// }
//
// // Test of length().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_length)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentLength(tag);
//
//     String<Dna, Alloc<> > const constTag;
//     testSegmentLength(constTag);
// }
//
// // Test of moveValue().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_move_value)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentMoveValue(tag);
// }
//
// // Test of value().
// SEQAN_DEFINE_TEST(test_sequence_alloc_segment_dna_value)
// {
//     using namespace seqan;
//
//     String<Dna, Alloc<> > tag;
//     testSegmentValue(tag);
// }

// TODO (singer): add tests to ensure that the underlying string of a segment is not modified.

#endif // TESTS_SEQUENCE_TEST_INFIX_H_
