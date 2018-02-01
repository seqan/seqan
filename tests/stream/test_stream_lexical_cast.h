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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// String <=> Numerical conversion tests
// ==========================================================================

#include <cmath>

#ifndef TEST_STREAM_TEST_LEXICAL_CAST_H_

using namespace seqan;

// ==========================================================================
// Types
// ==========================================================================

// --------------------------------------------------------------------------
// String Types
// --------------------------------------------------------------------------

typedef
    TagList<CharString,
    TagList<std::string
    > >
    StringTypes;

// --------------------------------------------------------------------------
// Signed Integer Types
// --------------------------------------------------------------------------

typedef
    TagList<short,
    TagList<int,
    TagList<long,
    TagList<long long
    > > > >
    SignedIntegerTypes;

// --------------------------------------------------------------------------
// Unsigned Integer Types
// --------------------------------------------------------------------------

typedef
    TagList<unsigned short,
    TagList<unsigned int,
    TagList<unsigned long,
    TagList<unsigned long long
    > > > >
    UnsignedIntegerTypes;

// --------------------------------------------------------------------------
// Floating Point Types
// --------------------------------------------------------------------------

typedef
    TagList<float,
    TagList<double
    > >
    FloatingPointTypes;

// --------------------------------------------------------------------------
// Aggregated Numerical Types
// --------------------------------------------------------------------------

typedef Sum<SignedIntegerTypes, UnsignedIntegerTypes>::Type    IntegerTypes;
typedef Sum<SignedIntegerTypes, FloatingPointTypes>::Type      SignedTypes;
typedef Sum<IntegerTypes, FloatingPointTypes>::Type            NumericalTypes;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class LexicalTest
// --------------------------------------------------------------------------

template <typename TTargetSourcePair>
class LexicalTest : public Test
{
public:
    typedef typename Value<TTargetSourcePair, 1>::Type  TTarget;
    typedef typename Value<TTargetSourcePair, 2>::Type  TSource;

    TTarget target;
    TSource source;
};

// --------------------------------------------------------------------------
// Class LexicalCastTest
// --------------------------------------------------------------------------

template <typename TTargetSourcePair>
class LexicalCastTest : public LexicalTest<TTargetSourcePair> {};

typedef Product<NumericalTypes, StringTypes>::Type LexicalCastTypes;

SEQAN_TYPED_TEST_CASE(LexicalCastTest, LexicalCastTypes);

// --------------------------------------------------------------------------
// Class AppendUnsignedTest
// --------------------------------------------------------------------------

template <typename TTargetSourcePair>
class AppendUnsignedTest : public LexicalTest<TTargetSourcePair> {};

typedef Product<CharString, NumericalTypes>::Type AppendNumericalTypes;

SEQAN_TYPED_TEST_CASE(AppendUnsignedTest, AppendNumericalTypes);

// --------------------------------------------------------------------------
// Class AppendSignedTest
// --------------------------------------------------------------------------

template <typename TTargetSourcePair>
class AppendSignedTest : public LexicalTest<TTargetSourcePair> {};

typedef Product<CharString, SignedTypes>::Type AppendSignedTypes;

SEQAN_TYPED_TEST_CASE(AppendSignedTest, AppendSignedTypes);

// --------------------------------------------------------------------------
// Class AppendFloatingPointTest
// --------------------------------------------------------------------------

template <typename TTargetSourcePair>
class AppendFloatingPointTest : public LexicalTest<TTargetSourcePair> {};

typedef Product<CharString, FloatingPointTypes>::Type AppendFloatingPointTypes;

SEQAN_TYPED_TEST_CASE(AppendFloatingPointTest, AppendFloatingPointTypes);

// ==========================================================================
// Tests
// ==========================================================================

// --------------------------------------------------------------------------
// Test lexicalCast(TTarget, UnsignedSource)
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(LexicalCastTest, UnsignedSource)
{
    assign(this->source, "12345");
    SEQAN_ASSERT(lexicalCast(this->target, this->source));
    SEQAN_ASSERT_EQ(this->target, static_cast<typename TestFixture::TTarget>(12345));
}

// --------------------------------------------------------------------------
// Test lexicalCast(TTarget, SignedSource)
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(LexicalCastTest, SignedSource)
{
    assign(this->source, "-12345");
    bool success = lexicalCast(this->target, this->source);
    SEQAN_ASSERT(success ^ IsUnsignedInteger<typename TestFixture::TTarget>::VALUE);
    if (success) SEQAN_ASSERT_EQ(this->target, static_cast<typename TestFixture::TTarget>(-12345));
}

// --------------------------------------------------------------------------
// Test lexicalCast(TTarget, FloatingPointSource)
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(LexicalCastTest, FloatingPointSource)
{
    typename TestFixture::TTarget reciprocal = (typename TestFixture::TTarget)123.45;
    typename TestFixture::TTarget epsilon = std::numeric_limits<typename TestFixture::TTarget>::epsilon();

    assign(this->source, "-123.45");
    bool success = lexicalCast(this->target, this->source);
    SEQAN_ASSERT(success ^ Is<IntegerConcept<typename TestFixture::TTarget> >::Type::VALUE);
    if (success) SEQAN_ASSERT_LT(this->target + reciprocal, epsilon);
}

// --------------------------------------------------------------------------
// Test lexicalCast(TTarget, WrongPrefixSource)
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(LexicalCastTest, WrongPrefixSource)
{
    assign(this->source, "foo123");
    SEQAN_ASSERT_NOT(lexicalCast(this->target, this->source));
}

// --------------------------------------------------------------------------
// Test lexicalCast(TTarget, WrongSuffixSource)
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(LexicalCastTest, WrongSuffixSource)
{
    assign(this->source, "123foo");
    SEQAN_ASSERT_NOT(lexicalCast(this->target, this->source));
}

// --------------------------------------------------------------------------
// Test lexicalCast<TTarget>(WrongSource)
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(LexicalCastTest, Exception)
{
    assign(this->source, "foo");
    try
    {
        this->target = lexicalCast<typename TestFixture::TTarget>(this->source);
    }
    catch (BadLexicalCast &)
    {
        return;
    }
    SEQAN_FAIL("The expected exception was not catched.");
}

// --------------------------------------------------------------------------
// Test appendNumber(TTarget, UnsignedSource)
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(AppendUnsignedTest, AppendNumber)
{
    assign(this->target, "foo");
    appendNumber(this->target, static_cast<typename TestFixture::TSource>(12345));
    SEQAN_ASSERT_EQ(this->target, "foo12345");
}

// --------------------------------------------------------------------------
// Test appendNumber(TTarget, SignedSource)
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(AppendSignedTest, AppendNumber)
{
    assign(this->target, "foo");
    appendNumber(this->target, static_cast<typename TestFixture::TSource>(-12345));
    SEQAN_ASSERT_EQ(this->target, "foo-12345");
}

// --------------------------------------------------------------------------
// Test appendNumber(TTarget, FloatingPointSource)
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(AppendFloatingPointTest, AppendNumber)
{
    assign(this->target, "foo");
    appendNumber(this->target, static_cast<typename TestFixture::TSource>(-123.45));
    SEQAN_ASSERT_EQ(this->target, "foo-123.45");
}

#endif // ifndef TEST_STREAM_TEST_LEXICAL_CAST_H_
