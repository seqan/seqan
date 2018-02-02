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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
//         David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Tests for seqan/stream/tokenize.h
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_TOKENIZATION_H_
#define TEST_STREAM_TEST_STREAM_TOKENIZATION_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

char const * EXAMPLE_STR1 = "This is a string ....foobar AAAACCCGGGTTT\n\
 TCG > foo.. SJAUDOF78456kapLP345LPL AAAAAAAAAAAAA\n\
\a 123gogogo \r\n\
Tetstetetstststetststetststetsstetetstetststetstdtetstestst    \r\n\
etstetetstststetststetststV           tststetstdtetsteststetstedetstet\r\n\
etstetetstststetststetststetsstetetstetststetstdtetsteststetstedetst|t\r\n\
\n\
\v\n\
\r\n\
ACGTACGTACGTACGATACGATCTnn\n\
\n\
ACGT*";

void createFile(CharString &fileName, char const *text)
{
    fileName = SEQAN_TEMP_FILENAME();
    std::fstream file(toCString(fileName), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file.is_open());
    file.write(text, strlen(text));
}

// ==========================================================================
// Types
// ==========================================================================

// --------------------------------------------------------------------------
// Input Stream Types
// --------------------------------------------------------------------------

typedef
    TagList<CharString,
    TagList<std::ifstream,
    TagList<std::fstream,
    TagList<Iter<std::ifstream, StreamIterator<Input> >
    > > > >
    TokenizationInputStreamTypes;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class TokenizationTest
// --------------------------------------------------------------------------

template <typename TStream_>
class TokenizationTest : public Test
{
public:
    typedef TStream_ TStream;
};

SEQAN_TYPED_TEST_CASE(TokenizationTest, TokenizationInputStreamTypes);

// --------------------------------------------------------------------------
// Class TokenizationContext
// --------------------------------------------------------------------------

template <typename TStream>
struct TokenizationContext
{};
//{
//    typedef typename Iterator<TStream, Standard>::Type  TIter;
//
//    TStream file;
//    TIter   iter;
//
//    TokenizationContext(const char *text)
//    {
//        CharString fileName;
//        createFile(fileName, text);
//        SEQAN_ASSERT(open(file, toCString(fileName)));
//        iter = TIter(file);
//    }
//};

template <>
struct TokenizationContext<std::fstream>
{
    std::fstream file;
    std::istreambuf_iterator<char> iter;

    TokenizationContext(const char *text)
    {
        CharString fileName;
        createFile(fileName, text);
        open(file, toCString(fileName));
        SEQAN_ASSERT(file.is_open());
        iter = std::istreambuf_iterator<char>(file);
    }
};

template <>
struct TokenizationContext<std::ifstream>
{
    std::ifstream file;
    std::istreambuf_iterator<char> iter;

    TokenizationContext(const char *text)
    {
        CharString fileName;
        createFile(fileName, text);
        open(file, toCString(fileName));
        SEQAN_ASSERT(file.is_open());
        iter = std::istreambuf_iterator<char>(file);
    }
};

template <typename TStream>
struct TokenizationContext<Iter<TStream, StreamIterator<Input> > >
{
    TStream file;
    Iter<TStream, StreamIterator<Input> > iter;

    TokenizationContext(const char *text)
    {
        CharString fileName;
        createFile(fileName, text);
        open(file, toCString(fileName));
        SEQAN_ASSERT(file.is_open());
        iter = Iter<TStream, StreamIterator<Input> >(file);
    }
};

// --------------------------------------------------------------------------
// Class TokenizationContext<String>
// --------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct TokenizationContext<String<TValue, TSpec> >
{
    String<TValue, TSpec> str;
    Iterator<CharString, Rooted>::Type iter;

    TokenizationContext(const char *text)
    {
        assign(str, text);
        iter = begin(str);
    }
};

// ==========================================================================
// Tests
// ==========================================================================

// readNChars
SEQAN_TYPED_TEST(TokenizationTest, read)
{
    TokenizationContext<typename TestFixture::TStream> ctx(EXAMPLE_STR1);
    SEQAN_ASSERT_NOT(atEnd(ctx.iter));

    CharString buf;

    SEQAN_ASSERT_EQ(read(buf, ctx.iter, 16), 16);
    SEQAN_ASSERT_EQ(buf, "This is a string");

    // check EOF-handling
    clear(buf);
    readUntil(buf, ctx.iter, EqualsChar<'*'>());
    SEQAN_ASSERT_EQ(*(ctx.iter), '*');

    clear(buf);
    SEQAN_ASSERT_EQ(read(buf,ctx.iter, 3), 1);
    SEQAN_ASSERT_EQ(buf, "*");
}


// _readHelper function with ignored characters
SEQAN_TYPED_TEST(TokenizationTest, ReadIgnoring)
{
    TokenizationContext<typename TestFixture::TStream> ctx("fefew AAAACC CGGGTT  TTCG ffe ATTAAaggtc agT");
    SEQAN_ASSERT_NOT(atEnd(ctx.iter));

    CharString buf;
    Dna5String buf2;

    // skip to where we want to go
    readUntil(buf, ctx.iter, IsInAlphabet<Dna5>());
    SEQAN_ASSERT_EQ(*(ctx.iter), 'A');

    readUntil(buf2, ctx.iter, NotFunctor<OrFunctor<IsInAlphabet<Dna5>, IsWhitespace> >(), IsWhitespace());
    SEQAN_ASSERT_EQ(buf2, "AAAACCCGGGTTTTCG");
    SEQAN_ASSERT_EQ(*(ctx.iter), 'f');

    readUntil(buf, ctx.iter, IsInAlphabet<Dna5>());
    SEQAN_ASSERT_EQ(*(ctx.iter), 'A');

    clear(buf2);
    readUntil(buf2, ctx.iter, EqualsChar<'_'>(), NotFunctor<IsInAlphabet<Dna5> >());
    SEQAN_ASSERT_EQ(buf2, "ATTAAAGGTCAGT");
    SEQAN_ASSERT(atEnd(ctx.iter));
}


SEQAN_TYPED_TEST(TokenizationTest, ReadUntil_ReadLine)
{
    TokenizationContext<typename TestFixture::TStream> ctx(EXAMPLE_STR1);
    SEQAN_ASSERT_NOT(atEnd(ctx.iter));

    CharString buf;
    // skip to where we want to go
    readUntil(buf, ctx.iter, EqualsChar<'A'>());
    SEQAN_ASSERT_EQ(*(ctx.iter), 'A');

    clear(buf);
    readLine(buf, ctx.iter);
    SEQAN_ASSERT_EQ(buf, "AAAACCCGGGTTT");
    SEQAN_ASSERT_EQ(*(ctx.iter), ' ');

    clear(buf);
    readLine(buf, ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), '\a');

    clear(buf);
    readUntil(buf, ctx.iter, EqualsChar<'1'>());

    clear(buf);
    readLine(buf, ctx.iter);
    SEQAN_ASSERT_EQ(buf, "123gogogo ");
    SEQAN_ASSERT_EQ(*(ctx.iter), 'T');

    // check EOF-handling
    clear(buf);
    readUntil(buf, ctx.iter, EqualsChar<'*'>());
    SEQAN_ASSERT_EQ(*(ctx.iter), '*');

    clear(buf);
    readLine(buf, ctx.iter);
    SEQAN_ASSERT_EQ(buf, "*");
    SEQAN_ASSERT(atEnd(ctx.iter));
}

// _skipHelper function "skipUntil"
SEQAN_TYPED_TEST(TokenizationTest, SkipUntil)
{
    TokenizationContext<typename TestFixture::TStream> ctx(EXAMPLE_STR1);
    SEQAN_ASSERT_NOT(atEnd(ctx.iter));

    skipUntil(ctx.iter, EqualsChar<'.'>());
    SEQAN_ASSERT_EQ(*(ctx.iter), '.');

    skipUntil(ctx.iter, IsBlank());
    SEQAN_ASSERT_EQ(*(ctx.iter), ' ');

    skipUntil(ctx.iter, EqualsChar<'A'>());
    SEQAN_ASSERT_EQ(*(ctx.iter), 'A');

    skipUntil(ctx.iter, IsWhitespace());
    SEQAN_ASSERT_EQ(*(ctx.iter), '\n');

    skipUntil(ctx.iter, EqualsChar<'A'>());
    SEQAN_ASSERT_EQ(*(ctx.iter), 'A');

    skipUntil(ctx.iter, IsBlank());
    SEQAN_ASSERT_EQ(*(ctx.iter), ' ');

    skipUntil(ctx.iter, EqualsChar<'A'>());
    SEQAN_ASSERT_EQ(*(ctx.iter), 'A');

    skipUntil(ctx.iter, IsWhitespace());
    SEQAN_ASSERT_EQ(*(ctx.iter), '\n');

    skipUntil(ctx.iter, IsGraph()); // skip over \a
    SEQAN_ASSERT_EQ(*(ctx.iter), '1');

    // check EOF-handling
    skipUntil(ctx.iter, EqualsChar<'*'>());
    SEQAN_ASSERT_EQ(*(ctx.iter), '*');

    skipUntil(ctx.iter, IsWhitespace());
    SEQAN_ASSERT(atEnd(ctx.iter));
}

// skipLine
SEQAN_TYPED_TEST(TokenizationTest, SkipLine)
{
    TokenizationContext<typename TestFixture::TStream> ctx(EXAMPLE_STR1);
    SEQAN_ASSERT_NOT(atEnd(ctx.iter));

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), ' ');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), '\a');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), 'T');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), 'e');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), 'e');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), '\n');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), '\v');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), '\r');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), 'A');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), '\n');

    skipLine(ctx.iter);
    SEQAN_ASSERT_EQ(*(ctx.iter), 'A');

    skipLine(ctx.iter);
    SEQAN_ASSERT(atEnd(ctx.iter));
}

#endif // ifndef TEST_STREAM_TEST_STREAM_TOKENIZATION_H_
