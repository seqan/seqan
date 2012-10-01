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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for char arrays to the stream concept adaption.
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_CHAR_ARRAY_H_
#define TEST_STREAM_TEST_STREAM_CHAR_ARRAY_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

SEQAN_DEFINE_TEST(test_stream_char_array_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature<Stream<CharArray<char *> >, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
        b = HasStreamFeature<Stream<CharArray<char const *> >, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<CharArray<char *> >, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
        b = HasStreamFeature<Stream<CharArray<char const *> >, IsOutput>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<CharArray<char *> >, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
        b = HasStreamFeature<Stream<CharArray<char const *> >, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<CharArray<char *> >, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
        b = HasStreamFeature<Stream<CharArray<char const *> >, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<CharArray<char *> >, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
        b = HasStreamFeature<Stream<CharArray<char const *> >, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<CharArray<char *> >, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
        b = HasStreamFeature<Stream<CharArray<char const *> >, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<CharArray<char *> >, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
        b = HasStreamFeature<Stream<CharArray<char const *> >, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<CharArray<char *> >, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
        b = HasStreamFeature<Stream<CharArray<char const *> >, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple usage of pointers as stream for reading.
SEQAN_DEFINE_TEST(test_stream_char_array_read_simple_usage)
{
    using namespace seqan;

    char const * STR = "This is a string!\nWith two lines.";
    char buffer[100];
    strncpy(buffer, STR, 99);
    char * ptr = &(buffer[0]);
    char const * constPtr = ptr;

    Stream<CharArray<char *> > stream(ptr, ptr + strlen(STR));
    Stream<CharArray<char const *> > stream2(constPtr, constPtr + strlen(STR));

    testStreamReadSimpleUsage(stream);
    testStreamReadSimpleUsage(stream2);
}

// A bit more complex usage of pointers as strea for reading.
SEQAN_DEFINE_TEST(test_stream_char_array_read_complex_usage)
{
    using namespace seqan;

    char const * STR = "This is a string!\nWith two lines.";
    char buffer[100];
    strncpy(buffer, STR, 99);
    char * ptr = &(buffer[0]);
    char const * constPtr = ptr;

    Stream<CharArray<char *> > stream(ptr, ptr + strlen(STR));
    Stream<CharArray<char const *> > stream2(constPtr, constPtr + strlen(STR));

    testStreamReadComplexUsage(stream);
    testStreamReadComplexUsage(stream2);
}

// Simple usage of pointers as stream for writing.
SEQAN_DEFINE_TEST(test_stream_char_array_write_simple_usage)
{
    using namespace seqan;

    char buffer[100];
    char *ptr = &(buffer[0]);
    Stream<CharArray<char *> > stream(ptr, ptr + 100);

    streamWriteChar(stream, '1');
    streamWriteChar(stream, '2');
    streamWriteChar(stream, '\0');
    SEQAN_ASSERT_EQ(strcmp(buffer, "12"), 0);
    SEQAN_ASSERT_EQ(streamTell(stream), 3u);
}

// A bit more complex usage of pointers as strea for writing.
SEQAN_DEFINE_TEST(test_stream_char_array_write_complex_usage)
{
   using namespace seqan;
   
   char buffer[100];
   char *ptr = &(buffer[0]);
   Stream<CharArray<char *> > stream(ptr, ptr + 100);

   for (int i = 0; i < 10; ++i)
       streamWriteChar(stream, '0' + i);
   streamWriteChar(stream, '\0');
   SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
   SEQAN_ASSERT_EQ(streamTell(stream), 11u);
}

// Test of streamEof().
SEQAN_DEFINE_TEST(test_stream_char_array_eof)
{
    using namespace seqan;

    char const * STR = "This is a test.";
    char buffer[100];
    strcpy(buffer, STR);
    char * ptr = &(buffer[0]);
    char const * constPtr = ptr;

    Stream<CharArray<char *> > stream(ptr, ptr + strlen(STR));
    Stream<CharArray<char const *> > stream2(constPtr, constPtr + strlen(STR));

    testStreamEof(stream);
    testStreamEof(stream2);
}

// Test of streamPeek().
SEQAN_DEFINE_TEST(test_stream_char_array_peek)
{
    using namespace seqan;

    char const * STR = "This is a test.";
    char buffer[100];
    strncpy(buffer, STR, 100);
    char * ptr = &(buffer[0]);
    char const * constPtr = ptr;

    Stream<CharArray<char *> > stream(ptr, ptr + strlen(STR));
    Stream<CharArray<char const *> > stream2(constPtr, constPtr + strlen(STR));

    testStreamPeek(stream);
    testStreamPeek(stream2);
}

// Test of streamReadChar().
SEQAN_DEFINE_TEST(test_stream_char_array_read_char)
{
    using namespace seqan;

    char const * STR = "123";
    char buffer[100];
    strncpy(buffer, STR, 100);
    char * ptr = &(buffer[0]);
    char const * constPtr = ptr;

    Stream<CharArray<char *> > stream(ptr, ptr + strlen(STR));
    Stream<CharArray<char const *> > stream2(constPtr, constPtr + strlen(STR));

    testStreamReadChar(stream);
    testStreamReadChar(stream2);
}

// Test of streamReadBlock().
SEQAN_DEFINE_TEST(test_stream_char_array_read_block)
{
    using namespace seqan;

    char buffer[100];
    char * ptr = &(buffer[0]);
    for (int i = 0; i < 10; ++i, ++ptr)
        *ptr = 'X';
    *ptr = '\0';

    ptr = &(buffer[0]);
    char const * constPtr = ptr;

    Stream<CharArray<char *> > stream(ptr, ptr + 10);
    Stream<CharArray<char const *> > stream2(constPtr, constPtr + 10);

    testStreamReadBlockHitLimit(stream);
    testStreamReadBlockHitLimit(stream2);

    Stream<CharArray<char *> > stream3(ptr, ptr + 10);
    Stream<CharArray<char const *> > stream4(constPtr, constPtr + 10);

    testStreamReadBlockHitNoLimit(stream3);
    testStreamReadBlockHitNoLimit(stream4);
}

// Test of streamWriteChar().
SEQAN_DEFINE_TEST(test_stream_char_array_write_char)
{
    using namespace seqan;

    char buffer[100];
    char * ptr = &(buffer[0]);
    Stream<CharArray<char *> > stream(ptr, ptr + 100);

    testStreamWriteChar(stream);
    streamWriteChar(stream, '\0');
    SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);
}

// Test of streamWriteBlock().
SEQAN_DEFINE_TEST(test_stream_char_array_write_block)
{
    using namespace seqan;

    char buffer[100];
    char * ptr = &(buffer[0]);
    Stream<CharArray<char *> > stream(ptr, ptr + 100);

    testStreamWriteBlock(stream);
    buffer[8] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);
}


// Test of streamPut().
SEQAN_DEFINE_TEST(test_stream_char_array_streamPut)
{
    using namespace seqan;

    char buffer[100];
    char * ptr = &(buffer[0]);
    Stream<CharArray<char *> > stream(ptr, ptr + 100);

    testStreamPut(stream);
    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";
    buffer[strlen(cmp)] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);
}


// Test of streamFlush().
SEQAN_DEFINE_TEST(test_stream_char_array_flush)
{
    using namespace seqan;
    // Only test that the function is there.

    char const * STR = "123";
    char buffer[100];
    strncpy(buffer, STR, 100);
    char * ptr = &(buffer[0]);
    Stream<CharArray<char *> > stream(ptr, ptr + strlen(STR));

    streamFlush(stream);
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_char_array_seek)
{
    using namespace seqan;

    char const * STR = "0123456789";
    char buffer[100];
    strncpy(buffer, STR, 100);
    char * ptr = &(buffer[0]);
    char const * constPtr = ptr;

    Stream<CharArray<char *> > stream(ptr, ptr + strlen(STR));
    Stream<CharArray<char const *> > stream2(constPtr, constPtr + strlen(STR));

    testStreamSeek(stream);
    testStreamSeek(stream2);
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_char_array_tell)
{
    using namespace seqan;

    char const * STR = "0123456789";
    char buffer[100];
    strcpy(buffer, STR);
    char * ptr = &(buffer[0]);
    char const * constPtr = ptr;

    Stream<CharArray<char *> > stream(ptr, ptr + strlen(STR));
    Stream<CharArray<char const *> > stream2(constPtr, constPtr + strlen(STR));

    char c;
    size_t pos = streamTell(stream);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamReadChar(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamReadChar(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(stream);
    SEQAN_ASSERT_EQ(pos, 2u);
}

#endif  // TEST_STREAM_TEST_STREAM_CHAR_ARRAY_H_
