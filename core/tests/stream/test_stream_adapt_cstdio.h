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
// Tests for FILE* to the stream concept adaption.
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_ADAPT_CSTDIO_H_
#define TEST_STREAM_TEST_STREAM_ADAPT_CSTDIO_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature<FILE *, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<FILE *, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<FILE *, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<FILE *, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<FILE *, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<FILE *, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<FILE *, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<FILE *, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple example of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_read_simple_usage)
{
    using namespace seqan;

    char const * STR = "This is a string!\nWith two lines.";
    FILE * stream = tmpfile();
    // LLVM warns about possible insecurities if we use STR as fmt string here.
    fprintf(stream, "%s", STR);
    rewind(stream);
    testStreamReadSimpleUsage(stream);
    fclose(stream);
}

// More complex example of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_read_complex_usage)
{
    using namespace seqan;

    char const * STR = "This is a string!\nWith two lines.";
    FILE * stream = tmpfile();
    // LLVM warns about possible insecurities if we use STR as fmt string here.
    fprintf(stream, "%s", STR);
    rewind(stream);
    testStreamReadComplexUsage(stream);
    fclose(stream);
}

// Simple example of reading from FILE *.
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_write_simple_usage)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    streamWriteChar(stream, '1');
    streamWriteChar(stream, '2');
    rewind(stream);

    char buffer[100];
    char * unused = fgets(buffer, 10, stream);
    (void) unused;
    SEQAN_ASSERT_EQ(strcmp(buffer, "12"), 0);
    SEQAN_ASSERT_EQ(streamTell(stream), 2u);

    fclose(stream);
}

// A bit more complex usage of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_write_complex_usage)
{
   using namespace seqan;
   
   FILE * stream = tmpfile();
   
   for (int i = 0; i < 10; ++i)
       streamWriteChar(stream, '0' + i);
   rewind(stream);

   char buffer[100];
   char * unused = fgets(buffer, 20, stream);
   (void) unused;
   SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
   SEQAN_ASSERT_EQ(streamTell(stream), 10u);

   fclose(stream);
}

// Test of streamEof().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_eof)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    fprintf(stream, "This is a test.");
    rewind(stream);
    testStreamEof(stream);
    fclose(stream);
}

// Test of streamPeek().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_peek)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    fprintf(stream, "This is a test.");
    rewind(stream);
    testStreamPeek(stream);
    fclose(stream);
}

// Test of streamReadChar().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_read_char)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    fprintf(stream, "123");
    rewind(stream);
    testStreamReadChar(stream);
    fclose(stream);
}

// Test of streamReadBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_read_block)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    fprintf(stream, "XXXXXXXXXX");
    rewind(stream);
    testStreamReadBlockHitLimit(stream);
    rewind(stream);
    testStreamReadBlockHitNoLimit(stream);
    fclose(stream);
}

// Test of streamWriteChar().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_write_char)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    testStreamWriteChar(stream);

    rewind(stream);
    char buffer[100];
    char * unused = fgets(buffer, 99, stream);
    (void) unused;
    SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);
    fclose(stream);
}

// Test of streamWriteBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_write_block)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    testStreamWriteBlock(stream);

    rewind(stream);
    char buffer[100];
    char * unused = fgets(buffer, 99, stream);
    (void) unused;
    SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);
    fclose(stream);
}

// Test of streamPut().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_streamPut)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    testStreamPut(stream);

    rewind(stream);
    char buffer[1000];
    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";

    int i = 0;

    while (!feof(stream) && i <998) buffer[i++] = fgetc(stream);
    buffer[i-1] = 0;
    SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);
    fclose(stream);
}

// Test of streamFlush().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_flush)
{
    using namespace seqan;
    // Only test that the function is there.

    FILE * stream = tmpfile();
    streamFlush(stream);
    fclose(stream);
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_seek)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    fprintf(stream, "0123456789");
    rewind(stream);
    testStreamSeek(stream);
    fclose(stream);
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_adapt_cstdio_tell)
{
    using namespace seqan;

    FILE * stream = tmpfile();
    fprintf(stream, "0123456789");
    rewind(stream);

    char c;
    size_t pos = streamTell(stream);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamReadChar(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamReadChar(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(stream);
    SEQAN_ASSERT_EQ(pos, 2u);

    fclose(stream);
}

#endif  // TEST_STREAM_TEST_STREAM_ADAPT_CSTDIO_H_
