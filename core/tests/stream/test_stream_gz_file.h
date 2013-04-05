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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the GZ File Stream.
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_GZ_FILE_H_
#define TEST_STREAM_TEST_STREAM_GZ_FILE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

SEQAN_DEFINE_TEST(test_stream_gz_file_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature<Stream<GZFile>, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<GZFile>, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<GZFile>, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<GZFile>, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<GZFile>, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<GZFile>, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<GZFile>, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<GZFile>, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple example of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_gz_file_read_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    char const * STR = "This is a string!\nWith two lines.";
    gzFile gzOut = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(gzOut != NULL);
    gzwrite(gzOut, STR, strlen(STR));
    gzclose(gzOut);

    gzFile f = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    testStreamReadSimpleUsage(f2);
    gzclose(f);
}

// More complex example of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_gz_file_read_complex_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    char const * STR = "This is a string!\nWith two lines.";
    gzFile gzOut = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(gzOut != NULL);
    gzwrite(gzOut, STR, strlen(STR));
    gzclose(gzOut);

    gzFile f = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    testStreamReadComplexUsage(f2);
    gzclose(f);
}

// Simple example of reading from FILE *.
SEQAN_DEFINE_TEST(test_stream_gz_file_write_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    gzFile f = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    streamWriteChar(f2, '1');
    streamWriteChar(f2, '2');
    gzclose(f);

    // Read in data and compare.
    gzFile gzIn = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(gzIn != NULL);
    char buffer[100];
    int bytesRead = gzread(gzIn, buffer, 10);
    SEQAN_ASSERT_EQ(bytesRead, 2);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "12"), 0);
    gzclose(gzIn);
}

// A bit more complex usage of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_gz_file_write_complex_usage)
{
    using namespace seqan;
    
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);
   
    gzFile f = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    for (int i = 0; i < 10; ++i)
        streamWriteChar(f2, '0' + i);
    gzclose(f);

    // Read in data and compare.
    gzFile gzIn = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(gzIn != NULL);
    char buffer[100];
    int bytesRead = gzread(gzIn, buffer, 10);
    SEQAN_ASSERT_EQ(bytesRead, 10);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
    gzclose(gzIn);
}

// Test of streamEof().
SEQAN_DEFINE_TEST(test_stream_gz_file_eof)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    char const * STR = "This is a test.";
    gzFile gzOut = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(gzOut != NULL);
    gzwrite(gzOut, STR, strlen(STR));
    gzclose(gzOut);

    gzFile f = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    testStreamEof(f2);
    gzclose(f);
}

// Test of streamPeek().
SEQAN_DEFINE_TEST(test_stream_gz_file_peek)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    char const * STR = "This is a test.";
    gzFile gzOut = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(gzOut != NULL);
    gzwrite(gzOut, STR, strlen(STR));
    gzclose(gzOut);

    gzFile f = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    testStreamPeek(f2);
    gzclose(f);
}

// Test of streamReadChar().
SEQAN_DEFINE_TEST(test_stream_gz_file_read_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    char const * STR = "123";
    gzFile gzOut = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(gzOut != NULL);
    gzwrite(gzOut, STR, strlen(STR));
    gzclose(gzOut);

    gzFile f = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    testStreamReadChar(f2);
    gzclose(f);
}

// Test of streamReadBlock().
SEQAN_DEFINE_TEST(test_stream_gz_file_read_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    char const * STR = "XXXXXXXXXX";
    gzFile gzOut = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(gzOut != NULL);
    gzwrite(gzOut, STR, strlen(STR));
    gzclose(gzOut);

    {
        gzFile f = gzopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        Stream<GZFile> f2(f);
        testStreamReadBlockHitLimit(f2);
        gzclose(f);
    }
    {
        gzFile f = gzopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        Stream<GZFile> f2(f);
        testStreamReadBlockHitNoLimit(f2);
        gzclose(f);
    }
}

// Test of streamWriteChar().
SEQAN_DEFINE_TEST(test_stream_gz_file_write_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    gzFile f = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    testStreamWriteChar(f2);
    gzclose(f);

    // Read in data and compare.
    gzFile gzIn = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(gzIn != NULL);
    char buffer[100];
    int bytesRead = gzread(gzIn, buffer, 99);
    SEQAN_ASSERT_EQ(bytesRead, 3);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);
    gzclose(gzIn);
}

// Test of streamWrite().
SEQAN_DEFINE_TEST(test_stream_gz_file_write_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    gzFile f = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    testStreamWriteBlock(f2);
    gzclose(f);

    // Read in data and compare.
    gzFile gzIn = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(gzIn != NULL);
    char buffer[100];
    int bytesRead = gzread(gzIn, buffer, 99);
    SEQAN_ASSERT_EQ(bytesRead, 8);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);
    gzclose(gzIn);
}

// Test of streamWrite().
SEQAN_DEFINE_TEST(test_stream_gz_file_streamPut)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    gzFile f = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    testStreamPut(f2);
    gzclose(f);

    // Read in data and compare.
    gzFile gzIn = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(gzIn != NULL);
    char buffer[100];
    int bytesRead = gzread(gzIn, buffer, 99);
    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(bytesRead, int(sizeof(cmp) - sizeof(char)));
    SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);
    gzclose(gzIn);
}


// Test of streamFlush().
SEQAN_DEFINE_TEST(test_stream_gz_file_flush)
{
    using namespace seqan;

    // Only test that the function is there.
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    gzFile f = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    streamFlush(f2);
    gzclose(f);
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_gz_file_seek)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    char const * STR = "0123456789";
    gzFile gzOut = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(gzOut != NULL);
    gzwrite(gzOut, STR, strlen(STR));
    gzclose(gzOut);

    gzFile f = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);
    testStreamSeek(f2);
    gzclose(f);
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_gz_file_tell)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    char const * STR = "0123456789";
    gzFile gzOut = gzopen(filenameBuffer, "wb");
    SEQAN_ASSERT(gzOut != NULL);
    gzwrite(gzOut, STR, strlen(STR));
    gzclose(gzOut);

    gzFile f = gzopen(filenameBuffer, "rb");
    SEQAN_ASSERT(f != NULL);
    Stream<GZFile> f2(f);

    char c;
    size_t pos = streamTell(f2);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamReadChar(c, f2);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamReadChar(c, f2);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(f2);
    SEQAN_ASSERT_EQ(pos, 2u);

    gzclose(f);
}

#endif  // TEST_STREAM_TEST_STREAM_GZ_FILE_H_
