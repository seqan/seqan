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
// Tests for the BZ2 File Stream.
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_BZ2_FILE_H_
#define TEST_STREAM_TEST_STREAM_BZ2_FILE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

SEQAN_DEFINE_TEST(test_stream_bz2_file_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature<Stream<BZ2File>, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<BZ2File>, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<BZ2File>, HasPeek>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<BZ2File>, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<BZ2File>, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<BZ2File>, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<BZ2File>, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<BZ2File>, Tell>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
}

// Simple example of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_bz2_file_read_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    {
        char const * STR = "This is a string!\nWith two lines.";
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWrite(&err, f2, const_cast<char *>(STR), strlen(STR));
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    FILE * f = fopen(filenameBuffer, "rb");
    SEQAN_ASSERT(f != NULL);
    int err = BZ_OK;
    BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
    SEQAN_ASSERT_EQ(err, BZ_OK);
    Stream<BZ2File> f3(f2);
    testStreamReadSimpleUsage(f3);
    BZ2_bzReadClose(&err, f2);
    SEQAN_ASSERT_EQ(err, BZ_OK);
    fclose(f);
}

// More complex example of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_bz2_file_read_complex_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    {
        char const * STR = "This is a string!\nWith two lines.";
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWrite(&err, f2, const_cast<char *>(STR), strlen(STR));
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    FILE * f = fopen(filenameBuffer, "rb");
    SEQAN_ASSERT(f != NULL);
    int err = BZ_OK;
    BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
    SEQAN_ASSERT_EQ(err, BZ_OK);
    Stream<BZ2File> f3(f2);
    testStreamReadComplexUsage(f3);
    BZ2_bzReadClose(&err, f2);
    SEQAN_ASSERT_EQ(err, BZ_OK);
    fclose(f);
}

// Simple example of reading from FILE *.
SEQAN_DEFINE_TEST(test_stream_bz2_file_write_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Perform test with writing.
    {
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        Stream<BZ2File> f3(f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);

        streamWriteChar(f3, '1');
        streamWriteChar(f3, '2');

        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    // Compare resulting file.
    {
        FILE * f = fopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);

        char buffer[100];
        int bytesRead = BZ2_bzRead(&err, f2, buffer, 100);
        SEQAN_ASSERT_EQ(err, BZ_STREAM_END);
        SEQAN_ASSERT_EQ(bytesRead, 2);
        buffer[bytesRead] = '\0';
        SEQAN_ASSERT_EQ(strcmp(buffer, "12"), 0);

        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }
}

// A bit more complex usage of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_bz2_file_write_complex_usage)
{
    using namespace seqan;
    
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);
   
    // Perform test with writing.
    {
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        Stream<BZ2File> f3(f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);

        for (int i = 0; i < 10; ++i)
            streamWriteChar(f3, '0' + i);

        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    // Compare resulting file.
    {
        FILE * f = fopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);

        char buffer[100];
        int bytesRead = BZ2_bzRead(&err, f2, buffer, 100);
        SEQAN_ASSERT_EQ(err, BZ_STREAM_END);
        SEQAN_ASSERT_EQ(bytesRead, 10);
        buffer[bytesRead] = '\0';
        SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);

        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }
}

// Test of streamEof().
SEQAN_DEFINE_TEST(test_stream_bz2_file_eof)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    {
        char const * STR = "This is a test.";
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWrite(&err, f2, const_cast<char *>(STR), strlen(STR));
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    // Open file for reading and run test.
    {
        FILE * f = fopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        Stream<BZ2File> f3(f2);
        testStreamEof(f3);
        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }
}

// Test of streamReadChar().
SEQAN_DEFINE_TEST(test_stream_bz2_file_read_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    {
        char const * STR = "123";
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWrite(&err, f2, const_cast<char *>(STR), strlen(STR));
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    // Open file for reading and run test.
    {
        FILE * f = fopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        Stream<BZ2File> f3(f2);
        testStreamReadChar(f3);
        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }
}

// Test of streamReadBlock().
SEQAN_DEFINE_TEST(test_stream_bz2_file_read_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    {
        char const * STR = "XXXXXXXXXX";
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWrite(&err, f2, const_cast<char *>(STR), strlen(STR));
        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    // Open file for reading and run test.
    {
        FILE * f = fopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        Stream<BZ2File> f3(f2);
        testStreamReadBlockHitLimit(f3);
        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }
    {
        FILE * f = fopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        Stream<BZ2File> f3(f2);
        testStreamReadBlockHitNoLimit(f3);
        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }
}

// Test of streamWriteChar().
SEQAN_DEFINE_TEST(test_stream_bz2_file_write_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Run write test.
    {
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        Stream<BZ2File> f3(f2);

        testStreamWriteChar(f3);

        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    // Compare resulting file.
    {
        FILE * f = fopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);

        char buffer[100];
        int bytesRead = BZ2_bzRead(&err, f2, buffer, 100);
        SEQAN_ASSERT_EQ(bytesRead, 3);
        buffer[bytesRead] = '\0';
        SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);

        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }
}

// Test of streamWriteBlock().
SEQAN_DEFINE_TEST(test_stream_bz2_file_write_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Run write test.
    {
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        Stream<BZ2File> f3(f2);

        testStreamWriteBlock(f3);

        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    // Compare resulting file.
    {
        FILE * f = fopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);

        char buffer[100];
        int bytesRead = BZ2_bzRead(&err, f2, buffer, 99);
        SEQAN_ASSERT_EQ(bytesRead, 8);
        buffer[bytesRead] = '\0';
        SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);

        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }
}

// Test of streamPut()
SEQAN_DEFINE_TEST(test_stream_bz2_file_streamPut)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Run write test.
    {
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT_NOT(f == NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzWriteOpen(&err, f, 9, 0, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        Stream<BZ2File> f3(f2);

        testStreamPut(f3);

        SEQAN_ASSERT_EQ(err, BZ_OK);
        BZ2_bzWriteClose(&err, f2, 0, NULL, NULL);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

    // Compare resulting file.
    {
        FILE * f = fopen(filenameBuffer, "rb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);

        char buffer[100];
        char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";

        int bytesRead = BZ2_bzRead(&err, f2, buffer, 99);
        SEQAN_ASSERT_EQ(err, BZ_STREAM_END);
        buffer[bytesRead] = '\0';
        SEQAN_ASSERT_EQ(bytesRead, int(sizeof(cmp) - sizeof(char)));

        SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);

        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }

}

// Test of streamFlush().
SEQAN_DEFINE_TEST(test_stream_bz2_file_flush)
{
    using namespace seqan;

    // Only test that the function is there.
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Perform flush command
    {
        FILE * f = fopen(filenameBuffer, "wb");
        SEQAN_ASSERT(f != NULL);
        int err = BZ_OK;
        BZFILE * f2 = BZ2_bzReadOpen(&err, f, 0, 0, NULL, 0);
        SEQAN_ASSERT_EQ(err, BZ_OK);

        Stream<BZ2File> f3(f2);
        streamFlush(f3);

        BZ2_bzReadClose(&err, f2);
        SEQAN_ASSERT_EQ(err, BZ_OK);
        fclose(f);
    }
}

#endif  // TEST_STREAM_TEST_STREAM_BZ2_FILE_H_
