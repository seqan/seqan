// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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

#ifndef CORE_TESTS_STREAM_TEST_STREAM_BGZF_H_
#define CORE_TESTS_STREAM_TEST_STREAM_BGZF_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

SEQAN_DEFINE_TEST(test_stream_bgzf_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature<Stream<Bgzf>, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<Bgzf>, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<Bgzf>, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<Bgzf>, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<Bgzf>, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<Bgzf>, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<Bgzf>, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<Bgzf>, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple example of reading from BGZF Stream.
SEQAN_DEFINE_TEST(test_stream_bgzf_read_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    char const * STR = "This is a string!\nWith two lines.";
    streamWriteBlock(bgzfOut, STR, strlen(STR));
    close(bgzfOut);

    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));
    SEQAN_ASSERT_EQ(f._fileSize, 87);
    testStreamReadSimpleUsage(f);
    close(f);
}

// More complex example of reading from BGZF Stream.
SEQAN_DEFINE_TEST(test_stream_bgzf_read_complex_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    char const * STR = "This is a string!\nWith two lines.";
    streamWriteBlock(bgzfOut, STR, strlen(STR));
    close(bgzfOut);

    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));
    SEQAN_ASSERT_EQ(f._fileSize, 87);
    testStreamReadComplexUsage(f);
    close(f);
}

// Simple example of reading from FILE *.
SEQAN_DEFINE_TEST(test_stream_bgzf_write_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "wb"));
    streamWriteChar(f, '1');
    streamWriteChar(f, '2');
    close(f);

    // Read in data and compare.
    Stream<Bgzf> bgzfIn;
    SEQAN_ASSERT(open(bgzfIn, filenameBuffer, "rb"));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, bgzfIn, 10);
    SEQAN_ASSERT_EQ(bytesRead, 2);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "12"), 0);
    close(bgzfIn);
}

// A bit more complex usage of writing to FILE *.
SEQAN_DEFINE_TEST(test_stream_bgzf_write_complex_usage)
{
    using namespace seqan;
    
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);
   
    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "wb"));
    for (int i = 0; i < 10; ++i)
        streamWriteChar(f, '0' + i);
    close(f);

    // Read in data and compare.
    Stream<Bgzf> bgzfIn;
    SEQAN_ASSERT(open(bgzfIn, filenameBuffer, "rb"));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, bgzfIn, 10);
    SEQAN_ASSERT_EQ(bytesRead, 10);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
    close(bgzfIn);
}

// Test of streamEof().
SEQAN_DEFINE_TEST(test_stream_bgzf_eof)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    char const * STR = "This is a test.";
    streamWriteBlock(bgzfOut, STR, strlen(STR));
    close(bgzfOut);

    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));
    testStreamEof(f, false);
    close(f);
}

// Test of streamPeek().
SEQAN_DEFINE_TEST(test_stream_bgzf_peek)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    char const * STR = "This is a test.";
    streamWriteBlock(bgzfOut, STR, strlen(STR));
    close(bgzfOut);

    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));
    testStreamPeek(f);
    close(f);
}

// Test of streamReadChar().
SEQAN_DEFINE_TEST(test_stream_bgzf_read_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    char const * STR = "123";
    streamWriteBlock(bgzfOut, STR, strlen(STR));
    close(bgzfOut);

    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));
    testStreamReadChar(f, false);
    close(f);
}

// Test of streamReadBlock().
SEQAN_DEFINE_TEST(test_stream_bgzf_read_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    char const * STR = "XXXXXXXXXX";
    streamWriteBlock(bgzfOut, STR, strlen(STR));
    close(bgzfOut);

    {
        Stream<Bgzf> f;
        SEQAN_ASSERT(open(f, filenameBuffer, "r"));
        testStreamReadBlockHitLimit(f, false);
        close(f);
    }
    {
        Stream<Bgzf> f;
        SEQAN_ASSERT(open(f, filenameBuffer, "r"));
        testStreamReadBlockHitNoLimit(f, false);
        close(f);
    }
}

// Test of streamWriteChar().
SEQAN_DEFINE_TEST(test_stream_bgzf_write_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    testStreamWriteChar(bgzfOut);
    close(bgzfOut);

    // Read in data and compare.
    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, f, 99);
    SEQAN_ASSERT_EQ(bytesRead, 3);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);
    close(f);
}

// Test of streamWrite().
SEQAN_DEFINE_TEST(test_stream_bgzf_write_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    testStreamWriteBlock(bgzfOut);
    close(bgzfOut);

    // Read in data and compare.
    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, f, 99);
    SEQAN_ASSERT_EQ(bytesRead, 8);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);
    close(f);
}

// Test of streamWrite().
SEQAN_DEFINE_TEST(test_stream_bgzf_streamPut)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    testStreamPut(bgzfOut);
    close(bgzfOut);

    // Read in data and compare.
    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, f, 99);
    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(bytesRead, int(sizeof(cmp) - sizeof(char)));
    SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);
    close(f);
}


// Test of streamFlush().
SEQAN_DEFINE_TEST(test_stream_bgzf_flush)
{
    using namespace seqan;

    // Only test that the function is there.
    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    streamFlush(bgzfOut);
    close(bgzfOut);
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_bgzf_seek)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    char const * STR = "0123456789";
    streamWriteBlock(bgzfOut, STR, strlen(STR));
    close(bgzfOut);

    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));

    streamSeek(f, 2, SEEK_SET);
    char c = '\0';
    streamReadChar(c, f);
    SEQAN_ASSERT_EQ(c, '2');
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_bgzf_tell)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    // Write out test data.
    Stream<Bgzf> bgzfOut;
    SEQAN_ASSERT(open(bgzfOut, filenameBuffer, "w"));
    char const * STR = "0123456789";
    streamWriteBlock(bgzfOut, STR, strlen(STR));
    close(bgzfOut);

    Stream<Bgzf> f;
    SEQAN_ASSERT(open(f, filenameBuffer, "r"));

    char c;
    size_t pos = streamTell(f);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamReadChar(c, f);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamReadChar(c, f);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(f);
    SEQAN_ASSERT_EQ(pos, 2u);
}

SEQAN_DEFINE_TEST(test_stream_bgzf_write_large_and_compare_with_file)
{
    using namespace seqan;

    // -----------------------------------------------------------------------
    // Generate Paths, Open Files
    // -----------------------------------------------------------------------
    
    // Open test file for reading.
    const char * r = SEQAN_PATH_TO_ROOT();
    char tempPath[1000];
    strcpy(tempPath, r);
    strcat(tempPath, "/core/tests/stream/SRR067601_1.1k.fasta");
    FILE * fp = fopen(tempPath, "rb");

    // Open BGZF stream for writing.
    const char * p = SEQAN_TEMP_FILENAME();
    char outFilename[1000];
    strcpy(outFilename, p);

    Stream<Bgzf> stream;
    open(stream, outFilename, "w");

    // -----------------------------------------------------------------------
    // Write Data.
    // -----------------------------------------------------------------------

    // Copy from fp to stream.
    String<char> buffer;
    resize(buffer, 765);
    while (!feof(fp))
    {
        int len = fread(&buffer[0], 1, 765, fp);
        streamWriteBlock(stream, &buffer[0], len);
    }

    // Close Stream.
    close(stream);

    // -----------------------------------------------------------------------
    // Compare Data
    // -----------------------------------------------------------------------
    char inPath1[1000];
    strcpy(inPath1, SEQAN_PATH_TO_ROOT());
    strcat(inPath1, "/core/tests/stream/SRR067601_1.1k.fasta.gz");
    FILE * fin1 = fopen(inPath1, "rb");
    SEQAN_ASSERT(fin1 != NULL);
    // printf("inpath:%s\n", inPath1);

    FILE * fin2 = fopen(outFilename, "rb");
    SEQAN_ASSERT(fin2 != NULL);
    // printf("inpath2:%s\n", outFilename);

    int i = 0;
    while (!feof(fin1) && !feof(fin2))
    {
        int i1 = fgetc(fin1);
        int i2 = fgetc(fin2);
        SEQAN_ASSERT_EQ_MSG(i1, i2, "At character pos %d", i);
        ++i;
    }

    SEQAN_ASSERT(feof(fin1));
    SEQAN_ASSERT(feof(fin2));
}

SEQAN_DEFINE_TEST(test_stream_bgzf_from_file_and_compare)
{
    using namespace seqan;

    // Define paths to BGZF and FASTA file.
    char gzPath[1000];
    strcpy(gzPath, SEQAN_PATH_TO_ROOT());
    strcat(gzPath, "/core/tests/stream/SRR067601_1.1k.fasta.gz");
    char fastaPath[1000];
    strcpy(fastaPath, SEQAN_PATH_TO_ROOT());
    strcat(fastaPath, "/core/tests/stream/SRR067601_1.1k.fasta");

    // Open files.
    Stream<Bgzf> inBgzf;
    SEQAN_ASSERT(open(inBgzf, gzPath, "r"));

    FILE * inFasta = fopen(fastaPath, "rb");
    SEQAN_ASSERT(inFasta != NULL);

    // Read files, expect both to be of the same length and have the same
    // (uncompressed) content.
    int i = 0;
    while (!streamEof(inBgzf) && !streamEof(inFasta))
    {
        char c = '\0';
        int res = streamReadChar(c, inBgzf);
        SEQAN_ASSERT_EQ_MSG(res, 0, "Failed at character pos %d", i);
        int i1 = c;
        int i2 = fgetc(inFasta);
        SEQAN_ASSERT_EQ_MSG(i1, i2, "At character pos %d", i);
        ++i;
    }

    SEQAN_ASSERT(streamEof(inBgzf));
    // Note that BGZF files know that they are EOF before reading the last char.  For normal files we have to read
    // beyond the end of the file.
    SEQAN_ASSERT_LT(fgetc(inFasta), 0);
    SEQAN_ASSERT(feof(inFasta));
}

#endif // #ifndef CORE_TESTS_STREAM_TEST_STREAM_BGZF_H_
