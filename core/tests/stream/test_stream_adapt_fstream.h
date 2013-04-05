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
// Tests for FILE* to the stream concept adaption.
// ==========================================================================

// TODO(holtgrew): Fix comments above tests.

#ifndef TEST_STREAM_TEST_STREAM_ADAPT_FSTREAM_H_
#define TEST_STREAM_TEST_STREAM_ADAPT_FSTREAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

// ==========================================================================
// For std::fstream
// ==========================================================================

SEQAN_DEFINE_TEST(test_stream_adapt_fstream_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature< ::std::fstream, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::fstream, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::fstream, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::fstream, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::fstream, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::fstream, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::fstream, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::fstream, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple example reading with fstream.
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_read_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    char const * STR = "This is a string!\nWith two lines.";
    fstream.write(STR, strlen(STR));
    fstream.seekg(0);
    fstream.seekp(0);
    testStreamReadSimpleUsage(fstream);
}

// A bit more complex example reading with fstream.
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_read_complex_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    char const * STR = "This is a string!\nWith two lines.";
    fstream << STR;
    fstream.seekg(0);
    fstream.seekp(0);
    testStreamReadComplexUsage(fstream);
}

// Simple example of read with fstream.
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_write_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    streamWriteChar(fstream, '1');
    streamWriteChar(fstream, '2');

    fstream.seekg(0);
    fstream.seekp(0);
    char buffer[10];
    fstream.read(buffer, 2);
    buffer[2] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "12"), 0);
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(fstream)), 2);
}

// A bit more complex example of writing with fstream.
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_write_complex_usage)
{
   using namespace seqan;

   const char * tempFilename = SEQAN_TEMP_FILENAME();
   char filenameBuffer[1000];
   strcpy(filenameBuffer, tempFilename);
   
   std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
   SEQAN_ASSERT(fstream.is_open());

   for (int i = 0; i < 10; ++i)
       streamWriteChar(fstream, '0' + i);
    
   fstream.seekg(0);
   fstream.seekp(0);
   char buffer[100];
   fstream.read(buffer, 20);
   buffer[10] = '\0';
   SEQAN_ASSERT_EQ(fstream.gcount(), 10);
   SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
   SEQAN_ASSERT_EQ(static_cast<int>(streamTell(fstream)), -1);
}

// Test of streamEof().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_eof)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);
    
    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());
    fstream.write("This is a test.", 15);
    fstream.seekg(0);
    fstream.seekp(0);
    testStreamEof(fstream);
}

// Test of streamPeek().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_peek)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());
    fstream.write("This is a test.", 15);
    fstream.seekg(0);
    fstream.seekp(0);
    testStreamPeek(fstream);
}

// Test of streamReadChar().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_read_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());
    fstream.write("123", 3);
    fstream.seekg(0);
    fstream.seekp(0);

    testStreamReadChar(fstream);
}

// Test of streamReadBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_read_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);
    
    {
        std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        fstream.write("XXXXXXXXXX", 10);
        fstream.seekg(0);
        fstream.seekp(0);

        testStreamReadBlockHitLimit(fstream);
    }
    {
        std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        fstream.write("XXXXXXXXXX", 10);
        fstream.seekg(0);
        fstream.seekp(0);
        
        testStreamReadBlockHitNoLimit(fstream);
    }
}

// Test of streamWriteChar().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_write_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());
    testStreamWriteChar(fstream);

    fstream.seekg(0);
    fstream.seekp(0);

    char buffer[100];
    fstream.read(buffer, 99);
    SEQAN_ASSERT_EQ(fstream.gcount(), 3);
    buffer[3] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);
}

// Test of streamWriteBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_write_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    testStreamWriteBlock(fstream);

    fstream.seekg(0);
    fstream.seekp(0);

    char buffer[100];
    fstream.read(buffer, 99);
    SEQAN_ASSERT_EQ(fstream.gcount(), 8);
    buffer[8] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);
}

// Test of streamPut().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_streamPut)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());

    testStreamPut(fstream);

    fstream.seekg(0);
    fstream.seekp(0);

    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";

    char buffer[100];
    fstream.read(buffer, 99);
    SEQAN_ASSERT_EQ(fstream.gcount(), 43);
    buffer[43] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);
    fstream.close();
}


// Test of streamFlush().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_flush)
{
    using namespace seqan;
    // Only test that the function is there.

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());
    streamFlush(fstream);
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_seek)
{
    using namespace seqan;
    // Only test that the function is there.

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());
    fstream.write("0123456789", 10);
    fstream.seekg(0);
    fstream.seekp(0);

    testStreamSeek(fstream);
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_adapt_fstream_tell)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(fstream.is_open());
    fstream.write("0123456789", 10);
    fstream.seekg(0);
    fstream.seekp(0);

    char c;
    size_t pos = streamTell(fstream);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamReadChar(c, fstream);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamReadChar(c, fstream);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(fstream);
    SEQAN_ASSERT_EQ(pos, 2u);
}

// ==========================================================================
// For std::ifstream
// ==========================================================================

SEQAN_DEFINE_TEST(test_stream_adapt_ifstream_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature< ::std::ifstream, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ifstream, IsOutput>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ifstream, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ifstream, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ifstream, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ifstream, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ifstream, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ifstream, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple usage of pointers reading with ifstream.
SEQAN_DEFINE_TEST(test_stream_adapt_ifstream_read_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    { // Write data.
        std::fstream fstream(filenameBuffer, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        char const * STR = "This is a string!\nWith two lines.";
        fstream.write(STR, strlen(STR));
    }

    std::ifstream ifstream(filenameBuffer);
    SEQAN_ASSERT(ifstream.is_open());
    testStreamReadSimpleUsage(ifstream);
}

// A bit more complex usage writing to ifstream.
SEQAN_DEFINE_TEST(test_stream_adapt_ifstream_read_complex_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    { // Write data.
        std::fstream fstream(filenameBuffer, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        char const * STR = "This is a string!\nWith two lines.";
        fstream.write(STR, strlen(STR));
    }

    std::ifstream ifstream(filenameBuffer);
    SEQAN_ASSERT(ifstream.is_open());
    testStreamReadComplexUsage(ifstream);
}

// Test of streamEof().
SEQAN_DEFINE_TEST(test_stream_adapt_ifstream_eof)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    { // Write data.
        std::fstream fstream(filenameBuffer, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        fstream.write("This is a test.", 15);
    }
   
    std::ifstream ifstream(filenameBuffer);
    SEQAN_ASSERT(ifstream.is_open());
    testStreamEof(ifstream);
}

// Test of streamPeek().
SEQAN_DEFINE_TEST(test_stream_adapt_ifstream_peek)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    { // Write data.
        std::fstream fstream(filenameBuffer, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        fstream.write("This is a test.", 15);
    }

    std::ifstream ifstream(filenameBuffer);
    testStreamPeek(ifstream);
}

// Test of streamReadChar().
SEQAN_DEFINE_TEST(test_stream_adapt_ifstream_read_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    { // Write data.
        std::fstream fstream(filenameBuffer, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        fstream.write("123", 3);
    }

    std::ifstream ifstream(filenameBuffer);
    testStreamReadChar(ifstream);
}

// Test of streamReadBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_ifstream_read_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);
    

    { // Write data.
        std::fstream fstream(filenameBuffer, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        fstream.write("XXXXXXXXXX", 10);
    }

    {
        std::ifstream ifstream(filenameBuffer);
        SEQAN_ASSERT(ifstream.is_open());
        testStreamReadBlockHitLimit(ifstream);
    }
    {
        std::ifstream ifstream(filenameBuffer);
        SEQAN_ASSERT(ifstream.is_open());
        testStreamReadBlockHitNoLimit(ifstream);
    }
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_adapt_ifstream_seek)
{
    using namespace seqan;
    // Only test that the function is there.

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    { // Write data.
        std::fstream fstream(filenameBuffer, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        fstream.write("0123456789", 10);
    }

    std::ifstream ifstream(filenameBuffer);
    SEQAN_ASSERT(ifstream.is_open());
    testStreamSeek(ifstream);
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_adapt_ifstream_tell)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    { // Write data.
        std::fstream fstream(filenameBuffer, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        SEQAN_ASSERT(fstream.is_open());
        fstream.write("0123456789", 10);
    }

    std::ifstream ifstream(filenameBuffer);
    SEQAN_ASSERT(ifstream.is_open());

    char c;
    size_t pos = streamTell(ifstream);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamReadChar(c, ifstream);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamReadChar(c, ifstream);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(ifstream);
    SEQAN_ASSERT_EQ(pos, 2u);
}

// ==========================================================================
// For std::ofstream
// ==========================================================================

SEQAN_DEFINE_TEST(test_stream_adapt_ofstream_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature< ::std::ofstream, IsInput>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ofstream, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ofstream, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ofstream, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ofstream, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ofstream, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ofstream, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ofstream, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple usage of writing to ofstream.
SEQAN_DEFINE_TEST(test_stream_adapt_ofstream_write_simple_usage)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::ofstream ofstream(filenameBuffer);
    SEQAN_ASSERT(ofstream.is_open());
    streamWriteChar(ofstream, '1');
    streamWriteChar(ofstream, '2');
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(ofstream)), 2);
    ofstream.close();

    // Check result.
    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(fstream.is_open());
    char buffer[10];
    fstream.read(buffer, 2);
    buffer[2] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "12"), 0);
}

// A bit more complex usage of writing to ofstream.
SEQAN_DEFINE_TEST(test_stream_adapt_ofstream_write_complex_usage)
{
   using namespace seqan;

   const char * tempFilename = SEQAN_TEMP_FILENAME();
   char filenameBuffer[1000];
   strcpy(filenameBuffer, tempFilename);
   
   std::ofstream ofstream(filenameBuffer);
   SEQAN_ASSERT(ofstream.is_open());
   for (int i = 0; i < 10; ++i)
       streamWriteChar(ofstream, '0' + i);
   ofstream.close();

   // Check result.
   std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::binary);
   SEQAN_ASSERT(fstream.is_open());
   char buffer[100];
   fstream.read(buffer, 20);
   buffer[10] = '\0';
   SEQAN_ASSERT_EQ(fstream.gcount(), 10);
   SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
   SEQAN_ASSERT_EQ(static_cast<int>(streamTell(fstream)), -1);
}

// Test of streamWriteChar().
SEQAN_DEFINE_TEST(test_stream_adapt_ofstream_write_char)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::ofstream ofstream(filenameBuffer);
    SEQAN_ASSERT(ofstream.is_open());
    testStreamWriteChar(ofstream);
    ofstream.close();

   // Check result.
   std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::binary);
   SEQAN_ASSERT(fstream.is_open());
   char buffer[100];
   fstream.read(buffer, 99);
   SEQAN_ASSERT_EQ(fstream.gcount(), 3);
   buffer[3] = '\0';
   SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);
}

// Test of streamWriteBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_ofstream_write_block)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::ofstream ofstream(filenameBuffer);
    SEQAN_ASSERT(ofstream.is_open());
    testStreamWriteBlock(ofstream);
    ofstream.close();

    // Check result.
    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(fstream.is_open());
    char buffer[100];
    fstream.read(buffer, 99);
    SEQAN_ASSERT_EQ(fstream.gcount(), 8);
    buffer[8] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);
}

// Test of streamPut().
SEQAN_DEFINE_TEST(test_stream_adapt_ofstream_streamPut)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::ofstream ofstream(filenameBuffer, std::ios_base::out | std::ios_base::binary);
    SEQAN_ASSERT(ofstream.is_open());

    testStreamPut(ofstream);
    ofstream.close();

    // Check result.
    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(fstream.is_open());
    fstream.seekg(0);
    fstream.seekp(0);

    char buffer[1000];
    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";

    fstream.read(buffer, 99);
    SEQAN_ASSERT_EQ(fstream.gcount(), 43);
    buffer[43] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);
    fstream.close();
}

// Test of streamFlush().
SEQAN_DEFINE_TEST(test_stream_adapt_ofstream_flush)
{
    using namespace seqan;
    // Only test that the function is there.

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::ofstream ofstream(filenameBuffer);
    SEQAN_ASSERT(ofstream.is_open());
    streamFlush(ofstream);
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_adapt_ofstream_seek)
{
    using namespace seqan;
    // Only test that the function is there.

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::ofstream ofstream(filenameBuffer);
    SEQAN_ASSERT(ofstream.is_open());
    ofstream.write("0123456789", 10);
    ofstream.seekp(1);
    ofstream.write("0123456789", 10);
    ofstream.close();

    // Check result.
    std::fstream fstream(filenameBuffer, std::ios_base::in | std::ios_base::binary);
    SEQAN_ASSERT(fstream.is_open());
    char buffer[100];
    fstream.read(buffer, 99);
    SEQAN_ASSERT_EQ(fstream.gcount(), 11);
    buffer[11] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "00123456789"), 0);
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_adapt_ofstream_tell)
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::ofstream ofstream(filenameBuffer);
    SEQAN_ASSERT(ofstream.is_open());
    ofstream.write("0123456789", 10);
    ofstream.seekp(0);

    char c;
    size_t pos = streamTell(ofstream);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamWriteChar(ofstream, c);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamWriteChar(ofstream, c);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(ofstream);
    SEQAN_ASSERT_EQ(pos, 2u);
}

#endif  // TEST_STREAM_TEST_STREAM_ADAPT_FSTREAM_H_
