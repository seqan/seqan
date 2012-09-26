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
// Tests for FILE* to the stream concept adaption.
// ==========================================================================

// TODO(holtgrew): Fix comments above tests.

#ifndef TEST_STREAM_TEST_STREAM_ADAPT_SSTREAM_H_
#define TEST_STREAM_TEST_STREAM_ADAPT_SSTREAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

// ==========================================================================
// For std::stringstream
// ==========================================================================

SEQAN_DEFINE_TEST(test_stream_adapt_sstream_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature< ::std::stringstream, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::stringstream, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::stringstream, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::stringstream, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::stringstream, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::stringstream, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::stringstream, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::stringstream, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple example reading with sstream.
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_read_simple_usage)
{
    using namespace seqan;

    std::stringstream sstream;

    char const * STR = "This is a string!\nWith two lines.";
    sstream.write(STR, strlen(STR));
    sstream.seekg(0);
    sstream.seekp(0);
    testStreamReadSimpleUsage(sstream);
}

// A bit more complex example reading with sstream.
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_read_complex_usage)
{
    using namespace seqan;

    std::stringstream sstream;

    char const * STR = "This is a string!\nWith two lines.";
    sstream << STR;
    sstream.seekg(0);
    sstream.seekp(0);
    testStreamReadComplexUsage(sstream);
}

// Simple example of read with sstream.
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_write_simple_usage)
{
    using namespace seqan;

    std::stringstream sstream;

    streamWriteChar(sstream, '1');
    streamWriteChar(sstream, '2');

    sstream.seekg(0);
    char buffer[10];
    sstream.read(buffer, 2);
    sstream.seekp(sstream.tellp());
    buffer[2] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "12"), 0);
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(sstream)), 2);
}

// A bit more complex example of writing with sstream.
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_write_complex_usage)
{
   using namespace seqan;

   std::stringstream sstream;

   for (int i = 0; i < 10; ++i)
       streamWriteChar(sstream, '0' + i);
    
   sstream.seekg(0);
   sstream.seekp(0);
   char buffer[100];
   sstream.read(buffer, 20);
   buffer[10] = '\0';
   SEQAN_ASSERT_EQ(sstream.gcount(), 10);
   SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
   SEQAN_ASSERT_EQ(static_cast<int>(streamTell(sstream)), -1);
}

// Test of streamEof().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_eof)
{
    using namespace seqan;

    std::stringstream sstream;

    sstream.write("This is a test.", 15);
    sstream.seekg(0);
    sstream.seekp(0);
    testStreamEof(sstream);
}

// Test of streamPeek().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_peek)
{
    using namespace seqan;

    std::stringstream sstream;

    sstream.write("This is a test.", 15);
    sstream.seekg(0);
    sstream.seekp(0);
    testStreamPeek(sstream);
}

// Test of streamReadChar().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_read_char)
{
    using namespace seqan;

    std::stringstream sstream;

    sstream.write("123", 3);
    sstream.seekg(0);
    sstream.seekp(0);

    testStreamReadChar(sstream);
}

// Test of streamReadBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_read_block)
{
    using namespace seqan;

    {
        std::stringstream sstream;

        sstream.write("XXXXXXXXXX", 10);
        sstream.seekg(0);
        sstream.seekp(0);

        testStreamReadBlockHitLimit(sstream);
    }
    {
        std::stringstream sstream;

        sstream.write("XXXXXXXXXX", 10);
        sstream.seekg(0);
        sstream.seekp(0);
        
        testStreamReadBlockHitNoLimit(sstream);
    }
}

// Test of streamWriteChar().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_write_char)
{
    using namespace seqan;

    std::stringstream sstream;

    testStreamWriteChar(sstream);

    sstream.seekg(0);
    sstream.seekp(0);

    char buffer[100];
    sstream.read(buffer, 99);
    SEQAN_ASSERT_EQ(sstream.gcount(), 3);
    buffer[3] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);
}

// Test of streamWriteBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_write_block)
{
    using namespace seqan;

    std::stringstream sstream;

    testStreamWriteBlock(sstream);

    sstream.seekg(0);
    sstream.seekp(0);

    char buffer[100];
    sstream.read(buffer, 99);
    SEQAN_ASSERT_EQ(sstream.gcount(), 8);
    buffer[8] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);
}

// Test of streamPut().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_streamPut)
{
    using namespace seqan;

    std::stringstream sstream;

    testStreamPut(sstream);

    sstream.seekg(0);
    sstream.seekp(0);

    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";

    char buffer[100];
    sstream.read(buffer, 99);
    SEQAN_ASSERT_EQ(sstream.gcount(), 43);
    buffer[43] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);
}


// Test of streamFlush().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_flush)
{
    using namespace seqan;
    // Only test that the function is there.

    std::stringstream sstream;

    streamFlush(sstream);
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_seek)
{
    using namespace seqan;
    // Only test that the function is there.

    std::stringstream sstream("0123456789");

    testStreamSeek(sstream);
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_adapt_sstream_tell)
{
    using namespace seqan;

    std::stringstream sstream;

    sstream.write("0123456789", 10);
    sstream.seekg(0);
    sstream.seekp(0);

    char c;
    size_t pos = streamTell(sstream);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamReadChar(c, sstream);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamReadChar(c, sstream);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(sstream);
    SEQAN_ASSERT_EQ(pos, 2u);
}

// ==========================================================================
// For std::isstream
// ==========================================================================

SEQAN_DEFINE_TEST(test_stream_adapt_istringstream_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature< ::std::istringstream, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::istringstream, IsOutput>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::istringstream, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::istringstream, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::istringstream, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::istringstream, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::istringstream, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::istringstream, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple usage of pointers reading with istringstream.
SEQAN_DEFINE_TEST(test_stream_adapt_istringstream_read_simple_usage)
{
    using namespace seqan;

    char const * STR = "This is a string!\nWith two lines.";
    std::istringstream istringstream(STR);

    testStreamReadSimpleUsage(istringstream);
}

// A bit more complex usage writing to istringstream.
SEQAN_DEFINE_TEST(test_stream_adapt_istringstream_read_complex_usage)
{
    using namespace seqan;

    char const * STR = "This is a string!\nWith two lines.";
    std::istringstream istringstream(STR);

    testStreamReadComplexUsage(istringstream);
}

// Test of streamEof().
SEQAN_DEFINE_TEST(test_stream_adapt_istringstream_eof)
{
    using namespace seqan;

    char const * STR = "This is a test.";
    std::istringstream istringstream(STR);

    testStreamEof(istringstream);
}

// Test of streamPeek().
SEQAN_DEFINE_TEST(test_stream_adapt_istringstream_peek)
{
    using namespace seqan;

    char const * STR = "This is a test.";
    std::istringstream istringstream(STR);

    testStreamPeek(istringstream);
}

// Test of streamReadChar().
SEQAN_DEFINE_TEST(test_stream_adapt_istringstream_read_char)
{
    using namespace seqan;

    char const * STR = "123";
    std::istringstream istringstream(STR);

    testStreamReadChar(istringstream);
}

// Test of streamReadBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_istringstream_read_block)
{
    using namespace seqan;

    char const * STR = "XXXXXXXXXX";
    std::istringstream istringstream(STR);

    testStreamReadBlockHitLimit(istringstream);
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_adapt_istringstream_seek)
{
    using namespace seqan;
    // Only test that the function is there.

    char const * STR = "0123456789";
    std::istringstream istringstream(STR);

    testStreamSeek(istringstream);
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_adapt_istringstream_tell)
{
    using namespace seqan;

    char const * STR = "0123456789";
    std::istringstream istringstream(STR);

    char c;
    size_t pos = streamTell(istringstream);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamReadChar(c, istringstream);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamReadChar(c, istringstream);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(istringstream);
    SEQAN_ASSERT_EQ(pos, 2u);
}

// ==========================================================================
// For std::ostringstream
// ==========================================================================

SEQAN_DEFINE_TEST(test_stream_adapt_ostringstream_metafunctions)
{
    using namespace seqan;

    {
        bool b = HasStreamFeature< ::std::ostringstream, IsInput>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ostringstream, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ostringstream, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ostringstream, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ostringstream, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ostringstream, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ostringstream, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature< ::std::ostringstream, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple usage of writing to ostringstream.
SEQAN_DEFINE_TEST(test_stream_adapt_ostringstream_write_simple_usage)
{
    using namespace seqan;

    std::ostringstream ostringstream;
    streamWriteChar(ostringstream, '1');
    streamWriteChar(ostringstream, '2');
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(ostringstream)), 2);

    // Check result.
    SEQAN_ASSERT_EQ(CharString("12"), CharString(ostringstream.str()));
}

// A bit more complex usage of writing to ostringstream.
SEQAN_DEFINE_TEST(test_stream_adapt_ostringstream_write_complex_usage)
{
   using namespace seqan;

   std::ostringstream ostringstream;
   for (int i = 0; i < 10; ++i)
       streamWriteChar(ostringstream, '0' + i);

   // Check result.
   SEQAN_ASSERT_EQ(CharString("0123456789"), CharString(ostringstream.str()));
}

// Test of streamWriteChar().
SEQAN_DEFINE_TEST(test_stream_adapt_ostringstream_write_char)
{
    using namespace seqan;

    std::ostringstream ostringstream;

    testStreamWriteChar(ostringstream);

    // Check result.
    SEQAN_ASSERT_EQ(CharString("345"), CharString(ostringstream.str()));
}

// Test of streamWriteBlock().
SEQAN_DEFINE_TEST(test_stream_adapt_ostringstream_write_block)
{
    using namespace seqan;

    std::ostringstream ostringstream;

    testStreamWriteBlock(ostringstream);

    // Check result.
    SEQAN_ASSERT_EQ(CharString("ABCDEFGH"), CharString(ostringstream.str()));
}

// Test of streamPut().
SEQAN_DEFINE_TEST(test_stream_adapt_ostringstream_streamPut)
{
    using namespace seqan;

    std::ostringstream ostringstream;

    testStreamPut(ostringstream);

    // Check result.
    SEQAN_ASSERT_EQ(CharString("c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n"), CharString(ostringstream.str()));
}

// Test of streamFlush().
SEQAN_DEFINE_TEST(test_stream_adapt_ostringstream_flush)
{
    using namespace seqan;
    // Only test that the function is there.

    std::ostringstream ostringstream;

    streamFlush(ostringstream);
}

// Test of streamSeek().
SEQAN_DEFINE_TEST(test_stream_adapt_ostringstream_seek)
{
    using namespace seqan;
    // Only test that the function is there.

    std::ostringstream ostringstream;

    ostringstream.write("0123456789", 10);
    ostringstream.seekp(1);
    ostringstream.write("0123456789", 10);

    // Check result.
    SEQAN_ASSERT_EQ(CharString("00123456789"), CharString(ostringstream.str()));
}

// Test of streamTell().
SEQAN_DEFINE_TEST(test_stream_adapt_ostringstream_tell)
{
    using namespace seqan;

    std::ostringstream ostringstream;

    ostringstream.write("0123456789", 10);
    ostringstream.seekp(0);

    char c = 'x';
    size_t pos = streamTell(ostringstream);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamWriteChar(ostringstream, c);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamWriteChar(ostringstream, c);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(ostringstream);
    SEQAN_ASSERT_EQ(pos, 2u);
}

#endif  // TEST_STREAM_TEST_STREAM_ADAPT_SSTREAM_H_
