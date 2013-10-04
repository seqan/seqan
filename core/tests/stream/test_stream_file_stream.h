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

#ifndef TEST_STREAM_TEST_STREAM_FILE_STREAM_H_
#define TEST_STREAM_TEST_STREAM_FILE_STREAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

template <typename TSpec>
void runTestStreamFileStreamMetafunctions()
{
    using namespace seqan;

    {
        bool b = HasStreamFeature<Stream<TSpec>, IsInput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<TSpec>, IsOutput>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<TSpec>, HasPeek>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<TSpec>, HasFilename>::Type::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    {
        bool b = HasStreamFeature<Stream<TSpec>, Seek<OriginBegin> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<TSpec>, Seek<OriginEnd> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<TSpec>, Seek<OriginCurrent> >::Type::VALUE;
        SEQAN_ASSERT(b);
    }
    {
        bool b = HasStreamFeature<Stream<TSpec>, Tell>::Type::VALUE;
        SEQAN_ASSERT(b);
    }
}

// Simple example of writing to FILE *.
template <typename TSpec>
void runTestStreamFileStreamReadSimpleUsage()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    char const * STR = "This is a string!\nWith two lines.";
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    int len = strlen(STR);
    SEQAN_ASSERT_EQ(streamWriteBlock(stream, STR, len), len);
    close(stream);

    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamReadSimpleUsage(stream2);
}

// More complex example of writing to FILE *.
template <typename TSpec>
void runTestStreamFileStreamReadComplexUsage()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    char const * STR = "This is a string!\nWith two lines.";
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    int len = strlen(STR);
    SEQAN_ASSERT_EQ(streamWriteBlock(stream, STR, len), len);
    close(stream);

    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamReadComplexUsage(stream2);
}

// Simple example of reading from FILE *.
template <typename TSpec>
void runTestStreamFileStreamWriteSimpleUsage()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    SEQAN_ASSERT_EQ(streamWriteChar(stream, '1'), 0);
    SEQAN_ASSERT_EQ(streamWriteChar(stream, '2'), 0);
    close(stream);

    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, stream2, 10);
    SEQAN_ASSERT_EQ(bytesRead, 2);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "12"), 0);
    close(stream2);
}

// A bit more complex usage of writing to FILE *.
template <typename TSpec>
void runTestStreamFileStreamWriteComplexUsage()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    for (int i = 0; i < 10; ++i)
        SEQAN_ASSERT_EQ(streamWriteChar(stream, (char)('0' + i)), 0);
    close(stream);

    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, stream2, 10);
    SEQAN_ASSERT_EQ(bytesRead, 10);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
    close(stream2);
}

// Test of streamEof().
template <typename TSpec>
void runTestStreamFileStreamEof()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "This is a test.";
    SEQAN_ASSERT_EQ(streamWriteBlock(stream, STR, strlen(STR)), strlen(STR));
    close(stream);

    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamEof(stream2);
}

// Test of streamPeek().
template <typename TSpec>
void runTestStreamFileStreamPeek()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "This is a test.";
    SEQAN_ASSERT_EQ(streamWriteBlock(stream, STR, strlen(STR)), strlen(STR));
    close(stream);

    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamPeek(stream2);
}

// Test of streamReadChar().
template <typename TSpec>
void runTestStreamFileStreamReadChar()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "123";
    SEQAN_ASSERT_EQ(streamWriteBlock(stream, STR, strlen(STR)), strlen(STR));
    close(stream);

    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamReadChar(stream2);
}

// Test of streamReadBlock().
template <typename TSpec>
void runTestStreamFileStreamReadBlock()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "XXXXXXXXXX";
    SEQAN_ASSERT_EQ(streamWriteBlock(stream, STR, strlen(STR)), strlen(STR));
    close(stream);

    {
        Stream<TSpec> stream2;
        SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
        testStreamReadBlockHitLimit(stream2);
    }
    {
        Stream<TSpec> stream2;
        SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
        testStreamReadBlockHitNoLimit(stream2);
    }
}

// Test of streamWriteChar().
template <typename TSpec>
void runTestStreamFileStreamWriteChar()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamWriteChar(stream);
    close(stream);

    // Read in data and compare.
    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, stream2, 99);
    SEQAN_ASSERT_EQ(bytesRead, 3);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);
}

// Test of streamWrite().
template <typename TSpec>
void runTestStreamFileStreamWriteBlock()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamWriteBlock(stream);
    close(stream);

    // Read in data and compare.
    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, stream2, 99);
    SEQAN_ASSERT_EQ(bytesRead, 8);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);
}

// Test of streamWrite().
template <typename TSpec>
void runTestStreamFileStreamStreamPut()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamPut(stream);
    close(stream);

    // Read in data and compare.
    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char buffer[100];
    int bytesRead = streamReadBlock(buffer, stream2, 99);
    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(bytesRead, int(sizeof(cmp) - sizeof(char)));
    SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);
}


// Test of streamFlush().
template <typename TSpec>
void runTestStreamFileStreamFlush()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    streamFlush(stream);
}

// Test of streamSeek().
template <typename TSpec>
void runTestStreamFileStreamSeek()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "0123456789";
    SEQAN_ASSERT_EQ(streamWriteBlock(stream, STR, strlen(STR)), strlen(STR));
    close(stream);

    // Read in data and compare.
    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamSeek(stream2);
}

// Test of streamTell().
template <typename TSpec>
void runTestStreamFileStreamTell()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "0123456789";
    SEQAN_ASSERT_EQ(streamWriteBlock(stream, STR, strlen(STR)), strlen(STR));
    close(stream);

    // Test tell.
    Stream<TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));

    char c = '\0';
    size_t pos = streamTell(stream2);
    SEQAN_ASSERT_EQ(pos, 0u);
    int res = streamReadChar(c, stream2);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamReadChar(c, stream2);
    SEQAN_ASSERT_EQ(res, 0);
    pos = streamTell(stream2);
    SEQAN_ASSERT_EQ(pos, 2u);
}

// Test for reading a file of size 5MB, thus more than one page.
template <typename TSpec>
void runTestStreamFileStreamReadLarge()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();
    unsigned FILE_SIZE = 5 * 1000 * 1000;

    // Write out test data.
    std::fstream f(toCString(tempFilename), std::ios::binary | std::ios::out);
    for (unsigned i = 0; i < FILE_SIZE; ++i)
        f.put('!' + i % 53);  // 53 is prime
    f.flush();
    f.close();

    // Test reading of 5MB file.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    CharString buffer;
    resize(buffer, FILE_SIZE, '\0');
    SEQAN_ASSERT_EQ(streamReadBlock(&buffer[0], stream, FILE_SIZE), FILE_SIZE);

    for (unsigned i = 0; i < FILE_SIZE; ++i)
        SEQAN_ASSERT_EQ((unsigned)buffer[i], '!' + i % 53);
}

// Test for writing a file of size 5MB, thus more than one page.
template <typename TSpec>
void runTestStreamFileStreamWriteLarge()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();
    unsigned FILE_SIZE = 5 * 1000 * 1000;

    // Prepare buffer.
    CharString buffer;
    resize(buffer, FILE_SIZE, '\0');
    for (unsigned i = 0; i < FILE_SIZE; ++i)
        buffer[i] = ('!' + i % 53);  // 53 is prime

    // Write out.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    SEQAN_ASSERT_EQ(streamWriteBlock(stream, &buffer[0], (int)FILE_SIZE), (int)FILE_SIZE);
    close(stream);

    // Read into buffer again and compare.
    clear(buffer);
    resize(buffer, FILE_SIZE, '\0');
    std::fstream f(toCString(tempFilename), std::ios::binary | std::ios::in);
    f.read(&buffer[0], FILE_SIZE);
    SEQAN_ASSERT_EQ((unsigned)f.gcount(), FILE_SIZE);

    for (unsigned i = 0; i < FILE_SIZE; ++i)
        SEQAN_ASSERT_EQ((unsigned)buffer[i], '!' + i % 53);
}

// Test for seeking in a file of size 5MB, thus more than one page.
template <typename TSpec>
void runTestStreamFileStreamSeekLarge()
{
    using namespace seqan;

    CharString tempFilename = SEQAN_TEMP_FILENAME();
    unsigned FILE_SIZE = 5 * 1000 * 1000;

    // Write out test data.
    std::fstream f(toCString(tempFilename), std::ios::binary | std::ios::out);
    for (unsigned i = 0; i < FILE_SIZE; ++i)
        f.put('!' + i % 53);  // 53 is prime
    f.flush();
    f.close();

    // Test seeking in 5MB file.
    Stream<TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    CharString buffer;
    resize(buffer, 13, '\0');

    for (unsigned pos = 0; pos + 13 < FILE_SIZE; pos += 1235)
    {
        SEQAN_ASSERT_EQ(streamSeek(stream, pos, SEEK_SET), 0);
        SEQAN_ASSERT_EQ(streamReadBlock(&buffer[0], stream, 13), 13);
        for (unsigned i = 0; i < 13; ++i)
            SEQAN_ASSERT_EQ((unsigned)buffer[i], '!' + (pos + i) % 53);
    }

    SEQAN_ASSERT_EQ(streamSeek(stream, -5, SEEK_END), 0);
    SEQAN_ASSERT_EQ(streamReadBlock(&buffer[0], stream, 13), 5);
    for (unsigned i = 0; i < 5; ++i)
        SEQAN_ASSERT_EQ((unsigned)buffer[i], '!' + (FILE_SIZE - 5 + i) % 53);

}

// ===========================================================================
// Calls to the generic test functions
// ===========================================================================

SEQAN_DEFINE_TEST(test_stream_file_stream_metafunctions_file)
{
    runTestStreamFileStreamMetafunctions<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_metafunctions_mmap)
{
    runTestStreamFileStreamMetafunctions<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_simple_usage_file)
{
    runTestStreamFileStreamReadSimpleUsage<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_simple_usage_mmap)
{
    runTestStreamFileStreamReadSimpleUsage<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_complex_usage_file)
{
    runTestStreamFileStreamReadComplexUsage<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_complex_usage_mmap)
{
    runTestStreamFileStreamReadComplexUsage<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_simple_usage_file)
{
    runTestStreamFileStreamWriteSimpleUsage<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_simple_usage_mmap)
{
    runTestStreamFileStreamWriteSimpleUsage<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_complex_usage_file)
{
    runTestStreamFileStreamWriteComplexUsage<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_complex_usage_mmap)
{
    runTestStreamFileStreamWriteComplexUsage<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_eof_file)
{
    runTestStreamFileStreamEof<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_eof_mmap)
{
    runTestStreamFileStreamEof<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_peek_file)
{
    runTestStreamFileStreamPeek<seqan::FileStream<> >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_peek_mmap)
{
    runTestStreamFileStreamPeek<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_char_file)
{
    runTestStreamFileStreamReadChar<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_char_mmap)
{
    runTestStreamFileStreamReadChar<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_block_file)
{
    runTestStreamFileStreamReadBlock<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_block_mmap)
{
    runTestStreamFileStreamReadBlock<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_char_file)
{
    runTestStreamFileStreamWriteChar<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_char_mmap)
{
    runTestStreamFileStreamWriteChar<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_block_file)
{
    runTestStreamFileStreamWriteBlock<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_block_mmap)
{
    runTestStreamFileStreamWriteBlock<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_streamPut_file)
{
    runTestStreamFileStreamStreamPut<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_streamPut_mmap)
{
    runTestStreamFileStreamStreamPut<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_flush_file)
{
    runTestStreamFileStreamFlush<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_flush_mmap)
{
    runTestStreamFileStreamFlush<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_seek_file)
{
    runTestStreamFileStreamSeek<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_seek_mmap)
{
    runTestStreamFileStreamSeek<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_tell_file)
{
    runTestStreamFileStreamTell<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_tell_mmap)
{
    runTestStreamFileStreamTell<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_large_file)
{
    runTestStreamFileStreamReadLarge<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_read_large_mmap)
{
    runTestStreamFileStreamReadLarge<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_large_file)
{
    runTestStreamFileStreamWriteLarge<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_write_large_mmap)
{
    runTestStreamFileStreamWriteLarge<seqan::FileStream<seqan::MMap<> > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_seek_large_file)
{
    runTestStreamFileStreamSeekLarge<seqan::FileStream<seqan::File<seqan::Async<> > > >();
}

SEQAN_DEFINE_TEST(test_stream_file_stream_seek_large_mmap)
{
    runTestStreamFileStreamSeekLarge<seqan::FileStream<seqan::MMap<> > >();
}

#endif  // TEST_STREAM_TEST_STREAM_FILE_STREAM_H_
