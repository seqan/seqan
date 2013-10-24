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
//         David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Tests for the File Stream.
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_FILE_STREAM_H_
#define TEST_STREAM_TEST_STREAM_FILE_STREAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

using namespace seqan;

// ==========================================================================
// Types
// ========================================================================== 

template <typename TSpec_>
class FileStreamTest : public Test
{
public:
    typedef TSpec_ TSpec;
};

// --------------------------------------------------------------------------
// FileStream Specs
// --------------------------------------------------------------------------

typedef
    TagList<Async<>,
    TagList<MMap<>
    > >
    FileStreamSpecs;

SEQAN_TYPED_TEST_CASE(FileStreamTest, FileStreamSpecs);

// --------------------------------------------------------------------------
// FileStream Tests
// --------------------------------------------------------------------------

// Simple example of writing to FILE *.
SEQAN_TYPED_TEST(FileStreamTest, ReadSimpleUsage)
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    char const * STR = "This is a string!\nWith two lines.";
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename)));
    stream.write(STR, strlen(STR));
    SEQAN_ASSERT(stream.good());
    SEQAN_ASSERT(close(stream));

    FileStream<char, Input, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename)));
    testStreamReadSimpleUsage(stream2);
}

// More complex example of writing to FILE *.
SEQAN_TYPED_TEST(FileStreamTest, ReadComplexUsage)
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    char const * STR = "This is a string!\nWith two lines.";
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename)));
    stream.write(STR, strlen(STR));
    SEQAN_ASSERT(stream.good());
    SEQAN_ASSERT(close(stream));

    FileStream<char, Input, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename)));
    testStreamReadComplexUsage(stream2);
}

// Simple example of reading from FILE *.
SEQAN_TYPED_TEST(FileStreamTest, WriteSimpleUsage)
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename)));
    stream.put('1');
    SEQAN_ASSERT(stream.good());
    stream.put('2');
    SEQAN_ASSERT(stream.good());
    SEQAN_ASSERT(close(stream));

    FileStream<char, Input, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename)));
    char buffer[100] = "__________";
    stream2.read(buffer, 10);
    SEQAN_ASSERT(stream2.eof());
    SEQAN_ASSERT_EQ(buffer[0], '1');
    SEQAN_ASSERT_EQ(buffer[1], '2');
    SEQAN_ASSERT(close(stream2));
}

// A bit more complex usage of writing to FILE *.
SEQAN_TYPED_TEST(FileStreamTest, WriteComplexUsage)
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename)));
    for (char c = '0'; c <= '9'; ++c)
        stream.put(c);
    SEQAN_ASSERT(stream.good());
    SEQAN_ASSERT(close(stream));

    FileStream<char, Input, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename)));
    char buffer[100] = "__________";
    stream2.read(buffer, 10);
    SEQAN_ASSERT_NOT(stream2.eof());
    SEQAN_ASSERT(stream.good());
    buffer[10] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
    SEQAN_ASSERT(close(stream2));
}

/*
// Test of streamEof().
template <typename TSpec>
void runTestStreamFileStreamEof()
SEQAN_TYPED_TEST(FileStreamTest, WriteComplexUsage)
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "This is a test.";
    SEQAN_ASSERT_EQ(stream.write(STR, strlen(STR)), strlen(STR));
    close(stream);

    FileStream<char, Output, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamEof(stream2);
}

// Test of streamPeek().
template <typename TSpec>
void runTestStreamFileStreamPeek()
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "This is a test.";
    SEQAN_ASSERT_EQ(stream.write(STR, strlen(STR)), strlen(STR));
    close(stream);

    FileStream<char, Output, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamPeek(stream2);
}

// Test of streamReadChar().
template <typename TSpec>
void runTestStreamFileStreamReadChar()
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "123";
    SEQAN_ASSERT_EQ(stream.write(STR, strlen(STR)), strlen(STR));
    close(stream);

    FileStream<char, Output, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamReadChar(stream2);
}

// Test of streamReadBlock().
template <typename TSpec>
void runTestStreamFileStreamReadBlock()
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "XXXXXXXXXX";
    SEQAN_ASSERT_EQ(stream.write(STR, strlen(STR)), strlen(STR));
    close(stream);

    {
        FileStream<char, Output, typename TestFixture::TSpec> stream2;
        SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
        testStreamReadBlockHitLimit(stream2);
    }
    {
        FileStream<char, Output, typename TestFixture::TSpec> stream2;
        SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
        testStreamReadBlockHitNoLimit(stream2);
    }
}

// Test of streamWriteChar().
template <typename TSpec>
void runTestStreamFileStreamWriteChar()
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamWriteChar(stream);
    close(stream);

    // Read in data and compare.
    FileStream<char, Output, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char buffer[100];
    int bytesRead = stream2.read(buffer, 99);
    SEQAN_ASSERT_EQ(bytesRead, 3);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "345"), 0);
}

// Test of streamWrite().
template <typename TSpec>
void runTestStreamFileStreamWriteBlock()
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamWriteBlock(stream);
    close(stream);

    // Read in data and compare.
    FileStream<char, Output, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char buffer[100];
    int bytesRead = stream2.read(buffer, 99);
    SEQAN_ASSERT_EQ(bytesRead, 8);
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "ABCDEFGH"), 0);
}

// Test of streamWrite().
template <typename TSpec>
void runTestStreamFileStreamStreamPut()
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamPut(stream);
    close(stream);

    // Read in data and compare.
    FileStream<char, Output, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char buffer[100];
    int bytesRead = stream2.read(buffer, 99);
    char cmp[] = "c\nseq\nsss\n12\n34\n56\n78\n5.4\n6.5\nA\nACGT\nACGTN\n";
    buffer[bytesRead] = '\0';
    SEQAN_ASSERT_EQ(bytesRead, int(sizeof(cmp) - sizeof(char)));
    SEQAN_ASSERT_EQ(strcmp(buffer, cmp), 0);
}


// Test of streamFlush().
template <typename TSpec>
void runTestStreamFileStreamFlush()
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    streamFlush(stream);
}

// Test of streamSeek().
template <typename TSpec>
void runTestStreamFileStreamSeek()
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "0123456789";
    SEQAN_ASSERT_EQ(stream.write(STR, strlen(STR)), strlen(STR));
    close(stream);

    // Read in data and compare.
    FileStream<char, Output, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    testStreamSeek(stream2);
}

// Test of streamTell().
template <typename TSpec>
void runTestStreamFileStreamTell()
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "0123456789";
    SEQAN_ASSERT_EQ(stream.write(STR, strlen(STR)), strlen(STR));
    close(stream);

    // Test tell.
    FileStream<char, Output, typename TestFixture::TSpec> stream2;
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
    CharString tempFilename = SEQAN_TEMP_FILENAME();
    unsigned FILE_SIZE = 5 * 1000 * 1000;

    // Write out test data.
    std::fstream f(toCString(tempFilename), std::ios::binary | std::ios::out);
    for (unsigned i = 0; i < FILE_SIZE; ++i)
        f.put('!' + i % 53);  // 53 is prime
    f.flush();
    f.close();

    // Test reading of 5MB file.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
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
    CharString tempFilename = SEQAN_TEMP_FILENAME();
    unsigned FILE_SIZE = 5 * 1000 * 1000;

    // Prepare buffer.
    CharString buffer;
    resize(buffer, FILE_SIZE, '\0');
    for (unsigned i = 0; i < FILE_SIZE; ++i)
        buffer[i] = ('!' + i % 53);  // 53 is prime

    // Write out.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    SEQAN_ASSERT_EQ(stream.write(&buffer[0], (int)FILE_SIZE), (int)FILE_SIZE);
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
    CharString tempFilename = SEQAN_TEMP_FILENAME();
    unsigned FILE_SIZE = 5 * 1000 * 1000;

    // Write out test data.
    std::fstream f(toCString(tempFilename), std::ios::binary | std::ios::out);
    for (unsigned i = 0; i < FILE_SIZE; ++i)
        f.put('!' + i % 53);  // 53 is prime
    f.flush();
    f.close();

    // Test seeking in 5MB file.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
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
*/
#endif  // TEST_STREAM_TEST_STREAM_FILE_STREAM_H_
