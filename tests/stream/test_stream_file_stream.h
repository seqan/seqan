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

#ifdef STDLIB_VS
typedef TagList<Async<> > FileStreamSpecs;
#else
typedef
    TagList<Async<>,
    TagList<MMap<>
    > >
    FileStreamSpecs;
#endif

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
    close(stream);

    File<> file;
    open(file, toCString(tempFilename));
    SEQAN_ASSERT_GT((size_t)length(file), 0u);
    close(file);

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
    close(stream);

    File<> file;
    open(file, toCString(tempFilename));
    SEQAN_ASSERT_GT((size_t)length(file), 0u);
    close(file);

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
    close(stream);

    FileStream<char, Input, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename)));
    char buffer[100] = "__________";
    stream2.read(buffer, 10);
    SEQAN_ASSERT_EQ(buffer[0], '1');
    SEQAN_ASSERT_EQ(buffer[1], '2');
    SEQAN_ASSERT(stream2.eof());
    close(stream2);
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
    close(stream);

    FileStream<char, Input, typename TestFixture::TSpec> stream2;
    SEQAN_ASSERT(open(stream2, toCString(tempFilename)));
    char buffer[100] = "__________";
    stream2.read(buffer, 10);
    SEQAN_ASSERT_NOT(stream2.eof());
    SEQAN_ASSERT(stream.good());
    buffer[10] = '\0';
    SEQAN_ASSERT_EQ(strcmp(buffer, "0123456789"), 0);
    close(stream2);
}

// Test eof().
SEQAN_TYPED_TEST(FileStreamTest, Eof)
{
    char const * STR = "This is a test.";
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    stream.exceptions(std::ios::failbit | std::ios::badbit);
    SEQAN_ASSERT(open(stream, toCString(tempFilename)));
    stream.write(STR, strlen(STR));
    close(stream);

    char buffer[100];
    FileStream<char, Input, typename TestFixture::TSpec> stream2;
    stream2.exceptions(std::ios::badbit);
    SEQAN_ASSERT(open(stream2, toCString(tempFilename)));
    stream2.read(buffer, strlen(STR));
    SEQAN_ASSERT_NOT(stream2.eof());
    int i = stream2.get();
    SEQAN_ASSERT_EQ(i, EOF);
    SEQAN_ASSERT(stream2.eof());
}

//// Test of streamFlush().
//template <typename TSpec>
//void runTestStreamFileStreamFlush()
//{
//    CharString tempFilename = SEQAN_TEMP_FILENAME();
//
//    // Write out test data.
//    FileStream<char, Output, typename TestFixture::TSpec> stream;
//    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
//    streamFlush(stream);
//}

// Test of seek()
SEQAN_TYPED_TEST(FileStreamTest, Seek)
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    FileStream<char, Output, typename TestFixture::TSpec> stream;
    stream.exceptions(std::ios::failbit | std::ios::badbit);
    SEQAN_ASSERT(open(stream, toCString(tempFilename)));
    char const * STR = "0123456789";
    stream.write(STR, strlen(STR));
    close(stream);

    // Read in data and compare.
    FileStream<char, Input, typename TestFixture::TSpec> stream2;
    stream2.exceptions(std::ios::failbit | std::ios::badbit);
    SEQAN_ASSERT(open(stream2, toCString(tempFilename)));
    testStreamSeek(stream2);
}

// Test of streamTell().
SEQAN_TYPED_TEST(FileStreamTest, Tell)
{
    typedef FileStream<char, Output, typename TestFixture::TSpec> TStream;
    typedef FileStream<char, Input, typename TestFixture::TSpec> TStream2;

    CharString tempFilename = SEQAN_TEMP_FILENAME();

    // Write out test data.
    TStream stream;
    stream.exceptions(std::ios::failbit | std::ios::badbit);
    SEQAN_ASSERT(open(stream, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));
    char const * STR = "0123456789";
    stream.write(STR, strlen(STR));;
    close(stream);

    // Test tell.
    TStream2 stream2;
    stream2.exceptions(std::ios::failbit | std::ios::badbit);
    SEQAN_ASSERT(open(stream2, toCString(tempFilename), OPEN_RDWR | OPEN_CREATE | OPEN_APPEND));

    typename DirectionIterator<TStream2, Input>::Type iter = directionIterator(stream2, Input());

    char c = '\0';
    size_t pos = position(iter);
    SEQAN_ASSERT_EQ(pos, 0u);
    c = *iter++;
    c = *iter++;
    ignoreUnusedVariableWarning(c);
    pos = position(iter);
    SEQAN_ASSERT_EQ(pos, 2u);
}

// Test for reading a file of size 5MB, thus more than one page.
SEQAN_TYPED_TEST(FileStreamTest, ReadLarge)
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();
    unsigned FILE_SIZE = 5 * 1000 * 1000;

    // Write out test data.
    std::fstream f(toCString(tempFilename), std::ios::binary | std::ios::out);
    f.exceptions(std::ios::failbit | std::ios::badbit);
    for (unsigned i = 0; i < FILE_SIZE; ++i)
        f.put('!' + i % 53);  // 53 is prime
    f.close();

    // Test reading of 5MB file.
    FileStream<char, Input, typename TestFixture::TSpec> stream;
    stream.exceptions(std::ios::failbit | std::ios::badbit);
    SEQAN_ASSERT(open(stream, toCString(tempFilename)));
    CharString buffer;
    resize(buffer, FILE_SIZE, '\0');
    stream.read(begin(buffer, Standard()), FILE_SIZE);

    for (unsigned i = 0; i < FILE_SIZE; ++i)
        SEQAN_ASSERT_EQ((unsigned)buffer[i], '!' + i % 53);
}

// Test for writing a file of size 5MB, thus more than one page.
SEQAN_TYPED_TEST(FileStreamTest, WriteLarge)
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
    stream.exceptions(std::ios::failbit | std::ios::badbit);
    SEQAN_ASSERT(open(stream, toCString(tempFilename)));
    stream.write(begin(buffer, Standard()), FILE_SIZE);
    close(stream);

    // Read into buffer again and compare.
    clear(buffer);
    resize(buffer, FILE_SIZE, '\0');
    std::fstream f(toCString(tempFilename), std::ios::binary | std::ios::in);
    f.exceptions(std::ios::failbit | std::ios::badbit);

    f.seekg(0, std::ios::end);
    SEQAN_ASSERT_EQ((unsigned)f.tellg(), FILE_SIZE);
    f.seekg(0);

    f.read(begin(buffer, Standard()), FILE_SIZE);
    SEQAN_ASSERT_EQ((unsigned)f.gcount(), FILE_SIZE);

    for (unsigned i = 0; i < FILE_SIZE; ++i)
        SEQAN_ASSERT_EQ((unsigned)buffer[i], '!' + i % 53);
}

// Test for seeking in a file of size 5MB, thus more than one page.
SEQAN_TYPED_TEST(FileStreamTest, SeekLarge)
{
    CharString tempFilename = SEQAN_TEMP_FILENAME();
    unsigned FILE_SIZE = 5 * 1000 * 1000;

    // Write out test data.
    std::fstream f(toCString(tempFilename), std::ios::binary | std::ios::out);
    f.exceptions(std::ios::failbit | std::ios::badbit);
    for (unsigned i = 0; i < FILE_SIZE; ++i)
        f.put('!' + i % 53);  // 53 is prime
    f.flush();
    f.close();

    // Test seeking in 5MB file.
    FileStream<char, Input, typename TestFixture::TSpec> stream;
    stream.exceptions(std::ios::failbit | std::ios::badbit);

    SEQAN_ASSERT(open(stream, toCString(tempFilename)));
    char buffer[13];

    for (unsigned pos = 0; pos + 13 < FILE_SIZE; pos += 1235)
    {
        stream.seekg(pos);
        stream.read(buffer, 13);
        SEQAN_ASSERT_NOT(stream.eof());
        for (unsigned i = 0; i < 13; ++i)
            SEQAN_ASSERT_EQ((unsigned)buffer[i], '!' + (pos + i) % 53);
    }

    stream.exceptions(std::ios::badbit);

    stream.seekg(-5, std::ios::end);
    stream.read(buffer, 13);
    for (unsigned i = 0; i < 5; ++i)
        SEQAN_ASSERT_EQ((unsigned)buffer[i], '!' + (FILE_SIZE - 5 + i) % 53);
    SEQAN_ASSERT(stream.fail());
}

#endif  // TEST_STREAM_TEST_STREAM_FILE_STREAM_H_
