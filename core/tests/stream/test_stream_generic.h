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
// Generic test code for stream module.
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_GENERIC_H_
#define TEST_STREAM_TEST_STREAM_GENERIC_H_

// Simple usage of the stream.
template <typename TStream>
void testStreamReadSimpleUsage(TStream & stream)
{
    using namespace seqan;

    // Read first line.
    char buffer[100];
    char *ptr = buffer;
    while (!streamEof(stream)) {
        int res = streamReadChar(*ptr, stream);
        SEQAN_ASSERT_EQ(res, 0);
        SEQAN_ASSERT_LT(ptr, buffer + 100);
        if (*(ptr++) == '\n')
            break;
    }
    char c;
    streamReadChar(c, stream);
    SEQAN_ASSERT_EQ(c, 'W');
}

// A bit more complex usage of the stream.
template <typename TStream>
void testStreamReadComplexUsage(TStream & stream)
{
    using namespace seqan;

    // Read lines and check whether the last one is terminated with a '\n'.
    char c = '\0';
    bool wasEol = false;
    while (!streamEof(stream)) {
        int res = streamReadChar(c, stream);
        if (res == EOF)
            break;
        wasEol = (c == '\n');
    }
    SEQAN_ASSERT_NOT(wasEol);
}

// Helper for testStreamEof() below, to conditionally run streamTell, or
// return defaultValue.
template <typename TStream>
int _runStreamTell(TStream & stream, int /*defaultValue*/, seqan::True const &)
{
    using namespace seqan;
    return streamTell(stream);
}

template <typename TStream>
int _runStreamTell(TStream & /*stream*/, int defaultValue, seqan::False const &)
{
    using namespace seqan;
    return defaultValue;
}

// Test of streamEof().
template <typename TStream>
void testStreamEof(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    const char buffer[] = "This is a test.";
    
    SEQAN_ASSERT_NOT(streamEof(stream));
    char c;
    for (int i = 0; i < 15; ++i) {
        SEQAN_ASSERT_NOT(streamEof(stream));
        int res = streamReadChar(c, stream);
        SEQAN_ASSERT_EQ(res, 0);
        SEQAN_ASSERT_EQ(c, buffer[i]);
        if (checkTell)
            SEQAN_ASSERT_EQ(_runStreamTell(stream, i + 1, typename HasStreamFeature<TStream, Tell>::Type()), i + 1);
    }
    // SEQAN_ASSERT_NOT(streamEof(stream));  // Not portable, inconsistent behaviour.
    int res = streamReadChar(c, stream);
    SEQAN_ASSERT_NEQ(res, 0);
    SEQAN_ASSERT(streamEof(stream));
    if (checkTell)
        SEQAN_ASSERT(_runStreamTell(stream, 15, typename HasStreamFeature<TStream, Tell>::Type()) == 15 ||
                     _runStreamTell(stream, -1, typename HasStreamFeature<TStream, Tell>::Type()) == -1);
}

// Test of streamPeek().
template <typename TStream>
void testStreamPeek(TStream & stream)
{
    using namespace seqan;

    char c;
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(stream)), 0);
    int res = streamPeek(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(c, 'T');
    res = streamPeek(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(c, 'T');
    res = streamPeek(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(c, 'T');
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(stream)), 0);
}

// Test of testStreamReadBlock(), limit is longer than stream.
template <typename TStream>
void testStreamReadBlockHitLimit(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    String<char> buffer;
    resize(buffer, 20);
    size_t charsRead = streamReadBlock(begin(buffer, Standard()), stream, 20);
    resize(buffer, charsRead);
    SEQAN_ASSERT_EQ(charsRead, 10u);
    SEQAN_ASSERT_EQ(strcmp(toCString(buffer), "XXXXXXXXXX"), 0);
    if (checkTell)
        SEQAN_ASSERT(_runStreamTell(stream, 10, typename HasStreamFeature<TStream, Tell>::Type()) == 10 ||
                     _runStreamTell(stream, -1, typename HasStreamFeature<TStream, Tell>::Type()) == -1);
}

// Test of testStreamReadBlock(), limit is shorter than stream.
template <typename TStream>
void testStreamReadBlockHitNoLimit(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    String<char> buffer;
    resize(buffer, 20);
    size_t charsRead = streamReadBlock(begin(buffer, Standard()), stream, 5);
    resize(buffer, charsRead);
    SEQAN_ASSERT_EQ(charsRead, 5u);
    SEQAN_ASSERT(buffer == "XXXXX");
    if (checkTell)
        SEQAN_ASSERT(_runStreamTell(stream, 5, typename HasStreamFeature<TStream, Tell>::Type()) == 5);
}

// Test of testStreamWriteBlock(), iterator interface.
template <typename TStream>
void testStreamWriteBlock(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    String<char> buffer;
    append(buffer, "ABCDEFGH");
    size_t charsWritten = streamWriteBlock(stream, begin(buffer, Standard()), length(buffer));
    SEQAN_ASSERT_EQ(charsWritten, 8u);
    if (checkTell)
        SEQAN_ASSERT_EQ(_runStreamTell(stream, 8, typename HasStreamFeature<TStream, Tell>::Type()), 8);
}

// Test of streamWriteChar().
template <typename TStream>
void testStreamWriteChar(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    int res = streamWriteChar(stream, '3');
    SEQAN_ASSERT_EQ(res, 0);
    res = streamWriteChar(stream, '4');
    SEQAN_ASSERT_EQ(res, 0);
    res = streamWriteChar(stream, '5');
    SEQAN_ASSERT_EQ(res, 0);
    if (checkTell)
        SEQAN_ASSERT_EQ(_runStreamTell(stream, 3, typename HasStreamFeature<TStream, Tell>::Type()), 3);
}

template <typename TStream>
void testStreamPut(TStream & stream)
{
    using namespace seqan;

    int res = streamPut(stream, 'c');
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, CharString("seq"));
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, "sss");
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, 12);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, 34u);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, 56l);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, 78ul);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, 5.4f);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, 6.5);
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, Dna('A'));
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, DnaString("ACGT"));
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);

    res = streamPut(stream, Dna5String("ACGTN"));
    SEQAN_ASSERT_EQ(res, 0);
    res = streamPut(stream, '\n');
    SEQAN_ASSERT_EQ(res, 0);
}


// Test of streamReadChar().
template <typename TStream>
void testStreamReadChar(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    char c;
    int res = streamReadChar(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(c, '1');
    res = streamReadChar(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(c, '2');
    res = streamReadChar(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(c, '3');
    // SEQAN_ASSERT_NOT(streamEof(stream));  // Non-portable, different behaviour.
    res = streamReadChar(c, stream);
    SEQAN_ASSERT_NEQ(res, 0);
    SEQAN_ASSERT(streamEof(stream));
    if (checkTell)
        SEQAN_ASSERT(_runStreamTell(stream, 3, typename HasStreamFeature<TStream, Tell>::Type()) == 3 ||
                     _runStreamTell(stream, -1, typename HasStreamFeature<TStream, Tell>::Type()) == -1);
}

// Test of streamSeek().
template <typename TStream>
void testStreamSeek(TStream & stream)
{
    using namespace seqan;

    // TODO(holtgrew): Do real tests for all variants.
    // TODO(holtgrew): Use constant for origin in the future.
    char c;
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(stream)), 0);
    int res = streamPeek(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(c, '0');
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(stream)), 0);

    res = streamSeek(stream, 4, SEEK_CUR);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(stream)), 4);

    res = streamPeek(c, stream);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(c, '4');
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(stream)), 4);

    res = streamSeek(stream, -2, SEEK_CUR);
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(stream)), 2);

    res = streamPeek(c, stream);
    SEQAN_ASSERT_EQ(c, '2');
    SEQAN_ASSERT_EQ(static_cast<int>(streamTell(stream)), 2);
}

#endif  // TEST_STREAM_TEST_STREAM_GENERIC_H_
