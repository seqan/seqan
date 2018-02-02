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

    typename DirectionIterator<TStream, Input>::Type iter = directionIterator(stream, Input());

    // Read first line.
    char buffer[100];
    char *ptr = buffer;
    while (!atEnd(iter))
    {
        *ptr = *iter++;
        SEQAN_ASSERT(stream.good());
        SEQAN_ASSERT_LT(ptr, buffer + 100);
        if (*(ptr++) == '\n')
            break;
    }
    char c = *iter++;
    SEQAN_ASSERT_EQ(c, 'W');
}

// A bit more complex usage of the stream.
template <typename TStream>
void testStreamReadComplexUsage(TStream & stream)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Input>::Type iter = directionIterator(stream, Input());

    // Read lines and check whether the last one is terminated with a '\n'.
    bool wasEol = false;
    while (!atEnd(iter))
    {
        char c = *iter++;
        wasEol = (c == '\n');
    }
    SEQAN_ASSERT_NOT(wasEol);
}

// Test of streamEof().
template <typename TStream>
void testStreamEof(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Output>::Type iter = directionIterator(stream, Output());

    const char buffer[] = "This is a test.";

    SEQAN_ASSERT_NOT(atEnd(iter));
    char c;
    for (int i = 0; i < 15; ++i)
    {
        SEQAN_ASSERT_NOT(atEnd(iter));
        readOne(c, iter);
        SEQAN_ASSERT_EQ(c, buffer[i]);
        if (checkTell)
            SEQAN_ASSERT(position(iter) == i + 1 || position(iter) == -1);
    }
    // SEQAN_ASSERT_NOT(streamEof(stream));  // Not portable, inconsistent behaviour.
    readOne(c, iter);
    SEQAN_ASSERT_NOT(atEnd(iter));
    if (checkTell)
        SEQAN_ASSERT(position(iter) == 15 || position(iter) == -1);
}

// Test of streamPeek().
template <typename TStream>
void testStreamPeek(TStream & stream)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Output>::Type iter = directionIterator(stream, Output());

    SEQAN_ASSERT_EQ(position(iter), 0);
    SEQAN_ASSERT_EQ(*iter, 'T');
    SEQAN_ASSERT_EQ(*iter, 'T');
    SEQAN_ASSERT_EQ(*iter, 'T');
    SEQAN_ASSERT_EQ(position(iter), 0);
}

// Test of testStreamReadBlock(), limit is longer than stream.
template <typename TStream>
void testStreamReadBlockHitLimit(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Output>::Type iter = directionIterator(stream, Output());

    String<char> buffer;
    size_t charsRead = read(buffer, iter, 20);
    SEQAN_ASSERT_EQ(charsRead, 10u);
    SEQAN_ASSERT_EQ(buffer, "XXXXXXXXXX");
    if (checkTell)
        SEQAN_ASSERT(position(iter) == 10 || position(iter) == -1);
}

// Test of testStreamReadBlock(), limit is shorter than stream.
template <typename TStream>
void testStreamReadBlockHitNoLimit(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Input>::Type iter = directionIterator(stream, Input());

    String<char> buffer;
    SEQAN_ASSERT_EQ(read(buffer, iter, 5), 5);
    SEQAN_ASSERT(buffer == "XXXXX");
    if (checkTell)
        SEQAN_ASSERT(position(iter) == 5 || position(iter) == -1);
}

// Test of testStreamWriteBlock(), iterator interface.
template <typename TStream>
void testStreamWriteBlock(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Output>::Type iter = directionIterator(stream, Output());

    String<char> buffer;
    append(buffer, "ABCDEFGH");
    write(iter, buffer);
    if (checkTell)
        SEQAN_ASSERT(position(iter) == 8 || position(iter) == -1);
}

// Test of streamWriteChar().
template <typename TStream>
void testStreamWriteChar(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Output>::Type iter = directionIterator(stream, Output());

    writeValue(iter, '3');
    writeValue(iter, '4');
    writeValue(iter, '5');
    if (checkTell)
        SEQAN_ASSERT(position(iter) == 3 || position(iter) == -1);
}

template <typename TStream>
void testStreamPut(TStream & stream)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Output>::Type iter = directionIterator(stream, Output());

    writeValue(iter, 'c');
    write(iter, '\n');

    write(iter, CharString("seq"));
    writeValue(iter, '\n');

    write(iter, "sss");
    writeValue(iter, '\n');

    writeValue(iter, 12);
    writeValue(iter, '\n');

    writeValue(iter, 34u);
    writeValue(iter, '\n');

    writeValue(iter, 56l);
    writeValue(iter, '\n');

    writeValue(iter, 78ul);
    writeValue(iter, '\n');

    writeValue(iter, 5.4f);
    writeValue(iter, '\n');

    writeValue(iter, 6.5);
    writeValue(iter, '\n');

    writeValue(iter, Dna('A'));
    writeValue(iter, '\n');

    write(iter, DnaString("ACGT"));
    writeValue(iter, '\n');

    write(iter, Dna5String("ACGTN"));
    writeValue(iter, '\n');
}


// Test of readOne().
template <typename TStream>
void testStreamReadChar(TStream & stream, bool checkTell = true)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Input>::Type iter = directionIterator(stream, Input());

    char c;
    readOne(c, iter);
    SEQAN_ASSERT_EQ(c, '1');
    readOne(c, iter);
    SEQAN_ASSERT_EQ(c, '2');
    readOne(c, iter);
    SEQAN_ASSERT_EQ(c, '3');
    // SEQAN_ASSERT_NOT(streamEof(stream));  // Non-portable, different behaviour.
    readOne(c, iter);
    SEQAN_ASSERT(atEnd(iter));
    if (checkTell)
        SEQAN_ASSERT(position(iter) == 3 || position(iter) == -1);
}

// Test of streamSeek().
template <typename TStream>
void testStreamSeek(TStream & stream)
{
    using namespace seqan;

    typename DirectionIterator<TStream, Input>::Type iter = directionIterator(stream, Input());

    // TODO(holtgrew): Do real tests for all variants.
    // TODO(holtgrew): Use constant for origin in the future.
    SEQAN_ASSERT_EQ(static_cast<int>(position(iter)), 0);
    SEQAN_ASSERT_EQ(*iter, '0');
    SEQAN_ASSERT_EQ(static_cast<int>(position(iter)), 0);

    setPosition(iter, 4);
    SEQAN_ASSERT_EQ(static_cast<int>(position(iter)), 4);

    SEQAN_ASSERT_EQ(*iter, '4');
    SEQAN_ASSERT_EQ(static_cast<int>(position(iter)), 4);

    setPosition(iter, 2);
    SEQAN_ASSERT_EQ(static_cast<int>(position(iter)), 2);

    SEQAN_ASSERT_EQ(*iter, '2');
    SEQAN_ASSERT_EQ(static_cast<int>(position(iter)), 2);
}

#endif  // TEST_STREAM_TEST_STREAM_GENERIC_H_
