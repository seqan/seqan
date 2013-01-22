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

#include <sstream>

#include <seqan/stream.h>

SEQAN_DEFINE_TEST(test_stream_record_reader_single_pass_position)
{
    // Prepare input stream.
    std::stringstream ss;
    ss << "1234567890";
    ss.seekg(0);
    ss.seekp(0);
    ss.clear();

    // Create RecordReader object.
    typedef seqan::RecordReader<std::stringstream, seqan::SinglePass<> > TRecordReader;
    TRecordReader reader(ss, /*bufferSize=*/ 3);

    typedef typename seqan::Position<TRecordReader>::Type TPos;

    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 3), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)4);
    SEQAN_ASSERT_EQ(skipNChars(reader, 6), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
    SEQAN_ASSERT(atEnd(reader));
}

SEQAN_DEFINE_TEST(test_stream_record_reader_single_pass_set_position)
{
    // Prepare input stream.
    std::stringstream ss;
    ss << "1234567890";
    ss.seekg(0);
    ss.seekp(0);
    ss.clear();

    // Create RecordReader object.
    typedef seqan::RecordReader<std::stringstream, seqan::SinglePass<> > TRecordReader;
    TRecordReader reader(ss, /*bufferSize=*/ 3);

    typedef typename seqan::Position<TRecordReader>::Type TPos;

    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 3), 0);
    SEQAN_ASSERT_EQ(setPosition(reader, 1), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 6), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)7);
    SEQAN_ASSERT_EQ(setPosition(reader, 10), 0);
    SEQAN_ASSERT(atEnd(reader));
}

SEQAN_DEFINE_TEST(test_stream_record_reader_double_pass_position)
{
    // Prepare input stream.
    std::stringstream ss;
    ss << "1234567890";
    ss.seekg(0);
    ss.seekp(0);
    ss.clear();

    // Create RecordReader object.
    typedef seqan::RecordReader<std::stringstream, seqan::DoublePass<> > TRecordReader;
    TRecordReader reader(ss, /*bufferSize=*/ 3);

    typedef typename seqan::Position<TRecordReader>::Type TPos;

    startFirstPass(reader);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 3), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)4);

    startSecondPass(reader);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 3), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)4);

    startFirstPass(reader);
    SEQAN_ASSERT_EQ(skipNChars(reader, 6), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
    SEQAN_ASSERT(atEnd(reader));

    startSecondPass(reader);
    SEQAN_ASSERT_EQ(skipNChars(reader, 6), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
    SEQAN_ASSERT(atEnd(reader));
}

SEQAN_DEFINE_TEST(test_stream_record_reader_double_pass_set_position)
{
    // Prepare input stream.
    std::stringstream ss;
    ss << "1234567890";
    ss.seekg(0);
    ss.seekp(0);
    ss.clear();

    // Create RecordReader object.
    typedef seqan::RecordReader<std::stringstream, seqan::DoublePass<> > TRecordReader;
    TRecordReader reader(ss, /*bufferSize=*/ 3);

    typedef typename seqan::Position<TRecordReader>::Type TPos;

    startFirstPass(reader);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    startSecondPass(reader);

    setPosition(reader, 0);
    SEQAN_ASSERT_EQ(reader._passNo, 1);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_EQ(value(reader), '1');

    SEQAN_ASSERT_EQ(skipNChars(reader, 10), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
    SEQAN_ASSERT(atEnd(reader));

    startSecondPass(reader);
    SEQAN_ASSERT_EQ(reader._passNo, 2);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_EQ(skipNChars(reader, 10), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_single_pass_mmap_position)
{
    // Prepare input stream.
    seqan::CharString s = "1234567890";

    // Create RecordReader object.
    typedef seqan::RecordReader<seqan::CharString, seqan::SinglePass<seqan::Mapped> > TRecordReader;
    TRecordReader reader(s);

    typedef typename seqan::Position<TRecordReader>::Type TPos;

    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 3), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)4);
    SEQAN_ASSERT_EQ(skipNChars(reader, 6), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
    SEQAN_ASSERT(atEnd(reader));
}

SEQAN_DEFINE_TEST(test_stream_record_reader_single_pass_mmap_set_position)
{
    // Prepare input stream.
    seqan::CharString s = "1234567890";

    // Create RecordReader object.
    typedef seqan::RecordReader<seqan::CharString, seqan::SinglePass<seqan::Mapped> > TRecordReader;
    TRecordReader reader(s);

    typedef typename seqan::Position<TRecordReader>::Type TPos;

    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 3), 0);
    SEQAN_ASSERT_EQ(setPosition(reader, 1), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 6), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)7);
    SEQAN_ASSERT_EQ(setPosition(reader, 10), 0);
    SEQAN_ASSERT(atEnd(reader));
}

SEQAN_DEFINE_TEST(test_stream_record_reader_double_pass_mmap_position)
{
    // Prepare input stream.
    seqan::CharString s = "1234567890";

    // Create RecordReader object.
    typedef seqan::RecordReader<seqan::CharString, seqan::DoublePass<seqan::Mapped> > TRecordReader;
    TRecordReader reader(s);

    typedef typename seqan::Position<TRecordReader>::Type TPos;

    startFirstPass(reader);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 3), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)4);

    startSecondPass(reader);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    SEQAN_ASSERT_EQ(skipNChars(reader, 3), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)4);

    startFirstPass(reader);
    SEQAN_ASSERT_EQ(skipNChars(reader, 6), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
    SEQAN_ASSERT(atEnd(reader));

    startSecondPass(reader);
    SEQAN_ASSERT_EQ(skipNChars(reader, 6), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
    SEQAN_ASSERT(atEnd(reader));
}

SEQAN_DEFINE_TEST(test_stream_record_reader_double_pass_mmap_set_position)
{
    // Prepare input stream.
    seqan::CharString s = "1234567890";

    // Create RecordReader object.
    typedef seqan::RecordReader<seqan::CharString, seqan::DoublePass<seqan::Mapped> > TRecordReader;
    TRecordReader reader(s);

    typedef typename seqan::Position<TRecordReader>::Type TPos;

    startFirstPass(reader);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_NOT(goNext(reader));
    SEQAN_ASSERT_EQ(position(reader), (TPos)1);
    startSecondPass(reader);

    setPosition(reader, 0);
    SEQAN_ASSERT_EQ(reader._passNo, 1);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_EQ(value(reader), '1');

    SEQAN_ASSERT_EQ(skipNChars(reader, 10), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
    SEQAN_ASSERT(atEnd(reader));

    startSecondPass(reader);
    SEQAN_ASSERT_EQ(reader._passNo, 2);
    SEQAN_ASSERT_EQ(position(reader), (TPos)0);
    SEQAN_ASSERT_EQ(skipNChars(reader, 10), 0);
    SEQAN_ASSERT_EQ(position(reader), (TPos)10);
}
