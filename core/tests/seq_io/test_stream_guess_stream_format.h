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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Tests for seqan/stream/guess_stream_format.h
// NOTE that tests for guessStreamFormat() specific to each FileFormat are
// located with the other tests for that file format, e.g.
// tests/stream/test_stream_record_reader_fasta.h for FASTA
// ==========================================================================


#ifndef TEST_STREAM_TEST_STREAM_GUESS_STREAM_FORMAT_H_
#define TEST_STREAM_TEST_STREAM_GUESS_STREAM_FORMAT_H_

SEQAN_DEFINE_TEST(test_stream_guess_stream_format_auto_fasta)
{
    using namespace seqan;

    CharString filename;
    std::fstream *file = createFastAFile(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);
    AutoSeqStreamFormat tagSelect;

    SEQAN_ASSERT(guessStreamFormat(reader, tagSelect));
    SEQAN_ASSERT_EQ(value(tagSelect), +(Find<AutoSeqStreamFormat, Fasta>::VALUE));

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_guess_stream_format_auto_fastq)
{
    using namespace seqan;

    CharString filename;
    std::fstream *file = createFastQFile(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);
    AutoSeqStreamFormat tagSelect;

    SEQAN_ASSERT(guessStreamFormat(reader, tagSelect));
    SEQAN_ASSERT_EQ(value(tagSelect), +(Find<AutoSeqStreamFormat, Fastq>::VALUE));

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_guess_stream_format_auto_bogus)
{
    using namespace seqan;

    CharString filename;
    std::fstream *file = createBogusFile(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);
    AutoSeqStreamFormat tagSelect;

    SEQAN_ASSERT_NOT(guessStreamFormat(reader, tagSelect));
    SEQAN_ASSERT_EQ(tagSelect.tagId, -1);

    file->close();
}

#endif // def TEST_STREAM_TEST_STREAM_GUESS_STREAM_FORMAT_H_
