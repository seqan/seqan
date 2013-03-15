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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// FASTQ-Tests for seqan/stream/read_fasta_fastq.h
// ==========================================================================


#ifndef TEST_STREAM_TEST_STREAM_READ_FASTQ_H_
#define TEST_STREAM_TEST_STREAM_READ_FASTQ_H_

template <typename TRecordReader>
void FASTQ_TEST(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(guessStreamFormat(reader, Fastq()));

    CharString meta;
    DnaString seq;
    CharString qual;

    // qualities are not asked for, but should be skipped internally
    int res = readRecord(meta, seq, reader, Fastq());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, "SEQ_ID");
    SEQAN_ASSERT_EQ(seq,
                "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(value(reader), '@');

    // now qualities are also requested
    res = readRecord(meta, seq, qual, reader, Fastq());
    SEQAN_ASSERT_EQ(res, 0);
    SEQAN_ASSERT_EQ(meta, " 2ndSequence with formatting obscurities");
    SEQAN_ASSERT_EQ(seq,
                "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(qual, "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");
    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void FASTQ_TEST_BATCH(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(guessStreamFormat(reader, Fastq()));

    StringSet<CharString> metas;
    StringSet<String<DnaQ> > seqs;

    int res = read2(metas, seqs, reader, Fastq());
    SEQAN_ASSERT_EQ(res, 0);

    CharString quals;

    SEQAN_ASSERT_EQ(metas[0], "SEQ_ID");
    SEQAN_ASSERT_EQ(seqs[0], "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(metas[1], " 2ndSequence with formatting obscurities");
    SEQAN_ASSERT_EQ(seqs[1], "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT(atEnd(reader));

    resize(quals, length(seqs[0]));
    assignQualities(quals, seqs[0]);
    SEQAN_ASSERT_EQ(quals, "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");
    resize(quals, length(seqs[1]));
    assignQualities(quals, seqs[1]);
    SEQAN_ASSERT_EQ(quals, "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");

    SEQAN_ASSERT_EQ(length(metas), 2u);
    SEQAN_ASSERT_EQ(length(seqs), 2u);
}

template <typename TRecordReader>
void FASTQ_TEST_BATCH_CONCAT(TRecordReader & reader)
{
    using namespace seqan;

    SEQAN_ASSERT(guessStreamFormat(reader, Fastq()));

    StringSet<CharString, Owner<ConcatDirect<> > > metas;
    StringSet<String<Dna5Q>, Owner<ConcatDirect<> > > seqs;
    // StringSet<CharString, Owner<ConcatDirect<> > > quals;

    
    int res = read2(metas, seqs, reader, Fastq());
    SEQAN_ASSERT_EQ(res, 0);

    CharString quals;

    SEQAN_ASSERT_EQ(metas[0], "SEQ_ID");
    SEQAN_ASSERT_EQ(seqs[0], "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT_EQ(metas[1], " 2ndSequence with formatting obscurities");
    SEQAN_ASSERT_EQ(seqs[1], "GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT");
    SEQAN_ASSERT(atEnd(reader));

    resize(quals, length(seqs[0]));
    assignQualities(quals, seqs[0]);
    SEQAN_ASSERT_EQ(quals, "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");
    resize(quals, length(seqs[1]));
    assignQualities(quals, seqs[1]);
    SEQAN_ASSERT_EQ(quals, "!''*((((***+))%%%++)(%%%%).@***-+*''))**55CCF>>>>>>CCCCCCC65");

    SEQAN_ASSERT_EQ(length(metas), 2u);
    SEQAN_ASSERT_EQ(length(seqs), 2u);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_single_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    typedef RecordReader<std::fstream, SinglePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTQ_TEST(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_double_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    typedef RecordReader<std::fstream, DoublePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTQ_TEST(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_batch_single_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    typedef RecordReader<std::fstream, DoublePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTQ_TEST_BATCH(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_batch_double_fstream)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    typedef RecordReader<std::fstream, DoublePass<void> > TRecordReader;
    TRecordReader reader(*file);

    FASTQ_TEST_BATCH(reader);

    file->close();
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_single_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, SinglePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_double_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_batch_single_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, SinglePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST_BATCH(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_batch_single_concat_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, SinglePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST_BATCH_CONCAT(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_batch_double_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST_BATCH(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_batch_double_concat_mmap)
{
    using namespace seqan;
    CharString filename;
    std::fstream *file = createFastQFile(filename);

    file->close();
    String<char, MMap<> > mmapString;
    open(mmapString, toCString(filename));

    typedef RecordReader<String<char, MMap<> >, DoublePass<Mapped> > TRecordReader;
    TRecordReader reader(mmapString);

    FASTQ_TEST_BATCH_CONCAT(reader);

    close(mmapString);
}

SEQAN_DEFINE_TEST(test_stream_record_reader_fastq_check_stream_format)
{
    using namespace seqan;

    CharString fileName = SEQAN_PATH_TO_ROOT();
    append(fileName, "/core/tests/stream/files/lane_5_p1.fastq");
    std::fstream inFile(toCString(fileName), std::ios::in | std::ios::binary);

    RecordReader<std::fstream, SinglePass<> > reader(inFile);
    AutoSeqStreamFormat tagSelector;
    bool b = guessStreamFormat(reader, tagSelector);

    SEQAN_ASSERT_MSG(b, "File format detection must have been successful.");
    SEQAN_ASSERT_EQ_MSG(tagSelector.tagId, 2, "Format must be FASTQ.");
}

#endif // def TEST_STREAM_TEST_STREAM_READ_FASTQ_H_
