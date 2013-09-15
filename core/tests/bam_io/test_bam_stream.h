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

#ifndef CORE_TESTS_BAM_IO_TEST_EASY_BAM_IO_H_
#define CORE_TESTS_BAM_IO_TEST_EASY_BAM_IO_H_

#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

// ---------------------------------------------------------------------------
// Read Header
// ---------------------------------------------------------------------------

void testBamIOBamStreamReadHeader(char const * pathFragment)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragment);

    seqan::BamStream bamIO(toCString(filePath), seqan::BamStream::READ);
    SEQAN_ASSERT(isGood(bamIO));

    SEQAN_ASSERT_EQ(length(bamIO.header.records), 2u);
    SEQAN_ASSERT_EQ(bamIO.header.records[0].type, seqan::BAM_HEADER_FIRST);
    SEQAN_ASSERT_EQ(bamIO.header.records[0].tags[0].i1, "VN");
    SEQAN_ASSERT_EQ(bamIO.header.records[0].tags[0].i2, "1.3");
    SEQAN_ASSERT_EQ(bamIO.header.records[0].tags[1].i1, "SO");
    SEQAN_ASSERT_EQ(bamIO.header.records[0].tags[1].i2, "coordinate");
    SEQAN_ASSERT_EQ(bamIO.header.records[1].type, seqan::BAM_HEADER_REFERENCE);
    SEQAN_ASSERT_EQ(bamIO.header.records[1].tags[0].i1, "SN");
    SEQAN_ASSERT_EQ(bamIO.header.records[1].tags[0].i2, "REFERENCE");
    SEQAN_ASSERT_EQ(bamIO.header.records[1].tags[1].i1, "LN");
    SEQAN_ASSERT_EQ(bamIO.header.records[1].tags[1].i2, "10000");
    SEQAN_ASSERT_EQ(length(bamIO.header.sequenceInfos), 1u);
    SEQAN_ASSERT_EQ(bamIO.header.sequenceInfos[0].i1, "REFERENCE");
    SEQAN_ASSERT_EQ(bamIO.header.sequenceInfos[0].i2, 10000);
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_sam_read_header)
{
    testBamIOBamStreamReadHeader("/core/tests/bam_io/small.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_read_header)
{
    testBamIOBamStreamReadHeader("/core/tests/bam_io/small.bam");
}

// ---------------------------------------------------------------------------
// Read Records
// ---------------------------------------------------------------------------

void testBamIOBamStreamReadRecords(char const * pathFragment)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragment);

    seqan::BamStream bamIO(toCString(filePath));
    SEQAN_ASSERT(isGood(bamIO));

    seqan::BamAlignmentRecord record;

    seqan::String<seqan::BamAlignmentRecord> alignments;
    while (!atEnd(bamIO))
    {
        resize(alignments, length(alignments) + 1);
        SEQAN_ASSERT_EQ(readRecord(back(alignments), bamIO), 0);
        SEQAN_ASSERT(isGood(bamIO));
    }
    SEQAN_ASSERT_EQ(length(alignments), 3u);

    SEQAN_ASSERT_EQ(alignments[0].qName, "READ0");
    SEQAN_ASSERT_EQ(alignments[0].flag, 2);
    SEQAN_ASSERT_EQ(alignments[0].rID, 0);
    SEQAN_ASSERT_EQ(alignments[0].beginPos, 0);
    SEQAN_ASSERT_EQ(alignments[0].mapQ, 8);
    SEQAN_ASSERT_EQ(length(alignments[0].cigar), 3u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[0].count, 5u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[0].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[0].cigar[1].count, 1u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[1].operation, 'I');
    SEQAN_ASSERT_EQ(alignments[0].cigar[2].count, 4u);
    SEQAN_ASSERT_EQ(alignments[0].cigar[2].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[0].rNextId, 0);
    SEQAN_ASSERT_EQ(alignments[0].pNext, 30);
    SEQAN_ASSERT_EQ(alignments[0].tLen, 40);
    SEQAN_ASSERT_EQ(alignments[0].seq, "AAAAAAAAAA");
    SEQAN_ASSERT_EQ(alignments[0].qual, "!!!!!!!!!!");
    SEQAN_ASSERT_EQ(length(alignments[0].tags), 0u);

    SEQAN_ASSERT_EQ(alignments[1].qName, "READ0");
    SEQAN_ASSERT_EQ(alignments[1].flag, 1);
    SEQAN_ASSERT_EQ(alignments[1].rID, 0);
    SEQAN_ASSERT_EQ(alignments[1].beginPos, 1);
    SEQAN_ASSERT_EQ(alignments[1].mapQ, 8);
    SEQAN_ASSERT_EQ(length(alignments[1].cigar), 3u);
    SEQAN_ASSERT_EQ(alignments[1].cigar[0].count, 5u);
    SEQAN_ASSERT_EQ(alignments[1].cigar[0].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[1].cigar[1].count, 1u);
    SEQAN_ASSERT_EQ(alignments[1].cigar[1].operation, 'I');
    SEQAN_ASSERT_EQ(alignments[1].cigar[2].count, 4u);
    SEQAN_ASSERT_EQ(alignments[1].cigar[2].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[1].rNextId, 0);
    SEQAN_ASSERT_EQ(alignments[1].pNext, 30);
    SEQAN_ASSERT_EQ(alignments[1].tLen, 40);
    SEQAN_ASSERT_EQ(alignments[1].seq, "AAAAAAAAAA");
    SEQAN_ASSERT_EQ(alignments[1].qual, "!!!!!!!!!!");
    SEQAN_ASSERT_EQ(length(alignments[1].tags), 0u);

    SEQAN_ASSERT_EQ(alignments[2].qName, "READ0");
    SEQAN_ASSERT_EQ(alignments[2].flag, 3);
    SEQAN_ASSERT_EQ(alignments[2].rID, 0);
    SEQAN_ASSERT_EQ(alignments[2].beginPos, 2);
    SEQAN_ASSERT_EQ(alignments[2].mapQ, 8);
    SEQAN_ASSERT_EQ(length(alignments[2].cigar), 3u);
    SEQAN_ASSERT_EQ(alignments[2].cigar[0].count, 5u);
    SEQAN_ASSERT_EQ(alignments[2].cigar[0].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[2].cigar[1].count, 1u);
    SEQAN_ASSERT_EQ(alignments[2].cigar[1].operation, 'I');
    SEQAN_ASSERT_EQ(alignments[2].cigar[2].count, 4u);
    SEQAN_ASSERT_EQ(alignments[2].cigar[2].operation, 'M');
    SEQAN_ASSERT_EQ(alignments[2].rNextId, -1);
    SEQAN_ASSERT_EQ(alignments[2].pNext, 2147483647);
    SEQAN_ASSERT_EQ(alignments[2].tLen, 0);
    SEQAN_ASSERT_EQ(alignments[2].seq, "AAAAAAAAAA");
    SEQAN_ASSERT_EQ(alignments[2].qual, "!!!!!!!!!!");
    SEQAN_ASSERT_EQ(length(alignments[2].tags), 0u);

    SEQAN_ASSERT_EQ(bamIO._nameStore[0], "REFERENCE");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_sam_read_records)
{
    testBamIOBamStreamReadRecords("/core/tests/bam_io/small.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_read_records)
{
    testBamIOBamStreamReadRecords("/core/tests/bam_io/small.bam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_read_ex1)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/core/tests/bam_io/ex1.bam");

    seqan::BamStream bamIO(toCString(filePath));
    SEQAN_ASSERT(isGood(bamIO));

    SEQAN_ASSERT_EQ(nameStore(bamIO.bamIOContext)[0], "seq1");
    SEQAN_ASSERT_EQ(nameStore(bamIO.bamIOContext)[1], "seq2");
    SEQAN_ASSERT_EQ(bamIO.bamIOContext.translateFile2GlobalRefId[0], 0u);
    SEQAN_ASSERT_EQ(bamIO.bamIOContext.translateFile2GlobalRefId[1], 1u);

    seqan::BamAlignmentRecord record;
    seqan::String<int> counts;
    resize(counts, 2, 0);
    while (!atEnd(bamIO))
    {
        SEQAN_ASSERT_EQ(readRecord(record, bamIO), 0);
        ++counts[record.rID];
//        seqan::CharString name = nameStore(bamIO.bamIOContext)[record.rID];
//        std::cout << "Chrom: " << name << " (" << record.rID << ")" << std::endl;
    }
    SEQAN_ASSERT_EQ(counts[0], 1501);
    SEQAN_ASSERT_EQ(counts[1], 1806);
}

// ---------------------------------------------------------------------------
// Write Header
// ---------------------------------------------------------------------------

void testBamIOBamStreamWriteHeader(char const * pathFragmentExpected)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragmentExpected);

    seqan::CharString tmpPath = SEQAN_TEMP_FILENAME();
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        append(tmpPath, ".bam");
    else
        append(tmpPath, ".sam");

    // Initialize BamStream, build header.
    seqan::BamStream bamIO(toCString(tmpPath), seqan::BamStream::WRITE);
    resize(bamIO.header.sequenceInfos, 1);
    bamIO.header.sequenceInfos[0].i1 = "REFERENCE";
    bamIO.header.sequenceInfos[0].i2 = 10000;
    resize(bamIO.header.records, 2);
    resize(bamIO.header.records[0].tags, 2);
    bamIO.header.records[0].type = seqan::BAM_HEADER_FIRST;
    bamIO.header.records[0].tags[0].i1 = "VN";
    bamIO.header.records[0].tags[0].i2 = "1.3";
    bamIO.header.records[0].tags[1].i1 = "SO";
    bamIO.header.records[0].tags[1].i2 = "coordinate";
    resize(bamIO.header.records[1].tags, 2);
    bamIO.header.records[1].type = seqan::BAM_HEADER_REFERENCE;
    bamIO.header.records[1].tags[0].i1 = "SN";
    bamIO.header.records[1].tags[0].i2 = "REFERENCE";
    bamIO.header.records[1].tags[1].i1 = "LN";
    bamIO.header.records[1].tags[1].i2 = "10000";

    // Force writing of header on flush.
    flush(bamIO);
    close(bamIO);

    // Compare results.
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        SEQAN_ASSERT(seqan::_compareBinaryFiles(toCString(tmpPath), toCString(filePath)));
    else
        SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(filePath)));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_sam_write_header)
{
    testBamIOBamStreamWriteHeader("/core/tests/bam_io/header.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_write_header)
{
    testBamIOBamStreamWriteHeader("/core/tests/bam_io/header.bam");
}

// ---------------------------------------------------------------------------
// Write Records
// ---------------------------------------------------------------------------

void testBamIOBamStreamWriteRecords(char const * pathFragmentExpected)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragmentExpected);

    seqan::CharString tmpPath = SEQAN_TEMP_FILENAME();
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        append(tmpPath, ".bam");
    else
        append(tmpPath, ".sam");

    // Initialize BamStream, build header.
    seqan::BamStream bamIO(toCString(tmpPath), seqan::BamStream::WRITE);
    resize(bamIO.header.sequenceInfos, 1);
    bamIO.header.sequenceInfos[0].i1 = "REFERENCE";
    bamIO.header.sequenceInfos[0].i2 = 10000;
    resize(bamIO.header.records, 2);
    resize(bamIO.header.records[0].tags, 2);
    bamIO.header.records[0].type = seqan::BAM_HEADER_FIRST;
    bamIO.header.records[0].tags[0].i1 = "VN";
    bamIO.header.records[0].tags[0].i2 = "1.3";
    bamIO.header.records[0].tags[1].i1 = "SO";
    bamIO.header.records[0].tags[1].i2 = "coordinate";
    resize(bamIO.header.records[1].tags, 2);
    bamIO.header.records[1].type = seqan::BAM_HEADER_REFERENCE;
    bamIO.header.records[1].tags[0].i1 = "SN";
    bamIO.header.records[1].tags[0].i2 = "REFERENCE";
    bamIO.header.records[1].tags[1].i1 = "LN";
    bamIO.header.records[1].tags[1].i2 = "10000";

    // Construct first records.
    seqan::BamAlignmentRecord record;

    record.qName = "READ0";
    record.flag = 2;
    record.rID = 0;
    record.beginPos = 0;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = 0;
    record.pNext = 30;
    record.tLen = 40;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    SEQAN_ASSERT_EQ(writeRecord(bamIO, record), 0);

    record.qName = "READ0";
    record.flag = 1;
    record.rID = 0;
    record.beginPos = 1;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = 0;
    record.pNext = 30;
    record.tLen = 40;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    SEQAN_ASSERT_EQ(writeRecord(bamIO, record), 0);

    record.qName = "READ0";
    record.flag = 3;
    record.rID = 0;
    record.beginPos = 2;
    record.mapQ = 8;
    resize(record.cigar, 3);
    record.cigar[0].count = 5;
    record.cigar[0].operation = 'M';
    record.cigar[1].count = 1;
    record.cigar[1].operation = 'I';
    record.cigar[2].count = 4;
    record.cigar[2].operation = 'M';
    record.rNextId = seqan::BamAlignmentRecord::INVALID_REFID;
    record.pNext = seqan::BamAlignmentRecord::INVALID_POS;
    record.tLen = seqan::BamAlignmentRecord::INVALID_LEN;
    record.seq = "AAAAAAAAAA";
    record.qual = "!!!!!!!!!!";
    SEQAN_ASSERT_EQ(writeRecord(bamIO, record), 0);

    // Force writing of everything.
    close(bamIO);

    // Compare results.
    if (seqan::endsWith(pathFragmentExpected, ".bam"))
        SEQAN_ASSERT(seqan::_compareBinaryFiles(toCString(tmpPath), toCString(filePath)));
    else
        SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(filePath)));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_sam_write_records)
{
    testBamIOBamStreamWriteRecords("/core/tests/bam_io/small.sam");
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_write_records)
{
    testBamIOBamStreamWriteRecords("/core/tests/bam_io/small.bam");
}

// ---------------------------------------------------------------------------
// File Size / Byte Position In File
// ---------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_sam_file_size)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/core/tests/bam_io/small.sam");

    seqan::BamStream bamStream(toCString(filePath));
    SEQAN_ASSERT(isGood(bamStream));

    // TODO(holtgrew): tellg() on std::istream does not work correctly for some reason.

    SEQAN_ASSERT_EQ(fileSize(bamStream), 226u);
    // SEQAN_ASSERT_EQ(positionInFile(bamStream), 0u);

    seqan::BamAlignmentRecord record;
    SEQAN_ASSERT_EQ(readRecord(record, bamStream), 0);

    // SEQAN_ASSERT_EQ(positionInFile(bamStream), 0u);
}

SEQAN_DEFINE_TEST(test_bam_io_bam_stream_bam_file_size)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/core/tests/bam_io/small.bam");

    seqan::BamStream bamStream(toCString(filePath));
    SEQAN_ASSERT(isGood(bamStream));

    SEQAN_ASSERT_EQ(fileSize(bamStream), 181u);
    SEQAN_ASSERT_EQ(positionInFile(bamStream), 0u);

    seqan::BamAlignmentRecord record;
    SEQAN_ASSERT_EQ(readRecord(record, bamStream), 0);

    SEQAN_ASSERT_EQ(positionInFile(bamStream), 0u);  // Is block position.
}

#endif  // #ifndef CORE_TESTS_BAM_IO_TEST_EASY_BAM_IO_H_
