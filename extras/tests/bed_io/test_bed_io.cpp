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

#include <sstream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/bed_io.h>

SEQAN_DEFINE_TEST(test_bed_read_bed3_record)
{
    // Prepare in-memory data.
    std::stringstream ss;
    ss << "I\t123\t456\tsome data that is \tignored\n"
       << "II\t999\t1000\tdata again!";
    ss.seekg(0);

    // RecordReader to use.
    seqan::RecordReader<std::stringstream, seqan::SinglePass<> > reader(ss);

    // The record to load into.
    seqan::BedRecord<seqan::Bed3> record;

    // Perform tests.

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 122);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 998);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_read_bed3_record_with_context)
{
    // Prepare in-memory data.
    std::stringstream ss;
    ss << "I\t123\t456\tsome data that is \tignored\n"
       << "II\t999\t1000\tdata again!";
    ss.seekg(0);

    // RecordReader to use.
    seqan::RecordReader<std::stringstream, seqan::SinglePass<> > reader(ss);

    // The record to load into.
    seqan::BedRecord<seqan::Bed3> record;

    // The IO Context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    TNameStore refNames;
    seqan::NameStoreCache<TNameStore> refNamesCache(refNames);
    seqan::BedIOContext<TNameStore> bedIOContext(refNames, refNamesCache);

    // Perform tests.

    SEQAN_ASSERT_EQ(readRecord(record, reader, bedIOContext, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.rID, 0);
    SEQAN_ASSERT_EQ(record.beginPos, 122);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    SEQAN_ASSERT_EQ(length(refNames), 1u);
    SEQAN_ASSERT_EQ(refNames[0], "I");

    SEQAN_ASSERT_EQ(readRecord(record, reader, bedIOContext, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.rID, 1);
    SEQAN_ASSERT_EQ(record.beginPos, 998);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.data, "data again!");

    SEQAN_ASSERT_EQ(length(refNames), 2u);
    SEQAN_ASSERT_EQ(refNames[1], "II");
}

SEQAN_DEFINE_TEST(test_bed_read_bed4_record)
{
    // Prepare in-memory data.
    std::stringstream ss;
    ss << "I\t123\t456\tNAME\tsome data that is \tignored\n"
       << "II\t999\t1000\tNAME2\tdata again!";
    ss.seekg(0);

    // RecordReader to use.
    seqan::RecordReader<std::stringstream, seqan::SinglePass<> > reader(ss);

    // The record to load into.
    seqan::BedRecord<seqan::Bed4> record;

    // Perform tests.

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 122);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.name, "NAME");
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 998);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.name, "NAME2");
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_read_bed5_record)
{
    // Prepare in-memory data.
    std::stringstream ss;
    ss << "I\t123\t456\tNAME\t3\tsome data that is \tignored\n"
       << "II\t999\t1000\tNAME2\t2e5\tdata again!";
    ss.seekg(0);

    // RecordReader to use.
    seqan::RecordReader<std::stringstream, seqan::SinglePass<> > reader(ss);

    // The record to load into.
    seqan::BedRecord<seqan::Bed5> record;

    // Perform tests.

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 122);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.name, "NAME");
    SEQAN_ASSERT_EQ(record.score, "3");
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 998);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.name, "NAME2");
    SEQAN_ASSERT_EQ(record.score, "2e5");
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_read_bed6_record)
{
    // Prepare in-memory data.
    std::stringstream ss;
    ss << "I\t123\t456\tNAME\t3\t-\tsome data that is \tignored\n"
       << "II\t999\t1000\tNAME2\t2e5\t.\tdata again!";
    ss.seekg(0);

    // RecordReader to use.
    seqan::RecordReader<std::stringstream, seqan::SinglePass<> > reader(ss);

    // The record to load into.
    seqan::BedRecord<seqan::Bed6> record;

    // Perform tests.

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 122);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.name, "NAME");
    SEQAN_ASSERT_EQ(record.score, "3");
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 998);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.name, "NAME2");
    SEQAN_ASSERT_EQ(record.score, "2e5");
    SEQAN_ASSERT_EQ(record.strand, '.');
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_read_bed12_record)
{
    // Prepare in-memory data.
    std::stringstream ss;
    ss << "I\t123\t456\tNAME\t3\t-\t33\t66\t255,0,0\t3\t10,11,12\t1,2,3\tsome data that is \tignored\n"
       << "II\t999\t1000\tNAME2\t2e5\t.\t44\t55\t0,0,0\t3\t3,4,5\t4,5,6\tdata again!";
    ss.seekg(0);

    // RecordReader to use.
    seqan::RecordReader<std::stringstream, seqan::SinglePass<> > reader(ss);

    // The record to load into.
    seqan::BedRecord<seqan::Bed12> record;

    // Perform tests.

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 122);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.name, "NAME");
    SEQAN_ASSERT_EQ(record.score, "3");
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.thickBegin, 32);
    SEQAN_ASSERT_EQ(record.thickEnd, 66);
    SEQAN_ASSERT(record.itemRgb == seqan::BedRgb(255, 0, 0));
    SEQAN_ASSERT_EQ(record.blockCount, 3);
    SEQAN_ASSERT_EQ(length(record.blockSizes), 3u);
    SEQAN_ASSERT_EQ(record.blockSizes[0], 10);
    SEQAN_ASSERT_EQ(record.blockSizes[1], 11);
    SEQAN_ASSERT_EQ(record.blockSizes[2], 12);
    SEQAN_ASSERT_EQ(length(record.blockBegins), 3u);
    SEQAN_ASSERT_EQ(record.blockBegins[0], 0);
    SEQAN_ASSERT_EQ(record.blockBegins[1], 1);
    SEQAN_ASSERT_EQ(record.blockBegins[2], 2);
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 998);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.name, "NAME2");
    SEQAN_ASSERT_EQ(record.score, "2e5");
    SEQAN_ASSERT_EQ(record.strand, '.');
    SEQAN_ASSERT_EQ(record.thickBegin, 43);
    SEQAN_ASSERT_EQ(record.thickEnd, 55);
    SEQAN_ASSERT(record.itemRgb == seqan::BedRgb(0, 0, 0));
    SEQAN_ASSERT_EQ(record.blockCount, 3);
    SEQAN_ASSERT_EQ(length(record.blockSizes), 3u);
    SEQAN_ASSERT_EQ(record.blockSizes[0], 3);
    SEQAN_ASSERT_EQ(record.blockSizes[1], 4);
    SEQAN_ASSERT_EQ(record.blockSizes[2], 5);
    SEQAN_ASSERT_EQ(length(record.blockBegins), 3u);
    SEQAN_ASSERT_EQ(record.blockBegins[0], 3);
    SEQAN_ASSERT_EQ(record.blockBegins[1], 4);
    SEQAN_ASSERT_EQ(record.blockBegins[2], 5);
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_write_bed3_record)
{
    seqan::BedRecord<seqan::Bed3> record1;
    record1.ref = "I";
    record1.beginPos = 122;
    record1.endPos = 456;
    record1.data = "some data that is \tignored";

    seqan::BedRecord<seqan::Bed3> record2;
    record2.ref = "II";
    record2.beginPos = 999;
    record2.endPos = 1000;
    record2.data = "data again!";

    // Write BED records to string stream.
    std::stringstream ss;
    SEQAN_ASSERT_EQ(writeRecord(ss, record1, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(writeRecord(ss, record2, seqan::Bed()), 0);

    // Compar string stream to expected value.
    char const * EXPECTED =
            "I\t123\t456\tsome data that is \tignored\n"
            "II\t1000\t1000\tdata again!\n";
    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

SEQAN_DEFINE_TEST(test_bed_write_bed3_record_with_context)
{
    // The IO Context.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    TNameStore refNames;
    appendValue(refNames, "0");
    appendValue(refNames, "1");
    seqan::NameStoreCache<TNameStore> refNamesCache(refNames);
    seqan::BedIOContext<TNameStore> bedIOContext(refNames, refNamesCache);

    seqan::BedRecord<seqan::Bed3> record1;
    record1.ref = "I";
    record1.rID = 1;
    record1.beginPos = 122;
    record1.endPos = 456;
    record1.data = "some data that is \tignored";

    seqan::BedRecord<seqan::Bed3> record2;
    record2.ref = "II";
    record2.rID = 0;
    record2.beginPos = 999;
    record2.endPos = 1000;
    record2.data = "data again!";

    // Write BED records to string stream.
    std::stringstream ss;
    SEQAN_ASSERT_EQ(writeRecord(ss, record1, bedIOContext, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(writeRecord(ss, record2, bedIOContext, seqan::Bed()), 0);

    // Compar string stream to expected value.
    char const * EXPECTED =
            "1\t123\t456\tsome data that is \tignored\n"
            "0\t1000\t1000\tdata again!\n";
    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

SEQAN_DEFINE_TEST(test_bed_write_bed4_record)
{
    seqan::BedRecord<seqan::Bed4> record1;
    record1.ref = "I";
    record1.beginPos = 122;
    record1.endPos = 456;
    record1.name = "NAME1";
    record1.data = "some data that is \tignored";

    seqan::BedRecord<seqan::Bed4> record2;
    record2.ref = "II";
    record2.beginPos = 999;
    record2.endPos = 1000;
    record2.name = "NAME2";
    record2.data = "data again!";

    // Write BED records to string stream.
    std::stringstream ss;
    SEQAN_ASSERT_EQ(writeRecord(ss, record1, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(writeRecord(ss, record2, seqan::Bed()), 0);

    // Compar string stream to expected value.
    char const * EXPECTED =
            "I\t123\t456\tNAME1\tsome data that is \tignored\n"
            "II\t1000\t1000\tNAME2\tdata again!\n";
    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

SEQAN_DEFINE_TEST(test_bed_write_bed5_record)
{
    seqan::BedRecord<seqan::Bed5> record1;
    record1.ref = "I";
    record1.beginPos = 122;
    record1.endPos = 456;
    record1.name = "NAME1";
    record1.score = "5";
    record1.data = "some data that is \tignored";

    seqan::BedRecord<seqan::Bed5> record2;
    record2.ref = "II";
    record2.beginPos = 999;
    record2.endPos = 1000;
    record2.name = "NAME2";
    record2.score = "3e5";
    record2.data = "data again!";

    // Write BED records to string stream.
    std::stringstream ss;
    SEQAN_ASSERT_EQ(writeRecord(ss, record1, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(writeRecord(ss, record2, seqan::Bed()), 0);

    // Compar string stream to expected value.
    char const * EXPECTED =
            "I\t123\t456\tNAME1\t5\tsome data that is \tignored\n"
            "II\t1000\t1000\tNAME2\t3e5\tdata again!\n";
    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

SEQAN_DEFINE_TEST(test_bed_write_bed6_record)
{
    seqan::BedRecord<seqan::Bed6> record1;
    record1.ref = "I";
    record1.beginPos = 122;
    record1.endPos = 456;
    record1.name = "NAME1";
    record1.score = "5";
    record1.strand = '-';
    record1.data = "some data that is \tignored";

    seqan::BedRecord<seqan::Bed6> record2;
    record2.ref = "II";
    record2.beginPos = 999;
    record2.endPos = 1000;
    record2.name = "NAME2";
    record2.score = "3e5";
    record2.strand = '.';
    record2.data = "data again!";

    // Write BED records to string stream.
    std::stringstream ss;
    SEQAN_ASSERT_EQ(writeRecord(ss, record1, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(writeRecord(ss, record2, seqan::Bed()), 0);

    // Compar string stream to expected value.
    char const * EXPECTED =
            "I\t123\t456\tNAME1\t5\t-\tsome data that is \tignored\n"
            "II\t1000\t1000\tNAME2\t3e5\t.\tdata again!\n";
    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

SEQAN_DEFINE_TEST(test_bed_write_bed12_record)
{
    seqan::BedRecord<seqan::Bed12> record1;
    record1.ref = "I";
    record1.beginPos = 122;
    record1.endPos = 456;
    record1.name = "NAME1";
    record1.score = "5";
    record1.strand = '-';
    record1.thickBegin = 123;
    record1.thickEnd = 234;
    record1.itemRgb.red = 10;
    record1.itemRgb.green = 20;
    record1.itemRgb.blue = 30;
    record1.blockCount = 2;
    appendValue(record1.blockSizes, 10);
    appendValue(record1.blockSizes, 20);
    appendValue(record1.blockBegins, 3);
    appendValue(record1.blockBegins, 15);
    record1.data = "some data that is \tignored";

    seqan::BedRecord<seqan::Bed12> record2;
    record2.ref = "II";
    record2.beginPos = 999;
    record2.endPos = 1000;
    record2.name = "NAME2";
    record2.score = "3e5";
    record2.strand = '.';
    record2.thickBegin = 123;
    record2.thickEnd = 234;
    record2.itemRgb.red = 10;
    record2.itemRgb.green = 20;
    record2.itemRgb.blue = 30;
    record2.blockCount = 2;
    appendValue(record2.blockSizes, 10);
    appendValue(record2.blockSizes, 20);
    appendValue(record2.blockBegins, 3);
    appendValue(record2.blockBegins, 15);
    record2.data = "data again!";

    // Write BED records to string stream.
    std::stringstream ss;
    SEQAN_ASSERT_EQ(writeRecord(ss, record1, seqan::Bed()), 0);
    SEQAN_ASSERT_EQ(writeRecord(ss, record2, seqan::Bed()), 0);

    // Compar string stream to expected value.
    char const * EXPECTED =
            "I\t123\t456\tNAME1\t5\t-\t124\t234\t10,20,30\t2\t10,20\t2,14\tsome data that is \tignored\n"
            "II\t1000\t1000\tNAME2\t3e5\t.\t124\t234\t10,20,30\t2\t10,20\t2,14\tdata again!\n";
    SEQAN_ASSERT_EQ(ss.str(), EXPECTED);
}

SEQAN_BEGIN_TESTSUITE(test_bed_io)
{
    // Reading of BED records.
    SEQAN_CALL_TEST(test_bed_read_bed3_record);
    SEQAN_CALL_TEST(test_bed_read_bed3_record_with_context);
    SEQAN_CALL_TEST(test_bed_read_bed4_record);
    SEQAN_CALL_TEST(test_bed_read_bed5_record);
    SEQAN_CALL_TEST(test_bed_read_bed6_record);
    SEQAN_CALL_TEST(test_bed_read_bed12_record);

    // Writing of BED records.
    SEQAN_CALL_TEST(test_bed_write_bed3_record);
    SEQAN_CALL_TEST(test_bed_write_bed3_record_with_context);
    SEQAN_CALL_TEST(test_bed_write_bed4_record);
    SEQAN_CALL_TEST(test_bed_write_bed5_record);
    SEQAN_CALL_TEST(test_bed_write_bed6_record);
    SEQAN_CALL_TEST(test_bed_write_bed12_record);
}
SEQAN_END_TESTSUITE
