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

#include <sstream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/bed_io.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_bed_read_bed3_record)
{
    // Prepare in-memory data.
    String<char> test = "I\t123\t456\tsome data that is \tignored\n";
    append(test, "II\t999\t1000\tdata again!");

    DirectionIterator<String<char>, Input>::Type iter = directionIterator(test, Input());

    // The record to load into.
    seqan::BedRecord<seqan::Bed3> record;
    seqan::CharString buffer;

    // Perform tests.
    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 123);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 999);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_read_bed4_record)
{
    String<char> test = "I\t123\t456\tNAME\tsome data that is \tignored\n";
    append(test, "II\t999\t1000\tNAME2\tdata again!");

    // Iterator to use
    DirectionIterator<String<char>, Input>::Type iter = directionIterator(test, Input());

    // The record to load into.
    seqan::BedRecord<seqan::Bed4> record;
    seqan::CharString buffer;

    // Perform tests.

    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 123);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.name, "NAME");
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 999);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.name, "NAME2");
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_read_bed5_record)
{
    // Prepare in-memory data.
    String<char> test = "I\t123\t456\tNAME\t3\tsome data that is \tignored\n";
    append(test, "II\t999\t1000\tNAME2\t2e5\tdata again!");

    // Iterator to use
    DirectionIterator<String<char>, Input>::Type iter = directionIterator(test, Input());

    // The record to load into.
    seqan::BedRecord<seqan::Bed5> record;
    CharString buffer;

    // Perform tests.

    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 123);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.name, "NAME");
    SEQAN_ASSERT_EQ(record.score, "3");
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 999);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.name, "NAME2");
    SEQAN_ASSERT_EQ(record.score, "2e5");
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_read_bed6_record)
{
    // Prepare in-memory data.
    String<char> test = "I\t123\t456\tNAME\t3\t-\tsome data that is \tignored\n";
    append(test, "II\t999\t1000\tNAME2\t2e5\t.\tdata again!");

    // Iterator to use
    DirectionIterator<String<char>, Input>::Type iter = directionIterator(test, Input());

    // The record to load into.
    seqan::BedRecord<seqan::Bed6> record;
    CharString buffer;

    // Perform tests.

    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 123);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.name, "NAME");
    SEQAN_ASSERT_EQ(record.score, "3");
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 999);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.name, "NAME2");
    SEQAN_ASSERT_EQ(record.score, "2e5");
    SEQAN_ASSERT_EQ(record.strand, '.');
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_read_bed12_record)
{
    // Prepare in-memory data.
    String<char> test = "I\t123\t456\tNAME\t3\t-\t33\t66\t255,0,0\t3\t10,11,12\t1,2,3\tsome data that is \tignored\n";
    append(test, "II\t999\t1000\tNAME2\t2e5\t.\t44\t55\t0,0,0\t3\t3,4,5\t4,5,6\tdata again!");

    // Iterator to use
    DirectionIterator<String<char>, Input>::Type iter = directionIterator(test, Input());

    // The record to load into.
    seqan::BedRecord<seqan::Bed12> record;
    CharString buffer;

    // Perform tests.

    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 123);
    SEQAN_ASSERT_EQ(record.endPos, 456);
    SEQAN_ASSERT_EQ(record.name, "NAME");
    SEQAN_ASSERT_EQ(record.score, "3");
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.thickBegin, 33);
    SEQAN_ASSERT_EQ(record.thickEnd, 66);
    SEQAN_ASSERT(record.itemRgb == seqan::BedRgb(255, 0, 0));
    SEQAN_ASSERT_EQ(record.blockCount, 3);
    SEQAN_ASSERT_EQ(length(record.blockSizes), 3u);
    SEQAN_ASSERT_EQ(record.blockSizes[0], 10);
    SEQAN_ASSERT_EQ(record.blockSizes[1], 11);
    SEQAN_ASSERT_EQ(record.blockSizes[2], 12);
    SEQAN_ASSERT_EQ(length(record.blockBegins), 3u);
    SEQAN_ASSERT_EQ(record.blockBegins[0], 1);
    SEQAN_ASSERT_EQ(record.blockBegins[1], 2);
    SEQAN_ASSERT_EQ(record.blockBegins[2], 3);
    SEQAN_ASSERT_EQ(record.data, "some data that is \tignored");

    readRecord(record, buffer, iter, seqan::Bed());
    SEQAN_ASSERT_EQ(record.ref, "II");
    SEQAN_ASSERT_EQ(record.beginPos, 999);
    SEQAN_ASSERT_EQ(record.endPos, 1000);
    SEQAN_ASSERT_EQ(record.name, "NAME2");
    SEQAN_ASSERT_EQ(record.score, "2e5");
    SEQAN_ASSERT_EQ(record.strand, '.');
    SEQAN_ASSERT_EQ(record.thickBegin, 44);
    SEQAN_ASSERT_EQ(record.thickEnd, 55);
    SEQAN_ASSERT(record.itemRgb == seqan::BedRgb(0, 0, 0));
    SEQAN_ASSERT_EQ(record.blockCount, 3);
    SEQAN_ASSERT_EQ(length(record.blockSizes), 3u);
    SEQAN_ASSERT_EQ(record.blockSizes[0], 3);
    SEQAN_ASSERT_EQ(record.blockSizes[1], 4);
    SEQAN_ASSERT_EQ(record.blockSizes[2], 5);
    SEQAN_ASSERT_EQ(length(record.blockBegins), 3u);
    SEQAN_ASSERT_EQ(record.blockBegins[0], 4);
    SEQAN_ASSERT_EQ(record.blockBegins[1], 5);
    SEQAN_ASSERT_EQ(record.blockBegins[2], 6);
    SEQAN_ASSERT_EQ(record.data, "data again!");
}

SEQAN_DEFINE_TEST(test_bed_write_bed3_record)
{
    seqan::BedRecord<seqan::Bed3> record1;
    record1.ref = "I";
    record1.beginPos = 123;
    record1.endPos = 456;
    record1.data = "some data that is \tignored";

    seqan::BedRecord<seqan::Bed3> record2;
    record2.ref = "II";
    record2.beginPos = 999;
    record2.endPos = 1000;
    record2.data = "data again!";

    // Write BED records to string stream.
    String<char> out;
    writeRecord(out, record1, seqan::Bed());
    writeRecord(out, record2, seqan::Bed());

    // Compar string stream to expected value.
    String<char> expected = "I\t123\t456\tsome data that is \tignored\n";
    append(expected, "II\t999\t1000\tdata again!\n");
    SEQAN_ASSERT_EQ(out, expected);
}

SEQAN_DEFINE_TEST(test_bed_write_bed4_record)
{
    seqan::BedRecord<seqan::Bed4> record1;
    record1.ref = "I";
    record1.beginPos = 123;
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
    String<char> out;
    writeRecord(out, record1, seqan::Bed());
    writeRecord(out, record2, seqan::Bed());

    // Compar string stream to expected value.
    String<char> expected = "I\t123\t456\tNAME1\tsome data that is \tignored\n";
    append(expected, "II\t999\t1000\tNAME2\tdata again!\n");
    SEQAN_ASSERT_EQ(out, expected);
}

SEQAN_DEFINE_TEST(test_bed_write_bed5_record)
{
    seqan::BedRecord<seqan::Bed5> record1;
    record1.ref = "I";
    record1.beginPos = 123;
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

    // Write BED records to string stream.String<char> out;
    String<char> out;
    writeRecord(out, record1, seqan::Bed());
    writeRecord(out, record2, seqan::Bed());

    // Compar string stream to expected value.
    String<char> expected = "I\t123\t456\tNAME1\t5\tsome data that is \tignored\n";
    append(expected, "II\t999\t1000\tNAME2\t3e5\tdata again!\n");
    SEQAN_ASSERT_EQ(out, expected);
}

SEQAN_DEFINE_TEST(test_bed_write_bed6_record)
{
    seqan::BedRecord<seqan::Bed6> record1;
    record1.ref = "I";
    record1.beginPos = 123;
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
    String<char> out;
    writeRecord(out, record1, seqan::Bed());
    writeRecord(out, record2, seqan::Bed());

    // Compar string stream to expected value.
    String<char> expected = "I\t123\t456\tNAME1\t5\t-\tsome data that is \tignored\n";
    append(expected, "II\t999\t1000\tNAME2\t3e5\t.\tdata again!\n");
    SEQAN_ASSERT_EQ(out, expected);
}

SEQAN_DEFINE_TEST(test_bed_write_bed12_record)
{
    seqan::BedRecord<seqan::Bed12> record1;
    record1.ref = "I";
    record1.beginPos = 123;
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
    String<char> out;
    writeRecord(out, record1, seqan::Bed());
    writeRecord(out, record2, seqan::Bed());

    // Compar string stream to expected value.
    String<char> expected =
            "I\t123\t456\tNAME1\t5\t-\t123\t234\t10,20,30\t2\t10,20\t3,15\tsome data that is \tignored\n"
            "II\t999\t1000\tNAME2\t3e5\t.\t123\t234\t10,20,30\t2\t10,20\t3,15\tdata again!\n";
    SEQAN_ASSERT_EQ(out, expected);
}

SEQAN_DEFINE_TEST(test_bed_bed_file_read)
{
    seqan::CharString inPath = getAbsolutePath("/tests/bed_io/example.bed");

    seqan::BedFileIn bedStream(toCString(inPath));

    seqan::BedRecord<seqan::Bed3> record1;
    readRecord(record1, bedStream);

    SEQAN_ASSERT_EQ(record1.ref, "chr7");
    SEQAN_ASSERT_EQ(record1.beginPos, 127471196);
    SEQAN_ASSERT_EQ(record1.endPos, 127472363);
    SEQAN_ASSERT_EQ(record1.data, "Pos1\t0\t+\t127471196\t127472363\t255,0,0");

    seqan::BedRecord<seqan::Bed3> record2;
    readRecord(record2, bedStream);
    SEQAN_ASSERT(atEnd(bedStream));

    SEQAN_ASSERT_EQ(record2.ref, "chr8");
    SEQAN_ASSERT_EQ(record2.beginPos, 127472363);
    SEQAN_ASSERT_EQ(record2.endPos, 127473530);
    SEQAN_ASSERT_EQ(record2.data, "Pos2\t0\t+\t127472363\t127473530\t255,0,0");
}

SEQAN_DEFINE_TEST(test_bed_bed_file_write)
{
    seqan::CharString tmpPath = SEQAN_TEMP_FILENAME();
    append(tmpPath, ".bed");

    seqan::BedFileOut bedStream(toCString(tmpPath));

    seqan::BedRecord<seqan::Bed3> record1;
    record1.ref = "chr7";
    record1.beginPos = 127471196;
    record1.endPos = 127472363;
    record1.data = "Pos1\t0\t+\t127471196\t127472363\t255,0,0";
    writeRecord(bedStream, record1);

    seqan::BedRecord<seqan::Bed3> record2;
    record2.ref = "chr8";
    record2.beginPos = 127472363;
    record2.endPos = 127473530;
    record2.data = "Pos2\t0\t+\t127472363\t127473530\t255,0,0";
    writeRecord(bedStream, record2);

    close(bedStream);

    seqan::CharString goldPath(getAbsolutePath("/tests/bed_io/example.bed"));
    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(goldPath)));
}

SEQAN_DEFINE_TEST(test_bed_io_isOpen_fileIn)
{
    // Build path to file.
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/tests/bed_io/example.bed");

    // Create SequenceStream object.
    seqan::BedFileIn bedI;
    SEQAN_ASSERT(!isOpen(bedI));

    // open file
    open(bedI, toCString(filePath));
    SEQAN_ASSERT(isOpen(bedI));

    // close file
    close(bedI);
    SEQAN_ASSERT(!isOpen(bedI));
}

SEQAN_DEFINE_TEST(test_bed_io_isOpen_fileOut)
{
    // Build path to file.
    seqan::CharString filePath = SEQAN_TEMP_FILENAME();
    append(filePath, ".bed");

    // Create SequenceStream object.
    seqan::BedFileOut  bedO;
    SEQAN_ASSERT(!isOpen(bedO));

    // open files
    open(bedO, toCString(filePath));
    SEQAN_ASSERT(isOpen(bedO));

    // close files
    close(bedO);
    SEQAN_ASSERT(!isOpen(bedO));
}

SEQAN_BEGIN_TESTSUITE(test_bed_io)
{
    // Reading of BED records.
    SEQAN_CALL_TEST(test_bed_read_bed3_record);
    SEQAN_CALL_TEST(test_bed_read_bed4_record);
    SEQAN_CALL_TEST(test_bed_read_bed5_record);
    SEQAN_CALL_TEST(test_bed_read_bed6_record);
    SEQAN_CALL_TEST(test_bed_read_bed12_record);

    // Writing of BED records.
    SEQAN_CALL_TEST(test_bed_write_bed3_record);
    SEQAN_CALL_TEST(test_bed_write_bed4_record);
    SEQAN_CALL_TEST(test_bed_write_bed5_record);
    SEQAN_CALL_TEST(test_bed_write_bed6_record);
    SEQAN_CALL_TEST(test_bed_write_bed12_record);

    // BED Stream
    SEQAN_CALL_TEST(test_bed_bed_file_read);
    SEQAN_CALL_TEST(test_bed_bed_file_write);

    // isOpen
    SEQAN_CALL_TEST(test_bed_io_isOpen_fileIn);
    SEQAN_CALL_TEST(test_bed_io_isOpen_fileOut);
}

SEQAN_END_TESTSUITE
