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
#include <seqan/seq_io.h>
#include <seqan/roi_io.h>

SEQAN_DEFINE_TEST(test_roi_read_roi_record)
{
    seqan::String<char> inString = "I\t1\t3\tregion0\t3\t+\t4\t0.55\t33\t1,2,4\n";
    seqan::DirectionIterator<seqan::String<char>, seqan::Input>::Type iter = begin(inString);

    seqan::RoiRecord record;
    seqan::RoiIOContext context;

    readRecord(record, context, iter, seqan::Roi());
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 0);
    SEQAN_ASSERT_EQ(record.endPos, 3);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.len, 3u);
    SEQAN_ASSERT_EQ(record.name, "region0");
    SEQAN_ASSERT_EQ(record.countMax, 4u);
    SEQAN_ASSERT_EQ(length(record.data), 2u);
    SEQAN_ASSERT_EQ(record.data[0], "0.55");
    SEQAN_ASSERT_EQ(record.data[1], "33");
    SEQAN_ASSERT_EQ(length(record.count), 3u);
    SEQAN_ASSERT_EQ(record.count[0], 1u);
    SEQAN_ASSERT_EQ(record.count[1], 2u);
    SEQAN_ASSERT_EQ(record.count[2], 4u);
}

SEQAN_DEFINE_TEST(test_roi_write_roi_record)
{
    seqan::RoiRecord record;
    record.ref = "I";
    record.beginPos = 0;
    record.endPos = 3;
    record.strand = '+';
    record.len = 3;
    record.name = "region0";
    appendValue(record.data, "0.55");
    appendValue(record.data, "33");
    record.countMax = 4;
    appendValue(record.count, 1);
    appendValue(record.count, 2);
    appendValue(record.count, 4);

    seqan::String<char> outString;
    writeRecord(outString, record, seqan::Roi());

    seqan::String<char> expected = "I\t1\t3\tregion0\t3\t+\t4\t0.55\t33\t1,2,4\n";

    SEQAN_ASSERT_EQ(expected, outString);

}

SEQAN_DEFINE_TEST(test_roi_roi_file_read)
{
    seqan::CharString inPath = seqan::getAbsolutePath("/tests/roi_io/example.roi");

    seqan::RoiFileIn roiFileIn(toCString(inPath));

    seqan::RoiRecord record1;
    readRecord(record1, roiFileIn);

    SEQAN_ASSERT_EQ(record1.ref, "I");
    SEQAN_ASSERT_EQ(record1.beginPos, 0);
    SEQAN_ASSERT_EQ(record1.endPos, 3);
    SEQAN_ASSERT_EQ(record1.strand, '+');
    SEQAN_ASSERT_EQ(record1.len, 3u);
    SEQAN_ASSERT_EQ(record1.name, "region0");
    SEQAN_ASSERT_EQ(record1.countMax, 4u);
    SEQAN_ASSERT_EQ(length(record1.count), 3u);
    SEQAN_ASSERT_EQ(record1.count[0], 1u);
    SEQAN_ASSERT_EQ(record1.count[1], 2u);
    SEQAN_ASSERT_EQ(record1.count[2], 4u);

    seqan::RoiRecord record2;
    readRecord(record2, roiFileIn);
    SEQAN_ASSERT(atEnd(roiFileIn));

    SEQAN_ASSERT_EQ(record2.ref, "II");
    SEQAN_ASSERT_EQ(record2.beginPos, 1);
    SEQAN_ASSERT_EQ(record2.endPos, 4);
    SEQAN_ASSERT_EQ(record2.strand, '+');
    SEQAN_ASSERT_EQ(record2.len, 3u);
    SEQAN_ASSERT_EQ(record2.name, "region1");
    SEQAN_ASSERT_EQ(record2.countMax, 10u);
    SEQAN_ASSERT_EQ(length(record2.count), 3u);
    SEQAN_ASSERT_EQ(record2.count[0], 8u);
    SEQAN_ASSERT_EQ(record2.count[1], 9u);
    SEQAN_ASSERT_EQ(record2.count[2], 10u);
}

SEQAN_DEFINE_TEST(test_roi_roi_file_write)
{
    seqan::CharString tmpPath = SEQAN_TEMP_FILENAME();
    append(tmpPath, ".roi");

    seqan::RoiFileOut roiFileOut(toCString(tmpPath));

    seqan::RoiRecord record1;
    record1.ref = "I";
    record1.beginPos = 0;
    record1.endPos = 3;
    record1.strand = '+';
    record1.len = 3;
    record1.name = "region0";
    record1.countMax = 4;
    appendValue(record1.count, 1);
    appendValue(record1.count, 2);
    appendValue(record1.count, 4);
    writeRecord(roiFileOut, record1);

    seqan::RoiRecord record2;
    record2.ref = "II";
    record2.beginPos = 1;
    record2.endPos = 4;
    record2.strand = '+';
    record2.len = 3;
    record2.name = "region1";
    record2.countMax = 10;
    appendValue(record2.count, 8);
    appendValue(record2.count, 9);
    appendValue(record2.count, 10);
    writeRecord(roiFileOut, record2);

    close(roiFileOut);

    seqan::CharString goldPath(seqan::getAbsolutePath("/tests/roi_io/example.roi"));
    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(goldPath)));
}

SEQAN_DEFINE_TEST(test_roi_io_isOpen_fileIn)
{
    // Build path to file.
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/tests/roi_io/example.roi");

    // Create SequenceStream object.
    seqan::RoiFileIn roiI;
    SEQAN_ASSERT(!isOpen(roiI));

    // open file
    open(roiI, toCString(filePath));
    SEQAN_ASSERT(isOpen(roiI));

    // close file
    close(roiI);
    SEQAN_ASSERT(!isOpen(roiI));
}

SEQAN_DEFINE_TEST(test_roi_io_isOpen_fileOut)
{
    // Build path to file.
    seqan::CharString filePath = SEQAN_TEMP_FILENAME();
    append(filePath, ".roi");

    // Create SequenceStream object.
    seqan::RoiFileOut  roiO;
    SEQAN_ASSERT(!isOpen(roiO));

    // open files
    open(roiO, toCString(filePath));
    SEQAN_ASSERT(isOpen(roiO));

    // close files
    close(roiO);
    SEQAN_ASSERT(!isOpen(roiO));
}

SEQAN_BEGIN_TESTSUITE(test_roi_io)
{
    // Reading of ROI records.
    SEQAN_CALL_TEST(test_roi_read_roi_record);

    // Writing of ROI records.
    SEQAN_CALL_TEST(test_roi_write_roi_record);

    // RoiFile
    SEQAN_CALL_TEST(test_roi_roi_file_read);
    SEQAN_CALL_TEST(test_roi_roi_file_write);

    // isOpen
    SEQAN_CALL_TEST(test_roi_io_isOpen_fileIn);
    SEQAN_CALL_TEST(test_roi_io_isOpen_fileOut);
}
SEQAN_END_TESTSUITE
