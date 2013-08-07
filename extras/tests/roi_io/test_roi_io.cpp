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
#include <seqan/roi_io.h>

SEQAN_DEFINE_TEST(test_roi_read_roi_record)
{
    std::stringstream ss;
    ss << "I\t1\t3\tregion0\t3\t+\t4\t1,2,4";
    ss.seekg(0);

    seqan::RecordReader<std::stringstream, seqan::SinglePass<> > reader(ss);

    seqan::RoiRecord record;

    SEQAN_ASSERT_EQ(readRecord(record, reader, seqan::Roi()), 0);
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 0);
    SEQAN_ASSERT_EQ(record.endPos, 3);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.len, 3u);
    SEQAN_ASSERT_EQ(record.name, "region0");
    SEQAN_ASSERT_EQ(record.countMax, 4u);
    SEQAN_ASSERT_EQ(length(record.count), 3u);
    SEQAN_ASSERT_EQ(record.count[0], 1u);
    SEQAN_ASSERT_EQ(record.count[1], 2u);
    SEQAN_ASSERT_EQ(record.count[2], 4u);
}

SEQAN_DEFINE_TEST(test_roi_read_roi_record_context)
{
    std::stringstream ss;
    ss << "I\t1\t3\tregion0\t3\t+\t4\t1,2,4";
    ss.seekg(0);

    seqan::RecordReader<std::stringstream, seqan::SinglePass<> > reader(ss);

    typedef seqan::StringSet<seqan::CharString> TNameStore;
    TNameStore refNames;
    seqan::NameStoreCache<TNameStore> refNamesCache(refNames);
    seqan::RoiIOContext<TNameStore> roiIOContext(refNames, refNamesCache);
    
    seqan::RoiRecord record;

    SEQAN_ASSERT_EQ(readRecord(record, reader, roiIOContext, seqan::Roi()), 0);
    SEQAN_ASSERT_EQ(record.ref, "I");
    SEQAN_ASSERT_EQ(record.beginPos, 0);
    SEQAN_ASSERT_EQ(record.endPos, 3);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.len, 3u);
    SEQAN_ASSERT_EQ(record.name, "region0");
    SEQAN_ASSERT_EQ(record.countMax, 4u);
    SEQAN_ASSERT_EQ(length(record.count), 3u);
    SEQAN_ASSERT_EQ(record.count[0], 1u);
    SEQAN_ASSERT_EQ(record.count[1], 2u);
    SEQAN_ASSERT_EQ(record.count[2], 4u);

    SEQAN_ASSERT_EQ(length(refNames), 1u);
    SEQAN_ASSERT_EQ(refNames[0], "I");
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
    record.countMax = 4;
    appendValue(record.count, 1);
    appendValue(record.count, 2);
    appendValue(record.count, 4);

    std::stringstream ss;
    SEQAN_ASSERT_EQ(writeRecord(ss, record, seqan::Roi()), 0);

    char const * EXPECTED = "I\t1\t3\tregion0\t3\t+\t4\t1,2,4\n";
    
    SEQAN_ASSERT_EQ(EXPECTED, ss.str());
}

SEQAN_DEFINE_TEST(test_roi_write_roi_record_context)
{
    seqan::RoiRecord record;
    // NO record.ref = "I", comes from cache!
    record.rID = 0;
    record.beginPos = 0;
    record.endPos = 3;
    record.strand = '+';
    record.len = 3;
    record.name = "region0";
    record.countMax = 4;
    appendValue(record.count, 1);
    appendValue(record.count, 2);
    appendValue(record.count, 4);


    typedef seqan::StringSet<seqan::CharString> TNameStore;
    TNameStore refNames;
    appendValue(refNames, "I");
    seqan::NameStoreCache<TNameStore> refNamesCache(refNames);
    seqan::RoiIOContext<TNameStore> roiIOContext(refNames, refNamesCache);

    std::stringstream ss;
    SEQAN_ASSERT_EQ(writeRecord(ss, record, roiIOContext, seqan::Roi()), 0);

    char const * EXPECTED = "I\t1\t3\tregion0\t3\t+\t4\t1,2,4\n";
    
    SEQAN_ASSERT_EQ(EXPECTED, ss.str());
}

SEQAN_BEGIN_TESTSUITE(test_roi_io)
{
	SEQAN_CALL_TEST(test_roi_read_roi_record);
	SEQAN_CALL_TEST(test_roi_read_roi_record_context);
	SEQAN_CALL_TEST(test_roi_write_roi_record);
	SEQAN_CALL_TEST(test_roi_write_roi_record_context);
}
SEQAN_END_TESTSUITE
