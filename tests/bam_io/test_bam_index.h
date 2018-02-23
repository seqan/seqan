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

#ifndef TESTS_BAM_IO_TEST_BAM_INDEX_H_
#define TESTS_BAM_IO_TEST_BAM_INDEX_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_bam_io_bam_index_build)
{
    CharString expectedBaiFilename = getAbsolutePath("/tests/bam_io/small.bam.bai");

    CharString bamFilename = getAbsolutePath("/tests/bam_io/small.bam");

    CharString tmpOutPath = SEQAN_TEMP_FILENAME();
    append(tmpOutPath, ".bai");

    BamIndex<Bai> baiIndex;
    SEQAN_ASSERT(build(baiIndex, toCString(bamFilename)));
    SEQAN_ASSERT(save(baiIndex, toCString(tmpOutPath)));

    SEQAN_ASSERT(_compareBinaryFiles(toCString(tmpOutPath), toCString(expectedBaiFilename)));
}


SEQAN_DEFINE_TEST(test_bam_io_bam_index_open)
{
    CharString baiFilename = getAbsolutePath("/tests/bam_io/small.bam.bai");

    BamIndex<Bai> baiIndex;
    SEQAN_ASSERT(open(baiIndex, toCString(baiFilename)));

    SEQAN_ASSERT_EQ(length(baiIndex._binIndices), 1u);
    SEQAN_ASSERT_EQ(baiIndex._binIndices[0].size(), 2u);
    SEQAN_ASSERT(baiIndex._binIndices[0].find(4681) != baiIndex._binIndices[0].end());
    SEQAN_ASSERT(baiIndex._binIndices[0].find(37450) != baiIndex._binIndices[0].end());

    SEQAN_ASSERT_EQ(length(baiIndex._linearIndices), 1u);
    SEQAN_ASSERT_EQ(length(baiIndex._linearIndices[0]), 1u);

    SEQAN_ASSERT_EQ(getUnalignedCount(baiIndex), 0u);
}

SEQAN_DEFINE_TEST(test_bam_io_bam_index_save)
{
    CharString baiFilename = getAbsolutePath("/tests/bam_io/small.bam.bai");

    CharString tmpOutPath = SEQAN_TEMP_FILENAME();
    append(tmpOutPath, ".bai");

    BamIndex<Bai> baiIndex;
    SEQAN_ASSERT(open(baiIndex, toCString(baiFilename)));
    SEQAN_ASSERT(save(baiIndex, toCString(tmpOutPath)));

    SEQAN_ASSERT(_compareBinaryFiles(toCString(tmpOutPath), toCString(baiFilename)));
}

SEQAN_DEFINE_TEST(test_bam_io_bam_index_jump_to_region)
{
    CharString bamFileName(getAbsolutePath("/tests/bam_io/ex1.bam"));
    CharString baiFileName(getAbsolutePath("/tests/bam_io/ex1.bam.bai"));

    BamFileIn bamFile;  // Open BamFileIn for reading
    SEQAN_ASSERT(open(bamFile, toCString(bamFileName)));

    BamIndex<Bai> baiFile;  // Open BAI index file
    open(baiFile, toCString(baiFileName));

    BamHeader header;
    readHeader(header, bamFile);

    bool hasAlignments = true;
    BamAlignmentRecord record;

    // 0 is not allowed for region positions
    SEQAN_TEST_EXCEPTION(std::logic_error, jumpToRegion(bamFile, hasAlignments, 0, 0, 12, baiFile));
    SEQAN_ASSERT_NOT_MSG(hasAlignments, "There should be no alignments in region seq1:[1,0].");

    // begin > end
    SEQAN_TEST_EXCEPTION(std::logic_error, jumpToRegion(bamFile, hasAlignments, 0, 2, 1, baiFile));
    SEQAN_ASSERT_NOT_MSG(hasAlignments, "There should be no alignments in region seq1:[1,0].");

    // region is out of bounds
    SEQAN_TEST_EXCEPTION(std::logic_error, jumpToRegion(bamFile, hasAlignments, 0, 2000, 100000, baiFile));
    SEQAN_ASSERT_NOT_MSG(hasAlignments, "No alignments should be in region seq1:[2000,100000].");

    // invalid reference
    SEQAN_TEST_EXCEPTION(std::logic_error, jumpToRegion(bamFile, hasAlignments, 2, 1, 10, baiFile));
    SEQAN_ASSERT_NOT_MSG(hasAlignments, "No alignments should be on contig 2.");

    SEQAN_ASSERT(jumpToRegion(bamFile, hasAlignments, 0, 1, 10, baiFile));
    SEQAN_ASSERT_MSG(hasAlignments, "There should be alignments in region seq1:[1,10].");
    readRecord(record, bamFile);
    SEQAN_ASSERT_EQ_MSG(record.beginPos, 0, "Jumping to region seq1:[1,10] should go to the first read.");

    SEQAN_ASSERT(jumpToRegion(bamFile, hasAlignments, 0, 2, 100, baiFile));
    SEQAN_ASSERT_MSG(hasAlignments, "There should be alignments in region seq1:[2,100].");
    readRecord(record, bamFile);
    SEQAN_ASSERT_EQ_MSG(record.beginPos, 2, "Jumping to region seq1:[2,10] should go to the second read.");

    SEQAN_ASSERT(jumpToRegion(bamFile, hasAlignments, 0, 3, 100, baiFile));
    SEQAN_ASSERT_MSG(hasAlignments, "There should be alignments in region seq1:[3,100].");
    readRecord(record, bamFile);
    SEQAN_ASSERT_EQ_MSG(record.beginPos, 2, "Jumping to region seq1:[3,10] should go to the second read.");

    SEQAN_ASSERT(jumpToRegion(bamFile, hasAlignments, 1, 98, 98, baiFile));
    SEQAN_ASSERT_MSG(hasAlignments, "There should be alignments in region seq2:[99,99].");
    readRecord(record, bamFile);
    SEQAN_ASSERT_EQ(record.beginPos, 97);

    SEQAN_ASSERT(jumpToRegion(bamFile, hasAlignments, 1, 99, 103, baiFile));
    SEQAN_ASSERT_NOT_MSG(hasAlignments, "There should be no alignments in region seq2:[99,103].");

}

SEQAN_DEFINE_TEST(test_bam_io_bam_index_view_records)
{
    CharString bamFileName(getAbsolutePath("/tests/bam_io/ex1.bam"));
    CharString baiFileName(getAbsolutePath("/tests/bam_io/ex1.bam.bai"));

    BamFileIn bamFile;  // Open BamFileIn for reading
    SEQAN_ASSERT(open(bamFile, toCString(bamFileName)));

    BamIndex<Bai> baiFile;  // Open BAI index file
    open(baiFile, toCString(baiFileName));

    BamHeader header;
    readHeader(header, bamFile);

    std::vector<BamAlignmentRecord> records;

    // 0 is not allowed for region positions
    SEQAN_TEST_EXCEPTION(std::logic_error, viewRecords(records, bamFile, baiFile, 0, 0, 12));
    SEQAN_ASSERT_EQ_MSG(length(records), 0u, "No reads expected in region seq1:1-0");

    // begin > end
    SEQAN_TEST_EXCEPTION(std::logic_error, viewRecords(records, bamFile, baiFile, 0, 5, 4));
    SEQAN_ASSERT_EQ_MSG(length(records), 0u, "No reads expected in region seq1:1-0");

    // region is out of bounds
    SEQAN_TEST_EXCEPTION(std::logic_error, viewRecords(records, bamFile, baiFile, 0, 2000, 100000));
    SEQAN_ASSERT_EQ_MSG(length(records), 0u, "No reads expected in region seq1:2000-100000");

    // invalid reference
    SEQAN_TEST_EXCEPTION(std::logic_error, viewRecords(records, bamFile, baiFile, 2, 1, 2000));
    SEQAN_ASSERT_EQ_MSG(length(records), 0u, "No reads expected on invalid contigs.");

    viewRecords(records, bamFile, baiFile, 0, 1, 1);
    SEQAN_ASSERT_EQ_MSG(length(records), 1u, "1 read expected in region seq1:1-1[");
    SEQAN_ASSERT_EQ(front(records).beginPos, 0);
    clear(records);

    viewRecords(records, bamFile, baiFile, 0, 1, 10);
    SEQAN_ASSERT_EQ_MSG(length(records), 5u, "5 reads expected in region seq1:1-10");
    SEQAN_ASSERT_EQ(front(records).beginPos, 0);
    SEQAN_ASSERT_EQ(back(records).beginPos, 8);
    clear(records);

    viewRecords(records, bamFile, baiFile, 0, 1, 100);
    SEQAN_ASSERT_EQ_MSG(length(records), 39u, "39 reads expected in region seq1:1-100");
    SEQAN_ASSERT_EQ(front(records).beginPos, 0);
    SEQAN_ASSERT_EQ(back(records).beginPos, 99);
    clear(records);

    viewRecords(records, bamFile, baiFile, 0, 1, 1575);
    SEQAN_ASSERT_EQ_MSG(length(records), 1501u, "1501 reads expected in region seq1:1-1575");
    SEQAN_ASSERT_EQ(front(records).beginPos, 0);
    SEQAN_ASSERT_EQ(back(records).beginPos, 1534);
    clear(records);

    viewRecords(records, bamFile, baiFile, 0, 100, 1575);
    SEQAN_ASSERT_EQ_MSG(length(records), 1472u, "1472 reads expected in region seq1:100-1575");
    SEQAN_ASSERT_EQ(front(records).beginPos, 66);
    SEQAN_ASSERT_EQ(back(records).beginPos, 1534);
    clear(records);

    // There is an unaligned read at position seq1:100 that should be missed in region [101,1575[
    viewRecords(records, bamFile, baiFile, 0, 101, 1575);
    SEQAN_ASSERT_EQ_MSG(length(records), 1471u, "1471 reads expected in region seq1:101-1575");
    SEQAN_ASSERT_EQ(front(records).beginPos, 66);
    SEQAN_ASSERT_EQ(back(records).beginPos, 1534);
    clear(records);

    viewRecords(records, bamFile, baiFile, 0, 102, 1575);
    SEQAN_ASSERT_EQ_MSG(length(records), 1470u, "1470 reads expected in region seq1:101-1575");
    SEQAN_ASSERT_EQ(front(records).beginPos, 69);
    SEQAN_ASSERT_EQ(back(records).beginPos, 1534);
    clear(records);

    viewRecords(records, bamFile, baiFile, 1, 300, 400);
    SEQAN_ASSERT_EQ_MSG(length(records), 196u, "196 reads expected in region seq2:300-400");
    SEQAN_ASSERT_EQ(front(records).beginPos, 265);
    SEQAN_ASSERT_EQ(back(records).beginPos, 398);
    clear(records);
}

#endif  // TESTS_BAM_IO_TEST_BAM_INDEX_H_
