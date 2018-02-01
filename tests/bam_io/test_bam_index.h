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

    // File has same contents as in the SAM test.
    CharString bamFilename = getAbsolutePath("/tests/bam_io/small.bam");

    BamFileIn bamFile(toCString(bamFilename));

    BamHeader header;
    readHeader(header, bamFile);

    bool found = true;
    SEQAN_ASSERT(jumpToRegion(bamFile, found, 0, 1, 10, baiIndex));
    SEQAN_ASSERT(found);
    SEQAN_ASSERT(jumpToRegion(bamFile, found, 0, 2, 100, baiIndex));
    SEQAN_ASSERT(found);
    SEQAN_ASSERT_NOT(jumpToRegion(bamFile, found, 1, 1, 10, baiIndex));
    SEQAN_ASSERT_NOT(found);
    SEQAN_ASSERT(jumpToRegion(bamFile, found, 0, 20, 100, baiIndex));
    SEQAN_ASSERT_NOT(found);
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


#endif  // TESTS_BAM_IO_TEST_BAM_INDEX_H_
