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

#ifndef CORE_TESTS_BAM_IO_TEST_BAM_INDEX_H_
#define CORE_TESTS_BAM_IO_TEST_BAM_INDEX_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

SEQAN_DEFINE_TEST(test_bam_io_bam_index_bai)
{
    using namespace seqan;

    CharString baiFilename;
    append(baiFilename, SEQAN_PATH_TO_ROOT());
    append(baiFilename, "/core/tests/bam_io/small.bam.bai");

    BamIndex<Bai> baiIndex;
    SEQAN_ASSERT_EQ(read(baiIndex, toCString(baiFilename)), 0);

    SEQAN_ASSERT_EQ(length(baiIndex._binIndices), 1u);
    SEQAN_ASSERT_EQ(baiIndex._binIndices[0].size(), 2u);
    SEQAN_ASSERT(baiIndex._binIndices[0].find(4681) != baiIndex._binIndices[0].end());
    SEQAN_ASSERT(baiIndex._binIndices[0].find(37450) != baiIndex._binIndices[0].end());

    SEQAN_ASSERT_EQ(length(baiIndex._linearIndices), 1u);
    SEQAN_ASSERT_EQ(length(baiIndex._linearIndices[0]), 1u);

    SEQAN_ASSERT_EQ(getUnalignedCount(baiIndex), 0u);

    // File has same contents as in the SAM test.
    CharString bamFilename;
    append(bamFilename, SEQAN_PATH_TO_ROOT());
    append(bamFilename, "/core/tests/bam_io/small.bam");

    Stream<Bgzf> stream;
    open(stream, toCString(bamFilename), "r");

    StringSet<CharString> nameStore;
    NameStoreCache<StringSet<CharString> > nameStoreCache(nameStore);
    BamIOContext<StringSet<CharString> > bamIOContext(nameStore, nameStoreCache);
    
    BamHeader header;
    SEQAN_ASSERT_EQ(readRecord(header, bamIOContext, stream, Bam()), 0);

    bool found = true;
    SEQAN_ASSERT(jumpToRegion(stream, found, bamIOContext, 0, 1, 10, baiIndex));
    SEQAN_ASSERT(found);
    SEQAN_ASSERT(jumpToRegion(stream, found, bamIOContext, 0, 2, 100, baiIndex));
    SEQAN_ASSERT(found);
    SEQAN_ASSERT_NOT(jumpToRegion(stream, found, bamIOContext, 1, 1, 10, baiIndex));
    SEQAN_ASSERT_NOT(found);
}

#endif  // CORE_TESTS_BAM_IO_TEST_BAM_INDEX_H_
