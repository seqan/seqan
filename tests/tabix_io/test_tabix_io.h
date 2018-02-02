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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_TESTS_TABIX_TEST_TABIX_IO_H_
#define SEQAN_TESTS_TABIX_TEST_TABIX_IO_H_

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>
#include <seqan/tabix_io.h>


#if SEQAN_HAS_ZLIB
SEQAN_DEFINE_TEST(test_tabix_io_read_indexed_vcf)
{
    // Open TABIX file
    seqan::CharString vcfPath = seqan::getAbsolutePath("/tests/tabix_io/test.vcf.gz");
    seqan::VcfFileIn vcfFile(toCString(vcfPath));

    // Read header (to get the contig names)
    seqan::VcfHeader header;
    readHeader(header, vcfFile);
    
    // Open Tabix index
    seqan::CharString tbiPath = vcfPath;
    append(tbiPath, ".tbi");
    seqan::TabixIndex tabixIndex(toCString(tbiPath));

    // Search overlapping variants

    // 1st
    bool hasEntries = false;
    SEQAN_ASSERT(jumpToRegion(vcfFile,
                              hasEntries,
                              "chr1",
                              66441,
                              66442,
                              tabixIndex));
    SEQAN_ASSERT(hasEntries);

    seqan::VcfRecord record;
    SEQAN_ASSERT_NOT(atEnd(vcfFile));

    readRecord(record, vcfFile);
    SEQAN_ASSERT_EQ(record.beginPos, 66441);

    readRecord(record, vcfFile);
    SEQAN_ASSERT_EQ(record.beginPos, 66441);

    readRecord(record, vcfFile);
    SEQAN_ASSERT_EQ(record.beginPos, 66479);

    // 2nd
    SEQAN_ASSERT(jumpToRegion(vcfFile,
                              hasEntries,
                              "chr7",
                              62368,
                              62370,
                              tabixIndex));
    SEQAN_ASSERT(hasEntries);

    readRecord(record, vcfFile);
    SEQAN_ASSERT_EQ(record.beginPos, 62369);

    // 3rd - test failures
    SEQAN_ASSERT(jumpToRegion(vcfFile,
                                  hasEntries,
                                  "chr7",
                                  62368,
                                  62369,
                                  tabixIndex));
    SEQAN_ASSERT_NOT(hasEntries);

    SEQAN_ASSERT_NOT(jumpToRegion(vcfFile,
                                  hasEntries,
                                  "chr8",
                                  62368,
                                  62370,
                                  tabixIndex));
    SEQAN_ASSERT_NOT(hasEntries);


    SEQAN_ASSERT_NOT(atEnd(vcfFile));
}

#else // SEQAN_HAS_ZLIB
SEQAN_DEFINE_TEST(test_tabix_io_read_indexed_vcf)
{
    SEQAN_SKIP_TEST;
}
#endif // SEQAN_HAS_ZLIB


#endif  // SEQAN_TESTS_TABIX_TEST_TABIX_IO_H_
