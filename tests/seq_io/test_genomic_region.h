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

#ifndef TESTS_SEQ_IO_TEST_GENOMIC_REGION_H_
#define TESTS_SEQ_IO_TEST_GENOMIC_REGION_H_

#include <seqan/seq_io.h>

SEQAN_DEFINE_TEST(test_seq_io_genomic_region_default_constructed)
{
    seqan::GenomicRegion region;

    SEQAN_ASSERT(empty(region.seqName));
    SEQAN_ASSERT(region.rID == region.INVALID_ID);
    SEQAN_ASSERT(region.beginPos == region.INVALID_POS);
    SEQAN_ASSERT(region.endPos == region.INVALID_POS);
}

SEQAN_DEFINE_TEST(test_seq_io_genomic_region_from_string)
{
    seqan::GenomicRegion region("chr1:1,000-2,000");

    SEQAN_ASSERT_EQ(region.seqName, "chr1");
    SEQAN_ASSERT(region.rID == region.INVALID_ID);
    SEQAN_ASSERT_EQ(region.beginPos, 999);
    SEQAN_ASSERT_EQ(region.endPos, 2000);
}

SEQAN_DEFINE_TEST(test_seq_io_genomic_region_clear)
{
    seqan::GenomicRegion region;
    region.seqName = "chr1";
    region.rID = 10;
    region.beginPos = 100;
    region.endPos = 1000;

    clear(region);
    SEQAN_ASSERT(empty(region.seqName));
    SEQAN_ASSERT(region.rID == region.INVALID_ID);
    SEQAN_ASSERT(region.beginPos == region.INVALID_POS);
    SEQAN_ASSERT(region.endPos == region.INVALID_POS);
}

SEQAN_DEFINE_TEST(test_seq_io_genomic_region_parse_chrom)
{
    seqan::GenomicRegion region;
    seqan::parse(region, "chr1");

    SEQAN_ASSERT_EQ(region.seqName, "chr1");
    SEQAN_ASSERT(region.rID == region.INVALID_ID);
    SEQAN_ASSERT(region.beginPos == region.INVALID_POS);
    SEQAN_ASSERT(region.endPos == region.INVALID_POS);
}

SEQAN_DEFINE_TEST(test_seq_io_genomic_region_parse_chrom_begin)
{
    seqan::GenomicRegion region;
    seqan::parse(region, "chr1:1,000");

    SEQAN_ASSERT_EQ(region.seqName, "chr1");
    SEQAN_ASSERT(region.rID == region.INVALID_ID);
    SEQAN_ASSERT_EQ(region.beginPos, 999);
    SEQAN_ASSERT(region.endPos == region.INVALID_POS);
}

SEQAN_DEFINE_TEST(test_seq_io_genomic_region_parse_chrom_begin_end)
{
    seqan::GenomicRegion region;
    seqan::parse(region, "chr1:1,000-2,000");

    SEQAN_ASSERT_EQ(region.seqName, "chr1");
    SEQAN_ASSERT(region.rID == region.INVALID_ID);
    SEQAN_ASSERT_EQ(region.beginPos, 999);
    SEQAN_ASSERT_EQ(region.endPos, 2000);
}

SEQAN_DEFINE_TEST(test_seq_io_genomic_region_to_string_interval)
{
    seqan::GenomicRegion region;
    region.seqName = "chr1";
    region.beginPos = 1000;
    region.endPos = 2000;

    seqan::CharString buffer;
    region.toString(buffer);
    SEQAN_ASSERT_EQ(buffer, "chr1:1001-2000");
}

SEQAN_DEFINE_TEST(test_seq_io_genomic_region_to_string_point)
{
    seqan::GenomicRegion region;
    region.seqName = "chr1";
    region.beginPos = 1000;
    region.endPos = 1001;

    seqan::CharString buffer;
    region.toString(buffer);
    SEQAN_ASSERT_EQ(buffer, "chr1:1001");
}

#endif  // #ifndef TESTS_SEQ_IO_TEST_GENOMIC_REGION_H_
