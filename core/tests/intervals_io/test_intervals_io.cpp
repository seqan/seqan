// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

#include <seqan/basic.h>
#include <seqan/stream.h>

#include <seqan/intervals_io.h>

SEQAN_DEFINE_TEST(test_intervals_io_read_records)
{
    // Input file contents and input iterator.
    std::string input = "chr1:2489281-2489282\n"
                        "chr4:2489841\n";
    typename seqan::DirectionIterator<std::string, seqan::Input>::Type it =
        directionIterator(input, seqan::Input());

    // Buffers and record.
    seqan::CharString buffer;
    seqan::GenomicRegion region;

    // Read and check first record.
    readRecord(region, buffer, it, seqan::Intervals());
    SEQAN_ASSERT_EQ(region.seqName, "chr1");
    SEQAN_ASSERT_EQ(region.beginPos, 2489280);
    SEQAN_ASSERT_EQ(region.endPos, 2489282);

    // Read and check second record.
    readRecord(region, buffer, it, seqan::Intervals());
    SEQAN_ASSERT_EQ(region.seqName, "chr4");
    SEQAN_ASSERT_EQ(region.beginPos, 2489840);
    SEQAN_ASSERT_EQ(region.endPos, 2489841);

    // We should be at the end of the file.
    SEQAN_ASSERT(atEnd(it));
}

SEQAN_DEFINE_TEST(test_intervals_io_write_records)
{
    // Output buffer and output iterator.
    seqan::CharString out;
    typename seqan::DirectionIterator<seqan::CharString, seqan::Output>::Type it =
        directionIterator(out, seqan::Output());

    // Buffer and record.
    seqan::GenomicRegion region;
    seqan::CharString buffer;

    // Write out first record.
    region.seqName = "chr1";
    region.beginPos = 2489280;
    region.endPos = 2489282;
    writeRecord(it, buffer, region, seqan::Intervals());

    // Write out second record.
    region.seqName = "chr4";
    region.beginPos = 2489840;
    region.endPos = 2489841;
    writeRecord(it, buffer, region, seqan::Intervals());

    seqan::CharString expected = "chr1:2489281-2489282\n"
                                 "chr4:2489841\n";
    SEQAN_ASSERT_EQ(expected, out);
}

SEQAN_BEGIN_TESTSUITE(test_intervals_io)
{
	SEQAN_CALL_TEST(test_intervals_io_read_records);
	SEQAN_CALL_TEST(test_intervals_io_write_records);
}
SEQAN_END_TESTSUITE
