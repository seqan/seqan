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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
// Tests for the SeqAn module ucsc_io.
// ==========================================================================

#define SEQAN_ENABLE_CHECKPOINTS 0

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/ucsc_io.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_store_io_read_record_ucsc_known_genes)
{
    // The file contains 13 annotations in total which will be checked line
    // after line.
    String<char> ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_known_genes.tsv");

    String<char, MMap<> > mmapString;
    open(mmapString, toCString(ucscPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(mmapString);

    UcscRecord record;
    UcscContext ucscContext;

    readRecord(record, iter, ucscContext);
    SEQAN_ASSERT_EQ(record.transName, "uc002yoz.1");
    SEQAN_ASSERT_EQ(record.contigName, "chr21");
    SEQAN_ASSERT(record.format == record.KNOWN_GENE);
    SEQAN_ASSERT_EQ(record.cdsBegin, 33026870u);
    SEQAN_ASSERT_EQ(record.cdsEnd, 33026870u);
    SEQAN_ASSERT_EQ(record.exonBegin[0], 33027740u);
    SEQAN_ASSERT_EQ(record.exonBegin[1], 33030540u);
    SEQAN_ASSERT_EQ(record.exonBegin[2], 33031813u);
    SEQAN_ASSERT_EQ(record.exonEnd[0], 33026870u);
    SEQAN_ASSERT_EQ(record.exonEnd[1], 33030246u);
    SEQAN_ASSERT_EQ(record.exonEnd[2], 33031709u);
    SEQAN_ASSERT_EQ(record.proteinName, "");
    SEQAN_ASSERT_EQ(record.annotationBeginPos, 33031813u);
    SEQAN_ASSERT_EQ(record.annotationEndPos, 33026870u);

    readRecord(record, iter, ucscContext);
    SEQAN_ASSERT_EQ(record.transName, "uc002ypa.3");
    SEQAN_ASSERT_EQ(record.contigName, "chr21");
    SEQAN_ASSERT(record.format == record.KNOWN_GENE);
    SEQAN_ASSERT_EQ(record.cdsBegin, 33032082u);
    SEQAN_ASSERT_EQ(record.cdsEnd, 33040891u);
    SEQAN_ASSERT_EQ(record.exonBegin[0], 33031934u);
    SEQAN_ASSERT_EQ(record.exonBegin[1], 33036102u);
    SEQAN_ASSERT_EQ(record.exonBegin[2], 33038761u);
    SEQAN_ASSERT_EQ(record.exonBegin[3], 33039570u);
    SEQAN_ASSERT_EQ(record.exonBegin[4], 33040783u);
    SEQAN_ASSERT_EQ(record.exonEnd[0], 33032154u);
    SEQAN_ASSERT_EQ(record.exonEnd[1], 33036199u);
    SEQAN_ASSERT_EQ(record.exonEnd[2], 33038831u);
    SEQAN_ASSERT_EQ(record.exonEnd[3], 33039688u);
    SEQAN_ASSERT_EQ(record.exonEnd[4], 33041243u);
    SEQAN_ASSERT_EQ(record.proteinName, "P00441");
    SEQAN_ASSERT_EQ(record.annotationBeginPos, 33031934u);
    SEQAN_ASSERT_EQ(record.annotationEndPos, 33041243u);

    for (unsigned i = 0; i < 20; ++i)
    {
        SEQAN_TEST_EXCEPTION(ParseError,
                             seqan::readRecord(record, iter, ucscContext));
        skipLine(iter);
    }

}

SEQAN_DEFINE_TEST(test_store_io_read_ucsc_known_genes)
{
    // The file contains 13 annotations in total which will be checked line
    // after line.
    String<char> ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/ucsc_io/example_known_genes.tsv");

    String<char, MMap<> > mmapString;
    open(mmapString, toCString(ucscPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(mmapString);

    String<UcscRecord> records;
    read(records, iter);

    SEQAN_ASSERT_EQ(records[0].transName, "uc002yoz.1");
    SEQAN_ASSERT_EQ(records[0].contigName, "chr21");
    SEQAN_ASSERT(records[0].format == records[0].KNOWN_GENE);
    SEQAN_ASSERT_EQ(records[0].cdsBegin, 33026870u);
    SEQAN_ASSERT_EQ(records[0].cdsEnd, 33026870u);
    SEQAN_ASSERT_EQ(records[0].exonBegin[0], 33027740u);
    SEQAN_ASSERT_EQ(records[0].exonBegin[1], 33030540u);
    SEQAN_ASSERT_EQ(records[0].exonBegin[2], 33031813u);
    SEQAN_ASSERT_EQ(records[0].exonEnd[0], 33026870u);
    SEQAN_ASSERT_EQ(records[0].exonEnd[1], 33030246u);
    SEQAN_ASSERT_EQ(records[0].exonEnd[2], 33031709u);
    SEQAN_ASSERT_EQ(records[0].proteinName, "");
    SEQAN_ASSERT_EQ(records[0].annotationBeginPos, 33031813u);
    SEQAN_ASSERT_EQ(records[0].annotationEndPos, 33026870u);

    SEQAN_ASSERT_EQ(records[1].transName, "uc002ypa.3");
    SEQAN_ASSERT_EQ(records[1].contigName, "chr21");
    SEQAN_ASSERT(records[1].format == records[1].KNOWN_GENE);
    SEQAN_ASSERT_EQ(records[1].cdsBegin, 33032082u);
    SEQAN_ASSERT_EQ(records[1].cdsEnd, 33040891u);
    SEQAN_ASSERT_EQ(records[1].exonBegin[0], 33031934u);
    SEQAN_ASSERT_EQ(records[1].exonBegin[1], 33036102u);
    SEQAN_ASSERT_EQ(records[1].exonBegin[2], 33038761u);
    SEQAN_ASSERT_EQ(records[1].exonBegin[3], 33039570u);
    SEQAN_ASSERT_EQ(records[1].exonBegin[4], 33040783u);
    SEQAN_ASSERT_EQ(records[1].exonEnd[0], 33032154u);
    SEQAN_ASSERT_EQ(records[1].exonEnd[1], 33036199u);
    SEQAN_ASSERT_EQ(records[1].exonEnd[2], 33038831u);
    SEQAN_ASSERT_EQ(records[1].exonEnd[3], 33039688u);
    SEQAN_ASSERT_EQ(records[1].exonEnd[4], 33041243u);
    SEQAN_ASSERT_EQ(records[1].proteinName, "P00441");
    SEQAN_ASSERT_EQ(records[1].annotationBeginPos, 33031934u);
    SEQAN_ASSERT_EQ(records[1].annotationEndPos, 33041243u);
}

SEQAN_DEFINE_TEST(test_store_io_write_record_ucsc_known_genes)
{
    String<char> ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_known_genes.tsv");

    String<char, MMap<> > mmapString;
    open(mmapString, toCString(ucscPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(mmapString);

    String<char> outString;
    while (!atEnd(iter))
    {
        UcscRecord record;
        UcscContext ucscContext;
        readRecord(record, iter, ucscContext);
        writeRecord(outString, record);
    }

    SEQAN_ASSERT_EQ(mmapString, outString);
}

SEQAN_DEFINE_TEST(test_store_io_write_ucsc_known_genes)
{
    String<char> ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_known_genes.tsv");

    String<char, MMap<> > mmapString;
    open(mmapString, toCString(ucscPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(mmapString);

    String<UcscRecord> records;
    read(records, iter);

    String<char> outString;
    write(outString, records);

    SEQAN_ASSERT_EQ(mmapString, outString);
}
