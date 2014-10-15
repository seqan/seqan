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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the SeqAn module ucsc_io.
// ==========================================================================

#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/ucsc_io.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_ucsc_io_read_record_ucsc_known_genes)
{
    // The file contains 13 annotations in total which will be checked line
    // after line.
    CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/ucsc_io/example_known_genes_with_errors.tsv");

    String<char, MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(ucscPath)));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(mmapString);

    UcscRecord record;
    UcscIOContext ucscIOContext;

    SEQAN_ASSERT_DOES_NOT_THROW(readRecord(record, ucscIOContext, iter, Ucsc()));
    SEQAN_ASSERT_EQ(record.transName, "uc002yoz.1");
    SEQAN_ASSERT_EQ(record.contigName, "chr21");
    SEQAN_ASSERT(record.format == record.KNOWN_GENE);
    SEQAN_ASSERT_EQ(record.cdsBegin, 33026870);
    SEQAN_ASSERT_EQ(record.cdsEnd, 33026870);
    SEQAN_ASSERT_EQ(record.exonBegin[0], 33027740);
    SEQAN_ASSERT_EQ(record.exonBegin[1], 33030540);
    SEQAN_ASSERT_EQ(record.exonBegin[2], 33031813);
    SEQAN_ASSERT_EQ(record.exonEnds[0], 33026870);
    SEQAN_ASSERT_EQ(record.exonEnds[1], 33030246);
    SEQAN_ASSERT_EQ(record.exonEnds[2], 33031709);
    SEQAN_ASSERT_EQ(record.proteinName, "");
    SEQAN_ASSERT_EQ(record.annotationBeginPos, 33031813u);
    SEQAN_ASSERT_EQ(record.annotationEndPos, 33026870u);

    SEQAN_ASSERT_DOES_NOT_THROW(readRecord(record, ucscIOContext, iter, Ucsc()));
    SEQAN_ASSERT_EQ(record.transName, "uc002ypa.3");
    SEQAN_ASSERT_EQ(record.contigName, "chr21");
    SEQAN_ASSERT(record.format == record.KNOWN_GENE);
    SEQAN_ASSERT_EQ(record.cdsBegin, 33032082);
    SEQAN_ASSERT_EQ(record.cdsEnd, 33040891);
    SEQAN_ASSERT_EQ(record.exonBegin[0], 33031934);
    SEQAN_ASSERT_EQ(record.exonBegin[1], 33036102);
    SEQAN_ASSERT_EQ(record.exonBegin[2], 33038761);
    SEQAN_ASSERT_EQ(record.exonBegin[3], 33039570);
    SEQAN_ASSERT_EQ(record.exonBegin[4], 33040783);
    SEQAN_ASSERT_EQ(record.exonEnds[0], 33032154);
    SEQAN_ASSERT_EQ(record.exonEnds[1], 33036199);
    SEQAN_ASSERT_EQ(record.exonEnds[2], 33038831);
    SEQAN_ASSERT_EQ(record.exonEnds[3], 33039688);
    SEQAN_ASSERT_EQ(record.exonEnds[4], 33041243);
    SEQAN_ASSERT_EQ(record.proteinName, "P00441");
    SEQAN_ASSERT_EQ(record.annotationBeginPos, 33031934u);
    SEQAN_ASSERT_EQ(record.annotationEndPos, 33041243u);

    for (unsigned i = 0; i < 20; ++i)
    {
        SEQAN_ASSERT_THROWS(seqan::readRecord(record, ucscIOContext, iter, Ucsc()),
                            ParseError);
        skipLine(iter);
    }
}

SEQAN_DEFINE_TEST(test_ucsc_io_read_record_ucsc_known_isoforms)
{
    // The file contains 13 annotations in total which will be checked line
    // after line.
    CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/ucsc_io/example_known_isoforms.tsv");

    String<char, MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(ucscPath)));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(mmapString);

    UcscRecord record;
    UcscIOContext ucscIOContext;

    SEQAN_ASSERT_DOES_NOT_THROW(readRecord(record, ucscIOContext, iter, Ucsc()));
    SEQAN_ASSERT_EQ(record.transName, "GENE1");
    SEQAN_ASSERT_EQ(record.contigName, "NM_001025288");
    SEQAN_ASSERT(record.format == record.KNOWN_ISOFORMS);

    SEQAN_ASSERT_DOES_NOT_THROW(readRecord(record, ucscIOContext, iter, Ucsc()));
    SEQAN_ASSERT_EQ(record.transName, "GENE2");
    SEQAN_ASSERT_EQ(record.contigName, "NM_134386");
    SEQAN_ASSERT(record.format == record.KNOWN_ISOFORMS);

    SEQAN_ASSERT_DOES_NOT_THROW(readRecord(record, ucscIOContext, iter, Ucsc()));
    SEQAN_ASSERT_EQ(record.transName, "GENE3");
    SEQAN_ASSERT_EQ(record.contigName, "NM_001030033");
    SEQAN_ASSERT(record.format == record.KNOWN_ISOFORMS);

    for (unsigned i = 0; i < 20; ++i)
    {
        SEQAN_ASSERT_THROWS(
            seqan::readRecord(record, ucscIOContext, iter, Ucsc()),
            ParseError);
        skipLine(iter);
    }
}

SEQAN_DEFINE_TEST(test_ucsc_io_write_record_ucsc_known_genes)
{
    CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_known_genes.tsv");

    String<char, MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(ucscPath)));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(mmapString);

    CharString outString;
    while (!atEnd(iter))
    {
        UcscRecord record;
        UcscIOContext ucscIOContext;
        SEQAN_ASSERT_DOES_NOT_THROW(readRecord(record, ucscIOContext, iter, Ucsc()));
        SEQAN_ASSERT_DOES_NOT_THROW(writeRecord(outString, record, Ucsc()));
    }

    SEQAN_ASSERT_EQ(mmapString, outString);
}

SEQAN_DEFINE_TEST(test_ucsc_io_write_record_ucsc_known_isoforms)
{
    CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_known_isoforms.tsv");

    String<char, MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(ucscPath)));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(mmapString);

    UcscIOContext ucscIOContext;
    UcscRecord record;
    CharString outString;
    while (!atEnd(iter))
    {
        SEQAN_ASSERT_DOES_NOT_THROW(readRecord(record, ucscIOContext, iter, Ucsc()));
        SEQAN_ASSERT_DOES_NOT_THROW(writeRecord(outString, record, Ucsc()));
    }

    SEQAN_ASSERT_EQ(mmapString, outString);
}

SEQAN_DEFINE_TEST(test_ucsc_io_ucsc_file_in_read_record_ucsc_known_genes)
{
    // The file contains 13 annotations in total which will be checked line
    // after line.
    CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/ucsc_io/example_known_genes_with_errors.tsv");

    UcscFileIn ucscFileIn(toCString(ucscPath));
    UcscRecord record;

    readRecord(record, ucscFileIn);
    SEQAN_ASSERT_EQ(record.transName, "uc002yoz.1");
    SEQAN_ASSERT_EQ(record.contigName, "chr21");
    SEQAN_ASSERT(record.format == record.KNOWN_GENE);
    SEQAN_ASSERT_EQ(record.cdsBegin, 33026870);
    SEQAN_ASSERT_EQ(record.cdsEnd, 33026870);
    SEQAN_ASSERT_EQ(record.exonBegin[0], 33027740);
    SEQAN_ASSERT_EQ(record.exonBegin[1], 33030540);
    SEQAN_ASSERT_EQ(record.exonBegin[2], 33031813);
    SEQAN_ASSERT_EQ(record.exonEnds[0], 33026870);
    SEQAN_ASSERT_EQ(record.exonEnds[1], 33030246);
    SEQAN_ASSERT_EQ(record.exonEnds[2], 33031709);
    SEQAN_ASSERT_EQ(record.proteinName, "");
    SEQAN_ASSERT_EQ(record.annotationBeginPos, 33031813u);
    SEQAN_ASSERT_EQ(record.annotationEndPos, 33026870u);

    readRecord(record, ucscFileIn);
    SEQAN_ASSERT_EQ(record.transName, "uc002ypa.3");
    SEQAN_ASSERT_EQ(record.contigName, "chr21");
    SEQAN_ASSERT(record.format == record.KNOWN_GENE);
    SEQAN_ASSERT_EQ(record.cdsBegin, 33032082);
    SEQAN_ASSERT_EQ(record.cdsEnd, 33040891);
    SEQAN_ASSERT_EQ(record.exonBegin[0], 33031934);
    SEQAN_ASSERT_EQ(record.exonBegin[1], 33036102);
    SEQAN_ASSERT_EQ(record.exonBegin[2], 33038761);
    SEQAN_ASSERT_EQ(record.exonBegin[3], 33039570);
    SEQAN_ASSERT_EQ(record.exonBegin[4], 33040783);
    SEQAN_ASSERT_EQ(record.exonEnds[0], 33032154);
    SEQAN_ASSERT_EQ(record.exonEnds[1], 33036199);
    SEQAN_ASSERT_EQ(record.exonEnds[2], 33038831);
    SEQAN_ASSERT_EQ(record.exonEnds[3], 33039688);
    SEQAN_ASSERT_EQ(record.exonEnds[4], 33041243);
    SEQAN_ASSERT_EQ(record.proteinName, "P00441");
    SEQAN_ASSERT_EQ(record.annotationBeginPos, 33031934u);
    SEQAN_ASSERT_EQ(record.annotationEndPos, 33041243u);

    for (unsigned i = 0; i < 20; ++i)
    {
        SEQAN_ASSERT_THROWS(seqan::readRecord(record, ucscFileIn), ParseError);
        skipLine(ucscFileIn.iter);
    }
}

SEQAN_DEFINE_TEST(test_ucsc_io_ucsc_file_in_read_record_ucsc_known_isoforms)
{
    // The file contains 13 annotations in total which will be checked line
    // after line.
    CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/ucsc_io/example_known_isoforms.tsv");

    
    UcscFileIn ucscFileIn(toCString(ucscPath));
    UcscRecord record;

    readRecord(record, ucscFileIn);
    SEQAN_ASSERT_EQ(record.transName, "GENE1");
    SEQAN_ASSERT_EQ(record.contigName, "NM_001025288");
    SEQAN_ASSERT(record.format == record.KNOWN_ISOFORMS);

    readRecord(record, ucscFileIn);
    SEQAN_ASSERT_EQ(record.transName, "GENE2");
    SEQAN_ASSERT_EQ(record.contigName, "NM_134386");
    SEQAN_ASSERT(record.format == record.KNOWN_ISOFORMS);

    readRecord(record, ucscFileIn);
    SEQAN_ASSERT_EQ(record.transName, "GENE3");
    SEQAN_ASSERT_EQ(record.contigName, "NM_001030033");
    SEQAN_ASSERT(record.format == record.KNOWN_ISOFORMS);

    for (unsigned i = 0; i < 20; ++i)
    {
        SEQAN_ASSERT_THROWS(seqan::readRecord(record, ucscFileIn), ParseError);
        skipLine(ucscFileIn.iter);
    }
}

SEQAN_DEFINE_TEST(test_ucsc_io_ucsc_file_out_write_record_ucsc_known_genes)
{
    CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_known_genes.tsv");

    String<char, MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(ucscPath)));

    std::stringstream ss;
    UcscFileIn ucscFileIn(toCString(ucscPath));
    UcscFileOut ucscFileOut(ss, Ucsc());

    while (!atEnd(ucscFileIn))
    {
        UcscRecord record;
        SEQAN_ASSERT_DOES_NOT_THROW(readRecord(record, ucscFileIn));
        SEQAN_ASSERT_DOES_NOT_THROW(writeRecord(ucscFileOut, record));
    }

    SEQAN_ASSERT_EQ(mmapString, ss.str().c_str());
}

SEQAN_DEFINE_TEST(test_ucsc_io_ucsc_file_out_write_record_ucsc_known_isoforms)
{
    CharString ucscPath = SEQAN_PATH_TO_ROOT();
    append(ucscPath, "/core/tests/store/example_known_isoforms.tsv");

    String<char, MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(ucscPath)));

    std::stringstream ss;
    UcscFileIn ucscFileIn(toCString(ucscPath));
    UcscFileOut ucscFileOut(ss, Ucsc());

    while (!atEnd(ucscFileIn))
    {
        UcscRecord record;
        SEQAN_ASSERT_DOES_NOT_THROW(readRecord(record, ucscFileIn));
        SEQAN_ASSERT_DOES_NOT_THROW(writeRecord(ucscFileOut, record));
    }

    SEQAN_ASSERT_EQ(mmapString, ss.str().c_str());
}

SEQAN_BEGIN_TESTSUITE(test_ucsc_io)
{
    // Low-level reading of knownGenes and knownIsoforms format.
    SEQAN_CALL_TEST(test_ucsc_io_read_record_ucsc_known_genes);
    SEQAN_CALL_TEST(test_ucsc_io_read_record_ucsc_known_isoforms);

    // Low-level writing of knownGenes and knownIsoforms format.
    SEQAN_CALL_TEST(test_ucsc_io_write_record_ucsc_known_genes);
    SEQAN_CALL_TEST(test_ucsc_io_write_record_ucsc_known_isoforms);

    // Using UcscFileIn for reading of knownGenes and knownIsoforms format.
    SEQAN_CALL_TEST(test_ucsc_io_ucsc_file_in_read_record_ucsc_known_genes);
    SEQAN_CALL_TEST(test_ucsc_io_ucsc_file_in_read_record_ucsc_known_isoforms);

    // Using UcscFileOut for writing of knownGenes and knownIsoforms format.
    SEQAN_CALL_TEST(test_ucsc_io_ucsc_file_out_write_record_ucsc_known_genes);
    SEQAN_CALL_TEST(test_ucsc_io_ucsc_file_out_write_record_ucsc_known_isoforms);
}
SEQAN_END_TESTSUITE
