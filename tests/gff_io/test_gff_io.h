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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_GFF_IO_TEST_GFF_IO_H_
#define TESTS_GFF_IO_TEST_GFF_IO_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/gff_io.h>

using namespace seqan;


SEQAN_DEFINE_TEST(test_store_io_read_record_context_gff)
{
    CharString gffPath = getAbsolutePath("/tests/gff_io/example.gff");

    String<char, MMap<> > inString;
    open(inString, toCString(gffPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(inString);

    CharString buffer;

    GffRecord record;
    readRecord(record, buffer, iter);

    SEQAN_ASSERT_EQ(record.ref, "ctg123");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "mRNA");
    SEQAN_ASSERT_EQ(record.beginPos, 1299u);
    SEQAN_ASSERT_EQ(record.endPos, 9000u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "ID");
    SEQAN_ASSERT_EQ(record.tagValues[0], "mrna0001");
    SEQAN_ASSERT_EQ(record.tagNames[1], "Name");
    SEQAN_ASSERT_EQ(record.tagValues[1], "sonichedgehog;hehe");

    readRecord(record, buffer, iter);
    SEQAN_ASSERT_EQ(record.ref, "ctg123");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "exon");
    SEQAN_ASSERT_EQ(record.beginPos, 1299u);
    SEQAN_ASSERT_EQ(record.endPos, 1500u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "ID");
    SEQAN_ASSERT_EQ(record.tagValues[0], "exon00001");
    SEQAN_ASSERT_EQ(record.tagNames[1], "Parent");
    SEQAN_ASSERT_EQ(record.tagValues[1], "mrn a0001");

    readRecord(record, buffer, iter);
    SEQAN_ASSERT_EQ(record.ref, "ctg123");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "exon");
    SEQAN_ASSERT_EQ(record.beginPos, 1049u);
    SEQAN_ASSERT_EQ(record.endPos, 1500u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "ID");
    SEQAN_ASSERT_EQ(record.tagValues[0], "exon00002");
    SEQAN_ASSERT_EQ(record.tagNames[1], "Name");
    SEQAN_ASSERT_EQ(record.tagValues[1], "");
    SEQAN_ASSERT_EQ(record.tagNames[2], "Parent");
    SEQAN_ASSERT_EQ(record.tagValues[2], "mrna0001");
}

SEQAN_DEFINE_TEST(test_store_io_write_record_context_gff)
{
    CharString gffPath = getAbsolutePath("/tests/gff_io/example.gff");

    String<char, MMap<> > inString;
    open(inString, toCString(gffPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(inString);

    GffRecord record;
    CharString buffer;

    String<char> compare;

    String<char> temp = "A";
    unsigned count = 0;
    while (!atEnd(iter))
    {
        readRecord(record, buffer, iter);
        record.ref = temp;
        writeRecord(compare, record, Gff());
        appendValue(temp, 'A');
        ++count;
    }

    String<char> goldPath = getAbsolutePath("/tests/gff_io/example_context.gff");
    String<char, MMap<> > goldString;
    open(goldString, toCString(goldPath));

    SEQAN_ASSERT_EQ(goldString, compare);
}

SEQAN_DEFINE_TEST(test_store_io_comment_processing_context_gff)
{
    // slightly modified version of "example.gff". 3 lines of comments are included.
    CharString gffPath = getAbsolutePath("/tests/gff_io/example_with_comments.gff");

    String<char, MMap<> > inString;
    open(inString, toCString(gffPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(inString);

    CharString buffer;

    GffRecord record;
    readRecord(record, buffer, iter);

    SEQAN_ASSERT_EQ(record.ref, "ctg123");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "mRNA");
    SEQAN_ASSERT_EQ(record.beginPos, 1299u);
    SEQAN_ASSERT_EQ(record.endPos, 9000u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "ID");
    SEQAN_ASSERT_EQ(record.tagValues[0], "mrna0001");
    SEQAN_ASSERT_EQ(record.tagNames[1], "Name");
    SEQAN_ASSERT_EQ(record.tagValues[1], "sonichedgehog;hehe");

    readRecord(record, buffer, iter);
    SEQAN_ASSERT_EQ(record.ref, "ctg123");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "exon");
    SEQAN_ASSERT_EQ(record.beginPos, 1299u);
    SEQAN_ASSERT_EQ(record.endPos, 1500u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "ID");
    SEQAN_ASSERT_EQ(record.tagValues[0], "exon00001");
    SEQAN_ASSERT_EQ(record.tagNames[1], "Parent");
    SEQAN_ASSERT_EQ(record.tagValues[1], "mrn a0001");

    readRecord(record, buffer, iter);
    SEQAN_ASSERT_EQ(record.ref, "ctg123");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "exon");
    SEQAN_ASSERT_EQ(record.beginPos, 1049u);
    SEQAN_ASSERT_EQ(record.endPos, 1500u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "ID");
    SEQAN_ASSERT_EQ(record.tagValues[0], "exon00002");
    SEQAN_ASSERT_EQ(record.tagNames[1], "Name");
    SEQAN_ASSERT_EQ(record.tagValues[1], "");
    SEQAN_ASSERT_EQ(record.tagNames[2], "Parent");
}


// Complex GTF format, from pseudogenes.org
SEQAN_DEFINE_TEST(test_store_io_read_record_gtf_pseudogenes)
{
    CharString gtfPath = getAbsolutePath("/tests/gff_io/example_pseudogenes.gtf");

    String<char, MMap<> > inString;
    open(inString, toCString(gtfPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(inString);

    GffRecord record;
    CharString buffer;

    readRecord(record, buffer, iter);
    SEQAN_ASSERT_EQ(record.ref, "chrAL590842");
    SEQAN_ASSERT_EQ(record.source, "pgenes.org");
    SEQAN_ASSERT_EQ(record.type, "pseudogene (p)");
    SEQAN_ASSERT_EQ(record.beginPos, 34872u);
    SEQAN_ASSERT_EQ(record.endPos, 35130u);
    SEQAN_ASSERT_NEQ(record.score, record.score);  // NaN
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(length(record.tagNames), 7u);
    SEQAN_ASSERT_EQ(length(record.tagValues), 7u);
    SEQAN_ASSERT_EQ(record.tagNames[0], "pgene_type");
    SEQAN_ASSERT_EQ(record.tagValues[0], "p");
    SEQAN_ASSERT_EQ(record.tagNames[1], "protein");
    SEQAN_ASSERT_EQ(record.tagValues[1], "Q8ZKU6");
    SEQAN_ASSERT_EQ(record.tagNames[2], "exon_index");
    SEQAN_ASSERT_EQ(record.tagValues[2], "1");
    SEQAN_ASSERT_EQ(record.tagNames[3], "name");
    SEQAN_ASSERT_EQ(record.tagValues[3], "Q8ZKU6.Yersinia_pestis.chrAL590842.mb0");
    SEQAN_ASSERT_EQ(record.tagNames[4], "tax_id");
    SEQAN_ASSERT_EQ(record.tagValues[4], "632");
    SEQAN_ASSERT_EQ(record.tagNames[5], "gene_id");
    SEQAN_ASSERT_EQ(record.tagValues[5], "1");
    SEQAN_ASSERT_EQ(record.tagNames[6], "transcript_id");
    SEQAN_ASSERT_EQ(record.tagValues[6], "urn:lsid:pseudogene.org:632.Pseudogene:1");

    readRecord(record, buffer, iter);
    SEQAN_ASSERT_EQ(record.ref, "chrAL590842");
    SEQAN_ASSERT_EQ(record.source, "pgenes.org");
    SEQAN_ASSERT_EQ(record.type, "pseudogene (p)");
    SEQAN_ASSERT_EQ(record.beginPos, 72639u);
    SEQAN_ASSERT_EQ(record.endPos, 73496u);
    SEQAN_ASSERT_NEQ(record.score, record.score);  // NaN
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(length(record.tagNames), 7u);
    SEQAN_ASSERT_EQ(length(record.tagValues), 7u);
    SEQAN_ASSERT_EQ(record.tagNames[0], "pgene_type");
    SEQAN_ASSERT_EQ(record.tagValues[0], "p");
    SEQAN_ASSERT_EQ(record.tagNames[1], "protein");
    SEQAN_ASSERT_EQ(record.tagValues[1], "Q8ZL54");
    SEQAN_ASSERT_EQ(record.tagNames[2], "exon_index");
    SEQAN_ASSERT_EQ(record.tagValues[2], "1");
    SEQAN_ASSERT_EQ(record.tagNames[3], "name");
    SEQAN_ASSERT_EQ(record.tagValues[3], "Q8ZL54.Yersinia_pestis.chrAL590842.mb0");
    SEQAN_ASSERT_EQ(record.tagNames[4], "tax_id");
    SEQAN_ASSERT_EQ(record.tagValues[4], "632");
    SEQAN_ASSERT_EQ(record.tagNames[5], "gene_id");
    SEQAN_ASSERT_EQ(record.tagValues[5], "2");
    SEQAN_ASSERT_EQ(record.tagNames[6], "transcript_id");
    SEQAN_ASSERT_EQ(record.tagValues[6], "urn:lsid:pseudogene.org:632.Pseudogene:2");

    SEQAN_ASSERT(atEnd(iter));
}

SEQAN_DEFINE_TEST(test_store_io_read_record_context_gtf)
{
    CharString gtfPath = getAbsolutePath("/tests/gff_io/example.gtf");

    String<char, MMap<> > inString;
    open(inString, toCString(gtfPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(inString);

    CharString buffer;

    GffRecord record;
    readRecord(record, buffer, iter);
    SEQAN_ASSERT_EQ(record.ref, "140");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "inter");
    SEQAN_ASSERT_EQ(record.beginPos, 5140u);
    SEQAN_ASSERT_EQ(record.endPos, 8522u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "gene_id");
    SEQAN_ASSERT_EQ(record.tagValues[0], "gene1");
    SEQAN_ASSERT_EQ(record.tagNames[1], "transcript_id");
    SEQAN_ASSERT_EQ(record.tagValues[1], "trans2");
    SEQAN_ASSERT_EQ(record.tagNames[2], "position");
    SEQAN_ASSERT_EQ(record.tagValues[2], "43");

    readRecord(record, buffer, iter);
    SEQAN_ASSERT_EQ(record.ref, "240");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "CDS");
    SEQAN_ASSERT_EQ(record.beginPos, 66995u);
    SEQAN_ASSERT_EQ(record.endPos, 66999u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "gene_id");
    SEQAN_ASSERT_EQ(record.tagValues[0], "140.000");
    SEQAN_ASSERT_EQ(record.tagNames[1], "transcript_id");
    SEQAN_ASSERT_EQ(record.tagValues[1], "140.000.1");

    readRecord(record, buffer, iter);
    SEQAN_ASSERT_EQ(record.ref, "340");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "intron_CNS");
    SEQAN_ASSERT_EQ(record.beginPos, 70102u);
    SEQAN_ASSERT_EQ(record.endPos, 70151u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "gene_id");
    SEQAN_ASSERT_EQ(record.tagValues[0], "140.000");
    SEQAN_ASSERT_EQ(record.tagNames[1], "transcript_id");
    SEQAN_ASSERT_EQ(record.tagValues[1], "140.000.2");
}

SEQAN_DEFINE_TEST(test_store_io_write_record_context_gtf)
{
    CharString gtfPath = getAbsolutePath("/tests/gff_io/example.gtf");

    String<char, MMap<> > inString;
    open(inString, toCString(gtfPath));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(inString);

    GffRecord record;

//    StringSet<CharString> _nameStore;
//    NameStoreCache<StringSet<CharString> > _nameStoreCache(_nameStore);
//    GffIOContext<StringSet<CharString>, NameStoreCache<StringSet<CharString> >, Dependent<> > buffer(_nameStore, _nameStoreCache);

    CharString buffer;
    String<char> outString;

    String<char> temp = "A";
    unsigned count = 0;
    while (!atEnd(iter))
    {
        readRecord(record, buffer, iter);
        record.ref = temp;
        writeRecord(outString, record, Gtf());
        ++count;
        appendValue(temp, 'A');
    }

    String<char> goldPath = getAbsolutePath("/tests/gff_io/example_context.gtf");
    String<char, MMap<> > goldString;
    open(goldString, toCString(goldPath));

    SEQAN_ASSERT_EQ(goldString, outString);
}

SEQAN_DEFINE_TEST(test_store_io_gff_stream_read_record_gff)
{
    CharString gffPath = getAbsolutePath("/tests/gff_io/example.gff");

    GffFileIn gffStream(toCString(gffPath));

    GffRecord record;

    readRecord(record, gffStream);
    SEQAN_ASSERT_EQ(record.ref, "ctg123");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "mRNA");
    SEQAN_ASSERT_EQ(record.beginPos, 1299u);
    SEQAN_ASSERT_EQ(record.endPos, 9000u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "ID");
    SEQAN_ASSERT_EQ(record.tagValues[0], "mrna0001");
    SEQAN_ASSERT_EQ(record.tagNames[1], "Name");
    SEQAN_ASSERT_EQ(record.tagValues[1], "sonichedgehog;hehe");

    readRecord(record, gffStream);
    SEQAN_ASSERT_EQ(record.ref, "ctg123");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "exon");
    SEQAN_ASSERT_EQ(record.beginPos, 1299u);
    SEQAN_ASSERT_EQ(record.endPos, 1500u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "ID");
    SEQAN_ASSERT_EQ(record.tagValues[0], "exon00001");
    SEQAN_ASSERT_EQ(record.tagNames[1], "Parent");
    SEQAN_ASSERT_EQ(record.tagValues[1], "mrn a0001");

    readRecord(record, gffStream);
    SEQAN_ASSERT_EQ(record.ref, "ctg123");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "exon");
    SEQAN_ASSERT_EQ(record.beginPos, 1049u);
    SEQAN_ASSERT_EQ(record.endPos, 1500u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '+');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "ID");
    SEQAN_ASSERT_EQ(record.tagValues[0], "exon00002");
    SEQAN_ASSERT_EQ(record.tagNames[1], "Name");
    SEQAN_ASSERT_EQ(record.tagValues[1], "");
    SEQAN_ASSERT_EQ(record.tagNames[2], "Parent");
    SEQAN_ASSERT_EQ(record.tagValues[2], "mrna0001");
}

SEQAN_DEFINE_TEST(test_store_io_gff_stream_read_record_gtf)
{
    CharString gffPath = getAbsolutePath("/tests/gff_io/example.gtf");

    GffFileIn gffStream(toCString(gffPath));

    GffRecord record;

    readRecord(record, gffStream);
    SEQAN_ASSERT_EQ(record.ref, "140");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "inter");
    SEQAN_ASSERT_EQ(record.beginPos, 5140u);
    SEQAN_ASSERT_EQ(record.endPos, 8522u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "gene_id");
    SEQAN_ASSERT_EQ(record.tagValues[0], "gene1");
    SEQAN_ASSERT_EQ(record.tagNames[1], "transcript_id");
    SEQAN_ASSERT_EQ(record.tagValues[1], "trans2");

    readRecord(record, gffStream);
    SEQAN_ASSERT_EQ(record.ref, "240");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "CDS");
    SEQAN_ASSERT_EQ(record.beginPos, 66995u);
    SEQAN_ASSERT_EQ(record.endPos, 66999u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "gene_id");
    SEQAN_ASSERT_EQ(record.tagValues[0], "140.000");
    SEQAN_ASSERT_EQ(record.tagNames[1], "transcript_id");
    SEQAN_ASSERT_EQ(record.tagValues[1], "140.000.1");

    readRecord(record, gffStream);
    SEQAN_ASSERT_EQ(record.ref, "340");
    SEQAN_ASSERT_EQ(record.source, "");
    SEQAN_ASSERT_EQ(record.type, "intron_CNS");
    SEQAN_ASSERT_EQ(record.beginPos, 70102u);
    SEQAN_ASSERT_EQ(record.endPos, 70151u);
    SEQAN_ASSERT_NEQ(record.score, record.score);
    SEQAN_ASSERT_EQ(record.strand, '-');
    SEQAN_ASSERT_EQ(record.phase, '.');
    SEQAN_ASSERT_EQ(record.tagNames[0], "gene_id");
    SEQAN_ASSERT_EQ(record.tagValues[0], "140.000");
    SEQAN_ASSERT_EQ(record.tagNames[1], "transcript_id");
    SEQAN_ASSERT_EQ(record.tagValues[1], "140.000.2");
}

SEQAN_DEFINE_TEST(test_store_io_gff_stream_write_record_gff)
{
    CharString gtfPath = getAbsolutePath("/tests/gff_io/example.gff");

    GffFileIn inStream(toCString(gtfPath));

    CharString outPath  = SEQAN_TEMP_FILENAME();
    append(outPath, ".gff");

    GffFileOut outStream(toCString(outPath));

    GffRecord record;
    while (!atEnd(inStream))
    {
        readRecord(record, inStream);
        writeRecord(outStream, record);
    }

    close(outStream);

    SEQAN_ASSERT(_compareTextFilesAlt(toCString(outPath), toCString(gtfPath)));
}

SEQAN_DEFINE_TEST(test_store_io_gff_stream_write_record_gtf)
{
    CharString gtfPath = getAbsolutePath("/tests/gff_io/example.gtf");

    GffFileIn inStream(toCString(gtfPath));

    CharString outPath  = SEQAN_TEMP_FILENAME();
    append(outPath, ".gtf");

    GffFileOut outStream(toCString(outPath));

    GffRecord record;
    while (!atEnd(inStream))
    {
        readRecord(record, inStream);
        writeRecord(outStream, record);
    }

    close(outStream);

    SEQAN_ASSERT(_compareTextFilesAlt(toCString(outPath), toCString(gtfPath)));
}

SEQAN_DEFINE_TEST(test_gff_io_isOpen_fileIn)
{
    // Build path to file.
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/tests/gff_io/example.gff");

    // Create SequenceStream object.
    GffFileIn gffI;
    SEQAN_ASSERT(!isOpen(gffI));

    // open file
    open(gffI, toCString(filePath));
    SEQAN_ASSERT(isOpen(gffI));

    // close file
    close(gffI);
    SEQAN_ASSERT(!isOpen(gffI));
}

SEQAN_DEFINE_TEST(test_gff_io_isOpen_fileOut)
{
    // Build path to file.
    seqan::CharString filePath = SEQAN_TEMP_FILENAME();
    append(filePath, ".gff");

    // Create SequenceStream object.
    GffFileOut gffO;
    SEQAN_ASSERT(!isOpen(gffO));

    // open files
    open(gffO, toCString(filePath));
    SEQAN_ASSERT(isOpen(gffO));

    // close files
    close(gffO);
    SEQAN_ASSERT(!isOpen(gffO));
}

#endif  // TESTS_GFF_IO_TEST_GFF_IO_H_
