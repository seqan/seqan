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

#ifndef SEQAN_EXTRAS_TESTS_VCF_TEST_VCF_IO_H_
#define SEQAN_EXTRAS_TESTS_VCF_TEST_VCF_IO_H_

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

SEQAN_DEFINE_TEST(test_vcf_io_read_vcf_header)
{
    seqan::CharString vcfPath = SEQAN_PATH_TO_ROOT();
    append(vcfPath, "/extras/tests/vcf_io/example.vcf");

    std::fstream inF(toCString(vcfPath), std::ios::in | std::ios::binary);
    SEQAN_ASSERT(inF.good());

    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(inF);
    seqan::VcfHeader vcfHeader;
    seqan::VcfIOContext vcfIOContext(vcfHeader.sequenceNames, vcfHeader.sampleNames);

    SEQAN_ASSERT_EQ(read(vcfHeader, reader, vcfIOContext, seqan::Vcf()), 0);

    SEQAN_ASSERT_EQ(length(vcfHeader.headerRecords), 18u);
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[0].key, "fileformat");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[0].value, "VCFv4.1");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[1].key, "fileDate");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[1].value, "20090805");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[2].key, "source");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[2].value, "myImputationProgramV3.1");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[3].key, "reference");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[3].value, "file:///seq/references/1000GenomesPilot-NCBI36.fasta");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[4].key, "contig");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[4].value, "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[5].key, "phasing");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[5].value, "partial");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[6].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[6].value, "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[7].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[7].value, "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[8].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[8].value, "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[9].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[9].value, "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[10].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[10].value, "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[11].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[11].value, "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[12].key, "FILTER");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[12].value, "<ID=q10,Description=\"Quality below 10\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[13].key, "FILTER");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[13].value, "<ID=s50,Description=\"Less than 50% of samples have data\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[14].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[14].value, "<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[15].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[15].value, "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[16].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[16].value, "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[17].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfHeader.headerRecords[17].value, "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");

    SEQAN_ASSERT_EQ(length(vcfHeader.sequenceNames), 1u);
    SEQAN_ASSERT_EQ(vcfHeader.sequenceNames[0], "20");

    SEQAN_ASSERT_EQ(length(vcfHeader.sampleNames), 3u);
    SEQAN_ASSERT_EQ(vcfHeader.sampleNames[0], "NA00001");
    SEQAN_ASSERT_EQ(vcfHeader.sampleNames[1], "NA00002");
    SEQAN_ASSERT_EQ(vcfHeader.sampleNames[2], "NA00003");
}

SEQAN_DEFINE_TEST(test_vcf_io_read_vcf_record)
{
    seqan::CharString vcfPath = SEQAN_PATH_TO_ROOT();
    append(vcfPath, "/extras/tests/vcf_io/example.vcf");

    std::fstream inF(toCString(vcfPath), std::ios::in | std::ios::binary);
    SEQAN_ASSERT(inF.good());

    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(inF);
    seqan::VcfHeader vcfHeader;
    seqan::VcfIOContext vcfIOContext(vcfHeader.sequenceNames, vcfHeader.sampleNames);

    SEQAN_ASSERT_EQ(read(vcfHeader, reader, vcfIOContext, seqan::Vcf()), 0);

    seqan::String<seqan::VcfRecord> records;
    while (!atEnd(reader))
    {
        seqan::VcfRecord record;
        SEQAN_ASSERT_EQ(readRecord(record, reader, vcfIOContext, seqan::Vcf()), 0);
        appendValue(records, record);
    }

    SEQAN_ASSERT_EQ(length(records), 3u);

    SEQAN_ASSERT_EQ(records[0].rID, 0);
    SEQAN_ASSERT_EQ(records[0].beginPos, 14369);
    SEQAN_ASSERT_EQ(records[0].id, "rs6054257");
    SEQAN_ASSERT_EQ(records[0].ref, "G");
    SEQAN_ASSERT_EQ(records[0].alt, "A");
    SEQAN_ASSERT_EQ(records[0].qual, 29);
    SEQAN_ASSERT_EQ(records[0].filter, "PASS");
    SEQAN_ASSERT_EQ(records[0].info, "NS=3;DP=14;AF=0.5;DB;H2");
    SEQAN_ASSERT_EQ(records[0].format, "GT:GQ:DP:HQ");
    SEQAN_ASSERT_EQ(length(records[0].genotypeInfos), 3u);

    SEQAN_ASSERT_EQ(records[1].rID, 0);
    SEQAN_ASSERT_EQ(records[1].beginPos, 17329);
    SEQAN_ASSERT_EQ(records[1].id, ".");
    SEQAN_ASSERT_EQ(records[1].ref, "T");
    SEQAN_ASSERT_EQ(records[1].alt, "A");
    SEQAN_ASSERT_EQ(records[1].qual, 3);
    SEQAN_ASSERT_EQ(records[1].filter, "q10");
    SEQAN_ASSERT_EQ(records[1].info, "NS=3;DP=11;AF=0.017");
    SEQAN_ASSERT_EQ(records[1].format, "GT:GQ:DP:HQ");
    SEQAN_ASSERT_EQ(length(records[1].genotypeInfos), 3u);

    SEQAN_ASSERT_EQ(records[2].rID, 0);
    SEQAN_ASSERT_EQ(records[2].beginPos, 1110695);
    SEQAN_ASSERT_EQ(records[2].id, "rs6040355");
    SEQAN_ASSERT_EQ(records[2].ref, "A");
    SEQAN_ASSERT_EQ(records[2].alt, "G,T");
    SEQAN_ASSERT_EQ(records[2].qual, 67);
    SEQAN_ASSERT_EQ(records[2].filter, "PASS");
    SEQAN_ASSERT_EQ(records[2].info, "NS=2;DP=10;AF=0.333,0.667;AA=T;DB");
    SEQAN_ASSERT_EQ(records[2].format, "GT:GQ:DP:HQ");
    SEQAN_ASSERT_EQ(length(records[2].genotypeInfos), 3u);
}

SEQAN_DEFINE_TEST(test_vcf_io_vcf_stream_read_record)
{
    seqan::CharString vcfPath = SEQAN_PATH_TO_ROOT();
    append(vcfPath, "/extras/tests/vcf_io/example.vcf");

    seqan::VcfStream vcfStream(toCString(vcfPath));
    SEQAN_ASSERT(isGood(vcfStream));

    SEQAN_ASSERT_EQ(length(vcfStream.header.headerRecords), 18u);
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[0].key, "fileformat");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[0].value, "VCFv4.1");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[1].key, "fileDate");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[1].value, "20090805");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[2].key, "source");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[2].value, "myImputationProgramV3.1");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[3].key, "reference");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[3].value, "file:///seq/references/1000GenomesPilot-NCBI36.fasta");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[4].key, "contig");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[4].value, "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[5].key, "phasing");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[5].value, "partial");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[6].key, "INFO");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[6].value, "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[7].key, "INFO");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[7].value, "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[8].key, "INFO");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[8].value, "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[9].key, "INFO");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[9].value, "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[10].key, "INFO");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[10].value, "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[11].key, "INFO");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[11].value, "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[12].key, "FILTER");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[12].value, "<ID=q10,Description=\"Quality below 10\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[13].key, "FILTER");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[13].value, "<ID=s50,Description=\"Less than 50% of samples have data\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[14].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[14].value, "<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[15].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[15].value, "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[16].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[16].value, "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[17].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfStream.header.headerRecords[17].value, "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");

    SEQAN_ASSERT_EQ(length(vcfStream.header.sequenceNames), 1u);
    SEQAN_ASSERT_EQ(vcfStream.header.sequenceNames[0], "20");

    SEQAN_ASSERT_EQ(length(vcfStream.header.sampleNames), 3u);
    SEQAN_ASSERT_EQ(vcfStream.header.sampleNames[0], "NA00001");
    SEQAN_ASSERT_EQ(vcfStream.header.sampleNames[1], "NA00002");
    SEQAN_ASSERT_EQ(vcfStream.header.sampleNames[2], "NA00003");

    seqan::String<seqan::VcfRecord> records;
    while (!atEnd(vcfStream))
    {
        seqan::VcfRecord record;
        SEQAN_ASSERT_EQ(readRecord(record, vcfStream), 0);
        appendValue(records, record);
    }

    SEQAN_ASSERT_EQ(length(records), 3u);

    SEQAN_ASSERT_EQ(records[0].rID, 0);
    SEQAN_ASSERT_EQ(records[0].beginPos, 14369);
    SEQAN_ASSERT_EQ(records[0].id, "rs6054257");
    SEQAN_ASSERT_EQ(records[0].ref, "G");
    SEQAN_ASSERT_EQ(records[0].alt, "A");
    SEQAN_ASSERT_EQ(records[0].qual, 29);
    SEQAN_ASSERT_EQ(records[0].filter, "PASS");
    SEQAN_ASSERT_EQ(records[0].info, "NS=3;DP=14;AF=0.5;DB;H2");
    SEQAN_ASSERT_EQ(records[0].format, "GT:GQ:DP:HQ");
    SEQAN_ASSERT_EQ(length(records[0].genotypeInfos), 3u);

    SEQAN_ASSERT_EQ(records[1].rID, 0);
    SEQAN_ASSERT_EQ(records[1].beginPos, 17329);
    SEQAN_ASSERT_EQ(records[1].id, ".");
    SEQAN_ASSERT_EQ(records[1].ref, "T");
    SEQAN_ASSERT_EQ(records[1].alt, "A");
    SEQAN_ASSERT_EQ(records[1].qual, 3);
    SEQAN_ASSERT_EQ(records[1].filter, "q10");
    SEQAN_ASSERT_EQ(records[1].info, "NS=3;DP=11;AF=0.017");
    SEQAN_ASSERT_EQ(records[1].format, "GT:GQ:DP:HQ");
    SEQAN_ASSERT_EQ(length(records[1].genotypeInfos), 3u);

    SEQAN_ASSERT_EQ(records[2].rID, 0);
    SEQAN_ASSERT_EQ(records[2].beginPos, 1110695);
    SEQAN_ASSERT_EQ(records[2].id, "rs6040355");
    SEQAN_ASSERT_EQ(records[2].ref, "A");
    SEQAN_ASSERT_EQ(records[2].alt, "G,T");
    SEQAN_ASSERT_EQ(records[2].qual, 67);
    SEQAN_ASSERT_EQ(records[2].filter, "PASS");
    SEQAN_ASSERT_EQ(records[2].info, "NS=2;DP=10;AF=0.333,0.667;AA=T;DB");
    SEQAN_ASSERT_EQ(records[2].format, "GT:GQ:DP:HQ");
    SEQAN_ASSERT_EQ(length(records[2].genotypeInfos), 3u);
}

SEQAN_DEFINE_TEST(test_vcf_io_write_vcf_header)
{
    seqan::CharString tmpPath(SEQAN_TEMP_FILENAME());
    std::fstream outF(toCString(tmpPath), std::ios::out | std::ios::binary);
    SEQAN_ASSERT(outF.good());

    seqan::VcfHeader vcfHeader;
    appendValue(vcfHeader.sequenceNames, "20");
    appendValue(vcfHeader.sampleNames, "NA00001");
    appendValue(vcfHeader.sampleNames, "NA00002");
    appendValue(vcfHeader.sampleNames, "NA00003");

    resize(vcfHeader.headerRecords, 18);
    vcfHeader.headerRecords[0].key = "fileformat";
    vcfHeader.headerRecords[0].value = "VCFv4.1";
    vcfHeader.headerRecords[1].key = "fileDate";
    vcfHeader.headerRecords[1].value = "20090805";
    vcfHeader.headerRecords[2].key = "source";
    vcfHeader.headerRecords[2].value = "myImputationProgramV3.1";
    vcfHeader.headerRecords[3].key = "reference";
    vcfHeader.headerRecords[3].value = "file:///seq/references/1000GenomesPilot-NCBI36.fasta";
    vcfHeader.headerRecords[4].key = "contig";
    vcfHeader.headerRecords[4].value = "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>";
    vcfHeader.headerRecords[5].key = "phasing";
    vcfHeader.headerRecords[5].value = "partial";
    vcfHeader.headerRecords[6].key = "INFO";
    vcfHeader.headerRecords[6].value = "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">";
    vcfHeader.headerRecords[7].key = "INFO";
    vcfHeader.headerRecords[7].value = "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
    vcfHeader.headerRecords[8].key = "INFO";
    vcfHeader.headerRecords[8].value = "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">";
    vcfHeader.headerRecords[9].key = "INFO";
    vcfHeader.headerRecords[9].value = "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">";
    vcfHeader.headerRecords[10].key = "INFO";
    vcfHeader.headerRecords[10].value = "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">";
    vcfHeader.headerRecords[11].key = "INFO";
    vcfHeader.headerRecords[11].value = "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">";
    vcfHeader.headerRecords[12].key = "FILTER";
    vcfHeader.headerRecords[12].value = "<ID=q10,Description=\"Quality below 10\">";
    vcfHeader.headerRecords[13].key = "FILTER";
    vcfHeader.headerRecords[13].value = "<ID=s50,Description=\"Less than 50% of samples have data\">";
    vcfHeader.headerRecords[14].key = "FORMAT";
    vcfHeader.headerRecords[14].value = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    vcfHeader.headerRecords[15].key = "FORMAT";
    vcfHeader.headerRecords[15].value = "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
    vcfHeader.headerRecords[16].key = "FORMAT";
    vcfHeader.headerRecords[16].value = "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">";
    vcfHeader.headerRecords[17].key = "FORMAT";
    vcfHeader.headerRecords[17].value = "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">";

    seqan::VcfIOContext vcfIOContext(vcfHeader.sequenceNames, vcfHeader.sampleNames);
    SEQAN_ASSERT_EQ(write(outF, vcfHeader, vcfIOContext, seqan::Vcf()), 0);
    outF.close();

    seqan::CharString goldPath(SEQAN_PATH_TO_ROOT());
    append(goldPath, "/extras/tests/vcf_io/vcf_header.vcf");
    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(goldPath)));
}

SEQAN_DEFINE_TEST(test_vcf_io_write_vcf_record)
{
    seqan::CharString tmpPath(SEQAN_TEMP_FILENAME());
    std::fstream outF(toCString(tmpPath), std::ios::out | std::ios::binary);
    SEQAN_ASSERT(outF.good());

    seqan::VcfHeader vcfHeader;
    appendValue(vcfHeader.sequenceNames, "20");
    appendValue(vcfHeader.sampleNames, "NA00001");
    appendValue(vcfHeader.sampleNames, "NA00002");
    appendValue(vcfHeader.sampleNames, "NA00003");

    resize(vcfHeader.headerRecords, 18);
    vcfHeader.headerRecords[0].key = "fileformat";
    vcfHeader.headerRecords[0].value = "VCFv4.1";
    vcfHeader.headerRecords[1].key = "fileDate";
    vcfHeader.headerRecords[1].value = "20090805";
    vcfHeader.headerRecords[2].key = "source";
    vcfHeader.headerRecords[2].value = "myImputationProgramV3.1";
    vcfHeader.headerRecords[3].key = "reference";
    vcfHeader.headerRecords[3].value = "file:///seq/references/1000GenomesPilot-NCBI36.fasta";
    vcfHeader.headerRecords[4].key = "contig";
    vcfHeader.headerRecords[4].value = "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>";
    vcfHeader.headerRecords[5].key = "phasing";
    vcfHeader.headerRecords[5].value = "partial";
    vcfHeader.headerRecords[6].key = "INFO";
    vcfHeader.headerRecords[6].value = "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">";
    vcfHeader.headerRecords[7].key = "INFO";
    vcfHeader.headerRecords[7].value = "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
    vcfHeader.headerRecords[8].key = "INFO";
    vcfHeader.headerRecords[8].value = "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">";
    vcfHeader.headerRecords[9].key = "INFO";
    vcfHeader.headerRecords[9].value = "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">";
    vcfHeader.headerRecords[10].key = "INFO";
    vcfHeader.headerRecords[10].value = "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">";
    vcfHeader.headerRecords[11].key = "INFO";
    vcfHeader.headerRecords[11].value = "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">";
    vcfHeader.headerRecords[12].key = "FILTER";
    vcfHeader.headerRecords[12].value = "<ID=q10,Description=\"Quality below 10\">";
    vcfHeader.headerRecords[13].key = "FILTER";
    vcfHeader.headerRecords[13].value = "<ID=s50,Description=\"Less than 50% of samples have data\">";
    vcfHeader.headerRecords[14].key = "FORMAT";
    vcfHeader.headerRecords[14].value = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    vcfHeader.headerRecords[15].key = "FORMAT";
    vcfHeader.headerRecords[15].value = "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
    vcfHeader.headerRecords[16].key = "FORMAT";
    vcfHeader.headerRecords[16].value = "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">";
    vcfHeader.headerRecords[17].key = "FORMAT";
    vcfHeader.headerRecords[17].value = "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">";

    seqan::VcfIOContext vcfIOContext(vcfHeader.sequenceNames, vcfHeader.sampleNames);

    seqan::VcfRecord vcfRecord;
    vcfRecord.rID = 0;
    vcfRecord.beginPos = 14369;
    vcfRecord.id = "rs6054257";
    vcfRecord.ref = "G";
    vcfRecord.alt = "A";
    vcfRecord.qual = 29;
    vcfRecord.filter = "PASS";
    vcfRecord.info = "NS=3;DP=14;AF=0.5;DB;H2";
    vcfRecord.format = "GT:GQ:DP:HQ";
    appendValue(vcfRecord.genotypeInfos, "0|0:48:1:51,51");
    appendValue(vcfRecord.genotypeInfos, "1|0:48:8:51,51");
    appendValue(vcfRecord.genotypeInfos, "1/1:43:5:.,.");
    SEQAN_ASSERT_EQ(writeRecord(outF, vcfRecord, vcfIOContext, seqan::Vcf()), 0);
    outF.close();

    seqan::CharString goldPath(SEQAN_PATH_TO_ROOT());
    append(goldPath, "/extras/tests/vcf_io/vcf_record.vcf");
    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(goldPath)));
}

SEQAN_DEFINE_TEST(test_vcf_io_vcf_stream_write_record)
{
    seqan::CharString tmpPath(SEQAN_TEMP_FILENAME());
    seqan::VcfStream vcfStream(toCString(tmpPath), seqan::VcfStream::WRITE);
    SEQAN_ASSERT(isGood(vcfStream));

    // Build header.
    appendValue(vcfStream.header.sequenceNames, "20");
    appendValue(vcfStream.header.sampleNames, "NA00001");
    appendValue(vcfStream.header.sampleNames, "NA00002");
    appendValue(vcfStream.header.sampleNames, "NA00003");

    resize(vcfStream.header.headerRecords, 18);
    vcfStream.header.headerRecords[0].key = "fileformat";
    vcfStream.header.headerRecords[0].value = "VCFv4.1";
    vcfStream.header.headerRecords[1].key = "fileDate";
    vcfStream.header.headerRecords[1].value = "20090805";
    vcfStream.header.headerRecords[2].key = "source";
    vcfStream.header.headerRecords[2].value = "myImputationProgramV3.1";
    vcfStream.header.headerRecords[3].key = "reference";
    vcfStream.header.headerRecords[3].value = "file:///seq/references/1000GenomesPilot-NCBI36.fasta";
    vcfStream.header.headerRecords[4].key = "contig";
    vcfStream.header.headerRecords[4].value = "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>";
    vcfStream.header.headerRecords[5].key = "phasing";
    vcfStream.header.headerRecords[5].value = "partial";
    vcfStream.header.headerRecords[6].key = "INFO";
    vcfStream.header.headerRecords[6].value = "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">";
    vcfStream.header.headerRecords[7].key = "INFO";
    vcfStream.header.headerRecords[7].value = "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
    vcfStream.header.headerRecords[8].key = "INFO";
    vcfStream.header.headerRecords[8].value = "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">";
    vcfStream.header.headerRecords[9].key = "INFO";
    vcfStream.header.headerRecords[9].value = "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">";
    vcfStream.header.headerRecords[10].key = "INFO";
    vcfStream.header.headerRecords[10].value = "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">";
    vcfStream.header.headerRecords[11].key = "INFO";
    vcfStream.header.headerRecords[11].value = "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">";
    vcfStream.header.headerRecords[12].key = "FILTER";
    vcfStream.header.headerRecords[12].value = "<ID=q10,Description=\"Quality below 10\">";
    vcfStream.header.headerRecords[13].key = "FILTER";
    vcfStream.header.headerRecords[13].value = "<ID=s50,Description=\"Less than 50% of samples have data\">";
    vcfStream.header.headerRecords[14].key = "FORMAT";
    vcfStream.header.headerRecords[14].value = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    vcfStream.header.headerRecords[15].key = "FORMAT";
    vcfStream.header.headerRecords[15].value = "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
    vcfStream.header.headerRecords[16].key = "FORMAT";
    vcfStream.header.headerRecords[16].value = "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">";
    vcfStream.header.headerRecords[17].key = "FORMAT";
    vcfStream.header.headerRecords[17].value = "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">";

    // Write first record.
    {
        seqan::VcfRecord vcfRecord;
        vcfRecord.rID = 0;
        vcfRecord.beginPos = 14369;
        vcfRecord.id = "rs6054257";
        vcfRecord.ref = "G";
        vcfRecord.alt = "A";
        vcfRecord.qual = 29;
        vcfRecord.filter = "PASS";
        vcfRecord.info = "NS=3;DP=14;AF=0.5;DB;H2";
        vcfRecord.format = "GT:GQ:DP:HQ";
        appendValue(vcfRecord.genotypeInfos, "0|0:48:1:51,51");
        appendValue(vcfRecord.genotypeInfos, "1|0:48:8:51,51");
        appendValue(vcfRecord.genotypeInfos, "1/1:43:5:.,.");
        SEQAN_ASSERT_EQ(writeRecord(vcfStream, vcfRecord), 0);
    }

    // Write second record.
    {
        seqan::VcfRecord vcfRecord;
        vcfRecord.rID = 0;
        vcfRecord.beginPos = 17329;
        vcfRecord.id = ".";
        vcfRecord.ref = "T";
        vcfRecord.alt = "A";
        vcfRecord.qual = 3;
        vcfRecord.filter = "q10";
        vcfRecord.info = "NS=3;DP=11;AF=0.017";
        vcfRecord.format = "GT:GQ:DP:HQ";
        appendValue(vcfRecord.genotypeInfos, "0|0:49:3:58,50");
        appendValue(vcfRecord.genotypeInfos, "0|1:3:5:65,3");
        appendValue(vcfRecord.genotypeInfos, "0/0:41:3");
        SEQAN_ASSERT_EQ(writeRecord(vcfStream, vcfRecord), 0);
    }

    // Write third record.
    {
        seqan::VcfRecord vcfRecord;
        vcfRecord.rID = 0;
        vcfRecord.beginPos = 1110695;
        vcfRecord.id = "rs6040355";
        vcfRecord.ref = "A";
        vcfRecord.alt = "G,T";
        vcfRecord.qual = 67;
        vcfRecord.filter = "PASS";
        vcfRecord.info = "NS=2;DP=10;AF=0.333,0.667;AA=T;DB";
        vcfRecord.format = "GT:GQ:DP:HQ";
        appendValue(vcfRecord.genotypeInfos, "1|2:21:6:23,27");
        appendValue(vcfRecord.genotypeInfos, "2|1:2:0:18,2");
        appendValue(vcfRecord.genotypeInfos, "2/2:35:4");
        SEQAN_ASSERT_EQ(writeRecord(vcfStream, vcfRecord), 0);
    }

    close(vcfStream);

    seqan::CharString goldPath(SEQAN_PATH_TO_ROOT());
    append(goldPath, "/extras/tests/vcf_io/example.vcf");
    SEQAN_ASSERT(seqan::_compareTextFiles(toCString(tmpPath), toCString(goldPath)));
}

#endif  // SEQAN_EXTRAS_TESTS_VCF_TEST_VCF_IO_H_
