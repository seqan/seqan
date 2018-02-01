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

#ifndef SEQAN_TESTS_VCF_TEST_VCF_IO_H_
#define SEQAN_TESTS_VCF_TEST_VCF_IO_H_

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>


SEQAN_DEFINE_TEST(test_vcf_io_read_vcf_header)
{
    seqan::CharString vcfPath = seqan::getAbsolutePath("/tests/vcf_io/example.vcf");

    seqan::String<char, seqan::MMap<> > mmapString;
    SEQAN_ASSERT(open(mmapString, toCString(vcfPath)));
    seqan::Iterator<seqan::String<char, seqan::MMap<> >, seqan::Rooted>::Type iter = begin(mmapString);

    seqan::VcfIOContext<> vcfIOContext;
    seqan::VcfHeader vcfHeader;

    readHeader(vcfHeader, vcfIOContext, iter, seqan::Vcf());

    SEQAN_ASSERT_EQ(length(vcfHeader), 18u);
    SEQAN_ASSERT_EQ(vcfHeader[0].key, "fileformat");
    SEQAN_ASSERT_EQ(vcfHeader[0].value, "VCFv4.1");
    SEQAN_ASSERT_EQ(vcfHeader[1].key, "fileDate");
    SEQAN_ASSERT_EQ(vcfHeader[1].value, "20090805");
    SEQAN_ASSERT_EQ(vcfHeader[2].key, "source");
    SEQAN_ASSERT_EQ(vcfHeader[2].value, "myImputationProgramV3.1");
    SEQAN_ASSERT_EQ(vcfHeader[3].key, "reference");
    SEQAN_ASSERT_EQ(vcfHeader[3].value, "file:///seq/references/1000GenomesPilot-NCBI36.fasta");
    SEQAN_ASSERT_EQ(vcfHeader[4].key, "contig");
    SEQAN_ASSERT_EQ(vcfHeader[4].value, "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>");
    SEQAN_ASSERT_EQ(vcfHeader[5].key, "phasing");
    SEQAN_ASSERT_EQ(vcfHeader[5].value, "partial");
    SEQAN_ASSERT_EQ(vcfHeader[6].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader[6].value, "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    SEQAN_ASSERT_EQ(vcfHeader[7].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader[7].value, "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    SEQAN_ASSERT_EQ(vcfHeader[8].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader[8].value, "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
    SEQAN_ASSERT_EQ(vcfHeader[9].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader[9].value, "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
    SEQAN_ASSERT_EQ(vcfHeader[10].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader[10].value, "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">");
    SEQAN_ASSERT_EQ(vcfHeader[11].key, "INFO");
    SEQAN_ASSERT_EQ(vcfHeader[11].value, "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">");
    SEQAN_ASSERT_EQ(vcfHeader[12].key, "FILTER");
    SEQAN_ASSERT_EQ(vcfHeader[12].value, "<ID=q10,Description=\"Quality below 10\">");
    SEQAN_ASSERT_EQ(vcfHeader[13].key, "FILTER");
    SEQAN_ASSERT_EQ(vcfHeader[13].value, "<ID=s50,Description=\"Less than 50% of samples have data\">");
    SEQAN_ASSERT_EQ(vcfHeader[14].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfHeader[14].value, "<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    SEQAN_ASSERT_EQ(vcfHeader[15].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfHeader[15].value, "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    SEQAN_ASSERT_EQ(vcfHeader[16].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfHeader[16].value, "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    SEQAN_ASSERT_EQ(vcfHeader[17].key, "FORMAT");
    SEQAN_ASSERT_EQ(vcfHeader[17].value, "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");

    SEQAN_ASSERT_EQ(length(contigNames(vcfIOContext)), 1u);
    SEQAN_ASSERT_EQ(contigNames(vcfIOContext)[0], "20");

    SEQAN_ASSERT_EQ(length(sampleNames(vcfIOContext)), 3u);
    SEQAN_ASSERT_EQ(sampleNames(vcfIOContext)[0], "NA00001");
    SEQAN_ASSERT_EQ(sampleNames(vcfIOContext)[1], "NA00002");
    SEQAN_ASSERT_EQ(sampleNames(vcfIOContext)[2], "NA00003");
}


SEQAN_DEFINE_TEST(test_vcf_io_read_vcf_record)
{
    seqan::CharString vcfPath = seqan::getAbsolutePath("/tests/vcf_io/example_records_with_errors.vcf");

    seqan::String<char, seqan::MMap<> > mmapString;
    open(mmapString, toCString(vcfPath));
    seqan::Iterator<seqan::String<char, seqan::MMap<> >, seqan::Rooted>::Type iter = begin(mmapString);

    seqan::VcfIOContext<> vcfIOContext;
    seqan::VcfHeader vcfHeader;

    resize(sampleNames(vcfIOContext), 3);

    seqan::String<seqan::VcfRecord> records;
    seqan::VcfRecord record;
    for (unsigned i = 0; i < 5; ++i)
    {
        readRecord(record, vcfIOContext, iter, seqan::Vcf());
        appendValue(records, record);
    }

    SEQAN_ASSERT_EQ(length(records), 5u);

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

    // the next 2 recodrs are valid since vcf v4.2
    SEQAN_ASSERT_EQ(records[3].rID, 0);
    SEQAN_ASSERT_EQ(records[3].beginPos, 1110695);
    SEQAN_ASSERT_EQ(records[3].id, "rs6040355");
    SEQAN_ASSERT_EQ(records[3].ref, "A");
    SEQAN_ASSERT_EQ(records[3].alt, "G,T");
    SEQAN_ASSERT_EQ(records[3].qual, 67);
    SEQAN_ASSERT_EQ(records[3].filter, "PASS");
    SEQAN_ASSERT_EQ(records[3].info, "NS=2;DP=10;AF=0.333,0.667;AA=T;DB");
    SEQAN_ASSERT_EQ(records[3].format, ""); // empty formats are accepted since v4.2
    SEQAN_ASSERT_EQ(length(records[3].genotypeInfos), 3u);

    SEQAN_ASSERT_EQ(records[4].rID, 0);
    SEQAN_ASSERT_EQ(records[4].beginPos, 1110695);
    SEQAN_ASSERT_EQ(records[4].id, "rs6040355");
    SEQAN_ASSERT_EQ(records[4].ref, "A");
    SEQAN_ASSERT_EQ(records[4].alt, "G,T");
    SEQAN_ASSERT_EQ(records[4].qual, 67);
    SEQAN_ASSERT_EQ(records[4].filter, "PASS");
    SEQAN_ASSERT_EQ(records[4].info, "NS=2;DP=10;AF=0.333,0.667;AA=T;DB");
    SEQAN_ASSERT_EQ(records[4].format, ""); // empty formats are accepted since v4.2
    SEQAN_ASSERT_EQ(length(records[4].genotypeInfos), 3u);

    // the next 18 records are invalid and readRecord should throw ParseError
    // continuing to read after EOF file should also result in ParseError
    for (unsigned i = 0; i < 25; ++i)
    {
        SEQAN_TEST_EXCEPTION(seqan::ParseError,
                             seqan::readRecord(record, vcfIOContext, iter, seqan::Vcf()));
    }
}

SEQAN_DEFINE_TEST(test_vcf_io_vcf_file_read_record)
{
    seqan::CharString vcfPath = seqan::getAbsolutePath("/tests/vcf_io/example.vcf");

    seqan::VcfFileIn vcfStream(toCString(vcfPath));
    seqan::VcfHeader header;

    readHeader(header, vcfStream);

    SEQAN_ASSERT_EQ(length(header), 18u);
    SEQAN_ASSERT_EQ(header[0].key, "fileformat");
    SEQAN_ASSERT_EQ(header[0].value, "VCFv4.1");
    SEQAN_ASSERT_EQ(header[1].key, "fileDate");
    SEQAN_ASSERT_EQ(header[1].value, "20090805");
    SEQAN_ASSERT_EQ(header[2].key, "source");
    SEQAN_ASSERT_EQ(header[2].value, "myImputationProgramV3.1");
    SEQAN_ASSERT_EQ(header[3].key, "reference");
    SEQAN_ASSERT_EQ(header[3].value, "file:///seq/references/1000GenomesPilot-NCBI36.fasta");
    SEQAN_ASSERT_EQ(header[4].key, "contig");
    SEQAN_ASSERT_EQ(header[4].value, "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>");
    SEQAN_ASSERT_EQ(header[5].key, "phasing");
    SEQAN_ASSERT_EQ(header[5].value, "partial");
    SEQAN_ASSERT_EQ(header[6].key, "INFO");
    SEQAN_ASSERT_EQ(header[6].value, "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    SEQAN_ASSERT_EQ(header[7].key, "INFO");
    SEQAN_ASSERT_EQ(header[7].value, "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    SEQAN_ASSERT_EQ(header[8].key, "INFO");
    SEQAN_ASSERT_EQ(header[8].value, "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
    SEQAN_ASSERT_EQ(header[9].key, "INFO");
    SEQAN_ASSERT_EQ(header[9].value, "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
    SEQAN_ASSERT_EQ(header[10].key, "INFO");
    SEQAN_ASSERT_EQ(header[10].value, "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">");
    SEQAN_ASSERT_EQ(header[11].key, "INFO");
    SEQAN_ASSERT_EQ(header[11].value, "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">");
    SEQAN_ASSERT_EQ(header[12].key, "FILTER");
    SEQAN_ASSERT_EQ(header[12].value, "<ID=q10,Description=\"Quality below 10\">");
    SEQAN_ASSERT_EQ(header[13].key, "FILTER");
    SEQAN_ASSERT_EQ(header[13].value, "<ID=s50,Description=\"Less than 50% of samples have data\">");
    SEQAN_ASSERT_EQ(header[14].key, "FORMAT");
    SEQAN_ASSERT_EQ(header[14].value, "<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    SEQAN_ASSERT_EQ(header[15].key, "FORMAT");
    SEQAN_ASSERT_EQ(header[15].value, "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    SEQAN_ASSERT_EQ(header[16].key, "FORMAT");
    SEQAN_ASSERT_EQ(header[16].value, "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    SEQAN_ASSERT_EQ(header[17].key, "FORMAT");
    SEQAN_ASSERT_EQ(header[17].value, "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");

    SEQAN_ASSERT_EQ(length(contigNames(context(vcfStream))), 1u);
    SEQAN_ASSERT_EQ(contigNames(context(vcfStream))[0], "20");

    SEQAN_ASSERT_EQ(length(sampleNames(context(vcfStream))), 3u);
    SEQAN_ASSERT_EQ(sampleNames(context(vcfStream))[0], "NA00001");
    SEQAN_ASSERT_EQ(sampleNames(context(vcfStream))[1], "NA00002");
    SEQAN_ASSERT_EQ(sampleNames(context(vcfStream))[2], "NA00003");

    seqan::String<seqan::VcfRecord> records;
    while (!atEnd(vcfStream))
    {
        seqan::VcfRecord record;
        readRecord(record, vcfStream);
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
    seqan::VcfIOContext<> vcfIOContext;
    seqan::VcfHeader vcfHeader;
    appendName(contigNamesCache(vcfIOContext), "20");
    appendName(sampleNamesCache(vcfIOContext), "NA00001");
    appendName(sampleNamesCache(vcfIOContext), "NA00002");
    appendName(sampleNamesCache(vcfIOContext), "NA00003");

    resize(vcfHeader, 18);
    vcfHeader[0].key = "fileformat";
    vcfHeader[0].value = "VCFv4.1";
    vcfHeader[1].key = "fileDate";
    vcfHeader[1].value = "20090805";
    vcfHeader[2].key = "source";
    vcfHeader[2].value = "myImputationProgramV3.1";
    vcfHeader[3].key = "reference";
    vcfHeader[3].value = "file:///seq/references/1000GenomesPilot-NCBI36.fasta";
    vcfHeader[4].key = "contig";
    vcfHeader[4].value = "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>";
    vcfHeader[5].key = "phasing";
    vcfHeader[5].value = "partial";
    vcfHeader[6].key = "INFO";
    vcfHeader[6].value = "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">";
    vcfHeader[7].key = "INFO";
    vcfHeader[7].value = "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
    vcfHeader[8].key = "INFO";
    vcfHeader[8].value = "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">";
    vcfHeader[9].key = "INFO";
    vcfHeader[9].value = "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">";
    vcfHeader[10].key = "INFO";
    vcfHeader[10].value = "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">";
    vcfHeader[11].key = "INFO";
    vcfHeader[11].value = "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">";
    vcfHeader[12].key = "FILTER";
    vcfHeader[12].value = "<ID=q10,Description=\"Quality below 10\">";
    vcfHeader[13].key = "FILTER";
    vcfHeader[13].value = "<ID=s50,Description=\"Less than 50% of samples have data\">";
    vcfHeader[14].key = "FORMAT";
    vcfHeader[14].value = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    vcfHeader[15].key = "FORMAT";
    vcfHeader[15].value = "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
    vcfHeader[16].key = "FORMAT";
    vcfHeader[16].value = "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">";
    vcfHeader[17].key = "FORMAT";
    vcfHeader[17].value = "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">";

    std::string tmpPath = (std::string)SEQAN_TEMP_FILENAME() + ".vcf";
    std::ofstream file(tmpPath.c_str());
    seqan::DirectionIterator<std::ofstream, seqan::Output>::Type iter = directionIterator(file, seqan::Output());
    writeHeader(iter, vcfHeader, vcfIOContext, seqan::Vcf());
    file.close();

    std::string goldPath = seqan::getAbsolutePath("/tests/vcf_io/vcf_header.vcf");
    SEQAN_ASSERT(seqan::_compareTextFilesAlt(tmpPath.c_str(), goldPath.c_str()));
}


SEQAN_DEFINE_TEST(test_vcf_io_write_vcf_record)
{
    std::string goldPath = seqan::getAbsolutePath("/tests/vcf_io/example_records.vcf");
    std::string tmpPath = (std::string)SEQAN_TEMP_FILENAME() + ".vcf";

    std::ifstream file(goldPath.c_str());
    std::ofstream fileOut(tmpPath.c_str());

    seqan::DirectionIterator<std::ifstream, seqan::Input>::Type iter = directionIterator(file, seqan::Input());
    seqan::DirectionIterator<std::ofstream, seqan::Output>::Type iterOut = directionIterator(fileOut, seqan::Output());

    seqan::VcfHeader vcfHeader;
    seqan::VcfIOContext<> vcfIOContext;
    resize(sampleNames(vcfIOContext), 3);

    seqan::VcfRecord record;
    while (!atEnd(iter))
    {
        readRecord(record, vcfIOContext, iter, seqan::Vcf());
        writeRecord(iterOut, record, vcfIOContext, seqan::Vcf());
    }
    fileOut.close();

    SEQAN_ASSERT(seqan::_compareTextFilesAlt(tmpPath.c_str(), goldPath.c_str()));
}

SEQAN_DEFINE_TEST(test_vcf_io_vcf_file_write_record)
{
    std::string tmpPath = (std::string)SEQAN_TEMP_FILENAME() + ".vcf";
    seqan::VcfFileOut vcfStream(tmpPath.c_str());

    // Build header.
    appendName(contigNamesCache(context(vcfStream)), "20");
    appendName(sampleNamesCache(context(vcfStream)), "NA00001");
    appendName(sampleNamesCache(context(vcfStream)), "NA00002");
    appendName(sampleNamesCache(context(vcfStream)), "NA00003");

    seqan::VcfHeader header;
    resize(header, 18);
    header[0].key = "fileformat";
    header[0].value = "VCFv4.1";
    header[1].key = "fileDate";
    header[1].value = "20090805";
    header[2].key = "source";
    header[2].value = "myImputationProgramV3.1";
    header[3].key = "reference";
    header[3].value = "file:///seq/references/1000GenomesPilot-NCBI36.fasta";
    header[4].key = "contig";
    header[4].value = "<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>";
    header[5].key = "phasing";
    header[5].value = "partial";
    header[6].key = "INFO";
    header[6].value = "<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">";
    header[7].key = "INFO";
    header[7].value = "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
    header[8].key = "INFO";
    header[8].value = "<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">";
    header[9].key = "INFO";
    header[9].value = "<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">";
    header[10].key = "INFO";
    header[10].value = "<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">";
    header[11].key = "INFO";
    header[11].value = "<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">";
    header[12].key = "FILTER";
    header[12].value = "<ID=q10,Description=\"Quality below 10\">";
    header[13].key = "FILTER";
    header[13].value = "<ID=s50,Description=\"Less than 50% of samples have data\">";
    header[14].key = "FORMAT";
    header[14].value = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    header[15].key = "FORMAT";
    header[15].value = "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
    header[16].key = "FORMAT";
    header[16].value = "<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">";
    header[17].key = "FORMAT";
    header[17].value = "<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">";

    // Write header
    writeHeader(vcfStream, header);

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
        writeRecord(vcfStream, vcfRecord);
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
        writeRecord(vcfStream, vcfRecord);
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
        writeRecord(vcfStream, vcfRecord);
    }

    close(vcfStream);

    seqan::CharString goldPath(seqan::getAbsolutePath("/tests/vcf_io/example.vcf"));
    SEQAN_ASSERT(seqan::_compareTextFilesAlt(tmpPath.c_str(), toCString(goldPath)));
}

SEQAN_DEFINE_TEST(test_vcf_io_isOpen_fileIn)
{
    // Build path to file.
    seqan::CharString vcfPath = SEQAN_PATH_TO_ROOT();
    append(vcfPath, "/tests/vcf_io/example.vcf");

    // Create SequenceStream object.
    seqan::VcfFileIn vcfI;
    SEQAN_ASSERT(!isOpen(vcfI));

    // open file
    open(vcfI, toCString(vcfPath));
    SEQAN_ASSERT(isOpen(vcfI));

    // close file
    close(vcfI);
    SEQAN_ASSERT(!isOpen(vcfI));
}

SEQAN_DEFINE_TEST(test_vcf_io_isOpen_fileOut)
{
    // Build path to file.
    seqan::CharString vcfPath = SEQAN_TEMP_FILENAME();
    append(vcfPath, ".vcf");

    // Create SequenceStream object.
    seqan::VcfFileOut  vcfO;
    SEQAN_ASSERT(!isOpen(vcfO));

    // open files
    open(vcfO, toCString(vcfPath));
    SEQAN_ASSERT(isOpen(vcfO));

    // close files
    close(vcfO);
    SEQAN_ASSERT(!isOpen(vcfO));
}

#endif  // SEQAN_TESTS_VCF_TEST_VCF_IO_H_
