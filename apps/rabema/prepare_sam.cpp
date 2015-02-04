// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tool to prepare SAM file for Rabema.  Currently, this only consists of
// replacing all sequence and quality fields that store a "*" with the value
// from the primary alignment.
// ==========================================================================

#include <iostream>

#include <seqan/bam_io.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

#include "sorting.h"

struct Options
{
    // Disable checking for sortedness.  The hard requirement is on being clustered by read name and not being sorted.
    bool dontCheckSorting;

    // Input SAM file.
    seqan::CharString inputFile;
    // Output SAM file.
    seqan::CharString outputFile;

    Options() : dontCheckSorting(false)
    {}
};

int fixRecords(seqan::String<seqan::BamAlignmentRecord> & records)
{
    using namespace seqan;

    if (empty(records))
        return 0;  // OK to be empty.

    // Pick indices with sequences for first/second sequence.
    int idxSeqFirst = -1, idxSeqSecond = -1;
    for (unsigned i = 0; i < length(records); ++i)
    {
        if (hasFlagFirst(records[i]) && !empty(records[i].seq))
            idxSeqFirst = i;
        else if (hasFlagLast(records[i]) && !empty(records[i].seq))
            idxSeqSecond = i;
        else if (!hasFlagFirst(records[i]) && !hasFlagLast(records[i]) && !empty(records[i].seq))
            idxSeqFirst = i;
    }

    // Get sequences for first and second.
    Dna5String seqFirst;
    CharString qualFirst;
    if (idxSeqFirst != -1)
    {
        seqFirst = records[idxSeqFirst].seq;
        if (hasFlagRC(records[idxSeqFirst]))
            reverseComplement(seqFirst);
        qualFirst = records[idxSeqFirst].qual;
        if (hasFlagRC(records[idxSeqFirst]))
            reverse(qualFirst);
    }
    Dna5String seqFirstRC = seqFirst;
    reverseComplement(seqFirstRC);
    CharString qualFirstRC = qualFirst;
    reverse(qualFirstRC);
    Dna5String seqSecond;
    CharString qualSecond;
    if (idxSeqSecond != -1)
    {
        seqSecond = records[idxSeqSecond].seq;
        if (hasFlagRC(records[idxSeqSecond]))
            reverseComplement(seqSecond);
        qualSecond = records[idxSeqSecond].qual;
        if (hasFlagRC(records[idxSeqSecond]))
            reverse(qualSecond);
    }
    Dna5String seqSecondRC = seqSecond;
    reverseComplement(seqSecondRC);
    CharString qualSecondRC = qualSecond;
    reverse(qualSecond);

    // Actually assign sequences into records.
    for (unsigned i = 0; i < length(records); ++i)
    {
        if (hasFlagFirst(records[i]) || (!hasFlagFirst(records[i]) && !hasFlagLast(records[i])))
        {
            if (idxSeqFirst == -1)
            {
                std::cerr << "ERROR: No sequence for first mate of query name " << records[i].qName << ".\n";
                return 1;
            }
            if (hasFlagRC(records[i]))
            {
                if (!empty(records[i].seq) && records[i].seq != seqFirstRC)
                    SEQAN_FAIL("ERROR: Mismatching sequences for query name %s.", toCString(records[i].qName));
                records[i].seq = seqFirstRC;
                records[i].qual = qualFirstRC;
            }
            else
            {
                if (!empty(records[i].seq) && records[i].seq != seqFirst)
                    SEQAN_FAIL("ERROR: Mismatching sequences for query name %s.", toCString(records[i].qName));
                records[i].seq = seqFirst;
                records[i].qual = qualFirst;
            }
        }
        else if (hasFlagLast(records[i]))
        {
            if (idxSeqSecond == -1)
            {
                std::cerr << "ERROR: No sequence for second mate of query name " << records[i].qName << ".\n";
                return 1;
            }
            if (hasFlagRC(records[i]))
            {
                if (!empty(records[i].seq) && records[i].seq != seqSecondRC)
                    SEQAN_FAIL("ERROR: Mismatching sequences for query name %s.", toCString(records[i].qName));
                records[i].seq = seqSecondRC;
                records[i].qual = qualSecondRC;
            }
            else
            {
                if (!empty(records[i].seq) && records[i].seq != seqSecond)
                    SEQAN_FAIL("ERROR: Mismatching sequences for query name %s.", toCString(records[i].qName));
                records[i].seq = seqSecond;
                records[i].qual = qualSecond;
            }
        }
    }

    return 0;
}

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("rabema_prepare_sam");

    // Set short description, version, and date.
    setShortDescription(parser, "Prepare SAM For Rabema");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setCategory(parser, "Benchmarking");

    // Define usage line and long description.
    addUsageLine(parser, "\\fB-i\\fP \\fIIN.sam\\fP \\fB-o\\fP \\fIOUT.sam\\fP");
    addDescription(parser, "Prepare SAM file for usage with RABEMA.");

    // Define Options.
    addOption(parser, seqan::ArgParseOption(
            "i", "in-file", "Path to the input file.",
            seqan::ArgParseArgument::INPUT_FILE, "IN.sam"));
    setValidValues(parser, "in-file", "sam");
    setRequired(parser, "in-file");

    addOption(parser, seqan::ArgParseOption(
            "o", "out-file", "Path to the output file.",
            seqan::ArgParseArgument::OUTPUT_FILE, "OUT.sam"));
    setValidValues(parser, "out-file", "sam");
    setRequired(parser, "out-file");

    addOption(parser, seqan::ArgParseOption("", "dont-check-sorting", "Do not check sortedness."));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getOptionValue(options.inputFile, parser, "in-file");
    getOptionValue(options.outputFile, parser, "out-file");
    options.dontCheckSorting = isSet(parser, "dont-check-sorting");

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv)
{
    using namespace seqan;

    // Parse command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Open SAM file for reading.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(options.inputFile)))
    {
        std::cerr << "Could not open file " << argv[1] << " for reading.\n";
        return 1;
    }

    // Read header.
    BamHeader header;
    try
    {
        readHeader(header, bamFileIn);
    }
    catch (ParseError const & ioErr)
    {
        std::cerr << "ERROR: Could not read header from SAM file (" << ioErr.what() << ").\n";
        return 1;
    }

    BamFileOut bamFileOut(bamFileIn);
    if (!open(bamFileOut, toCString(options.outputFile)))
    {
        std::cerr << "Could not open file " << options.outputFile << " for writing.\n";
        return 1;
    }

    try
    {
        writeHeader(bamFileOut, header);
    }
    catch (ParseError const & ioErr)
    {
        std::cerr << "ERROR writing SAM header!\n";
        return 1;
    }

    // Read file in chunks, one for each query name.
    String<BamAlignmentRecord> records;
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        try
        {
            readRecord(record, bamFileIn);
        }
        catch (ParseError const & ioError)
        {
            std::cerr << "ERROR reading SAM record!\n";
            return 1;
        }

        if (!empty(records) && record.qName != back(records).qName)
        {
            // Sanity check for sorting.
            if (!options.dontCheckSorting && !empty(record) && lessThanSamtoolsQueryName(record.qName, back(records).qName))
            {
                std::cerr << "ERROR: " << record.qName << " succeeds " << back(records).qName << " in SAM file.\n"
                          << "File must be sorted by query name.\n"
                          << "If your input FASTQ file was not sorted by query name then you can disable this \n"
                          << "sanity check using the --dont-check-sorting parameter.\n";
                return 1;
            }
            if (fixRecords(records) != 0)
            {
                std::cerr << "Could not fix records!\n";
                return 1;
            }
            writeRecords(bamFileOut, records);
            clear(records);
        }

        appendValue(records, record);
    }
    if (fixRecords(records) != 0)
    {
        std::cerr << "Could not fix records!\n";
        return 1;
    }
    try
    {
        for (unsigned i = 0; i < length(records); ++i)
            writeRecord(bamFileOut, records[i]);
    }
    catch (IOError const & ioError)
    {
        std::cerr << "ERROR writing SAM record!\n";
        return 1;
    }

    return 0;
}
