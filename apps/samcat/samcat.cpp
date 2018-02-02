// ==========================================================================
//                                   samcat
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

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/arg_parse.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    StringSet<CharString> inFiles;
    CharString outFile;
    bool bamFormat;
    bool verbose;
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function mergeBamFiles()
// --------------------------------------------------------------------------

template <typename TWriter>
void catBamFiles(TWriter &writer, StringSet<CharString> &inFiles, AppOptions const &options)
{
    String<BamFileIn *> readerPtr;
    resize(readerPtr, length(inFiles));

    // Step 1: Merge all headers (if available)
    BamHeader header;
    for (unsigned i = 0; i < length(inFiles); ++i)
    {
        readerPtr[i] = new BamFileIn(writer);

        bool success;
        if (inFiles[i] != "-")
            success = open(*readerPtr[i], toCString(inFiles[i]));
        else
            // read from stdin (autodetect format from stream)
            success = open(*readerPtr[i], std::cin);
        if (!success)
        {
            std::cerr << "Couldn't open " << toCString(inFiles[i]) << " for reading." << std::endl;
            close(*readerPtr[i]);
            delete readerPtr[i];
            readerPtr[i] = NULL;
            continue;
        }

        readHeader(header, *(readerPtr[i]));
    }

    // Step 2: Remove duplicate header entries and write merged header
    if (length(inFiles) > 1)
        removeDuplicates(header);
    writeHeader(writer, header);

    // Step 3: Read and output alignment records
    BamAlignmentRecord record;
    String<BamAlignmentRecord> records;
    uint64_t numRecords = 0;
    double start = sysTime();
    for (unsigned i = 0; i != length(inFiles); ++i)
    {
        if (readerPtr[i] == NULL)
            continue;

        BamFileIn &reader = *readerPtr[i];

        // copy all alignment records
        if (isEqual(format(writer), Sam()))
        {
            // For Sam parallel batch processing is faster
            while (!atEnd(reader))
            {
                unsigned size = readRecords(records, reader, 100000);
                writeRecords(writer, prefix(records, size));
                numRecords += size;
            }
        }
        else
        {
            while (!atEnd(reader))
            {
                readRecord(record, reader);
                writeRecord(writer, record);
                ++numRecords;
            }
        }
        close(reader);
        delete readerPtr[i];
    }
    double stop = sysTime();
    if (options.verbose)
    {
        std::cerr << "Number of alignments: " << numRecords << std::endl;
        std::cerr << "Elapsed time:         " << stop - start << " seconds" << std::endl;
    }
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("samcat");
    // Set short description, version, and date.
    setCategory(parser, "Utilities");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIINFILE\\fP> [<\\fIINFILE\\fP> ...] [-o <\\fIOUTFILE\\fP>]");
#if SEQAN_HAS_ZLIB
    setShortDescription(parser, "SAM/BAM file concatenation and conversion");
    addDescription(parser, "This tool reads a set of input files in SAM or BAM format "
#else
    setShortDescription(parser, "SAM file concatenation and conversion");
    addDescription(parser, "This tool reads a set of input files in SAM format "
#endif
                           "and outputs the concatenation of them. "
                           "If the output file name is omitted the result is written to stdout.");

    addDescription(parser, "(c) Copyright in 2014 by David Weese.");

    addOption(parser, ArgParseOption("o", "output", "Output file name.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output", BamFileOut::getFileExtensions());

    // We require one argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "INFILE", true));
    setValidValues(parser, 0, BamFileIn::getFileExtensions());
#if SEQAN_HAS_ZLIB
    setHelpText(parser, 0, "Input SAM or BAM file (or - for stdin).");
    addOption(parser, ArgParseOption("b", "bam", "Use BAM format for standard output. Default: SAM."));
#else
    setHelpText(parser, 0, "Input SAM file (or - for stdin).");
#endif
    addOption(parser, ArgParseOption("v", "verbose", "Print some stats."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBsamcat\\fP \\fBmapped1.sam\\fP \\fBmapped2.sam\\fP \\fB-o\\fP \\fBmerged.sam\\fP",
                "Merge two SAM files.");
#if SEQAN_HAS_ZLIB
    addListItem(parser, "\\fBsamcat\\fP \\fBinput.sam\\fP \\fB-o\\fP \\fBoutput.bam\\fP",
                "Convert a SAM file into BAM format.");
#endif

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    options.inFiles = getArgumentValues(parser, 0);
    getOptionValue(options.outFile, parser, "output");
#if SEQAN_HAS_ZLIB
    getOptionValue(options.bamFormat, parser, "bam");
#endif
    getOptionValue(options.verbose, parser, "verbose");

    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    bool success;
    BamFileOut writer;
    if (!empty(options.outFile))
        success = open(writer, toCString(options.outFile));
    else
        // write to stdout
#if SEQAN_HAS_ZLIB
        if (options.bamFormat)
            success = open(writer, std::cout, Bam());
        else
#endif
            success = open(writer, std::cout, Sam());

    if (!success)
    {
        std::cerr << "Couldn't open " << options.outFile << " for writing." << std::endl;
        return 1;
    }

    catBamFiles(writer, options.inFiles, options);
    return 0;
}
