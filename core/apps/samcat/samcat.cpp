// ==========================================================================
//                                   samcat
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/misc/misc_name_store_cache.h>

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
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function mergeBamFiles()
// --------------------------------------------------------------------------

template <typename TWriter>
void mergeBamFiles(TWriter &writer, StringSet<CharString> &inFiles)
{
    // Step 1: Merge all headers (if available)
    BamHeader header;
    for (unsigned i = 0; i < length(inFiles); ++i)
    {
        BamFile<Input> reader(writer, toCString(inFiles[i]));
        readRecord(header, reader);
    }

    // Step 2: Remove duplicate header entries and write merged header
    removeDuplicates(header);
    write(writer, header);

    // Step 3: Read and output alignment records
    BamAlignmentRecord record;
    for (unsigned i = 0; i != length(inFiles); ++i)
    {
        BamFile<Input> reader(writer, toCString(inFiles[i]));

        // skip header
        readRecord(header, reader);

        while (!atEnd(reader))
        {
            readRecord(record, reader);
            write(writer, record);
        }
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
    setShortDescription(parser, "SAM/BAM file concatenation and conversion");
    setCategory(parser, "Utilities");
    setVersion(parser, "0.1");
    setDate(parser, "May 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIINFILE\\fP> [<\\fIINFILE\\fP> ...] [-o <\\fIOUTFILE\\fP>]");
    addDescription(parser, "This tool reads a set of input files in SAM or BAM format and outputs the concatenation of them. "
                           "If the output file name is ommitted the result is written to standard output in SAM format.");

    addDescription(parser, "(c) Copyright 2014 by David Weese.");

    // We require one argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "INFILE", true));
    setValidValues(parser, 0, ".sam .bam");
    setHelpText(parser, 0, "Input SAM or BAM file.");

    addOption(parser, ArgParseOption("o", "output", "Output file name", ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "output", ".sam .bam");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBsamcat\\fP \\fBmapped1.sam\\fP \\fBmapped2.sam\\fP \\fB-o\\fP \\fBmerged.sam\\fP",
                "Merge two SAM files.");
    addListItem(parser, "\\fBsamcat\\fP \\fBinput.sam\\fP \\fB-o\\fP \\fBouput.bam\\fP",
                "Convert a SAM file into BAM format.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    options.inFiles = getArgumentValues(parser, 0);
    getOptionValue(options.outFile, parser, "output");

    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    BamFile<Output> writer(toCString(options.outFile));
    mergeBamFiles(writer, options.inFiles);
//    else
//    {
//        // dump to standard out stream
//        return mergeBamFiles(std::cout, Sam(), options.inFiles);
//    }
}
