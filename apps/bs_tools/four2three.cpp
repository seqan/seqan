// ==========================================================================
//                              four2three
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
// Author: Sabrina Krakau <sabrina.krakau@fu-berlin.de>
// ==========================================================================

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/arg_parse.h>
#include "four2three.h"

using namespace std;
using namespace seqan;

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    CharString inputFileName;
    bool ctConversion;

    CharString outputFileName;

    AppOptions() :
        verbosity(1),
        ctConversion(true)
    {}
};

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("four2three");
    // Set short description, version, and date.
    setShortDescription(parser, "Four to three-letter alphabet reduction.");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setCategory(parser, "BS-Seq Analysis");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fISEQUENCE FILE\\fP\"");
    addDescription(parser, "This program converts four-letter sequences into three-letter sequences for bisulfite sequence analysis.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "SEQUENCES"));
    setHelpText(parser, 0, "A sequence file containing reads or genome.");
    setValidValues(parser, 0, SeqFileIn::getFileExtensions());

    addOption(parser, ArgParseOption("o", "output-file", "Name of output file.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "o", SeqFileOut::getFileExtensions());
    setRequired(parser, "output-file", true);
    addOption(parser, ArgParseOption("ga", "ga-conversion", "Convert Gs to As, instead of Cs to Ts."));

    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBfour2three\\fP \\fB-o\\fP \\fBreads_CT.fastq\\fP \\fBreads.fastq\\fP ",
                "Converts sequences into three-letter alphabet by converting Cs to Ts.");
    addListItem(parser, "\\fBfour2three\\fP \\fB-ga\\fP \\fB-o\\fP \\fBgenome_GA.fa\\fP \\fBgenome.fa\\fP ",
                "Converts sequences into three-letter alphabet by converting Gs to As.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.inputFileName, parser, 0);
    getOptionValue(options.outputFileName, parser, "output-file");

    if (isSet(parser, "ga-conversion"))
        options.ctConversion = false;

    CharString tmp1 = options.inputFileName;
    toLower(tmp1);
    CharString tmp2 = options.outputFileName;
    toLower(tmp2);

    if ( ( (endsWith(tmp1, ".fa") || endsWith(tmp1, ".fasta")) &&  (endsWith(tmp2, ".fastq") || endsWith(tmp2, ".fq")) ) ||
         ( (endsWith(tmp2, ".fa") || endsWith(tmp2, ".fasta")) &&  (endsWith(tmp1, ".fastq") || endsWith(tmp2, ".fq")) )  )
    {
        std::cerr << "ERROR: Output file must have the same format as input file!" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    return ArgumentParser::PARSE_OK;
}


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

    preProcess(options);
    return 0;
}



