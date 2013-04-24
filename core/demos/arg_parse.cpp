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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.

#include <seqan/arg_parse.h>

using seqan::ArgumentParser;
using seqan::ArgParseOption;
using seqan::ArgParseArgument;
using seqan::CharString;

int main(int argc, char const ** argv)
{
    ArgumentParser parser("arg_parse_demo");
    setShortDescription(parser, "Just a demo of the new seqan::ArgumentParser!");
    setVersion(parser, "0.1");
    setDate(parser, "Mar 2012");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIIN\\fP \\fIOUT\\fP ");

    addDescription(parser, "This is just a little demo to show what seqan::ArgumentParser is able to do.");
    addDescription(parser, "\\fIIN\\fP is a multi-FASTA input.");
    addDescription(parser, "\\fIOUT\\fP is a txt output.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "IN"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE, "OUT"));

    // allow only fasta files as input
    setValidValues(parser, 0, "FASTA fa");
    setValidValues(parser, 1, "txt");

    addSection(parser, "Important Tool Parameters");
    addOption(parser, ArgParseOption("", "id", "Sequence identity between [0.0:1.0]", ArgParseArgument::DOUBLE, "ID"));
    setRequired(parser, "id", true);
    setMinValue(parser, "id", "0.0");
    setMaxValue(parser, "id", "1.0");

    addSection(parser, "Miscellaneous");
    addOption(parser, ArgParseOption("v", "verbose", "Turn on verbose output."));
    addOption(parser, ArgParseOption("H", "hidden", "Super mysterious flag that will not be shown in the help screen or man-page."));
    hideOption(parser, "H");

    addTextSection(parser, "References");
    addText(parser, "http://www.seqan.de");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == ArgumentParser::PARSE_OK)
    {
        bool verbose = false;
        getOptionValue(verbose, parser, "verbose");
        std::cout << "Verbose:     " << (verbose ? "on" : "off") << std::endl;

        double identity = -1.0;
        getOptionValue(identity, parser, "id");
        std::cout << "Identity:    " << identity << std::endl;

        CharString inputFile, outputFile;
        getArgumentValue(inputFile, parser, 0);
        getArgumentValue(outputFile, parser, 1);

        std::cout << "Input-File:  " << inputFile << std::endl;
        std::cout << "Output-File: " << outputFile << std::endl;

        return 0;
    }
    else
    {
        return res == ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise
    }
}
