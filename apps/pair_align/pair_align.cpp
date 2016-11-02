/*==========================================================================
               SeqAn - The Library for Sequence Analysis
                         http://www.seqan.de
============================================================================
Copyright (C) 2007-2012
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
==========================================================================*/

#include "lib/pair_align_lib.h"

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("pair_align");
    // Set short description, version, and date.
    setShortDescription(parser, "Pairwise alignment");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-s\\fP \\fIIN\\fP");
    setCategory(parser, "Sequence Alignment");
    addDescription(parser,
                   "The program allows one to align two sequences using dyamic programming alignment algorithms while "
                   "tweaking various parameters.");

    addSection(parser, "Main Options");
    addOption(parser, seqan::ArgParseOption("s", "seq", "FASTA file with two sequences.", seqan::ArgParseOption::INPUT_FILE, "IN"));
    setRequired(parser, "seq");
    setValidValues(parser, "seq", getFileExtensions(Fasta()));
    addOption(parser, seqan::ArgParseOption("a", "alphabet", "Sequence alphabet.", seqan::ArgParseOption::STRING, "ALPHABET"));
    setValidValues(parser, "alphabet", "protein dna rna text");
    setDefaultValue(parser, "alphabet", "protein");
    addOption(parser, seqan::ArgParseOption("m", "method",
                                            "DP alignment method: Needleman-Wunsch, Gotoh, Smith-Waterman, "
                                            "Longest Common Subsequence",
                                            seqan::ArgParseOption::STRING, "METHOD"));
    setValidValues(parser, "method", "nw gotoh sw lcs");
    setDefaultValue(parser, "method", "gotoh");
    addOption(parser, seqan::ArgParseOption("o", "outfile", "Output filename.", seqan::ArgParseOption::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "outfile", "out.fasta");
    std::vector<std::string> outFileNames = getFileExtensions(Fasta());
    outFileNames.push_back(".msf");
    setValidValues(parser, "outfile", outFileNames);

    addSection(parser, "Scoring Options");
    addOption(parser, seqan::ArgParseOption("g", "gop", "Gap open penalty.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "gop", "-11");
    addOption(parser, seqan::ArgParseOption("e", "gex", "Gap extension penalty.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "gex", "-1");
    addOption(parser, seqan::ArgParseOption("ma", "matrix", "Score matrix.", seqan::ArgParseOption::STRING, "MATRIX_FILE"));
    addOption(parser, seqan::ArgParseOption("ms", "msc", "Match score.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "msc", "5");
    addOption(parser, seqan::ArgParseOption("mm", "mmsc", "Mismatch penalty.", seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "mmsc", "-4");

    addSection(parser, "Banded Alignment Options");
    addOption(parser, seqan::ArgParseOption("lo", "low", "Lower diagonal.", seqan::ArgParseOption::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("hi", "high", "Upper diagonal.", seqan::ArgParseOption::INTEGER, "INT"));

    addSection(parser, "DP Matrix Configuration Options");
    addOption(parser, seqan::ArgParseOption("c", "config", "Alignment configuration.", seqan::ArgParseOption::STRING, "CONF"));
    setValidValues(parser, "config", "ffff ffft fftf fftt ftff ftft fttf fttt tfff tfft tftf tftt ttff ttft tttf tttt");

    addTextSection(parser, "Alignment configuration");
    addText(parser,
            "The alignment configuration is a string of four characters, each being either t or f. All "
            "combinations are allowed. The meaning is as follows.");
    addListItem(parser, "tfff", "First row initialized with 0s.");
    addListItem(parser, "ftff", "First column initialized with 0s.");
    addListItem(parser, "fftf", "Search last column for maximum.");
    addListItem(parser, "ffft", "Search last row for maximum.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.inputFile, parser, "seq");
    getOptionValue(options.outputFile, parser, "outfile");

    getOptionValue(options.alphabet, parser, "alphabet");
    getOptionValue(options.method, parser, "method");

    //getOptionValue(options.format, parser, "format");
    getOptionValue(options.gop, parser, "gop");
    getOptionValue(options.gex, parser, "gex");
    getOptionValue(options.matrix, parser, "matrix");
    getOptionValue(options.msc, parser, "msc");
    getOptionValue(options.mmsc, parser, "mmsc");
    getOptionValue(options.low, parser, "low");
    getOptionValue(options.high, parser, "high");
    getOptionValue(options.config, parser, "config");

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, const char* argv[])
{
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Defer to the library code pair_align_lib.a.
    // For every major alignment configuration there is a single translation unit
    // breaking down the compiled source code to single sources. This reduces the
    // memory needed for the compilation to some MB. Note before it were several GB.
    // Also the library can built in parallel which also speeds up the compilation
    // process.

    pairAlignMain(options);

    return res;
}
