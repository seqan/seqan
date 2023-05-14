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

seqan2::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan2::ArgumentParser parser("pair_align");
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
    addOption(parser, seqan2::ArgParseOption("s", "seq", "FASTA file with two sequences.", seqan2::ArgParseOption::INPUT_FILE, "IN"));
    setRequired(parser, "seq");
    setValidValues(parser, "seq", getFileExtensions(Fasta()));
    addOption(parser, seqan2::ArgParseOption("a", "alphabet", "Sequence alphabet.", seqan2::ArgParseOption::STRING, "ALPHABET"));
    setValidValues(parser, "alphabet", "protein dna rna text");
    setDefaultValue(parser, "alphabet", "protein");
    addOption(parser, seqan2::ArgParseOption("m", "method",
                                            "DP alignment method: Needleman-Wunsch, Gotoh, Smith-Waterman, "
                                            "Longest Common Subsequence",
                                            seqan2::ArgParseOption::STRING, "METHOD"));
    setValidValues(parser, "method", "nw gotoh sw lcs");
    setDefaultValue(parser, "method", "gotoh");
    addOption(parser, seqan2::ArgParseOption("o", "outfile", "Output filename.", seqan2::ArgParseOption::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "outfile", "out.fasta");
    std::vector<std::string> outFileNames = getFileExtensions(Fasta());
    outFileNames.push_back(".msf");
    setValidValues(parser, "outfile", outFileNames);

    addSection(parser, "Scoring Options");
    addOption(parser, seqan2::ArgParseOption("g", "gop", "Gap open penalty.", seqan2::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "gop", "-11");
    addOption(parser, seqan2::ArgParseOption("e", "gex", "Gap extension penalty.", seqan2::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "gex", "-1");
    addOption(parser, seqan2::ArgParseOption("ma", "matrix", "Score matrix.", seqan2::ArgParseOption::STRING, "MATRIX_FILE"));
    addOption(parser, seqan2::ArgParseOption("ms", "msc", "Match score.", seqan2::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "msc", "5");
    addOption(parser, seqan2::ArgParseOption("mm", "mmsc", "Mismatch penalty.", seqan2::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "mmsc", "-4");

    addSection(parser, "Banded Alignment Options");
    addOption(parser, seqan2::ArgParseOption("lo", "low", "Lower diagonal.", seqan2::ArgParseOption::INTEGER, "INT"));
    addOption(parser, seqan2::ArgParseOption("hi", "high", "Upper diagonal.", seqan2::ArgParseOption::INTEGER, "INT"));

    addSection(parser, "DP Matrix Configuration Options");
    addOption(parser, seqan2::ArgParseOption("c", "config", "Alignment configuration.", seqan2::ArgParseOption::STRING, "CONF"));
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
    seqan2::ArgumentParser::ParseResult res = seqan2::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan2::ArgumentParser::PARSE_OK)
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

    return seqan2::ArgumentParser::PARSE_OK;
}

int main(int argc, const char* argv[])
{
    // Parse the command line.
    Options options;
    seqan2::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan2::ArgumentParser::PARSE_OK)
        return res == seqan2::ArgumentParser::PARSE_ERROR;

    // Defer to the library code pair_align_lib.a.
    // For every major alignment configuration there is a single translation unit
    // breaking down the compiled source code to single sources. This reduces the
    // memory needed for the compilation to some MB. Note before it were several GB.
    // Also the library can built in parallel which also speeds up the compilation
    // process.

    pairAlignMain(options);

    return res;
}
