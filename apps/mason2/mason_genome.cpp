// ==========================================================================
//                         Mason - A Read Simulator
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
// Simulate a random genome.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>

#include "simulate_genome.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MasonGenomeOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct MasonGenomeOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The output file name.
    seqan::CharString outputFilename;

    // Lengths of the contigs.
    seqan::String<int> contigLengths;

    // The seed to use for the RNG.
    int seed;

    MasonGenomeOptions() : verbosity(1), seed(0)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonGenomeOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_genome");
    // Set short description, version, and date.
    setShortDescription(parser, "Random Genome Simulation");
    setDateAndVersion(parser);
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] [\\fB-l\\fP \\fILEN\\fP]+ \\fB-o\\fP \\fIOUT.fa\\fP");
    addDescription(parser,
                   "Simulate a random genome to the output file.  For each \\fB-l\\fP/\\fB--contig-length\\fP "
                   "entry, a contig with the given length will be simulated.");

    // We require one argument.
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addSection(parser, "Simulation Configuration");
    addOption(parser, seqan::ArgParseOption("l", "contig-length",
                                            "Length of the contig to simulate. Give one \\fB-l\\fP "
                                            "value for each contig to simulate.",
                                            seqan::ArgParseOption::INTEGER, "LENGTH", true));
    setMinValue(parser, "contig-length", "1");
    setRequired(parser, "contig-length");

    addOption(parser, seqan::ArgParseOption("s", "seed", "The seed to use for the random number generator.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "seed", 0);

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o", "out-file", "Output file.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "FILE"));
    setValidValues(parser, "out-file", seqan::SeqFileOut::getFileExtensions());
    setRequired(parser, "out-file");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBmason_genome\\fP \\fB-l\\fP 1000 \\fB-l\\fP 4000 \\fB-o\\fP \\fIgenome.fa\\fP",
                "Simulate a genome with two contigs of lengths 1000 and 4000 and write it to genome.fa.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.outputFilename, parser, "out-file");
    getOptionValue(options.seed, parser, "seed");

    for (unsigned i = 0; i < getOptionValueCount(parser, "contig-length"); ++i)
    {
        int len = 0;
        getOptionValue(len, parser, "contig-length", i);
        appendValue(options.contigLengths, len);
    }

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    MasonGenomeOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "MASON GENOME SIMULATOR\n"
              << "======================\n\n";

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY  \t" << options.verbosity << '\n'
                  << "\n"
                  << "SEED       \t" << options.seed << '\n'
                  << "\n"
                  << "OUTPUT FILE\t" << options.outputFilename << "\n"
                  << "CONTIG LENS\t";
        for (unsigned i = 0; i < length(options.contigLengths); ++i)
        {
            if (i > 0)
                std::cout << ", ";
            std::cout << options.contigLengths[i];
        };
        std::cout << "\n\n";
    }

    // Perform genome simulation.
    std::cout << "__SIMULATING GENOME__________________________________________________________\n"
              << "\n";

    MasonSimulateGenomeOptions simOptions;
    simOptions.contigLengths = options.contigLengths;
    simOptions.seed = options.seed;
    if (simulateGenome(toCString(options.outputFilename), simOptions) != 0)
        return 1;

    std::cerr << "\nDone.\n";
    return 0;
}
