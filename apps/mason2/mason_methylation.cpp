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
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

#include "mason_options.h"
#include "methylation_levels.h"

// ==========================================================================
// Classes
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonMethylationOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_methylation");
    // Set short description, version, and date.
    setShortDescription(parser, "Methylation Level Simulation");
    setDateAndVersion(parser);
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser, "[OPTIONS] \\fB-i\\fP \\fIIN.fa\\fP \\fB-o\\fP \\fIOUT.fa\\fP");
    addDescription(parser, "Simulate methylation levels for \\fIIN.fa\\fP and write them to \\fIOUT.fa\\fP.");

    options.addOptions(parser);

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    options.getOptionValues(parser);
    options.methOptions.simulateMethylationLevels = true;

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    MasonMethylationOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cerr << "MASON METHYLATION SIMULATION\n"
              << "============================\n\n";

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
        options.print(std::cerr);

    std::cerr << "\n__PREPARATION_________________________________________________________________\n"
              << "\n";

    std::cerr << "Loading Reference Index " << options.fastaInFile << " ...";
    seqan::FaiIndex faiIndex;
    if (!open(faiIndex, toCString(options.fastaInFile)))
    {
        std::cerr << " FAILED (not fatal, we can just build it)\n";
        std::cerr << "Building Index        " << options.fastaInFile << ".fai ...";
        if (!build(faiIndex, toCString(options.fastaInFile)))
        {
            std::cerr << "Could not build FAI index.\n";
            return 1;
        }
        std::cerr << " OK\n";
        seqan::CharString faiPath = options.fastaInFile;
        append(faiPath, ".fai");
        std::cerr << "Reference Index       " << faiPath << " ...";
        if (!save(faiIndex, toCString(faiPath)))
        {
            std::cerr << "Could not write FAI index we just built.\n";
            return 1;
        }
        std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
    }
    else
    {
        std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
    }

    std::cerr << "Opening output File " << options.methFastaOutFile << " ...";
    seqan::SeqFileOut outStream;
    if (!open(outStream, toCString(options.methFastaOutFile)))
    {
        std::cerr << "\nERROR: Could not open output file.\n";
        return 1;
    }
    std::cerr << " OK\n";

    std::cerr << "\n__SIMULATION__________________________________________________________________\n"
              << "\n";

    TRng rng(options.seed);
    MethylationLevelSimulator methSim(rng, options.methOptions);

    MethylationLevels levels;

    seqan::Dna5String contig;
    for (unsigned i = 0; i < numSeqs(faiIndex); ++i)
    {
        levels.clear();

        std::cerr << "Simulating for " << sequenceName(faiIndex, i) << " ...";
        readSequence(contig, faiIndex, i);

        methSim.run(levels, contig);

        std::stringstream ssTop;
        ssTop << sequenceName(faiIndex, i) << "/TOP";
        writeRecord(outStream, ssTop.str().c_str(), levels.forward);

        std::stringstream ssBottom;
        ssBottom << sequenceName(faiIndex, i) << "/BOT";
        writeRecord(outStream, ssBottom.str().c_str(), levels.reverse);

        std::cerr << " OK\n";
    }
    std::cerr << "\nDone with methylation simulation.\n";

    return 0;
}
