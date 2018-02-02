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
// Simulate sequencing process on fragments.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/arg_parse.h>

#include "mason_types.h"
#include "mason_options.h"
#include "sequencing.h"

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
parseCommandLine(MasonFragmentSequencingOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_frag_sequencing");
    // Set short description, version, and date.
    setShortDescription(parser, "Fragment Sequencing Simulation");
    setDateAndVersion(parser);
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-i\\fP \\fIIN.fa\\fP \\fB-o\\fP \\fIOUT.{fa,fq}\\fP "
                 "[\\fB-or\\fP \\fIOUT2.{fa,fq}\\fP]");
    addDescription(parser,
                   "Given a FASTA file with fragments, simulate sequencing thereof.");
    addDescription(parser,
                   "This program is a more lightweight version of mason_sequencing without support for the "
                   "application of VCF and fragment sampling.  Output of SAM is also not available.  However, "
                   "it uses the same code for the simulation of the reads as the more powerful mason_simulator.");
    addDescription(parser,
                   "You can use mason_frag_sequencing if you want to implement you rown fragmentation behaviour, e.g. "
                   "if you have implemented your own bias models.");

    // Add option and text sections.
    options.addOptions(parser);
    options.addTextSections(parser);

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    options.getOptionValues(parser);

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function trimAfterSpace()
// --------------------------------------------------------------------------

void trimAfterSpace(seqan::CharString & s)
{
    unsigned i = 0;
    for (; i < length(s); ++i)
        if (isspace(s[i]))
            break;
    resize(s, i);
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    MasonFragmentSequencingOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cerr << "MASON SEQUENCING SIMULATOR\n"
              << "==========================\n\n";

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cerr << "__OPTIONS____________________________________________________________________\n"
                  << '\n';
        options.print(std::cerr);
    }

    std::cerr << "\n__PREPARATION________________________________________________________________\n"
              << "\n";

    // Open fragments FASTA file.
    std::cerr << "Opening fragments       " << options.inputFileName << " ...";
    seqan::SeqFileIn inFragments;
    if (!open(inFragments, toCString(options.inputFileName)))
    {
        std::cerr << " ERROR\n"
                  << "Could not open " << options.inputFileName << "\n";
        return 1;
    }
    std::cerr << " OK\n";

    // Open reads output file.
    std::cerr << "Opening output file (L) " << options.outFileNameLeft << " ...";
    seqan::SeqFileOut outReads;
    if (!open(outReads, toCString(options.outFileNameLeft)))
    {
        std::cerr << " ERROR\n"
                  << "Could not open " << options.outFileNameLeft << "\n";
        return 1;
    }
    std::cerr << " OK\n";

    // Open output file for the right reads.
    seqan::SeqFileOut outReadsRight;
    if (!empty(options.outFileNameRight))
    {
        std::cerr << "Opening output file (R) " << options.outFileNameRight << " ...";
        if (!open(outReadsRight, toCString(options.outFileNameRight)))
        {
            std::cerr << " ERROR\n"
                      << "Could not open " << options.outFileNameRight << "\n";
            return 1;
        }
        std::cerr << " OK\n";
    }

    // Configure output streams to write out each sequence in a single line.
    context(outReads).options.lineLength = 0;
    context(outReads).options.lineLength = 0;

    // Perform genome simulation.
    std::cerr << "\n__SIMULATING READS___________________________________________________________\n"
              << "\n"
              << "Simulating reads ...";

    TRng rng(options.seed);
    TRng ignoredMethRng(0);

    // Create sequencing simulator.
    SequencingSimulatorFactory simFactory(rng, ignoredMethRng, options.seqOptions, options.illuminaOptions,
                                          options.rocheOptions, options.sangerOptions);
    std::unique_ptr<SequencingSimulator> sim = simFactory.make();

    // Buffers for reading in fragments.
    seqan::CharString fragId;
    seqan::Dna5String fragSeq;
    // Buffers for simulated reads.
    seqan::Dna5String seqL, seqR;
    seqan::CharString qualsL, qualsR;
    // The information for storing the simulation info.
    SequencingSimulationInfo simInfoL, simInfoR;

    // We will use these string streams to generate the read identifier strings.
    std::stringstream ssL, ssR;

    for (unsigned readId = 1; !atEnd(inFragments); ++readId)
    {
        // Reset the string streams.
        ssL.str("");
        ssL.clear();
        ssR.str("");
        ssR.clear();

        // Read fragment to simulate from.
        readRecord(fragId, fragSeq, inFragments);

        // Trim fragment identifier after first whitespace.
        trimAfterSpace(fragId);

        if (empty(options.outFileNameRight))  // Single-end sequencing.
        {
            sim->simulateSingleEnd(seqL, qualsL, simInfoL, infix(fragSeq, 0, length(fragSeq)));
            ssL << options.seqOptions.readNamePrefix << readId;
            if (options.seqOptions.embedReadInfo)
            {
                ssL << ' ';
                simInfoL.serialize(ssL);
                ssL << " FRAG_ID=" << fragId;
            }
            writeRecord(outReads, ssL.str(), seqL, qualsL);
        }
        else  // Paired sequencing.
        {
            sim->simulatePairedEnd(seqL, qualsL, simInfoL, seqR, qualsR, simInfoR, infix(fragSeq, 0, length(fragSeq)));
            ssL << options.seqOptions.readNamePrefix << readId;
            ssR << options.seqOptions.readNamePrefix << readId;
            if (options.seqOptions.embedReadInfo)
            {
                ssL << ' ';
                simInfoL.serialize(ssL);
                ssL << " FRAG_ID=" << fragId;
                ssR << ' ';
                simInfoR.serialize(ssR);
                ssR << " FRAG_ID=" << fragId;
            }

            // std::cerr << seqL << "\t" << qualsL << "\n"
            //           << seqR << "\t" << qualsR << "\n\n";

            writeRecord(outReads, ssL.str(), seqL, qualsL);
            writeRecord(outReadsRight, ssR.str(), seqR, qualsR);
        }
    }

    std::cerr << " OK\n";

    std::cerr << "\nDONE.\n";
    return 0;
}
