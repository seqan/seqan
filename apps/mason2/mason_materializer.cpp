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
// Apply variants from a VCF file to a genomic sequence.
//
// The variants must be equivalent to the variants written by mason_variator.
// See the documentation of mason_materializer and mason_variator for details
// on this.
// ==========================================================================

// Note: We treat all given variants as phased.

#include <fstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

#include "vcf_materialization.h"
#include "mason_options.h"
#include "mason_types.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class MasonMaterializerApp
// --------------------------------------------------------------------------

class MasonMaterializerApp
{
public:
    // The configuration to use.
    MasonMaterializerOptions const & options;

    // The random number generation.
    TRng rng, methRng;

    // Materialization of VCF.
    VcfMaterializer vcfMat;

    // Output sequence stream.
    seqan::SeqFileOut outStream;
    // Output breakpoints file.
    std::fstream breakpointsOut;
    // Input and output for methylation.
    seqan::FaiIndex methFaiIndex;
    seqan::SeqFileOut outMethLevelStream;

    MasonMaterializerApp(MasonMaterializerOptions const & _options) :
            options(_options), rng(options.seed), methRng(options.methSeed),
            vcfMat(rng,
                   toCString(options.matOptions.fastaFileName),
                   toCString(options.matOptions.vcfFileName),
                   toCString(options.methFastaInFile),
                   &options.methOptions)
    {}

    int run()
    {
        // Intialization
        std::cerr << "__INITIALIZATION_____________________________________________________________\n"
                  << "\n";

        std::cerr << "Opening files...";
        try
        {
            vcfMat.init();

            if (!open(outStream, toCString(options.outputFileName)))
                throw MasonIOException("Could not open output file.");

            // Open output breakpoints TSV file.
            if (!empty(options.outputBreakpointFile))
            {
                breakpointsOut.open(toCString(options.outputBreakpointFile), std::ios::binary | std::ios::out);
                if (!breakpointsOut.good())
                    throw MasonIOException("Could not open breakpoints output file.");
                breakpointsOut << "#ref\tid\tpos\n";
            }

            if (options.methOptions.simulateMethylationLevels)
            {
                if (!open(outMethLevelStream, toCString(options.methFastaOutFile)))
                    throw MasonIOException("Could not open methylation output file.");
            }
        }
        catch (MasonIOException & e)
        {
            std::cerr << "\nERROR: " << e.what() << "\n";
            return 1;
        }
        std::cerr << " OK\n";

        // Perform genome simulation.
        std::cerr << "\n__MATERIALIZING______________________________________________________________\n"
                  << "\n";

        // The identifiers of the just materialized data.
        int rID = 0, hID = 0;
        seqan::Dna5String seq;
        std::cerr << "Materializing...";
        MethylationLevels levels;
        std::vector<SmallVarInfo> varInfos;  // small variants for counting in read alignments
        std::vector<std::pair<int, int> > breakpoints;
        if (options.methOptions.simulateMethylationLevels)  // methylation level simulation
            while (vcfMat.materializeNext(seq, levels, varInfos, breakpoints, rID, hID))
            {
                std::stringstream ssName;
                ssName << contigNames(context(vcfMat.vcfFileIn))[rID] << options.haplotypeNameSep << (hID + 1);
                std::cerr << " " << ssName.str();

                writeRecord(outStream, ssName.str(), seq);

                if (!empty(options.outputBreakpointFile))
                    for (std::vector<std::pair<int, int> >::const_iterator it = breakpoints.begin(); it != breakpoints.end(); ++it)
                        breakpointsOut << ssName.str() << "\t" << vcfMat.contigVariants.getVariantName(it->second)
                                       << "\t" << (it->first + 1) << "\n";

                std::stringstream ssTop;
                ssTop << ssName.str() << "/TOP";
                writeRecord(outMethLevelStream, ssTop.str(), levels.forward);
                std::stringstream ssBottom;
                ssBottom << ssName.str() << "/BOT";
                writeRecord(outMethLevelStream, ssBottom.str(), levels.reverse);
            }
        else  // NO methylation level simulation
            while (vcfMat.materializeNext(seq, varInfos, breakpoints, rID, hID))
            {
                std::stringstream ssName;
                ssName << contigNames(context(vcfMat.vcfFileIn))[rID] << options.haplotypeNameSep << (hID + 1);
                std::cerr << " " << ssName.str();

                writeRecord(outStream, ssName.str(), seq);

                if (!empty(options.outputBreakpointFile))
                    for (std::vector<std::pair<int, int> >::const_iterator it = breakpoints.begin(); it != breakpoints.end(); ++it)
                        breakpointsOut << ssName.str() << "\t" << vcfMat.contigVariants.getVariantName(it->second)
                                       << "\t" << (it->first + 1) << "\n";
            }
        std::cerr << " DONE\n";

        std::cerr << "\nDone materializing VCF file.\n";

        return 0;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(MasonMaterializerOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("mason_materializer");
    // Set short description, version, and date.
    setShortDescription(parser, "VCF Materialization");
    setDateAndVersion(parser);
    setCategory(parser, "Simulators");

    // Define usage line and long description.
    addUsageLine(parser,
                 "[OPTIONS] \\fB-ir\\fP \\fIIN.fa\\fP \\fB-iv\\fP \\fIIN.vcf\\fP \\fB-o\\fP \\fIOUT.fa\\fP ");
    addDescription(parser,
                   "Apply variants from \\fIIN.vcf\\fP to \\fIIN.fa\\fP and write the results to \\fIout.fa\\fP.");

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
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    MasonMaterializerOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cerr << "MASON VARIANT MATERIALIZER\n"
              << "==========================\n\n";

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cerr << "__OPTIONS____________________________________________________________________\n"
                  << "\n";
        options.print(std::cerr);
    }

    MasonMaterializerApp app(options);
    return app.run();
}
