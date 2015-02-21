// ==========================================================================
//                               FX Tools
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Computes density of variants given a VCF file, window-based.
// ==========================================================================

#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/vcf_io.h>
#include <seqan/seq_io.h>

// --------------------------------------------------------------------------
// Class FxVariantDensityOptions
// --------------------------------------------------------------------------

struct FxVariantDensityOptions
{
    // Verbosity level.  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;
    
    // Path to SAM file.
    seqan::CharString inVcfPath;

    // Path to output file.
    seqan::CharString outPath;

    // Window size to use for computation.
    __int32 windowSize;

    FxVariantDensityOptions() : verbosity(1), windowSize(10*1000)
    {}
};


// --------------------------------------------------------------------------
// Class BinData
// --------------------------------------------------------------------------

struct BinData
{
    // Number of reads whose alignment starts here.
    unsigned coverage;
    // Length of underlying sequence.
    unsigned length;
    
    BinData() : coverage(0), length(0)
    {}
};

// --------------------------------------------------------------------------
// Function parseArgs()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseArgs(FxVariantDensityOptions & options,
          int argc,
          char const ** argv)
{
    seqan::ArgumentParser parser("fx_variant_density");
    setShortDescription(parser, "Variant Density Computation.");
    setCategory(parser, "Utilities");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fB-o\\fP \\fIOUT.bed\\fP "
                 "\\fB-i\\fP \\fIVARIANTS.vcf\\fP");
    addDescription(parser, "Compute variant density for a genome.");

    // Input files: vcf.
    addOption(parser, seqan::ArgParseOption("i", "in-variants", "Path to the variant file to analyze.",
                                            seqan::ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "in-variants", seqan::VcfFileIn::getFileExtensions());
    setRequired(parser, "in-variants");

    // TODO(holtgrew): I want a custom help text!
    // addOption(parser, seqan::ArgParseOption("h", "help", "This helpful screen."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose, log to STDERR."));
    hideOption(parser, "verbose");
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose, log to STDERR."));
    hideOption(parser, "very-verbose");

    addSection(parser, "Main Options");
    addOption(parser, seqan::ArgParseOption("w", "window-size", "Set the size of the non-overlapping windows in base pairs.", seqan::ArgParseArgument::INTEGER, "NUM"));
    setDefaultValue(parser, "window-size", options.windowSize);

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o", "out-path", "Path to the resulting file.  If omitted, result is printed to stdout.", seqan::ArgParseArgument::OUTPUT_FILE, "BED"));
    setRequired(parser, "out-path");
    setValidValues(parser, "out-path", "bed");

    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
    {
        getOptionValue(options.inVcfPath, parser, "in-variants");
        getOptionValue(options.outPath, parser, "out-path");
        getOptionValue(options.windowSize, parser, "window-size");

        if (isSet(parser, "verbose"))
            options.verbosity = 2;
        if (isSet(parser, "very-verbose"))
            options.verbosity = 3;
    }

    return res;
}

// ---------------------------------------------------------------------------
// Function main()
// ---------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    double startTime = seqan::sysTime();
    
    // -----------------------------------------------------------------------
    // Parse command line.
    // -----------------------------------------------------------------------
    FxVariantDensityOptions options;
    seqan::ArgumentParser::ParseResult res = parseArgs(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // -----------------------------------------------------------------------
    // Show options.
    // -----------------------------------------------------------------------
    if (options.verbosity >= 1)
    {
        std::cerr << "____OPTIONS___________________________________________________________________\n"
                  << "\n"
                  << "VERBOSITY    " << options.verbosity << "\n"
                  << "SAM/BAM      " << options.inVcfPath << "\n"
                  << "OUT          " << options.outPath << "\n"
                  << "WINDOW SIZE  " << options.windowSize << "\n";
    }

    // -----------------------------------------------------------------------
    // Compute Coverage
    // -----------------------------------------------------------------------

    std::cerr << "\n"
              << "___DENSITY COMPUTATION________________________________________________________\n"
              << "\n"
              << "Computing Density...";

    seqan::VcfFileIn vcfFile;
    if (!open(vcfFile, toCString(options.inVcfPath)))
    {
        std::cerr << "Could not open " << options.inVcfPath << "!\n";
        return 1;
    }
    
    seqan::String<seqan::String<BinData> > bins;

    seqan::VcfHeader header;
    readHeader(header, vcfFile);

    seqan::VcfRecord record;
    while (!atEnd(vcfFile))
    {
        readRecord(record, vcfFile);

        unsigned binNo = record.beginPos / options.windowSize;
        if (length(bins) <= (unsigned)record.rID)
            resize(bins, record.rID + 1);
        if (length(bins[record.rID]) <= binNo)
            resize(bins[record.rID], binNo + 1);
        bins[record.rID][binNo].coverage++;
    }

    std::cerr << "DONE\n";

    // -----------------------------------------------------------------------
    // Write Output
    // -----------------------------------------------------------------------

    std::ostream * out = &std::cout;
    std::ofstream outFile;
    if (options.outPath != "-")
    {
        outFile.open(toCString(options.outPath), std::ios::binary | std::ios::out);
        if (!outFile.good())
        {
            std::cerr << "ERROR: Could not open output file " << options.outPath << "!\n";
            return 1;
        }
        out = &outFile;
    }

    for (unsigned i = 0, globalBin = 0; i < length(bins); ++i)
    {
        for (unsigned refBin = 0; refBin < length(bins[i]); ++refBin, ++globalBin)
        {
            (*out) << contigNames(context(vcfFile))[i] << '\t'
                   << refBin * options.windowSize << '\t'
                   << (refBin + 1) * options.windowSize << '\t'
                   << bins[i][refBin].coverage << '\n';
        }
    }

    if (options.verbosity >= 2)
        std::cerr << "Took " << (seqan::sysTime() - startTime) << " s\n";

    return 0;
}
