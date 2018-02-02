// ==========================================================================
//                              bisar
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

//#define POST_PRO_PROFILE

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/bam_io.h>
#include <seqan/score.h>

#include "bisar_score_data.h"
#include "bisar_score.h"
#include "bisar_base.h"
#include "bisar.h"


using namespace seqan;

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    CharString readFileName;
    CharString readFileName2;
    CharString samFileName;
    CharString refFileName;
    CharString outputFileName;
    int intervalOffset;
    double minMapq;
    double max4Error;        // max. allowed real error rate
    double max3Error;        // max. allowed error rate in 3 letter alphabet (corresponding to mapper settings)
    double maxScore;
    unsigned maxBasePenalty;    // limit the penalty for a single base

    int minScore;
    bool outputSingleMates;

    double scoreMatch;
    double scoreMismatch;

    bool simpleScore;
    bool nonSimpleSubstErrors;
    bool nonSimpleInsErrors;
    bool nonSimpleDelErrors;
    double delErrorRate;

    double lambda;
    double gapOpenScore;
    double gapExtendScore;
    double scalingFactorDelErrors;
    double scalingFactorInsErrors;

    double bsConversionRate;
    double globalMethRate;
    double seqIdentity;         // Used for substitution matrix construction [0.0-1.0]
    double refNRate;            // Used for substitution matrix construction
    double pseudoMatchScale;

    AppOptions() :
        verbosity(1),
        intervalOffset(3),
        minMapq(1),
        max4Error(4),
        max3Error(3),
        maxScore(1000000),  // TODO: what would be reasonable?
        maxBasePenalty(-3),  // scaled to single penalties
        minScore(0),
        outputSingleMates(true),    // Output also read whose mate didn't map & if no match mate pair found, output mates single
        scoreMatch(10.0),           // at the moment only used for pseudoWorstScore
        scoreMismatch(0.01),
        simpleScore(true),
        nonSimpleSubstErrors(false),
        nonSimpleInsErrors(false),
        nonSimpleDelErrors(false),
        delErrorRate(0.001),
        lambda(1.0),
        gapOpenScore(-4.5),
        gapExtendScore(-2.0),
        scalingFactorDelErrors(5.0),
        scalingFactorInsErrors(5.0),
        bsConversionRate(0.99),
        globalMethRate(0.5),
        seqIdentity(0.9),
        refNRate(0.01),
        pseudoMatchScale(0.9)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("bisar");
    // Set short description, version, and date.
    setShortDescription(parser, "Pairwise four-letter realignment computation for bisulfite reads");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setCategory(parser, "BS-Seq Analysis");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIALIGNMENT FILE\\fP> <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP>");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIALIGNMENT FILE\\fP> <\\fIGENOME FILE\\fP> <\\fIPE-READS FILE1\\fP> <\\fIPE-READS FILE2\\fP>");
    addDescription(parser, "This program reads three-letter mappings of bisulfite reads and computes local pairwise four-letter realignments using an advanced statistical alignment model.");

    // We require ... arguments.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "ALIGNMENTS"));
    setHelpText(parser, 0, "SAM input file containing three-letter read alignments (must be sorted by query names).");
    setValidValues(parser, 0, BamFileIn::getFileExtensions());
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "GENOME"));
    setHelpText(parser, 1, "A reference genome file.");
    setValidValues(parser, 1, SeqFileIn::getFileExtensions());
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READS", true));
    setHelpText(parser, 2, "Either one (single-end) or two (paired-end) read files.");
    setValidValues(parser, 2, SeqFileIn::getFileExtensions());

    addSection(parser, "Options");
    addOption(parser, ArgParseOption("o", "output-file", "Mapping output file.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "output-file", BamFileOut::getFileExtensions());
    setRequired(parser, "output-file", true);

    addOption(parser, ArgParseOption("e3", "max3-error", "Max. error rate in 3-letter alphabet.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "max3-error", "0");
    setMaxValue(parser, "max3-error", "100");
    setDefaultValue(parser, "max3-error", options.max3Error);
    addOption(parser, ArgParseOption("e4", "max4-error", "Max. error rate in 4-letter alphabet.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "max4-error", "0");
    setMaxValue(parser, "max4-error", "100");
    setDefaultValue(parser, "max4-error", options.max4Error);
    addOption(parser, ArgParseOption("mq", "min-mapq", "Min required mapping quality.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "min-mapq", "0");
    setDefaultValue(parser, "min-mapq", options.minMapq);
    addOption(parser, ArgParseOption("ns", "non-simple", "Use non-uniform SNP distributions."));
    hideOption(parser, "ns");
    addOption(parser, ArgParseOption("nse", "ns-subst-errors", "Use empirical substitution error frequencies of Illumina sequencing data for alignment scoring scheme (corresponding to Dohm et al. 2008)."));
    addOption(parser, ArgParseOption("nsi", "ns-ins-errors", "Use empirical insertion error frequencies of Illumina sequencing data for alignment scoring scheme (corresponding to Minoche et al. 2011)."));
    addOption(parser, ArgParseOption("nsd", "ns-del-errors", "Use empirical deletion error frequencies of Illumina sequencing data for alignment scoring scheme (corresponding to Minoche et al. 2011)."));
    addOption(parser, ArgParseOption("der", "del-error-rate", "Deletion error rate.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "del-error-rate", "0");
    setMaxValue(parser, "del-error-rate", "1");
    setDefaultValue(parser, "del-error-rate", options.delErrorRate);
    addOption(parser, ArgParseOption("gas", "gap-open-score", "Gap open score (original, must be proportional to mismatch scores).", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "gap-open-score", options.gapOpenScore);
    addOption(parser, ArgParseOption("ges", "gap-extend-score", "Gap extend score.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "gap-extend-score", options.gapExtendScore);
    addOption(parser, ArgParseOption("bsc", "bs-conversion-rate", "Bisulfite conversion rate.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "bs-conversion-rate", "0");
    setMaxValue(parser, "bs-conversion-rate", "1");
    setDefaultValue(parser, "bs-conversion-rate", options.bsConversionRate);
    addOption(parser, ArgParseOption("gmr", "global-meth-rate", "Global methylation rate.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "global-meth-rate", "0");
    setMaxValue(parser, "global-meth-rate", "1");
    setDefaultValue(parser, "global-meth-rate", options.globalMethRate);
    addOption(parser, ArgParseOption("i", "seq-identity", "Sequence identity used for substitution matrix.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "seq-identity", "0");
    setMaxValue(parser, "seq-identity", "1");
    setDefaultValue(parser, "seq-identity", options.seqIdentity);
    addOption(parser, ArgParseOption("rn", "ref-n", "Rate of Ns in reference sequence.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "ref-n", "0");
    setMaxValue(parser, "ref-n", "1");
    setDefaultValue(parser, "ref-n", options.refNRate);
    addOption(parser, ArgParseOption("pms", "pseudo-match-scale", "Scaling for pseudo match score. ", ArgParseArgument::DOUBLE));
    setMinValue(parser, "pseudo-match-scale", "0");
    setDefaultValue(parser, "pseudo-match-scale", options.pseudoMatchScale);
    hideOption(parser, "pms");

    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBbisar\\fP \\fB-e3\\fP \\fB4\\fP \\fB-e4\\fP \\fB5\\fP \\fB-o\\fP \\fBmapped_reads_verified.sam\\fP \\fBmapped_reads.sam\\fP \\fBgenome.fa\\fP \\fBreads.fastq\\fP",
                "Compute realignments for all reads with up to 4% errors in their three-letter alignment, while allowing up to 5% errors in four-letter alignments.");
    addListItem(parser, "\\fBbisar\\fP \\fB-e3\\fP \\fB4\\fP \\fB-e4\\fP \\fB5\\fP \\fB-o\\fP \\fBmapped_reads_verified.sam\\fP \\fBmapped_reads.sam\\fP \\fBgenome.fa\\fP \\fBreads_L.fastq\\fP \\fBreads_R.fastq \\fP",
                "Compute realignments for paired-end reads.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.samFileName, parser, 0);
    getArgumentValue(options.refFileName, parser, 1);

    if (1 == getArgumentValueCount(parser, 2))
        getArgumentValue(options.readFileName, parser, 2, 0);
    else if (2 == getArgumentValueCount(parser, 2))
    {
        getArgumentValue(options.readFileName, parser, 2, 0);
        getArgumentValue(options.readFileName2, parser, 2, 1);
    }
    else
    {
        std::cerr << "ERROR: " << getArgumentValueCount(parser, 2) << " read files specified (must be one or two)." << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(options.outputFileName, parser, "output-file");
    getOptionValue(options.max3Error, parser, "max3-error");
    getOptionValue(options.max4Error, parser, "max4-error");
    getOptionValue(options.minMapq, parser, "min-mapq");
    options.nonSimpleSubstErrors = isSet(parser, "ns-subst-errors");
    options.nonSimpleInsErrors = isSet(parser, "ns-ins-errors");
    options.nonSimpleDelErrors = isSet(parser, "ns-del-errors");
    getOptionValue(options.delErrorRate, parser, "del-error-rate");
    getOptionValue(options.gapOpenScore, parser, "gap-open-score");
    getOptionValue(options.gapExtendScore, parser, "gap-extend-score");
    getOptionValue(options.bsConversionRate, parser, "bs-conversion-rate");
    getOptionValue(options.globalMethRate, parser, "global-meth-rate");
    getOptionValue(options.seqIdentity, parser, "seq-identity");
    getOptionValue(options.refNRate, parser, "ref-n");
    getOptionValue(options.pseudoMatchScale, parser, "pseudo-match-scale");

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

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

#ifdef POST_PRO_PROFILE
    double timeStamp = sysTime();
#endif
    if (!options.simpleScore)
        postProcessMain(options, BsNonSimple());
    else
        postProcessMain(options, BsSimple());

#ifdef POST_PRO_PROFILE
    Times::instance().time_all = sysTime() - timeStamp;
    std::cout << "  Time needed for all: " << Times::instance().time_all/60.0 << "min" << std::endl;
    std::cout << "  Time needed for globalAlignment: " << Times::instance().time_globalAlignment/60.0 << "min" << std::endl;
    std::cout << "  Time needed for writeBsAlignment: " << Times::instance().time_writeBsAlignment/60.0 << "min" << std::endl;
#endif

    return 0;
}
