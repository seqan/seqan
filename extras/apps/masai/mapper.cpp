// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the masai_mapper application.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "options.h"
#include "mapper.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options : public MasaiOptions
{
    CharString  genomeFile;
    CharString  genomeIndexFile;
    IndexType   genomeIndexType;
    CharString  readsFile;

    CharString  mappedReadsFile;
    OutputFormat outputFormat;
    bool        outputCigar;

    MappingMode mappingMode;
    unsigned    errorsPerRead;
    unsigned    errorsLossy;
    unsigned    seedLength;
    bool        mismatchesOnly;

    bool        verifyHits;
    bool        dumpResults;
    bool        multipleBacktracking;

    Options() :
        MasaiOptions(),
        genomeIndexType(INDEX_ESA),
        outputFormat(RAW),
        outputCigar(true),
        mappingMode(ANY_BEST),
        errorsPerRead(5),
        errorsLossy(0),
        seedLength(33),
        mismatchesOnly(false),
        verifyHits(true),
        dumpResults(true),
        multipleBacktracking(true)
    {}
};

// ============================================================================

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "masai_mapper");
    setShortDescription(parser, "Masai Mapper");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta");
    setValidValues(parser, 1, "fastq");

    addSection(parser, "Mapping Options");

    addOption(parser, ArgParseOption("mm", "mapping-mode", "Select mapping mode.", ArgParseOption::STRING));
    setValidValues(parser, "mapping-mode", options.mappingModeList);
    setDefaultValue(parser, "mapping-mode", options.mappingModeList[options.mappingMode]);

    addOption(parser, ArgParseOption("e", "errors", "Maximum number of errors per read.", ArgParseOption::INTEGER));
    setMinValue(parser, "errors", "0");
    setMaxValue(parser, "errors", "32");
    setDefaultValue(parser, "errors", options.errorsPerRead);

    addOption(parser, ArgParseOption("el", "errors-lossy",
                                     "Maximum number of errors per read to report. For any-best mode only.",
                                     ArgParseOption::INTEGER));
    setMinValue(parser, "errors-lossy", "0");
    setMaxValue(parser, "errors-lossy", "32");
    setDefaultValue(parser, "errors-lossy", options.errorsLossy);

    addOption(parser, ArgParseOption("sl", "seed-length", "Minimum seed length.", ArgParseOption::INTEGER));
    setMinValue(parser, "seed-length", "10");
    setMaxValue(parser, "seed-length", "100");
    setDefaultValue(parser, "seed-length", options.seedLength);

    addOption(parser, ArgParseOption("ng", "no-gaps", "Do not align reads with gaps."));


    addSection(parser, "Genome Index Options");

    addOption(parser, ArgParseOption("x", "index", "Select genome index type.", ArgParseOption::STRING));
    setValidValues(parser, "index", options.indexTypeList);
    setDefaultValue(parser, "index", options.indexTypeList[options.genomeIndexType]);

    addOption(parser, ArgParseOption("xp", "index-prefix", "Specify genome index prefix name.", ArgParseOption::STRING));


    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-file", "Specify output file. Required.", ArgParseOption::OUTPUTFILE));
    setRequired(parser, "output-file");

    addOption(parser, ArgParseOption("of", "output-format", "Select output format.", ArgParseOption::STRING));
    setValidValues(parser, "output-format", options.outputFormatList);
    setDefaultValue(parser, "output-format", options.outputFormatList[options.outputFormat]);


    addSection(parser, "Debug Options");

    addOption(parser, ArgParseOption("nv", "no-verify", "Do not verify seed hits."));
    addOption(parser, ArgParseOption("nd", "no-dump", "Do not dump results."));
    addOption(parser, ArgParseOption("nm", "no-multiple", "Disable multiple backtracking."));
}

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Parse genome input file.
    getArgumentValue(options.genomeFile, parser, 0);

    // Parse reads input file.
    getArgumentValue(options.readsFile, parser, 1);

    // Parse mapping mode.
    getOptionValue(options.mappingMode, parser, "mapping-mode", options.mappingModeList);

    // Parse mapping options.
    getOptionValue(options.errorsPerRead, parser, "errors");
    getOptionValue(options.errorsLossy, parser, "errors-lossy");
    options.errorsLossy = std::max(options.errorsLossy, options.errorsPerRead);
    getOptionValue(options.seedLength, parser, "seed-length");
    options.mismatchesOnly = isSet(parser, "no-gaps");

    // Parse genome index prefix.
    getOptionValue(options.genomeIndexFile, parser, "index-prefix");
    if (!isSet(parser, "index-prefix"))
    {
        // TODO(esiragusa): Trim extension .fasta from genomeIndexFile.
        options.genomeIndexFile = options.genomeFile;
    }

    // Parse genome index type.
    getOptionValue(options.genomeIndexType, parser, "index", options.indexTypeList);

    // Parse output file.
    getOptionValue(options.mappedReadsFile, parser, "output-file");
    if (!isSet(parser, "output-file"))
    {
        options.mappedReadsFile = options.readsFile;
        // TODO(esiragusa): Guess output file extension.
        append(options.mappedReadsFile, ".out");
    }

    // Parse output format.
    getOptionValue(options.outputFormat, parser, "output-format", options.outputFormatList);

    // Parse debug options.
    options.verifyHits = !isSet(parser, "no-verify");
    options.dumpResults = !isSet(parser, "no-dump");
    options.multipleBacktracking = !isSet(parser, "no-multiple");

    return seqan::ArgumentParser::PARSE_OK;
}

// ============================================================================

template <typename TStrategy, typename TDistance, typename TFormat, typename TBacktracking, typename TIndex>
int runMapper(Options & options)
{
    // TODO(esiragusa): Use typename TGenomeIndex instead of TSpec.
    typedef Mapper<TIndex>    TMapper;

    // TODO(esiragusa): Remove writeCigar from Mapper members.
    TMapper mapper(options.seedLength, options.outputCigar, options.verifyHits, options.dumpResults);

    double start, finish;

    // Loading genome.
    std::cout << "Loading genome:\t\t\t" << std::flush;
    start = sysTime();
    if (!loadGenome(mapper.indexer, options.genomeFile))
    {
        std::cerr << "Error while loading genome" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;
//	std::cout << "Contigs count:\t\t\t" << mapper.indexer.contigsCount << std::endl;

    // Loading genome index.
    std::cout << "Loading genome index:\t\t" << std::flush;
    start = sysTime();
    if (!loadGenomeIndex(mapper.indexer, options.genomeIndexFile))
    {
        std::cout << "Error while loading genome index" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    // Loading reads.
    std::cout << "Loading reads:\t\t\t" << std::flush;
    start = sysTime();
    if (!loadReads(mapper, options.readsFile))
    {
        std::cerr << "Error while loading reads" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;
    std::cout << "Reads count:\t\t\t" << mapper.readsCount << std::endl;

    // Mapping reads.
    start = sysTime();
    mapReads(mapper, options.mappedReadsFile, options.errorsPerRead, options.errorsLossy,
             TDistance(), TStrategy(), TBacktracking(), TFormat());
    finish = sysTime();
    std::cout << "Mapping time:\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    return 0;
}

// ============================================================================

template <typename TStrategy, typename TDistance, typename TFormat, typename TBacktracking>
int configureGenomeIndex(Options & options)
{
    switch (options.genomeIndexType)
    {
    case Options::INDEX_ESA:
        return runMapper<TStrategy, TDistance, TFormat, TBacktracking, TGenomeEsa>(options);

    case Options::INDEX_SA:
        return runMapper<TStrategy, TDistance, TFormat, TBacktracking, TGenomeSa>(options);

//    case Options::INDEX_QGRAM:
//        return runMapper<TStrategy, TDistance, TFormat, TBacktracking, TGenomeQGram>(options);

        case Options::INDEX_FM:
            return runMapper<TStrategy, TDistance, TFormat, TBacktracking, TGenomeFM>(options);

    default:
        return 1;
    }
}

template <typename TStrategy, typename TDistance, typename TFormat>
int configureBacktracking(Options & options)
{
    if (options.multipleBacktracking)
        return configureGenomeIndex<TStrategy, TDistance, TFormat, MultipleBacktracking>(options);
    else
        return configureGenomeIndex<TStrategy, TDistance, TFormat, SingleBacktracking>(options);
}

template <typename TStrategy, typename TDistance>
int configureOutputFormat(Options & options)
{
    switch (options.outputFormat)
    {
    case Options::RAW:
        return configureBacktracking<TStrategy, TDistance, Raw>(options);

    case Options::SAM:
        return configureBacktracking<TStrategy, TDistance, Sam>(options);

    case Options::SAM_NO_CIGAR:
        options.outputCigar = false;
        return configureBacktracking<TStrategy, TDistance, Sam>(options);

    default:
        return 1;
    }
}

template <typename TStrategy>
int configureDistance(Options & options)
{
    if (options.mismatchesOnly)
        return configureOutputFormat<TStrategy, HammingDistance>(options);
    else
        return configureOutputFormat<TStrategy, EditDistance>(options);
}

int mainWithOptions(Options & options)
{
    switch (options.mappingMode)
    {
    case Options::ANY_BEST:
        return configureDistance<AnyBest>(options);

    case Options::ALL_BEST:
        return configureDistance<AllBest>(options);

    case Options::ALL:
        return configureDistance<All>(options);

    default:
        return 1;
    }
}

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
        return mainWithOptions(options);
    else
        return res;
}
