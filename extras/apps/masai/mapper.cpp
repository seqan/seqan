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
#include "index.h"
#include "mapper.h"
#include "writer.h"

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
    int         mappingBlock;

    CharString  mappedReadsFile;
    OutputFormat outputFormat;
    bool        outputCigar;

    MappingMode mappingMode;
    unsigned    errorsPerRead;
    unsigned    seedLength;
    bool        mismatchesOnly;

    bool        noVerify;
    bool        noDump;
    bool        noMultiple;

    Options() :
        MasaiOptions(),
        genomeIndexType(INDEX_SA),
        mappingBlock(MaxValue<int>::VALUE),
        outputFormat(RAW),
        outputCigar(true),
        mappingMode(ANY_BEST),
        errorsPerRead(5),
        seedLength(33),
        mismatchesOnly(false),
        noVerify(false),
        noDump(false),
        noMultiple(false)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()                              [ArgumentParser]
// ----------------------------------------------------------------------------

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
    setValidValues(parser, 0, "fasta fa");
    setValidValues(parser, 1, "fastq fasta fa");

    addSection(parser, "Mapping Options");

    addOption(parser, ArgParseOption("mm", "mapping-mode", "Select mapping mode.", ArgParseOption::STRING));
    setValidValues(parser, "mapping-mode", options.mappingModeList);
    setDefaultValue(parser, "mapping-mode", options.mappingModeList[options.mappingMode]);

    addOption(parser, ArgParseOption("mb", "mapping-block", "Maximum number of reads to be mapped at once.", ArgParseOption::INTEGER));
    setMinValue(parser, "mapping-block", "10000");
    setDefaultValue(parser, "mapping-block", options.mappingBlock);

    addOption(parser, ArgParseOption("e", "errors", "Maximum number of errors per read.", ArgParseOption::INTEGER));
    setMinValue(parser, "errors", "0");
    setMaxValue(parser, "errors", "32");
    setDefaultValue(parser, "errors", options.errorsPerRead);

    addOption(parser, ArgParseOption("sl", "seed-length", "Minimum seed length.", ArgParseOption::INTEGER));
    setMinValue(parser, "seed-length", "10");
    setMaxValue(parser, "seed-length", "100");
    setDefaultValue(parser, "seed-length", options.seedLength);

    addOption(parser, ArgParseOption("ng", "no-gaps", "Do not align reads with gaps."));


    addSection(parser, "Genome Index Options");

    setIndexType(parser, options);
    setIndexPrefix(parser);


    addSection(parser, "Output Options");

    setOutputFile(parser);
    setOutputFormat(parser, options);


    addSection(parser, "Debug Options");

    addOption(parser, ArgParseOption("nv", "no-verify", "Do not verify seed hits."));
    addOption(parser, ArgParseOption("nd", "no-dump", "Do not dump results."));
    addOption(parser, ArgParseOption("nm", "no-multiple", "Disable multiple backtracking."));
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()                                        [Options]
// ----------------------------------------------------------------------------

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

    // Parse mapping block.
    getOptionValue(options.mappingBlock, parser, "mapping-block");

    // Parse mapping options.
    getOptionValue(options.errorsPerRead, parser, "errors");
    getOptionValue(options.seedLength, parser, "seed-length");
    options.mismatchesOnly = isSet(parser, "no-gaps");

    // Parse genome index prefix.
    getIndexPrefix(options, parser);

    // Parse genome index type.
    getIndexType(options, parser);

    // Parse output format.
    getOutputFormat(options, parser);

    // Parse output file.
    getOutputFile(options.mappedReadsFile, options, parser, options.readsFile, "");

    // Parse debug options.
    options.noVerify = isSet(parser, "no-verify");
    options.noDump = isSet(parser, "no-dump");
    options.noMultiple = isSet(parser, "no-multiple");

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TFormat, typename TMapperConfig, typename TReadsConfig>
int runMapper(Options & options, TMapperConfig const & /* config */, TReadsConfig const & /* config */)
{
    typedef Genome<void>                                                            TGenome;
    typedef GenomeIndex<TGenome, TIndex>                                            TGenomeIndex;
    typedef Reads<void, TReadsConfig>                                               TReads;
    typedef ReadsLoader<void, TReadsConfig>                                         TReadsLoader;
    typedef Writer<TGenome, TReads, TFormat, typename TMapperConfig::TDistance>     TWriter;
    typedef Mapper<TReads, TWriter, void, TMapperConfig>                            TMapper;

    TFragmentStore      store;
    TGenome             genome(store);
    TGenomeIndex        genomeIndex(genome);
    TReads              reads(store);
    TReadsLoader        readsLoader(reads);
    TWriter             writer(genome, options.noDump);
//    TMapper             mapper(writer, options.noVerify);

    double start, finish;

    // Load genome.
    std::cout << "Loading genome:\t\t\t" << std::flush;
    start = sysTime();
    if (!load(genome, options.genomeFile))
    {
        std::cerr << "Error while loading genome" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    // Load genome index.
    std::cout << "Loading genome index:\t\t" << std::flush;
    start = sysTime();
    if (!load(genomeIndex, options.genomeIndexFile))
    {
        std::cout << "Error while loading genome index" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    // Open reads file.
    start = sysTime();
    if (!open(readsLoader, options.readsFile))
    {
        std::cerr << "Error while opening reads file" << std::endl;
        return 1;
    }

    // Open output file.
    if (!open(writer, options.mappedReadsFile))
    {
        std::cerr << "Error while opening output file" << std::endl;
        return 1;
    }

    // Configure writer.
    writeAlignments(writer, options.outputCigar);

    // Reserve space for reads.
    if (options.mappingBlock < MaxValue<int>::VALUE)
        reserve(reads, options.mappingBlock);
    else
        reserve(reads);

    // Process reads in blocks.
    while (!atEnd(readsLoader))
    {
        // Load reads.
        std::cout << "Loading reads:\t\t\t" << std::flush;
        if (!load(readsLoader, options.mappingBlock))
        {
            std::cerr << "Error while loading reads" << std::endl;
            return 1;
        }
        finish = sysTime();
        std::cout << finish - start << " sec" << std::endl;
        std::cout << "Reads count:\t\t\t" << reads.readsCount << std::endl;

        // Pass reads to writer.
        setReads(writer, reads);

        // Configure mapper.
        TMapper mapper(reads, writer, options.noVerify);
        setSeedLength(mapper, options.seedLength);
        setReads(mapper, reads);

        // Map reads.
        start = sysTime();
        mapReads(mapper, genomeIndex, options.errorsPerRead);
        finish = sysTime();
        std::cout << "Mapping time:\t\t\t" << std::flush;
        std::cout << finish - start << " sec" << std::endl;

        // Clear mapped reads.
        clear(reads);
    }

    // Close output file.
    close(writer);

    // Close reads file.
    close(readsLoader);

    return 0;
}

template <typename TIndex, typename TFormat, typename TDistance, typename TStrategy, typename TBacktracking>
int runMapper(Options & options)
{
    typedef ReadMapperConfig<TDistance, TStrategy, MultipleBacktracking>    TMapperConfig;

    typedef typename IsSameType<TFormat, Sam>::Type                         TUseReadStore;
    typedef typename IsSameType<TFormat, Sam>::Type                         TUseReadNameStore;
    typedef ReadsConfig<TUseReadStore, TUseReadNameStore>                   TReadsConfig;

    return runMapper<TIndex, TFormat>(options, TMapperConfig(), TReadsConfig());
}

// ----------------------------------------------------------------------------
// Functions configure*()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TFormat, typename TDistance, typename TStrategy>
int configureMapperBacktracking(Options & options)
{
    if (options.noMultiple)
        return runMapper<TIndex, TFormat, TDistance, TStrategy, SingleBacktracking>(options);
    else
        return runMapper<TIndex, TFormat, TDistance, TStrategy, MultipleBacktracking>(options);
}

template <typename TIndex, typename TFormat, typename TDistance>
int configureMapperStrategy(Options & options)
{
    switch (options.mappingMode)
    {
    case Options::ANY_BEST:
        return configureMapperBacktracking<TIndex, TFormat, TDistance, AnyBest>(options);

    case Options::ALL_BEST:
        return configureMapperBacktracking<TIndex, TFormat, TDistance, AllBest>(options);

    case Options::ALL:
        return configureMapperBacktracking<TIndex, TFormat, TDistance, All>(options);

    default:
        return 1;
    }
}

template <typename TIndex, typename TFormat>
int configureMapperDistance(Options & options)
{
    if (options.mismatchesOnly)
        return configureMapperStrategy<TIndex, TFormat, HammingDistance>(options);
    else
        return configureMapperStrategy<TIndex, TFormat, EditDistance>(options);
}

template <typename TIndex>
int configureOutputFormat(Options & options)
{
    switch (options.outputFormat)
    {
    case Options::RAW:
        return configureMapperDistance<TIndex, Raw>(options);

    case Options::SAM:
        return configureMapperDistance<TIndex, Sam>(options);

    case Options::SAM_NO_CIGAR:
        options.outputCigar = false;
        return configureMapperDistance<TIndex, Sam>(options);

    default:
        return 1;
    }
}

int configureGenomeIndex(Options & options)
{
    switch (options.genomeIndexType)
    {
        case Options::INDEX_ESA:
            return configureOutputFormat<TGenomeEsa>(options);

        case Options::INDEX_SA:
            return configureOutputFormat<TGenomeSa>(options);

//    case Options::INDEX_QGRAM:
//        return configureOutputFormat<TGenomeQGram>(options);

        case Options::INDEX_FM:
            return configureOutputFormat<TGenomeFM>(options);

        default:
            return 1;
    }
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    return configureGenomeIndex(options);
}
