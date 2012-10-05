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

#include <ctime>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "options.h"
#include "sorter.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Options : public MasaiOptions
{
    CharString  genomeFile;

    CharString  readsFile;
    CharString  mappedReadsFile;

    CharString  sortedReadsFile;
    OutputFormat outputFormat;
    bool        outputCigar;

//    unsigned    errorsPerRead;
    bool        mismatchesOnly;
    unsigned    matchesPerRead;

    bool        dumpResults;

    Options() :
        MasaiOptions(),
        outputFormat(SAM),
        outputCigar(true),
//        errorsPerRead(5),
        mismatchesOnly(false),
        matchesPerRead(MaxValue<unsigned>::VALUE),
        dumpResults(true)
    {}
};

// ============================================================================

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
//    CharString rev  = "$Revision$";
//    CharString date = "$Date$";

    setAppName(parser, "masai_output_se");
    setShortDescription(parser, "Masai Output - Single End Mode");
    setCategory(parser, "Read Mapping");

    setVersion(parser, "0.4");
    setDate(parser, "October 2012");

    addDescription(parser, "Masai is a fast and sensitive read mapper based on approximate seeds and multiple backtracking.");
    addDescription(parser, "See \\fIhttp://www.seqan.de/projects/masai\\fP for more information.");
    addDescription(parser, "(c) Copyright 2011-2012 by Enrico Siragusa.");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP> <\\fIRAW FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta");
    setValidValues(parser, 1, "fastq");
    setValidValues(parser, 2, "raw");


    addSection(parser, "Mapping Options");

//    addOption(parser, ArgParseOption("e",  "errors", "Maximum number of errors per read.", ArgParseOption::INTEGER));
//    setMinValue(parser, "errors", "0");
//    setMaxValue(parser, "errors", "32");
//    setDefaultValue(parser, "errors", options.errorsPerRead);

    addOption(parser, ArgParseOption("ng", "no-gaps", "Do not align reads with gaps."));

    addOption(parser, ArgParseOption("m", "matches", "Maximum number of matches per read.", ArgParseOption::INTEGER));
//    setDefaultValue(parser, "matches", options.matchesPerRead);


    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-file", "Specify output file.", ArgParseOption::OUTPUTFILE));
    setRequired(parser, "output-file");

    addOption(parser, ArgParseOption("of", "output-format", "Select output format.", ArgParseOption::STRING));
    setValidValues(parser, "output-format", options.outputFormatList);
    setDefaultValue(parser, "output-format", options.outputFormatList[options.outputFormat]);


    addSection(parser, "Debug Options");

    addOption(parser, ArgParseOption("nd", "no-dump", "Do not dump results."));
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

    // Parse raw input files.
    getArgumentValue(options.mappedReadsFile, parser, 2);

    // Parse mapping options.
//    getOptionValue(options.errorsPerRead, parser, "errors");
    options.mismatchesOnly = isSet(parser, "no-gaps");
    getOptionValue(options.matchesPerRead, parser, "matches");

    // Parse output file.
    getOptionValue(options.sortedReadsFile, parser, "output-file");
    if (!isSet(parser, "output-file"))
    {
        options.sortedReadsFile = options.readsFile;
        // TODO(esiragusa): Guess output file extension.
        append(options.sortedReadsFile, ".out");
    }

    // Parse output format.
    getOptionValue(options.outputFormat, parser, "output-format", options.outputFormatList);

    // Parse debug options.
    options.dumpResults = !isSet(parser, "no-dump");

    return seqan::ArgumentParser::PARSE_OK;
}

// ============================================================================

template <typename TDistance, typename TFormat>
int runSorter(Options & options)
{
    typedef Sorter<>    TSorter;

    // TODO(esiragusa): Remove outputCigar from Sorter members.
    TSorter sorter(options.matchesPerRead, options.outputCigar, options.dumpResults);

    clock_t start, end;

    // Loading genome.
    std::cout << "Loading genome:\t\t\t" << std::flush;
    start = clock();
    if (!loadGenome(sorter.indexer, options.genomeFile))
    {
        std::cerr << "Error while loading genome" << std::endl;
        return 1;
    }
    end = clock();
    std::cout << ((end - start) / (double)CLOCKS_PER_SEC) << " sec" << std::endl;

    // Loading reads.
    std::cout << "Loading reads:\t\t\t" << std::flush;
    start = clock();
    if (!loadReads(sorter, options.readsFile))
    {
        std::cerr << "Error while loading reads" << std::endl;
        return 1;
    }
    end = clock();
    std::cout << ((end - start) / (double)CLOCKS_PER_SEC) << " sec" << std::endl;
    std::cout << "Reads count:\t\t\t" << sorter.readsCount << std::endl;

    // Sorting mapped reads.
    std::cout << "Sorting matches:\t\t" << std::flush;
    start = clock();
    sortMappedReads(sorter, options.mappedReadsFile, options.sortedReadsFile, TDistance(), TFormat());
    end = clock();
    std::cout << ((end - start) / (double)CLOCKS_PER_SEC) << " sec" << std::endl;

    std::cout << "Matches count:\t\t\t" << sorter.matchesCount << std::endl;

    return 0;
}

// ============================================================================

template <typename TDistance>
int configureOutputFormat(Options & options)
{
    switch (options.outputFormat)
    {
    case Options::RAW:
        return runSorter<TDistance, Raw>(options);

    case Options::SAM:
        return runSorter<TDistance, Sam>(options);

    case Options::SAM_NO_CIGAR:
        options.outputCigar = false;
        return runSorter<TDistance, Sam>(options);

    default:
        return 1;
    }
}

int mainWithOptions(Options & options)
{
    if (options.mismatchesOnly)
        return configureOutputFormat<HammingDistance>(options);
    else
        return configureOutputFormat<EditDistance>(options);
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
