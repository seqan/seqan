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
// This file contains the masai_output_se application.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "options.h"
#include "sorter.h"
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

    CharString  readsFile;
    CharString  mappedReadsFile;

    CharString  sortedReadsFile;
    OutputFormat outputFormat;
    bool        outputCigar;

//    unsigned    errorsPerRead;
    bool        mismatchesOnly;
    unsigned    matchesPerRead;

    bool        noDump;

    Options() :
        MasaiOptions(),
        outputFormat(SAM),
        outputCigar(true),
//        errorsPerRead(5),
        mismatchesOnly(false),
        matchesPerRead(MaxValue<unsigned>::VALUE),
        noDump(false)
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
    setAppName(parser, "masai_output_se");
    setShortDescription(parser, "Masai Output - Single End Mode");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP> <\\fIREADS FILE\\fP> <\\fIRAW FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta fa");
    setValidValues(parser, 1, "fastq fasta fa");
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

    setTmpFolder(parser);
    setOutputFile(parser, options);
    addOption(parser, ArgParseOption("nc", "no-cigar", "Do not output CIGAR string. This only affects SAM output."));


    addSection(parser, "Debug Options");

    addOption(parser, ArgParseOption("nd", "no-dump", "Do not dump results."));
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

    // Parse raw input files.
    getArgumentValue(options.mappedReadsFile, parser, 2);

    // Parse mapping options.
//    getOptionValue(options.errorsPerRead, parser, "errors");
    options.mismatchesOnly = isSet(parser, "no-gaps");
    getOptionValue(options.matchesPerRead, parser, "matches");

    // Parse tmp folder.
    getTmpFolder(options, parser);

    // Parse output file.
    getOutputFile(options.sortedReadsFile, options, parser, options.readsFile, "_se");

    // Parse output format.
    getOutputFormat(options, options.sortedReadsFile);
    options.outputCigar = !isSet(parser, "no-cigar");

    // Parse debug options.
    options.noDump = isSet(parser, "no-dump");

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function runSorter()
// ----------------------------------------------------------------------------

template <typename TDistance, typename TFormat, typename TReadsConfig>
int runSorter(Options & options, TReadsConfig const & /* config */)
{
    typedef Genome<void>                                                        TGenome;
    typedef Reads<void, TReadsConfig>                                           TReads;
    typedef ReadsLoader<void, TReadsConfig>                                     TReadsLoader;
    typedef Writer<TGenome, TReads, TFormat, TDistance>                         TWriter;
    typedef Sorter<TWriter>                                                     TSorter;

    TFragmentStore      store;
    TGenome             genome(store);
    TReads              reads(store);
    TReadsLoader        readsLoader(reads);
    TWriter             writer(genome, options.noDump);
    TSorter             sorter(writer);

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

    // Open reads file.
    start = sysTime();
    if (!open(readsLoader, options.readsFile))
    {
        std::cerr << "Error while opening reads file" << std::endl;
        return 1;
    }

    // Open matches file.
    if (!open(sorter, options.mappedReadsFile))
    {
        std::cerr << "Error while opening matches file" << std::endl;
        return 1;
    }

    // Open output file.
    if (!open(writer, options.sortedReadsFile))
    {
        std::cerr << "Error while opening output file" << std::endl;
        return 1;
    }

    // Configure writer.
    writeAlignments(writer, options.outputCigar);

    // Reserve space for reads.
    reserve(reads);

    // Load all reads.
    std::cout << "Loading reads:\t\t\t" << std::flush;
    start = sysTime();
    load(readsLoader);
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;
    std::cout << "Reads count:\t\t\t" << reads.readsCount << std::endl;

    // Pass reads to writer.
    setReads(writer, reads);

    // Sort matches.
    std::cout << "Sorting matches:\t\t" << std::flush;
    start = sysTime();
    sort(sorter, options.matchesPerRead);
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;
    std::cout << "Matches count:\t\t\t" << sorter.matchesCount << std::endl;

    // Close output file.
    close(writer);

    // Close matches file.
    close(sorter);

    // Close reads file.
    close(readsLoader);

    return 0;
}

template <typename TDistance, typename TFormat>
int runSorter(Options & options)
{
    typedef typename IsSameType<TFormat, Sam>::Type                         TUseReadStore;
    typedef typename IsSameType<TFormat, Sam>::Type                         TUseReadNameStore;
    typedef ReadsConfig<TUseReadStore, TUseReadNameStore>                   TReadsConfig;

    return runSorter<TDistance, TFormat>(options, TReadsConfig());
}

// ----------------------------------------------------------------------------
// Functions configure*()
// ----------------------------------------------------------------------------

template <typename TDistance>
int configureOutputFormat(Options & options)
{
    switch (options.outputFormat)
    {
    case Options::RAW:
        return runSorter<TDistance, Raw>(options);

    case Options::SAM:
        return runSorter<TDistance, Sam>(options);

    default:
        return 1;
    }
}

int configureDistance(Options & options)
{
    if (options.mismatchesOnly)
        return configureOutputFormat<HammingDistance>(options);
    else
        return configureOutputFormat<EditDistance>(options);
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

    return configureDistance(options);
}
