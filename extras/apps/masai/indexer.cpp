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
// This file contains the masai_indexer application.
// ==========================================================================

#define SEQAN_EXTRAS_MASAI_DISABLE_MMAP

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "options.h"
#include "index.h"

using namespace seqan;


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options : public MasaiOptions
{
    CharString genomeFile;
    CharString genomeIndexFile;
    IndexType  genomeIndexType;

    Options() :
        MasaiOptions(),
        genomeIndexType(INDEX_SA)
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
    setAppName(parser, "masai_indexer");
    setShortDescription(parser, "Masai Indexer");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta fa");

    addSection(parser, "Genome Index Options");

    setIndexType(parser, options);
    setIndexPrefix(parser);
    
    addSection(parser, "Output Options");

    setTmpFolder(parser);
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

    // Parse genome index prefix.
    getIndexPrefix(options, parser);

    // Parse genome index type.
    getIndexType(options, parser);

    // Parse tmp folder.
    getTmpFolder(options, parser);

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function runIndexer()
// ----------------------------------------------------------------------------

template <typename TIndex>
int runIndexer(Options & options)
{
    typedef Genome<>                        TGenome;
    typedef GenomeIndex<TGenome, TIndex>    TGenomeIndex;

    TFragmentStore      store;
    TGenome             genome(store);
    TGenomeIndex        genomeIndex(genome);

    double start, finish;

    std::cout << "Loading genome:\t\t\t" << std::flush;
    start = sysTime();
    if (!load(genome, options.genomeFile))
    {
        std::cerr << "Error while loading genome" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    std::cout << "Building genome index:\t\t" << std::flush;
    start = sysTime();
    build(genomeIndex);
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    std::cout << "Dumping genome index:\t\t" << std::flush;
    start = sysTime();
    if (!dump(genomeIndex, options.genomeIndexFile))
    {
        std::cerr << "Error while dumping genome index" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    return 0;
}

// ----------------------------------------------------------------------------
// Functions configure*()
// ----------------------------------------------------------------------------

int configureIndex(Options & options)
{
    switch (options.genomeIndexType)
    {
    case Options::INDEX_ESA:
        return runIndexer<TGenomeEsa>(options);

    case Options::INDEX_SA:
        return runIndexer<TGenomeSa>(options);

    case Options::INDEX_QGRAM:
        return runIndexer<TGenomeQGram>(options);

    case Options::INDEX_FM:
        return runIndexer<TGenomeFM>(options);

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

    return configureIndex(options);
}
