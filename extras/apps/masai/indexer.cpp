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

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "options.h"
#include "indexer.h"

using namespace seqan;


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options()
// ----------------------------------------------------------------------------

struct Options : public MasaiOptions
{
    CharString genomeFile;
    CharString genomeIndexFile;
    IndexType  genomeIndexType;

    Options() :
        MasaiOptions(),
        genomeIndexType(INDEX_ESA)
    {}
};

// ============================================================================
// Functions
// ============================================================================

void setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setAppName(parser, "masai_indexer");
    setShortDescription(parser, "Masai Indexer");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta");

    addSection(parser, "Genome Index Options");

    addOption(parser, ArgParseOption("x", "index", "Select genome index type.", ArgParseOption::STRING));
    setValidValues(parser, "index", options.indexTypeList);
    setDefaultValue(parser, "index", options.indexTypeList[options.genomeIndexType]);

    addOption(parser, ArgParseOption("xp", "index-prefix", "Specify genome index prefix name.", ArgParseOption::STRING));
}

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Parse genome input file.
    getArgumentValue(options.genomeFile, parser, 0);

    // Parse genome index prefix.
    getOptionValue(options.genomeIndexFile, parser, "index-prefix");
    if (!isSet(parser, "index-prefix"))
    {
        // TODO(esiragusa): Trim extension .fasta from genomeIndexFile.
        options.genomeIndexFile = options.genomeFile;
    }

    // Parse genome index type.
    getOptionValue(options.genomeIndexType, parser, "index", options.indexTypeList);

    return seqan::ArgumentParser::PARSE_OK;
}

// ============================================================================

template <typename TIndex, typename TDepth>
void visit(TIndex & index, TDepth depth)
{
    typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type  TIndexIterator;
    TIndexIterator readsIt(index);

    do
    {
        std::cout << representative(readsIt) << std::endl;
        if (repLength(readsIt) >= depth || !goDown(readsIt))
            if (!goRight(readsIt))
                while (goUp(readsIt) && !goRight(readsIt))
                    ;
    }
    while (!isRoot(readsIt));
}

template <typename TGenomeIndex>
int executeIndexer(Options & options)
{
    TFragmentStore          store;
    Indexer<TGenomeIndex>   indexer(store);

    double start, finish;

    std::cout << "Loading genome:\t\t\t" << std::flush;
    start = sysTime();
    if (!loadGenome(indexer, options.genomeFile))
    {
        std::cerr << "Error while loading genome" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    std::cout << "Building genome index:\t\t" << std::flush;
    start = sysTime();
    indexGenome(indexer);
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    std::cout << "Dumping genome index:\t\t" << std::flush;
    start = sysTime();
    if (!dumpIndexedGenome(indexer, options.genomeIndexFile))
    {
        std::cerr << "Error while dumping genome index" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    return 0;
}

int mainWithOptions(Options & options)
{
    switch (options.genomeIndexType)
    {
    case Options::INDEX_ESA:
        return executeIndexer<TGenomeFM>(options);

    case Options::INDEX_SA:
        return executeIndexer<TGenomeSa>(options);

    case Options::INDEX_QGRAM:
        return executeIndexer<TGenomeQGram>(options);

//    case Options::INDEX_FM:
//        return executeIndexer<TGenomeFM>(options);

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
