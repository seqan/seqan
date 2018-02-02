// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This appication performs String Similarity Search over Edit Distance
// on a database of DNA reads or Geographical names.
//
// Siragusa, E., Weese D., & Reinert, K. (2013).
// Scalable String Similarity Search/Join with Approximate Seeds and Multiple Backtracking.
// EDBT/ICDT ’13, March 18 – 22 2013, Genoa, Italy
// ==========================================================================

#include <seqan/platform.h>

#if defined(_OPENMP)
    #include <omp.h>

    #if defined(STDLIB_GNU)
        #include <parallel/algorithm>
    #endif
#else
    #if !defined(SEQAN_IGNORE_MISSING_OPENMP) || (SEQAN_IGNORE_MISSING_OPENMP == 0)
        #pragma message("OpenMP not found! Shared-memory parallelization will be disabled in search tool.")
    #endif  // #if !defined(SEQAN_IGNORE_MISSING_OPENMP) || (SEQAN_IGNORE_MISSING_OPENMP == 0)
#endif  // #ifdef _OPENMP

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/parallel.h>

#include "db.h"
#include "finder.h"
#include "writer.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString      databaseFile;
    CharString      queryFile;
    CharString      inputType;
    CharString      resultsFile;

    unsigned        threadsCount;
    unsigned        seedLength;

    bool            online;
    bool            noWait;
    bool            hugeDb;

    Options() :
        threadsCount(8),
        seedLength(0),
        online(false),
        noWait(false),
        hugeDb(false)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _waitForQueryFile()
// ----------------------------------------------------------------------------

void _waitForQueryFile()
{
    typedef std::fstream        TStream;

    static const char * INDEX_GEN_FINISHED  = "/competition/track1/index_generation_finished.txt";
    static const char * QUERIES_AVAILABLE   = "/competition/track1/queries_available.txt";

    // Signal that index was built.
    TStream indexGenFinished(INDEX_GEN_FINISHED, std::ios::out);
    indexGenFinished.close();

    // Wait for queries.
    TStream queries_available;
    do
    {
        queries_available.open(QUERIES_AVAILABLE, std::ios::in);
    }
    while (!queries_available.is_open());
}

// ----------------------------------------------------------------------------
// Function setupArgumentParser()                              [ArgumentParser]
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser)
{
    setAppName(parser, "search");
    setShortDescription(parser, "EDBT/ICDT 2013 Search");
    setCategory(parser, "Databases");

    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIDATABASE FILE\\fP> <\\fIQUERY FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE));

    addSection(parser, "Options");

    // Add threads option.
    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
    setMaxValue(parser, "threads", "64");
    setDefaultValue(parser, "threads", "8");

    // Add input type option.
    addOption(parser, ArgParseOption("i", "input-type", "Specify the type of input.", ArgParseOption::STRING));
    setValidValues(parser, "input-type", "dna geo");
    setRequired(parser, "input-type", true);

    // Add huge db option.
#ifdef SEARCHJOIN_HUGEDB
    addOption(parser, ArgParseOption("g", "huge", "Required if the db contains more than 16M entries."));
#endif // SEARCHJOIN_HUGEDB

    // Add output file option.
    addOption(parser, ArgParseOption("o", "output-file", "Specify an output file.", ArgParseOption::STRING));
    setDefaultValue(parser, "output-file", "result_track1.out");

    // Add online search option.
    addOption(parser, ArgParseOption("l", "online", "Perform online search."));

    // Add seed length option.
    addOption(parser, ArgParseOption("sl", "seed-length", "Minimum seed length.", ArgParseOption::INTEGER));
    setMinValue(parser, "seed-length", "0");
    setMaxValue(parser, "seed-length", "100");
    setDefaultValue(parser, "seed-length", 0);

    // Add no wait option.
    addOption(parser, ArgParseOption("n", "no-wait", "Do not wait for query file."));
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()                                        [Options]
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse database input file.
    getArgumentValue(options.databaseFile, parser, 0);

    // Parse query input file.
    getArgumentValue(options.queryFile, parser, 1);

    // Parse the number of threads.
    getOptionValue(options.threadsCount, parser, "threads");

    // Parse input type.
    getOptionValue(options.inputType, parser, "input-type");

    // Parse huge db option.
#ifdef SEARCHJOIN_HUGEDB
    options.hugeDb = isSet(parser, "huge");
#endif // SEARCHJOIN_HUGEDB

    // Parse output file.
    getOptionValue(options.resultsFile, parser, "output-file");

    // Parse online search option.
    options.online = isSet(parser, "online");

    // Parse seed length option.
    getOptionValue(options.seedLength, parser, "seed-length");

    // Parse no wait option.
    options.noWait = isSet(parser, "no-wait");

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function runSearcher()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TSpec>
int runSearcher(Options & options, TText const & /* tag */, TIndex const & /* tag */, TSpec /* tag */, Nothing const & /* tag */)
{
    typedef Db<TText>                                       TDb;
    typedef Db<TText, Query>                                TDbQuery;
    typedef Writer<TDb, TDbQuery, Search>                   TWriter;
    typedef DbFinder<TText, TIndex, Query, TWriter, TSpec>  TSearcher;

    double start, finish;

#ifdef _OPENMP
    // Set the number of threads that OpenMP can spawn.
    omp_set_num_threads(options.threadsCount);
#endif

    // Instantiate objects.
    TDb db;
    TDbQuery query;
    TWriter writer(db, query);
    TSearcher searcher(db, writer);

    // Load database.
    std::cout << "Loading database:\t\t\t" << std::flush;
    start = sysTime();
    if (!load(db, options.databaseFile))
    {
        std::cerr << "Error while loading database" << std::endl;
        return 1;
    }
    finish = sysTime();

    // Check database.
    if (IsSameType<TDbDnaSaSmall, TIndex>::VALUE || IsSameType<TDbGeoSaSmall, TIndex>::VALUE)
    {
        if (length(db.text) >= Power<2, 24>::VALUE)
        {
            std::cerr << "Please specify the option '--huge'" << std::endl;
            return 1;
        }
        if (db.maxLength >= Power<2, 8>::VALUE)
        {
            std::cerr << "Database strings are too long" << std::endl;
            return 1;
        }
    }

    std::cout << finish - start << " sec" << std::endl;
    std::cout << "Database entries:\t\t\t" << length(db.text) << std::endl;
    std::cout << "Min length:\t\t\t\t" << db.minLength << std::endl;
    std::cout << "Avg length:\t\t\t\t" << db.avgLength << std::endl;
    std::cout << "Max length:\t\t\t\t" << db.maxLength << std::endl;

    // Index database.
    start = sysTime();
    setMinSeedLength(searcher, options.seedLength);
    index(searcher);
    finish = sysTime();
    std::cout << "Indexing time:\t\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Wait for query file.
    if (!options.noWait)
    {
        std::cout << "Waiting queries..." << std::endl;
        _waitForQueryFile();
    }

    // Load queries.
    std::cout << "Loading queries:\t\t\t" << std::flush;
    start = sysTime();
    if (!load(query, options.queryFile))
    {
        std::cerr << "Error while loading queries" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    std::cout << "Queries count:\t\t\t\t" << length(query.text) << std::endl;
    std::cout << "Min length:\t\t\t\t" << query.minLength << std::endl;
    std::cout << "Avg length:\t\t\t\t" << query.avgLength << std::endl;
    std::cout << "Max length:\t\t\t\t" << query.maxLength << std::endl;

    std::cout << "Min errors:\t\t\t\t" << static_cast<unsigned>(getMinErrors(query)) << std::endl;
    std::cout << "Avg errors:\t\t\t\t" << static_cast<unsigned>(getAvgErrors(query)) << std::endl;
    std::cout << "Max errors:\t\t\t\t" << static_cast<unsigned>(getMaxErrors(query)) << std::endl;

    // Open results file.
    if (!open(writer, options.resultsFile))
    {
        std::cerr << "Error while opening results file" << std::endl;
        return 1;
    }

    // Prepare queries.
    start = sysTime();
    prepare(searcher, query);
    finish = sysTime();
    std::cout << "Preparation time:\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Search database.
    start = sysTime();
    execute(searcher);
    finish = sysTime();
    std::cout << "Search time:\t\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Close results file.
    close(writer);

    // Print some additional stats.
    std::cout << "Length filter:\t\t\t\t" << searcher.verifier.lengthFilter << std::endl;
    std::cout << "Verifications:\t\t\t\t" << searcher.verifier.verifications << std::endl;
    std::cout << "Results:\t\t\t\t" << writer.recordsCount << std::endl;

    return 0;
}

template <typename TText, typename TSpec, typename TIndex>
int runSearcher(Options & options, TText const & /* tag */, TIndex const & /* tag */, TSpec /* tag */, Online /* tag */)
{
    typedef Db<TText>                                       TDb;
    typedef Db<TText, Query>                                TDbQuery;
    typedef Writer<TDb, TDbQuery, Search>                   TWriter;
    typedef DbFinder<TText, TIndex, Query, TWriter, TSpec>  TSearcher;
    typedef DbFinder<TText, void, Query, TWriter, Online>   TSearcherOnline;

    double start, finish;

#ifdef _OPENMP
    // Set the number of threads that OpenMP can spawn.
    omp_set_num_threads(options.threadsCount);
#endif

    // Instantiate objects.
    TDb db;
    TDbQuery query;
    TDbQuery queryShort;
    TDbQuery queryLong;
    TWriter writer(db, queryLong);
    TSearcher searcher(db, writer);
    TSearcherOnline searcherOnline(db, writer);

    // Load database.
    std::cout << "Loading database:\t\t\t" << std::flush;
    start = sysTime();
    if (!load(db, options.databaseFile))
    {
        std::cerr << "Error while loading database" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    // Check database.
    if (IsSameType<TDbDnaSaSmall, TIndex>::VALUE || IsSameType<TDbGeoSaSmall, TIndex>::VALUE)
    {
        if (length(db.text) >= Power<2, 24>::VALUE)
        {
            std::cerr << "Please specify the option '--huge'" << std::endl;
            return 1;
        }
        if (db.maxLength >= Power<2, 8>::VALUE)
        {
            std::cerr << "Database strings are too long." << std::endl;
            return 1;
        }
    }

    std::cout << "Database entries:\t\t\t" << length(db.text) << std::endl;
    std::cout << "Min length:\t\t\t\t" << db.minLength << std::endl;
    std::cout << "Avg length:\t\t\t\t" << db.avgLength << std::endl;
    std::cout << "Max length:\t\t\t\t" << db.maxLength << std::endl;

    // Index database.
    start = sysTime();
    setMinSeedLength(searcher, options.seedLength);
    index(searcher);
    finish = sysTime();
    std::cout << "Indexing time:\t\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Wait for query file.
    if (!options.noWait)
    {
        std::cout << "Waiting queries..." << std::endl;
        _waitForQueryFile();
    }

    // Load queries.
    std::cout << "Loading queries:\t\t\t" << std::flush;
    start = sysTime();
    if (!load(query, options.queryFile))
    {
        std::cerr << "Error while loading queries" << std::endl;
        return 1;
    }
    finish = sysTime();
    std::cout << finish - start << " sec" << std::endl;

    // Split queries by length.
    split(queryShort, queryLong, query, options.seedLength);

    // Open results file.
    if (!open(writer, options.resultsFile))
    {
        std::cerr << "Error while opening results file" << std::endl;
        return 1;
    }

    std::cout << "Queries count:\t\t\t\t" << length(queryLong.text) << std::endl;
    std::cout << "Min length:\t\t\t\t" << queryLong.minLength << std::endl;
    std::cout << "Avg length:\t\t\t\t" << queryLong.avgLength << std::endl;
    std::cout << "Max length:\t\t\t\t" << queryLong.maxLength << std::endl;

    std::cout << "Min errors:\t\t\t\t" << static_cast<unsigned>(getMinErrors(queryLong)) << std::endl;
    std::cout << "Avg errors:\t\t\t\t" << static_cast<unsigned>(getAvgErrors(queryLong)) << std::endl;
    std::cout << "Max errors:\t\t\t\t" << static_cast<unsigned>(getMaxErrors(queryLong)) << std::endl;

    // Prepare queries.
    start = sysTime();
    prepare(searcher, queryLong);
    finish = sysTime();
    std::cout << "Preparation time:\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Search database.
    start = sysTime();
    execute(searcher);
    finish = sysTime();
    std::cout << "Search time:\t\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Switch query in writer.
    writer.query = queryShort;

    std::cout << "Queries count:\t\t\t\t" << length(queryShort.text) << std::endl;
    std::cout << "Min length:\t\t\t\t" << queryShort.minLength << std::endl;
    std::cout << "Avg length:\t\t\t\t" << queryShort.avgLength << std::endl;
    std::cout << "Max length:\t\t\t\t" << queryShort.maxLength << std::endl;

    std::cout << "Min errors:\t\t\t\t" << static_cast<unsigned>(getMinErrors(queryShort)) << std::endl;
    std::cout << "Avg errors:\t\t\t\t" << static_cast<unsigned>(getAvgErrors(queryShort)) << std::endl;
    std::cout << "Max errors:\t\t\t\t" << static_cast<unsigned>(getMaxErrors(queryShort)) << std::endl;

    // Prepare queries.
    start = sysTime();
    prepare(searcherOnline, queryShort);
    finish = sysTime();
    std::cout << "Preparation time:\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Search database.
    start = sysTime();
    execute(searcherOnline);
    finish = sysTime();
    std::cout << "Search time:\t\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Close results file.
    close(writer);

    // Print some additional stats.
    std::cout << "Length filter:\t\t\t\t" << searcher.verifier.lengthFilter + searcherOnline.verifier.lengthFilter << std::endl;
    std::cout << "Verifications:\t\t\t\t" << searcher.verifier.verifications + searcherOnline.verifier.verifications << std::endl;
    std::cout << "Results:\t\t\t\t" << writer.recordsCount << std::endl;

    return 0;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int mainWithOptions(Options & options)
{
    if (isEqual(options.inputType, "dna"))
    {
        if (options.online)
        {
            return runSearcher(options, TDbDna(), Nothing(), Online(), Nothing());
        }
        else
        {
#ifdef SEARCHJOIN_HUGEDB
            if (options.hugeDb)
            {
                if (options.threadsCount > 1)
                    return runSearcher(options, TDbDna(), TDbDnaSaHuge(), Parallel(), Nothing());
                else
                    return runSearcher(options, TDbDna(), TDbDnaSaHuge(), Nothing(), Nothing());
            }
            else
            {
#endif // SEARCHJOIN_HUGEDB
                if (options.threadsCount > 1)
                    return runSearcher(options, TDbDna(), TDbDnaSaSmall(), Parallel(), Nothing());
                else
                    return runSearcher(options, TDbDna(), TDbDnaSaSmall(), Nothing(), Nothing());
            }
#ifdef SEARCHJOIN_HUGEDB
        }
#endif // SEARCHJOIN_HUGEDB
    }
    else if (isEqual(options.inputType, "geo"))
    {
        if (options.online)
        {
            return runSearcher(options, TDbGeo(), Nothing(), Online(), Nothing());
        }
        else
        {
#ifdef SEARCHJOIN_HUGEDB
            if (options.hugeDb)
            {
                if (options.threadsCount > 1)
                    return runSearcher(options, TDbGeo(), TDbGeoSaHuge(), Parallel(), Online());
                else
                    return runSearcher(options, TDbGeo(), TDbGeoSaHuge(), Nothing(), Online());
            }
            else
            {
#endif // SEARCHJOIN_HUGEDB
                if (options.threadsCount > 1)
                    return runSearcher(options, TDbGeo(), TDbGeoSaSmall(), Parallel(), Online());
                else
                    return runSearcher(options, TDbGeo(), TDbGeoSaSmall(), Nothing(), Online());
#ifdef SEARCHJOIN_HUGEDB
            }
#endif // SEARCHJOIN_HUGEDB
        }
    }
    else
        return 1;
}

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res == ArgumentParser::PARSE_OK)
        return mainWithOptions(options);
    else
        return res;
}

