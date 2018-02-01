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
// This appication performs String Similarity Join over Edit Distance
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
        #pragma message("OpenMP not found! Shared-memory parallelization will be disabled in join tool.")
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
    CharString      inputType;
    CharString      resultsFile;

    unsigned        threadsCount;
    unsigned        maxErrors;
    unsigned        seedLength;

    bool            online;
    bool            hugeDb;

    Options() :
        threadsCount(8),
        maxErrors(0),
        seedLength(0),
        online(false),
        hugeDb(false)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()                              [ArgumentParser]
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser)
{
    setAppName(parser, "join");
    setShortDescription(parser, "EDBT/ICDT 2013 Join");
    setCategory(parser, "Databases");

    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIDATABASE FILE\\fP> <\\fIERRORS\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INTEGER));

    // Set errors range to [0,32].
    setMinValue(parser, 1, "0");
    setMaxValue(parser, 1, "32");

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
#endif

    // Add output file option.
    addOption(parser, ArgParseOption("o", "output-file", "Specify an output file.", ArgParseOption::STRING));
    setDefaultValue(parser, "output-file", "result_track2.out");

    // Add online join option.
    addOption(parser, ArgParseOption("l", "online", "Perform online join."));

    // Add seed length option.
    addOption(parser, ArgParseOption("sl", "seed-length", "Minimum seed length.", ArgParseOption::INTEGER));
    setMinValue(parser, "seed-length", "0");
    setMaxValue(parser, "seed-length", "100");
    setDefaultValue(parser, "seed-length", 0);
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

    // Parse errors.
    getArgumentValue(options.maxErrors, parser, 1);

    // Parse the number of threads.
    getOptionValue(options.threadsCount, parser, "threads");

    // Parse input type.
    getOptionValue(options.inputType, parser, "input-type");

    // Parse huge db option.
#ifdef SEARCHJOIN_HUGEDB
    options.hugeDb = isSet(parser, "huge");
#endif

    // Parse output file.
    getOptionValue(options.resultsFile, parser, "output-file");

    // Parse online join option.
    options.online = isSet(parser, "online");

    // Parse seed length option.
    getOptionValue(options.seedLength, parser, "seed-length");

    return ArgumentParser::PARSE_OK;
}

template <typename TJoiner>
inline void printStats(TJoiner & joiner)
{
    std::cout << "Length filter:\t\t\t\t" << joiner.verifier.lengthFilter << std::endl;
    std::cout << "Verifications:\t\t\t\t" << joiner.verifier.verifications << std::endl;
}

template <typename TText, typename TIndex, typename TDbQuerySpec, typename TDelegate>
inline void printStats(DbFinder<TText, TIndex, TDbQuerySpec, TDelegate, Exact> &)
{}

// ----------------------------------------------------------------------------
// Function runJoiner()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndex, typename TSpec>
int runJoiner(Options & options, TText const & /* tag */, TIndex const & /* tag */, TSpec /* tag */)
{
    typedef Db<TText>                                       TDb;
    typedef Writer<TDb, TDb, Join>                          TWriter;
    typedef DbFinder<TText, TIndex, void, TWriter, TSpec>   TJoiner;

    double start, finish;

#ifdef _OPENMP
    // Set the number of threads that OpenMP can spawn.
    omp_set_num_threads(options.threadsCount);
#endif

    // Instantiate objects.
    TDb db;
    TWriter writer(db);
    TJoiner joiner(db, writer);

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
    setMinSeedLength(joiner, options.seedLength);
    index(joiner);
    finish = sysTime();
    std::cout << "Indexing time:\t\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Open results file.
    if (!open(writer, options.resultsFile))
    {
        std::cerr << "Error while opening results file" << std::endl;
        return 1;
    }

    // Prepare database.
    start = sysTime();
    db.errors = options.maxErrors;
    prepare(joiner, db);
    finish = sysTime();
    std::cout << "Preparation time:\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;

    // Join database.
    start = sysTime();
    execute(joiner);
    finish = sysTime();
    std::cout << "Join time:\t\t\t\t" << std::flush;
    std::cout << finish - start << " sec" << std::endl;
    printStats(joiner);

    // Close results file.
    close(writer);

    // Print some additional stats.
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
            return runJoiner(options, TDbDna(), Nothing(), Online());
        }
        else
        {
#ifdef SEARCHJOIN_HUGEDB
            if (options.hugeDb)
            {
                if (options.maxErrors == 0)
                    return runJoiner(options, TDbDna(), TDbDnaSaHuge(), Exact());
                else
                {
                    if (options.threadsCount > 1)
                        return runJoiner(options, TDbDna(), TDbDnaSaHuge(), Parallel());
                    else
                        return runJoiner(options, TDbDna(), TDbDnaSaHuge(), Nothing());
                }
            }
            else
            {
#endif // SEARCHJOIN_HUGEDB
                if (options.maxErrors == 0)
                    return runJoiner(options, TDbDna(), TDbDnaSaSmall(), Exact());
                else
                {
                    if (options.threadsCount > 1)
                        return runJoiner(options, TDbDna(), TDbDnaSaSmall(), Parallel());
                    else
                        return runJoiner(options, TDbDna(), TDbDnaSaSmall(), Nothing());
                }
#ifdef SEARCHJOIN_HUGEDB
            }
#endif // SEARCHJOIN_HUGEDB
        }
    }
    else if (isEqual(options.inputType, "geo"))
    {
        if (options.online)
        {
            return runJoiner(options, TDbGeo(), Nothing(), Online());
        }
        else
        {
#ifdef SEARCHJOIN_HUGEDB
            if (options.hugeDb)
            {
                if (options.maxErrors == 0)
                    return runJoiner(options, TDbGeo(), TDbGeoSaHuge(), Exact());
                else
                {
                    if (options.threadsCount > 1)
                        return runJoiner(options, TDbGeo(), TDbGeoSaHuge(), Parallel());
                    else
                        return runJoiner(options, TDbGeo(), TDbGeoSaHuge(), Nothing());
                }
            }
            else
            {
#endif // SEARCHJOIN_HUGEDB
                if (options.maxErrors == 0)
                    return runJoiner(options, TDbGeo(), TDbGeoSaSmall(), Exact());
                else
                {
                    if (options.threadsCount > 1)
                        return runJoiner(options, TDbGeo(), TDbGeoSaSmall(), Parallel());
                    else
                        return runJoiner(options, TDbGeo(), TDbGeoSaSmall(), Nothing());
                }
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

