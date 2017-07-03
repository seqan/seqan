// ===========================================================================
//                 SGIP - Solution of Graph Isomorphism Problem
// ===========================================================================
// Copyright (C) 2012 by Jialu Hu
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ===========================================================================
// Author: Jialu Hu <Jialu.Hu@fu-berlin.de>
// ===========================================================================
// This application is used to determine whether two given graphs are
// isomorphic or not by using heuristic approach.
// ===========================================================================

#include <time.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>
#include <seqan/arg_parse.h>

#include "sgip.h"
#include "sgip_base.h"
#include "sgip_output.h"

using namespace seqan;

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Enum FileOption
// --------------------------------------------------------------------------

enum FileOption
{
    ZERO,                  // graph given through parameter
    FIRST,                 // graph created from  orignalFile
    SECOND                 // graph created from  comparisive file
};

// --------------------------------------------------------------------------
// Enum SearchingType
// --------------------------------------------------------------------------

enum SearchingType
{
    HEURISTIC,            // heuristic searching approach
    BRUTEFORTH
};

// --------------------------------------------------------------------------
// Enum SgipOption
// --------------------------------------------------------------------------

// Option for sgip.
struct SgipOption
{
    // I/O options.
    CharString originalFile;        // name of original file(first graph)
    CharString comparFile;         // name of comparisive file(second graph)
    CharString outputFile;         // name of result file
    CharString compareFolder;

    // More options.
    bool autoMetric;               // search Metric dimension if true
    FileOption activeFile;
    SearchingType searchingType;
    CharString algorithm;          // Search strategy for metric dimension,e.g. Greedy, genetic etc.
    unsigned odimension;           // metric dimension of original graph specified by user
    unsigned cdimension;           // metric dimension of comparative graph specified by user
    bool showHelp;
    bool showVersion;
    bool isoCheck;                //to check whether two input graphs are isomorphic
    bool isPrintFile;             //print result if outputFile is specified
    bool isAllatOnce;
    int verbose;

    SgipOption()
    {
        algorithm = "greedy";
        autoMetric = false;
        activeFile = FIRST;
        searchingType = HEURISTIC;
        showHelp = 0;
        showVersion = 0;
        isoCheck = 0;
        isPrintFile = false;
        odimension = 3;
        cdimension = 3;
        verbose = 0;
        isAllatOnce = false;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function _sgip()
// --------------------------------------------------------------------------

// Test sgip with Options.
template <typename TOption>
int _sgip(TOption & options)
{
    typedef Graph<Directed<> > TGraph;
    typedef String<bool>       TMat;
    TMat leastmat1, leastmat2;
    TGraph g1, g2;
    char const * file1, * file2;
    if (!options.isoCheck)
    {
        file1 = toCString(options.originalFile);
        if (!_createGraph(g1, SivaLab(), file1))
            return 1;

        getCanonicalLabel(leastmat1, g1);
        if (options.verbose > 1)
            outputLabel(leastmat1, FFFF());
    }
    else
    {
        file1 = toCString(options.originalFile);
        file2 = toCString(options.comparFile);
        if (!_createGraph(g1, SivaLab(), file1))
            return 1;

        if (!_createGraph(g2, SivaLab(), file2))
            return 1;

        if (checkIsomorphic(g1, g2))
        {
            if (options.verbose > 0)
                std::cout << "They are isomorphic!" << std::endl;
        }
        else
        {
            if (options.verbose > 0)
                std::cout << "They are not isomorphic!" << std::endl;
        }
    }
    return 0;
}

// --------------------------------------------------------------------------
// Function setupParser()
// --------------------------------------------------------------------------

// Set up parser.
template <typename TParser>
void _setupParser(TParser & parser)
{
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setShortDescription(parser, "Solution of Graph Isomorphism Problem");
    addUsageLine(parser, "-o <original graph> [Option]");    
    addSection(parser, "Mandatory Options");
    addOption(parser, ArgParseOption("o", "original", "File containing original graph", ArgParseArgument::INPUT_FILE,"IN"));
    setRequired(parser, "o");
    addSection(parser, "Main Options");
    addOption(parser, ArgParseOption("a", "algorithm", "Algorithm used for searching metric dimension", ArgParseArgument::STRING));
    setDefaultValue(parser, "algorithm", "greedy");
    addOption(parser, ArgParseOption("s", "searching", "Searching algorithm used for detecting resolving set,heuristic 0 bruteforce 1.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "searching", "0");
    addOption(parser, ArgParseOption("i", "isomorphism", "To check whether two given graphs are isomorphic"));
    addOption(parser, ArgParseOption("c", "comparative", "File containing comparative graph", ArgParseArgument::INPUT_FILE,"IN"));
    addOption(parser, ArgParseOption("od", "odimension", "Specified initial dimension of original graph by user", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "odimension", "3");
    addOption(parser, ArgParseOption("cd", "cdimension", "Specified initial dimension of comparative graph by user", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "cdimension", "3");
    addOption(parser, ArgParseOption("ad", "directory", "test for all graphs in a specified directory", ArgParseArgument::STRING));
    addOption(parser, ArgParseOption("v", "verbose", "control the level of output files ", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "verbose", "0");
    addSection(parser, "Output Option");
    addOption(parser, ArgParseOption("r", "result", "File containing canonical label of input graphs", ArgParseArgument::OUTPUT_FILE,"OUT"));
}

// --------------------------------------------------------------------------
// Function _parseOptions()
// --------------------------------------------------------------------------

template <typename TOption>
ArgumentParser::ParseResult _parseOptions(TOption & options,
                                          int argc,
                                          char const ** argv)
{
    // Setup ArgumentParsre.
    ArgumentParser parser("sgip");
    _setupParser(parser);

    // Parse commandline.
    ArgumentParser::ParseResult res=parse(parser,argc,argv);

    // Only extract options if the program will continue after parse().
    if(res != ArgumentParser::PARSE_OK)
        return res;
        
    // Extract option value.
    getOptionValue(options.originalFile, parser, "o");
    if (isSet(parser, "i"))
    {
        if (!isSet(parser, "c"))
        {
            if (!isSet(parser, "ad"))
            {
                std::cerr << "sgip" << ":comparative file has not been specified!" << std::endl;
                printShortHelp(parser, std::cerr);
                return ArgumentParser::PARSE_ERROR;
            }
        }
        options.isoCheck = true;
        getOptionValue(options.comparFile, parser, "c");
    }
    if (isSet(parser, "r"))
    {
        options.isPrintFile = true;
        getOptionValue(options.outputFile, parser, "r");
    }
    if (isSet(parser, "od"))
        getOptionValue(options.odimension, parser, "od");
    if (isSet(parser, "cd"))
        getOptionValue(options.cdimension, parser, "cd");
    if (isSet(parser, "a"))
        getOptionValue(options.algorithm, parser, "a");
    if (isSet(parser, "v"))
        getOptionValue(options.verbose, parser, "v");
    if (isSet(parser, "s"))
    {
        int s;
        getOptionValue(s, parser, "s");
        options.searchingType = static_cast<SearchingType>(s);
    }
    printVersion(parser, std::cerr);
    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.
int main(int argc, char const ** argv)
{
    // Initial options and time variables.
    SgipOption options;
    time_t start, end;
    double dif;
    int ret = 0;
    time(&start);

    // Setup parser and parse the commandline.
    ArgumentParser::ParseResult res = _parseOptions(options, argc, argv);
    if(res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Finally, launch the program.
    if (!options.isAllatOnce)
        ret = _sgip(options);
    time(&end);
    dif = difftime(end, start);
    if (options.verbose > 2)
        std::cout << std::setprecision(10) << "escaped time:" << dif << "seconds" << std::endl;
    return ret;
}
