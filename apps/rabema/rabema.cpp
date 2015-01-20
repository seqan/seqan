// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Read mapper benchmark tool.  This is the top level program that
// parses the command line parameters and then calls the main routine
// from the headers "build_gold_standard.h" and "evaluation.h".  These
// headers contain the code for building the gold standard and
// comparing the read mapper results against the gold standard.  They
// also contain the command line parsing code used here.
//
// Usage: rabema build_standard -o GOLD_STANDARD.gsi GENOME.fasta GOLD.sam
//        rabema compare GENOME.fasta GOLD_STANDARD.gsi READ_MAPPER_OUTPUT.sam
// ==========================================================================

#include "rabema.h"

#include <iostream>

#include "evaluation.h"
#include "build_gold_standard.h"

/* Print global help.  This function is called when the user specifies
   a wrong simulation command.
 */
void printHelpGlobal()
{
    std::cerr << "RABEMA - Read Alignment Benchmark" << std::endl
              << "(c) 2010 by Manuel Holtgrewe" << std::endl
              << std::endl
              << "Usage: rabema build_standard -o GOLD_STANDARD.gsi GENOME.fasta GOLD.sam" << std::endl
              << "       rabema compare GENOME.fasta GOLD_STANDARD.gsi READ_MAPPER_OUTPUT.sam" << std::endl
              << std::endl
              << "Call with 'rabema COMMAND --help' to get detailed help." << std::endl;
}

int parseOptions(Options<Global> & options, const int argc, const char * argv[])
{
    if (argc == 1 ||
        (argc >= 2 && CharString(argv[1]) == "--help") ||
        (argc >= 2 && CharString(argv[1]) == "-h"))
    {
        printHelpGlobal();
        return 0;
    }
    else if (argc >= 2 && CharString(argv[1]) != "build_standard" && CharString(argv[1]) != "compare")
    {
        printHelpGlobal();
        return 1;
    }

    if (CharString(argv[1]) == "build_standard")
    {
        options.selectedCommand = COMMAND_BUILD_STANDARD;
    }
    else if (CharString(argv[1]) == "compare")
    {
        options.selectedCommand = COMMAND_EVALUATE;
    }
    else
    {
        printHelpGlobal();
        return 1;
    }

    return 0;
}

int main(int argc, char const ** argv)
{
    Options<Global> globalOptions;
    int ret = parseOptions(globalOptions, argc, argv);
    if (ret)
        return ret;

    if (globalOptions.selectedCommand == COMMAND_BUILD_STANDARD)
    {
        CommandLineParser parser;
        setUpCommandLineParser(parser, BuildGoldStandard());
        Options<BuildGoldStandard> options;
        int ret = parseCommandLineAndCheck(options, parser, argc, argv);
        if (options.showHelp)
            return 0;

        if (ret != 0)
            return ret;

        return buildGoldStandard(options);
    }
    else if (globalOptions.selectedCommand == COMMAND_EVALUATE)
    {
        CommandLineParser parser;
        setUpCommandLineParser(parser, EvaluateResults());
        Options<EvaluateResults> options;
        int ret = parseCommandLineAndCheck(options, parser, argc, argv);
        if (options.showHelp)
            return 0;

        if (ret != 0)
            return ret;

        return evaluateReadMapperResult(options);
    }

    return 0;
}
