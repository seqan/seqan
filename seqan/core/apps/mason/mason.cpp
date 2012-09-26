// ==========================================================================
//                          Mason - A Read Simulator
// ==========================================================================
// Copyright (C) 2010 Manuel Holtgrewe, FU Berlin
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
// Mason - A read simulator.
//
// Usage: mason --help
//       mason illumina [options] source file
//           Simulation of Illumina reads.
//       mason 454 [options] source file
//           Simulation of 454 reads.
//       mason sanger [options] source file
//           Simulation of Sanger reads.
// ==========================================================================
// This file only contains code to parse the first command line
// argument and delegates the actual simulation to the simulate_*.h
// headers.
// ==========================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/arg_parse.h>
#include <seqan/misc/misc_random.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/store.h>

#include "mason.h"
#include "simulate_illumina.h"
#include "simulate_454.h"
#include "simulate_sanger.h"

using namespace seqan;

/* Print global help.  This function is called when the user specifies
 * a wrong simulation command.
 */
void printHelpGlobal() {
    std::cerr << "Mason - A Read Simulator" << std::endl
              << "(c) 2010 by Manuel Holtgrewe" << std::endl
              << std::endl
              << "Usage: mason illumina [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << "       mason 454 [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << "       mason sanger [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << std::endl
              << "Call with 'mason READS-TYPE --help' to get detailed help." << std::endl;
}

int parseOptions(Options<Global> & options, const int argc, const char * argv[]) {
    if (argc == 1 ||
        (argc >= 2 && CharString(argv[1]) == "--help") ||
        (argc >= 2 && CharString(argv[1]) == "-h") ||
        (CharString(argv[1]) != "illumina" && CharString(argv[1]) != "454" &&
         CharString(argv[1]) != "sanger")) {
        printHelpGlobal();
        return 1;
    }

    if (CharString(argv[1]) == "illumina") {
        options.readsType = READS_TYPE_ILLUMINA;
    } else if (CharString(argv[1]) == "454") {
        options.readsType = READS_TYPE_454;
    } else if (CharString(argv[1]) == "sanger") {
        options.readsType = READS_TYPE_SANGER;
    } else {
        printHelpGlobal();
        return 1;
    }

    return 0;
}

int main(const int argc, const char * argv[]) {
    // Switch command type (which reads to simulate) or show global help.
    Options<Global> globalOptions;
    int ret = parseOptions(globalOptions, argc, argv);
    if (ret != 0)
        return ret;

    // Kick off read simulation, depending on the chosen command.  We
    // use a type-based compile-time dispatch for selecting the
    // appropriate function for simulation.  This if/then/else switch
    // selects the right code path.
    if (globalOptions.readsType == READS_TYPE_ILLUMINA) {
        ArgumentParser parser("mason");
        setUpArgumentParser(parser, IlluminaReads());
        Options<IlluminaReads> options;
        CharString referenceFilename;
        ArgumentParser::ParseResult res = parseArgumentsAndCheck(options, referenceFilename, parser, argc, argv);
        // If parsing was not successful then exit with code 1 if there were errors.  Otherwise, exit with code 0
        // (e.g. help was printed).
        if (res != seqan::ArgumentParser::PARSE_OK)
            return res == seqan::ArgumentParser::PARSE_ERROR;
        return simulateReads(options, referenceFilename, IlluminaReads());
    } else if (globalOptions.readsType == READS_TYPE_454) {
        ArgumentParser parser("mason");
        setUpArgumentParser(parser, LS454Reads());
        Options<LS454Reads> options;
        CharString referenceFilename;
        ArgumentParser::ParseResult res = parseArgumentsAndCheck(options, referenceFilename, parser, argc, argv);
        // If parsing was not successful then exit with code 1 if there were errors.  Otherwise, exit with code 0
        // (e.g. help was printed).
        if (res != seqan::ArgumentParser::PARSE_OK)
            return res == seqan::ArgumentParser::PARSE_ERROR;
        return simulateReads(options, referenceFilename, LS454Reads());
    } else if (globalOptions.readsType == READS_TYPE_SANGER) {
        ArgumentParser parser("mason");
        setUpArgumentParser(parser, SangerReads());
        Options<SangerReads> options;
        CharString referenceFilename;
        ArgumentParser::ParseResult res = parseArgumentsAndCheck(options, referenceFilename, parser, argc, argv);
        // If parsing was not successful then exit with code 1 if there were errors.  Otherwise, exit with code 0
        // (e.g. help was printed).
        if (res != seqan::ArgumentParser::PARSE_OK)
            return res == seqan::ArgumentParser::PARSE_ERROR;
        return simulateReads(options, referenceFilename, SangerReads());
    } else {
        SEQAN_ASSERT_FAIL("Invalid reads type!");
    }

    return 0;
}
