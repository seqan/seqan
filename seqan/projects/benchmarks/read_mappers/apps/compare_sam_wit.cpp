/*==========================================================================
   SeqAn - The Library for Sequence Analysis
   http://www.seqan.de 
  ===========================================================================
   Copyright (C) 2010
  
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.
  
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.
  
  ===========================================================================
   Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
   Usage: compare_sam_wit [options] <contig.fasta> <result.sam> <golden.wit>

   Call as "compare_sam_wit --help" for a full list of options.
  
   This program is used to compare the read mapper results in a Sam file
   with the golden standard in a WIT file.  Log messages are printed to
   stderr, the result is written to the first line of stdout.

   The first line of stdout consists of a JSON encoded dictionary that has
   the following entries:

   * total_intervals       Total number of intervals with the given error
                           rate in the WIT file.
   * found_intervals       Number of intervals in the WIT file that were
                           found by the read mapper.
   * superflous_intervals  Number of alignments in the Sam file that do not
                           have their end position in an interval from the
                           WIT file and the alignment at this position has
                           a higher error rate than the specified one.
   * additional_intervals  Same as superflous_intervals, but these alignments
                           have an error rate that is as low as the specified
                           one or lower.  If this happens in the non-weighted
                           case then this is a bug in the benchmark tools.
                           For the weighted case, this happens if RazerS does
                           not find the alignment with low weighted but too
                           high unweighted error rate and the read mapper
                           generating the Sam file finds it.  Read the paper
                           and/or manual for more details.
                           In the non-weighte case, the program will exit
                           with a non-zero exit code and print a PANIC message.

   This .cpp file only contains the code to parse the command line arguments
   and kicks off the flesh of the program which is defined in
   compare_sam_wit.h.
  ===========================================================================*/

#include "compare_sam_wit.h"

#include <cassert>
#include <cstdlib>    // exit()

#include <algorithm>  // sort()
#include <map>
#include <set>

#include <seqan/basic.h>
#include <seqan/store.h>

#include "wit_store.h"
#include "return_codes.h"
#include "find_myers_ukkonen_reads.h"
#include "wit_store.h"

using namespace seqan;  // Remove some syntatic noise.


// The revision string of the program.
const char * kRevision = "0.0alpha";


// Set up the CommandLineParser options with options.
void setUpCommandLineParser(CommandLineParser &parser) {
    // Add a version string.
    std::string versionString = "compare_sam_wit ";
    versionString += kRevision;
    addVersionLine(parser, versionString);

    // Add usage lines.
    addTitleLine(parser, "Compare Sam hits against WIT file.");
    addUsageLine(parser, "[OPTIONS] <SEQUENCE FILE> <Sam FILE> <WIT FILE>");

    // Set options.
    addOption(parser, CommandLineOption("e", "max-error-rate", "the maximal error rate in percent, default: 0", OptionType::Int));
    addOption(parser, CommandLineOption("d", "distance-function", "the distance function to use, default: hamming", OptionType::String | OptionType::Label));
    addHelpLine(parser, "hamming          = Hamming distance");
    addHelpLine(parser, "edit             = Edit distance");
    addOption(parser, CommandLineOption("o", "out-file", "Path to the output file.  Use \"-\" for stdout, default is \"-\"", OptionType::String));
    addOption(parser, CommandLineOption("sm", "show-missed-intervals", "the missed intervals are printed to stderr for debugging if this option is set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("sh", "show-hit-intervals", "the hit intervals are printed to stderr for debugging if this option is set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("st", "show-try-hit-intervals", "The last positions tried to hit against an intervals are printd to stderr if set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("ss", "show-superflous-intervals", "The intervals that are in the Sam file but have a too bad score are printed to stderr if set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("sa", "show-additional-intervals", "The intervals that are in the Sam with good score but not in WIT file to stderr if set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("mN", "match-N", "If set, N matches all as a wildcard character, otherwise it never matches.", OptionType::Boolean));
    addOption(parser, CommandLineOption("wm", "weighted-distances", "If set, use weighted distances instead of unit ones.", OptionType::Boolean));
    addOption(parser, CommandLineOption("c", "benchmark-category", "The benchmark category to compare for.  One of {all, any-best, all-best}.  Default: all.", OptionType::String));
    addOption(parser, CommandLineOption("DP", "DONT-PANIC", "Don't panic on additional hits in non-weightedmode.  Default: False.", OptionType::Boolean));
    addOption(parser, CommandLineOption("ow", "oracle-wit-mode", "Enable 'oracle wit mode', alignment distances are ignored, assumed to be 0.  Default: False.", OptionType::Boolean));
    
    // We require 4 command line options.
    requiredArguments(parser, 3);
}


// Actually parse the arguments.  Return kRetOk or the return value
// for main.
int parseCommandLineAndCheck(Options &options,
                             CommandLineParser &parser,
                             CharString &outFile,
                             const int argc, const char *argv[]) {
    // Show short help on invalid arguments and long help if help
    // argument was given.
    if (not parse(parser, argc, argv)) {
        return kRetArgsErr;
    } else if (isSetShort(parser, 'h')) {
        exit(kRetOk);
    }

    // Get arguments.
    if (isSetLong(parser, "DONT-PANIC"))
        options.dontPanic = true;
    if (isSetLong(parser, "oracle-wit-mode"))
        options.oracleWitMode = true;
    if (isSetLong(parser, "match-N"))
        options.matchN = true;
    if (isSetLong(parser, "weighted-distances"))
        options.weightedDistances = true;
    if (isSetLong(parser, "show-missed-intervals"))
        options.showMissedIntervals = true;
    if (isSetLong(parser, "show-hit-intervals"))
        options.showHitIntervals = true;
    if (isSetLong(parser, "show-try-hit-intervals"))
        options.showTryHitIntervals = true;
    if (isSetLong(parser, "show-superflous-intervals"))
        options.showSuperflousIntervals = true;
    if (isSetLong(parser, "show-additional-intervals"))
        options.showAdditionalIntervals = true;
    if (isSetLong(parser, "max-error-rate"))
        getOptionValueLong(parser, "max-error-rate", options.maxError);
    if (isSetLong(parser, "distance-function"))
        getOptionValueLong(parser, "distance-function", options.distanceFunction);
    if (isSetLong(parser, "out-file"))
        getOptionValueLong(parser, "out-file", outFile);
    if (isSetLong(parser, "benchmark-category"))
        getOptionValueLong(parser, "benchmark-category", options.benchmarkCategory);

    // Validate values.
    if (options.maxError < 0) {
        std::cerr << "ERROR: Invalid maximum error value: " << options.maxError << std::endl;
        return kRetArgsErr;
    }
    if (!options.validDistanceFunction()) {      
      std::cerr << "ERROR: Invalid distance function: " << options.distanceFunction << std::endl;
      return kRetArgsErr;
    }
    if (!options.validBenchmarkCategory()) {
        std::cerr << "ERROR: Invalid benchmark category: " << options.benchmarkCategory << std::endl;
        return kRetArgsErr;
    }
    
    // Get positional arguments.
    options.seqFileName = getArgumentValue(parser, 0);
    options.samFileName = getArgumentValue(parser, 1);
    options.witFileName = getArgumentValue(parser, 2);

    return kRetOk;
}


int main(int argc, const char *argv[]) {
    // =================================================================
    // Parse Options.
    // =================================================================
    // Initialize Options object with defaults.
    Options options;
    options.maxError = 0;
    options.distanceFunction = "edit";
    options.showMissedIntervals = false;
    options.showHitIntervals = false;
    options.showTryHitIntervals = false;
    options.matchN = false;
    options.weightedDistances = false;
    options.showSuperflousIntervals = false;
    options.showAdditionalIntervals = false;
    options.benchmarkCategory = "all";
    options.dontPanic = false;
    options.oracleWitMode = false;
    CharString outFile = "-";

    // Setup the parser, parse command line and return if an error occured.
    CommandLineParser parser;
    setUpCommandLineParser(parser);
    int ret = parseCommandLineAndCheck(options, parser, outFile, argc, argv);
    if (ret != kRetOk)
        return ret;

    // =================================================================
    // Load FASTA Sequence And Sam File Into FragmentStore.
    // =================================================================
    typedef FragmentStore<> TFragmentStore;
    TFragmentStore fragments;

    // Load Contigs.
    double startTime = sysTime();
    std::cerr << "Reading FASTA contigs sequence file " << options.seqFileName << " ..." << std::endl;
    if (!loadContigs(fragments, options.seqFileName)) {
        std::cerr << "Could not read contigs." << std::endl;
        return kRetIoErr;
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // Load Sam File.
    std::cerr << "Reading Sam file file " << options.samFileName << " ..." << std::endl;
    startTime = sysTime();
    {
        std::fstream fstrm(toCString(options.samFileName),
                           std::ios_base::in | std::ios_base::binary);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open Sam file." << std::endl;
            return kRetIoErr;
        }
        read(fstrm, fragments, Sam());
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;
    //for (unsigned i = 0; i < length(fragments.readNameStore); ++i) {
      //std::cerr << ">>>" << fragments.readNameStore[i] << " " << i << std::endl;
    //}

    // =================================================================
    // Load WIT file.
    // =================================================================
    std::cerr << "Loading intervals from " << options.witFileName << std::endl;
    startTime = sysTime();
    WitStore witStore;
    loadWitFile(witStore, fragments, options.witFileName);
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;
    
    // =================================================================
    // Perform Interval Score Lowering.
    // =================================================================
    performIntervalScoreLowering(witStore, options.maxError);

    // =================================================================
    // Compare The Sam Hits Against WIT Intervals.
    // =================================================================
    std::cerr << "Compare reader hits from Sam file against WIT file." << std::endl;
    startTime = sysTime();
    typedef Position<WitStore::TIntervalStore>::Type TPos;
    // The result will be a list of ids to entries in witStore.
    String<size_t> result;
    if (options.distanceFunction == "edit")
        compareAlignedReadsToReference(result, fragments, witStore, options, MyersUkkonenReads());
    else  // options.distanceFunction == "hamming"
        compareAlignedReadsToReference(result, fragments, witStore, options, HammingSimple());
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // =================================================================
    // Perform Counting On Result Indices, Yield ComparisonResult.
    // =================================================================
    ComparisonResult comparisonResult;
    evaluateFoundIntervals(comparisonResult, witStore, result, fragments, options);

    // =================================================================
    // Write Output.
    // =================================================================
    startTime = sysTime();
    // The output consists of one line that describes the total and
    // found intervals as a JSON record with the entries
    // "total_intervals", "found_itervals", "superflous_intervals",
    // "additional_intervals".
    if (outFile == "-") {
        // Print to stdout.
        std::cout << comparisonResult << std::endl;
    } else {
        // Write output to file.
        std::fstream fstrm(toCString(outFile), std::ios_base::out);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open output JSON file." << std::endl;
            return kRetIoErr;
        }
        fstrm << comparisonResult << std::endl;
    }

    return kRetOk;
}
