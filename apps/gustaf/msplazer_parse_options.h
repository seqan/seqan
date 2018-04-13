// ==========================================================================
//                                  Gustaf
// ==========================================================================
// Copyright (c) 2011-2018, Kathrin Trappe, FU Berlin
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
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_APPS_GUSTAF_MSPLAZER_PARSEOPTIONS_H_
#define SEQAN_APPS_GUSTAF_MSPLAZER_PARSEOPTIONS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include "msplazer.h"

using namespace seqan;

// /////////////////////////////////////////////////////////////////////////////
// Parses options from command line parser and writes them into options object
bool
_parseOptions(ArgumentParser & parser, StellarOptions & options, MSplazerOptions & msplazerOptions)
{
// IOREV _notio_
    // i/o options
    getArgumentValue(msplazerOptions.databaseFile, parser, 0);
    // getArgumentValue(msplazerOptions.queryFile, parser, 1);
    resize(msplazerOptions.queryFile, getArgumentValueCount(parser, 1), Exact());
    for (unsigned i = 0; i < length(msplazerOptions.queryFile); ++i)
        getArgumentValue(msplazerOptions.queryFile[i], parser, 1, i);
    for (unsigned i = 0; i < length(msplazerOptions.queryFile); ++i)
        std::cout << msplazerOptions.queryFile[i] << std::endl;
    if (length(msplazerOptions.queryFile) > 1 && !isSet(parser, "m"))
    {
        std::cerr << "Please provide an input (Stellar) match file when using paired-end mode!" << std::endl;
        return false;
    }
    if (length(msplazerOptions.queryFile) > 1 && !isSet(parser, "ll"))
    {
        std::cerr << "Please provide library length when using paired-end mode!" << std::endl;
        return false;
    }
    if (length(msplazerOptions.queryFile) > 1 && !isSet(parser, "le"))
    {
        std::cerr << "Please provide library error when using paired-end mode!" << std::endl;
        return false;
    }
    // getOptionValue(msplazerOptions.queryFile2, parser, "q2");
    // getOptionValue(msplazerOptions.outDir, parser, "i");
    getOptionValue(msplazerOptions.vcfOutFile, parser, "vcf");
    getOptionValue(msplazerOptions.gffOutFile, parser, "gff");
    getOptionValue(msplazerOptions.jobName, parser, "j");
    getOptionValue(msplazerOptions.stellarInputFile, parser, "m");
    getOptionValue(msplazerOptions.dotOut, parser, "do");

    // main options
    getOptionValue(msplazerOptions.diffDBPen, parser, "tp");
    getOptionValue(msplazerOptions.diffStrandPen, parser, "ip");
    getOptionValue(msplazerOptions.diffOrderPen, parser, "op");
    getOptionValue(msplazerOptions.simThresh, parser, "oth");
    getOptionValue(msplazerOptions.gapThresh, parser, "gth");
    getOptionValue(msplazerOptions.initGapThresh, parser, "ith");
    getOptionValue(msplazerOptions.breakendThresh, parser, "bth");
    getOptionValue(msplazerOptions.tandemThresh, parser, "tth");
    getOptionValue(msplazerOptions.breakpointPosRange, parser, "pth");
    if (isSet(parser, "cbp"))    
        msplazerOptions.inferComplexBP = false;
    getOptionValue(msplazerOptions.support, parser, "st");
    getOptionValue(msplazerOptions.mateSupport, parser, "mst");
    getOptionValue(msplazerOptions.libSize, parser, "ll");
    getOptionValue(msplazerOptions.libError, parser, "le");
    if (isSet(parser, "rc"))
        msplazerOptions.revCompl = false;

    if (length(msplazerOptions.queryFile) > 1)
        msplazerOptions.pairedEndMode = true;
    getOptionValue(msplazerOptions.numThreads, parser, "nth");

    // Parsing Stellar options
    getArgumentValue(options.databaseFile, parser, 0);
    //getArgumentValue(options.queryFile, parser, 1);
    options.queryFile = msplazerOptions.queryFile[0];

    getOptionValue(options.qGram, parser, "kmer");
    getOptionValue(options.minLength, parser, "minLength");
    getOptionValue(options.epsilon, parser, "epsilon");
    getOptionValue(options.xDrop, parser, "xDrop");
    getOptionValue(options.alphabet, parser, "alphabet");

    if (isSet(parser, "forward") && !isSet(parser, "reverse"))
        options.reverse = false;
    if (!isSet(parser, "forward") && isSet(parser, "reverse"))
        options.forward = false;

    getOptionValue(options.fastOption, parser, "verification");
    getOptionValue(options.disableThresh, parser, "disableThresh");
    getOptionValue(options.numMatches, parser, "numMatches");
    getOptionValue(options.compactThresh, parser, "sortThresh");
    getOptionValue(options.maxRepeatPeriod, parser, "repeatPeriod");
    getOptionValue(options.minRepeatLength, parser, "repeatLength");
    getOptionValue(options.qgramAbundanceCut, parser, "abundanceCut");

    getOptionValue(options.verbose, parser, "verbose");

    if (isSet(parser, "kmer") && options.qGram >= 1 / options.epsilon)
    {
        std::cerr << "Invalid parameter value: Please choose q-gram length lower than 1/epsilon." << std::endl;
        return false;
    }

    if (options.numMatches > options.compactThresh)
    {
        std::cerr << "Invalid parameter values: Please choose numMatches <= sortThresh." << std::endl;
        return false;
    }
    return true;
}

// /////////////////////////////////////////////////////////////////////////////
// Set-Up of Argument Parser
void _setupArgumentParser(ArgumentParser & parser)
{
    setShortDescription(
        parser,
        "Gustaf - Generic mUlti-SpliT Alignment Finder: Tool for split-read mapping allowing multiple splits.");
    setVersion(parser, "1.0.0");
    setDate(parser, "August 2014");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FASTA FILE\\fP> <\\fIREAD FASTA FILE\\fP> \n "
                 );
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIGENOME FASTA FILE\\fP> <\\fIREAD FASTA FILE\\fP> <\\fIREAD FASTA FILE 2\\fP> \n "
                 );
    // addUsageLine(parser, "-d <FASTA sequence file> -q <FASTA sequence file> [\\fIOPTIONS\\fP]");
    addDescription(parser,
                   "GUSTAF uses SeqAns STELLAR to find splits as local matches on different strands or "
                   "chromosomes. Criteria and penalties to chain these matches can be specified. "
                   "Output file contains the breakpoints along the best chain.");
    addDescription(parser, "The genome file is used as database input, the read file as query input.");
    addDescription(parser,
                   "All STELLAR options are supported. See STELLAR documentation for STELLAR parameters and options.");
    addDescription(parser, "(c) 2011-2012 by Kathrin Trappe");

    addSection(parser, "GUSTAF Options");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILE 1"));
    setValidValues(parser, 0, "fa fasta fq fastq");  // allow only fasta/q files as input
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILE 2", true));
    setValidValues(parser, 1, "fa fasta fq fastq");  // allow only fasta/q files as input
    setHelpText(parser, 1, "Either one (single-end) or two (paired-end) read files.");


    addSection(parser, "Main Options");  // addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILE 3"));
    // setValidValues(parser, 2, "fa fasta");  // allow only fasta files as input

    addOption(parser,
              ArgParseOption("tp", "transPen", "Interchromosomal translocation penalty", ArgParseArgument::INTEGER,
                             "INT"));
    setDefaultValue(parser, "tp", "5");
    addOption(parser, ArgParseOption("ip", "invPen", "Inversion penalty", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "ip", "5");
    addOption(parser,
              ArgParseOption("op", "orderPen", "Intrachromosomal order change penalty", ArgParseArgument::INTEGER,
                             "INT"));
    setDefaultValue(parser, "op", "0");
    addOption(parser,
              ArgParseOption("oth", "overlapThresh", "Allowed overlap between matches", ArgParseArgument::DOUBLE,
                             "DOUBLE"));
    setDefaultValue(parser, "oth", "0.5");
    addOption(parser,
              ArgParseOption("gth", "gapThresh", "Allowed gap length between matches, default value corresponse to expected size of microindels (5 bp)", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "gth", "5");
    addOption(parser, ArgParseOption(
                  "ith", "initGapThresh", "Allowed initial or ending gap length at begin and end of read with no breakpoint (e.g. due to sequencing errors at the end)",
                  ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "ith", "15");
    addOption(parser, ArgParseOption(
                  "bth", "breakendThresh", "Allowed initial or ending gap length at begin and end of read that creates a breakend/breakpoint (e.g. for reads extending into insertions)",
                  ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "bth", "30");
    addOption(parser, ArgParseOption(
                  "tth", "tandemThresh", "Minimal length of (small) insertion/duplication with double overlap to be considered tandem repeat",
                  ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "tth", "50");
    addOption(parser, ArgParseOption(
                  "pth", "breakpoint-pos-range", "Allowed difference in breakpoint position", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "pth", "5");
    addOption(parser, ArgParseOption("cbp", "complex-breakpoints", "Disable inferring complex SVs"));
    addOption(parser, ArgParseOption("st", "support", "Number of supporting reads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "st", "2");
    addOption(parser, ArgParseOption("mst", "mate-support", "Number of supporting concordant mates", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "mst", "2");
    addOption(parser, ArgParseOption("ll", "library-size", "Library size of paired-end reads", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("le", "library-error", "Library error (sd) of paired-end reads", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("rc", "revcompl", "Disable reverse complementing second mate pair input file."));
    // set min values?

    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("m", "matchfile", "File of (stellar) matches", ArgParseArgument::INPUT_FILE, "FILE"));
    setValidValues(parser, "m", "gff GFF");

    addSection(parser, "Output Options");
    addOption(parser,
              ArgParseOption("gff", "gffOut", "Name of gff breakpoint output file.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "gff", "gff txt");
    setDefaultValue(parser, "gff", "breakpoints.gff");
    addOption(parser,
              ArgParseOption("vcf", "vcfOut", "Name of vcf breakpoint output file.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "vcf", "vcf txt");
    setDefaultValue(parser, "vcf", "breakpoints.vcf");

    addOption(parser, ArgParseOption("j", "jobName", "Job/Queue name", ArgParseArgument::STRING, "STR"));
    setDefaultValue(parser, "j", "");

    addOption(parser, ArgParseOption("do", "dots", "Enable graph output in dot format"));

    addSection(parser, "Parallelization Options");
    addOption(parser,
              ArgParseOption("nth", "numThreads", "Number of threads for parallelization of I/O.", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "nth", "1");

    addSection(parser, "Stellar Options");

    addSection(parser, "Main Options");

    addOption(parser, ArgParseOption("e", "epsilon", "Maximal error rate (max 0.25).", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "e", "0.05");
    setMinValue(parser, "e", "0.0000001");
    setMaxValue(parser, "e", "0.25");
    addOption(parser, ArgParseOption("l", "minLength", "Minimal length of epsilon-matches.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "l", "100");
    setMinValue(parser, "l", "0");
    addOption(parser, ArgParseOption("f", "forward", "Search only in forward strand of database."));
    addOption(parser, ArgParseOption("r", "reverse", "Search only in reverse complement of database."));
    addOption(parser, ArgParseOption("a", "alphabet",
                                     "Alphabet type of input sequences (dna, rna, dna5, rna5, protein, char).",
                                     ArgParseArgument::STRING));
    setValidValues(parser, "a", "dna dna5 rna rna5 protein char");
    addOption(parser, ArgParseOption("v", "verbose", "Set verbosity mode."));

    addSection(parser, "Filtering Options");

    addOption(parser, ArgParseOption("k", "kmer", "Length of the q-grams (max 32).", ArgParseArgument::INTEGER));
    setMinValue(parser, "k", "1");
    setMaxValue(parser, "k", "32");
    addOption(parser, ArgParseOption("rp", "repeatPeriod",
                                     "Maximal period of low complexity repeats to be filtered.",
                                     ArgParseArgument::INTEGER));
    setDefaultValue(parser, "rp", "1");
    addOption(parser, ArgParseOption("rl", "repeatLength",
                                     "Minimal length of low complexity repeats to be filtered.",
                                     ArgParseArgument::INTEGER));
    setDefaultValue(parser, "rl", "1000");
    addOption(parser, ArgParseOption("c", "abundanceCut", "k-mer overabundance cut ratio.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "c", "1");
    setMinValue(parser, "c", "0");
    setMaxValue(parser, "c", "1");

    addSection(parser, "Verification Options");

    addOption(parser, ArgParseOption("x", "xDrop", "Maximal x-drop for extension.", ArgParseArgument::DOUBLE));
    setDefaultValue(parser, "x", "5");
    addOption(parser, ArgParseOption("vs", "verification", "Verification strategy: exact or bestLocal or bandedGlobal",
                                     ArgParseArgument::STRING));
    // addHelpLine(parser, "exact        = compute and extend all local alignments in SWIFT hits");
    // addHelpLine(parser, "bestLocal    = compute and extend only best local alignment in SWIFT hits");
    // addHelpLine(parser, "bandedGlobal = banded global alignment on SWIFT hits");
    setDefaultValue(parser, "vs", "exact");
    setValidValues(parser, "vs", "exact bestLocal bandedGlobal");
    addOption(parser, ArgParseOption("dt", "disableThresh",
                                     "Maximal number of verified matches before disabling verification for one query "
                                     "sequence (default infinity).", ArgParseArgument::INTEGER));
    setMinValue(parser, "dt", "0");
    addOption(parser, ArgParseOption("n", "numMatches",
                                     "Maximal number of kept matches per query and database. If STELLAR finds more matches, "
                                     "only the longest ones are kept.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "n", "50");
    addOption(parser, ArgParseOption("s", "sortThresh",
                                     "Number of matches triggering removal of duplicates. Choose a smaller value for saving "
                                     "space.", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "s", "500");

    /*
     * Stellar output options are not supported bc. no Stellar output is supported
    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "out", "Name of output file.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "o", "gff txt");
    setDefaultValue(parser, "o", "stellar.gff");
    addOption(parser, ArgParseOption("od", "outDisabled",
                                     "Name of output file for disabled query sequences.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "outDisabled", "fa FASTA");
    setDefaultValue(parser, "od", "stellar.disabled.fasta");
    addOption(parser, ArgParseOption("t", "no-rt", "Suppress printing running time."));
    hideOption(parser, "t");
    */
}

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_MSPLAZER_PARSEOPTIONS_
