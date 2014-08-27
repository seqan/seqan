// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "anise_app.h"

#include <chrono>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/vcf_io.h>

#include "shared/progress_indicator.h"

#include "anise/anise_options.h"
#include "anise/app_state.h"
#include "anise/assembly_substep.h"
#include "anise/finishing_step.h"
#include "anise/preparation_step.h"
#include "anise/read_mapping_substep.h"
#include "anise/site_state.h"
#include "anise/temporary_file_manager.h"
#include "anise/time_log.h"

namespace  // anonymous namespace
{

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

// Parse options from command line.  Result can be one of OK, ERROR, OTHER.

seqan::ArgumentParser::ParseResult parseCommandLine(AniseOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("anise");

    // Set short description, version, and date

    setShortDescription(parser, "Assembly of Novel Inserted SEquence");
#ifdef SEQAN_REVISION
        setVersion(parser, "0.2.0-beta.1 [" + std::string(SEQAN_REVISION) + "]");
#else
        setVersion(parser, "0.2.0-beta.1");
#endif
#ifdef SEQAN_DATE
        setDate(parser, SEQAN_DATE);
#endif

    // Define usage line and long description

    addUsageLine(parser,
                 "\\fB-ir\\fP \\fIIN.fa\\fP \\fB-iv\\fP \\fIIN.vcf\\fP \\fB-im\\fP \\fIIN.{sam,bam}\\fP "
                 "\\fB-of\\fP \\fIOUT.fa\\fP [\\fB-om\\fP \\fIOUT.sam\\fP]");
    addDescription(parser,
                   "Read insertion site candidates from \\fIIN.vcf\\fP and the mapping from \\fIIN.{sam,bam}\\fP.  "
                   "\\fIANISE\\fP will then try to assemble the inserted sequences at the sites in \\fIIN.vcf\\fP and "
                   "write the assembled sequences to \\fIOUT.vcf\\fP.");
    addDescription(parser,
                   "The reference sequence is taken from \\fIIN.fa\\fP");

    // General Options

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addOption(parser, seqan::ArgParseOption("", "num-threads",
                                            "Number of threads to use.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "num-threads", "1");
    setDefaultValue(parser, "num-threads", "1");

    addOption(parser, seqan::ArgParseOption("", "debug-site-id",
                                            "Debug site ID (-1 to disable).",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "debug-site-id", "-1");
    setDefaultValue(parser, "debug-site-id", "-1");

    addOption(parser, seqan::ArgParseOption("", "debug-step-no",
                                            "Debug step no (-1 to disable).",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "debug-step-no", "-1");
    setDefaultValue(parser, "debug-step-no", "-1");

    addOption(parser, seqan::ArgParseOption("", "no-auto-tuning", "Disable auto-tuning (see below)."));

    // I/O Related Options

    addSection(parser, "Input / Output");

    addOption(parser, seqan::ArgParseOption("ir", "input-reference", "Input FASTA file with reference.",
                                            seqan::ArgParseOption::INPUTFILE, "IN.fa"));
    setValidValues(parser, "input-reference", "fa fasta");
    setRequired(parser, "input-reference");

    addOption(parser, seqan::ArgParseOption("iv", "input-vcf", "Input VCF file with insert site candidates.",
                                            seqan::ArgParseOption::INPUTFILE, "IN.vcf"));
    setValidValues(parser, "input-vcf", "vcf");
    setRequired(parser, "input-vcf");

    addOption(parser, seqan::ArgParseOption("im", "input-mapping", "Input SAM/BAM mapping file.",
                                            seqan::ArgParseOption::INPUTFILE, "IN.{sam,bam}"));
    setValidValues(parser, "input-mapping", "sam bam");
    setRequired(parser, "input-mapping");

    addOption(parser, seqan::ArgParseOption("of", "output-fasta", "Output FASTA with contigs",
                                            seqan::ArgParseOption::OUTPUTFILE, "OUT.fa"));
    setValidValues(parser, "output-fasta", "fa fasta");
    setRequired(parser, "output-fasta");

    addOption(parser, seqan::ArgParseOption("om", "output-mapping", "Output SAM/BAM file with mapping fo reads "
                                            "to contigs in \\fB--output-fasta\\fP.",
                                            seqan::ArgParseOption::OUTPUTFILE, "OUT.sam"));
    setValidValues(parser, "output-mapping", "sam bam");

    addOption(parser, seqan::ArgParseOption("", "output-debug-dir",
                                            "Directory for debug output.  Leave empty for no such output.",
                                            seqan::ArgParseOption::STRING, "OUT.vcf"));

    addOption(parser, seqan::ArgParseOption("", "clean-up-tmp-files", "Clean up temporary files."));

    // Algorithm Options

    addSection(parser, "Algorithm");

    addOption(parser, seqan::ArgParseOption("", "recursion-max-steps",
                                            "Maximal recursion depth.  0 for infinity.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "recursion-max-steps", "0");
    setDefaultValue(parser, "recursion-max-steps", "50");

    addOption(parser, seqan::ArgParseOption("", "no-realign-assembly",
                                            "Do not realign the reads after assembly."));

    addOption(parser, seqan::ArgParseOption("", "max-reads-factor",
                                            "Factor to use for the maximal read computation.  ANISE stops for a site "
                                            "if more than the number of reads expected from the expected coverage "
                                            "times the max reads factor are assigned to the site.",
                                            seqan::ArgParseOption::DOUBLE, "FLOAT"));
    setMinValue(parser, "max-reads-factor", "1");
    setDefaultValue(parser, "max-reads-factor", "2");

    addOption(parser, seqan::ArgParseOption("", "stop-initial-read-count",
                                            "If there are more than this number of reads for a site in the initial round "
                                            "then no assembly is performed.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "stop-initial-read-count", "0");
    setDefaultValue(parser, "stop-initial-read-count", "4000");

    addOption(parser, seqan::ArgParseOption("", "stop-tex-read-count",
                                            "If there are more than this number of reads for a site in a later round "
                                            "then no triplet library extension is performed.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "stop-tex-read-count", "0");
    setDefaultValue(parser, "stop-tex-read-count", "3000");

    addOption(parser, seqan::ArgParseOption("", "stop-read-count",
                                            "If there are more than this number of reads for a site in a later round "
                                            "then no assembly is performed.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "stop-read-count", "0");
    setDefaultValue(parser, "stop-read-count", "30000");

    addOption(parser, seqan::ArgParseOption("", "stop-coverage",
                                            "If the length sum of all reads for a site divided by the length sum of "
                                            "its contigs is higher than this value before assembly then the site is "
                                            "deactivated.  Set to 0 to deactivate check.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "stop-coverage", "0");
    setDefaultValue(parser, "stop-coverage", "100");

    addOption(parser, seqan::ArgParseOption("", "realignment-bandwidth", "The bandwidth to use in the realignment step.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "realignment-bandwidth", "0");
    setDefaultValue(parser, "realignment-bandwidth", "40");

    addOption(parser, seqan::ArgParseOption("", "realignment-border", "The border from the profile to extract around alignments.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "realignment-border", "0");
    setDefaultValue(parser, "realignment-border", "30");

    // Read Separation Options

    addSection(parser, "Repeat Separation");

    addOption(parser, seqan::ArgParseOption("", "no-separate-repeats", "Dont' repeat separation algorithm after realignment."));

    addOption(parser, seqan::ArgParseOption("", "repsep-tammi-method", "Variant of the Tammi method to use for repeat "
                                            "separation (simple or phred).", seqan::ArgParseOption::STRING, "STR"));
    setValidValues(parser, "repsep-tammi-method", "phred simple");
    setDefaultValue(parser, "repsep-tammi-method", "simple");

    addOption(parser, seqan::ArgParseOption("", "repsep-p-err", "Repeat separation per-base error for simple Tammi method.",
                                            seqan::ArgParseOption::DOUBLE, "FLOAT"));
    setMinValue(parser, "repsep-p-err", "0.0");
    setMaxValue(parser, "repsep-p-err", "1.0");
    setDefaultValue(parser, "repsep-p-err", "0.01");

    addOption(parser, seqan::ArgParseOption("", "repsep-max-random-correlation", "Repeat separation maximal random correlation.",
                                            seqan::ArgParseOption::DOUBLE, "FLOAT"));
    setMinValue(parser, "repsep-max-random-correlation", "0.0");
    setMaxValue(parser, "repsep-max-random-correlation", "1.0");
    setDefaultValue(parser, "repsep-max-random-correlation", "0.00001");

    addOption(parser, seqan::ArgParseOption("", "repsep-tau-min", "Repeat separation tau_min value.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "repsep-tau-min", "0");
    setDefaultValue(parser, "repsep-tau-min", "2");
    setDefaultValue(parser, "repsep-tau-min", "100000");  // single-column separation sites

    addOption(parser, seqan::ArgParseOption("", "repsep-r-min", "Repeat separation r_min value.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "repsep-r-min", "0");
    setDefaultValue(parser, "repsep-r-min", "2");
    setDefaultValue(parser, "repsep-r-min", "100000");  // single-column separation sites

    addOption(parser, seqan::ArgParseOption("", "repsep-min-overlap", "Repeat separation minimal overlap value.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "repsep-min-overlap", "0");
    setDefaultValue(parser, "repsep-min-overlap", "2");

    addOption(parser, seqan::ArgParseOption("", "repsep-start-compression-at", "Repeat separation start compression.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "repsep-start-compression-at", "2");
    setDefaultValue(parser, "repsep-start-compression-at", "100");

    addOption(parser, seqan::ArgParseOption("", "repsep-split-d-min", "Repeat separation split at d_min deviations."));

    // Paired-End Options

    addSection(parser, "Library Info");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-factor",
                                            "Factor to multiple fragment size stddev with to get allowed error.",
                                            seqan::ArgParseOption::DOUBLE, "FACTOR"));
    setMinValue(parser, "fragment-size-factor", "0");
    setDefaultValue(parser, "fragment-size-factor", "8");

    addOption(parser, seqan::ArgParseOption("", "auto-library-num-records",
                                            "Number of records to use for automatic library evaluation.  Set to 0 to "
                                            "evaluate all.", seqan::ArgParseOption::INTEGER, "NUM"));
    setMinValue(parser, "auto-library-num-records", "0");
    setDefaultValue(parser, "auto-library-num-records", "100000");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-median",
                                            "Median fragment size.", seqan::ArgParseOption::DOUBLE, "SIZE"));
    setMinValue(parser, "fragment-size-median", "0");
    setDefaultValue(parser, "fragment-size-median", "250");

    addOption(parser, seqan::ArgParseOption("", "fragment-size-std-dev", "Fragment size standard deviation.",
                                            seqan::ArgParseOption::DOUBLE, "STDDEV"));
    setMinValue(parser, "fragment-size-std-dev", "0");
    setDefaultValue(parser, "fragment-size-std-dev", "30");

    addOption(parser, seqan::ArgParseOption("", "fragment-default-orientation", "Default orientation.",
                                            seqan::ArgParseOption::STRING, "FACTOR"));
    setValidValues(parser, "fragment-default-orientation", "F+ F- R+ R-");
    setDefaultValue(parser, "fragment-default-orientation", "R+");

    // Assembly Options

    addSection(parser, "Assembly");

    addOption(parser, seqan::ArgParseOption("", "assembly-site-window-radius",
                                            "Radius around insert site to cut for initial contigs.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "assembly-site-window-radius", "100");
    setDefaultValue(parser, "assembly-site-window-radius", "1000");
    addOption(parser, seqan::ArgParseOption("", "assembly-site-fringe-radius",
                                            "Radius around insert site to cut for collecting clippings.  Set to -1 "
                                            "(default) to consider all records with >= 15 clipped bases.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "assembly-site-fringe-radius", "-1");
    setDefaultValue(parser, "assembly-site-fringe-radius", "-1");

    // Overlapper Read Mapping

    addSection(parser, "Read Mapping");

    addOption(parser, seqan::ArgParseOption("", "read-mapping-error-rate",
                                            "Error rate of internal read mapping step in percent.",
                                            seqan::ArgParseOption::DOUBLE, "PERCENT"));
    setMinValue(parser, "read-mapping-error-rate", "0");
    setMaxValue(parser, "read-mapping-error-rate", "20");
    setDefaultValue(parser, "read-mapping-error-rate", "5");

    addOption(parser, seqan::ArgParseOption("", "read-mapping-batch-size",
                                            "Batch size for read mapping.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "read-mapping-batch-size", "1");
    setDefaultValue(parser, "read-mapping-batch-size", "10000");

    // Overlapper Options

    addSection(parser, "Overlapper");

    addOption(parser, seqan::ArgParseOption("", "overlapper-min-overlap-ratio",
                                            "Overlapper min overlap rate in percent of the longer read.",
                                            seqan::ArgParseOption::INTEGER, "PERCENT"));
    setMinValue(parser, "overlapper-min-overlap-ratio", "0");
    setDefaultValue(parser, "overlapper-min-overlap-ratio", "40");

    addOption(parser, seqan::ArgParseOption("", "overlapper-max-error-rate",
                                            "Overlapper maximum error rate in percent.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setMinValue(parser, "overlapper-max-error-rate", "0");
    setMaxValue(parser, "overlapper-max-error-rate", "30");
    setDefaultValue(parser, "overlapper-max-error-rate", "5");

    // MSA Options

    addSection(parser, "Multiple Sequence Alignment");

    addOption(parser, seqan::ArgParseOption("", "no-read-correction", "Whether or not to perform read correction"));

    addOption(parser, seqan::ArgParseOption("", "msa-score-match", "PW match score in MSA.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "msa-score-match", "2");

    addOption(parser, seqan::ArgParseOption("", "msa-score-mismatch", "PW mismatch score in MSA.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "msa-score-mismatch", "-6");

    addOption(parser, seqan::ArgParseOption("", "msa-score-gap-open", "PW gap open score in MSA.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "msa-score-gap-open", "-4");

    addOption(parser, seqan::ArgParseOption("", "msa-score-gap-extend", "PW gap extension score in MSA.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "msa-score-gap-extend", "-9");

    // Consensus Calling Options

    addSection(parser, "Consensus Calling");

    addOption(parser, seqan::ArgParseOption("", "consensus-min-base-support",
                                            "Minimal base support for non-N call in consensus calling.",
                                            seqan::ArgParseOption::INTEGER, "INT"));
    setDefaultValue(parser, "consensus-min-base-support", "2");

    addOption(parser, seqan::ArgParseOption("", "consensus-min-contig-length-rate",
                                            "Minimal contig length in percent of average read length.",
                                            seqan::ArgParseOption::INTEGER, "PERCENT"));
    setDefaultValue(parser, "consensus-min-contig-length-rate", "150");

    // Adding section on library property detection.
    addTextSection(parser, "Library Properties");
    addText(parser,
            "The terms insert size, fragment, and template length all denote the length of the physical fragment "
            "that was extracted and is then sequenced from both sides to yield paired reads.");
    addText(parser,
            "Note that if you set \\fB--fragment-size-mean\\fP or \\fB--fragment-size-std-dev\\fP then you have "
            "to set both.");

    addTextSection(parser, "Repeat Separation");
    addText(parser,
            "In each step, ANISE can try to separate repeats after the realignment step.  This can reduce "
            "problems such as \"chimeric assemblies\", i.e. when an insert is assembled from reads that come "
            "from different copies of the repeat.");

    addTextSection(parser, "Auto Tuning");
    addText(parser,
            "ANISE expresses the minimal overlap and maximal allowed error in MSA computation and the allowed "
            "errors in terms of rates (%) of the average read length.  This works well for longer Illumina reads "
            "but for short reads, one should require an overlap of 14 bp and allow up to 2 errors in the overlap. "
            "Likewise, read mapping should allow at least up to 2 errors.");
    addText(parser,
            "Unless \\fB--no-auto-tuning\\fP is specified, ANISE will adjust the settings to these values.");

    // References Text Section

    addTextSection(parser, "References");
    addText(parser,
            "Hajirasouliha, I., Hormozdiari, F., Alkan, C., Kidd, J.M., Birol, I., Eichler, E.E., Sahinalp, S.C.  "
            "Detection and characterization of sequence insertions using paired-end next-generation sequencing.  "
            "Bioinformatics 2010 May; 15;26(10):1277-83.");

    // Parse command line
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // General Options

    options.verbosity = 1;
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.numThreads, parser, "num-threads");
    getOptionValue(options.debugSiteID, parser, "debug-site-id");
    getOptionValue(options.debugStepNo, parser, "debug-step-no");
    bool noAutoTuning = options.autoTuning;
    getOptionValue(noAutoTuning, parser, "no-auto-tuning");
    options.autoTuning = !noAutoTuning;

    // I/O Related Options

    getOptionValue(options.inputReference, parser, "input-reference");
    getOptionValue(options.inputVcf, parser, "input-vcf");
    getOptionValue(options.inputMapping, parser, "input-mapping");
    getOptionValue(options.outputFasta, parser, "output-fasta");
    getOptionValue(options.outputMapping, parser, "output-mapping");

    getOptionValue(options.outputDebugDir, parser, "output-debug-dir");
    getOptionValue(options.cleanUpTemporaryFiles, parser, "clean-up-tmp-files");
    // TODO(holtgrew): Allow specification of tmpdir or use TEMPDIR environment variable.
    getOptionValue(options.tmpDir, parser, "output-fasta");
    append(options.tmpDir, ".tmp");

    // Algorithm Related Options

    getOptionValue(options.recursionMaxSteps, parser, "recursion-max-steps");
    options.realignAssembly = !isSet(parser, "no-realign-assembly");
    getOptionValue(options.stopInitialReadCount, parser, "stop-initial-read-count");
    getOptionValue(options.stopReadCount, parser, "stop-read-count");
    getOptionValue(options.stopTripletExtensionReadCount, parser, "stop-tex-read-count");
    getOptionValue(options.stopCoverage, parser, "stop-coverage");
    getOptionValue(options.realignmentBandwidth, parser, "realignment-bandwidth");
    getOptionValue(options.realignmentBorder, parser, "realignment-border");

    // Repeat separation options.

    bool noSeparateRepeats = false;
    getOptionValue(noSeparateRepeats, parser, "no-separate-repeats");
    options.separateRepeats = !noSeparateRepeats;
    {
        seqan::CharString tmp;
        getOptionValue(tmp, parser, "repsep-tammi-method");
        if (tmp == "simple")
            options.readSepOptions.tammiMethod = rep_sep::ReadSeparatorOptions::TAMMI_SIMPLE;
        else
            options.readSepOptions.tammiMethod = rep_sep::ReadSeparatorOptions::TAMMI_PHRED;
        getOptionValue(options.readSepOptions.pErr, parser, "repsep-p-err");
        getOptionValue(options.readSepOptions.maxRandomCorrelation, parser, "repsep-max-random-correlation");
        getOptionValue(options.readSepOptions.tauMin, parser, "repsep-tau-min");
        getOptionValue(options.readSepOptions.rMin, parser, "repsep-r-min");
        getOptionValue(options.readSepOptions.minOverlap, parser, "repsep-min-overlap");
        getOptionValue(options.readSepOptions.startCompressionAt, parser, "repsep-start-compression-at");
        getOptionValue(options.readSepOptions.splitDMin, parser, "repsep-split-d-min");  // TODO(holtgrew): Remove?
    }

    // Library Options

    getOptionValue(options.fragmentSizeFactor, parser, "fragment-size-factor");
    getOptionValue(options.autoLibraryNumRecords, parser, "auto-library-num-records");
    options.autoLibraryInfo = !!options.autoLibraryNumRecords;
    getOptionValue(options.libraryInfo.median, parser, "fragment-size-median");
    getOptionValue(options.libraryInfo.stdDev, parser, "fragment-size-std-dev");
    {
        seqan::CharString tmp;
        getOptionValue(tmp, parser, "fragment-default-orientation");
        if (tmp == "R-")
            options.libraryInfo.defaultOrient = BamLibraryInfo::R_MINUS;
        else if (tmp == "R+")
            options.libraryInfo.defaultOrient = BamLibraryInfo::R_PLUS;
        else if (tmp == "F-")
            options.libraryInfo.defaultOrient = BamLibraryInfo::F_MINUS;
        else // if (tmp == "F+")
            options.libraryInfo.defaultOrient = BamLibraryInfo::F_PLUS;
    }

    // Assembly Options

    getOptionValue(options.assemblySiteWindowRadius, parser, "assembly-site-window-radius");
    getOptionValue(options.assemblySiteFringeRadius, parser, "assembly-site-fringe-radius");

    // Read Mapping Options

    getOptionValue(options.readMappingErrorRate, parser, "read-mapping-error-rate");
    getOptionValue(options.readMappingBatchSize, parser, "read-mapping-batch-size");

    // Overlapper Options

    getOptionValue(options.overlapperMinOverlapRatio, parser, "overlapper-min-overlap-ratio");
    getOptionValue(options.overlapperMaxErrorRate, parser, "overlapper-max-error-rate");

    // MSA Options

    getOptionValue(options.msaScoreMatch, parser, "msa-score-match");
    getOptionValue(options.msaScoreMismatch, parser, "msa-score-mismatch");
    getOptionValue(options.msaScoreGapOpen, parser, "msa-score-gap-open");
    getOptionValue(options.msaScoreGapExtend, parser, "msa-score-gap-extend");

    // Consensus Calling Options

    getOptionValue(options.readCorrection, parser, "no-read-correction");
    options.readCorrection = !options.readCorrection;
    
    getOptionValue(options.consensusMinBaseSupport, parser, "consensus-min-base-support");
    getOptionValue(options.consensusMinContigLengthRate, parser, "consensus-min-contig-length-rate");

    // Copy out command line.
    options.commandLine = argv[0];
    for (int i = 1; i < argc; ++i)
    {
        options.commandLine.append(" ");
        options.commandLine.append(argv[i]);
    }

    return seqan::ArgumentParser::PARSE_OK;
}

}  // anonymous namespace

// --------------------------------------------------------------------------
// Class AniseAppImpl
// --------------------------------------------------------------------------

class AniseAppImpl
{
public:
    // Entry point for the applications.
    void run(int argc, char const ** argv);

private:
    // Parse the command line.
    seqan::ArgumentParser::ParseResult parseCommandLine(int argc, char const ** argv);
    // Load or initial the global state.
    void loadGlobalState();
    // Check the input for consistency.
    void checkInput();
    // Estimate library size.
    void estimateLibrarySize();
    // Perform the preparation (step 0).
    void performPreparationStep();
    // Perform the assembly by alternating assembly and mapping.
    void performAssembly();
    // Delete unused site state files.
    void reapOldFiles();
    // Delete all remaining temporary files.
    void reapTemporaryFiles();
    // Perform the finishing that tries to maximize the usability by the end user.
    void performFinishing();
    // Check whether there are any active sites left.  If there are none then advance the state to FINISHING.
    void checkActiveSites();

    // The program configuration.
    AniseOptions options;
    // The temporary file manager.
    TemporaryFileManager tempFileManager;
    // The application state.
    AppState appState;
};

seqan::ArgumentParser::ParseResult AniseAppImpl::parseCommandLine(int argc, char const ** argv)
{
    std::cerr << "Parsing command line...\n";
    return ::parseCommandLine(options, argc, argv);
}

void AniseAppImpl::loadGlobalState()
{
    std::cerr << "Loading global state...\n";
    // Load the application state.
    appState.loadOrInit(tempFileManager);
    if (options.debugStepNo != -1)
    {
        appState.superStep = AppState::SuperStep::ASSEMBLING;
        appState.stepNo = options.debugStepNo;
        appState.save(tempFileManager);
        if (options.debugSiteID != -1)
        {
            AssemblySiteState siteState;
            siteState.load(tempFileManager, options.debugSiteID);
            siteState.active = true;
            siteState.stepNo = appState.stepNo;
            siteState.save(tempFileManager);
        }
    }
}

void AniseAppImpl::checkInput()
{
    std::cerr << "Checking input...\n";

    // Read / build FAI Index.
    seqan::FaiIndex faiIndex;
    if (read(faiIndex, toCString(options.inputReference)) != 0)
    {
        if (build(faiIndex, toCString(options.inputReference)) != 0)
            throw std::runtime_error("Could not build FAI index.");
        seqan::CharString faiPath = options.inputReference;
        append(faiPath, ".fai");
        if (write(faiIndex, toCString(faiPath)) != 0)
            throw std::runtime_error("Could not write FAI index.");
    }

    std::set<std::string> faiSeqIds;
    for (unsigned i = 0; i < numSeqs(faiIndex); ++i)
        faiSeqIds.insert(toCString(sequenceName(faiIndex, i)));

    // Open BAM Input File.
    seqan::BamStream bamStream;
    open(bamStream, toCString(options.inputMapping));
    if (!isGood(bamStream))
        throw std::runtime_error("Could not open BAM file.");

    std::set<std::string> bamSeqIds;
    for (unsigned i = 0; i < length(bamStream.header.sequenceInfos); ++i)
        bamSeqIds.insert(toCString(bamStream.header.sequenceInfos[i].i1));

    // Open VCF input file.
    seqan::VcfStream vcfStream;
    open(vcfStream, toCString(options.inputVcf));
    if (!isGood(vcfStream))
        throw std::runtime_error("Could not open input VCF file.");

    std::set<std::string> vcfSeqIds;
    for (unsigned i = 0; i < length(vcfStream.header.sequenceNames); ++i)
        vcfSeqIds.insert(toCString(vcfStream.header.sequenceNames[i]));

    if (faiSeqIds.size() != bamSeqIds.size() || vcfSeqIds.size() != bamSeqIds.size() ||
        !std::equal(faiSeqIds.begin(), faiSeqIds.end(), bamSeqIds.begin()) ||
        !std::equal(faiSeqIds.begin(), faiSeqIds.end(), vcfSeqIds.begin()))
    {
        std::cerr << "ERROR: sequences in BAM file do not match sequences in FASTA/VCF file!\n";

        std::cerr << "BAM sequence IDs:\n";
        for (std::set<std::string>::const_iterator it = bamSeqIds.begin(); it != bamSeqIds.end(); ++it)
            std::cerr << "  " << *it << "\n";

        std::cerr << "\nFASTA sequence IDs:\n";
        for (std::set<std::string>::const_iterator it = faiSeqIds.begin(); it != faiSeqIds.end(); ++it)
            std::cerr << "  " << *it << "\n";

        std::cerr << "\nVCF sequence IDs:\n";
        for (std::set<std::string>::const_iterator it = vcfSeqIds.begin(); it != vcfSeqIds.end(); ++it)
            std::cerr << "  " << *it << "\n";

        throw std::runtime_error("Mismatching sequence names!");
    }

    // Everything is fine.
}

void AniseAppImpl::estimateLibrarySize()
{
    if (appState.superStep == AppState::SuperStep::INITIAL)
    {
        if (options.autoLibraryInfo)
        {
            std::cerr << "Estimating library size...\n";

            ProgressBar pb(std::cerr, 0, options.autoLibraryNumRecords, (options.verbosity == AniseOptions::NORMAL));
            pb.setLabel("  estimation complete");
            pb.updateDisplay();
            auto updateLambda = [&pb](int done) { pb.advanceTo(done); };
            BamLibraryEstimator estimator(options.autoLibraryNumRecords);
            if (estimator.run(options.libraryInfo, toCString(options.inputMapping), updateLambda) != 0)
                throw std::runtime_error("Estimating library size failed!");
            pb.finish();  // finish in case there were too few records
        }
        appState.libraryInfos[0] = options.libraryInfo;
        appState.save(tempFileManager);
    }
    else
    {
        std::cerr << "Getting library size estimate from app state.\n";
        options.libraryInfo = appState.libraryInfos[0];
    }

    if (options.verbosity >= AniseOptions::NORMAL)
        std::cerr << "  library size: median=" << options.libraryInfo.median << ", stdDev="
                  << options.libraryInfo.stdDev << "\n";

    if (false && options.autoTuning)  // TODO(holtgrew): Remove.
    {
        int const MIN_OVERLAP = 14;
        int const MIN_OVL_ERR = 2;
        int const MIN_MAP_ERR = 2;
        if (options.verbosity >= AniseOptions::NORMAL)
            std::cerr << "Auto-tuning...\n";

        if (options.libraryInfo.avgReadLen / 100.0 * options.overlapperMinOverlapRatio < MIN_OVERLAP)
        {
            options.overlapperMinOverlapRatio = ceil(MIN_OVERLAP * 100.0 / options.libraryInfo.avgReadLen);
            if (options.verbosity >= AniseOptions::NORMAL)
                std::cerr << "  min overlap ratio: " << options.overlapperMinOverlapRatio << "\n";
        }
        if (options.libraryInfo.avgReadLen / 100.0 * options.overlapperMaxErrorRate < MIN_OVL_ERR)
        {
            options.overlapperMaxErrorRate = ceil(MIN_OVL_ERR * 100.0 / options.libraryInfo.avgReadLen);
            if (options.verbosity >= AniseOptions::NORMAL)
                std::cerr << "  overlapper max error rate: " << options.overlapperMaxErrorRate << "\n";
        }
        if (options.libraryInfo.avgReadLen / 100.0 * options.readMappingErrorRate < MIN_MAP_ERR)
        {
            options.readMappingErrorRate = ceil(MIN_MAP_ERR * 100.0 / options.libraryInfo.avgReadLen);
            if (options.verbosity >= AniseOptions::NORMAL)
                std::cerr << "  read mapping max error rate: " << options.readMappingErrorRate << "\n";
        }

        if (options.verbosity >= AniseOptions::NORMAL)
            std::cerr << "  DONE\n";
    }
}

void AniseAppImpl::performPreparationStep()
{
    if (!appState.isInitial())
        return;

    std::cerr << "Perform preparation step...\n";
    // Perform the preparation step.
    PreparationStep prepStep(tempFileManager, options, appState);
    prepStep.run();

    // Application step was updated on disk, update from there.
    appState.load(tempFileManager);
    SEQAN_CHECK(appState.isAssembling(), "Incorrect state.");
}

void AniseAppImpl::reapOldFiles()
{
    if (!options.cleanUpTemporaryFiles)
        return;
    for (int siteID = 0; siteID < appState.numSites; ++siteID)
    {
        AssemblySiteState siteState;
        siteState.load(tempFileManager, siteID);

        // Delete previous state if contains data for the current step.
        if (siteState.stepNo == appState.stepNo)
            tempFileManager.reapSiteFiles(siteState.stepNo - 1, siteID);
    }
}

void AniseAppImpl::reapTemporaryFiles()
{
    if (!options.cleanUpTemporaryFiles)
        return;
    tempFileManager.cleanup();
}

void AniseAppImpl::performAssembly()
{
    if (!appState.isAssembling())
        return;  // already done assembling
    std::cerr << "Performing assembly...\n";
    SEQAN_CHECK(appState.stepNo >= 1, "Invalid step no, should not be 0 (initial).");

    for (int i = appState.stepNo; i <= options.recursionMaxSteps && appState.isAssembling(); ++i)
    {
        std::cerr << "Performing step #" << i << "\n";
        // Perform the assembly substep for the current step.
        AssemblySubstep assembly(tempFileManager, options);
        assembly.run();

        // Check whether there are any active sites left.  If this is not then case then we switch to finishing.
        checkActiveSites();

        // Perform read mapping substep for the current step if the recursion is not forced to stop here by a configured
        // limit.
        if (appState.isAssembling() && i < options.recursionMaxSteps)
        {
            ReadMappingSubstep readMapping(tempFileManager, i, options);
            readMapping.run();

            // Update the application state from the disk again.  The read mapping can switch the step to FINISHING if
            // no or too few new read alignments were found.
            appState.load(tempFileManager);
        }
        // Remove the files of the previous step for all sites that were still active in this step.
        reapOldFiles();

        // Advance step in app state.  Change application state to FINISHING if (i == options.recursionMaxSteps).
        appState.load(tempFileManager);
        if (appState.isAssembling() && i == options.recursionMaxSteps)
            appState.superStep = AppState::SuperStep::FINISHING;
        else
            appState.stepNo = (i + 1);
        appState.save(tempFileManager);
    }

    // Application step was updated on disk, update from there.
    appState.load(tempFileManager);
    SEQAN_CHECK(appState.isFinishing(), "State should be finishing after assembly.");
}

void AniseAppImpl::checkActiveSites()
{
    if (!appState.isAssembling())
        return;
    for (int siteID = 0; siteID < appState.numSites; ++siteID)
    {
        AssemblySiteState siteState;
        siteState.load(tempFileManager, siteID);
        if (siteState.active)
            return;  // one is enough, no need to advance to finishing
    }

    // If we reach here, we should start finishing.
    appState.superStep = AppState::SuperStep::FINISHING;
    appState.save(tempFileManager);
}

void AniseAppImpl::performFinishing()
{
    if (!appState.isFinishing())
        return;  // already done finishing
    std::cerr << "Performing finishing...\n";

    FinishingStep finishing(tempFileManager, options);
    finishing.run();
}

void AniseAppImpl::run(int argc, char const ** argv)
{
    // Parse command line and handle error vs. special action (such as --help and --write-ctd).
    auto res = parseCommandLine(argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
    {
        if (res == seqan::ArgumentParser::PARSE_ERROR)
            throw std::runtime_error("Could not parse command line.");
        else
            return;
    }

    // Initialize the path manager object, creates temporary directory.  This enables us to load the global application
    // state.
    tempFileManager.init(toCString(options.tmpDir));
    loadGlobalState();

    // Initialize the TimeLog singleton.
    TimeLog::instance().open(tempFileManager);

    // Check the input and estimate library size.
    checkInput();
    estimateLibrarySize();

    // Perform preparation step, assembly, and finishing.
    withTimeLog("PREPARATION", [&]() {
            performPreparationStep();
        });
    withTimeLog("ASSEMBLY", [&]() {
            performAssembly();
        });
    withTimeLog("FINISHING", [&]() {
            performFinishing();
        });

    reapTemporaryFiles();
}

// --------------------------------------------------------------------------
// Class AniseApp
// --------------------------------------------------------------------------

AniseApp::AniseApp() : impl(new AniseAppImpl())
{}

AniseApp::~AniseApp()
{}

int AniseApp::run(int argc, char const ** argv)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    int result = 0;
    try
    {
        impl->run(argc, argv);
    }
    catch (std::exception const & e)
    {
        std::cerr << "\nERROR: " << e.what() << "\n";
        result = 1;
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cerr << "\nTook " << elapsed_seconds.count() << " s\n";
    return result;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    AniseApp app;
    return app.run(argc, argv);
}
