// ==========================================================================
//                                 SeqCons
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "seqcons_options.h"

namespace {  // anonymous namespace

// --------------------------------------------------------------------------
// Function strToMethod()
// --------------------------------------------------------------------------

SeqConsOptions::Operation strToMethod(std::string opStr)
{
    if (opStr == "align_consensus")
        return SeqConsOptions::ALN_CONSENSUS;
    else if (opStr == "overlap_consensus")
        return SeqConsOptions::OVL_CONSENSUS;
    else if (opStr == "contig_consensus")
        return SeqConsOptions::CTG_CONSENSUS;
    else if (opStr == "pos_consensus")
        return SeqConsOptions::POS_CONSENSUS;
    else if (opStr == "realign")
        return SeqConsOptions::REALIGN;
    else  // opStr == "nop"
        return SeqConsOptions::NOP;
}

// --------------------------------------------------------------------------
// Function operationStr()
// --------------------------------------------------------------------------

char const * operationStr(SeqConsOptions::Operation op)
{
    static char const * labels[] = {
        "convert only",
        "realign only",
        "consensus with positions",
        "contig-wise MSA + consensus",
        "overlap MSA + consensus",
        "global MSA + consensus"
    };
    return labels[op];
}

// --------------------------------------------------------------------------
// Function verbosityStr()
// --------------------------------------------------------------------------

char const * verbosityStr(int verbosity)
{
    static char const * labels[] = {
        "QUIET",
        "NORMAL",
        "VERBOSE",
        "VERY VERBOSE"
    };
    verbosity = std::max(std::min(verbosity, 3), 0);
    return labels[verbosity];
}

}  // anonymous namespace

// --------------------------------------------------------------------------
// Class SeqConsOptions
// --------------------------------------------------------------------------

void SeqConsOptions::checkConsistency()
{
    // Position information is only available when input is in SAM format and is required for
    // position-/contig-based consensus and realignment.
    seqan::CharString inFileLowerCase = inputFile;
    seqan::toLower(inFileLowerCase);
    if ((operation == POS_CONSENSUS || operation == CTG_CONSENSUS || operation == REALIGN) &&
        (!endsWith(inFileLowerCase, ".sam")))
        throw std::runtime_error("SAM input required for coordinates.  Either specify SAM file for the "
                                 "input or use \"--method overlap_consensus\" or \"--method align_consensus\".");
}

void SeqConsOptions::print(std::ostream & out) const
{
    out << "VERBOSITY              \t" << verbosityStr(verbosity) << "\n"
        << "OPERATION              \t" << operationStr(operation) << "\n"
        << "\n"
        << "INPUT FILE             \t" << inputFile << "\n"
        << "OUTPUT READS           \t" << outputFileAlignment << "\n"
        << "OUTPUT CONSENSUS       \t" << outputFileConsensus << "\n"
        << "\n"
        << "OVERLAP MIN LENGTH     \t" << overlapMinLength << "\n"
        << "OVERLAP MAX ERR-RATE   \t" << overlapMaxErrorRate << "\n"
        << "OVERLAP MIN COUNT      \t" << overlapMinCount << "\n"
        << "OVERLAP WINDOW SIZE    \t" << overlapWindowSize << "\n"
        << "\n"
        << "K-MER SIZE             \t" << kMerSize << "\n"
        << "K-MER MAX OCCURRENCES  \t" << kMerMaxOcc << "\n"
        << "\n"
        << "REALIGNMENT BANDWIDTH  \t" << reAlignmentBandwidth << "\n"
        << "REALIGNMENT ENVIRONMENT\t" << reAlignmentEnvironment << "\n";
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(SeqConsOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("seqcons2");
    // Set short description, version, and date.
    setShortDescription(parser, "Compute consensus from sequences.");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Define usage line and long description.
    addUsageLine(parser,
                 "\\fB-i\\fP \\fIINPUT.{fa,sam}\\fP [\\fB-oa\\fP \\fIOUT_ALIGN.{fa,sam}\\fP] "
                 "[\\fB-oc\\fP \\fIOUT_CONSENSUS.fa\\fP]");
    addDescription(parser,
                   "Compute consensus from sequences with and without approximate alignment information.");

    // Overall Program Options
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addOption(parser, seqan::ArgParseOption("m", "method", "Method to perform.  See section \\fIMethods\\fP "
                                            "below for details.", seqan::ArgParseOption::STRING, "METHOD"));
    setValidValues(parser, "method", "nop realign align_consensus overlap_consensus contig_consensus pos_consensus");
    setDefaultValue(parser, "method", "pos_consensus");

    // I/O Options
    addSection(parser, "I/O Options");

    addOption(parser, seqan::ArgParseOption("i", "input-file", "Input file.", seqan::ArgParseOption::INPUT_FILE,
                                            "INPUT"));
    setRequired(parser, "input-file", true);
    setValidValues(parser, "input-file", "sam fa fasta");

    addOption(parser, seqan::ArgParseOption("oa", "output-alignment-file", "Output file with alignment.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "OUT_ALIGNMENT"));
    setRequired(parser, "output-alignment-file", false);
    setValidValues(parser, "output-alignment-file", "sam txt");

    addOption(parser, seqan::ArgParseOption("oc", "output-consensus-file", "Output file with consensus sequence.",
                                            seqan::ArgParseOption::OUTPUT_FILE, "OUT_CONSENSUS"));
    setRequired(parser, "output-consensus-file", false);
    setValidValues(parser, "output-consensus-file", "fa fasta");

    // Alignment Quality Filter Options
    addSection(parser, "Alignment Quality Filter Options");

    addOption(parser, seqan::ArgParseOption("", "overlap-min-length", "Minimal overlap length.",
                                            seqan::ArgParseOption::INTEGER, "LENGTH"));
    setMinValue(parser, "overlap-min-length", "0");
    setDefaultValue(parser, "overlap-min-length", "20");

    addOption(parser, seqan::ArgParseOption("", "overlap-max-error", "Maximal error rate in overlap as percentage.",
                                            seqan::ArgParseOption::DOUBLE, "RATE"));
    setMinValue(parser, "overlap-max-error", "0.0");
    setDefaultValue(parser, "overlap-max-error", "5.0");

    addOption(parser, seqan::ArgParseOption("", "overlap-min-count", "Minimal overlap count.",
                                            seqan::ArgParseOption::INTEGER, "COUNT"));
    setMinValue(parser, "overlap-min-count", "0");
    setDefaultValue(parser, "overlap-min-count", "3");

    addOption(parser, seqan::ArgParseOption("", "overlap-window-size", "Window size to look for alignments.",
                                            seqan::ArgParseOption::INTEGER, "SIZE"));
    setMinValue(parser, "overlap-window-size", "0");
    setDefaultValue(parser, "overlap-window-size", "20");

    // K-mer Filter Options
    addSection(parser, "K-Mer Filter Options");

    addOption(parser, seqan::ArgParseOption("", "k-mer-size", "The k-mer size to use.",
                                            seqan::ArgParseOption::INTEGER, "LENGTH"));
    setMinValue(parser, "k-mer-size", "5");
    setDefaultValue(parser, "k-mer-size", "20");

    addOption(parser, seqan::ArgParseOption("", "k-mer-max-occ", "Ignore k-mer with higher occurrence count, 0 to disable.",
                                            seqan::ArgParseOption::INTEGER, "COUNT"));
    setMinValue(parser, "k-mer-max-occ", "0");
    setDefaultValue(parser, "k-mer-max-occ", "200");

    // Realignment Options
    addSection(parser, "Realignment Options");

    addOption(parser, seqan::ArgParseOption("", "realign-bandwidth",
                                            "Bandwidth to use for pairwise alignments in realignment.",
                                            seqan::ArgParseOption::INTEGER, "LENGTH"));
    setMinValue(parser, "realign-bandwidth", "5");
    setDefaultValue(parser, "realign-bandwidth", "10");

    addOption(parser, seqan::ArgParseOption("", "realign-environment",
                                            "Environment for extraction in realignment.",
                                            seqan::ArgParseOption::INTEGER, "COUNT"));
    setMinValue(parser, "realign-environment", "5");
    setDefaultValue(parser, "realign-environment", "20");

    // Add Methods Section
    addTextSection(parser, "Methods");
    addListItem(parser, "\\fBnop\\fP",
                "Perform no action, just perform file conversion if possible.");
    addListItem(parser, "\\fBrealign\\fP",
                "Perform realignment, requires input to be a SAM file to provide approximate position "
                "information, creates consensus sequence after realignment.");
    addListItem(parser, "\\fBoverlap_consensus\\fP",
                "Perform MSA with overlap alignments of the input ignoring any given coordinates, then realign. "
                "This is most suited when computing the consensus of reads where the underlying sequence is very "
                "similar and most differences stem from sequencing errors and not genomic variation. All "
                "pairwise alignments computed here are banded.");
    addListItem(parser, "\\fBalign_consensus\\fP",
                "Perform MSA with global alignments of the input ignoring any given coordinates, then realign. "
                "This will computed unbanded global ends-gap free pairwise alignments.  This is also suitable "
                "when aligning different sequences, e.g. clustered transcripts.  Using this method, seqcons "
                "will be similar to calling SeqAn::T-Coffee, followed by realignment and consensus computation.");
    addListItem(parser, "\\fBcontig_consensus\\fP",
                "Perform MSA of the input, contig by contig, requires contig information, then realign. Input "
                "must be SAM.");
    addListItem(parser, "\\fBpos_consensus\\fP",
                "Perform consensus of the input, then realign. Requires approximate coordinate information in "
                "SAM file.");

    // Add Output Section
    addTextSection(parser, "Output Formats");
    addText(parser,
            "The program can write out the consensus sequence in FASTA format and optionally the alignment of the "
            "input sequences against the consensus in SAM/BAM format.  When using the extension \\fI.txt\\fP, seqcons "
            "will write out the MSA as a plain text visualization.");

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBseqcons\\fP \\fB-m\\fP \\fIovl_consensus\\fP \\fB-i\\fP \\fIreads.fa\\fP \\fB-oa\\fP "
                "\\fIout.sam\\fP \\fB-oc\\fP \\fIcons.fa\\fP",
                "Compute MSA of the sequences in \\fIreads.fa\\fP.  The consensus sequence is written to "
                "\\fIcons.fa\\fP and the alignment of the sequences in \\fIreads.fa\\fP is written to "
                "\\fIout.sam\\fP.");
    addListItem(parser,
                "\\fBseqcons\\fP \\fB-m\\fP \\fIrealign\\fP \\fB-i\\fP \\fIin.sam\\fP \\fB-oa\\fP \\fIout.sam\\fP",
                "Read in multi-read alignment from \\fIin.sam\\fP, refine it using Anson-Myers realignment and "
                "write out the refined alignment to \\fIout.sam\\fP");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.

    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    std::string opStr;
    getOptionValue(opStr, parser, "method");
    options.operation = strToMethod(opStr.c_str());

    getOptionValue(options.inputFile, parser, "input-file");
    getOptionValue(options.outputFileAlignment, parser, "output-alignment-file");
    getOptionValue(options.outputFileConsensus, parser, "output-consensus-file");

    getOptionValue(options.overlapMinLength, parser, "overlap-min-length");
    getOptionValue(options.overlapMaxErrorRate, parser, "overlap-max-error");
    getOptionValue(options.overlapWindowSize, parser, "overlap-window-size");

    getOptionValue(options.kMerSize, parser, "k-mer-size");
    getOptionValue(options.kMerMaxOcc, parser, "k-mer-max-occ");

    getOptionValue(options.reAlignmentBandwidth, parser, "realign-bandwidth");
    getOptionValue(options.reAlignmentEnvironment, parser, "realign-environment");

    return seqan::ArgumentParser::PARSE_OK;
}
