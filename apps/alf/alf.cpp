// ==========================================================================
//                  ALF - Alignment free sequence comparison
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
// Author: Jonathan Goeke <goeke@molgen.mpg.de>
// ==========================================================================
// Alignment free sequence comparison.
//
// This application can be used to calculate pairwise scores of DNA Sequences
// without alignments.
//
// The following scores are implemented: N2, D2, D2Star, D2z.
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/misc/edit_environment.h>
#include <seqan/arg_parse.h>
#include <seqan/alignment_free.h>

using namespace seqan;
using namespace std;

// TODO(holtgrew): Adapt parameters to naming conventions, i.e. use --parameter-name.

int main(int argc, const char * argv[])
{
    // -----------------------------------------------------------------------
    // Setup argument parser
    // -----------------------------------------------------------------------
    seqan::ArgumentParser parser("alf");

    // Set short description, version, date.
    setShortDescription(parser, "Alignment free sequence comparison");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setCategory(parser, "Sequence Comparison");

    // Usage line and description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-i\\fP \\fIIN.FASTA\\fP [\\fB-o\\fP \\fIOUT.TXT\\fP]");
    addDescription(parser, "Compute pairwise similarity of sequences using alignment-free methods in \\fIIN.FASTA\\fP and write out tab-delimited matrix with pairwise scores to \\fIOUT.TXT\\fP.");
    addOption(parser, seqan::ArgParseOption("v", "verbose", "When given, details about the progress are printed to the screen."));

    // Options Section: Input / Output parameters.
    addSection(parser, "Input / Output");
    addOption(parser, seqan::ArgParseOption("i", "input-file", "Name of the multi-FASTA input file.", seqan::ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "input-file", seqan::SeqFileIn::getFileExtensions());
    setRequired(parser, "input-file");
    addOption(parser, seqan::ArgParseOption("o", "output-file", "Name of the file to which the tab-delimtied matrix with pairwise scores will be written to.  Default is to write to stdout.", seqan::ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "output-file", "alf.tsv");

    addSection(parser, "General Algorithm Parameters");
    addOption(parser, seqan::ArgParseOption("m", "method", "Select method to use.", seqan::ArgParseArgument::STRING, "METHOD"));
    setValidValues(parser, "method", "N2 D2 D2Star D2z");
    setDefaultValue(parser, "method", "N2");
    addOption(parser, seqan::ArgParseOption("k", "k-mer-size", "Size of the k-mers.", seqan::ArgParseArgument::INTEGER, "K"));
    setDefaultValue(parser, "k-mer-size", "4");
    addOption(parser, seqan::ArgParseOption("mo", "bg-model-order", "Order of background Markov Model.", seqan::ArgParseArgument::INTEGER, "ORDER"));
    setDefaultValue(parser, "bg-model-order", "1");

    addSection(parser, "N2 Algorithm Parameters");
    // addText(parser, "The following parameters are only used in the N2 algorithm.");
    addOption(parser, seqan::ArgParseOption("rc", "reverse-complement", "Which strand to score.  Use \\fIboth_strands\\fP to score both strands simultaneously.", seqan::ArgParseArgument::STRING, "MODE"));
    setValidValues(parser, "reverse-complement", "input both_strands mean min max");
    setDefaultValue(parser, "reverse-complement", "input");
    addOption(parser, seqan::ArgParseOption("mm", "mismatches", "Number of mismatches, one of \\fI0\\fP and \\fI1\\fP.  When \\fI1\\fP is used, N2 uses the k-mer-neighbour with one mismatch.", seqan::ArgParseArgument::INTEGER, "MISMATCHES"));
    setDefaultValue(parser, "mismatches", "0");
    addOption(parser, seqan::ArgParseOption("mmw", "mismatch-weight", "Real-valued weight of counts for words with mismatches.", seqan::ArgParseArgument::DOUBLE, "WEIGHT"));
    setDefaultValue(parser, "mismatch-weight", "0.1");
    addOption(parser, seqan::ArgParseOption("kwf", "k-mer-weights-file", "Print k-mer weights for every sequence to this file if given.", seqan::ArgParseArgument::OUTPUT_FILE, "FILE.TXT"));
    setValidValues(parser, "k-mer-weights-file", "txt");

    addTextSection(parser, "Contact and References");
    addListItem(parser, "For questions or comments, contact:", "Jonathan Goeke <goeke@molgen.mpg.de>");
    addListItem(parser, "Please reference the following publication if you used ALF or the N2 method for your analysis:", "Jonathan Goeke, Marcel H. Schulz, Julia Lasserre, and Martin Vingron. Estimation of Pairwise Sequence Similarity of Mammalian Enhancers with Word Neighbourhood Counts. Bioinformatics (2012).");
    addListItem(parser, "Project Homepage:", "http://www.seqan.de/projects/alf");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return (res == seqan::ArgumentParser::PARSE_ERROR);

    // Declare all parameters
    String<char> kmerWeightsFileTmp;
    String<char> inFileTmp;
    String<char> outFileTmp;

    getOptionValue(inFileTmp, parser, "input-file");
    String<char, CStyle> inFile = inFileTmp;

    getOptionValue(outFileTmp, parser, "output-file");
    String<char, CStyle> outFile = outFileTmp;

    String<char> method;
    getOptionValue(method, parser, "method");

    int kmerSize = 0;
    getOptionValue(kmerSize, parser, "k-mer-size");

    int bgModelOrder = 1;
    getOptionValue(bgModelOrder, parser, "bg-model-order");

    String<char> revComp;
    getOptionValue(revComp, parser, "reverse-complement");
    if (revComp == "input")
        clear(revComp);

    unsigned mismatches = 0;
    getOptionValue(mismatches, parser, "mm");

    double mismatchWeight = 0.1;
    getOptionValue(mismatchWeight, parser, "mmw");

    String<char, CStyle> kmerWeightsFile;
    getOptionValue(kmerWeightsFileTmp, parser, "k-mer-weights-file");
    kmerWeightsFile = kmerWeightsFileTmp;

    bool verbose = isSet(parser, "verbose");

    // Definition of type DNA string sets
    typedef String<Dna5>        TText;
    typedef StringSet<TText>    TStringSet;

    // Definition of mxn two-dimensional matrix
    typedef Matrix<double, 2> TMatrix;

    TMatrix myMatrix;  // myMatrix stores pairwise kmerScores

    TStringSet mySequenceSet;  // mySequenceSet stores all sequences from the multi-fasta file
    if (inFile != "")  // read in file
    {
        SeqFileIn seqFile;
        open(seqFile, toCString(inFile));
        StringSet<CharString> seqIDs;
        readRecords(seqIDs, mySequenceSet, seqFile);
    }

    // Dispatch to alignment free comparisons with different scores.
    if (method == "D2")
    {
        AFScore<D2> myScoreD2(kmerSize, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreD2);
    }
    else if (method == "D2z")
    {
        AFScore<D2z> myScoreD2z(kmerSize, bgModelOrder, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreD2z);
    }
    else if (method == "D2Star")
    {
        AFScore<D2Star> myScoreD2Star(kmerSize, bgModelOrder, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreD2Star);
    }
    else if (method == "N2")
    {
        AFScore<N2> myScoreN2(kmerSize, bgModelOrder, revComp, mismatches, mismatchWeight, kmerWeightsFile, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreN2);
    }

    // Write out resulting matrix; to file if file name was given, to stdout otherwise.
    if (!empty(outFile))
    {
        ofstream myfile(outFile, std::ios::binary | std::ios::out);
        myfile << myMatrix;
        myfile.close();
    }
    else
    {
        std::cout << "\n" << myMatrix;
    }
    return 0;
}
