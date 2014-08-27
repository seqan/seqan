// ==========================================================================
//                           LINEAR BLOCK ALIGNER
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

#include <algorithm>
#include <iostream>
#include <string>

#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

// ==========================================================================
// Class LbaOptions
// ==========================================================================

struct LbaOptions
{
    int verbosity { 1 };

    // -----------------------------------------------------------------------
    // Input And Output Files
    // -----------------------------------------------------------------------

    std::string inSeq0, inSeq1;
    std::string out;

    // -----------------------------------------------------------------------
    // Alignment Algorithms
    // -----------------------------------------------------------------------

    // K-mer Size.
    int kMerSize       { 20 };  // ~5% errors

    // Local Alignment Scores.
    int scoreMatch     { 3 };
    int scoreMismatch  { -5 };
    int scoreGapOpen   { -8 };
    int scoreGapExtend { -2 };

    // Gaps in seqH with greater length are considered "uncovered" instead of errors.
    int longGapMinLen  { 5 };

    // Alignment bandwidth.
    int band           { 15 };

    // Smallest score for a block.
    int minScore       { 20 };

    // X-dropoff
    int xDrop          { 30 };

    // Ignore overabundant q-grams.
    unsigned maxOcc    { 200 };

    void print(std::ostream & out) const;
};

void LbaOptions::print(std::ostream & out) const
{
    out << "__OPTIONS_____________________________________________________________________\n"
        << "\n"
        << "INPUT 0\t" << inSeq0 << "\n"
        << "INPUT 1\t" << inSeq1 << "\n"
        << "OUTPUT \t" << this->out << "\n"
        << "\n"
        << "SCORE\n"
        << "\tMATCH    \t" << scoreMatch << "\n"
        << "\tMISMATCH \t" << scoreMismatch << "\n"
        << "\tGAP OPEN \t" << scoreGapOpen << "\n"
        << "\tGAP EXT  \t" << scoreGapExtend << "\n"
        << "\n"
        << "BAND       \t" << band << "\n"
        << "MIN SCORE  \t" << minScore << "\n"
        << "X-DROP     \t" << xDrop << "\n"
        << "\n";
}

// ==========================================================================
// Function parseCommandLine()
// ==========================================================================

seqan::ArgumentParser::ParseResult
parseCommandLine(LbaOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("lba");
    // Set short description, version, and date.
    setShortDescription(parser, "Local Block Aligner");
#ifdef SEQAN_REVISION
    setVersion(parser, "0.1 [" + std::string(SEQAN_REVISION) + "]");
#else
    setVersion(parser, "0.1");
#endif
#ifdef SEQAN_DATE
    setDate(parser, SEQAN_DATE);
#else
    setDate(parser, "May 2014");
#endif
    setCategory(parser, "Pairwise Alignment");

    // Define usage line and long description.
    addUsageLine(parser, "[OPTIONS] -o OUT.txt -i0 IN0.fa -i1 IN1.fa");

    addDescription(parser, "Perform local pairwise alignment based on seeds.");

    // ----------------------------------------------------------------------
    // General Options
    // ----------------------------------------------------------------------

    addSection(parser, "General Options");

    // We require one argument.
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // ----------------------------------------------------------------------
    // Input / Output Options
    // ----------------------------------------------------------------------

    addSection(parser, "Input / Output");

    addOption(parser, seqan::ArgParseOption("i0", "in-first", "First FASTA file with reference.",
                                            seqan::ArgParseOption::INPUTFILE, "FASTA"));
    setValidValues(parser, "in-first", "fasta fa");
    setRequired(parser, "in-first");

    addOption(parser, seqan::ArgParseOption("i1", "in-second", "Second FASTA file with reference.",
                                            seqan::ArgParseOption::INPUTFILE, "FASTA"));
    setValidValues(parser, "in-second", "fasta fa");
    setRequired(parser, "in-second");

    addOption(parser, seqan::ArgParseOption("o", "out", "Alignment out file.",
                                            seqan::ArgParseOption::INPUTFILE, "TXT"));
    setValidValues(parser, "out", "txt");

    // ----------------------------------------------------------------------
    // Alignment Options
    // ----------------------------------------------------------------------

    addSection(parser, "Alignment");

    addOption(parser, seqan::ArgParseOption("", "k-mer-size", "K-mer size.",
                                            seqan::ArgParseOption::INTEGER, "K"));
    setMinValue(parser, "k-mer-size", "5");
    setDefaultValue(parser, "k-mer-size", "20");

    addOption(parser, seqan::ArgParseOption("", "score-match", "Match score.",
                                            seqan::ArgParseOption::INTEGER, "SCORE"));
    setMinValue(parser, "score-match", "1");
    setDefaultValue(parser, "score-match", "10");

    addOption(parser, seqan::ArgParseOption("", "score-mismatch", "Mismatch score.",
                                            seqan::ArgParseOption::INTEGER, "SCORE"));
    setMaxValue(parser, "score-mismatch", "0");
    setDefaultValue(parser, "score-mismatch", "-8");

    addOption(parser, seqan::ArgParseOption("", "score-gap-open", "Gap open score.",
                                            seqan::ArgParseOption::INTEGER, "SCORE"));
    setMaxValue(parser, "score-gap-open", "0");
    setDefaultValue(parser, "score-gap-open", "-20");

    addOption(parser, seqan::ArgParseOption("", "score-gap-extend", "Gap extension score.",
                                            seqan::ArgParseOption::INTEGER, "SCORE"));
    setMaxValue(parser, "score-gap-extend", "0");
    setDefaultValue(parser, "score-gap-extend", "-1");

    addOption(parser, seqan::ArgParseOption("", "min-block-score", "Minimal local score.",
                                            seqan::ArgParseOption::INTEGER, "SCORE"));
    setMinValue(parser, "min-block-score", "0");
    setDefaultValue(parser, "min-block-score", "120");

    addOption(parser, seqan::ArgParseOption("", "x-drop", "Maximal x-drop.",
                                            seqan::ArgParseOption::INTEGER, "SCORE"));
    setMinValue(parser, "x-drop", "0");
    setDefaultValue(parser, "x-drop", "20");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    options.verbosity = 1;
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.inSeq0, parser, "in-first");
    getOptionValue(options.inSeq1, parser, "in-second");
    seqan::CharString out;
    getOptionValue(out, parser, "out");
    options.out = toCString(out);

    getOptionValue(options.kMerSize, parser, "k-mer-size");
    getOptionValue(options.scoreMatch, parser, "score-match");
    getOptionValue(options.scoreMismatch, parser, "score-mismatch");
    getOptionValue(options.scoreGapOpen, parser, "score-gap-open");
    getOptionValue(options.scoreGapExtend, parser, "score-gap-extend");

    getOptionValue(options.minScore, parser, "min-block-score");
    getOptionValue(options.xDrop, parser, "x-drop");

    return seqan::ArgumentParser::PARSE_OK;
}

// ===========================================================================
// Function trimAfterSpace()
// ===========================================================================

// Trim after the first whitespace.
void trimAfterSpace(seqan::CharString & s)
{
    unsigned i = 0;
    for (; i < length(s); ++i)
        if (isspace(s[i]))
            break;
    resize(s, i);
}

// ===========================================================================
// Class Counter
// ===========================================================================

// Counters for horizontal and vertical direction.

struct Counter
{
    unsigned h { 0 };
    unsigned v { 0 };

    unsigned max() const { return std::max(h, v); }
    unsigned sum() const { return h + v; }
};

// ===========================================================================
// Function countLeading()
// ===========================================================================

Counter countLeading(seqan::Align<seqan::Dna5String const> const & align)
{
    Counter result;

    auto itH = begin(row(align, 0)), itV = begin(row(align, 1));
    while ((isGap(itH) || isGap(itV)) && !atEnd(itH) && !atEnd(itV))
    {
        unsigned h = 0, v = 0;
        if (isGap(itH))
            h = countGaps(itH);
        else if (isGap(itV))
            v = countGaps(itV);
        result.h += h;
        result.v += v;
        itH += std::max(h, v);
        itV += std::max(h, v);
    }

    return result;
}

// ===========================================================================
// Function countTrailing()
// ===========================================================================

Counter countTrailing(seqan::Align<seqan::Dna5String const> const & align)
{
    Counter result;
    auto itH = end(row(align, 0)), itV = end(row(align, 1));

    enum { NONE, GAP_H, GAP_V } state = NONE;
    do {
        if (state == NONE && !isGap(itV) && !isGap(itH))
            break;  // no trailing gap
        --itH, --itV;
        if (isGap(itH))
        {
            if (state == GAP_V)
                break;
            state = GAP_H;
            ++result.h;
        }
        if (state == GAP_H && !isGap(itH))
            break;
        if (isGap(itV))
        {
            if (state == GAP_H)
                break;
            state = GAP_V;
            ++result.v;
        }
        if (state == GAP_V && !isGap(itV))
            break;
    } while (itH != begin(row(align, 0)) && itV != begin(row(align, 1)));

    return result;
}

// ===========================================================================
// Function computeChunks()
// ===========================================================================

// Return list of intervals in seqH/0 that are covered.

std::vector<std::pair<unsigned, unsigned>>
computeChunks(seqan::Align<seqan::Dna5String const> const & align,
              int longGapMinLen)
{
    std::vector<std::pair<unsigned, unsigned> > result;
    SEQAN_CHECK(!isGap(begin(row(align, 0))), "Must not start with gap.");
    SEQAN_CHECK(!isGap(begin(row(align, 1))), "Must not start with gap.");

    std::pair<unsigned, unsigned> current(0, 0);

    for (auto itH = begin(row(align, 0)), itV = begin(row(align, 1)); !atEnd(itH); /*see below*/)
    {
        current.second = position(itV) + 1;
        if (isGap(itV))
        {
            if ((int)countGaps(itV) >= longGapMinLen)
            {
                result.push_back(current);
                int num = countGaps(itV);
                itH += num;
                itV += num;
                current.first = current.second = position(itV);
            }
            else
            {
                int num = countGaps(itV);
                itH += num;
                itV += num;
            }
        }
        else
        {
            ++itH;
            ++itV;
        }
    }
    if (current.first != current.second)
        result.push_back(current);

    // Remove Ns and gaps from borders of chunks.
    auto & rowV = row(align, 1);
    for (auto & chunk : result)
    {
        while (chunk.first <= chunk.second)
        {
            if (isGap(rowV, chunk.first))
                ++chunk.first;
            else if (rowV[chunk.first] == 'N')
                ++chunk.first;
            else
                break;
        }
        while (chunk.first <= chunk.second)
        {
            if (isGap(rowV, chunk.second - 1))
                --chunk.second;
            else if (rowV[chunk.second - 1] == 'N')
                --chunk.second;
            else
                break;
        }
    }

    // Remove empty chunks.
    result.erase(std::remove_if(result.begin(), result.end(),
                                [](std::pair<unsigned, unsigned> const & p) { return (p.first == p.second); }),
                 result.end());

    // std::cerr << "CHUNKS\n";
    // for (auto pair : result)
    //     std::cerr << pair.first << " -- " << pair.second << "\n";

    return result;
}

// ==========================================================================
// Class LocalBlockAlignerResult
// ==========================================================================

struct LocalBlockAlignerResult
{
    struct Stats
    {
        int alignmentLength { 0 };  // length of alignment
        int numGaps         { 0 };  // number of gaps in blocks
        int numMatches      { 0 };  // number of matches in blocks
        int numMismatches   { 0 };  // number of mismatches in blocks

        int missedBlockSeqH { 0 };  // unmatched in block from seqH
        int missedBlockSeqV { 0 };  // unmatched in block from seqV
        int inBlockSeqH     { 0 };  // number of seqH chars in blocks
        int inBlockSeqV     { 0 };  // number of seqV chars in blocks
    };

    // The resulting alignments from the local block aligner.
    seqan::Align<seqan::Dna5String const> align;
    // Alignment statistics.
    Stats stats;

    void print(std::ostream & out) const;
};

void LocalBlockAlignerResult::print(std::ostream & out) const
{
    int seqHLen = 0, seqVLen = 0;
    seqHLen = length(source(row(align, 0)));
    seqVLen = length(source(row(align, 1)));

    out << "#ali_len\tnum_gap\tnum_match\tnum_mismatch\tcov_seq0\tcov_seq1\tuncov_num_seq0\tuncov_num_seq1\t"
        << "cov_perc_seq0\tcov_perc_seq1\terr_rate\tseq_h_len\tseq_v_len\n"
        << stats.alignmentLength << "\t"
        << stats.numGaps << "\t"
        << stats.numMatches << "\t"
        << stats.numMismatches << "\t"
        << stats.inBlockSeqH << "\t"
        << stats.inBlockSeqH << "\t"
        << stats.missedBlockSeqH << "\t"
        << stats.missedBlockSeqV << "\t"
        << ((stats.inBlockSeqH + stats.missedBlockSeqH) ? (100.0 * stats.inBlockSeqH / (stats.inBlockSeqH + stats.missedBlockSeqH)) : 0) << "\t"
        << ((stats.inBlockSeqV + stats.missedBlockSeqV) ? (100.0 * stats.inBlockSeqV / (stats.inBlockSeqV + stats.missedBlockSeqV)) : 0) << "\t"
        << (stats.inBlockSeqV ? (100.0 * (stats.numGaps + stats.numMismatches) / stats.inBlockSeqV) : 0) << "\t"
        << seqHLen << "\t"
        << seqVLen << "\n\n";

    out << "# alignment length:  " << stats.alignmentLength << "\n"
        << "# num gaps:          " << stats.numGaps << "\n"
        << "# num matches:       " << stats.numMatches << "\n"
        << "# num mismatches:    " << stats.numMismatches << "\n"
        << "# seq0 in block:     " << stats.inBlockSeqH << "\n"
        << "# seq1 in block:     " << stats.inBlockSeqV << "\n"
        << "# seq0 not in block: " << stats.missedBlockSeqH << "\n"
        << "# seq1 not in block: " << stats.missedBlockSeqV << "\n"
        << "\n"
        << "# covered 0 [%]:  " << ((stats.inBlockSeqH + stats.missedBlockSeqH) ? (100.0 * stats.inBlockSeqH / (stats.inBlockSeqH + stats.missedBlockSeqH)) : 0) << "\t"
        << "# covered 1 [%]:  " << ((stats.inBlockSeqV + stats.missedBlockSeqV) ? (100.0 * stats.inBlockSeqV / (stats.inBlockSeqV + stats.missedBlockSeqV)) : 0) << "\t"
        << "# error rate [%]: " << (stats.inBlockSeqV ? (100.0 * (stats.numGaps + stats.numMismatches) / stats.inBlockSeqV) : 0) << "\n"
        << "#\n"
        << "# seq 0 len:      " << seqHLen << "\n"
        << "# seq 1 len:      " << seqVLen << "\n"
        << "\n"
        << align << "\n";
}

// ===========================================================================
// Function computeStats()
// ===========================================================================

void computeStats(LocalBlockAlignerResult::Stats & result,
                  seqan::Align<seqan::Dna5String const> const & align,
                  std::vector<std::pair<unsigned, unsigned>> const & chunks)
{
    if (chunks.empty())
        return;

    // Compute statistics for each chunk.
    for (auto chunk : chunks)
    {
        result.alignmentLength += chunk.second - chunk.first;

        for (auto itH = iter(row(align, 0), chunk.first, seqan::Standard()),
                     itV = iter(row(align, 1), chunk.first, seqan::Standard()),
                     itHEnd = iter(row(align, 0), chunk.second, seqan::Standard());
             itH != itHEnd; ++itH, ++itV)
        {
            result.inBlockSeqH += !isGap(itH);
            result.inBlockSeqV += !isGap(itV);
            if (isGap(itH) || isGap(itV))
                result.numGaps += 1;
            else if (*itH == *itV)
                result.numMatches += 1;
            else
                result.numMismatches += 1;
        }
    }

    // Compute statistic between the chunks.
    for (auto it = chunks.begin(), itNext = std::next(chunks.begin()); itNext != chunks.end(); ++it, ++itNext)
    {
        result.missedBlockSeqH += toSourcePosition(row(align, 0), itNext->first) -
                toSourcePosition(row(align, 0), it->second);
        result.missedBlockSeqV += toSourcePosition(row(align, 1), itNext->first) -
                toSourcePosition(row(align, 1), it->second);
    }
}

// ===========================================================================
// Function main()
// ===========================================================================

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    LbaOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality was triggered then we exit the program.
    // The return code is 1 if there were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    if (options.verbosity >= 1)
        options.print(std::cout);

    // Load both sequences.
    seqan::CharString idH, idV;
    seqan::Dna5String seqH, seqV;
    seqan::SequenceStream seqIn(options.inSeq0.c_str());
    if (!isGood(seqIn))
        throw std::runtime_error("Could not open input file -i0");
    if (readRecord(idH, seqH, seqIn) != 0)
        throw std::runtime_error("Problem reading from input file -i0");

    open(seqIn, options.inSeq1.c_str());
    if (!isGood(seqIn))
        throw std::runtime_error("Could not open input file -i1");
    if (readRecord(idV, seqV, seqIn) != 0)
        throw std::runtime_error("Problem reading from input file -i1");

    // Get sequence ids from meta line and print in verbose mode.
    trimAfterSpace(idH);
    trimAfterSpace(idV);

    if (options.verbosity >= 2)
        std::cerr << "SEQ H\t" << idH << "\t" << seqH << "\n"
                  << "SEQ V\t" << idV << "\t" << seqV << "\n";

    LocalBlockAlignerResult result;
    auto & align = result.align;
    resize(rows(align), 2);
    setSource(row(align, 0), seqH);
    setSource(row(align, 1), seqV);

    seqan::AlignConfig<true, true, true, true> alignConfig;

    seqan::Score<int, seqan::Simple> scoringScheme(options.scoreMatch,
                                                   options.scoreMismatch,
                                                   options.scoreGapExtend,
                                                   options.scoreGapOpen);

    globalAlignment(align, scoringScheme, alignConfig);

    Counter leading = countLeading(align);
    Counter trailing = countTrailing(align);
    bool noAlignment = (leading.h + trailing.h + leading.v + trailing.v == length(seqH) + length(seqV));

    // Skip evaluation if seqV is a marker sequence.
    if (idV == "marker" || noAlignment)
    {
        result.stats.missedBlockSeqH = length(seqH);
        if (idV != "marker")
            result.stats.missedBlockSeqV = length(seqV);
    }
    else
    {
        if (options.verbosity >= 2)
            std::cerr << "INITIAL ALIGNMENT\n"
                      << align << "\n";

        // Clip away leading and traling gap.
        setClippedEndPosition(row(align, 0), length(row(align, 0)) - trailing.sum());
        setClippedBeginPosition(row(align, 0), leading.sum());
        setClippedEndPosition(row(align, 1), length(row(align, 1)) - trailing.sum());
        setClippedBeginPosition(row(align, 1), leading.sum());

        if (options.verbosity >= 2)
            std::cerr << "CLIPPED ALIGNMENT\n"
                      << align << "\n";

        auto chunks = computeChunks(align, options.longGapMinLen);
        computeStats(result.stats, align, chunks);

        clearClipping(row(align, 0));
        clearClipping(row(align, 1));

        // Add leading and trailing sequence to result.
        result.stats.missedBlockSeqH += leading.h + trailing.h;
        result.stats.missedBlockSeqV += leading.v + trailing.v;
    }

    // Print output.
    std::fstream f;
    if (!options.out.empty())
    {
        f.open(options.out.c_str(), std::ios::binary | std::ios::out);
        if (!f.good())
            throw std::runtime_error("Could not open output file.");
    }
    std::ostream & out = options.out.empty() ? std::cout : f;
    result.print(out);

    return 0;
}
