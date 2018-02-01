// ==========================================================================
//                                 FX Tools
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

#include <map>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The input file name is a string.
    seqan::CharString inFilename;

    // The out file name is an out file.
    seqan::CharString outFilename;

    AppOptions() :
        verbosity(1)
    {}
};

// --------------------------------------------------------------------------
// Class FastqStats
// --------------------------------------------------------------------------

struct FastqStats
{
    // -----------------------------------------------------------------------
    // Members with Results
    // -----------------------------------------------------------------------

    unsigned maxLength;
    
    // Number of bases in column i.
    seqan::String<int64_t> numBases;
    // Smallest score in column i.
    seqan::String<int32_t> minScores;
    // Largest score in column i.
    seqan::String<int32_t> maxScores;
    // Sum of scores in column i.
    seqan::String<int64_t> sumScores;
    // Mean of scores in column i.
    seqan::String<double> meanScores;
    // First quartile quality score (Q1) for column i.
    seqan::String<double> firstQuartiles;
    // Median quality score for column i.
    seqan::String<double> medianScores;
    // Third quartile quality score (Q3) for column i.
    seqan::String<double> thirdQuartiles;
    // Inter-quartile range  (Q3-Q1) for column i.
    seqan::String<double> interQuartileRanges;
    // Left-whisker value for boxplotting for column i.
    seqan::String<int32_t> leftWhiskers;
    // Right-whisker value for boxplotting for column i.
    seqan::String<int32_t> rightWhiskers;
    // Number of nucleotides A, C, G, T, N for column i.
    seqan::String<seqan::String<int64_t> > nucleotideCounts;

    // -----------------------------------------------------------------------
    // Histogram Members
    // -----------------------------------------------------------------------

    // Quality histogram.
    seqan::String<std::map<int32_t, int32_t> > qualHistos;

    // -----------------------------------------------------------------------
    // Constructor
    // -----------------------------------------------------------------------

    FastqStats() : maxLength(0)
    {}

    // -----------------------------------------------------------------------
    // Member Functions
    // -----------------------------------------------------------------------

    // Resize members to read of length.
    void resizeToReadLength(unsigned n)
    {
        if (maxLength == n)
            return;
        maxLength = n;

        resize(numBases, n, 0);
        resize(minScores, n, 80);
        resize(maxScores, n, 0);
        resize(sumScores, n, 0);
        resize(meanScores, n, 0);
        resize(firstQuartiles, n, 0);
        resize(medianScores, n, 0);
        resize(thirdQuartiles, n, 0);
        resize(interQuartileRanges, n, 0);
        resize(leftWhiskers, n, 0);
        resize(rightWhiskers, n, 0);
        resize(nucleotideCounts, n);
        for (unsigned i = 0; i < n; ++i)
            resize(nucleotideCounts[i], 5, 0);

        resize(qualHistos, n);
    }

    // Update histogram and statistics for the given read.
    void registerRead(seqan::Dna5String const & seq, seqan::CharString const & quals)
    {
        resizeToReadLength(length(seq));

        // Update nucleotide counts.
        for (unsigned i = 0; i < length(seq); ++i)
        {
            numBases[i] += 1;
            nucleotideCounts[i][ordValue(seq[i])] += 1;
        }

        // Update quality histograms and statistics.
        for (unsigned i = 0; i < length(quals); ++i)
        {
            int qual = quals[i] - '!';  // PHRED scores.
            if (numBases[i] == 0u || minScores[i] > qual)
                minScores[i] = qual;
            if (numBases[i] == 0u || maxScores[i] < qual)
                maxScores[i] = qual;
            qualHistos[i][qual] += 1;
            sumScores[i] += qual;
        }
    }

    // Compute statistics after updating for the last read.
    void finalizeStats()
    {
        // Compute score means.
        for (unsigned i = 0; i < length(meanScores); ++i)
            meanScores[i] = (1.0 * sumScores[i]) / numBases[i];
        // Compute score medians and quartiles.
        for (unsigned i = 0; i < length(qualHistos); ++i)
        {
            if (qualHistos[i].size() == 0u)
                continue;  // Skip if empty.

            // TODO(holtgrew): Quartile/median computation not mathematically correct yet.
            unsigned n = numBases[i];  // Number of bases.
            unsigned firstQN = n / 4;
            unsigned medianN = n / 2;
            unsigned thirdQN = (3 * n) / 4;
            unsigned count = 0;  // Number of bases up to here.
            for (std::map<int32_t, int32_t>::const_iterator it = qualHistos[i].begin(); it != qualHistos[i].end(); ++it)
            {
                if (count < firstQN && count + it->second >= firstQN)
                    firstQuartiles[i] = it->first;
                if (count < medianN && count + it->second >= medianN)
                    medianScores[i] = it->first;
                if (count < thirdQN && count + it->second >= thirdQN)
                    thirdQuartiles[i] = it->first;

                count += it->second;
            }
            interQuartileRanges[i] = (thirdQuartiles[i] - firstQuartiles[i]);
        }

        // Compute whiskers as the data point that is still within 1.5 IQR of the lower quartile.
        for (unsigned i = 0; i < length(qualHistos); ++i)
        {
            if (qualHistos[i].size() == 0u)
                continue;  // Skip if empty.
            leftWhiskers[i] = static_cast<int>(firstQuartiles[i]);
            rightWhiskers[i] = static_cast<int>(thirdQuartiles[i]);
            double leftWhiskerBound = ((double)firstQuartiles[i]) - 1.5 * interQuartileRanges[i];
            double rightWhiskerBound = ((double)thirdQuartiles[i]) + 1.5 * interQuartileRanges[i];
            for (std::map<int32_t, int32_t>::const_iterator it = qualHistos[i].begin(); it != qualHistos[i].end(); ++it)
                if (leftWhiskers[i] > it->first && it->first >= leftWhiskerBound)
                    leftWhiskers[i] = it->first;
            for (std::map<int32_t, int32_t>::const_iterator it = qualHistos[i].begin(); it != qualHistos[i].end(); ++it)
                if (rightWhiskers[i] < it->first && it->first <= rightWhiskerBound)
                    rightWhiskers[i] = it->first;
        }
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("fx_fastq_stats");
    // Set short description, version, and date.
    setShortDescription(parser, "Compute FASTQ statistics.");
    setCategory(parser, "NGS Quality Control");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    // Define usage line and long description.
    addUsageLine(parser, "\\fB-i\\fP \\fIINPUT.fq\\fP \\fB-o\\fP \\fIOUTPUT.fq.stats.tsv\\fP");
    addDescription(parser,
                   "Read a FASTQ file with reads sequences of the same length.  Writes out a TSV file with one record for "
                   "each column/position with statistics on the nucleotides and qualities.");

    addSection(parser, "Input / Output");
    addOption(parser, seqan::ArgParseOption("i", "input", "Input FASTQ file.", seqan::ArgParseOption::INPUT_FILE, "INPUT"));
    setValidValues(parser, "input", "fastq fq");
    setRequired(parser, "input");
    addOption(parser, seqan::ArgParseOption("o", "output", "Output TSV file.", seqan::ArgParseOption::OUTPUT_FILE, "OUTPUT"));
    setRequired(parser, "output");
    setValidValues(parser, "output", "fq_stats_tsv");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    seqan::getOptionValue(options.inFilename, parser, "input");
    seqan::getOptionValue(options.outFilename, parser, "output");

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Open input file.
    seqan::SeqFileIn seqFile;
    if (!open(seqFile, toCString(options.inFilename)))
    {
        std::cerr << "ERROR: Could not open file " << options.inFilename << " for reading.\n";
        return 1;
    }

    // Read sequences and build result.
    FastqStats stats;
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString quals;
    while (!atEnd(seqFile))
    {
        readRecord(id, seq, quals, seqFile);
        if (empty(quals))  // Fill with Q40 if there are no qualities.
            resize(quals, length(seq), '!' + 40);

        // Update statistics.
        stats.registerRead(seq, quals);
    }

    // Finalize statistics and write to output.
    stats.finalizeStats();

    std::ostream * out = &std::cout;
    std::ofstream outStream;
    if (options.outFilename != "-")
    {
        outStream.open(toCString(options.outFilename), std::ios::binary | std::ios::out);
        if (!outStream.good())
        {
            std::cerr << "ERROR: Could not open file " << options.outFilename << " for writing.\n";
            return 1;
        }
        out = &outStream;
    }

    *out << "#column\tcount\tmin\tmax\tsum\tmean\tQ1\tmedian\tQ3\tIQR\tlW\trW\tA_count\tC_count\tG_count\tT_count\tN_count\n";
    for (unsigned i = 0; i < stats.maxLength; ++i)
    {
        *out << i << "\t"
             << stats.numBases[i] << "\t"
             << stats.minScores[i] << "\t"
             << stats.maxScores[i] << "\t"
             << stats.sumScores[i] << "\t"
             << stats.meanScores[i] << "\t"
             << stats.firstQuartiles[i] << "\t"
             << stats.medianScores[i] << "\t"
             << stats.thirdQuartiles[i] << "\t"
             << stats.interQuartileRanges[i] << "\t"
             << stats.leftWhiskers[i] << "\t"
             << stats.rightWhiskers[i] << "\t"
             << stats.nucleotideCounts[i][0] << "\t"
             << stats.nucleotideCounts[i][1] << "\t"
             << stats.nucleotideCounts[i][2] << "\t"
             << stats.nucleotideCounts[i][3] << "\t"
             << stats.nucleotideCounts[i][4] << "\n";
    }

    return 0;
}
