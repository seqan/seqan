// ==========================================================================
//                               FX Tools
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Computes coverage a genome given a SAM or BAM file, window-based.
// ==========================================================================

#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

// --------------------------------------------------------------------------
// Class FxBamCoverageOptions
// --------------------------------------------------------------------------

struct FxBamCoverageOptions
{
    // Verbosity level.  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Path to SAM file.
    seqan::CharString inBamPath;

    // Path to output file.
    seqan::CharString outPath;

    // Window size to use for computation.
    __int32 windowSize;

    FxBamCoverageOptions() : verbosity(1), windowSize(0)
    {}
};


// --------------------------------------------------------------------------
// Class Range
// --------------------------------------------------------------------------

struct Range
{
    // Reference id.
    __int32 rID;
    // Coordinates.
    __int32 beginPos;
    __int32 endPos;
    // Number of reads whose alignment cover this range.
    double coverage;
    
    seqan::CharString contigName;
    
    Range() :
        rID(-1),
        beginPos(0),
        endPos(0),
        coverage(0)
    {}

    Range(int32_t rID, int32_t beginPos, int32_t endPos, double coverage) :
        rID(rID),
        beginPos(beginPos),
        endPos(endPos),
        coverage(coverage)
    {}
};

void clear(Range & range, int rID)
{
    range.rID = rID;
    range.beginPos = 0;
    range.endPos = 0;
    range.coverage = 0;
}

std::ostream & operator<< (std::ostream & stream, Range const & range)
{
    stream << range.beginPos << '\t'
           << range.endPos << '\t'
           << range.coverage << '\n';
    return stream;
}

// --------------------------------------------------------------------------
// Function parseArgs()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseArgs(FxBamCoverageOptions & options,
          int argc,
          char const ** argv)
{
    seqan::ArgumentParser parser("fx_read_depth");
    setShortDescription(parser, "Read Coverage Computation.");
    setCategory(parser, "Utilities");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fB-o\\fP \\fIOUT.bed\\fP "
                 "\\fB-m\\fP \\fIMAPPING.bam\\fP");
    addDescription(parser, "Compute read depth, i.e. the number of reads covering a genomic position. "
                           "The result is a bed file with intervals and the maximal read depth in each interval. "
                           "Intervals start at multiples of the windows size, if given. "
                           "If no window size is given, positions with identical read depth are aggregated as intervals.");

    addOption(parser, seqan::ArgParseOption("m", "in-mapping", "Path to the mapping file to analyze.",
                                            seqan::ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "in-mapping", "sam bam");
    setRequired(parser, "in-mapping");

    // TODO(holtgrew): I want a custom help text!
    // addOption(parser, seqan::ArgParseOption("h", "help", "This helpful screen."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose, log to STDERR."));
    hideOption(parser, "verbose");
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose, log to STDERR."));
    hideOption(parser, "very-verbose");

    addSection(parser, "Main Options");
    addOption(parser, seqan::ArgParseOption("w", "window-size", "Set the size of the non-overlapping windows in base pairs.", seqan::ArgParseArgument::INTEGER, "NUM"));
    setDefaultValue(parser, "window-size", options.windowSize);

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o", "out-path", "Path to the resulting file.  If omitted, result is printed to stdout.", seqan::ArgParseArgument::OUTPUT_FILE, "OUTFILE"));
    setValidValues(parser, "out-path", ".bed");

    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
    {
        getOptionValue(options.inBamPath, parser, "in-mapping");
        getOptionValue(options.outPath, parser, "out-path");
        getOptionValue(options.windowSize, parser, "window-size");

        if (isSet(parser, "verbose"))
            options.verbosity = 2;
        if (isSet(parser, "very-verbose"))
            options.verbosity = 3;
    }

    return res;
}

template <typename TBamContext>
void writeRange(std::ostream & stream, Range & lastBin, Range const & range, FxBamCoverageOptions & options, TBamContext & bamContext, bool lastOfContig)
{
    SEQAN_ASSERT_EQ(lastBin.rID, range.rID);
    
    if (options.windowSize > 0)
    {
        unsigned binNo     = lastBin.endPos / options.windowSize;
        unsigned binNoEnd  = range.endPos / options.windowSize;
        unsigned binOfsEnd = range.endPos % options.windowSize;
        
        for (; binNo < binNoEnd; ++binNo)
        {
            // aggregate the max coverage
            if (lastBin.coverage < range.coverage)
                lastBin.coverage = range.coverage;
            lastBin.endPos = (binNo + 1) * options.windowSize;
            stream << contigNames(bamContext)[range.rID] << '\t' << lastBin;

            lastBin.beginPos = lastBin.endPos;
            lastBin.coverage = 0;
        }
        if (binOfsEnd != 0)
        {
            // aggregate the max coverage
            if (lastBin.coverage < range.coverage)
                lastBin.coverage = range.coverage;
            lastBin.endPos = range.endPos;
            
            if (lastOfContig)
                stream << contigNames(bamContext)[lastBin.rID] << '\t' << lastBin;
        }
    }
    else
    {
        stream << contigNames(bamContext)[range.rID] << '\t' << range;
    }
}


// ---------------------------------------------------------------------------
// Function main()
// ---------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    double startTime = seqan::sysTime();
    
    // -----------------------------------------------------------------------
    // Parse command line.
    // -----------------------------------------------------------------------
    FxBamCoverageOptions options;
    seqan::ArgumentParser::ParseResult res = parseArgs(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // -----------------------------------------------------------------------
    // Show options.
    // -----------------------------------------------------------------------
    if (options.verbosity >= 1)
    {
        std::cerr << "____OPTIONS___________________________________________________________________\n"
                  << "\n"
                  << "VERBOSITY    " << options.verbosity << "\n"
                  << "SAM/BAM      " << options.inBamPath << "\n"
                  << "OUT          " << options.outPath << "\n"
                  << "WINDOW SIZE  " << options.windowSize << "\n";
    }

    // -----------------------------------------------------------------------
    // Open Output
    // -----------------------------------------------------------------------

    std::ostream * out = &std::cout;
    std::ofstream outFile;
    if (!empty(options.outPath))
    {
        outFile.open(toCString(options.outPath), std::ios::binary | std::ios::out);
        if (!outFile.good())
        {
            std::cerr << "ERROR: Could not open output file " << options.outPath << "!\n";
            return 1;
        }
        out = &outFile;
    }

    // -----------------------------------------------------------------------
    // Compute Coverage
    // -----------------------------------------------------------------------

    std::cerr << "\n"
              << "___COVERAGE COMPUTATION________________________________________________________\n"
              << "\n"
              << "Computing Coverage...";

    seqan::BamFileIn bamFile;
    if (!open(bamFile, toCString(options.inBamPath)))
    {
        std::cerr << "Could not open " << options.inBamPath << "!\n";
        return 1;
    }

    seqan::BamHeader header;
    readHeader(header, bamFile);

    // Prepare bins.
    seqan::PriorityType<__int32, std::greater<__int32> > heap;
    
//    seqan::String<seqan::String<Range> > bins;
//    resize(bins, length(contigNames(context(bamFile))));
//    for (unsigned i = 0; i < length(bins); ++i)
//    {
//        resize(bins[i], (contigLengths(context(bamFile))[i] + options.windowSize - 1) / options.windowSize, options.windowSize);
//        back(bins[i]).length = (contigLengths(context(bamFile))[i] + options.windowSize - 1) % options.windowSize + 1;
//    }
//
    typedef int32_t TPos;
    
    int lastID = -1;
    TPos lastBreak = 0;
    int coverage = 0;
    seqan::BamAlignmentRecord record;
    Range lastBin;

    if (!atEnd(bamFile))
    {
        readRecord(record, bamFile);
        clear(lastBin, record.rID);
        lastID = record.rID;
    }
    
    while (!atEnd(bamFile))
    {
        while (!empty(heap) && top(heap) <= record.beginPos)
        {
            TPos endPos = top(heap);
            Range range(record.rID, lastBreak, endPos, coverage);
            writeRange(*out, lastBin, range, options, context(bamFile), lastID != record.rID);
            lastBreak = endPos;

            // remove all reads ending here and update coverage
            while (!empty(heap) && top(heap) == endPos)
            {
                --coverage;
                pop(heap);
            }
        }

        if (record.rID == lastID)
        {
            // output next range
            Range range(record.rID, lastBreak, record.beginPos, coverage);
            writeRange(*out, lastBin, range, options, context(bamFile), false);
        }
        else
        {
            SEQAN_ASSERT_EQ(coverage, 0u);
            
            // output last range
            if (lastBreak != contigLengths(context(bamFile))[record.rID])
            {
                Range range(lastID, lastBreak, contigLengths(context(bamFile))[lastID], 0u);
                writeRange(*out, lastBin, range, options, context(bamFile), true);
            }
            clear(lastBin, record.rID);
            if (record.beginPos != 0)
            {
                Range range(record.rID, 0, record.beginPos, 0u);
                writeRange(*out, lastBin, range, options, context(bamFile), false);
            }
                
            coverage = 0;
            lastID = record.rID;
        }
        lastBreak = record.beginPos;

        // consume all read at the current position until the next position
        while (true)
        {
            ++coverage;
            
            TPos len = 0;
            _getLengthInRef(len, record.cigar);
            push(heap, record.beginPos + len);

            while (!atEnd(bamFile))
            {
                readRecord(record, bamFile);

                if (!hasFlagUnmapped(record) && !hasFlagSecondary(record) && record.rID != seqan::BamAlignmentRecord::INVALID_REFID)
                    break;  // Only use aligned records.
            }
            
            if (atEnd(bamFile) || record.beginPos != lastBreak || record.rID != lastID)
                break;
            
        }
    }

    std::cerr << "DONE\n";

    if (options.verbosity >= 2)
        std::cerr << "Took " << (seqan::sysTime() - startTime) << " s\n";

    return 0;
}
