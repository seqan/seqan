// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// A simple converter between SAM and BAM.  BAM files can be dumped region-
// wise, using a .bai index.
// ==========================================================================

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/stream.h>

using namespace seqan;

struct Options
{
    bool showHelp;
    bool showVersion;
    CharString inFile;
    CharString baiFile;
    CharString outFile;
    unsigned verbosity;
    bool bamFormat;
    String<GenomicRegion> regions;  // <(refName, beginPos, endPos)>

    Options()
    {
        showHelp = false;
        showVersion = false;
        verbosity = 1;
        bamFormat = false;
    }

};

void
setupArgumentParser(ArgumentParser & parser, Options const & options)
{
    setVersion(parser, "1.0");

    addDescription(parser, "***********");
    addDescription(parser, "* BAMUTIL *");
    addDescription(parser, "***********");
    addDescription(parser, "");
    addDescription(parser, "SeqAn SAM/BAM Utility.");
    addDescription(parser, "");
    addDescription(parser, "Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>");

    addUsageLine(parser, "bamutil [OPTIONS] <IN >OUT");

    addSection(parser, "General Options");
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose mode."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose mode."));

    addSection(parser, "Input Specification");
    addOption(parser, ArgParseOption("i", "input-file", "Path to input, '-' for stdin.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("bi", "bai-index-file", "Path to BAI index, default: input + '.bai'.", ArgParseOption::STRING));

    addSection(parser, "Range Specification");
    addOption(parser, ArgParseOption("r", "region", "Regions to dump (in which aligments start). REF:FROM-TO, e.g. IV:1,000-2,000 will alignments with ref IV in range [1000,2000) (zero based).", ArgParseOption::STRING));

    addSection(parser, "Output Specification");
    addOption(parser, ArgParseOption("o", "output-file", "Path to output, '-' for stdout.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("b", "output-bam", "Output file is BAM."));
}

// Parse the command line and check for any syntatical errors.

ArgumentParser::ParseResult
parseArguments(Options & options,
               ArgumentParser & parser,
               int argc,
               char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValue(options.inFile, parser, "input-file");
    getOptionValue(options.baiFile, parser, "bai-index-file");
    if (empty(options.baiFile) && !empty(options.inFile))
    {
        options.baiFile = options.inFile;
        append(options.baiFile, ".bai");
    }

    getOptionValue(options.outFile, parser, "output-file");
    options.bamFormat = isSet(parser, "output-bam");
    if (isSet(parser, "region"))
    {
        String<CharString> regions = getOptionValues(parser, "region");
        GenomicRegion region;
        for (unsigned i = 0; i < length(regions); ++i)
        {
            try
            {
                parse(region, regions[i]);
                region.beginPos++;
                appendValue(options.regions, region);
                if (options.verbosity >= 2)
                    std::cerr << "[VERBOSE] Region " << region.seqName << ":"
                              << region.beginPos << "-" << region.endPos << std::endl;
            }
            catch (ParseError & e)
            {
                std::cerr << "[WARNING] could not parse region \"" << regions[i] << "\". IGNORING." << std::endl;
            }
        }
    }

    return ArgumentParser::PARSE_OK;
}

template <typename TOptions>
int _dumpRegion(BamFileIn & in, BamFileOut & out, TOptions const & options)
{
    // TOOD(holtgrew): The index is loaded for each region.  This should probably not be the case!

    // Dump each region after loading the index.
    BamIndex<Bai> bamIndex;
    if (!open(bamIndex, toCString(options.baiFile)))
    {
        std::cerr << "[ERROR] Could not open index file " << options.baiFile << ", required when specifying regions." << std::endl;
        return 1;
    }

    for (unsigned i = 0; i < length(options.regions); ++i)
    {
        // Jump near range.
        CharString refName = options.regions[i].seqName;
        unsigned refId = 0;
        if (!getIdByName(refId, contigNamesCache(context(in)), refName))
        {
            std::cerr << "[ERROR] Unknown reference " << refName << std::endl;
            return 1;
        }
        bool hasAlignments = false;
        int beginPos = options.regions[i].beginPos;
        int endPos = options.regions[i].endPos;
        if (!jumpToRegion(in, hasAlignments, refId, beginPos, endPos, bamIndex))
        {
            std::cerr << "[ERROR] Could not jump to " << refName << ":" << beginPos << "-" << endPos << std::endl;
            return 1;
        }
        if (!hasAlignments)
            continue;
        // Dump range.
        BamAlignmentRecord record;
        while (!atEnd(in))
        {
            readRecord(record, in);
            if (record.beginPos < options.regions[i].beginPos)
                continue;  // Skip, before region.
            if (record.beginPos >= options.regions[i].endPos)
                return 0;

            writeRecord(out, record);
        }
    }
    return 0;
}

template <typename TOptions>
int performConversion(BamFileIn & in, BamFileOut & out, TOptions const & options)
{
    BamHeader header;
    readHeader(header, in);
    writeHeader(out, header);

    if (length(options.regions) == 0u)
    {
        // Simply dump whole file.
        BamAlignmentRecord record;
        while (!atEnd(in))
        {
            readRecord(record, in);
            writeRecord(out, record);
        }
        return 0;
    }
    else
    {
        return _dumpRegion(in, out, options);
    }
}

// The main function.
//
// Don't be intimidated by its longness, most of the code is for branching the different options to types.

int main(int argc, char const * argv[])
{
    // -----------------------------------------------------------------------
    // Handle Command Line
    // -----------------------------------------------------------------------

    // Setup command line parser.
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    ArgumentParser::ParseResult res = parseArguments(options, parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // -----------------------------------------------------------------------
    // Open Input / Output Files.
    // -----------------------------------------------------------------------

    BamFileIn reader;
    BamFileOut writer(reader);
    bool success = true;

    if (!empty(options.inFile))
        success = open(reader, toCString(options.inFile));
    else
        success = open(reader, std::cin);

    if (!success)
    {
        std::cerr << "Couldn't open " << options.inFile << " for reading." << std::endl;
        return 1;
    }

    if (!empty(options.outFile))
        success = open(writer, toCString(options.outFile));
    else
    {
        // write to stdout
#if SEQAN_HAS_ZLIB
        if (options.bamFormat)
            success = open(writer, std::cout, Bam());
        else
#endif
        success = open(writer, std::cout, Sam());
    }

    if (!success)
    {
        std::cerr << "Couldn't open " << options.outFile << " for writing." << std::endl;
        return 1;
    }

    // (weese:) We go through no pain here anymore.

    if (performConversion(reader, writer, options) != 0)
    {
        std::cerr << "Error during conversion!" << std::endl;
        return 1;
    }

    return 0;
}
