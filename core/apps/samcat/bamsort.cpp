// ==========================================================================
//                                   bamsort
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#define _GLIBCXX_PARALLEL

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/arg_parse.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

namespace seqan {

template <typename TDirection>
struct FileFormat<SmartFile<Bam, TDirection> >
{
    typedef Bam Type;
};

}

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    CharString  inFile;
    CharString  outPrefix;
    CharString  outFile;
    size_t      maxMem;

    bool        verbose;
    __uint64    numRecords;

    AppOptions() :
        maxMem(512 * 1024 * 1024),      // default memory to use is 512MB
        numRecords(0)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

struct LessCoord
{
    char const *start;

    LessCoord(CharString &buffer) :
        start(begin(buffer, Standard()) + 4)
    {}

    bool operator() (size_t a, size_t b)
    {
        return *reinterpret_cast<const __uint64 *>(start + a) < *reinterpret_cast<const __uint64 *>(start + b);
    }
};

struct LessCoordMerge
{
    String<char const *> &buffers;

    LessCoord(String<char const *> &buffers) :
        buffers(buffers)
    {}

    bool operator() (size_t a, size_t b)
    {
        return *reinterpret_cast<const __uint64 *>(buffers[a] + 4) < *reinterpret_cast<const __uint64 *>(buffers[b] + 4);
    }
};

// --------------------------------------------------------------------------
// Function sortBamInChunks()
// --------------------------------------------------------------------------

template <typename TLess>
bool sortChunks(String<CharString> &outFiles, BamFileIn &bamFileIn, AppOptions const &options)
{
    BamHeader header;
    readRecord(header, bamFileIn);

    CharString buffer;
    String<size_t> ofs;
    reserve(buffer, options.maxMem, Exact());

    __int32 recordLen = -1;
    bool needMerge = false;
    bool doReadRecordLen = true;

    while (!atEnd(bamFileIn))
    {
        clear(ofs);
        clear(buffer);

        while (!atEnd(bamFileIn))
        {
            appendValue(ofs, length(buffer));                                   // Remember begin position of record in buffer.
            if (doReadRecordLen)
                readRawPod(recordLen, bamFileIn.iter);                          // Read size of the remaining block.

            if (length(buffer) + 4 + recordLen > options.maxMem)
            {
                needMerge = true;
                doReadRecordLen = false;
                break;
            }

            appendRawPod(buffer, recordLen);                                    // Append record size.
            write(buffer, bamFileIn.iter, recordLen);                           // Read remaining block in one chunk (fastest).

            doReadRecordLen = true;
            options.numRecords++;
        }

        std::sort(begin(ofs, Standard()), end(ofs, Standard()), TLess(buffer)); // Parallely sort records

        CharString fname = options.outFile;
        if (needMerge)                                                          // Get output file name.
        {
            fname = options.outPrefix;
            appendValue(fname, '.');
            appendNumber(fname, length(outFiles));
            append(fname, ".bam");
            appendValue(outFiles, fname);
        }

        BamFileOut bamFileOut(bamFileIn);                                       // Open output/temporary bam files.
        bool success;
        if (!empty(options.outFile) || needMerge)
            success = open(bamFileOut, toCString(fname));
        else
            success = open(bamFileOut, std::cout, Bam());
        if (!success)
        {
            std::cerr << "Couldn't open " << fname << " for writing." << std::endl;
            return false;
        }

        writeRecord(bamFileOut, header);                                        // Write header.

        typedef Iterator<String<size_t>, Standard>::Type TOfsIter;              // Write records in sorted order.
        TOfsIter it = begin(ofs, Standard());
        TOfsIter itEnd = end(ofs, Standard());
        char const *start = begin(buffer, Standard());
        for (; it != itEnd; ++it)
            write(bamFileOut.iter, start + *it, *reinterpret_cast<const __uint32 *>(start + *it));
    }

    return true;
}

// --------------------------------------------------------------------------
// Function mergeBamFiles()
// --------------------------------------------------------------------------

template <typename TLess>
bool mergeBamFiles(BamFileOut &bamFileOut, String<CharString> &chunkFiles)
{
    String<BamFileIn *> readerPtr;
    resize(readerPtr, length(chunkFiles));

    // Step 1: Merge all headers (if available)
    BamHeader header;
    for (unsigned i = 0; i < length(chunkFiles); ++i)
    {
        readerPtr[i] = new BamFileIn(bamFileOut);
        if (!open(*readerPtr[i], toCString(chunkFiles[i])))
        {
            std::cerr << "Couldn't open " << toCString(chunkFiles[i]) << " for reading." << std::endl;
            return false;
        }
        readRecord(header, *(readerPtr[i]));
    }

    // Step 2: Remove duplicate header entries and write merged header
    removeDuplicates(header);
    writeRecord(bamFileOut, header);

    // Step 3: Read and output alignment records
    String<CharString> buffers;
    std::priority_queue<size_t> queue(TLess(buffers));
    for (unsigned i = 0; i != length(chunkFiles); ++i)
        queue.push_back(i);

    BamAlignmentRecord record;
    while (!queue.empty())
    {

        readRawPod(recordLen, bamFileIn.iter);                          // Read size of the remaining block.


        write(buffer, bamFileIn.iter, recordLen);                           // Read remaining block in one chunk (fastest).

        // copy all alignment records
        while (!atEnd(*readerPtr[i]))
        {

            readRecord(record, *readerPtr[i]);
            writeRecord(bamFileOut, record);
            ++numRecords;
        }
    }

    for (unsigned i = 0; i < length(chunkFiles); ++i)
    {
        close(*readerPtr[i]);
        delete readerPtr[i];
        unlink(toCString(chunkFiles[i]));
    }
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("bamsort");
    // Set short description, version, and date.
    setCategory(parser, "Utilities");
    setVersion(parser, "0.1");
    setDate(parser, "September 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIINFILE\\fP> [-o <\\fIOUTFILE\\fP>]");
    setShortDescription(parser, "BAM file sorter");
    addDescription(parser, "This tool reads a BAM file and outputs a sorted BAM file "
                           "and outputs the concatenation of them. "
                           "If the output file name is ommitted the result is written to stdout.");

    addDescription(parser, "(c) Copyright 2014 by David Weese.");

    addOption(parser, ArgParseOption("o", "output", "Output file name", ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "output", BamFileOut::getFileFormatExtensions());

    // We require one argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "INFILE"));
    setValidValues(parser, 0, BamFileIn::getFileFormatExtensions());
    setHelpText(parser, 0, "Input BAM file (or - for stdin).");
    addOption(parser, ArgParseOption("v", "verbose", "Print some stats."));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBbamsort\\fP \\fBunsorted.bam\\fP \\fB-o\\fP \\fBsorted.bam\\fP",
                "Sort a bam file.");
    addListItem(parser, "\\fBbamsort\\fP \\fB-\\fP",
                "Sort a bam stream from stdin and write to stdout.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.inFile, parser, 0);
    getOptionValue(options.outFile, parser, "output");
    getOptionValue(options.verbose, parser, "verbose");

    if (length(options.outFile) > 4)
        options.outPrefix = prefix(options.outFile, length(options.outFile) - 4);
    else if (length(options.inFile) > 4)
        options.outPrefix = prefix(options.inFile, length(options.inFile) - 4);
    else
        options.outPrefix = "temp";

    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    double start = sysTime();

    // Step 1: Open input file
    BamFileIn bamFileIn;
    bool success;
    if (!empty(options.inFile))
        success = open(bamFileIn, toCString(options.inFile));
    else
        success = open(bamFileIn, std::cin);

    if (!success)
    {
        std::cerr << "Couldn't open " << options.inFile << " for reading." << std::endl;
        return 1;
    }

    // Step 2: Sort chunks in memory and write to disk
    String<CharString> chunkFiles;
    if (!sortChunks<LessCoord>(chunkFiles, bamFileIn, options))
        return 1;

    // Step 3: Merge chunks and write to output file
    if (!empty(chunkFiles))
    {
        BamFileOut bamFileOut;
        if (!empty(options.outFile))
            success = open(bamFileOut, toCString(options.outFile));
        else
            success = open(bamFileOut, std::cout);

        if (!success)
        {
            std::cerr << "Couldn't open " << options.outFile << " for writing." << std::endl;
            return 1;
        }

        if (!mergeBamFiles<LessCoordMerge>(bamFileOut, chunkFiles))
            return 1;
    }

    double stop = sysTime();
    if (options.verbose)
    {
        std::cerr << "Number of alignments: " << options.numRecords << std::endl;
        std::cerr << "Elapsed time:         " << stop - start << " seconds" << std::endl;
    }
    return 0;
}
