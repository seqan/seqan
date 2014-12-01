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

#ifdef _OPENMP
#define _GLIBCXX_PARALLEL
#else
#if SEQAN_IGNORE_MISSING_OPENMP != 1
#pragma message("OpenMP not found! Parallelization will be limited in bamsort.")
#endif
#endif

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
struct FileFormat<SmartFile<Bam, TDirection, int> >
{
    typedef Bam Type;
};

typedef SmartFile<Bam, Input, int> BamOnlyFileIn;
typedef SmartFile<Bam, Output, int> BamOnlyFileOut;

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
    CharString  order;
    __uint64    maxMem;

    __uint64    numRecords;
    bool        verbose;

    AppOptions() :
        order("coord"),
        maxMem(512 << 20),      // default memory to use is 512MB
        numRecords(0),
        verbose(false)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

struct BamInStream
{
    std::ifstream           file;
    BamOnlyFileIn           reader;
    String<unsigned char>   buffer;

    template <typename TObject>
    BamInStream(TObject &object) :
        reader(object)
    {}
};

// compare strings lexicographically and runs of digits as numerically
inline int
numericalStrCmp(const unsigned char * SEQAN_RESTRICT a,
                const unsigned char * SEQAN_RESTRICT b)
{
    IsDigit isDigit;
    while (*a && *b)
    {
        if (isDigit(*a) && isDigit(*b))
        {
            while (*a == '0') ++a;                          // skip leading zeros
            while (*b == '0') ++b;                          // dito
            while (isDigit(*a) && *a == *b) { ++a; ++b; };  // skip until first difference

            int tmp = 0;
            if (isDigit(*a) && isDigit(*b))
            {
                tmp = (int)*a - (int)*b;
                while (isDigit(*a) && isDigit(*b)) { ++a; ++b; };
            }
            if (isDigit(*a))
                return 1;                                   // a has more non-zero digits than b
            else
                if (isDigit(*b))
                    return -1;                              // b has more non-zero digits than a
                else if (tmp != 0)
                    return tmp;                             // same amount of digits, the first difference (if exists) decides
        }
        else if (*a != *b)
            return (int)*a - (int)*b;
        else
            ++a, ++b;
    }
    return (int)*a - (int)*b;
}

// sort by genomic coordinate, query name
struct LessCoord
{
    unsigned char const *start;

    LessCoord() :
        start(NULL)
    {}

    LessCoord(String<unsigned char> &buffer) :
        start(begin(buffer, Standard()) + 4)
    {}

    bool operator() (size_t a, size_t b)
    {
        __uint64 x = *reinterpret_cast<const __uint64 *>(start + a);
        __uint64 y = *reinterpret_cast<const __uint64 *>(start + b);
        x = x << 32 | x >> 32;
        y = y << 32 | y >> 32;
        if (x != y)
            return x < y;

        return numericalStrCmp(start + a + sizeof(BamAlignmentRecordCore),
                               start + b + sizeof(BamAlignmentRecordCore)) < 0;
    }

    bool operator() (BamInStream *a, BamInStream *b)
    {
        __uint64 x = *reinterpret_cast<const __uint64 *>(&a->buffer[0]);
        __uint64 y = *reinterpret_cast<const __uint64 *>(&b->buffer[0]);
        x = x << 32 | x >> 32;
        y = y << 32 | y >> 32;
        if (x != y)
            return x > y;

        return numericalStrCmp(&a->buffer[sizeof(BamAlignmentRecordCore)],
                               &b->buffer[sizeof(BamAlignmentRecordCore)]) > 0;
    }
};

// sort by query name, genomic coordinate
struct LessQName
{
    unsigned char const *start;

    LessQName() :
        start(NULL)
    {}

    LessQName(String<unsigned char> &buffer) :
        start(begin(buffer, Standard()) + 4)
    {}

    bool operator() (size_t a, size_t b)
    {
        int res = numericalStrCmp(start + a + sizeof(BamAlignmentRecordCore),
                                  start + b + sizeof(BamAlignmentRecordCore));
        if (res != 0)
            return res < 0;

        __uint64 x = *reinterpret_cast<const __uint64 *>(start + a);
        __uint64 y = *reinterpret_cast<const __uint64 *>(start + b);
        x = x << 32 | x >> 32;
        y = y << 32 | y >> 32;
        return x < y;
    }

    bool operator() (BamInStream *a, BamInStream *b)
    {
        int res = numericalStrCmp(&a->buffer[sizeof(BamAlignmentRecordCore)],
                                  &b->buffer[sizeof(BamAlignmentRecordCore)]);
        if (res != 0)
            return res > 0;

        __uint64 x = *reinterpret_cast<const __uint64 *>(&a->buffer[0]);
        __uint64 y = *reinterpret_cast<const __uint64 *>(&b->buffer[0]);
        x = x << 32 | x >> 32;
        y = y << 32 | y >> 32;
        return x > y;
    }
};

// --------------------------------------------------------------------------
// Function sortBamInChunks()
// --------------------------------------------------------------------------

template <typename TLess>
bool sortChunks(String<CharString> &outFiles, BamOnlyFileIn &bamFileIn, AppOptions &options)
{
    BamHeader header;
    readRecord(header, bamFileIn);

    if (options.order == "coord")
        setSortOrder(header, BAM_SORT_COORDINATE);
    else
        setSortOrder(header, BAM_SORT_QUERYNAME);

    String<unsigned char> buffer;
    String<size_t> ofs;
    reserve(buffer, options.maxMem, Exact());

    __uint32 recordLen = 0;
    bool needMerge = false;
    bool doReadRecordLen = true;

    while (!atEnd(bamFileIn))
    {
        clear(ofs);
        clear(buffer);
        double start = sysTime();

        if (options.verbose)
            std::cerr << "Reading chunk #" << length(outFiles) << std::flush;

        while (!atEnd(bamFileIn))
        {
            if (doReadRecordLen)
                readRawPod(recordLen, bamFileIn.iter);                          // Read size of the remaining block.

            if (length(buffer) + 4 + recordLen > options.maxMem)
            {
                needMerge = true;
                doReadRecordLen = false;
                break;
            }

            appendValue(ofs, length(buffer));                                   // Remember begin position of record in buffer.
            appendRawPod(buffer, recordLen);                                    // Append record size.
            write(buffer, bamFileIn.iter, recordLen);                           // Read remaining block in one chunk (fastest).

            doReadRecordLen = true;
            options.numRecords++;
        }

        if (options.verbose)
        {
            std::cerr << " took\t" << sysTime() - start << " seconds" << std::endl;
            std::cerr << "Sorting chunk #" << length(outFiles) << std::flush;
            start = sysTime();
        }

        std::sort(begin(ofs, Standard()), end(ofs, Standard()), TLess(buffer)); // Parallely sort records

        if (options.verbose)
        {
            std::cerr << " took\t" << sysTime() - start << " seconds" << std::endl;
            std::cerr << "Writing chunk #" << length(outFiles) << std::flush;
            start = sysTime();
        }

        CharString fname = options.outFile;
        if (needMerge)                                                          // Get output file name.
        {
            fname = options.outPrefix;
            appendValue(fname, '.');
            appendNumber(fname, length(outFiles));
            append(fname, ".bam");
            appendValue(outFiles, fname);
        }

        std::ofstream rawFile;
        BamOnlyFileOut bamFileOut(bamFileIn);                                       // Open output/temporary bam files.
        bool success;
/*        if (needMerge)
        {
            rawFile.open(toCString(fname));
            success = rawFile.is_open() && _open(bamFileOut, rawFile, Nothing(), False());
        }
        else*/ if (!empty(options.outFile) || needMerge)
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
        unsigned char const *beg = begin(buffer, Standard());
        for (; it != itEnd; ++it)
            write(bamFileOut.iter, beg + *it, 4 + *reinterpret_cast<const __uint32 *>(beg + *it));
        close(bamFileOut);

        if (options.verbose)
            std::cerr << " took\t" << sysTime() - start << " seconds" << std::endl;
    }

    return true;
}

// --------------------------------------------------------------------------
// Function mergeBamFiles()
// --------------------------------------------------------------------------

template <typename TLess>
bool mergeBamFiles(BamOnlyFileOut &bamFileOut, String<CharString> &chunkFiles, AppOptions &options)
{
    String<BamInStream *> streamPtr;
    resize(streamPtr, length(chunkFiles));

    double start = sysTime();
    if (options.verbose)
        std::cerr << "Merging chunks" << std::flush;

    // Step 1: Merge all headers (if available)
    BamHeader header;
    for (unsigned i = 0; i < length(chunkFiles); ++i)
    {
        streamPtr[i] = new BamInStream(bamFileOut);
        streamPtr[i]->file.open(toCString(chunkFiles[i]));
        if (!streamPtr[i]->file.is_open() || !open(streamPtr[i]->reader, streamPtr[i]->file))
        {
            std::cerr << "Couldn't open " << toCString(chunkFiles[i]) << " for reading." << std::endl;
            return false;
        }
        readRecord(header, streamPtr[i]->reader);
    }

    // Step 2: Remove duplicate header entries and write merged header
    removeDuplicates(header);
    writeRecord(bamFileOut, header);

    // Step 3: Read first records and fill priority queue
    String<CharString> buffers;
    std::priority_queue<BamInStream *, std::vector<BamInStream *>, TLess> queue;
    for (unsigned i = 0; i != length(chunkFiles); ++i)
        if (!atEnd(streamPtr[i]->reader))
        {
            _readBamRecordWithoutSize(streamPtr[i]->buffer, streamPtr[i]->reader.iter);
            queue.push(streamPtr[i]);
        }

    // Step 4: Extract smallest element, write it and read next from the same stream
    BamAlignmentRecord record;
    while (!queue.empty())
    {
        BamInStream *top = queue.top();
        appendRawPod(bamFileOut.iter, (__uint32)length(top->buffer));
        write(bamFileOut.iter, top->buffer);
        queue.pop();
        if (!atEnd(top->reader))
        {
            _readBamRecordWithoutSize(top->buffer, top->reader.iter);
            queue.push(top);
        }
    }

    for (unsigned i = 0; i < length(chunkFiles); ++i)
    {
        close(streamPtr[i]->reader);
        delete streamPtr[i];
    }

    if (options.verbose)
        std::cerr << " took\t" << sysTime() - start << " seconds" << std::endl;
    return true;
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

    addDescription(parser, "(c) Copyright in 2014 by David Weese.");

    addOption(parser, ArgParseOption("o", "output", "Output file name.", ArgParseOption::OUTPUT_FILE));
    setValidValues(parser, "output", BamOnlyFileOut::getFileExtensions());

    // We require one argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "INFILE"));
    setValidValues(parser, 0, BamOnlyFileIn::getFileExtensions());
    setHelpText(parser, 0, "Input BAM file (or - for stdin).");
	addOption(parser, ArgParseOption("s", "sort-order", "Sort by either reference coordinate or query name.", ArgParseOption::STRING));
    setValidValues(parser, "sort-order", "coord qname");
    setDefaultValue(parser, "sort-order", options.order);
	addOption(parser, ArgParseOption("m", "max-memory", "Set maximal amount of memory (in MB) used for buffering.", ArgParseOption::INTEGER));
    setMinValue(parser, "max-memory", "1");
    setDefaultValue(parser, "max-memory", options.maxMem >> 20);
    addOption(parser, ArgParseOption("v", "verbose", "Print some stats."));

    // Add explanation for sort order
    addTextSection(parser, "Sort Order");
    addText(parser,
            "The three sort orders \\fIcoord\\fP, \\fIqname\\fP, and \\fIqname_lex\\fP have the following meaning.");
    addListItem(parser, "coord",
                "Sort by coordinate, unmapped reads come after all others. Secondary key is qname.");
    addListItem(parser, "qname",
                "Sort by query names lexicographically and included number numerically. Secondary key is coordinate. ");
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
    getOptionValue(options.order, parser, "sort-order");
	getOptionValue(options.maxMem, parser, "max-memory");
    options.maxMem <<= 20;
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
    BamOnlyFileIn bamFileIn;
    bool success;
    if (options.inFile != "-")
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
    if (options.order == "coord")
        success = sortChunks<LessCoord>(chunkFiles, bamFileIn, options);
    else
        success = sortChunks<LessQName>(chunkFiles, bamFileIn, options);

    // Step 3: Merge chunks and write to output file
    if (!empty(chunkFiles) && success)
    {
        BamOnlyFileOut bamFileOut;
        if (!empty(options.outFile))
            success = open(bamFileOut, toCString(options.outFile));
        else
            success = open(bamFileOut, std::cout);

        if (success)
        {
            if (options.order == "coord")
                success = mergeBamFiles<LessCoord>(bamFileOut, chunkFiles, options);
            else
                success = mergeBamFiles<LessQName>(bamFileOut, chunkFiles, options);
        }
        else
            std::cerr << "Couldn't open " << options.outFile << " for writing." << std::endl;
    }

    // remove temp files
    for (unsigned i = 0; i < length(chunkFiles); ++i)
        unlink(toCString(chunkFiles[i]));

    if (!success)
        return 1;

    double stop = sysTime();
    if (options.verbose)
    {
        std::cerr << "Number of alignments:\t" << options.numRecords << std::endl;
        std::cerr << "Elapsed total time:  \t" << stop - start << " seconds" << std::endl;
    }
    return 0;
}
