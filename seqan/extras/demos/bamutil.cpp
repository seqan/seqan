// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// A simple converter between SAM and BAM.  BAM files can be dumped region-
// wise, using a .bai index.
// ==========================================================================

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/bam_io.h>
#include <seqan/stream.h>

#if SEQAN_HAS_ZLIB

using namespace seqan;

enum Format
{
    FORMAT_SAM,
    FORMAT_BAM
};

struct Options
{
    bool showHelp;
    bool showVersion;
    CharString inFile;
    CharString baiFile;
    CharString outFile;
    Format inFormat;
    Format outFormat;
    unsigned verbosity;
    String<Triple<CharString, int, int> > regions;  // <(refName, beginPos, endPos)>

    Options()
    {
        showHelp = false;
        showVersion = false;
        inFile = "-";
        outFile = "-";
        inFormat = FORMAT_SAM;
        outFormat = FORMAT_SAM;
        verbosity = 1;
    }
};

void
setupCommandLineParser(CommandLineParser & parser, Options const & options)
{
    addVersionLine(parser, "1.0");
    
    addTitleLine(parser, "***********");
    addTitleLine(parser, "* BAMUTIL *");
    addTitleLine(parser, "***********");
    addTitleLine(parser, "");
    addTitleLine(parser, "SeqAn SAM/BAM Utility.");
    addTitleLine(parser, "");
    addTitleLine(parser, "Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>");

    addUsageLine(parser, "bamutil [OPTIONS] <IN >OUT");
    
	addSection(parser, "General Options");
    addOption(parser, CommandLineOption("v", "verbose", "Enable verbose mode.", OptionType::Bool));
    addOption(parser, CommandLineOption("vv", "very-verbose", "Enable very verbose mode.", OptionType::Bool));

	addSection(parser, "Input Specification");
    addOption(parser, CommandLineOption("i", "input-file", "Path to input, '-' for stdin.", OptionType::String, options.inFile));
    addOption(parser, CommandLineOption("S", "input-sam", "Input file is SAM.", OptionType::Bool, options.inFormat == FORMAT_SAM));
    addOption(parser, CommandLineOption("B", "input-bam", "Input file is BAM.", OptionType::Bool, options.inFormat == FORMAT_BAM));
    addOption(parser, CommandLineOption("bi", "bai-index-file", "Path to BAI index, default: input + '.bai'.", OptionType::String));

	addSection(parser, "Range Specification");
    addOption(parser, CommandLineOption("r", "region", "Regions to dump (in which aligments start). REF:FROM-TO, e.g. IV:1,000-2,000 will alignments with ref IV in range [1000,2000) (zero based).", OptionType::String | OptionType::List, options.inFile));

	addSection(parser, "Output Specification");
    addOption(parser, CommandLineOption("o", "output-file", "Path to output, '-' for stdout.", OptionType::String, options.outFile));
    addOption(parser, CommandLineOption("s", "output-sam", "Output file is SAM.", OptionType::Bool, options.outFormat == FORMAT_SAM));
    addOption(parser, CommandLineOption("b", "output-bam", "Output file is BAM.", OptionType::Bool, options.outFormat == FORMAT_BAM));
}

// Parse region string CHR:FROM-TO, e.g. "X:10,000-100,000" and write the result to a triple.  Return true on success,
// false on failure.

bool
parseRegion(Triple<CharString, int, int> & result, CharString regionString)
{
    Stream<CharArray<char const *> > stream(&regionString[0], &regionString[0] + length(regionString));
    RecordReader<Stream<CharArray<char const *> >, SinglePass<> > reader(stream);

    clear(result.i1);
    int res = readUntilChar(result.i1, reader, ':');
    if (res != 0)
        return false;
    if (goNext(reader))
        return false;
    CharString buf1;
    res = readUntilChar(buf1, reader, '-');
    if (res != 0)
        return false;
    if (goNext(reader))
        return false;
    CharString castBuf;
    for (unsigned i = 0; i < length(buf1); ++i)
        if (isdigit(buf1[i]))
            appendValue(castBuf, buf1[i]);
    if (!lexicalCast2(result.i2, castBuf))
        return false;

    clear(buf1);
    res = readLine(buf1, reader);
    if (res != 0 && res != EOF_BEFORE_SUCCESS)
        return false;
    clear(castBuf);
    for (unsigned i = 0; i < length(buf1); ++i)
        if (isdigit(buf1[i]))
            appendValue(castBuf, buf1[i]);
    if (!lexicalCast2(result.i3, castBuf))
        return false;
    return true;
}

// Parse the command line and check for any syntatical errors.

int parseCommandLineAndCheck(Options & options,
                             CommandLineParser & parser,
                             int argc,
                             char const ** argv)
{
    bool stop = !parse(parser, argc, argv);
    if (stop)
        return 1;
    if (isSetLong(parser, "help"))
    {
        options.showHelp = true;
        return 0;
    }
    if (isSetLong(parser, "version"))
    {
        options.showVersion = true;
        return 0;
    }

    if (isSetLong(parser, "verbose"))
        options.verbosity = 2;
    if (isSetLong(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValueLong(parser, "input-file", options.inFile);
    if (isSetLong(parser, "input-sam"))
        options.inFormat = FORMAT_SAM;
    if (isSetLong(parser, "input-bam"))
        options.inFormat = FORMAT_BAM;
    getOptionValueLong(parser, "bai-index-file", options.baiFile);
    if (length(options.baiFile) == 0 && options.inFormat == FORMAT_BAM && options.inFile != "-")
    {
        options.baiFile = options.inFile;
        append(options.baiFile, ".bai");
    }

    getOptionValueLong(parser, "output-file", options.outFile);
    if (isSetLong(parser, "output-sam"))
        options.outFormat = FORMAT_SAM;
    if (isSetLong(parser, "output-bam"))
        options.outFormat = FORMAT_BAM;
    if (isSetLong(parser, "region"))
    {
        String<CharString> regions = getOptionValuesLong(parser, "region");
        for (unsigned i = 0; i < length(regions); ++i)
        {
            Triple<CharString, int, int> t;
            if (parseRegion(t, regions[i]))
            {
                appendValue(options.regions, t);
                if (options.verbosity >= 2)
                    std::cerr << "[VERBOSE] Region " << t.i1 << ":" << t.i2 << "-" << t.i3 << std::endl;
            }
            else
            {
                std::cerr << "[WARNING] could not parse region \"" << regions[i] << "\". IGNORING." << std::endl;
            }
        }
    }

	return 0;
}

template <typename TInStreamOrRecordReader, typename TOutStream, typename TOutTag, typename TContext, typename TOptions>
int _dumpRegion(TInStreamOrRecordReader &, TOutStream &, Sam const & /*sam*/, TOutTag const &, TContext &, TOptions const &)
{
    std::cerr << "Dumping of regions not supported in SAM!" << std::endl;
    return 0;
}

template <typename TInStreamOrRecordReader, typename TOutStream, typename TOutTag, typename TContext, typename TOptions>
int _dumpRegion(TInStreamOrRecordReader & in, TOutStream & out, Bam const & /*bam*/, TOutTag const & outTag, TContext & context, TOptions const & options)
{
    // TOOD(holtgrew): The index is loaded for each region.  This should probably not be the case!

    // Dump each region after loading the index.
    BamIndex<Bai> bamIndex;
    if (read(bamIndex, toCString(options.baiFile)) != 0)
    {
        std::cerr << "[ERROR] Could not open index file " << options.baiFile << ", required when specifying regions." << std::endl;
        return 1;
    }

    for (unsigned i = 0; i < length(options.regions); ++i)
    {
        // Jump near range.
        CharString refName = options.regions[i].i1;
        unsigned refId = 0;
        if (!getIdByName(nameStore(context), refName, refId, nameStoreCache(context)))
        {
            std::cerr << "[ERROR] Unknown reference " << refName << std::endl;
            return 1;
        }
        bool hasAlignments = false;
        if (!jumpToPos(in, hasAlignments, context, refId, options.regions[i].i2, bamIndex))
        {
            std::cerr << "[ERROR] Could not jump to " << refName << ":" << options.regions[i].i2 << std::endl;
            return 1;
        }
        if (!hasAlignments)
            continue;
        // Dump range.
        BamAlignmentRecord record;
        bool stop = false;
        while (!atEnd(in) && !stop)
        {
            if (readRecord(record, context, in, Bam()) != 0)
            {
                std::cerr << "Could not read alignment record!" << std::endl;
                return 1;
            }
            if (record.pos < options.regions[i].i2)
                continue;  // Skip, before region.
            if (record.pos >= options.regions[i].i3)
            {
                stop = true;
                continue;  // Stop, wrote out all required ones.
            }
            if (write2(out, record, context, outTag) != 0)
            {
                std::cerr << "Could not write alignment record!" << std::endl;
                return 1;
            }
        }
    }
    return 0;
}

template <typename TInStreamOrRecordReader, typename TOutStream, typename TInTag, typename TOutTag, typename TOptions>
int performConversion(TInStreamOrRecordReader & in, TOutStream & out, TOptions const & options, TInTag const & inTag, TOutTag const & outTag)
{
    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache;

    TNameStore nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    BamIOContext<TNameStore> context(nameStore, nameStoreCache);

    BamHeader header;
    if (readRecord(header, context, in, inTag) != 0)
    {
        std::cerr << "Could not read header!" << std::endl;
        return 1;
    }
    if (write2(out, header, context, outTag) != 0)
    {
        std::cerr << "Could not write header!" << std::endl;
        return 1;
    }

    if (length(options.regions) == 0u)
    {
        // Simply dump whole file.
        BamAlignmentRecord record;
        while (!atEnd(in))
        {
            if (readRecord(record, context, in, inTag) != 0)
            {
                std::cerr << "Could not read alignment record!" << std::endl;
                return 1;
            }
            if (write2(out, record, context, outTag) != 0)
            {
                std::cerr << "Could not write alignment record!" << std::endl;
                return 1;
            }
        }
    }
    else
    {
        return _dumpRegion(in, out, inTag, outTag, context, options);
    }
    return 0;
}
#endif  // #if SEQAN_HAS_ZLIB

// The main function.
//
// Don't be intimidated by its longness, most of the code is for branching the different options to types.

int main(int argc, char const * argv[])
{
    (void)argc;
    (void)argv;
#if SEQAN_HAS_ZLIB
    using namespace seqan;

    // -----------------------------------------------------------------------
    // Handle Command Line
    // -----------------------------------------------------------------------

    // Setup command line parser.
    CommandLineParser parser;
    Options options;
    setupCommandLineParser(parser, options);
    
    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    int res = parseCommandLineAndCheck(options, parser, argc, argv);
    if (res != 0)
        return 1;
    if (options.showHelp || options.showVersion)
        return 0;

    // -----------------------------------------------------------------------
    // Open Input / Output Files.
    // -----------------------------------------------------------------------

    // We go through quite some pain here since we want to support both stdin/stdout and file writing for both SAM and
    // BAM.  The only way to do this for BGZF/BAM is to use the raw bgzf API.

    // We will use 
    std::istream * inS = &std::cin;
    std::ostream * outS = &std::cout;
    int inF = 0; //STDIN_FILENO;
    int outF = 1; //STDOUT_FILENO;

    bool ioGood = true;
    if (options.inFormat == FORMAT_SAM)
    {
        if (options.inFile != "-")
        {
            inS = new std::fstream(toCString(options.inFile), std::ios::binary | std::ios::in);
            if (!inS->good())
            {
                std::cerr << "Could not open file " << options.inFile << std::endl;
                ioGood = false;
            }
        }
    }
    else
    {
        if (options.inFile != "-")
        {
            inF = open(toCString(options.inFile), O_RDONLY);
            if (inF == -1)
            {
                std::cerr << "Could not open file " << options.inFile << std::endl;
                ioGood = false;
            }
        }
    }
    if (options.outFormat == FORMAT_SAM)
    {
        if (options.outFile != "-")
        {
            outS = new std::fstream(toCString(options.outFile), std::ios::binary | std::ios::out);
            if (!outS->good())
            {
                std::cerr << "Could not open file " << options.inFile << std::endl;
                ioGood = false;
            }
        }
    }
    else
    {
        if (options.outFile != "-")
        {
#ifdef PLATFORM_WINDOWS
            outF = open(toCString(options.outFile), O_CREAT | O_WRONLY, 00666);
#else  // #ifdef PLATFORM_WINDOWS
            outF = open(toCString(options.outFile), O_CREAT | O_WRONLY | O_DIRECT, 00666);
#endif   // #ifdef PLATFORM_WINDOWS
            if (outF == -1)
            {
                std::cerr << "Could not open file " << options.outFile << std::endl;
                ioGood = false;
            }
        }
    }

    if (!ioGood)
        goto main_end;

    // -----------------------------------------------------------------------
    // 4-way branch to SAM/BAM formats to read and write.
    // -----------------------------------------------------------------------

    if (options.inFormat == FORMAT_SAM)
    {
        RecordReader<std::istream, SinglePass<> > reader(*inS);
        if (options.outFormat == FORMAT_SAM)  // SAM -> SAM
        {
            res = performConversion(reader, *outS, options, Sam(), Sam());
        }
        else  // SAM -> BAM
        {
            Stream<Bgzf> bamOutStream;
            attachToFile(bamOutStream, outF, OPEN_CREATE | OPEN_WRONLY);
            res = performConversion(reader, bamOutStream, options, Sam(), Bam());
        }
    }
    else
    {
        if (options.outFormat == FORMAT_SAM)  // BAM -> SAM
        {
            Stream<Bgzf> bamInStream;
            attachToFile(bamInStream, inF, OPEN_RDONLY);
            res = performConversion(bamInStream, *outS, options, Bam(), Sam());
        }
        else  // BAM -> BAM
        {
            Stream<Bgzf> bamInStream;
            Stream<Bgzf> bamOutStream;
            attachToFile(bamInStream, inF, OPEN_RDONLY);
            attachToFile(bamOutStream, outF, OPEN_CREATE | OPEN_WRONLY);
            res = performConversion(bamInStream, bamOutStream, options, Bam(), Bam());
        }
    }
    if (res != 0)
    {
        std::cerr << "Error during conversion!" << std::endl;
        return 1;
    }

    // -----------------------------------------------------------------------
    // Cleanup.
    // -----------------------------------------------------------------------

main_end:

    if (options.outFile != "-" && options.outFormat == FORMAT_SAM)
        delete outS;
    if (options.outFile != "-" && options.outFormat != FORMAT_SAM)
        close(outF);
    if (options.inFile != "-" && options.inFormat == FORMAT_SAM)
        delete inS;
    if (options.inFile != "-" && options.inFormat != FORMAT_SAM)
        close(inF);
    
    return ioGood ? 0 : 1;
#else  // #if SEQAN_HAS_ZLIB
    return 0;
#endif  // #if SEQAN_HAS_ZLIB
}
