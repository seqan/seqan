// ==========================================================================
//                                   samcat
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

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/misc/misc_name_store_cache.h>

#include <seqan/arg_parse.h>

using namespace seqan;

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
    std::vector<std::string> inFiles;
    std::string outFile;
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function mergeBamFiles()
// --------------------------------------------------------------------------

template <typename TWriter, typename TFormat>
int mergeBamFiles(TWriter &writer, Tag<TFormat> writeFormat, std::vector<std::string> &inFiles)
{
    typedef StringSet<CharString>       TNameStore;
    typedef NameStoreCache<TNameStore>  TNameStoreCache;

    TNameStore contigNameStore;
    TNameStoreCache contigNameStoreCache(contigNameStore);

    BamIOContext<TNameStore> bamIOContext(contigNameStore, contigNameStoreCache);
    int returnedError = 0;
    int error;


    // Step 1: Merge all headers (if available)
    BamHeader header;
    for (unsigned i = 0; i < length(inFiles); ++i)
    {
        if (guessFormatFromFilename(inFiles[i], Bam()))
        {
            Stream<Bgzf> reader;
            if (!open(reader, inFiles[i].c_str(), "r"))
            {
                std::cout << "ERROR: Couldn't open " << inFiles[i] << " for reading." << std::endl;
                inFiles[i].clear();
                returnedError = -1;
                continue;
            }

            error = readRecord(header, bamIOContext, reader, Bam());
        }
        else
        {
            // Construct a RecordReader from the input file.
            std::ifstream inStream(inFiles[i].c_str());
            if (!inStream.is_open())
            {
                std::cout << "ERROR: Couldn't open " << inFiles[i] << " for reading." << std::endl;
                clear(inFiles[i]);
                returnedError = -1;
                continue;
            }

            RecordReader<std::ifstream, SinglePass<> > reader(inStream);
            if (atEnd(reader))
                continue;

            error = readRecord(header, bamIOContext, reader, Sam());
        }

        if (error != 0)
        {
            std::cerr << "ERROR: Problem reading header from file " << inFiles[i] << std::endl;
            returnedError = -1;
        }
    }


    // Step 2: Remove duplicate header entries and write merged header
    removeDuplicates(header);
    write2(writer, header, bamIOContext, writeFormat);


    // Step 3: Read
    BamAlignmentRecord record;
    for (unsigned i = 0; i != length(inFiles); ++i)
    {
        // avoid to open failed files twice
        if (inFiles[i].empty())
            continue;

        if (guessFormatFromFilename(inFiles[i], Bam()))
        {
            Stream<Bgzf> reader;
            if (!open(reader, inFiles[i].c_str(), "r"))
            {
                std::cout << "ERROR: Couldn't open " << inFiles[i] << " for reading." << std::endl;
                returnedError = -1;
                continue;
            }

            // skip header
            BamHeader tmp;
            readRecord(tmp, bamIOContext, reader, Bam());

            while (!atEnd(reader))
            {
                error = readRecord(record, bamIOContext, reader, Bam());
                if (error)
                {
                    std::cerr << "ERROR: Problem reading record from file " << inFiles[i] << std::endl;
                    returnedError = -1;
                    break;
                }
                write2(writer, record, bamIOContext, writeFormat);
            }
        }
        else
        {
            // Construct a RecordReader from the input file.
            std::ifstream inStream(inFiles[i].c_str());
            if (!inStream.is_open())
            {
                std::cout << "ERROR: Couldn't open " << inFiles[i] << " for reading." << std::endl;
                returnedError = -1;
                continue;
            }

            RecordReader<std::ifstream, SinglePass<> > reader(inStream);
            if (atEnd(reader))
                continue;

            // skip header
            BamHeader tmp;
            readRecord(tmp, bamIOContext, reader, Sam());

            while (!atEnd(reader))
            {
                error = readRecord(record, bamIOContext, reader, Sam());
                if (error)
                {
                    std::cerr << "ERROR: Problem reading record from file " << inFiles[i] << std::endl;
                    returnedError = -1;
                    break;
                }
                write2(writer, record, bamIOContext, writeFormat);
            }
        }
    }
    return returnedError;
}


// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("samcat");
    // Set short description, version, and date.
    setShortDescription(parser, "SAM/BAM file concatenation and conversion");
    setCategory(parser, "Utilities");
    setVersion(parser, "0.1");
    setDate(parser, "May 2014");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIINFILE\\fP> [<\\fIINFILE\\fP> ...] [-o <\\fIOUTFILE\\fP>]");
    addDescription(parser, "This tool reads a set of input files in SAM or BAM format and outputs the concatenation of them. "
                           "If the output file name is ommitted the result is written to standard output in SAM format.");

    addDescription(parser, "(c) Copyright 2014 by David Weese.");

    // We require one argument.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, "INFILE", true));
    setValidValues(parser, 0, ".sam .bam");
    setHelpText(parser, 0, "Input SAM or BAM file.");

    addOption(parser, ArgParseOption("o", "output", "Output file name", ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "output", ".sam .bam");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBsamcat\\fP mapped1.sam mapped2.sam \\fI-o merged.sam\\fP",
                "Merge two SAM files.");
    addListItem(parser, "\\fBsamcat\\fP input.sam \\fI-o ouput.bam\\fP",
                "Convert a SAM file into BAM format.");

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    options.inFiles = getArgumentValues(parser, 0);
    getOptionValue(options.outFile, parser, "output");

    return ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    if (guessFormatFromFilename(options.outFile, Bam()))
    {
        Stream<Bgzf> writer;
        if (!open(writer, options.outFile.c_str(), "w"))
        {
            std::cout << "ERROR: Couldn't open " << options.outFile << " for writing." << std::endl;
            return -1;
        }
        return mergeBamFiles(writer, Bam(), options.inFiles);
    }
    else if (!options.outFile.empty())
    {
        std::ofstream writer(options.outFile.c_str());
        if (!writer.is_open())
        {
            std::cout << "ERROR: Couldn't open " << options.outFile << " for writing." << std::endl;
            return -1;
        }
        return mergeBamFiles(writer, Sam(), options.inFiles);
    }
    else
    {
        // dump to standard out stream
        return mergeBamFiles(std::cout, Sam(), options.inFiles);
    }
}
