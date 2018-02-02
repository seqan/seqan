// ==========================================================================
//                              sam2matrix
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
// This program determines for each read in the reference file if it has an 
// entry in the provided sam files stating that it mapped.
//
// ==========================================================================
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
#include <map>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class SamToGasicOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
using namespace seqan;

struct SamToGasicOptions
{
    String<String<char> > samFileNames;
    String<String<char> > genomeFileNames;
    String<char> readNameFileName;
    String<char> outPutFileName;
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(SamToGasicOptions& options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("sam2matrix");
    // Set short description, version, and date.
    setShortDescription(parser, "This program outputs for each read the ids of references it maps to.");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setCategory(parser, "Metagenomics");

    // Define usage line and long description.
    addUsageLine(parser, "\\fB-m\\fP \\fIa.sam\\fP \\fB-m\\fP \\fIb.sam\\fP \\fB-r\\fP \\fIreads\\fP \\fB-rf\\fP "
                         "\\fIref_a.fasta\\fP \\fB-rf\\fP \\fIref_b.fasta\\fP \\fB-o\\fP \\fIout.tsv\\fP");
    addDescription(parser, "This program determines for each read in the reference file if it has an entry in the "
                           "provided sam files stating that it mapped. Afterwards a file is generated containing a row"
                           " for each read which contains the read ID and the index of the mapped references.");

    addOption(parser, ArgParseOption("m", "mapping", "File containing the mappings.", ArgParseOption::INPUT_FILE, 
                                     "FILE", true));
    setValidValues(parser, "mapping", BamFileIn::getFileExtensions());
    setRequired(parser, "mapping");
    addOption(parser, ArgParseOption("r", "reads", "File containing the reads contained in the mapping file(s).", 
                                     ArgParseOption::INPUT_FILE, "FILE"));
    setValidValues(parser, "reads", SeqFileIn::getFileExtensions());
    setRequired(parser, "reads");

    addOption(parser, ArgParseOption("rf", "reference", "Name of the file used as reference of the corresponding sam "
                                     "file. If not specified the names of the mapping files are taken", 
                                     ArgParseOption::STRING, "STRING", true));

    addOption(parser, ArgParseOption("o", "out", "Output file.", ArgParseOption::OUTPUT_FILE));
    setRequired(parser, "o");
    setValidValues(parser, "o", ".tsv");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fB./sam2matrix\\fP \\fB-m\\fP \\fIa.sam\\fP \\fB-m\\fP \\fIb.sam\\fP \\fB-r\\fP "
                        "\\fIreads.fasta\\fP \\fB-rf\\fP \\fIref_a.fasta\\fP \\fB-rf\\fP \\fIref_b.fasta\\fP "
                        "\\fB-o\\fP \\fIout.tsv\\fP",
                        "Storing in \\fIout.tsv\\fP for each read contained in \\fIreads.fasta\\fP if it mapped to "
                        "the references in \\fIref_a.fasta\\fP or \\fIref_b.fasta\\fP using the corresponding sam "
                        "files \\fIa.sam\\fP and \\fIb.sam\\fP.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    options.samFileNames = getOptionValues(parser, "m");
    getOptionValue(options.readNameFileName, parser, "r");
    options.genomeFileNames = getOptionValues(parser, "rf");
    getOptionValue(options.outPutFileName, parser, "o");

    if (length(options.samFileNames) > length(options.genomeFileNames))
        for (unsigned i = length(options.genomeFileNames); i < length(options.samFileNames); ++i)
            appendValue(options.genomeFileNames, options.samFileNames[i]);

    return seqan::ArgumentParser::PARSE_OK;
}

// A std::map is used to quickly find a read id when parsing the sam
// files.
bool _initializeMap(std::map<String<char>, unsigned> & nameToPos, SamToGasicOptions const & options)
{
    seqan::CharString id;
    seqan::Dna5String seq;

    SeqFileIn seqStream(toCString(options.readNameFileName));

    seqan::CharString fixedId;
    for (unsigned i = 0; !atEnd(seqStream); ++i)
    {
        readRecord(id, seq, seqStream);

        clear(fixedId);
        for(unsigned s = 0; s < length(id); ++s)
        {
            if(id[s] == ' ') break;
            appendValue(fixedId, id[s]);
        }

        nameToPos[fixedId] = i;
    }

    return true;
}

bool _parseSamFiles(StringSet<String<unsigned> > & mappedReads, 
                    std::map<String<char>, unsigned> & nameToPos, 
                    SamToGasicOptions const & options)
{
    resize(mappedReads, nameToPos.size());

    for (unsigned i = 0; i < length(options.samFileNames); ++i)
    {
        BamFileIn bamIO(toCString(options.samFileNames[i]));
        BamHeader header;
        BamAlignmentRecord record;

        readHeader(header, bamIO);

        for (unsigned j = 0; !atEnd(bamIO); ++j)
        {
            readRecord(record, bamIO);
            if (!hasFlagUnmapped(record))
            {
                if (nameToPos.find(record.qName) != nameToPos.end())
                {
                    if (length(mappedReads[nameToPos[record.qName]]) == 0 || back(mappedReads[nameToPos[record.qName]]) != i)
                        appendValue(mappedReads[nameToPos[record.qName]], i);
                }
            }
        }
    }
    return true;
}

void _writeFile(StringSet<String<unsigned> > const & result,
                std::map<String<char>, unsigned> & nameToPos,
                SamToGasicOptions const & options)
{
    std::fstream fout(toCString(options.outPutFileName), std::ios::binary | std::ios::out);
    Iter<std::fstream, StreamIterator<Output> > target(fout);

    appendNumber(target, length(options.samFileNames));
    writeValue(target, '\t');
    appendNumber(target, length(result));
    writeValue(target, '\n');

    for (unsigned i = 0; i < length(options.genomeFileNames); ++i)
    {
        writeValue(target, '>');
        write(target, options.genomeFileNames[i]);
        writeValue(target, '\n');
    }

    std::map<String<char>, unsigned>::iterator it = nameToPos.begin();
    for (unsigned i = 0; i < length(result); ++i)
    {
        write(target, it->first);
        writeValue(target, '\t');

        for (unsigned j = 0; j < length(result[i]); ++j)
        {
            appendNumber(target, result[i][j]);
            writeValue(target, '\t');
        }
        writeValue(target, '\n');
        ++it;
    }
}



// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    SamToGasicOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    StringSet<String<unsigned> > mappedReads;
    std::map<String<char>, unsigned> nameToPos;

    try
    {
        // A std::map is used to quickly find a read id when parsing the sam
        if (!_initializeMap(nameToPos, options))
            std::cerr << "Problem extracting the read names." << std::endl;
        if (!_parseSamFiles(mappedReads, nameToPos, options))
            std::cerr << "Problem parsing the sam files." << std::endl;
        _writeFile(mappedReads, nameToPos, options);
    }
    catch (IOError &e)
    {
        std::cerr << "Problem occurred: " << e.what() << std::endl;
    }

    return 0;
}

