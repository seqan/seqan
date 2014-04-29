// ==========================================================================
//                              sam2Matrix
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
    seqan::ArgumentParser parser("sam2Matrix");
    // Set short description, version, and date.
    setShortDescription(parser, "This program determines for each read in the reference file if it has an entry in"
                                " the provided sam files stating that it mapped.");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "\\fB-sf\\fP \\fISAMFILE\\fP \\fB-rf\\fP \\fIREADSFILE\\fP \\fB-gf\\fP \\fIGENOMEFILE\\fP"
                         "\\fB-o\\fP \\fIOUTFILE\\fP");
    addDescription(parser, "This program determines for each read in the reference file if it has an entry in the "
                           "provided sam files stating that it mapped. Afterwards a file is generated containing a row"
                           " for each read which contains the read ID and the index of the mapped references.");

    addOption(parser, ArgParseOption("sf", "samFiles", "Sam files.", ArgParseOption::INPUTFILE, "SAMFILES", true));
    setValidValues(parser, "sf", "sam");
    setRequired(parser, "sf");
    addOption(parser, ArgParseOption("rf", "referenceFile", "Read name file.", ArgParseOption::INPUTFILE));
    setValidValues(parser, "referenceFile", "fa fasta fq fastq");
    setRequired(parser, "rf");

    addOption(parser, ArgParseOption("gf", "genomeFiles", "Names of the references of the corresponding sam file..", ArgParseOption::INPUTFILE, "GENOMES", true));
    setValidValues(parser, "genomeFiles", "fa fasta fq fastq");
    setRequired(parser, "gf");

    addOption(parser, ArgParseOption("o", "outputFile", "Output file.", ArgParseOption::OUTPUTFILE));
    setRequired(parser, "o");
    setValidValues(parser, "outputFile", "csv");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBsam2Gasic\\fP \\fB-sf\\fP \\fIsam.sam\\fP \\fB-rf\\fP \\fIreads.fa\\fP \\fB-gf\\fP "
                        "\\fIecoli.fa\\fP \\fB-o\\fP \\fIout.csv\\fP",
                        "Call with \\fISAMFILE\\fP set to \"sam.sam\", with  \\fIREADSFILE\\fP set to \"reads.fa\", with "
                        "\\fIGENOMEFILE\\fP set to \"ecoli.fa\" and \\fIOUTFILE\\fP set to \"out.csv\", where \"ecoli.fa\""
                        " is the file used as reference to generate \"sam.sam\" ");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    options.samFileNames = getOptionValues(parser, "samFiles");
    getOptionValue(options.readNameFileName, parser, "referenceFile");
    options.genomeFileNames = getOptionValues(parser, "genomeFiles");
    getOptionValue(options.outPutFileName, parser, "outputFile");

    return seqan::ArgumentParser::PARSE_OK;
}

// A std::map is used to quickly find a read id when parsing the sam
// files.
bool _initializeMap(std::map<String<char>, unsigned> & nameToPos, SamToGasicOptions const & options)
{
    seqan::CharString id;
    seqan::Dna5String seq;

    SequenceStream seqStream(toCString(options.readNameFileName));
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    seqan::CharString fixedId;
    for (unsigned i = 0 ;!atEnd(seqStream); ++i)
    {
        if (readRecord(id, seq, seqStream) != 0)
        {
            std::cerr << "ERROR: Could not read from example.fa!\n";
            return 1;
        }

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
        BamStream bamIO(toCString(options.samFileNames[i])); 
        BamAlignmentRecord record;

        for (unsigned j = 0 ;!atEnd(bamIO); ++j)
        {
            if (readRecord(record, bamIO) != 0)
                return false;
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

bool _writeFile(StringSet<String<unsigned> > const & result,
                std::map<String<char>, unsigned> & nameToPos,
                SamToGasicOptions const & options)
{
    std::fstream fout(toCString(options.outPutFileName), std::ios::binary | std::ios::out);

    if (streamPut(fout, length(options.samFileNames)))    
        return 1;
    if (streamWriteChar(fout, '\t'))
        return 1;
    if (streamPut(fout, length(result)))
        return 1;
    if (streamWriteChar(fout, '\n'))
        return 1;

    for (unsigned i = 0; i < length(options.genomeFileNames); ++i)
    {
        if (streamWriteChar(fout, '>'))
            return 1;
        if (streamWriteBlock(fout, &(options.genomeFileNames[i])[0], length(options.genomeFileNames[i])) != length(options.genomeFileNames[i]))
            return 1;
        if (streamWriteChar(fout, '\n'))
            return 1;
    }

    std::map<String<char>, unsigned>::iterator it = nameToPos.begin();
    for (unsigned i = 0; i < length(result); ++i)
    {
        if (streamWriteBlock(fout, &(it->first)[0], length(it->first)) != length(it->first))
            return 1;
        if (streamWriteChar(fout, '\t'))
            return 1;

        for (unsigned j = 0; j < length(result[i]); ++j)
        {
            if (streamPut(fout, result[i][j]))
                return 1;
            if (streamWriteChar(fout, '\t'))
                return 1;
        }
        if (streamWriteChar(fout, '\n'))
            return 1;

        ++it;
    }

    fout.close();

    return 0;
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

    // A std::map is used to quickly find a read id when parsing the sam
    if (!_initializeMap(nameToPos, options))
        std::cerr << "Problem extracting the read names." << std::endl;
    if (!_parseSamFiles(mappedReads, nameToPos, options))
        std::cerr << "Problem parsing the sam files." << std::endl;
    if (!_writeFile(mappedReads, nameToPos, options))
        std::cerr << "Problem writing the result." << std::endl;

    return 0;
}

