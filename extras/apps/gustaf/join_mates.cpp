// ==========================================================================
//                                 join_mates
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#include <sstream>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>

#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

// using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.

struct JoinMatesOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Path to Fasta/Fastq input files
    seqan::CharString inPath1;
    seqan::CharString inPath2;

    // Path to Fasta output file
    seqan::CharString outPath;
    // CharString outPathPos;

    // Whether or not to rev-compl the second input file
    bool revCompl;

    JoinMatesOptions() :
        verbosity(1),
        revCompl(true)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(JoinMatesOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("join_mates");
    // Set short description, version, and date.
    setShortDescription(parser, "Joining paired-end files.");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIMATES1 FASTA/FASTQ FILE\\fP\" \"\\fIMATES2 FASTA/FASTQ FILE\\fP\"");
    addDescription(parser, "Joining two paired-end files into one file with joined (single-end) reads. Automatically reverse-complements reads of the second input file.");

    // We require two arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "FASTA/FASTQ FILE 1"));
    setValidValues(parser, 0, "fa fasta fq fastq");
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "FASTA/FASTQ FILE 2"));
    setValidValues(parser, 1, "fasta fa fastq fq");

    addOption(parser, seqan::ArgParseOption("o", "outPath", "Set name of output FASTA file.", seqan::ArgParseOption::OUTPUTFILE, "FASTA"));
    setValidValues(parser, "o", "fasta fa");
    addOption(parser, seqan::ArgParseOption("rc", "revcompl", "Disable reverse complementing second input file."));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    /*
    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBjoin_mates\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");
    */

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.inPath1, parser, 0);
    getArgumentValue(options.inPath2, parser, 1);

    getOptionValue(options.outPath, parser, "o");
    options.revCompl = !isSet(parser, "rc");

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function _getShortId()
// ----------------------------------------------------------------------------

// Creates a short Id out of a long one (i.e. it takes the prefix til the first white space)
template <typename TId>
void _getShortId(TId & shortId, TId const & longId)
{
    clear(shortId);
    for (typename seqan::Position<TId>::Type i = 0; i < length(longId) && isgraph(value(longId, i)); ++i)
    {
        appendValue(shortId, value(longId, i));
    }
}

// ----------------------------------------------------------------------------
// Function checkUniqueId()
// ----------------------------------------------------------------------------

// Checks whether the short ID (sId) from a long ID (id) is uniq regarding all short IDs of the given
// read set (sQueryIds). Prints IDs to cerr if not.
template <typename TId>
bool _checkUniqueId(TId const & sId, TId const & id, seqan::StringSet<TId> & ids, seqan::StringSet<TId> & sQueryIds)
{
    bool unique = true;
    for (unsigned j = 0; j < length(sQueryIds); ++j)
    {
        if (sId == sQueryIds[j])
        {
            std::cerr << "Found nonunique sequence ID!" << std::endl;
            std::cerr << ids[j] << std::endl;
            std::cerr << id << std::endl;
            std::cerr << "###########################" << std::endl;
            unique = false;
        }
    }
    return unique;
}

// --------------------------------------------------------------------------
// Function _importSequences()
// --------------------------------------------------------------------------

// Imports mate pairs from two files, joins them, stores the joining position in
// String readJoinPositions, and stores the sequences in the StringSet seqs and
// their identifiers in the StringSet ids
template <typename TSequence, typename TId>
inline bool
_importSequences(seqan::CharString const & fileNameL,
                 seqan::CharString const & fileNameR,
                 seqan::StringSet<TSequence> & seqs,
                 seqan::StringSet<TId> & ids,
                 seqan::StringSet<TId> & sIds,
                 seqan::String<unsigned> & readJoinPositions
                 )
{
    seqan::MultiSeqFile leftMates;
    seqan::MultiSeqFile rightMates;

    if (!open(leftMates.concat, toCString(fileNameL), seqan::OPEN_RDONLY))
    {
        std::cerr << "Failed to open " << fileNameL << " file." << std::endl;
        return false;
    }
    if (!open(rightMates.concat, toCString(fileNameR), seqan::OPEN_RDONLY))
    {
        std::cerr << "Failed to open " << fileNameR << " file." << std::endl;
        return false;
    }

    seqan::AutoSeqFormat formatL;
    guessFormat(leftMates.concat, formatL);
    split(leftMates, formatL);

    seqan::AutoSeqFormat formatR;
    guessFormat(rightMates.concat, formatR);
    split(rightMates, formatR);

    unsigned seqCount = length(leftMates);

    resize(readJoinPositions, seqCount);
    reserve(seqs, seqCount, seqan::Exact());
    reserve(ids, seqCount, seqan::Exact());
    reserve(sIds, seqCount, seqan::Exact());

    TSequence seq;
    TSequence seqL;
    TSequence seqR;
    TId id;
    TId sId;
    unsigned counter = 0;
    for (unsigned i = 0; i < seqCount; ++i)
    {
        assignSeq(seqL, leftMates[i], formatL);
        assignSeq(seqR, rightMates[i], formatR);
        reverseComplement(seqR);
        readJoinPositions[i] = length(seq);
        append(seq, seqL);
        append(seq, seqR);
        assignSeqId(id, leftMates[i], formatL);
        appendValue(seqs, seq, seqan::Generous());
        appendValue(ids, id, seqan::Generous());

        _getShortId(sId, id);
        if (!_checkUniqueId(sId, id, ids, sIds))
            ++counter;
        appendValue(sIds, sId);
        clear(seq);
    }

    std::cout << "Loaded " << seqCount << " mate pair sequence" << ((seqCount > 1) ? "s." : ".") << std::endl;
    if (counter > 0)
        std::cout << "Found " << counter << " nonunique sequence IDs" << std::endl;

    return true;
}

// --------------------------------------------------------------------------
// Function _writeSequences()
// --------------------------------------------------------------------------

// Writes out sequences and ids in FASTA format


template <typename TSequence>
int _writeSequences(seqan::CharString & outPath,
                seqan::StringSet<TSequence> const & seqs,
                seqan::StringSet<seqan::CharString> const & sIds)
{
    std::ofstream f(toCString(outPath));
    // seqan::SequenceStream outStream(toCString(outPath), seqan::SequenceStream::WRITE);
    if (!f.good())
    {
        std::cerr << "Error: Could not open output file!" << std::endl;
        return 1;
    }
    for (unsigned i = 0; i < length(seqs); ++i)
    {
        if (writeRecord(f, sIds[i], seqs[i], seqan::Fasta()) != 0)
        {
            std::cerr << "Error: Could not write to file!" << std::endl;
            return 1;
        }

    }
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
    JoinMatesOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << '\n'
                  << "INPUT FILE 1     \t" << options.inPath1 << '\n'
                  << "INPUT FILE 2     \t" << options.inPath2 << '\n'
                  << "OUTPUT FILE     \t" << options.outPath << '\n'
                  << "OPTION REV COMPL     \t" << options.revCompl << "\n\n";
    }

    typedef seqan::String<seqan::Dna5> TSequence;
    seqan::StringSet<TSequence> seqs;
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> sIds;
    seqan::String<unsigned> joinPos;
    // Read in paired-end reads
    _importSequences(options.inPath1, options.inPath2, seqs, ids, sIds, joinPos);

    // Write out one file with joined sequences in FASTA format
    _writeSequences(options.outPath, seqs, sIds);


    return 0;
}
