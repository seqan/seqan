// ==========================================================================
//                                 gustaf_mate_joining
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

// This struct stores the options from the command line.
struct JoinMatesOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Path to Fasta/Fastq input files
    seqan::String<seqan::CharString> inPaths;

    // Path to Fasta output file
    seqan::String<seqan::CharString> outPaths;
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
    seqan::ArgumentParser parser("gustaf_mate_joining");
    // Set short description, version, and date.
    setShortDescription(parser, "Joining paired-end files.");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIMATES1 FASTA/FASTQ FILE\\fP\" \"\\fIMATES2 FASTA/FASTQ FILE\\fP\"");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fIMATES FASTA/FASTQ FILE\\fP\" \"\\fI-o MATES1 FASTA/FASTQ FILE\\fP\" \"\\fI-o MATES2 FASTA/FASTQ FILE\\fP\"");
    addDescription(parser, "Joining two paired-end files into one file with joined (single-end) reads. Automatically reverse-complements reads of the second input file.");
    addDescription(parser, "This simple program takes as input two mate pair or paired-end files and outputs a file "
                            "where both mate sequences have been joined together. The FASTA file with joined mates "
                            "is an required input file for the paired-end mode of Gustaf. "
                            "The tool assumes the mates in the second file to be reverse complemented compared to the "
                            "first file. This behaviour can be turned off using the command line argument \"-rc\".");

    addDescription(parser, "Given only one input file and two output files, the program will split the reads from "
                            "the input files at half length, and write the first half of each sequence as mates1 into "
                            "the first output file and the reversed complemented second half of each sequence as "
                            "mates2 into the second output file. Reverse complementing the sequences can again be "
                            "turned off using \"-rc\".");

    addDescription(parser, "To prepare the joined mate file for the Gustaf paired-end example, call \n "

    "./gustaf_mate_joining adeno_modified_reads_mates1.fa adeno_modified_reads_mates2.fa "
        "-rc -o adeno_modified_reads_joinedMates.fa");

    // We require two arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "FASTA/FASTQ FILE(S)", true));
    setValidValues(parser, 0, "fa fasta fq fastq");
    /*
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, "FASTA/FASTQ FILE 2"));
    setValidValues(parser, 1, "fasta fa fastq fq");
    */

    addOption(parser, seqan::ArgParseOption("o", "outPath", "Set name of output FASTA/FASTQ file(s).", seqan::ArgParseOption::OUTPUTFILE, "FASTA/FASTQ", true));
    setValidValues(parser, "o", "fasta fa fq fastq");
    addOption(parser, seqan::ArgParseOption("rc", "revcompl", "Disable reverse complementing second input file."));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    resize(options.inPaths, getArgumentValueCount(parser, 0), seqan::Exact());
    for (unsigned i = 0; i < length(options.inPaths); ++i)
        getArgumentValue(options.inPaths[i], parser, 0, i);
    /*
    getArgumentValue(options.inPath1, parser, 0);
    getArgumentValue(options.inPath2, parser, 1);
    */

    resize(options.outPaths, getOptionValueCount(parser, "o"), seqan::Exact());
    for (unsigned i = 0; i < length(options.outPaths); ++i)
        getOptionValue(options.outPaths[i], parser, "o", i);
    //getOptionValue(options.outPath, parser, "o");
    if (isSet(parser, "rc"))
        options.revCompl = false;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    if (length(options.inPaths) == 1 && length(options.outPaths) == 1)
        return seqan::ArgumentParser::PARSE_ERROR;
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
                 bool revCompl,
                 seqan::StringSet<TSequence> & seqs,
                 seqan::StringSet<TId> & ids,
                 seqan::StringSet<TId> & sIds,
                 seqan::StringSet<seqan::CharString> & quals,
                 seqan::String<unsigned> & readJoinPositions
                 )
{
    seqan::SequenceStream l(toCString(fileNameL));
    seqan::SequenceStream r(toCString(fileNameR));
    if (!isGood(l) || !isGood(r))
    {
        std::cerr << "Failed to open file." << std::endl;
        return false;
    }

    TSequence seq;
    TSequence seqL;
    TSequence seqR;
    TId id;
    TId sId;
    seqan::CharString qual;
    seqan::CharString qualL;
    seqan::CharString qualR;
    unsigned counter = 0;
    while (!atEnd(l) || !atEnd(r))
    {
        if (readRecord(id, seqL, qualL, l) != 0)
        {
            std::cerr << "Problem reading from first input file." << std::endl;
            return false;
        }
        if (readRecord(id, seqR, qualR, r) != 0)
        {
            std::cerr << "Problem reading from first input file." << std::endl;
            return false;
        }

        appendValue(readJoinPositions, length(seqL));
        if (revCompl)
        {
            reverseComplement(seqR);
            reverse(qualR);
        }
        append(seq, seqL);
        append(seq, seqR);
        append(qual, qualL);
        append(qual, qualR);
        appendValue(seqs, seq, seqan::Generous());
        appendValue(quals, qual, seqan::Generous());
        appendValue(ids, id, seqan::Generous());

        _getShortId(sId, id);
        if (!_checkUniqueId(sId, id, ids, sIds))
            ++counter;
        appendValue(sIds, sId);
        clear(seq);
        clear(qual);
    }
    return true;
}

// Imports mate pairs from one file, separates them, and stores the sequences in the StringSets seqs and and mateSeqs
// and their identifiers in the StringSet ids
// Note: Assumes equally long mates, i.e. splits in the middle of the read
template <typename TSequence, typename TId>
inline bool
_importSequences(seqan::CharString const & fileName,
                 seqan::StringSet<TSequence> & seqs,
                 seqan::StringSet<TSequence> & mateSeqs,
                 seqan::StringSet<TId> & ids,
                 seqan::StringSet<TId> & sIds,
                 seqan::StringSet<seqan::CharString> & quals,
                 seqan::StringSet<seqan::CharString> & mateQuals
                 )
{
    typedef typename seqan::Position<TSequence>::Type TPos;
    seqan::SequenceStream f(toCString(fileName));
    if (!isGood(f))
    {
        std::cerr << "Failed to open file." << std::endl;
        return false;
    }

    TSequence seq;
    TSequence seqL;
    TSequence seqR;
    TId id;
    TId sId;
    seqan::CharString qual;
    seqan::CharString qualL;
    seqan::CharString qualR;
    unsigned counter = 0;
    TPos splitPos;
    while (!atEnd(f))
    {
        if (readRecord(id, seq, qual, f) != 0)
        {
            std::cerr << "Problem reading from first input file." << std::endl;
            return false;
        }

        splitPos = static_cast<TPos>(length(seq)/2);
        append(seqL, prefix(seq, splitPos));
        append(seqR, suffix(seq, splitPos));
        appendValue(seqs, seqL, seqan::Generous());
        appendValue(mateSeqs, seqR, seqan::Generous());
        append(qualL, prefix(qual, splitPos));
        append(qualR, suffix(qual, splitPos));
        appendValue(quals, qualL, seqan::Generous());
        appendValue(mateQuals, qualR, seqan::Generous());
        appendValue(ids, id, seqan::Generous());


        _getShortId(sId, id);
        if (!_checkUniqueId(sId, id, ids, sIds))
            ++counter;
        appendValue(sIds, sId);
        clear(seqL);
        clear(seqR);
        clear(qualL);
        clear(qualR);
    }
    return true;
}

// --------------------------------------------------------------------------
// Function _writeSequences()
// --------------------------------------------------------------------------

// Writes out sequences and ids in FASTA format
template <typename TSequence>
int _writeSequences(seqan::CharString & outPath,
                seqan::StringSet<TSequence> const & seqs,
                seqan::StringSet<seqan::CharString> const & sIds,
                seqan::StringSet<seqan::CharString> const & quals)
{
    seqan::SequenceStream seqStream(toCString(outPath), seqan::SequenceStream::WRITE);
    if (!isGood(seqStream))
    {
        std::cerr << "Error: Could not open output file!" << std::endl;
        return 1;
    }
    if (writeAll(seqStream, sIds, seqs, quals) != 0)
    {
        std::cerr << "Error: Could not write to file!" << std::endl;
        return 1;
    }
    return 0;
}

// Writes out sequences and ids in FASTA format
template <typename TSequence>
int _writeSequences(seqan::CharString & outPath1,
                seqan::CharString & outPath2,
                seqan::StringSet<TSequence> const & seqs,
                seqan::StringSet<TSequence> const & mateSeqs,
                seqan::StringSet<seqan::CharString> const & sIds,
                seqan::StringSet<seqan::CharString> const & quals,
                seqan::StringSet<seqan::CharString> const & mateQuals
                )
{
    seqan::SequenceStream f1(toCString(outPath1), seqan::SequenceStream::WRITE);
    seqan::SequenceStream f2(toCString(outPath2), seqan::SequenceStream::WRITE);
    // seqan::SequenceStream outStream(toCString(outPath), seqan::SequenceStream::WRITE);
    if (!isGood(f1) || !isGood(f2))
    {
        std::cerr << "Error: Could not open output file!" << std::endl;
        return 1;
    }
    for (unsigned i = 0; i < length(seqs); ++i)
    {
        if (writeRecord(f1, sIds[i], seqs[i], quals[i]) != 0 || writeRecord(f2, sIds[i], mateSeqs[i], mateQuals[i]) != 0)
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
                  << "INPUT FILE 1     \t" << options.inPaths[0] << '\n';
                  if (length(options.inPaths) > 1)
                      std::cout << "INPUT FILE 2     \t" << options.inPaths[1] << '\n';
                  std::cout << "OUTPUT FILE     \t" << options.outPaths[0] << '\n';
                  if (length(options.outPaths) > 1)
                      std::cout << "OUTPUT FILE 2     \t" << options.outPaths[1] << '\n';
                  std::cout << "OPTION REV COMPL     \t" << options.revCompl << "\n\n";
    }

    typedef seqan::String<seqan::Dna5> TSequence;
    seqan::StringSet<TSequence> seqs;
    seqan::StringSet<TSequence> mateSeqs;
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> sIds;
    seqan::StringSet<seqan::CharString> quals;
    seqan::StringSet<seqan::CharString> mateQuals;
    seqan::String<unsigned> joinPos;
    if (length(options.inPaths) > 1)
    {
        // Read in paired-end reads
        _importSequences(options.inPaths[0], options.inPaths[1], options.revCompl, seqs, ids, sIds, quals, joinPos);
        // Write out one file with joined sequences in FASTA format
        _writeSequences(options.outPaths[0], seqs, sIds, quals);
    } else
    {
        // Read in joined reads and output two FASTA files
        _importSequences(options.inPaths[0], seqs, mateSeqs, ids, sIds, quals, mateQuals);
        // Write out one file with joined sequences in FASTA format
        _writeSequences(options.outPaths[0], options.outPaths[1], seqs, mateSeqs, sIds, quals, mateQuals);
    }

    return 0;
}
