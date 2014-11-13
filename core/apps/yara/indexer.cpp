// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the yara_indexer application.
// ==========================================================================

#define YARA_INDEXER

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "store_seqs.h"
#include "misc_timer.h"
#include "misc_tags.h"
#include "misc_types.h"
#include "misc_options.h"
#include "bits_matches.h"
#include "index_fm.h"

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString contigsFile;
    CharString contigsIndexFile;

    bool    verbose;

    Options() :
        verbose(false)
    {}
};

// ----------------------------------------------------------------------------
// Class Indexer
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec = void>
struct Indexer
{
    typedef SeqStore<void, YaraContigsConfig>       TContigs;
    typedef Index<YaraContigsFM, TIndexSpec>        TIndex;

    TContigs            contigs;
    SeqFileIn           contigsFile;
    TIndex              index;
    Timer<double>       timer;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & /* options */)
{
    setAppName(parser, "yara_indexer");
    setShortDescription(parser, "Yara Indexer");
    setCategory(parser, "Read Mapping");

    setDateAndVersion(parser);
    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] <\\fIREFERENCE FILE\\fP>");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "REFERENCE FILE"));
    setValidValues(parser, 0, SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "A reference genome file.");

    addOption(parser, ArgParseOption("v", "verbose", "Displays verbose output."));

    addSection(parser, "Output Options");

    addOption(parser, ArgParseOption("o", "output-prefix", "Specify a filename prefix for the reference genome index. \
                                     Default: use the filename prefix of the reference genome.", ArgParseOption::OUTPUT_PREFIX));

    addOption(parser, ArgParseOption("td", "tmp-dir", "Specify a temporary directory where to construct the index. \
                                     Default: use the output directory.", ArgParseOption::STRING));
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Parse verbose output option.
    getOptionValue(options.verbose, parser, "verbose");

    // Parse contigs input file.
    getArgumentValue(options.contigsFile, parser, 0);

    // Parse contigs index prefix.
    getOptionValue(options.contigsIndexFile, parser, "output-prefix");
    if (!isSet(parser, "output-prefix"))
        options.contigsIndexFile = trimExtension(options.contigsFile);

    // Parse and set temp dir.
    CharString tmpDir;
    getOptionValue(tmpDir, parser, "tmp-dir");
    if (!isSet(parser, "tmp-dir"))
        tmpDir = getPath(options.contigsIndexFile);
    setEnv("TMPDIR", tmpDir);

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function loadContigs()
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec>
void loadContigs(Indexer<TIndexSpec, TSpec> & me, Options const & options)
{
    if (options.verbose)
        std::cerr << "Loading reference:\t\t\t" << std::flush;

    start(me.timer);
    open(me.contigsFile, toCString(options.contigsFile));

    // Compute file size.
//    me.contigsFile.seekg(0, std::ios::end);
//    unsigned fileSize = me.contigsFile.tellg();
//    me.contigsFile.seekg(0, std::ios::beg);

    // Reserve space for contigs.
//    reserve(me.contigs, fileSize);

    try
    {
        readRecords(me.contigs, me.contigsFile, Iupac());
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }
    stop(me.timer);

    if (length(me.contigs.seqs) > MemberLimits<Match<TSpec>, ContigId>::VALUE)
        throw RuntimeError("Maximum number of contigs exceeded.");

    if (maxLength(me.contigs.seqs) > MemberLimits<Match<TSpec>, ContigSize>::VALUE)
        throw RuntimeError("Maximum contig length exceeded.");

    if (options.verbose)
        std::cerr << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function saveContigs()
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec>
void saveContigs(Indexer<TIndexSpec, TSpec> & me, Options const & options)
{
    if (options.verbose)
        std::cerr << "Dumping reference:\t\t\t" << std::flush;

    start(me.timer);
    if (!save(me.contigs, toCString(options.contigsIndexFile)))
        throw RuntimeError("Error while dumping reference file.");
    stop(me.timer);

    if (options.verbose)
        std::cerr << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function buildIndex()
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec>
void buildIndex(Indexer<TIndexSpec, TSpec> & me, Options const & options)
{
    typedef typename Indexer<TIndexSpec, TSpec>::TIndex TIndex;

    if (options.verbose)
        std::cerr << "Building reference index:\t\t" << std::flush;

    start(me.timer);

    try
    {
        // Randomize Ns in contigs.
        randomizeNs(me.contigs);

        // IndexFM is built on the reversed contigs.
        reverse(me.contigs);

        // Set the index text.
        // NOTE(esiragusa): this assignment implicitly assigns and converts the contigs to the index contigs.
        setValue(me.index.text, me.contigs.seqs);

        // Clears the contigs.
        // NOTE(esiragusa): the index now owns its own contigs.
        shrinkToFit(me.contigs.seqs);

        // Iterator instantiation calls automatic index construction.
        typename Iterator<TIndex, TopDown<> >::Type it(me.index);
        ignoreUnusedVariableWarning(it);
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to index the reference.");
    }
    catch (IOError const & /* e */)
//    catch (PageFrameError const & /* e */)
    {
        throw RuntimeError("Insufficient disk space to index the reference. \
                            Specify a bigger temporary folder using the options --tmp-dir.");
    }

    // NOTE(esiragusa): the contigs should be already saved.
//    reverse(me.contigs);
    stop(me.timer);

    if (options.verbose)
        std::cerr << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function saveIndex()
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec>
void saveIndex(Indexer<TIndexSpec, TSpec> & me, Options const & options)
{
    if (options.verbose)
        std::cerr << "Dumping reference index:\t\t" << std::flush;

    start(me.timer);
    if (!save(me.index, toCString(options.contigsIndexFile)))
        throw RuntimeError("Error while dumping reference index file.");
    stop(me.timer);

    if (options.verbose)
        std::cerr << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function runIndexer()
// ----------------------------------------------------------------------------

template <typename TIndexSpec, typename TSpec>
void runIndexer(Indexer<TIndexSpec, TSpec> & me, Options const & options)
{
    loadContigs(me, options);
    saveContigs(me, options);
    buildIndex(me, options);
    saveIndex(me, options);
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    try
    {
        Indexer<YaraIndexSpec, YaraStringSpec> indexer;
        runIndexer(indexer, options);
    }
    catch (BadAlloc const & /* e */)
    {
        std::cerr << "Insufficient memory." << std::endl;
        return 1;
    }
    catch (Exception const & e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
