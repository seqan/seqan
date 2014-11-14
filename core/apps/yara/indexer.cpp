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
#include "bits_matches.h"
#include "misc_options.h"
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
    CharString      contigsFile;
    CharString      contigsIndexFile;

    __uint64        maxContigLength;
    __uint64        maxContigSetLength;
    __uint64        maxContigSetLengthSum;

    bool            verbose;

    Options() :
        maxContigLength(),
        maxContigSetLength(),
        maxContigSetLengthSum(),
        verbose(false)
    {}
};

// ----------------------------------------------------------------------------
// Class YaraIndexer
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = void>
struct YaraIndexer
{
    typedef SeqStore<TSpec, YaraContigsConfig<> >   TContigs;

    Options const &     options;
    TContigs            contigs;
    SeqFileIn           contigsFile;
    Timer<double>       timer;

    YaraIndexer(Options const & options) :
        options(options)
    {}
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

    addOption(parser, ArgParseOption("", "max-contigs-length", "Enforce the maximum length of any contig.", ArgParseOption::INT64));
    setMinValue(parser, "max-contigs-length", "0");
    hideOption(getOption(parser, "max-contigs-length"));

    addOption(parser, ArgParseOption("", "max-contigs-count", "Enforce the maximum number of contigs.", ArgParseOption::INT64));
    setMinValue(parser, "max-contigs-count", "0");
    hideOption(getOption(parser, "max-contigs-count"));

    addOption(parser, ArgParseOption("", "max-contigs-lengthsum", "Enforce the maximum length of all contigs.", ArgParseOption::INT64));
    setMinValue(parser, "max-contigs-lengthsum", "0");
    hideOption(getOption(parser, "max-contigs-lengthsum"));

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

    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Parse verbose output option.
    getOptionValue(options.verbose, parser, "verbose");

    // Parse contigs input file.
    getArgumentValue(options.contigsFile, parser, 0);

    // Parse contigs index prefix.
    getOptionValue(options.contigsIndexFile, parser, "output-prefix");
    if (!isSet(parser, "output-prefix"))
        options.contigsIndexFile = trimExtension(options.contigsFile);

    // Parse hidden limits.
    getOptionValue(options.maxContigLength, parser, "max-contigs-length");
    getOptionValue(options.maxContigSetLength, parser, "max-contigs-count");
    getOptionValue(options.maxContigSetLengthSum, parser, "max-contigs-lengthsum");
    options.maxContigSetLengthSum = std::max(options.maxContigSetLengthSum, options.maxContigLength);

    // Parse and set temp dir.
    CharString tmpDir;
    getOptionValue(tmpDir, parser, "tmp-dir");
    if (!isSet(parser, "tmp-dir"))
        tmpDir = getPath(options.contigsIndexFile);
    setEnv("TMPDIR", tmpDir);

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function loadContigs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void loadContigs(YaraIndexer<TSpec, TConfig> & me)
{
    if (me.options.verbose)
        std::cerr << "Loading reference:\t\t\t" << std::flush;

    start(me.timer);

    if (!open(me.contigsFile, toCString(me.options.contigsFile)))
        throw RuntimeError("Error while opening the reference file.");

    try
    {
        readRecords(me.contigs, me.contigsFile, Iupac());
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }

    stop(me.timer);

    if (me.options.verbose)
        std::cerr << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function saveContigs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void saveContigs(YaraIndexer<TSpec, TConfig> & me)
{
    if (me.options.verbose)
        std::cerr << "Saving reference:\t\t\t" << std::flush;

    start(me.timer);
    if (!saveContigsLimits(me.options) || !save(me.contigs, toCString(me.options.contigsIndexFile)))
        throw RuntimeError("Error while saving the reference.");
    stop(me.timer);

    if (me.options.verbose)
        std::cerr << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function saveIndex()
// ----------------------------------------------------------------------------

template <typename TContigsLen, typename TContigsSize, typename TContigsSum, typename TSpec, typename TConfig>
void saveIndex(YaraIndexer<TSpec, TConfig> & me)
{
    typedef typename YaraFMIndexContigs<TContigsLen, TContigsSize, TContigsSum>::Type   TIndexContigs;
    typedef YaraFMIndexConfig<TContigsLen, TContigsSize, TContigsSum>                   TIndexConfig;
    typedef FMIndex<void, TIndexConfig>                                                 TIndexSpec;
    typedef Index<TIndexContigs, TIndexSpec>                                            TIndex;

    if (me.options.verbose)
        std::cerr << "Building reference index:\t\t" << std::flush;

    start(me.timer);

    // Randomly replace Ns with A, C, G, T.
    randomizeNs(me.contigs);

    // IndexFM is built on the reversed contigs.
    reverse(me.contigs);

    TIndex index;

    // This assignment *copies* the contigs to the index as the types differ.
    setValue(index.text, me.contigs.seqs);

    // Clear the contigs - the index now owns its own copy.
    clear(me.contigs);
    shrinkToFit(me.contigs);

    try
    {
        // Iterator instantiation triggers index construction.
        typename Iterator<TIndex, TopDown<> >::Type it(index);
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

    stop(me.timer);

    if (me.options.verbose)
        std::cerr << me.timer << std::endl;

    if (me.options.verbose)
        std::cerr << "Saving reference index:\t\t\t" << std::flush;

    start(me.timer);
    if (!save(index, toCString(me.options.contigsIndexFile)))
        throw RuntimeError("Error while saving the reference index file.");
    stop(me.timer);

    if (me.options.verbose)
        std::cerr << me.timer << std::endl;
}

template <typename TContigsLen, typename TSpec, typename TConfig>
void saveIndex(YaraIndexer<TSpec, TConfig> & me)
{
    if (me.options.maxContigSetLengthSum < MaxValue<__uint32>::VALUE)
    {
        if (me.options.maxContigSetLength < MaxValue<__uint8>::VALUE)
        {
            saveIndex<__uint8, TContigsLen, __uint32>(me);
        }
        else if (me.options.maxContigSetLength < MaxValue<__uint16>::VALUE)
        {
            saveIndex<__uint16, TContigsLen, __uint32>(me);
        }
        else
        {
            throw RuntimeError("Maximum number of contigs exceeded.");
        }
    }
#ifdef YARA_LARGE_CONTIGS
    else
    {
        if (me.options.maxContigSetLength < MaxValue<__uint8>::VALUE)
        {
            saveIndex<__uint8, TContigsLen, __uint64>(me);
        }
        else if (me.options.maxContigSetLength < MaxValue<__uint16>::VALUE)
        {
            saveIndex<__uint16, TContigsLen, __uint64>(me);
        }
        else
        {
            throw RuntimeError("Maximum number of contigs exceeded.");
        }
    }
#else
    else
    {
        throw RuntimeError("Maximum contigs lengthsum exceeded.");
    }
#endif
}

template <typename TSpec, typename TConfig>
void saveIndex(YaraIndexer<TSpec, TConfig> & me)
{
    if (me.options.maxContigLength < MaxValue<__uint32>::VALUE)
    {
        saveIndex<__uint32>(me);
    }
    else
    {
#ifdef YARA_LARGE_CONTIGS
        saveIndex<__uint64>(me);
#else
        throw RuntimeError("Maximum contig length exceeded.");
#endif
    }
}

// ----------------------------------------------------------------------------
// Function runYaraIndexer()
// ----------------------------------------------------------------------------

void runYaraIndexer(Options & options)
{
    YaraIndexer<> indexer(options);

    loadContigs(indexer);
    saveContigs(indexer);
    setContigsLimits(options, indexer.contigs.seqs);
    saveIndex(indexer);
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

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    try
    {
        runYaraIndexer(options);
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
