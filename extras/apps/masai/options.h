// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the common Option class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_OPTIONS_H_
#define SEQAN_EXTRAS_MASAI_OPTIONS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class MasaiOptions
// ----------------------------------------------------------------------------

struct MasaiOptions
{
    typedef std::string             TString;
    typedef std::vector<TString>    TList;

    enum IndexType
    {
        INDEX_ESA, INDEX_SA, INDEX_QGRAM, INDEX_FM
    };

    enum MappingMode
    {
        ALL, ALL_BEST, ANY_BEST
    };

    enum OutputFormat
    {
        RAW, SAM, SAM_NO_CIGAR
    };

    TList       indexTypeList;
    TList       mappingModeList;
    TList       outputFormatList;
    TList       outputFormatExtensions;

    MasaiOptions()
    {
        indexTypeList.push_back("esa");
        indexTypeList.push_back("sa");
        indexTypeList.push_back("qgram");
        indexTypeList.push_back("fm");

        mappingModeList.push_back("all");
        mappingModeList.push_back("all-best");
        mappingModeList.push_back("any-best");

        outputFormatList.push_back("raw");
        outputFormatList.push_back("sam");
        outputFormatList.push_back("sam-no-cigar");

        outputFormatExtensions.push_back("raw");
        outputFormatExtensions.push_back("sam");
        outputFormatExtensions.push_back("sam");
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setEnv()
// ----------------------------------------------------------------------------

template <typename TString, typename TValue>
bool setEnv(TString & key, TValue & value)
{
#ifdef PLATFORM_WINDOWS
    CharString env(key);
    appendValue(env, '=');
    append(env, value);
    return !putenv(toCString(env));
#else
    return !setenv(toCString(key), toCString(value), true);
#endif
}

// ----------------------------------------------------------------------------
// Function lastOf()
// ----------------------------------------------------------------------------

template <typename TString, typename TToken>
typename Iterator<TString, Standard>::Type
lastOf(TString const & string, TToken const & token)
{
    typedef typename Iterator<TString, Standard>::Type    TIterator;

    TIterator it = end(string, Standard()) - length(token);

    for (TIterator itBegin = begin(string, Standard());
         it != itBegin && !isEqual(infix(string, it, it + length(token)), token);
         goPrevious(it)) ;

    return it;
}

// ----------------------------------------------------------------------------
// Function trimExtension()
// ----------------------------------------------------------------------------

template <typename TString>
Segment<TString, PrefixSegment>
trimExtension(TString & string)
{
    return prefix(string, lastOf(string, '.'));
}

// ----------------------------------------------------------------------------
// Function getPath()
// ----------------------------------------------------------------------------

template <typename TString>
Segment<TString, PrefixSegment>
getPath(TString & string)
{
#ifdef PLATFORM_WINDOWS
    return prefix(string, lastOf(string, '\\'));
#else
    return prefix(string, lastOf(string, '/'));
#endif
}

// ----------------------------------------------------------------------------
// Function getFilename()
// ----------------------------------------------------------------------------

template <typename TString>
Segment<TString, SuffixSegment>
getFilename(TString & string)
{
#ifdef PLATFORM_WINDOWS
    return suffix(string, lastOf(string, '\\'));
#else
    return suffix(string, lastOf(string, '/'));
#endif
}

// ----------------------------------------------------------------------------
// Function getOptionValue()
// ----------------------------------------------------------------------------

template <typename TOption, typename TString, typename TOptionsList>
void getOptionValue(TOption & option,
                    ArgumentParser const & parser,
                    TString const & optionName,
                    TOptionsList & optionsList)
{
    typedef typename Iterator<TOptionsList, Standard>::Type TOptionsIterator;
    typedef typename Value<TOptionsList>::Type              TOptionString;

    TOptionsIterator optionsBegin = begin(optionsList, Standard());
    TOptionsIterator optionsEnd = end(optionsList, Standard());

    TOptionString optionStr;
    getOptionValue(optionStr, parser, optionName);

    TOptionsIterator optionType = std::find(optionsBegin, optionsEnd, optionStr);

    SEQAN_ASSERT(optionType != optionsEnd);
    option = TOption(optionType - optionsBegin);
}

// ----------------------------------------------------------------------------
// Function setDateAndVersion()
// ----------------------------------------------------------------------------

void setDateAndVersion(ArgumentParser & parser)
{
    std::string rev  = "$Revision$";
    std::string date = "$Date$";

    setCategory(parser, "Read Mapping");
    setVersion(parser, "0.7 [" + rev.substr(11, rev.size() - 13) + "]");
    setDate(parser, date.substr(7, std::min((int)date.size() - 8, 10)));
}

// ----------------------------------------------------------------------------
// Function setDescription()
// ----------------------------------------------------------------------------

void setDescription(ArgumentParser & parser)
{
    addDescription(parser, "Masai is a fast and accurate read mapper based on approximate seeds and multiple backtracking.");
    addDescription(parser, "See \\fIhttp://www.seqan.de/projects/masai\\fP for more information.");
    addDescription(parser, "(c) Copyright 2011-2012 by Enrico Siragusa.");
}

// ----------------------------------------------------------------------------
// Function setIndexType()
// ----------------------------------------------------------------------------

template <typename TOptions>
void setIndexType(ArgumentParser & parser, TOptions const & options)
{
    addOption(parser, ArgParseOption("x", "index", "Select the genome index type.", ArgParseOption::STRING));
    setValidValues(parser, "index", options.indexTypeList);
    setDefaultValue(parser, "index", options.indexTypeList[options.genomeIndexType]);
}

// ----------------------------------------------------------------------------
// Function getIndexType()
// ----------------------------------------------------------------------------

template <typename TOptions>
void getIndexType(TOptions & options, ArgumentParser const & parser)
{
    getOptionValue(options.genomeIndexType, parser, "index", options.indexTypeList);
}

// ----------------------------------------------------------------------------
// Function setIndexPrefix()
// ----------------------------------------------------------------------------

void setIndexPrefix(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("xp", "index-prefix", "Specify an genome index prefix name. \
                                     Default: use the genome filename prefix.", ArgParseOption::STRING));
}

// ----------------------------------------------------------------------------
// Function getIndexPrefix()
// ----------------------------------------------------------------------------

template <typename TOptions>
void getIndexPrefix(TOptions & options, ArgumentParser const & parser)
{
    getOptionValue(options.genomeIndexFile, parser, "index-prefix");
    if (!isSet(parser, "index-prefix"))
        options.genomeIndexFile = trimExtension(options.genomeFile);
}

// ----------------------------------------------------------------------------
// Function setOutputFormat()
// ----------------------------------------------------------------------------

template <typename TOptions>
void setOutputFormat(ArgumentParser & parser, TOptions const & options)
{
    addOption(parser, ArgParseOption("of", "output-format", "Select the output format.", ArgParseOption::STRING));
    setValidValues(parser, "output-format", options.outputFormatList);
    setDefaultValue(parser, "output-format", options.outputFormatList[options.outputFormat]);
}

// ----------------------------------------------------------------------------
// Function getOutputFormat()
// ----------------------------------------------------------------------------

template <typename TOptions>
void getOutputFormat(TOptions & options, ArgumentParser const & parser)
{
    getOptionValue(options.outputFormat, parser, "output-format", options.outputFormatList);
}

// ----------------------------------------------------------------------------
// Function setOutputFile()
// ----------------------------------------------------------------------------

void setOutputFile(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("o", "output-file", "Specify an output file. \
                                     Default: use the reads filename and guess the extension.",
                                     ArgParseOption::OUTPUTFILE));
//    setValidValues(parser, "output-file", "raw sam");
}

// ----------------------------------------------------------------------------
// Function getOutputFile()
// ----------------------------------------------------------------------------

template <typename TString, typename TOptions, typename TSuffix>
void getOutputFile(TString & file,
                   TOptions const & options,
                   ArgumentParser const & parser,
                   TString const & from,
                   TSuffix const & suffix)
{
    getOptionValue(file, parser, "output-file");
    if (!isSet(parser, "output-file"))
    {
        file = trimExtension(from);
        append(file, suffix);
        appendValue(file, '.');
        append(file, options.outputFormatExtensions[options.outputFormat]);
    }
}

// ----------------------------------------------------------------------------
// Function setTmpFolder()
// ----------------------------------------------------------------------------

void setTmpFolder(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("t", "tmp-folder", "Specify a huge temporary folder. \
                                     Default: use the genome folder.", ArgParseOption::STRING));
}

// ----------------------------------------------------------------------------
// Function getTmpFolder()
// ----------------------------------------------------------------------------

template <typename TOptions>
void getTmpFolder(TOptions const & options, ArgumentParser const & parser)
{
    CharString tmpFolder;
    getOptionValue(tmpFolder, parser, "tmp-folder");
    if (!isSet(parser, "tmp-folder"))
        tmpFolder = getPath(options.genomeFile);
    setEnv("TMPDIR", tmpFolder);
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_OPTIONS_H_
