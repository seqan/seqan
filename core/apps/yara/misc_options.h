// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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

#ifndef APP_YARA_MISC_OPTIONS_H_
#define APP_YARA_MISC_OPTIONS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Functor EqualsChar
// ----------------------------------------------------------------------------
// TODO(esiragusa): remove this when new tokenization gets into develop.

template <char VALUE>
struct EqualsChar
{
    template <typename TValue>
    bool operator() (TValue const & val) const
    {
        return val == VALUE;
    }
};

// ----------------------------------------------------------------------------
// Composite Functors
// ----------------------------------------------------------------------------
// TODO(esiragusa): remove this when new tokenization gets into develop.

typedef EqualsChar<'.'>        IsDot;
typedef EqualsChar<' '>        IsSpace;
typedef EqualsChar<'/'>        IsSlash;
typedef EqualsChar<'\\'>       IsBackSlash;

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

template <typename TString, typename TFunctor>
inline typename Iterator<TString const, Standard>::Type
lastOf(TString const & string, TFunctor const & f)
{
    typedef typename Iterator<TString const, Standard>::Type TIterator;

    TIterator itBegin = begin(string, Standard());
    TIterator itEnd = end(string, Standard());

    for (TIterator it = itEnd; it != itBegin; )
        if (f(value(--it)))
            return it;

    return itEnd;
}

// ----------------------------------------------------------------------------
// Function trimExtension()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Prefix<TString>::Type
trimExtension(TString & string)
{
    return prefix(string, lastOf(string, IsDot()));
}

// ----------------------------------------------------------------------------
// Function getExtension()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Suffix<TString>::Type
getExtension(TString & string)
{
    return suffix(string, lastOf(string, IsDot()) + 1);
}

// ----------------------------------------------------------------------------
// Function getPath()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Prefix<TString const>::Type
getPath(TString const & string)
{
#ifdef PLATFORM_WINDOWS
    return prefix(string, lastOf(string, IsBackSlash()));
#else
    return prefix(string, lastOf(string, IsSlash()));
#endif
}

// ----------------------------------------------------------------------------
// Function getFilename()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Suffix<TString>::Type
getFilename(TString & string)
{
#ifdef PLATFORM_WINDOWS
    return suffix(string, lastOf(string, IsBackSlash()));
#else
    return suffix(string, lastOf(string, IsSlash()));
#endif
}

// ----------------------------------------------------------------------------
// Function getOptionEnum()
// ----------------------------------------------------------------------------

template <typename TOption, typename TString, typename TOptionsList>
void getOptionEnum(TOption & option,
                   TString const & optionStr,
                   TOptionsList & optionsList)
{
    typedef typename Iterator<TOptionsList, Standard>::Type TOptionsIterator;

    TOptionsIterator optionsBegin = begin(optionsList, Standard());
    TOptionsIterator optionsEnd = end(optionsList, Standard());

    TOptionsIterator optionPos = std::find(optionsBegin, optionsEnd, optionStr);

    option = (optionPos != optionsEnd) ? TOption(optionPos - optionsBegin) : TOption();
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
    typedef typename Value<TOptionsList>::Type              TOptionString;

    TOptionString optionStr;
    getOptionValue(optionStr, parser, optionName);

    return getOptionEnum(option, optionStr, optionsList);
}

// ----------------------------------------------------------------------------
// Function setDateAndVersion()
// ----------------------------------------------------------------------------

void setDateAndVersion(ArgumentParser & parser)
{
    setCategory(parser, "Read Mapping");

#ifdef SEQAN_REVISION
    setVersion(parser, "0.8.0 [" + std::string(SEQAN_REVISION) + "]");
#else
    setVersion(parser, "0.8.0");
#endif
#ifdef SEQAN_DATE
    setDate(parser, SEQAN_DATE);
#endif
}

// ----------------------------------------------------------------------------
// Function setDescription()
// ----------------------------------------------------------------------------

void setDescription(ArgumentParser & parser)
{
    addDescription(parser, "Yara - Yet Another Read Aligner.");
    addDescription(parser, "See \\fIhttp://www.seqan.de/projects/yara\\fP for more information.");
    addDescription(parser, "(c) Copyright 2011-2014 by Enrico Siragusa <enrico.siragusa@fu-berlin.de>.");
    addDescription(parser, "(c) Copyright 2013 by NVIDIA Corporation.");
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
    addOption(parser, ArgParseOption("xp", "index-prefix", "Specify a filename prefix for the reference genome index. \
                                     Default: use the filename prefix of the reference genome.", ArgParseOption::STRING));
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
// Function getInputType()
// ----------------------------------------------------------------------------

template <typename TOptions, typename TString>
void getInputType(TOptions & options, TString const & inputFile)
{
    TString inputExtension = getExtension(inputFile);
    typename TOptions::TString inputTypeExtension(toCString(inputExtension));
    getOptionEnum(options.inputType, inputTypeExtension, options.inputTypeList);
}

// ----------------------------------------------------------------------------
// Function getOutputFormat()
// ----------------------------------------------------------------------------

template <typename TOptions, typename TString>
void getOutputFormat(TOptions & options, TString const & outputFile)
{
    TString outputExtension = getExtension(outputFile);
    typename TOptions::TString outputFormatExtension(toCString(outputExtension));
    getOptionEnum(options.outputFormat, outputFormatExtension, options.outputFormatList);
}

// ----------------------------------------------------------------------------
// Function setOutputFile()
// ----------------------------------------------------------------------------

template <typename TOptions>
void setOutputFile(ArgumentParser & parser, TOptions const & options)
{
    addOption(parser, ArgParseOption("o", "output-file", "Specify an output file. \
                                     Default: use the reads filename prefix.",
                                     ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "output-file", options.outputFormatList);
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
        append(file, options.outputFormatList[options.outputFormat]);
    }
}

// ----------------------------------------------------------------------------
// Function setTmpFolder()
// ----------------------------------------------------------------------------

void setTmpFolder(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("t", "tmp-folder", "Specify a temporary folder where to construct the index. \
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

#endif  // #ifndef APP_YARA_MISC_OPTIONS_H_
