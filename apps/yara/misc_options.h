// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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

#ifdef STDLIB_VS
#include <direct.h>
#endif

using namespace seqan;

// ============================================================================
// Functors
// ============================================================================

typedef EqualsChar<'.'>        IsDot;
typedef EqualsChar<'/'>        IsSlash;
typedef EqualsChar<'\\'>       IsBackSlash;

#ifdef STDLIB_VS
    typedef IsBackSlash        IsPathDelimited;
#else
    typedef IsSlash            IsPathDelimited;
#endif

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setEnv()
// ----------------------------------------------------------------------------

template <typename TString, typename TValue>
bool setEnv(TString const & key, TValue & value)
{
#ifdef STDLIB_VS
    return !_putenv_s(toCString(key), toCString(value));
#else
    return !setenv(toCString(key), toCString(value), true);
#endif
}

// ----------------------------------------------------------------------------
// Function getCwd()
// ----------------------------------------------------------------------------

template <typename TString>
void getCwd(TString & string)
{
    char cwd[1000];

#ifdef STDLIB_VS
    _getcwd(cwd, 1000);
#else
    ignoreUnusedVariableWarning(getcwd(cwd, 1000));
#endif

    assign(string, cwd);
}

// ----------------------------------------------------------------------------
// Function firstOf()
// ----------------------------------------------------------------------------

template <typename TString, typename TFunctor>
inline typename Iterator<TString const, Standard>::Type
firstOf(TString const & string, TFunctor const & f)
{
    typedef typename Iterator<TString const, Standard>::Type TIter;

    TIter it = begin(string, Standard());
    skipUntil(it, f);

    return it;
}

// ----------------------------------------------------------------------------
// Function lastOf()
// ----------------------------------------------------------------------------

template <typename TString, typename TFunctor>
inline typename Iterator<TString const, Standard>::Type
lastOf(TString const & string, TFunctor const & f)
{
    typedef ModifiedString<TString const, ModReverse>        TStringRev;
    typedef typename Iterator<TStringRev, Standard>::Type    TIterRev;

    TStringRev revString(string);
    TIterRev revIt = firstOf(revString, f);

    return end(string) - position(revIt, revString);
}

// ----------------------------------------------------------------------------
// Function trimExtension()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Prefix<TString const>::Type
trimExtension(TString const & string)
{
    return prefix(string, firstOf(string, IsDot()));
}

// ----------------------------------------------------------------------------
// Function getExtension()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Suffix<TString const>::Type
getExtension(TString const & string)
{
    return suffix(string, firstOf(string, IsDot()) + 1);
}

// ----------------------------------------------------------------------------
// Function getPath()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Prefix<TString const>::Type
getPath(TString const & string)
{
    typedef typename Iterator<TString const, Standard>::Type TIter;

    TIter it = lastOf(string, IsPathDelimited());

    if (it != begin(string, Standard())) --it;

    return prefix(string, it);
}

// ----------------------------------------------------------------------------
// Function getFilename()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Suffix<TString const>::Type
getFilename(TString const & string)
{
    return suffix(string, lastOf(string, IsPathDelimited()));
}

// ----------------------------------------------------------------------------
// Function addLeadingDot()
// ----------------------------------------------------------------------------

template <typename TString>
inline void addLeadingDot(TString & string)
{
    insert(string, 0, ".");
}

// ----------------------------------------------------------------------------
// Function stripLeadingDot()
// ----------------------------------------------------------------------------

template <typename TString>
inline void stripLeadingDot(TString & string)
{
    string.erase(0, 1);
}

// ----------------------------------------------------------------------------
// Function getExtensionsWithoutLeadingDot()
// ----------------------------------------------------------------------------

template <typename TStrings>
inline TStrings getExtensionsWithoutLeadingDot(TStrings const & strings)
{
    typedef typename Value<TStrings>::Type  TString;

    TStrings extensions = strings;
    forEach(extensions, [](TString & extension) { stripLeadingDot(extension); });

    return extensions;
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
// Function getReadErrors()
// ----------------------------------------------------------------------------
// Returns the absolute number of errors for a given read sequence.

template <typename TMatch, typename TOptions, typename TReadSeqSize>
inline TReadSeqSize getReadErrors(TOptions const & options, TReadSeqSize readSeqLength)
{
    return std::min((TReadSeqSize)(readSeqLength * options.errorRate),
                    (TReadSeqSize)MemberLimits<TMatch, Errors>::VALUE);
}

// ----------------------------------------------------------------------------
// Function getReadIndels()
// ----------------------------------------------------------------------------
// Returns the absolute number of indels for a given read sequence.

template <typename TMatch, typename TOptions, typename TReadSeqSize>
inline TReadSeqSize getReadIndels(TOptions const & options, TReadSeqSize readSeqLength)
{
    return std::min((TReadSeqSize)(readSeqLength * options.indelRate),
                    (TReadSeqSize)MemberLimits<TMatch, Errors>::VALUE);
}

// ----------------------------------------------------------------------------
// Function getReadStrata()
// ----------------------------------------------------------------------------
// Returns the absolute number of strata for a given read sequence.

template <typename TMatch, typename TOptions, typename TReadSeqSize>
inline TReadSeqSize getReadStrata(TOptions const & options, TReadSeqSize readSeqLength)
{
    return std::min((TReadSeqSize)(readSeqLength * options.strataRate),
                    (TReadSeqSize)MemberLimits<TMatch, Errors>::VALUE);
}

// ----------------------------------------------------------------------------
// Function saveContigsLimits()
// ----------------------------------------------------------------------------

template <typename TOptions>
bool saveContigsLimits(TOptions const & options)
{
    String<uint64_t> limits;

    appendValue(limits, options.contigsMaxLength);
    appendValue(limits, options.contigsSize);
    appendValue(limits, options.contigsSum);

    CharString contigsLimitFile(options.contigsIndexFile);
    append(contigsLimitFile, ".txt.size");

    return save(limits, toCString(contigsLimitFile));
}

// ----------------------------------------------------------------------------
// Function openContigsLimits()
// ----------------------------------------------------------------------------

template <typename TOptions>
bool openContigsLimits(TOptions & options)
{
    String<uint64_t> limits;

    CharString contigsLimitFile(options.contigsIndexFile);
    append(contigsLimitFile, ".txt.size");

    if (!open(limits, toCString(contigsLimitFile), OPEN_RDONLY))
        return false;

    if (length(limits) != 3)
        return false;

    options.contigsMaxLength = limits[0];
    options.contigsSize = limits[1];
    options.contigsSum = limits[2];

    return true;
}

// ----------------------------------------------------------------------------
// Function setContigsLimits()
// ----------------------------------------------------------------------------

template <typename TOptions, typename TSeqs>
void setContigsLimits(TOptions & options, TSeqs const & seqs)
{
    options.contigsMaxLength = maxLength(seqs);
    options.contigsSize = length(seqs);
    options.contigsSum = lengthSum(seqs);
}

// ----------------------------------------------------------------------------
// Function setDateAndVersion()
// ----------------------------------------------------------------------------

void setDateAndVersion(ArgumentParser & parser)
{
    setCategory(parser, "Read Mapping");

#if defined(SEQAN_APP_VERSION) && defined(SEQAN_REVISION)
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
#endif
#if defined(SEQAN_DATE)
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
    addDescription(parser, "(c) Copyright 2011-2014 by Enrico Siragusa.");
    addDescription(parser, "(c) Copyright 2013 by NVIDIA Corporation.");
}

#endif  // #ifndef APP_YARA_MISC_OPTIONS_H_
