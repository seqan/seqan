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
// Functors
// ============================================================================

typedef EqualsChar<'.'>        IsDot;
typedef EqualsChar<'/'>        IsSlash;
typedef EqualsChar<'\\'>       IsBackSlash;

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
    return !_putenv_s(toCString(key), toCString(value));
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
    setVersion(parser, "0.9.1 [" + std::string(SEQAN_REVISION) + "]");
#else
    setVersion(parser, "0.9.1");
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

#endif  // #ifndef APP_YARA_MISC_OPTIONS_H_
