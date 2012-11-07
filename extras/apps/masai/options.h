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

    MasaiOptions()
    {
        indexTypeList.push_back("esa");
        indexTypeList.push_back("sa");
//        indexTypeList.push_back("qgram");
        indexTypeList.push_back("fm");

        mappingModeList.push_back("all");
        mappingModeList.push_back("all-best");
        mappingModeList.push_back("any-best");

        outputFormatList.push_back("raw");
        outputFormatList.push_back("sam");
        outputFormatList.push_back("sam-no-cigar");
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getOptionValue()
// ----------------------------------------------------------------------------

template <typename TOption, typename TString, typename TOptionsList>
void getOptionValue(TOption & option, ArgumentParser & parser, TString const & optionName, TOptionsList & optionsList)
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
    setVersion(parser, "0.5 [" + rev.substr(11, rev.size() - 13) + "]");
    setDate(parser, date.substr(7, std::min((int)date.size() - 8, 10)));
}

// ----------------------------------------------------------------------------
// Function setDescription()
// ----------------------------------------------------------------------------

void setDescription(ArgumentParser & parser)
{
    addDescription(parser, "Masai is a fast and sensitive read mapper based on approximate seeds and multiple backtracking.");
    addDescription(parser, "See \\fIhttp://www.seqan.de/projects/masai\\fP for more information.");
    addDescription(parser, "(c) Copyright 2011-2012 by Enrico Siragusa.");
}

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

#endif  // #ifndef SEQAN_EXTRAS_MASAI_OPTIONS_H_
