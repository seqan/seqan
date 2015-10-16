// ==========================================================================
//                              generalProcessing.h
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// ==========================================================================

#ifndef GENERALPROCESSING_H
#define GENERALPROCESSING_H

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "helper_functions.h"
#include "general_stats.h"

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


// ============================================================================
// Functions
// ============================================================================

struct NoSubstitute {};

template<typename TSeqChar, typename TSub>
constexpr void replaceN(TSeqChar& seqChar, const TSub sub) noexcept
{
    seqChar = sub;
}

template<typename TSeqChar>
constexpr void replaceN(TSeqChar& seqChar, const NoSubstitute) noexcept
{
    (void)seqChar;
}

template<typename TSeq, typename TSub>
inline int findNUniversal(TSeq& seq, unsigned allowed, const TSub substitute) noexcept
{
    const TSeq wanted = 'N';
    unsigned c = 0;
    for (auto& seqChar : seq)
    {
        if (seqChar == wanted)
        {
            replaceN(seqChar, substitute);
            if (++c > allowed)
                return -1;       //sequence will be removed
        }
    }
    return c;                   //sequence not deleted, number of substitutions returned
}

template<template<typename> class TRead, typename TSeq, typename TSub,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>> ::value >>
constexpr inline int findNMultiplex(TRead<TSeq>& read, unsigned allowed, const TSub substitute, bool = false) noexcept
{
    (void)read;
    (void)allowed;
    (void)substitute;
    return 0;
}

template<template<typename> class TRead, typename TSeq, typename TSub,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq >> ::value >>
inline int findNMultiplex(TRead<TSeq>& read, unsigned allowed, const TSub substitute) noexcept
{
    return findNUniversal(read.demultiplex, allowed, substitute);
}

template<template<typename> class TRead, typename TSeq, typename TSub,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplex<TSeq >> ::value >>
constexpr inline int findNPairedEnd(TRead<TSeq>& read, unsigned allowed, const TSub substitute, bool = false) noexcept
{
    (void)read;
    (void)allowed;
    (void)substitute;
    return 0;
}

template<template<typename> class TRead, typename TSeq, typename TSub,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq >> ::value >>
inline int findNPairedEnd(TRead<TSeq>& read, unsigned allowed, const TSub substitute) noexcept
{
    return findNUniversal(read.seqRev, allowed, substitute);
}

template<template<typename> class TRead, typename TSeq, typename TSub>
int findN(TRead<TSeq>& read, unsigned allowed, const TSub substitute) noexcept
{
    auto c = findNUniversal(read.seq, allowed, substitute);
    if (c == -1)
        return -1;
    const auto c2 = findNMultiplex(read, allowed, substitute);
    if (c2 == -1)
        return -1;
    c += c2;
    const auto c3 = findNPairedEnd(read, allowed, substitute);
    if (c3 == -1)
        return -1;
    c += c3;
    return c;
}

//universal function for all combinations of options
template<template <typename> class TRead, typename TSeq, typename TSub, typename TStats>
void processN(std::vector<TRead<TSeq>>& reads, unsigned allowed, TSub substitute, TStats& stats) noexcept
{
    std::vector<int> res(length(reads));
    std::transform(reads.begin(), reads.end(), res.begin(),[allowed, substitute, &stats](auto& read) ->auto
    {
        auto const resElement = findN(read, allowed, substitute);
        if(resElement != -1)
            stats.uncalledBases += resElement;
        return resElement;
    });
    stats.removedN += _eraseSeqs(res, -1, reads);
}

template<template <typename> class TRead, typename TSeq>
unsigned int removeShortSeqs(std::vector<TRead<TSeq>>& reads, const unsigned min) noexcept
{
    const auto numReads = (int)length(reads);
    reads.erase(std::remove_if(reads.begin(), reads.end(), [min](const auto& read) {return read.minSeqLen() < min;}), reads.end());
    return numReads - length(reads);
}

// main preTrim function
template<template <typename> class TRead, typename TSeq, bool tagTrimming,
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplex < TSeq >> ::value >>
unsigned int _preTrim(std::vector<TRead<TSeq>>& reads, const unsigned head, const unsigned tail, const unsigned min, bool = false) noexcept(!tagTrimming)
{
    std::string insertToken;
    if(tagTrimming)
        insertToken.reserve(8 + head + tail);
    std::for_each(reads.begin(), reads.end(), [head, tail, &insertToken](auto& read)
    {
        auto seqLen = length(read.seq);
        if (seqLen > (head + tail))
        {
            if (head > 0)
            {
                if (tagTrimming)
                    insertToken = ":TL:" + std::string(prefix(read.seq, head));
                erase(read.seq, 0, head);
            }
            if (tail > 0)
            {
                seqLen = length(read.seq);
                if (tagTrimming)
                    insertToken += ":TR:" + std::string(suffix(read.seq, seqLen - tail));
                erase(read.seq, seqLen - tail, seqLen);
            }
            if (tagTrimming && insertToken.size() != 0)
                insertAfterFirstToken(read.id, std::move(insertToken));
        }
        else
            clear(read.seq);
    });

    return removeShortSeqs(reads, min);
}

template<template <typename> class TRead, typename TSeq, bool tagTrimming,
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplexPairedEnd < TSeq >> ::value >>
unsigned int _preTrim(std::vector<TRead<TSeq>>& reads, const unsigned head, const unsigned tail, const unsigned min) noexcept(!tagTrimming)
{
    std::string tempString;
    if(tagTrimming)
        tempString.reserve(8 + head + tail);
    std::for_each(reads.begin(), reads.end(), [head, tail, &tempString](auto& read)
    {
        if (read.minSeqLen() > (head + tail))
        {
            std::string insertToken;
            if (head > 0)
            {
                if (tagTrimming)
                {
                    tempString = ":TL:";
                    append(tempString, prefix(read.seq, std::min<int>(length(read.seq), head)));
                    insertAfterFirstToken(read.id, std::move(tempString));
                    tempString = ":TL:";
                    append(tempString, prefix(read.seqRev, std::min<int>(length(read.seqRev), head)));
                    insertAfterFirstToken(read.idRev, std::move(tempString));
                }
                erase(read.seq, 0, std::min<int>(length(read.seq),head));
                erase(read.seqRev, 0, std::min<int>(length(read.seqRev),head));
            }
            if (tail > 0)
            {
                const auto seqLen = length(read.seq);
                const auto seqLenRev = length(read.seqRev);
                if (tagTrimming)
                {
                    tempString = ":TR:";
                    append(tempString, suffix(read.seq, std::max<int>(0, seqLen - tail)));
                    insertAfterFirstToken(read.id, std::move(tempString));
                    tempString = ":TR:";
                    append(tempString, suffix(read.seqRev, std::max<int>(0, seqLenRev - tail)));
                    insertAfterFirstToken(read.idRev, std::move(tempString));
                }
                erase(read.seq, std::max<int>(0, seqLen - tail), seqLen);
                erase(read.seqRev, std::max<int>(0, seqLenRev - tail), seqLenRev);
            }
            if (tagTrimming && insertToken.size() != 0)
                insertAfterFirstToken(read.id, std::move(insertToken));
        }
        else
        {
            clear(read.seq);
            clear(read.seqRev);
        }
    });
    return removeShortSeqs(reads, min);
}

template <template<typename> class TRead, typename TSeq, typename TStats>
void preTrim(std::vector<TRead<TSeq>>& reads, const unsigned head, const unsigned tail, const unsigned min, const bool tagTrimming, TStats& stats)
{
    if(tagTrimming)
        stats.removedShort += _preTrim<TRead, TSeq, true>(reads, head, tail, min);
    else
        stats.removedShort += _preTrim<TRead, TSeq, false>(reads, head, tail, min);
}

//Trims sequences to specific length and deletes to short ones together with their IDs
template<template <typename> class TRead, typename TSeq, typename TStats,
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplex < TSeq >> ::value >>
void trimTo(std::vector<TRead<TSeq>>& reads, const unsigned len, TStats& stats, bool = true) noexcept
{
    for(auto& read : reads)
        if (read.minSeqLen() > len)
            erase(read.seq, len, length(read.seq));

    stats.removedShort += removeShortSeqs(reads, len);
}

template<template <typename> class TRead, typename TSeq, typename TStats,
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplexPairedEnd < TSeq >> ::value >>
void trimTo(std::vector<TRead<TSeq>>& reads, const unsigned len, TStats& stats) noexcept
{
    for (auto& read : reads)
    {
            if (read.minSeqLen() > len)
            {
                if (length(read.seq) > len)
                    erase(read.seq, len, length(read.seq));
                if (length(read.seqRev) > len)
                    erase(read.seqRev, len, length(read.seqRev));
            }
    }

    stats.removedShort += removeShortSeqs(reads, len);
}


#endif
