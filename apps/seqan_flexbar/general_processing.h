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
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================
// This file provides functions used by different parts of seqan-flexbar
// which is based in the implementation of the original flexbar program
// in [1].
// [1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, C.  FLEXBAR—Flexible
// Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
// Biology 2012, 1, 895-905.
// ==========================================================================



#ifndef SANDBOX_GROUP3_APPS_SEQDPT_GENERALPROCESSING_H_
#define SANDBOX_GROUP3_APPS_SEQDPT_GENERALPROCESSING_H_

#ifdef _OPENMP
#include <omp.h>
#endif
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>

#include <seqan/stream.h>

#include "helper_functions.h"
#include "general_stats.h"
#include <initializer_list>


// ============================================================================
// Tags, Classes, Enums
// ============================================================================


// ============================================================================
// Functions
// ============================================================================

struct NoSubstitute {};

template<typename TSeqChar, typename TSub>
void replaceN(TSeqChar& seqChar, const TSub sub) noexcept
{
    seqChar = sub;
}

template<typename TSeqChar>
void replaceN(TSeqChar& seqChar, const NoSubstitute) noexcept
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
            ++c;
            if (c > allowed)
                return -1;       //sequence will be removed
        }
    }
    return c;                   //sequence not deleted, number of substitutions returned
}

template<template<typename> typename TRead, typename TSeq, typename TSub,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>> ::value >>
inline int findNMultiplex(TRead<TSeq>& read, unsigned allowed, const TSub substitute, bool = false) noexcept
{
    return 0;
}

template<template<typename> typename TRead, typename TSeq, typename TSub,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadMultiplex<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq >> ::value >>
inline int findNMultiplex(TRead<TSeq>& read, unsigned allowed, const TSub substitute) noexcept
{
    return findNUniversal(read.demultiplex, allowed, substitute);
}

template<template<typename> typename TRead, typename TSeq, typename TSub,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplex<TSeq >> ::value >>
inline int findNPairedEnd(TRead<TSeq>& read, unsigned allowed, const TSub substitute, bool = false) noexcept
{
    return 0;
}

template<template<typename> typename TRead, typename TSeq, typename TSub,
    typename = std::enable_if_t<std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same<TRead<TSeq>, ReadMultiplexPairedEnd<TSeq >> ::value >>
inline int findNPairedEnd(TRead<TSeq>& read, unsigned allowed, const TSub substitute) noexcept
{
    return findNUniversal(read.seqRev, allowed, substitute);
}

template<template<typename> typename TRead, typename TSeq, typename TSub>
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
    return c;                   //sequence not deleted, number of substitutions returned
}

//universal function for all combinations of options
template<template <typename> typename TRead, typename TSeq, typename TSub>
void processN(std::vector<TRead<TSeq>>& reads, unsigned allowed, TSub substitute, GeneralStats& stats) noexcept
{
    int limit = length(reads);
    std::vector<int> res(limit);
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static))
        for (int i = 0; i < limit; ++i)
        {
            res[i] = findN(reads[i], allowed, substitute);
        }

    unsigned ex = 0;
    stats.removedN += _eraseSeqs(res, -1, reads);
    for (int i = length(res) - 1; i >= 0; --i)
    {
        if (res[i] != -1)
            stats.uncalledBases += res[i];
    }
}

template<template <typename> class TRead, typename TSeq>
    unsigned int postTrim(std::vector<TRead<TSeq>>& reads, const unsigned min) noexcept
{
    const auto numReads = (int)length(reads);
    reads.erase(std::remove_if(reads.begin(), reads.end(), [min](const auto& read) {return read.minSeqLen() < min;}), reads.end());
    return numReads - length(reads);
}

// main preTrim function
template<template <typename> class TRead, typename TSeq, bool tagTrimming,
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplex < TSeq >> ::value >>
void _preTrim(std::vector<TRead<TSeq>>& readSet, const unsigned head, const unsigned tail, const unsigned min, std::vector<bool>& rem, bool = false) noexcept(!tagTrimming)
{
    int i = 0;
    const auto limit = (int)length(readSet);
    assert(rem.size() == length(readSet));

    SEQAN_OMP_PRAGMA(parallel for default(shared) private(i)schedule(static))
        for (i = 0; i < limit; ++i)
        {
            const auto seqLen = length(readSet[i].seq);
            if (seqLen >(head + tail))
            {
                std::string insertToken;
                if (head > 0)
                {
                    if (tagTrimming)
                        insertToken = ":TL:" + std::string(prefix(readSet[i].seq, head));
                    erase(readSet[i].seq, 0, head);
                }
                if (tail > 0)
                {
                    const auto seqLen = length(readSet[i].seq);
                    if (tagTrimming)
                        insertToken += ":TR:" + std::string(prefix(readSet[i].seq, seqLen - tail));
                    erase(readSet[i].seq, seqLen - tail, seqLen);
                }
                if(!empty(insertToken))
                    insertAfterFirstToken(readSet[i].id, std::move(insertToken));

                // check if trimmed sequence is at least of length min
                // if not, remove it
                if (readSet[i].minSeqLen() >= min)
                    rem[i] = false;
                else
                    rem[i] = true;
            }
            // if the sequence was too short to be trimmed, remove it
            else
                rem[i] = true;
        }
}

template<template <typename> class TRead, typename TSeq, bool tagTrimming,
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplexPairedEnd < TSeq >> ::value >>
    void _preTrim(std::vector<TRead<TSeq>>& readSet, const unsigned head, const unsigned tail, const unsigned min, std::vector<bool>& rem) noexcept(!tagTrimming)
{
    int i = 0;
    const auto limit = (int)length(readSet);
    assert(rem.size() == length(readSet));

    SEQAN_OMP_PRAGMA(parallel for default(shared) private(i)schedule(static))
        for (i = 0; i < limit; ++i)
        {
            const auto seqLen = std::min(length(readSet[i].seq), length(readSet[i].seqRev));
            if (seqLen > (head + tail))
            {
                if (head > 0)
                {
                    if (tagTrimming)
                    {
                        std::string tempString = ":TL:";
                        append(tempString, prefix(readSet[i].seq, head));
                        insertAfterFirstToken(readSet[i].id, std::move(tempString));
                        tempString = ":TL:";
                        append(tempString, prefix(readSet[i].seqRev, head));
                        insertAfterFirstToken(readSet[i].idRev, std::move(tempString));
                    }
                    erase(readSet[i].seq, 0, head);
                    erase(readSet[i].seqRev, 0, head);
                }
                if (tail > 0)
                {
                    const auto seqLen = length(readSet[i].seq);
                    if (tagTrimming)
                    {
                        std::string tempString = ":TR:";
                        append(tempString, suffix(readSet[i].seq, seqLen - tail));
                        insertAfterFirstToken(readSet[i].id, std::move(tempString));
                        tempString = ":TR:";
                        append(tempString, suffix(readSet[i].seqRev, seqLen - tail));
                        insertAfterFirstToken(readSet[i].idRev, std::move(tempString));
                    }
                    erase(readSet[i].seq, seqLen - tail, seqLen);
                    erase(readSet[i].seqRev, seqLen - tail, seqLen);
                }

                // check if trimmed sequence is at least of length min
                // if not, remove it
                if (readSet[i].minSeqLen() >= min)
                    rem[i] = false;
                else
                    rem[i] = true;
            }
            // if the sequence was too short to be trimmed, remove it
            else
                rem[i] = true;
        }
}

template <template<typename> class TRead, typename TSeq>
void preTrim(std::vector<TRead<TSeq>>& reads, unsigned head, unsigned tail, unsigned min, const bool tagTrimming, GeneralStats& stats)
{
    std::vector<bool> rem(length(reads));
    if(tagTrimming)
        _preTrim<TRead, TSeq, true>(reads, head, tail, min, rem);
    else
        _preTrim<TRead, TSeq, false>(reads, head, tail, min, rem);
    stats.removedShort += _eraseSeqs(rem, true, reads);
}

//Trims sequences to specific length and deletes to short ones together with their IDs
template<template <typename> class TRead, typename TSeq, 
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplex < TSeq >> ::value >>
    void trimTo(std::vector<TRead<TSeq>>& reads, const unsigned len, GeneralStats& stats, bool = true) noexcept
{
    const auto limit = (int)length(reads);
    std::vector<bool> rem(limit);
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static))
        for (int i = 0; i < limit; ++i)
        {
            if (length(reads[i].seq) < len)
            {
                rem[i] = true;
            }
            else
            {
                rem[i] = false;
                if (length(reads[i].seq) > len)
                {
                    erase(reads[i].seq, len, length(reads[i].seq));
                }
            }
        }
    stats.removedShort += _eraseSeqs(rem, true, reads);
}

template<template <typename> class TRead, typename TSeq,
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplexPairedEnd < TSeq >> ::value >>
    void trimTo(std::vector<TRead<TSeq>>& reads, const unsigned len, GeneralStats& stats) noexcept
{
    const auto limit = (int)length(reads);
    std::vector<bool> rem(limit);
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static))
        for (int i = 0; i < limit; ++i)
        {
            if (std::min(length(reads[i].seq), length(reads[i].seqRev)) < len)
            {
                rem[i] = true;
            }
            else
            {
                rem[i] = false;
                if (length(reads[i].seq) > len)
                    erase(reads[i].seq, len, length(reads[i].seq));
                if (length(reads[i].seqRev) > len)
                    erase(reads[i].seqRev, len, length(reads[i].seqRev));
            }
        }
    stats.removedShort += _eraseSeqs(rem, true, reads);
}


#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_GENERALPROCESSING_H_
