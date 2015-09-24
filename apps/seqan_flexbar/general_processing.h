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
#include <initializer_list>


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct GeneralStats
{
	unsigned removedSeqs;       //Number of deleted sequences due to N's
    unsigned long uncalledBases;//Number of uncalled bases (evtl. Masked) in surviving sequences
    unsigned removedSeqsShort;  //Number of deleted sequences due to shortness.
    unsigned int readCount;
    double processTime;
    double ioTime;

    GeneralStats() : removedSeqs(0), uncalledBases(0), removedSeqsShort(0), readCount(0), processTime(0), ioTime(0) {};
};

// ============================================================================
// Functions
// ============================================================================



template<typename TAlpha, typename TSub>
int findN(seqan::String<TAlpha>& seq, unsigned allowed, TSub substitute)
{
    unsigned limit = length(seq);
    TAlpha wanted = 'N';
    unsigned c = 0;
    for (unsigned i = 0; i < limit; ++i)
    {
        if (seq[i] == wanted)
        {
            seq[i] = substitute;
            ++c;
            if (c > allowed)
            {
                return -1;       //sequence will be removed
            }
        }
    }
    return c;                   //sequence not deleted, number of substitutions returned
}

//Overload if no substitution shall be performed
template<typename TAlpha>
int findN(seqan::String<TAlpha>& seq, unsigned allowed)
{
    unsigned limit = length(seq);
    TAlpha wanted = 'N';
    unsigned c = 0;
    for (unsigned i = 0; i < limit; ++i)
    {
        if (seq[i] == wanted)
        {
            ++c;
            if (c > allowed)
            {
                return -1;       //sequence will be removed
            }
        }
    }
    return c;                   //sequence not deleted, number of N's returned
}

//single-end data with substitutions
template<typename TSeqs, typename TIds, typename TSub>
void processN(TSeqs& seqs, TIds& ids, unsigned allowed, TSub substitute, GeneralStats& stats)
{
    int limit = length(seqs);
    seqan::StringSet<int> res;
    resize(res, limit);
    SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static))
    for (int i = 0; i < limit; ++i)
    {
        res[i] = findN(seqs[i], allowed, substitute);
    }
    //Simply erasing a sequence would lead to unnecessary reallocations, therefore we use another method
    unsigned ex = 0;
    for (int i = length(res) - 1; i >= 0 ; --i)
    {                                   //integer necessary because unsigned would cause error after last iteration
        if (res[i] == -1)                //Placing all sequences/ids which shall be erased at the end of their container
        {
            SEQAN_OMP_PRAGMA(parallel default(shared))
            {
                SEQAN_OMP_PRAGMA(sections nowait)
                {
                    SEQAN_OMP_PRAGMA(section)     //distributing the swapping actions
                    swap(seqs[i], seqs[limit - ex - 1]);
                    SEQAN_OMP_PRAGMA(section)
                    swap(ids[i], ids[limit - ex - 1]);
                    SEQAN_OMP_PRAGMA(section)
                    ++ex;
                }
            }
        }
        else
        {
            stats.uncalledBases += res[i];
        } 
    }
    if (ex != 0)
    {
        SEQAN_OMP_PRAGMA(parallel default(shared))
        {
            SEQAN_OMP_PRAGMA(sections nowait)
            {   //Resizing the containers to erase all unwanted sequences/ids at once
                SEQAN_OMP_PRAGMA(section)
                resize(seqs, limit - ex);
                SEQAN_OMP_PRAGMA(section)
                resize(ids, limit - ex);
                SEQAN_OMP_PRAGMA(section)
                stats.removedSeqs += ex;
            }
        }
    }
}

//Overload for single-end data and no substitutions
template<typename TSeqs, typename TIds>
void processN(TSeqs& seqs, TIds& ids, unsigned allowed, GeneralStats& stats)
{
    int limit = length(seqs);
    seqan::StringSet<int> res;
    resize(res, limit);
    unsigned uncalled = 0;
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static) reduction(+:uncalled))
    for (int i = 0; i < limit; ++i)
    {
        res[i] = findN(seqs[i], allowed);
        if (res[i] != -1)
            uncalled += res[i];
    }
    stats.uncalledBases += uncalled;
    stats.removedSeqs += _eraseSeqs(res, -1, seqs, ids);
}

//paired-end data with substitutions
template<typename TSeqs, typename TIds, typename TSub>
void processN(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, unsigned allowed, TSub substitute, GeneralStats& stats)
{
    int limit = length(seqs);
    seqan::StringSet<int> res;
    resize(res, limit);
    unsigned uncalled = 0;
    SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static) reduction(+:uncalled))
    for (int i = 0; i < limit; ++i)
    {
        res[i] = findN(seqs[i], allowed, substitute);
        if (res[i] != -1)
        {
            int resTemp = res[i];
            res[i] = findN(seqsRev[i], allowed, substitute);
            if (res[i] != -1)
            {
                uncalled += (res[i] + resTemp);
            }
        }
    }
    stats.uncalledBases += uncalled;
    stats.removedSeqs += 2 * _eraseSeqs(res, -1, seqs, seqsRev, ids, idsRev);
}

//paired-end data without substitutions
template<typename TSeqs, typename TIds>
void processN(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, unsigned allowed, GeneralStats& stats)
{
    int limit = length(seqs);
    seqan::StringSet<int> res;
    resize(res, limit);
    unsigned uncalled = 0;
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static) reduction(+:uncalled))
    for (int i = 0; i < limit; ++i)
    {
        res[i] = findN(seqs[i], allowed);
        if (res[i] != -1)
        {
            int resTemp = res[i];
            res[i] = findN(seqsRev[i], allowed);
            if (res[i] != -1)
            {
                uncalled += (res[i] + resTemp);
            }
        }
    }
    stats.uncalledBases += uncalled;
    stats.removedSeqs += 2 * _eraseSeqs(res, -1, seqs, seqsRev, ids, idsRev);
}

//Overload single-end data with substitutions and multiplex barcodes
template<typename TSeqs, typename TIds, typename TMulti, typename TSub>
void processN(TSeqs& seqs, TIds& ids, TMulti& multiplex, unsigned allowed, TSub substitute, GeneralStats& stats)
{
    int limit = length(seqs);
    seqan::StringSet<int> res;
    resize(res, limit);
    unsigned uncalled = 0;
    SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static) reduction(+:uncalled))
    for (int i = 0; i < limit; ++i)
    {
        res[i] = findN(seqs[i], allowed, substitute);
        if (res[i] != -1)
        {
            int resTemp = res[i];
            res[i] = findN(multiplex[i], allowed, substitute);
            if (res[i] != -1)
            {
                uncalled += (res[i] + resTemp);
            }
        }
    }
    stats.uncalledBases += uncalled;
    stats.removedSeqs += _eraseSeqs(res, -1, seqs, ids, multiplex);
}

//Overload for single-end data, no substitutions and multiplex barcodes
template<typename TSeqs, typename TIds, typename TMulti>
void processN(TSeqs& seqs, TIds& ids, TMulti& multiplex, unsigned allowed, GeneralStats& stats)
{
    int limit = length(seqs);
    seqan::StringSet<int> res;
    resize(res, limit);
    unsigned uncalled = 0;
    SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static) reduction(+:uncalled))
    for (int i = 0; i < limit; ++i)
    {
        res[i] = findN(seqs[i], allowed);
        if (res[i] != -1)
        {
            int resTemp = res[i];
            res[i] = findN(multiplex[i], allowed);
            if (res[i] != -1)
            {
                uncalled += (res[i] + resTemp);
            }
        }
    }
    stats.uncalledBases += uncalled;
    stats.removedSeqs += _eraseSeqs(res, -1, seqs, ids, multiplex);
}

//paired-end data with substitutions and multiplex barcodes
template<typename TSeqs, typename TIds, typename TMulti, typename TSub>
void processN(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, TMulti& multiplex, unsigned allowed,
    TSub substitute, GeneralStats& stats)
{
    int limit = length(seqs);
    seqan::StringSet<int> res;
    resize(res, limit);
    unsigned uncalled = 0;
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static) reduction(+:uncalled))
    for (int i = 0; i < limit; ++i)
    {
        res[i] = findN(seqs[i], allowed, substitute);
        if (res[i] != -1)
        {
            int resTemp = res[i];
            res[i] = findN(seqsRev[i], allowed, substitute);
            if (res[i] != -1)
            {
                resTemp += res[i];
                res[i] = findN(multiplex[i], allowed, substitute);
                if (res[i] != -1)
                {
                    uncalled += (res[i] + resTemp);
                }
            }
        }
    }
    stats.uncalledBases += uncalled;
    stats.removedSeqs += 2 * _eraseSeqs(res, -1, seqs, seqsRev, ids, idsRev, multiplex);
}

//paired-end data without substitutions and multiplex barcodes
template<typename TSeqs, typename TIds, typename TMulti >
void processN(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, TMulti& multiplex, unsigned allowed, GeneralStats& stats)
{
    int limit = length(seqs);
    seqan::StringSet<int> res;
    resize(res, limit);
    unsigned uncalled = 0;
    SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static) reduction(+:uncalled))
    for (int i = 0; i < limit; ++i)
    {
        res[i] = findN(seqs[i], allowed);
        if (res[i] != -1)
        {
            int resTemp = res[i];
            res[i] = findN(seqsRev[i], allowed);
            if (res[i] != -1)
            {
                resTemp += res[i];
                res[i] = findN(multiplex[i], allowed);
                if (res[i] != -1)
                {
                    uncalled += (res[i] + resTemp);
                }
            }
        }
    }
    stats.uncalledBases += uncalled;
    stats.removedSeqs += 2 * _eraseSeqs(res, -1, seqs, seqsRev, ids, idsRev, multiplex);
}

template<typename TRead>    
unsigned int postTrim(std::vector<TRead>& reads, const unsigned min)
{
    const auto numReads = (int)length(reads);
    reads.erase(std::remove_if(reads.begin(), reads.end(), [min](const auto& read) {return length(read.seq) >= min;}), reads.end());
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
                if (length(readSet[i].seq) >= min)
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
                const auto seqLen = std::min(length(readSet[i].seq), length(readSet[i].seqRev));
                if (seqLen >= min)
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
    stats.removedSeqsShort += _eraseSeqs(rem, true, reads);
}

//Trims sequences to specific length and deletes to short ones together with their IDs
template<template <typename> class TRead, typename TSeq, 
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, Read<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplex < TSeq >> ::value >>
    void trimTo(std::vector<TRead<TSeq>>& reads, const unsigned len, GeneralStats& stats, bool = true) 
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
    stats.removedSeqsShort += _eraseSeqs(rem, true, reads);
}

template<template <typename> class TRead, typename TSeq,
    typename = std::enable_if_t < std::is_same<TRead<TSeq>, ReadPairedEnd<TSeq>>::value || std::is_same < TRead<TSeq>, ReadMultiplexPairedEnd < TSeq >> ::value >>
    void trimTo(std::vector<TRead<TSeq>>& reads, const unsigned len, GeneralStats& stats) 
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
    stats.removedSeqsShort += _eraseSeqs(rem, true, reads);
}


#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_GENERALPROCESSING_H_
