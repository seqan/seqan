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

    GeneralStats() : removedSeqs(0), uncalledBases(0), removedSeqsShort(0) {};
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
    StringSet<int> res;
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
    StringSet<int> res;
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
    StringSet<int> res;
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
    StringSet<int> res;
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
    StringSet<int> res;
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
    StringSet<int> res;
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
    stats.removedSeqs += _eraseSeqs(res, -1, seqs, ids);
}

//paired-end data with substitutions and multiplex barcodes
template<typename TSeqs, typename TIds, typename TMulti, typename TSub>
void processN(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, TMulti& multiplex, unsigned allowed,
    TSub substitute, GeneralStats& stats)
{
    int limit = length(seqs);
    StringSet<int> res;
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
    stats.removedSeqs += 2 * _eraseSeqs(res, -1, seqs, seqsRev, ids, idsRev);
}

//paired-end data without substitutions and multiplex barcodes
template<typename TSeqs, typename TIds, typename TMulti >
void processN(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, TMulti& multiplex, unsigned allowed, GeneralStats& stats)
{
    int limit = length(seqs);
    StringSet<int> res;
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
    stats.removedSeqs += 2 * _eraseSeqs(res, -1, seqs, seqsRev, ids, idsRev);
}



// main preTrim function
template<typename TRead, bool tagTrimming>
void _preTrim(std::vector<TRead>& readSet, const unsigned head, const unsigned tail, const unsigned min, std::vector<bool>& rem)
{
    int i = 0;
    const auto limit = (int)length(readSet);
    resize(rem, limit);

    SEQAN_OMP_PRAGMA(parallel for default(shared) private(i)schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if (length(readSet[i].seq) >(head + tail))
            {
                if (head > 0)
                {
                    if (tagTrimming)
                    {
                        std::string tempString = ":TL:";
                        append(tempString, prefix(readSet[i].seq, head));
                        insertAfterFirstToken(readSet[i].id, tempString);
                    }
                    erase(readSet[i].seq, 0, head);
                }
                if (tail > 0)
                {
                    if (tagTrimming)
                    {
                        std::string tempString = ":TR:";
                        append(tempString, suffix(readSet[i].seq, length(readSet[i].seq) - tail));
                        insertAfterFirstToken(readSet[i].id, std::move(tempString));
                    }
                    erase(readSet[i].seq, length(readSet[i].seq) - tail, length(readSet[i].seq));
                }

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



template<typename TSeqs, typename TIds, typename = std::enable_if_t < std::is_same<TIds, seqan::StringSet<seqan::CharString, seqan::Owner<seqan::Default>>>::value >>
void _preTrim(TSeqs& seqs, TIds& ids, const unsigned head, const bool tagTrimming, const unsigned tail, const unsigned min, seqan::String<bool>& rem)
{
	int i = 0;
	const auto limit = (int)length(seqs);
	resize(rem, limit);

	SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
	for (i = 0; i < limit; ++i)
	{
		if (length(seqs[i]) >(head + tail))
		{
			if (head > 0)
            {
                if (tagTrimming)
                {
                    seqan::CharString tempString = ":TL:";
                    append(tempString, prefix(seqs[i], head));
                    insertAfterFirstToken(ids[i], tempString);
                }
				erase(seqs[i], 0, head);
			}
            if (tail > 0)
            {
                if (tagTrimming)
                {
                    seqan::CharString  tempString = ":TR:";
                    append(tempString, suffix(seqs[i], length(seqs[i])-tail));
                    insertAfterFirstToken(ids[i], std::move(tempString));
                }
                erase(seqs[i], length(seqs[i]) - tail, length(seqs[i]));
            }

			// check if trimmed sequence is at least of length min
			// if not, remove it
			if (length(seqs[i]) >= min)
				rem[i] = false;
			else
				rem[i] = true;
		}
		// if the sequence was too short to be trimmed, remove it
		else
			rem[i] = true;
	}
}

// overload for single end multiplex
template <typename TRead>
void preTrim(std::vector<TRead>& reads, unsigned head, unsigned tail, unsigned min, const bool tagTrimming, GeneralStats& stats)
{
    std::vector<bool> rem;
    if(tagTrimming)
        _preTrim<TRead, true>(reads, head, tail, min, rem);
    else
        _preTrim<TRead, false>(reads, head, tail, min, rem);
    stats.removedSeqsShort += _eraseSeqs(rem, true, reads);
}


template<typename TSeqs, typename TIds, typename TDemultiplexingParams>
void preTrim(TSeqs& seqs, TIds& ids, TDemultiplexingParams&& demultiplexParams, unsigned head, const bool tagTrimming, unsigned tail, unsigned min, GeneralStats& stats)
{
	String<bool> rem;
	_preTrim(seqs, ids, head, tagTrimming, tail, min, rem);
    if(length(demultiplexParams.multiplexFile) > 0)
        stats.removedSeqsShort += _eraseSeqs(rem, true, seqs, ids, demultiplexParams.multiplex);
    else
        stats.removedSeqsShort += _eraseSeqs(rem, true, seqs, ids);
}

// overload for paired end multiplex
template<typename TSeqs, typename TIds, typename TDemultiplexingParams>
void preTrim(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, TDemultiplexingParams&& demultiplexParams, unsigned head, const bool tagTrimming, unsigned tail, unsigned min, GeneralStats& stats)
{
	String<bool> rem1, rem2;
	_preTrim(seqs, ids, head, tagTrimming, tail, min, rem1);
	_preTrim(seqsRev, idsRev, head, tagTrimming, tail, min, rem2);
	for (unsigned int i = 0; i < length(rem1); ++i)
		rem1[i] = rem1[i] | rem2[i];	// remove both strands if either is marked for removal (true = 1)
    if (length(demultiplexParams.multiplexFile) > 0)
        stats.removedSeqsShort += _eraseSeqs(rem1, true, seqs, seqsRev, ids, idsRev, demultiplexParams.multiplex);
    else
        stats.removedSeqsShort += _eraseSeqs(rem1, true, seqs, seqsRev, ids, idsRev);
}

// overload for single end 
template<typename TSeqs, typename TIds>
void preTrim(TSeqs& seqs, TIds& ids, unsigned head, const bool tagTrimming, unsigned tail, unsigned min, GeneralStats& stats)
{
	preTrim(seqs, ids, DemultiplexingParams(), head, tagTrimming, tail, min, stats);
}

// overload for paired end
template<typename TSeqs, typename TIds>
void preTrim(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, unsigned head, const bool tagTrimming, unsigned tail, unsigned min, GeneralStats& stats)
{
	preTrim(seqs, ids, seqsRev, idsRev, DemultiplexingParams(), head, tagTrimming, tail, min, stats);
}

//Trims sequences to specific length and deletes to short ones together with their IDs
template<typename TRead>
void trimTo(std::vector<TRead>& reads, const unsigned len, GeneralStats& stats)
{
    std::vector<bool> rem;
    const auto limit = length(reads);
    resize(rem, limit);
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

template<typename TSeqs,typename TIds>
void trimTo(TSeqs& seqs, TIds& ids, const unsigned len, GeneralStats& stats)
{
    StringSet<bool> rem;
    const auto limit = length(seqs);
    resize(rem, limit);
    SEQAN_OMP_PRAGMA(parallel for default(shared) schedule(static))
    for (int i = 0; i < limit; ++i)
    {
        if (length(seqs[i]) < len)
        {
            rem[i] = true;
        }
        else 
        {
            rem[i] = false;
            if (length(seqs[i]) > len)
            {
                erase(seqs[i], len, length(seqs[i]));
            }
        }    
    }
    stats.removedSeqsShort += _eraseSeqs(rem, true, seqs, ids);
}

//Overload for paired end-data
template<typename TSeqs,typename TIds>
void trimTo(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, const unsigned len, GeneralStats& stats)
{
    int i = 0;
    StringSet<bool> rem;
    int limit = length(seqs);
    resize(rem, limit);
    SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
    for (i = 0; i < limit; ++i)
    {
        if ((length(seqs[i]) < len) || (length(seqsRev[i]) < len))
        {
            rem[i] = true;
        }
        else 
        {
            rem[i] = false;
            if (length(seqs[i]) > len)
            {
                erase(seqs[i], len, length(seqs[i]));
            }
            if (length(seqsRev[i]) > len)
            {
                erase(seqsRev[i], len, length(seqsRev[i]));
            }
        }    
    }
    stats.removedSeqsShort += _eraseSeqs(rem, true, seqs, seqsRev, ids, idsRev);
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_GENERALPROCESSING_H_
