// ==========================================================================
//                              generalProcessing.h
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Sebastian Roskosch 
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

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.GeneralStats:
..summary:Struct holding information on general statistics.
..include:seqan/generalProcessing.h
.Memvar.GeneralStats#removedSeqs
..class:Class.GeneralStats
..summary:Number of deleted sequences due to N's
..type:nolink:unsigned
.Memvar.GeneralStats#uncalledBases
..class:Class.GeneralStats
..summary:Number of uncalled bases (possibly masked) in surviving sequences.
..type:nolink:unsigned long
.Memvar.GeneralStats#removedSeqsShort
..class:Class.GeneralStats
..summary:Number of deleted sequences due to shortness.
..type:nolink:unsigned

.Memfunc.GeneralStats#GeneralStats:
..class:Class.GeneralStats
..summary:Default Constructor
..signature:GeneralStats(void)
*/
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

/**
.Function.findN
..summary:Function for finding N's in a single sequence, substituting them and/or deleting it, if
 the number of found N's is too high.
..signature:findN(seq, allowed, substitute)
..param.seq:The sequence.
...type:nolink:String
..param.allowed:Number of allowed N's per sequence.
...type:nolink:unsigned
..param.substitute:The (optional) substitute base for N's.
...type:nolink:Dna5
...type:nolink:Dna5Q
...type:nolink:char
..remarks:Alphabet of the given sequence must be Dna5 or Dna5Q, otherwise results will be wrong.
..returns:The Number of found N's or substitutions. If too many N's were found, -1 is returned.
...type:nolink:int
..see:Function.processN
*/
template<typename TAlpha, typename TSub>
int findN(String<TAlpha>& seq, unsigned allowed, TSub substitute)
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
int findN(String<TAlpha>& seq, unsigned allowed)
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

/**
.Function.processN
..summary:Function for finding N's in all given sequences, substituting them and/or deleting them, if
 the number of found N's is too high.
..signature:processN(seqs, ids, seqsRev, idsRev, multiplex, allowed, substitute, stats)
..signature:processN(seqs, ids, seqsRev, idsRev, multiplex, allowed, stats)
..signature:processN(seqs, ids, seqsRev, idsRev, allowed, substitute, stats)
..signature:processN(seqs, ids, seqsRev, idsRev, allowed, stats)
..signature:processN(seqs, ids, multiplex, allowed, substitute, stats)
..signature:processN(seqs, ids, multiplex, allowed, stats)
..signature:processN(seqs, ids, allowed, substitute, stats)
..signature:processN(seqs, ids, allowed, stats)
..param.seqs:The forward reads.
...type:Class.StringSet
..param.ids:The ids associated with the forward reads.
...type:Class.StringSet
..param.seqsRev:The backward reads.
...type:Class.StringSet
..param.idsRev:The ids associated with the backward reads.
...type:Class.StringSet
..param.multiplex:The multiplex barcodes of the reads.
...type:Class.StringSet
..param.allowed:Number of allowed N's per sequence.
...type:nolink:unsigned
..param.substitute:The (optional) substitute base for N's.
...type:nolink:Dna5
...type:nolink:Dna5Q
...type:nolink:char
..param.stats:The GeneralStats object holding general statistic data.
...type:Class.GeneralStats
..remarks:This function simply calls findN and performs the deletion of ids/sequences/barcodes if necessary
..remarks:If the number of N's in one part of the sequence (e.g. the backward reads) is too hight, also the forward read
 will be deleted.
..returns:void
..see:Function.findN
*/
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
    SEQAN_OMP_PRAGMA(parallel for default(shared)schedule(static))
    for (int i = 0; i < limit; ++i)
    {
        res[i] = findN(seqs[i], allowed);
    }
    unsigned ex = 0;
    for (int  i = length(res) -1; i >= 0; --i)
    {
        if (res[i] == -1)
        {
            swap(seqs[i], seqs[limit - ex - 1]);
            swap(ids[i], ids[limit - ex - 1]);
            ++ex;
        }
        else
        {
            stats.uncalledBases += res[i];
        }
    }   
    if (ex != 0)
    {
        resize(seqs, limit - ex);
        resize(ids, limit - ex);
        stats.removedSeqs += ex;
    }
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
    unsigned ex = 0;
    for (int i = length(res) - 1; i >= 0 ; --i)
    {
        if (res[i] == -1)
        {
            swap(seqs[i], seqs[limit - ex - 1]);
            swap(ids[i], ids[limit - ex - 1]);
            swap(seqsRev[i], seqsRev[limit - ex - 1]);
            swap(idsRev[i], idsRev[limit - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit - ex);
        resize(ids, limit - ex);
        resize(seqsRev, limit - ex);
        resize(idsRev, limit - ex);
        stats.removedSeqs += (2 * ex);
    }
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
    unsigned ex = 0;
    for (int i = length(res) - 1; i >= 0 ; --i)
    {
        if (res[i] == -1)
        {
            swap(seqs[i], seqs[limit - ex - 1]);
            swap(ids[i], ids[limit - ex - 1]);
            swap(seqsRev[i], seqsRev[limit - ex - 1]);
            swap(idsRev[i], idsRev[limit - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit - ex);
        resize(ids, limit - ex);
        resize(seqsRev, limit - ex);
        resize(idsRev, limit - ex);
        stats.removedSeqs += (2 * ex);
    }
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
    unsigned ex = 0;
    for (int i = length(res) - 1; i >= 0 ; --i) 
    {
        if (res[i] == -1)
        {
            swap(seqs[i], seqs[limit - ex - 1]);
            swap(ids[i], ids[limit - ex - 1]);
            swap(multiplex[i], multiplex[limit - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit-ex);
        resize(ids, limit-ex);
        resize(multiplex, limit-ex);
        stats.removedSeqs += ex;
    }
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
    unsigned ex = 0;
    for (int i = length(res) - 1; i >= 0 ; --i) 
    {
        if (res[i] == -1)
        {                               
            swap(seqs[i], seqs[limit - ex - 1]);
            swap(ids[i], ids[limit - ex - 1]);
            swap(multiplex[i], multiplex[limit - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit-ex);
        resize(ids, limit-ex);
        resize(multiplex, limit-ex);
        stats.removedSeqs += ex;
    }
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
    unsigned ex = 0;
    for (int i = length(res) - 1; i >= 0 ; --i) 
    {
        if (res[i] == -1)
        {
            swap(seqs[i], seqs[limit - ex - 1]);
            swap(seqsRev[i], seqsRev[limit - ex - 1]);
            swap(ids[i], ids[limit - ex - 1]);
            swap(idsRev[i], idsRev[limit - ex - 1]);
            swap(multiplex[i], multiplex[limit - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit-ex);
        resize(seqsRev, limit-ex);
        resize(ids, limit-ex);
        resize(idsRev, limit-ex);
        resize(multiplex, limit-ex);
        stats.removedSeqs += (2 * ex);
    }
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
    unsigned ex = 0;
    for (int i = length(res) - 1; i >= 0 ; --i) 
    {
        if (res[i] == -1)
        {
            swap(seqs[i], seqs[limit - ex - 1]);
            swap(seqsRev[i], seqsRev[limit - ex - 1]);
            swap(ids[i], ids[limit - ex - 1]);
            swap(idsRev[i], idsRev[limit - ex - 1]);
            swap(multiplex[i], multiplex[limit - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit-ex);
        resize(seqsRev, limit-ex);
        resize(ids, limit-ex);
        resize(idsRev, limit-ex);
        resize(multiplex, limit-ex);
        stats.removedSeqs += (2 * ex);
    }
}


/**
.Function.preTrim
..summary:Function for trimming the sequences and deleting too short ones before further processing.
..signature:preTrim(seqs, ids, seqsRev, idsRev, multiplex, head, tail, min, stats)
..signature:preTrim(seqs, ids, seqsRev, idsRev, head, tail, min, stats)
..signature:preTrim(seqs, ids, multiplex, head, tail, min, stats)
..signature:preTrim(seqs, ids, head, tail, min, stats)
..param.seqs:The forward reads.
...type:Class.StringSet
..param.ids:The ids associated with the forward reads.
...type:Class.StringSet
..param.seqsRev:The backward reads.
...type:Class.StringSet
..param.idsRev:The ids associated with the backward reads.
...type:Class.StringSet
..param.multiplex:The multiplex barcodes of the reads.
...type:Class.StringSet
..param.head:Number of bases to be cut at the 5' end(s).
...type:nolink:unsigned
..param.tail:Number of bases to be cut at the 3' end(s).
...type:nolink:unsigned
..param.min:Minimum sequence length AFTER trimming.
...type:nolink:unsigned
..param.stats:The GeneralStats object holding general statistic data.
...type:Class.GeneralStats
..returns:void
*/
template<typename TSeqs, typename TIds> //Version for single end data
void preTrim(TSeqs& seqs, TIds& ids, unsigned head, unsigned tail, unsigned min, GeneralStats& stats)
{
    int i = 0;
    int limit = length(seqs);
    StringSet<bool> rem;

    resize(rem, limit);
    if (head > 0 && tail > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if (length(seqs[i]) > (head + tail))
            {
                erase(seqs[i], 0, head);
                erase(seqs[i], length(seqs[i]) - tail, length(seqs[i]));
                if (length(seqs[i]) >= min)
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else if (head > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
           if (length(seqs[i]) > head)
            {
                erase(seqs[i], 0, head);
                if (length(seqs[i]) >= min)
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else if (tail > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if (length(seqs[i]) > tail)
            {
                erase(seqs[i], length(seqs[i]) - tail, length(seqs[i]));
                if (length(seqs[i]) >= min)
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if (length(seqs[i]) >= min)
            {
                rem[i] = false;
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    unsigned  ex = 0;
    for (int j = length(rem) - 1; j >= 0; --j)
    {
        if (rem[j])
        {
            swap(seqs[j], seqs[limit - ex - 1]);
            swap(ids[j], ids[limit - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit - ex);
        resize(ids, limit - ex);
        stats.removedSeqsShort += ex;
    }
}

//Overload for single end data with multiplex barcodes
template<typename TSeqs, typename TIds, typename TMulti>
void preTrim(TSeqs& seqs, TIds& ids, TMulti& multiplex, unsigned head, unsigned tail, unsigned min, GeneralStats& stats)
{
    int i = 0;
    int limit = length(seqs);
    StringSet<bool> rem;

    resize(rem, limit);
    if (head > 0 && tail > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if (length(seqs[i]) > (head + tail))
            {
                erase(seqs[i], 0, head);
                erase(seqs[i], length(seqs[i]) - tail, length(seqs[i]));
                if (length(seqs[i]) >= min)
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else if (head > 0)
    {
       SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
           if (length(seqs[i]) > head)
            {
                erase(seqs[i], 0, head);
                if (length(seqs[i]) >= min)
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else if (tail > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if (length(seqs[i]) > tail)
            {
                erase(seqs[i], length(seqs[i]) - tail, length(seqs[i]));
                if (length(seqs[i]) >= min)
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if (length(seqs[i]) >= min)
            {
                rem[i] = false;
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    unsigned ex = 0;
    for (int j = length(rem) - 1; j >= 0; --j)
    {
        if (rem[j])
        {
            swap(seqs[j], seqs[limit - ex -1]);
            swap(ids[j], ids[limit - ex -1]);
            swap(multiplex[j], multiplex[limit - ex -1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit - ex);
        resize(ids, limit - ex);
        resize(multiplex, limit - ex);
        stats.removedSeqsShort += ex;
    }
}

//Overload for paired end data
template<typename TSeqs, typename TIds> 
void preTrim(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, unsigned head, unsigned tail, unsigned min,
    GeneralStats& stats)
{
    int i = 0;
    int limit = length(seqs);
    StringSet<bool> rem;
    resize(rem, limit);
    if (head > 0 && tail > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if ((length(seqs[i]) > (head + tail)) && (length(seqsRev[i]) > (head + tail)))
            {
                erase(seqs[i], 0, head);
                erase(seqs[i], length(seqs[i]) - tail, length(seqs[i]));
                erase(seqsRev[i], 0, head);
                erase(seqsRev[i], length(seqsRev[i]) - tail, length(seqsRev[i]));
                if ((length(seqs[i]) >= min) && (length(seqsRev[i]) >= min)) 
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else if (head > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if ((length(seqs[i]) > head) && (length(seqsRev[i]) > head))
            {
                erase(seqs[i], 0, head);
                erase(seqsRev[i], 0, head);
                if ((length(seqs[i]) >= min) && (length(seqsRev[i]) >= min))
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else if (tail > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if ((length(seqs[i]) > tail) && (length(seqsRev[i]) > tail))
            {
                erase(seqs[i], length(seqs[i]) - tail, length(seqs[i]));
                erase(seqsRev[i], length(seqsRev[i]) - tail, length(seqsRev[i]));
                if ((length(seqs[i]) >= min) && (length(seqsRev[i]) >= min))
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if ((length(seqs[i]) >= min) && (length(seqsRev[i]) >= min))
            {
                rem[i] = false;
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    unsigned ex = 0;
    for (int j = length(rem) - 1; j >= 0; --j)
    {
        if (rem[j])
        {
            swap(seqs[j], seqs[limit - ex - 1]);
            swap(seqsRev[j], seqsRev[limit - ex - 1]);
            swap(ids[j], ids[limit - ex - 1]);
            swap(idsRev[j], idsRev[limit - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit - ex);
        resize(seqsRev, limit - ex);
        resize(ids, limit - ex);
        resize(idsRev, limit - ex);
        stats.removedSeqsShort += (2 * ex);
    }

}

//Overload for paired end data with multiplex barcodes
template<typename TSeqs, typename TIds,  typename TMulti> 
void preTrim(TSeqs& seqs, TIds& ids, TSeqs& seqsRev, TIds& idsRev, TMulti& multiplex, unsigned head, unsigned tail, unsigned min,
    GeneralStats& stats)
{
    int i = 0;
    int limit = length(seqs);
    StringSet<bool> rem;
    resize(rem, limit);
    if (head > 0 && tail > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if ((length(seqs[i]) > (head + tail)) && (length(seqsRev[i]) > (head + tail)))
            {
                erase(seqs[i], 0, head);
                erase(seqs[i], length(seqs[i]) - tail, length(seqs[i]));
                erase(seqsRev[i], 0, head);
                erase(seqsRev[i], length(seqsRev[i]) - tail, length(seqsRev[i]));
                if ((length(seqs[i]) >= min) && (length(seqsRev[i]) >= min)) 
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else if (head > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if ((length(seqs[i]) > head) && (length(seqsRev[i]) > head))
            {
                erase(seqs[i], 0, head);
                erase(seqsRev[i], 0, head);
                if ((length(seqs[i]) >= min) && (length(seqsRev[i]) >= min))
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else if (tail > 0)
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if ((length(seqs[i]) > tail) && (length(seqsRev[i]) > tail))
            {
                erase(seqs[i], length(seqs[i]) - tail, length(seqs[i]));
                erase(seqsRev[i], length(seqsRev[i]) - tail, length(seqsRev[i]));
                if ((length(seqs[i]) >= min) && (length(seqsRev[i]) >= min))
                {
                    rem[i] = false;
                }
                else
                {
                    rem[i] = true;
                }
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    else
    {
        SEQAN_OMP_PRAGMA(parallel for default(shared) private(i) schedule(static))
        for (i = 0; i < limit; ++i)
        {
            if ((length(seqs[i]) >= min) && (length(seqsRev[i]) >= min))
            {
                rem[i] = false;
            }
            else
            {
                rem[i] = true;
            }
        }
    }
    unsigned ex = 0;
    for (int j = length(rem) - 1; j >= 0; --j)
    {
        if (rem[j])
        {
            swap(seqs[j], seqs[j - ex - 1]);
            swap(seqsRev[j], seqsRev[j - ex - 1]);
            swap(ids[j], ids[j - ex - 1]);
            swap(idsRev[j], idsRev[j - ex - 1]);
            swap(multiplex[j], multiplex[j - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit - ex);
        resize(seqsRev, limit - ex);
        resize(ids, limit - ex);
        resize(idsRev, limit - ex);
        resize(multiplex, limit - ex);
        stats.removedSeqsShort += (2 * ex);
    }
}

//Trims sequences to specific length and deletes to short ones together with their IDs
template<typename TSeqs,typename TIds>
void trimTo(TSeqs& seqs, TIds& ids, const unsigned len, GeneralStats& stats)
{
    StringSet<bool> rem;
    int limit = length(seqs);
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
    unsigned ex = 0;
    for (int j = limit - 1; j >= 0; --j)
    {
        if (rem[j])
        {
            swap(seqs[j], seqs[limit  - ex - 1]);
            swap(ids[j], ids[limit  - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit - ex);
        resize(ids, limit - ex);
        stats.removedSeqsShort += ex;
    }
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
    unsigned ex = 0;
    for (int j = length(rem) - 1; j >= 0; --j)
    {
        if (rem[j])
        {
            swap(seqs[j], seqs[limit - ex - 1]);            
            swap(seqsRev[j], seqsRev[limit - ex - 1]);
            swap(ids[j], ids[limit - ex - 1]);
            swap(idsRev[j], idsRev[limit - ex - 1]);
            ++ex;
        }
    }
    if (ex != 0)
    {
        resize(seqs, limit - ex);
        resize(seqsRev, limit - ex);
        resize(ids, limit - ex);
        resize(idsRev, limit - ex);
        stats.removedSeqsShort += (2 * ex);
    }
}

#endif  // #ifndef SANDBOX_GROUP3_APPS_SEQDPT_GENERALPROCESSING_H_
