// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// This file contains the SeqStore class.
// ==========================================================================

#ifndef APP_YARA_STORE_SEQS_H_
#define APP_YARA_STORE_SEQS_H_

#include <seqan/seq_io.h>
#include <random>

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeqConfig
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct SeqConfig
{
    typedef Dna5Q                   TAlphabet;
    typedef Alloc<>                 TSeqSpec;
    typedef Owner<ConcatDirect<> >  TSeqsSpec;
    typedef Owner<ConcatDirect<> >  TSeqNamesSpec;
};

// ----------------------------------------------------------------------------
// Class SeqStore
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = SeqConfig<TSpec> >
struct SeqStore
{
    typedef typename TConfig::TAlphabet                 TAlphabet;
    typedef typename TConfig::TSeqSpec                  TSeqSpec;
    typedef typename TConfig::TSeqsSpec                 TSeqsSpec;
    typedef typename TConfig::TSeqNamesSpec             TSeqNamesSpec;

    typedef String<TAlphabet, TSeqSpec>                 TSeq;
    typedef StringSet<TSeq, TSeqsSpec>                  TSeqs;
	typedef StringSet<CharString, TSeqNamesSpec>        TSeqNames;

    TSeqs           seqs;
    TSeqNames       names;

    SeqStore() :
        seqs(),
        names()
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clear(SeqStore<TSpec, TConfig> & me)
{
    clear(me.seqs);
    clear(me.names);
}

// ----------------------------------------------------------------------------
// Function shrinkToFit()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void shrinkToFit(SeqStore<TSpec, TConfig> & me)
{
    shrinkToFit(me.seqs);
    shrinkToFit(me.names);
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
inline void reserve(SeqStore<TSpec, TConfig> & me, TSize newCapacity)
{
    reserve(me.seqs, newCapacity, Exact());
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileName>
inline bool open(SeqStore<TSpec, TConfig> & me, TFileName const & fileName, int openMode)
{
    CharString name;

    name = fileName;    append(name, ".txt");
    if (!open(me.seqs, toCString(name), openMode)) return false;

    name = fileName;    append(name, ".rid");
    if (!open(me.names, toCString(name), openMode)) return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileName>
inline bool save(SeqStore<TSpec, TConfig> const & me, TFileName const & fileName)
{
    CharString name;

    name = fileName;    append(name, ".txt");
    if (!save(me.seqs, toCString(name))) return false;

    name = fileName;    append(name, ".rid");
    if (!save(me.names, toCString(name))) return false;

    return true;
}

// --------------------------------------------------------------------------
// Function swap()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void swap(SeqStore<TSpec, TConfig> & a, SeqStore<TSpec, TConfig> & b)
{
    std::swap(a.seqs, b.seqs);
    std::swap(a.names, b.names);
}

// ----------------------------------------------------------------------------
// Function readRecords()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void readRecords(SeqStore<TSpec, TConfig> & me,
                        FormattedFile<Fastq, Input, TFileSpec> & fileIn)
{
    readRecords(me.names, me.seqs, fileIn);
}

template <typename TSpec, typename TConfig, typename TFileSpec, typename TSize>
inline void readRecords(SeqStore<TSpec, TConfig> & me,
                        FormattedFile<Fastq, Input, TFileSpec> & fileIn,
                        TSize maxRecords)
{
    readRecords(me.names, me.seqs, fileIn, maxRecords);
}

// ----------------------------------------------------------------------------
// Function trimSeqNames()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void trimSeqNames(SeqStore<TSpec, TConfig> & me)
{
    typedef SeqStore<TSpec, TConfig>                    TSeqStore;
    typedef typename TSeqStore::TSeqNames               TSeqNames;
    typedef typename Value<TSeqNames>::Type const       TSeqName;
    typedef typename Iterator<TSeqName, Rooted>::Type   TSeqNameIt;
    typedef CharString                                  TSeqNameBuffer;

    TSeqNameBuffer name;
    TSeqNames names;
    reserve(names, length(me.names), Exact());

    for (unsigned nameId = 0; nameId < length(me.names); ++nameId)
    {
        TSeqNameIt nameIt = begin(me.names[nameId], Rooted());
        readUntil(name, nameIt, IsSpace());
        appendValue(names, name);
        clear(name);
    }

    assign(me.names, names);
}

// ----------------------------------------------------------------------------
// Function randomizeNs()
// ----------------------------------------------------------------------------

template <typename TString, typename TRng>
inline void randomizeNs(TString && str, TRng & rng)
{
    typedef typename Value<TString>::Type               TAlphabet;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    TIter it = begin(str, Standard());
    TIter itEnd = end(str, Standard());

    while (it != itEnd)
    {
        for (; it != itEnd && value(it) != TAlphabet('N'); ++it) ;

        if (it == itEnd) break;

        for (; it != itEnd && value(it) == TAlphabet('N'); ++it)
            value(it) = rng() % ValueSize<Dna>::VALUE;
    }
}

// ----------------------------------------------------------------------------
// Function randomizeNs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void randomizeNs(SeqStore<TSpec, TConfig> & me)
{
    std::mt19937 rng(0xDEADBEEF);

    for (unsigned seqId = 0; seqId < length(me.seqs); ++seqId)
        randomizeNs(me.seqs[seqId], rng);
}

// ----------------------------------------------------------------------------
// Function reverse()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void reverse(SeqStore<TSpec, TConfig> & me)
{
    for (unsigned seqId = 0; seqId < length(me.seqs); ++seqId)
        reverse(me.seqs[seqId]);
}

// ----------------------------------------------------------------------------
// Function appendReverseComplement()
// ----------------------------------------------------------------------------
// Append reverse complemented sequences.

template <typename TSpec, typename TConfig>
void appendReverseComplement(SeqStore<TSpec, TConfig> & me)
{
    typedef SeqStore<TSpec, TConfig>    TSeqStore;
    typedef typename TSeqStore::TSeqs   TSeqs;
    typedef typename Value<TSeqs>::Type TSeq;
    typedef typename Size<TSeqs>::Type  TSeqId;

    TSeqId seqsCount = length(me.seqs);

    reserve(me.seqs, 2 * seqsCount, Exact());
    reserve(concat(me.seqs), 2 * lengthSum(me.seqs), Exact());

    for (TSeqId seqId = 0; seqId < seqsCount; ++seqId)
    {
        TSeq const & seq = me.seqs[seqId];
        appendValue(me.seqs, seq);
        reverseComplement(back(me.seqs));
    }
}


#endif  // #ifndef APP_YARA_STORE_SEQS_H_
