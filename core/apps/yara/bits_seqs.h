// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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

#ifndef APP_YARA_BITS_SEQS_H_
#define APP_YARA_BITS_SEQS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/random.h>


namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline bool open(Pair<SmartFile<TFileType, TDirection, TSpec> > & me, Pair<const char *> const & fileName)
{
    return open(me.i1, fileName.i1) && open(me.i2, fileName.i2);
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline void close(Pair<SmartFile<TFileType, TDirection, TSpec> > & me)
{
    close(me.i1);
    close(me.i2);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TFileType, typename TDirection, typename TSpec>
inline bool atEnd(Pair<SmartFile<TFileType, TDirection, TSpec> > const & me)
{
    return atEnd(me.i1) || atEnd(me.i2);
}

// ----------------------------------------------------------------------------
// Function readRecords()
// ----------------------------------------------------------------------------

template <typename TRecords, typename TFileType, typename TDirection, typename TSpec>
inline void readRecords(TRecords & records,
                        Pair<SmartFile<TFileType, TDirection, TSpec> > & me,
                        __uint64 maxRecords = MaxValue<__uint64>::VALUE)
{
    readRecords(records, me.i1, maxRecords);
    readRecords(records, me.i2, maxRecords);
}

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class PrefetchedFile<Serial>
// ----------------------------------------------------------------------------
// Serial implies no prefetching.

template <typename TFile, typename TRecords, typename TThreading = Serial>
struct PrefetchedFile
{
    TFile       file;
    TRecords    records;
    __uint64    maxRecords;

    PrefetchedFile() :
        file(),
        records(),
        maxRecords()
    {}

    PrefetchedFile(__uint64 maxRecords) :
        file(),
        records(),
        maxRecords(maxRecords)
    {}

    void operator()()
    {
        readRecords(_toParameter(records), file, maxRecords);
    }
};

// ----------------------------------------------------------------------------
// Class PrefetchedFile<Parallel>
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords>
struct PrefetchedFile<TFile, TRecords, Parallel>
{
    typedef PrefetchedFile<TFile, TRecords *, Serial> TWorker;

    Pair<TRecords>      recordsQueue;
    TRecords *          records;
    TWorker             _worker;
    Thread<TWorker>     reader;

    PrefetchedFile(__uint64 maxRecords) :
        records(&recordsQueue.i1),
        _worker(maxRecords),
        reader(_worker)
    {
        reader.worker.records = &recordsQueue.i2;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open<Serial>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords, typename TThreading>
inline bool open(PrefetchedFile<TFile, TRecords, TThreading> & me, const char * fileName)
{
    return open(me.file, fileName);
}

// ----------------------------------------------------------------------------
// Function open<Parallel>()
// ----------------------------------------------------------------------------
// Prefetches the first batch of records.

template <typename TFile, typename TRecords>
inline bool open(PrefetchedFile<TFile, TRecords, Parallel> & me, const char * fileName)
{
    bool status = open(me.reader.worker, fileName);
    run(me.reader);
    return status;
}

// ----------------------------------------------------------------------------
// Function open<Pair<TFile>, Serial>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords, typename TThreading>
inline bool open(PrefetchedFile<Pair<TFile>, TRecords, TThreading> & me, Pair<const char *> const & fileName)
{
    return open(me.file, fileName);
}

// ----------------------------------------------------------------------------
// Function open<Pair<TFile>, Parallel>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords>
inline bool open(PrefetchedFile<Pair<TFile>, TRecords, Parallel> & me, Pair<const char *> const & fileName)
{
    bool status = open(me.reader.worker, fileName);
    run(me.reader);
    return status;
}

// ----------------------------------------------------------------------------
// Function close<Serial>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords, typename TThreading>
inline void close(PrefetchedFile<TFile, TRecords, TThreading> & me)
{
    close(me.file);
}

// ----------------------------------------------------------------------------
// Function close<Parallel>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords>
inline void close(PrefetchedFile<TFile, TRecords, Parallel> & me)
{
    close(me.reader);
    close(me.reader.worker);
}

// ----------------------------------------------------------------------------
// Function readRecords<Serial>()
// ----------------------------------------------------------------------------

//template <typename TFile, typename TRecords, typename TThreading>
//inline void readRecords(TRecords & records, PrefetchedFile<TFile, TRecords, TThreading> & me)
//{
//    readRecords(records, me.file, me.maxRecords);
//}

// ----------------------------------------------------------------------------
// Function readRecords<Parallel>()
// ----------------------------------------------------------------------------

//template <typename TFile, typename TRecords>
//inline void readRecords(TRecords & records, PrefetchedFile<TFile, TRecords, Parallel> & me)
//{
//    readRecords(me.records, me.file, me.maxRecords);
//    std::swap(records, me.records);
//}

// ----------------------------------------------------------------------------
// Function getRecords<Serial>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords, typename TThreading>
inline TRecords * getRecords(PrefetchedFile<TFile, TRecords, TThreading> & me)
{
    readRecords(me.records, me.file, me.maxRecords);

    return &me.records;
}

// ----------------------------------------------------------------------------
// Function getRecords<Parallel>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords>
inline TRecords * getRecords(PrefetchedFile<TFile, TRecords, Parallel> & me)
{
    // Wait the next batch of records.
    waitFor(me.reader);

    // Make the next batch of records.
    std::swap(me.records, me.reader.worker.records);

    // Read the next batch of records.
    run(me.reader);

    return me.records;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

}


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

template <>
struct SeqConfig<Nothing>
{
    typedef Dna5                            TAlphabet;
    typedef Packed<Alloc<> >                TSeqSpec;
    typedef Owner<ConcatDirect<> >          TSeqsSpec;
    typedef Owner<ConcatDirect<__uint32> >  TSeqNamesSpec;
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
	typedef NameStoreCache<TSeqNames, CharString>       TSeqNamesCache;

    TSeqs           seqs;
    TSeqNames       names;
    TSeqNamesCache  namesCache;

    SeqStore() :
        seqs(),
        names(),
		namesCache(names)
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
//    clear(me.namesCache);
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
inline bool open(SeqStore<TSpec, TConfig> & me, TFileName const & fileName)
{
    CharString name;

    name = fileName;    append(name, ".txt");
    if (!open(me.seqs, toCString(name))) return false;

    name = fileName;    append(name, ".rid");
    if (!open(me.names, toCString(name))) return false;

//    refresh(me.namesCache);

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

// ----------------------------------------------------------------------------
// Function readRecords()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileSpec>
inline void readRecords(SeqStore<TSpec, TConfig> & me,
                        SmartFile<Fastq, Input, TFileSpec> & fileIn,
                        __uint64 maxRecords = MaxValue<__uint64>::VALUE)
{
    readRecords(me.names, me.seqs, fileIn, maxRecords);

//    appendValue(reads.names, prefix(seqName, lastOf(seqName, IsSpace())), Generous());
//    refresh(me.namesCache);
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
            value(it) = pickRandomNumber(rng) % ValueSize<Dna>::VALUE;
    }
}

// ----------------------------------------------------------------------------
// Function randomizeNs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void randomizeNs(SeqStore<TSpec, TConfig> & me)
{
    Rng<MersenneTwister> rng(0xDEADBEEF);

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


#endif  // #ifndef APP_YARA_BITS_READS_H_
