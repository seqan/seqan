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
// This file contains the PrefetchedFile class.
// ==========================================================================

#ifndef APP_YARA_FILE_PREFETCHED_H_
#define APP_YARA_FILE_PREFETCHED_H_

namespace seqan {

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
inline bool open(PrefetchedFile<Pair<TFile>, TRecords, TThreading> & me, const char * fileName1, const char * fileName2)
{
    return open(me.file, fileName1, fileName2);
}

// ----------------------------------------------------------------------------
// Function open<Pair<TFile>, Parallel>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords>
inline bool open(PrefetchedFile<Pair<TFile>, TRecords, Parallel> & me, const char * fileName1, const char * fileName2)
{
    bool status = open(me.reader.worker, fileName1, fileName2);
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

#endif  // #ifndef APP_YARA_FILE_PREFETCHED_H_
