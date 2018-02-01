// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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
    uint64_t    maxRecords;

    PrefetchedFile(uint64_t maxRecords) :
        file(),
        maxRecords(maxRecords)
    {}
};

// ----------------------------------------------------------------------------
// Class PrefetchedFile<Parallel>
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords>
struct PrefetchedFile<TFile, TRecords, Parallel>
{
    TFile           file;
    TRecords        records;
    uint64_t        maxRecords;
    std::thread     reader;

    PrefetchedFile(uint64_t maxRecords) :
        file(),
        records(),
        maxRecords(maxRecords),
        reader()
    {}

    ~PrefetchedFile()
    {
        close(*this);
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
// Function open<Pair<TFile>, Serial>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords, typename TThreading>
inline bool open(PrefetchedFile<Pair<TFile>, TRecords, TThreading> & me, const char * fileName1, const char * fileName2)
{
    return open(me.file, fileName1, fileName2);
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
// Function readRecords<Serial>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords, typename TThreading>
inline void readRecords(TRecords & records, PrefetchedFile<TFile, TRecords, TThreading> & me)
{
    readRecords(records, me.file, me.maxRecords);
}

// ----------------------------------------------------------------------------
// Function _prefetchRecords<Parallel>()
// ----------------------------------------------------------------------------
// Prefetches the first batch of records.

template <typename TFile, typename TRecords>
inline void _prefetchRecords(PrefetchedFile<TFile, TRecords, Parallel> & me)
{
    me.reader = std::thread([&me]() { readRecords(me.records, me.file, me.maxRecords); });
}

// ----------------------------------------------------------------------------
// Function open<Parallel>()
// ----------------------------------------------------------------------------
// Prefetches the first batch of records.

template <typename TFile, typename TRecords>
inline bool open(PrefetchedFile<TFile, TRecords, Parallel> & me, const char * fileName)
{
    if (open(me.file, fileName))
    {
        _prefetchRecords(me);
        return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function open<Pair<TFile>, Parallel>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords>
inline bool open(PrefetchedFile<Pair<TFile>, TRecords, Parallel> & me, const char * fileName1, const char * fileName2)
{
    if (open(me.file, fileName1, fileName2))
    {
        _prefetchRecords(me);
        return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function close<Parallel>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords>
inline void close(PrefetchedFile<TFile, TRecords, Parallel> & me)
{
    if (me.reader.joinable())
        me.reader.join();

    close(me.file);
}

// ----------------------------------------------------------------------------
// Function readRecords<Parallel>()
// ----------------------------------------------------------------------------

template <typename TFile, typename TRecords>
inline void readRecords(TRecords & records, PrefetchedFile<TFile, TRecords, Parallel> & me)
{
    // Wait the current batch of records.
    if (me.reader.joinable())
        me.reader.join();

    // Return the current batch of records.
    swap(records, me.records);

    // Read the next batch of records.
    _prefetchRecords(me);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

}

#endif  // #ifndef APP_YARA_FILE_PREFETCHED_H_
