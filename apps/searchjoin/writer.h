// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// This file contains the Writer class.
// ==========================================================================

#ifndef SEQAN_APPS_SEARCHJOIN_WRITER_H_
#define SEQAN_APPS_SEARCHJOIN_WRITER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Search
// ----------------------------------------------------------------------------

struct Search_;
typedef Tag<Search_>            Search;

// ----------------------------------------------------------------------------
// Tag Join
// ----------------------------------------------------------------------------

struct Join_;
typedef Tag<Join_>              Join;

// ----------------------------------------------------------------------------
// Class Writer
// ----------------------------------------------------------------------------

template <typename TDb, typename TDbQuery, typename TSpec = Search>
struct Writer
{
    typedef std::ofstream                                           TOutputStream;
    typedef typename DirectionIterator<TOutputStream, Output>::Type TOutputIt;

    TDb const       & db;
    TDbQuery /* const */  & query;
    TOutputStream   outputFile;
    TOutputIt       outputIt;
    unsigned long   recordsCount;
    ReadWriteLock   writeLock;

    Writer(TDb /* const */ & db, TDbQuery /* const */ & query) :
        db(db),
        query(query),
        outputFile(),
        outputIt(),
        recordsCount(0)
    {}

    ~Writer()
    {
        close(*this);
    }

    inline void
    operator() (typename Size<TDb>::Type dbId, typename Size<TDbQuery>::Type queryId)
    {
        _write(*this, dbId, queryId);
    }
};

template <typename TDb>
struct Writer<TDb, TDb, Join> : public Writer<TDb, TDb>
{
    typedef Writer<TDb, TDb>    TBase;

    Writer(TDb /* const */ & db) :
        TBase(db, db)
    {}

    inline void
    operator() (typename Size<TDb>::Type id1, typename Size<TDb>::Type id2)
    {
        _write(*this, id1, id2);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
//  Function _writeRecord()                                            [Writer]
// ----------------------------------------------------------------------------

template <typename TDb, typename TDbQuery, typename TSpec, typename TString1, typename TString2>
inline void
_writeRecord(Writer<TDb, TDbQuery, TSpec> & writer, TString1 const & id1, TString2 const & id2)
{
    lockWriting(writer.writeLock);

    write(writer.outputIt, id1);
    writeValue(writer.outputIt, ':');
    write(writer.outputIt, id2);
    writeValue(writer.outputIt, '\n');

    unlockWriting(writer.writeLock);
}

// ----------------------------------------------------------------------------
//  Function _write()                                                  [Writer]
// ----------------------------------------------------------------------------

template <typename TDb, typename TDbQuery, typename TSpec>
inline void
_write(Writer<TDb, TDbQuery, TSpec> & writer, typename Size<TDb>::Type dbId, typename Size<TDbQuery>::Type queryId)
{
    atomicInc(writer.recordsCount);

    _writeRecord(writer, writer.db.ids[dbId], writer.query.ids[queryId]);
}

template <typename TDb>
inline void
_write(Writer<TDb, TDb, Join> & writer, typename Size<TDb>::Type id1, typename Size<TDb>::Type id2)
{
    atomicAdd(writer.recordsCount, 2);

    _writeRecord(writer, writer.db.ids[id1], writer.db.ids[id2]);
    _writeRecord(writer, writer.db.ids[id2], writer.db.ids[id1]);
}

// ----------------------------------------------------------------------------
// Function open()                                                     [Writer]
// ----------------------------------------------------------------------------

template <typename TDb, typename TDbQuery, typename TSpec, typename TFileName>
inline bool open(Writer<TDb, TDbQuery, TSpec> & writer, TFileName const & fileName)
{
    if (!open(writer.outputFile, toCString(fileName), OPEN_WRONLY | OPEN_CREATE))
        return false;

    writer.outputIt = directionIterator(writer.outputFile, Output());
    return true;
}

// ----------------------------------------------------------------------------
// Function close()                                                    [Writer]
// ----------------------------------------------------------------------------

template <typename TDb, typename TDbQuery, typename TSpec>
inline bool close(Writer<TDb, TDbQuery, TSpec> & writer)
{
    return close(writer.outputFile);
}

template <typename TDb>
inline bool close(Writer<TDb, TDb, Join> & writer)
{
    typedef typename Size<TDb>::Type            TDbSize;
    typedef typename MakeSigned<TDbSize>::Type  TDbSSize;

    TDbSSize dbSize = (TDbSSize)length(writer.db.text);

    // Add reflexive closures to the join results.
    for (TDbSSize dbId = 0; dbId < dbSize; ++dbId)
        _writeRecord(writer, writer.db.ids[dbId], writer.db.ids[dbId]);

    atomicAdd(writer.recordsCount, dbSize);

    return close(writer.outputFile);
}

#endif  // #ifndef SEQAN_APPS_SEARCHJOIN_WRITER_H_
