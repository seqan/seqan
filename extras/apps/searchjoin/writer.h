// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef SEQAN_EXTRAS_APPS_SEARCHJOIN_WRITER_H_
#define SEQAN_EXTRAS_APPS_SEARCHJOIN_WRITER_H_

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
    typedef Stream<FileStream<File<> > > TStream;

    TDb const       & db;
    TDbQuery /* const */  & query;
    TStream         stream;
    unsigned long   recordsCount;

    Writer(TDb /* const */ & db, TDbQuery /* const */ & query) :
        db(db),
        query(query),
        stream(),
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
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
//  Function _writeRecord()                                            [Stream]
// ----------------------------------------------------------------------------

template <typename TStream, typename TString1, typename TString2>
inline void
_writeRecord(TStream & stream, TString1 const & id1, TString2 const & id2)
{
    SEQAN_OMP_PRAGMA(critical(_searchjoin_writeRecord))
    {
        streamWriteBlock(stream, &id1[0], length(id1));
        streamWriteChar(stream, ':');
        streamWriteBlock(stream, &id2[0], length(id2));
        streamWriteChar(stream, '\n');
    }
}

// ----------------------------------------------------------------------------
//  Function _write()                                                  [Writer]
// ----------------------------------------------------------------------------

template <typename TDb, typename TDbQuery, typename TSpec>
inline void
_write(Writer<TDb, TDbQuery, TSpec> & writer, typename Size<TDb>::Type dbId, typename Size<TDbQuery>::Type queryId)
{
    SEQAN_OMP_PRAGMA(atomic)
    writer.recordsCount++;

    _writeRecord(writer.stream, writer.db.ids[dbId], writer.query.ids[queryId]);
}

template <typename TDb>
inline void
_write(Writer<TDb, TDb, Join> & writer, typename Size<TDb>::Type id1, typename Size<TDb>::Type id2)
{
    SEQAN_OMP_PRAGMA(atomic)
    writer.recordsCount += 2;

    _writeRecord(writer.stream, writer.db.ids[id1], writer.db.ids[id2]);
    _writeRecord(writer.stream, writer.db.ids[id2], writer.db.ids[id1]);
}

// ----------------------------------------------------------------------------
// Function open()                                                     [Writer]
// ----------------------------------------------------------------------------

template <typename TDb, typename TDbQuery, typename TSpec, typename TFileName>
inline bool open(Writer<TDb, TDbQuery, TSpec> & writer, TFileName const & fileName)
{
    return open(writer.stream, toCString(fileName), OPEN_RDWR | OPEN_CREATE);
}

// ----------------------------------------------------------------------------
// Function close()                                                    [Writer]
// ----------------------------------------------------------------------------

template <typename TDb, typename TDbQuery, typename TSpec>
inline bool close(Writer<TDb, TDbQuery, TSpec> & writer)
{
    return close(writer.stream);
}

template <typename TDb>
inline bool close(Writer<TDb, TDb, Join> & writer)
{
    typedef typename Size<TDb>::Type            TDbSize;
    typedef typename MakeSigned<TDbSize>::Type  TDbSSize;

    TDbSSize dbSize = (TDbSSize)length(writer.db.text);

    // Add reflexive closures to the join results.
    for (TDbSSize dbId = 0; dbId < dbSize; ++dbId)
        _writeRecord(writer.stream, writer.db.ids[dbId], writer.db.ids[dbId]);

    SEQAN_OMP_PRAGMA(atomic)
    writer.recordsCount += dbSize;

    return close(writer.stream);
}

#endif  // #ifndef SEQAN_EXTRAS_APPS_SEARCHJOIN_WRITER_H_
