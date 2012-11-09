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
// This file contains the Sorter class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_SORTER_H_
#define SEQAN_EXTRAS_MASAI_SORTER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#include "tags.h"
#include "store.h"
#include "indexer.h"
#include "matches.h"
#include "writer.h"
#include "stream.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Sorter
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Sorter
{
    typedef Indexer<Nothing>    TIndexer;

    TFragmentStore      store;
    TIndexer            indexer;

    unsigned            readsCount;
    unsigned            matchesCount;

    bool                writeCigar;
    bool                dumpResults;

    unsigned            matchesPerRead;

    // TODO(esiragusa): Remove writeCigar from Sorter members.
    Sorter(unsigned matchesPerRead, bool writeCigar = true, bool dumpResults = true) :
        indexer(store),
        readsCount(0),
        matchesCount(0),
        writeCigar(writeCigar),
        dumpResults(dumpResults),
        matchesPerRead(matchesPerRead)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function loadReads()                                                [Sorter]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TString>
bool loadReads(Sorter<TSpec> & sorter, TString const & readsFile)
{
    if (!loadReads(sorter.store, readsFile))
        return false;

    sorter.readsCount = length(sorter.store.readSeqStore);

    _loadReadsRC(sorter);

    return true;
}

template <typename TSpec>
bool _loadReadsRC(Sorter<TSpec> & sorter)
{
    reserve(sorter.store.readSeqStore, sorter.readsCount * 2, Exact());

    for (TReadSeqStoreSize readId = 0; readId < sorter.readsCount; ++readId)
    {
        TReadSeq const & read = sorter.store.readSeqStore[readId];
        appendValue(sorter.store.readSeqStore, read);
        reverseComplement(back(sorter.store.readSeqStore));
    }

    return true;
}

// ----------------------------------------------------------------------------
// Function sortMappedReads()                                          [Sorter]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TString, typename TDistance>
bool sortMappedReads(Sorter<TSpec> & sorter,
                     TString const & mappedReadsFile,
                     TString const & sortedReadsFile,
                     TDistance const & /*tag*/,
                     Raw const & /*tag*/)
{
    typedef Stream<FileStream<Match<>, MMapWriter> >    TWriterStream;
    typedef MatchWriter<TWriterStream, TDistance, Raw>  TMatchWriter;

    TWriterStream file;

    if (sorter.dumpResults)
        if (!open(file, toCString(sortedReadsFile), OPEN_RDWR | OPEN_CREATE))
            return false;

    TMatchWriter writer(file, sorter.store, sorter.readsCount, sorter.dumpResults);

    _sortMappedReads(sorter, mappedReadsFile, writer);

    return true;
}

template <typename TSpec, typename TString, typename TDistance>
bool sortMappedReads(Sorter<TSpec> & sorter,
                     TString const & mappedReadsFile,
                     TString const & sortedReadsFile,
                     TDistance const & /*tag*/,
                     Sam const & /*tag*/)
{
    typedef Stream<FileStream<char, MMapWriter> >       TWriterStream;
    typedef MatchWriter<TWriterStream, TDistance, Sam>  TMatchWriter;

    TWriterStream file;

    if (sorter.dumpResults)
        if (!open(file, toCString(sortedReadsFile), OPEN_RDWR | OPEN_CREATE))
            return false;

    TMatchWriter writer(file, sorter.store, sorter.readsCount, sorter.dumpResults);

    // TODO(esiragusa): Remove writeCigar from members.
    writer.writeCigar = sorter.writeCigar;

    _sortMappedReads(sorter, mappedReadsFile, writer);

    return true;
}

// ============================================================================

template <typename TSpec, typename TString, typename TMatchesDelegate>
bool _sortMappedReads(Sorter<TSpec> & sorter,
                      TString const & mappedReadsFile,
                      TMatchesDelegate & matchesDelegate)
{
    typedef Match<>                  TMatch;
    typedef MatchStore<TMatch>       TMatchStore;
    typedef String<TMatch>           TMatches;

    TMatchStore store;

    if (!open(store, mappedReadsFile))
        return false;

    TMatches matches;

    if (!getNext(store, matches))
        return false;

    do
    {
        removeDuplicateMatches(matches);

        if (sorter.matchesPerRead < MaxValue<unsigned>::VALUE)
            sortByErrors(matches);

        _processMatch(sorter, matches, matchesDelegate);
    }
    while (getNext(store, matches));

    close(store);

    return true;
}

template <typename TSpec, typename TRecordSpec, typename TStringSpec, typename TMatchesDelegate>
inline void _processMatch(Sorter<TSpec> & sorter,
                          String<Match<TRecordSpec>, TStringSpec> const & matches,
                          TMatchesDelegate & matchesDelegate)
{
    typedef Match<TRecordSpec>                              TMatch;
    typedef String<TMatch, TStringSpec>                     TMatches;
    typedef typename Iterator<TMatches, Standard>::Type     TIterator;

    TIterator matchesIt = begin(matches, Standard());
    TIterator matchesEnd = end(matches, Standard());

    unsigned matchesCount = 0;

    for (; matchesIt != matchesEnd && matchesCount < sorter.matchesPerRead; ++matchesIt, ++matchesCount)
        onMatch(matchesDelegate, *matchesIt);

    sorter.matchesCount += matchesCount;
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_SORTER_H_
