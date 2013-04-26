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

#include "store.h"
#include "matches.h"

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Sorter
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TSpec = void>
struct Sorter
{
    typedef MatchStore<>    TMatchStore;

    TMatchStore         _store;
    TDelegate           & delegate;
    unsigned long       matchesCount;

    Sorter(TDelegate & delegate) :
        _store(),
        delegate(delegate),
        matchesCount(0)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open()                                                     [Sorter]
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TSpec, typename TString>
bool open(Sorter<TDelegate, TSpec> & sorter, TString const & mappedReadsFile)
{
    return open(sorter._store, mappedReadsFile);
}

// ----------------------------------------------------------------------------
// Function close()                                                    [Sorter]
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TSpec>
bool close(Sorter<TDelegate, TSpec> & sorter)
{
    return close(sorter._store);
}

// ----------------------------------------------------------------------------
// Function sort()                                                     [Sorter]
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TSpec>
void sort(Sorter<TDelegate, TSpec> & sorter, unsigned matchesPerRead)
{
//    typedef Sorter<TDelegate, TSpec>        TSorter;
//    typedef typename TSorter::TMatchStore   TMatchStore;
//    typedef typename TMatchStore::TMatches  TMatches;
    typedef String<Match<> >                TMatches;

    TMatches matches;

    while (getNext(sorter._store, matches))
    {
        removeDuplicateMatches(matches);

        if (matchesPerRead < MaxValue<unsigned>::VALUE)
            sortByErrors(matches);

        _delegateMatches(sorter, matches, matchesPerRead);
    }
}

// ----------------------------------------------------------------------------
// Function _delegateMatches()                                         [Sorter]
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TSpec, typename TRecordSpec, typename TStringSpec>
inline void
_delegateMatches(Sorter<TDelegate, TSpec> & sorter,
                 String<Match<TRecordSpec>, TStringSpec> const & matches,
                 unsigned matchesPerRead)
{
    typedef Match<TRecordSpec>                              TMatch;
    typedef String<TMatch, TStringSpec> const               TMatches;
    typedef typename Iterator<TMatches, Standard>::Type     TIterator;

    TIterator matchesIt = begin(matches, Standard());
    TIterator matchesEnd = end(matches, Standard());

    unsigned matchesCount = 0;

    for (; matchesIt != matchesEnd && matchesCount < matchesPerRead; ++matchesIt, ++matchesCount)
        onMatch(sorter.delegate, *matchesIt);

    sorter.matchesCount += matchesCount;
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_SORTER_H_
