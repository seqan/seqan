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
// This file contains classes for storing and manipulating matches.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_MATCHES_H_
#define SEQAN_EXTRAS_MASAI_MATCHES_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/pipe.h>

#include "store.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Match
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Match
{
    unsigned        readId;
    unsigned        beginPos;
    short           endPosDelta;
    unsigned char   contigId;
    unsigned char   errors;
};

// ----------------------------------------------------------------------------
// Class MatchSorterByXXX
// ----------------------------------------------------------------------------

template <typename TMatch>
struct MatchSorterByReadId :
    std::binary_function<TMatch, TMatch, int>
{
    inline int operator()(TMatch const & a, TMatch const & b) const
    {
        return a.readId - b.readId;
    }

};

template <typename TMatch>
struct MatchSorterByBeginPos
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return (a.contigId < b.contigId) ||
               (a.contigId == b.contigId && (isForward(a) && isReverse(b))) ||
               (a.contigId == b.contigId && !(isReverse(a) && isForward(b)) && a.beginPos < b.beginPos);
    }

};

template <typename TMatch>
struct MatchSorterByEndPos
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return (a.contigId < b.contigId) ||
               (a.contigId == b.contigId && (isForward(a) && isReverse(b))) ||
               (a.contigId == b.contigId && !(isReverse(a) && isForward(b)) && endPos(a) < endPos(b));
    }

};

template <typename TMatch>
struct MatchSorterByErrors
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return a.errors < b.errors;
    }

};

// ----------------------------------------------------------------------------
// Class MatchIterator
// ----------------------------------------------------------------------------

template <typename TMatch, typename TMatchString, typename TSpec = void>
struct MatchIterator
{
    typedef String<TMatch, TMatchString> const              TMatches;
    typedef typename Iterator<TMatches, Standard>::Type     TIterator;

    TIterator it;
    TIterator const endIt;

    MatchIterator(String<TMatch, TMatchString> const & matches) :
        it(begin(matches)),
        endIt(end(matches))
    {}
};

// ----------------------------------------------------------------------------
// Class MatchStore
// ----------------------------------------------------------------------------

template <typename TMatchStoreString = External<>, typename TSpec = void, typename TMatch = Match<TSpec> >
struct MatchStore
{
    typedef String<TMatch, TMatchStoreString>       TMatches;
    typedef MatchSorterByReadId<TMatch>             TMatchSorter;
    typedef SorterSpec<SorterConfig<TMatchSorter> > TSorterSpec;
    typedef Pool<TMatch, TSorterSpec>               TSorterPool;

    TMatches    matches;
    TSorterPool sorterPool;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function fill()                                                      [Match]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void fill(Match<TSpec> & match,
                 TContigId contigId,
                 TContigPos beginPos,
                 TContigPos endPos,
                 TReadId readId,
                 TErrors errors,
                 bool reverseComplemented)
{
    match.readId = readId;
    match.beginPos = beginPos;
    match.endPosDelta = endPos - beginPos;
    if (reverseComplemented)
        match.endPosDelta = -match.endPosDelta;
    match.contigId = contigId;
    match.errors = errors;
}

// ----------------------------------------------------------------------------
// Function assign()                                                    [Match]
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void assign(Match<TSpec> & dest, Match<TSpec> const & source)
{
    dest.readId = source.readId;
    dest.beginPos = source.beginPos;
    dest.endPosDelta = source.endPosDelta;
    dest.contigId = source.contigId;
    dest.errors = source.errors;
}

// ----------------------------------------------------------------------------
// Function isForward()                                                 [Match]
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isForward(Match<TSpec> const & match)
{
    return match.endPosDelta > 0;
}

// ----------------------------------------------------------------------------
// Function isReverse()                                                 [Match]
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isReverse(Match<TSpec> const & match)
{
    return match.endPosDelta < 0;
}

// ----------------------------------------------------------------------------
// Function isConcordant()                                              [Match]
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isConcordant(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return (isForward(a) && isForward(b)) || (isReverse(a) && isReverse(b));
}

// ----------------------------------------------------------------------------
// Function endPos()                                                    [Match]
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned endPos(Match<TSpec> const & match)
{
    return match.beginPos + abs(match.endPosDelta);
}

// ----------------------------------------------------------------------------
// Function onMatch()                                        [TMatchesDelegate]
// ----------------------------------------------------------------------------

// NOTE(esiragusa): Syntactical sugar.
template <typename TMatchesDelegate, typename TSpec>
inline void onMatch(TMatchesDelegate & matchesDelegate, Match<TSpec> const & match)
{
    onMatch(matchesDelegate, match.contigId,
            match.beginPos, endPos(match), match.readId, match.errors, isReverse(match));
}

// NOTE(esiragusa): Syntactical sugar.
template <typename TMatchesDelegate, typename TSpec>
inline void onMatch(TMatchesDelegate & matchesDelegate, Match<TSpec> const & matchFwd, Match<TSpec> const & matchRev)
{
    onMatch(matchesDelegate, matchFwd.contigId,
            matchFwd.beginPos, endPos(matchFwd), matchFwd.readId, matchFwd.errors,
            matchRev.beginPos, endPos(matchRev), matchRev.readId, matchRev.errors);
}

// ----------------------------------------------------------------------------
// Function open()                                                 [MatchStore]
// ----------------------------------------------------------------------------

template <typename TMatchStoreString, typename TSpec, typename TMatch, typename TString>
inline bool open(MatchStore<TMatchStoreString, TSpec, TMatch> & store, TString const & file)
{
    if (!open(store.matches, toCString(file), OPEN_RDONLY))
        return false;

    store.sorterPool << store.matches;

    beginRead(store.sorterPool);

    return true;
}

// ----------------------------------------------------------------------------
// Function close()                                                [MatchStore]
// ----------------------------------------------------------------------------

template <typename TMatchStoreString, typename TSpec, typename TMatch>
inline bool close(MatchStore<TMatchStoreString, TSpec, TMatch> & store)
{
    endRead(store.sorterPool);

    return close(store.matches);
}

// ----------------------------------------------------------------------------
// Function getNext()                                              [MatchStore]
// ----------------------------------------------------------------------------

template <typename TMatchStoreString, typename TSpec, typename TMatch, typename TMatches>
inline bool getNext(MatchStore<TMatchStoreString, TSpec, TMatch> & store, TMatches & matches)
{
    TMatch match;

    clear(matches);

    if (eof(store.sorterPool))
        return false;

    do
    {
        pop(store.sorterPool, match);
        append(matches, match);
    }
    while (!eof(store.sorterPool) && front(store.sorterPool).readId == match.readId);

    return true;
}

// ----------------------------------------------------------------------------
// Function removeDuplicateMatches()
// ----------------------------------------------------------------------------

template <typename TRecordSpec, typename TStringSpec>
inline void removeDuplicateMatches(String<Match<TRecordSpec>, TStringSpec> & matches)
{
    typedef String<Match<TRecordSpec>, TStringSpec>             TMatches;
    typedef typename Iterator<TMatches, Standard>::Type         TMatchesIterator;

    TMatchesIterator matchesBegin;
    TMatchesIterator matchesEnd;
    TMatchesIterator newIt;
    TMatchesIterator oldIt;


    matchesBegin = begin(matches, Standard());
    matchesEnd   = end(matches, Standard());
    newIt = matchesBegin;
    oldIt = matchesBegin;

    std::sort(matchesBegin, matchesEnd, MatchSorterByEndPos<Match<TRecordSpec> >());

    while (oldIt != matchesEnd)
    {
        *newIt = *oldIt;

        ++oldIt;

        while (oldIt != matchesEnd &&
               (*newIt).contigId == (*oldIt).contigId &&
               isConcordant(*newIt, *oldIt) &&
               endPos(*newIt) == endPos(*oldIt))
        {
            ++oldIt;
        }

        ++newIt;
    }

    resize(matches, newIt - matchesBegin, Exact());


    matchesBegin = begin(matches, Standard());
    matchesEnd   = end(matches, Standard());
    newIt = matchesBegin;
    oldIt = matchesBegin;

    std::stable_sort(matchesBegin, matchesEnd, MatchSorterByBeginPos<Match<TRecordSpec> >());

    while (oldIt != matchesEnd)
    {
        *newIt = *oldIt;

        ++oldIt;

        while (oldIt != matchesEnd &&
               (*newIt).contigId == (*oldIt).contigId &&
               isConcordant(*newIt, *oldIt) &&
               (*newIt).beginPos == (*oldIt).beginPos)
        {
            ++oldIt;
        }

        ++newIt;
    }

    resize(matches, newIt - matchesBegin, Exact());
}

// ----------------------------------------------------------------------------
// Function sortByErrors()
// ----------------------------------------------------------------------------

template <typename TRecordSpec, typename TStringSpec>
inline void sortByErrors(String<Match<TRecordSpec>, TStringSpec> & matches)
{
    typedef String<Match<TRecordSpec>, TStringSpec>             TMatches;
    typedef typename Iterator<TMatches, Standard>::Type         TMatchesIterator;

    TMatchesIterator matchesBegin = begin(matches, Standard());
    TMatchesIterator matchesEnd = end(matches, Standard());

    std::sort(matchesBegin, matchesEnd, MatchSorterByErrors<Match<TRecordSpec> >());
}

// ----------------------------------------------------------------------------
// Function atEnd()                                             [MatchIterator]
// ----------------------------------------------------------------------------

template <typename TMatch, typename TMatchString, typename TSpec>
inline bool atEnd(MatchIterator<TMatch, TMatchString, TSpec> & matchIt)
{
    return matchIt.it == matchIt.endIt;
}

// ----------------------------------------------------------------------------
// Function getNext()                                           [MatchIterator]
// ----------------------------------------------------------------------------

template <typename TMatch, typename TMatchString, typename TSpec>
inline bool getNext(MatchIterator<TMatch, TMatchString, TSpec> & matchIt,
                    typename MatchIterator<TMatch, TMatchString, TSpec>::TIterator & matchesBegin,
                    typename MatchIterator<TMatch, TMatchString, TSpec>::TIterator & matchesEnd)
{
    if (atEnd(matchIt))
        return false;

    matchesBegin = matchIt.it;

    do
    {
        ++(matchIt.it);
    }
    while (!atEnd(matchIt) &&
           (*(matchIt.it)).contigId == (*(matchesBegin)).contigId &&
           isConcordant(*(matchIt.it), *matchesBegin));

    matchesEnd = matchIt.it;

    return true;
}

// ----------------------------------------------------------------------------
// Function getNextContig()                                     [MatchIterator]
// ----------------------------------------------------------------------------

template <typename TMatch, typename TMatchString, typename TSpec>
inline bool getNextContig(MatchIterator<TMatch, TMatchString, TSpec> & matchIt,
                          typename MatchIterator<TMatch, TMatchString, TSpec>::TIterator & matchesFwdBegin,
                          typename MatchIterator<TMatch, TMatchString, TSpec>::TIterator & matchesFwdEnd,
                          typename MatchIterator<TMatch, TMatchString, TSpec>::TIterator & matchesRevBegin,
                          typename MatchIterator<TMatch, TMatchString, TSpec>::TIterator & matchesRevEnd)
{
    typename MatchIterator<TMatch, TMatchString, TSpec>::TIterator matchesBegin;
    typename MatchIterator<TMatch, TMatchString, TSpec>::TIterator matchesEnd;

    if (getNext(matchIt, matchesBegin, matchesEnd))
    {
        if (isReverse(*matchesBegin))
        {
            matchesFwdBegin = matchesEnd;
            matchesFwdEnd = matchesEnd;

            matchesRevBegin = matchesBegin;
            matchesRevEnd = matchesEnd;
        }
        else
        {
            matchesFwdBegin = matchesBegin;
            matchesFwdEnd = matchesEnd;

            if (!atEnd(matchIt) && (*matchesBegin).contigId == (*matchesEnd).contigId)
            {
                getNext(matchIt, matchesBegin, matchesEnd);

                matchesRevBegin = matchesBegin;
                matchesRevEnd = matchesEnd;
            }
            else
            {
                matchesRevBegin = matchesEnd;
                matchesRevEnd = matchesEnd;
            }
        }

        return true;
    }

    return false;
}

// ============================================================================

// NOTE(esiragusa): Debug stuff.
template <typename TRecordSpec>
void printMatch(Match<TRecordSpec> const & match)
{
    std::cout << match.readId << " ";
    std::cout << (unsigned)match.contigId << " ";
    std::cout << isForward(match) << " ";
    std::cout << match.beginPos << " ";
    std::cout << endPos(match) << " ";
    std::cout << (unsigned)match.errors << std::endl;
}

//template <typename TRecordSpec>
//void printPair(Match<TRecordSpec> const & leftMate, Match<TRecordSpec> const & rightMate)
//{
//    std::cout << "Pair" << std::endl;
//    printMatch(leftMate);
//    printMatch(rightMate);
//    std::cout << "==========" << std::endl;
//}
//
//template <typename TRecordSpec, typename TStringSpec>
//inline void
//printMatches(String<Match<TRecordSpec>, TStringSpec> const &,
//             typename Iterator<String<Match<TRecordSpec>, TStringSpec>, Standard>::Type matchesBegin,
//             typename Iterator<String<Match<TRecordSpec>, TStringSpec>, Standard>::Type matchesEnd)
//{
//    for (; matchesBegin != matchesEnd; ++matchesBegin)
//        printMatch(*matchesBegin);
//}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_MATCHES_H_
