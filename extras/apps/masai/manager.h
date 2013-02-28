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
// This file contains the MatchManager class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_MANAGER_H_
#define SEQAN_EXTRAS_MASAI_MANAGER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find.h>

#include "tags.h"
#include "store.h"
#include "matches.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class MatchManager
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TSpec = void, typename TMatch = Match<TSpec> >
struct MatchManager
{
    TDelegate               & delegate;
    TReadSeqStoreSize       readsCount;
    unsigned long           matchesCount;

    MatchManager(TDelegate & delegate, TReadSeqStoreSize readsCount) :
        delegate(delegate),
        readsCount(readsCount),
        matchesCount(0)
    {}
};

template <typename TDelegate, typename TMatch>
struct MatchManager<TDelegate, AllBest, TMatch> :
    public MatchManager<TDelegate, void, TMatch>
{
    typedef MatchManager<TDelegate, void, TMatch>       TBase;

    unsigned                    errors;
    String<unsigned char>       minErrors;

    MatchManager(TDelegate & delegate, TReadSeqStoreSize readsCount) :
        TBase(delegate, readsCount),
        errors(0)
    {
        resize(minErrors, readsCount, MaxValue<unsigned char>::VALUE, Exact());
    }
};

template <typename TDelegate, typename TMatch>
struct MatchManager<TDelegate, KBest, TMatch>:
    public MatchManager<TDelegate, void, TMatch>
{
    typedef MatchManager<TDelegate, void, TMatch>       TBase;
    typedef String<TMatch>                              TMatches;

    unsigned                    errors;
    String<unsigned char>       minErrors;
    String<TMatches>            matches;

    MatchManager(TDelegate & delegate, TReadSeqStoreSize readsCount) :
        TBase(delegate, readsCount),
        errors(0)
    {
        resize(minErrors, readsCount, MaxValue<unsigned char>::VALUE, Exact());
        // TODO Change hardcoded size.
        resize(matches, 16, Exact());
    }
};

template <typename TDelegate, typename TMatch>
struct MatchManager<TDelegate, AnyBest, TMatch>:
    public MatchManager<TDelegate, void, TMatch>
{
    typedef MatchManager<TDelegate, void, TMatch>       TBase;
    typedef String<TMatch>                              TMatches;

    unsigned                    errors;
    TMatches                    matches;

    MatchManager(TDelegate & delegate, TReadSeqStoreSize readsCount) :
        TBase(delegate, readsCount),
        errors(0)
    {
        TMatch match;
        fill(match, 0, 0, 0, 0, MaxValue<unsigned char>::VALUE, false);
        resize(matches, readsCount, match, Exact());
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _fixReverseComplement()                              [MatchManager]
// ----------------------------------------------------------------------------

// TODO(esiragusa): Remove this.
template <typename TDelegate, typename TSpec, typename TMatch, typename TReadId>
inline TReadId _fixReverseComplement(MatchManager<TDelegate, TSpec, TMatch> const & manager, TReadId readId)
{
    // Deal with reverse complemented reads.
    return (readId < manager.readsCount) ? readId : readId - manager.readsCount;
}

// ----------------------------------------------------------------------------
// Function onMatch()                                            [MatchManager]
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TSpec, typename TMatch,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchManager<TDelegate, TSpec, TMatch> & manager,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // Call matches delegate.
    onMatch(manager.delegate, contigId, beginPos, endPos, readId, errors, reverseComplemented);

    // Increment matches counters.
    manager.matchesCount++;
}

template <typename TDelegate, typename TMatch,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchManager<TDelegate, AllBest, TMatch> & manager,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // This match is useful.
    if (errors <= manager.minErrors[readId])
    {
        // Call matches delegate.
        onMatch(manager.delegate, contigId, beginPos, endPos, readId, errors, reverseComplemented);

        // Increment matches counters.
        manager.matchesCount++;

        // Disable the read after current stratum.
        manager.minErrors[readId] = errors;
    }

    // Otherwise this match is superfluous.
}

template <typename TDelegate, typename TMatch,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchManager<TDelegate, KBest, TMatch> & manager,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // This match is useful.
    if (errors < manager.minErrors[readId])
    {
        if (errors <= manager.errors)
        {
            // Call matches delegate.
            onMatch(manager.delegate, contigId, beginPos, endPos, readId, errors, reverseComplemented);

            // Increment matches counters.
            manager.matchesCount++;
        }
        else
        {
            // Store match.
            TMatch match;
            fill(match, contigId, beginPos, endPos, readId, errors, reverseComplemented);
            appendValue(manager.matches[errors], match, Generous());
        }
        // Disable the read.
        manager.minErrors[readId] = errors;
    }

    // Otherwise this match is superfluous.
}

template <typename TDelegate, typename TMatch,
typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchManager<TDelegate, AnyBest, TMatch> & manager,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // NOTE(esiragusa): This match should always be useful...
    if (errors < manager.matches[readId].errors)
    {
        // Store match.
        fill(manager.matches[readId], contigId, beginPos, endPos, readId, errors, reverseComplemented);

        if (errors == manager.errors)
        {
            // Call matches delegate.
            onMatch(manager.delegate, manager.matches[readId]);

            // Increment matches counters.
            manager.matchesCount++;
        }
    }
}

// ----------------------------------------------------------------------------
// Function isDisabled()                                         [MatchManager]
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TSpec, typename TMatch, typename TReadId>
inline bool isDisabled(MatchManager<TDelegate, TSpec, TMatch> const &, TReadId)
{
    return false;
}

template <typename TDelegate, typename TMatch, typename TReadId>
inline bool isDisabled(MatchManager<TDelegate, AllBest, TMatch> const & manager, TReadId readId)
{
    // Reads with at least one match in lower strata are disabled.
    return manager.minErrors[_fixReverseComplement(manager, readId)] < manager.errors;
}

template <typename TDelegate, typename TMatch, typename TReadId>
inline bool isDisabled(MatchManager<TDelegate, KBest, TMatch> const & manager, TReadId readId)
{
    // Reads with at least one best match are disabled.
    return manager.minErrors[_fixReverseComplement(manager, readId)] <= manager.errors;
}

template <typename TDelegate, typename TMatch, typename TReadId>
inline bool isDisabled(MatchManager<TDelegate, AnyBest, TMatch> const & manager, TReadId readId)
{
    // Reads with at least one best match are disabled.
    return manager.matches[_fixReverseComplement(manager, readId)].errors <= manager.errors;
}

// ----------------------------------------------------------------------------
// Function minErrors()                                          [MatchManager]
// ----------------------------------------------------------------------------

//template <typename TDelegate, typename TSpec, typename TMatch, typename TReadId>
//inline unsigned char minErrors(MatchManager<TDelegate, TSpec, TMatch> const &, TReadId)
//{
//    // TODO
//    return 0;
//}
//
//template <typename TDelegate, typename TMatch, typename TReadId>
//inline unsigned char minErrors(MatchManager<TDelegate, AnyBest, TMatch> const & manager, TReadId readId)
//{
//    return manager.matches[_fixReverseComplement(manager, readId)].errors + 1;
//}

// ----------------------------------------------------------------------------
// Function maxErrors()                                          [MatchManager]
// ----------------------------------------------------------------------------

//template <typename TDelegate, typename TSpec, typename TMatch, typename TReadId>
//inline unsigned char maxErrors(MatchManager<TDelegate, TSpec, TMatch> const &, TReadId)
//{
//    // TODO
//    return 0;
//}
//
//template <typename TDelegate, typename TMatch, typename TReadId>
//inline unsigned char maxErrors(MatchManager<TDelegate, AnyBest, TMatch> const & manager, TReadId readId)
//{
//    // TODO
//    return 0;
//}

// ----------------------------------------------------------------------------
// Function raiseErrorThreshold()                                [MatchManager]
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TSpec, typename TMatch>
inline void raiseErrorThreshold(MatchManager<TDelegate, TSpec, TMatch> &)
{}

template <typename TDelegate, typename TMatch>
inline void raiseErrorThreshold(MatchManager<TDelegate, AllBest, TMatch> & manager)
{
    // Increment error threshold.
    manager.errors++;
}

template <typename TDelegate, typename TMatch>
inline void raiseErrorThreshold(MatchManager<TDelegate, KBest, TMatch> & manager)
{
    typedef String<TMatch>                                  TMatches;
    typedef typename Iterator<TMatches, Standard>::Type     TMatchesIterator;

    TMatchesIterator matchesEnd = end(manager.matches[manager.errors + 1], Standard());

    for (TMatchesIterator matchesIt = begin(manager.matches[manager.errors + 1], Standard()); matchesIt != matchesEnd; ++matchesIt)
    {
        // If a read is disabled we already wrote its best match.
        if (isDisabled(manager, (*matchesIt).readId))
            continue;

        // Call matches delegate.
        onMatch(manager.delegate, *matchesIt);

        // Increment matches counter.
        manager.matchesCount++;

        // Disable the read.
        manager.minErrors[(*matchesIt).readId] = manager.errors + 1;
    }

    // Increment error threshold.
    manager.errors++;

    // Forget matches below error threshold.
    clear(manager.matches[manager.errors]);
    shrinkToFit(manager.matches[manager.errors]);
}

template <typename TDelegate, typename TMatch>
inline void raiseErrorThreshold(MatchManager<TDelegate, AnyBest, TMatch> & manager)
{
    typedef String<TMatch>                                  TMatches;
    typedef typename Iterator<TMatches, Standard>::Type     TMatchesIterator;

    // Increment error threshold.
    manager.errors++;

    TMatchesIterator matchesEnd = end(manager.matches, Standard());

    for (TMatchesIterator matchesIt = begin(manager.matches, Standard()); matchesIt != matchesEnd; ++matchesIt)
    {
        // Write current best matches.
        if ((*matchesIt).errors == manager.errors)
        {
            // Call matches delegate.
            onMatch(manager.delegate, *matchesIt);

            // Increment matches counter.
            manager.matchesCount++;

        }
    }
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_MATCHES_H_
