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

template <typename TMatch, typename TMatchesDelegate, typename TSpec = void>
struct MatchManager
{
    TMatchesDelegate & matchesDelegate;

    TReadSeqStoreSize       readsCount;
    unsigned long           matchesCount;

    MatchManager(TMatchesDelegate & matchesDelegate, TReadSeqStoreSize readsCount) :
        matchesDelegate(matchesDelegate),
        readsCount(readsCount),
        matchesCount(0)
    {}
};

template <typename TMatch, typename TMatchesDelegate>
struct MatchManager<TMatch, TMatchesDelegate, AllBest>:
    public MatchManager<TMatch, TMatchesDelegate>
{
    typedef MatchManager<TMatch, TMatchesDelegate>  TBase;
    typedef String<TMatch>                          TMatches;

    unsigned                    errors;
    String<unsigned char>       minErrors;

    MatchManager(TMatchesDelegate & matchesDelegate, TReadSeqStoreSize readsCount) :
        TBase(matchesDelegate, readsCount),
        errors(0)
    {
        resize(minErrors, readsCount, MaxValue<unsigned char>::VALUE, Exact());
    }
};

template <typename TMatch, typename TMatchesDelegate>
struct MatchManager<TMatch, TMatchesDelegate, AnyBest>:
    public MatchManager<TMatch, TMatchesDelegate>
{
    typedef MatchManager<TMatch, TMatchesDelegate>  TBase;
    typedef String<TMatch>                          TMatches;

    unsigned                    errors;
    String<unsigned char>       minErrors;
    String<TMatches>            matches;

    MatchManager(TMatchesDelegate & matchesDelegate, TReadSeqStoreSize readsCount) :
        TBase(matchesDelegate, readsCount),
        errors(0)
    {
        resize(minErrors, readsCount, MaxValue<unsigned char>::VALUE, Exact());
        // TODO Change hardcoded size.
        resize(matches, 16, Exact());
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function onMatch()                                            [MatchManager]
// ----------------------------------------------------------------------------

template <typename TMatch, typename TMatchesDelegate, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchManager<TMatch, TMatchesDelegate, TSpec> & manager,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // Call matches delegate.
    onMatch(manager.matchesDelegate, contigId, beginPos, endPos, readId, errors, reverseComplemented);

    // Increment matches counters.
    manager.matchesCount++;
}

template <typename TMatch, typename TMatchesDelegate,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchManager<TMatch, TMatchesDelegate, AllBest> & manager,
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
        onMatch(manager.matchesDelegate, contigId, beginPos, endPos, readId, errors, reverseComplemented);

        // Increment matches counters.
        manager.matchesCount++;

        // Disable the read after current stratum.
        manager.minErrors[readId] = errors;
    }

    // Otherwise this match is superfluous.
}

template <typename TMatch, typename TMatchesDelegate,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchManager<TMatch, TMatchesDelegate, AnyBest> & manager,
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
            onMatch(manager.matchesDelegate, contigId, beginPos, endPos, readId, errors, reverseComplemented);

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

// ----------------------------------------------------------------------------
// Function isDisabled()                                         [MatchManager]
// ----------------------------------------------------------------------------

template <typename TMatch, typename TMatchesDelegate, typename TSpec, typename TReadId>
inline bool isDisabled(MatchManager<TMatch, TMatchesDelegate, TSpec> &, TReadId)
{
    return false;
}

template <typename TMatch, typename TMatchesDelegate, typename TReadId>
inline bool isDisabled(MatchManager<TMatch, TMatchesDelegate, AllBest> & manager, TReadId readId)
{
    // Deal with reverse complemented reads.
    // NOTE This should be done outside of Manager, i.e. in Seeder.
    if (readId >= manager.readsCount)
        readId -= manager.readsCount;

    // Reads with at least one match in lower strata are disabled.
    return manager.minErrors[readId] < manager.errors;
}

template <typename TMatch, typename TMatchesDelegate, typename TReadId>
inline bool isDisabled(MatchManager<TMatch, TMatchesDelegate, AnyBest> & manager, TReadId readId)
{
    // Deal with reverse complemented reads.
    // NOTE This should be done outside of Manager, i.e. in Seeder.
    if (readId >= manager.readsCount)
        readId -= manager.readsCount;

    // Reads with at least one best match are disabled.
    return manager.minErrors[readId] <= manager.errors;
}


// ----------------------------------------------------------------------------
// Function raiseErrorThreshold()                                [MatchManager]
// ----------------------------------------------------------------------------

template <typename TMatch, typename TMatchesDelegate, typename TSpec>
inline void raiseErrorThreshold(MatchManager<TMatch, TMatchesDelegate, TSpec> &)
{}

template <typename TMatch, typename TMatchesDelegate>
inline void raiseErrorThreshold(MatchManager<TMatch, TMatchesDelegate, AllBest> & manager)
{
    // Increment error threshold.
    manager.errors++;
}

template <typename TMatch, typename TMatchesDelegate>
inline void raiseErrorThreshold(MatchManager<TMatch, TMatchesDelegate, AnyBest> & manager)
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
        onMatch(manager.matchesDelegate, *matchesIt);

        // Increment matches counters.
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

#endif  // #ifndef SEQAN_EXTRAS_MASAI_MATCHES_H_
