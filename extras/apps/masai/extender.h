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
// This file contains the Extender class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_EXTENDER_H_
#define SEQAN_EXTRAS_MASAI_EXTENDER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/find.h>

#include "store.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Extender
// ----------------------------------------------------------------------------

template <typename TMatchesDelegate, typename TDistance = HammingDistance, typename TSpec = void>
struct Extender
{
    TFragmentStore &        store;
    TMatchesDelegate &      matchesDelegate;
    String<TContigSeqSize>  contigSizes;
    TReadSeqStoreSize       readsCount;

    TReadSeqSize            minErrorsPerRead;
    TReadSeqSize            maxErrorsPerRead;
    TReadSeqSize            seedLength;

    bool                    disabled;

    Extender(TFragmentStore & store,
             TMatchesDelegate & matchesDelegate,
             TReadSeqStoreSize readsCount,
             bool disabled = false) :
        store(store),
        matchesDelegate(matchesDelegate),
        readsCount(readsCount),
        minErrorsPerRead(0),
        maxErrorsPerRead(0),
        seedLength(0),
        disabled(disabled)
    {
        _init(*this);
    }

};

template <typename TMatchesDelegate, typename TSpec>
struct Extender<TMatchesDelegate, EditDistance, TSpec>:
    public Extender<TMatchesDelegate>
{
    typedef Extender<TMatchesDelegate>                      TBase;
    typedef Myers<AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_>, True, void> TAlgorithmSpec;

    typedef Segment<TReadSeq, InfixSegment>                 TReadInfix;
    typedef ModifiedString<TReadInfix, ModReverse>          TReadInfixRev;

    typedef PatternState_<TReadInfix, TAlgorithmSpec>       TPatternState;
    typedef PatternState_<TReadInfixRev, TAlgorithmSpec>    TPatternStateRev;

    TPatternState patternState;
    TPatternStateRev patternStateRev;

    Extender(TFragmentStore & store,
             TMatchesDelegate & matchesDelegate,
             TReadSeqStoreSize readsCount,
             bool disabled = false) :
        TBase(store, matchesDelegate, readsCount, disabled)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// TODO(esiragusa): Move this into Genome class.
template <typename TMatchesDelegate, typename TDistance, typename TSpec>
inline void _init(Extender<TMatchesDelegate, TDistance, TSpec> & extender)
{
    reserve(extender.contigSizes, length(extender.store.contigStore), Exact());
    for (TContigStoreSize contigId = 0; contigId < length(extender.store.contigStore); ++contigId)
        appendValue(extender.contigSizes, length(extender.store.contigStore[contigId].seq));
}

// TODO(esiragusa): Remove this.
template <typename TMatchesDelegate, typename TDistance, typename TSpec>
inline bool _fixReverseComplemented(Extender<TMatchesDelegate, TDistance, TSpec> & extender, TReadSeqStoreSize & readId)
{
    // Deal with reverse complemented reads.
    if (readId >= extender.readsCount)
    {
        readId -= extender.readsCount;
        return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function onSeedHit()                                              [Extender]
// ----------------------------------------------------------------------------

template <typename TMatchesDelegate, typename TSpec>
inline bool onSeedHit(Extender<TMatchesDelegate, HammingDistance, TSpec> & extender,
                      TContigStoreSize contigId,
                      TContigSeqSize contigBegin,
                      TReadSeqStoreSize readId,
                      TContigSeqSize seedBegin,
                      TReadSeqSize seedErrors)
{
    typedef Segment<TReadSeq, InfixSegment>                 TReadInfix;

    if (extender.disabled)
        return false;

    TReadSeqSize errors = seedErrors;

    TContigSeq & contig = extender.store.contigStore[contigId].seq;
    TReadSeq & read = extender.store.readSeqStore[readId];
    TReadSeqSize readLength = length(read);

    // Extend left.
    TContigSeqSize matchBegin = contigBegin;

    if (seedBegin > 0)
    {
        TContigSeqSize contigLeftBegin = 0;
        if (contigBegin > seedBegin)
            contigLeftBegin = contigBegin - seedBegin;

        TContigInfix contigLeft(contig, contigLeftBegin, contigBegin);
        TReadInfix readLeft(read, 0, seedBegin);

        if (!_extend(extender, contigLeft, readLeft, errors))
            return false;

        matchBegin = contigLeftBegin;
    }

    // This removes some duplicates.
    if (errors - seedErrors < extender.minErrorsPerRead)
        return false;

    // Extend right.
    TContigSeqSize matchEnd = contigBegin + extender.seedLength;

    if (seedBegin + extender.seedLength < readLength)
    {
        TContigSeqSize contigRightEnd = extender.contigSizes[contigId];
        if (contigRightEnd > contigBegin + readLength - seedBegin)
            contigRightEnd = contigBegin + readLength - seedBegin;

        TContigInfix contigRight(contig, contigBegin + extender.seedLength, contigRightEnd);
        TReadInfix readRight(read, seedBegin + extender.seedLength, readLength);

        if (!_extend(extender, contigRight, readRight, errors))
            return false;

        matchEnd = contigRightEnd;
    }

    // This removes some duplicates.
    if (errors < extender.minErrorsPerRead)
        return false;

    bool reverseComplemented = _fixReverseComplemented(extender, readId);
    onMatch(extender.matchesDelegate, contigId, matchBegin, matchEnd, readId, errors, reverseComplemented);

    return true;
}

template <typename TMatchesDelegate, typename TSpec, typename TContigInfix, typename TReadInfix>
inline bool _extend(Extender<TMatchesDelegate, HammingDistance, TSpec> & extender,
                    TContigInfix & contigInfix,
                    TReadInfix & readInfix,
                    TReadSeqSize & errors)
{
    typedef typename Iterator<TContigInfix, Standard>::Type          TContigInfixIterator;
    typedef typename Iterator<TReadInfix, Standard>::Type            TReadInfixIterator;

    if (length(contigInfix) != length(readInfix))
        return false;

    TContigInfixIterator contigIt = begin(contigInfix, Standard());
    TReadInfixIterator readBegin = begin(readInfix, Standard());
    TReadInfixIterator readEnd = end(readInfix, Standard());

    for (TReadInfixIterator readIt = readBegin; readIt != readEnd; ++readIt, ++contigIt)
        if (*readIt != *contigIt)
            if (++errors > extender.maxErrorsPerRead)
                return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function onSeedHit()                                              [Extender]
// ----------------------------------------------------------------------------

template <typename TMatchesDelegate, typename TSpec>
inline bool onSeedHit(Extender<TMatchesDelegate, EditDistance, TSpec> & extender,
                      TContigStoreSize contigId,
                      TContigSeqSize contigBegin,
                      TReadSeqStoreSize readId,
                      TContigSeqSize seedBegin,
                      TReadSeqSize seedErrors)
{
    typedef Segment<TReadSeq, InfixSegment> TReadInfix;

    if (extender.disabled)
        return false;

    TReadSeqSize errors = seedErrors;

    TContigSeq & contig = extender.store.contigStore[contigId].seq;
    TReadSeq & read = extender.store.readSeqStore[readId];
    TReadSeqSize readLength = length(read);

    // Extend left.
    TContigSeqSize matchBegin = contigBegin;

    if (seedBegin > 0)
    {
        TContigSeqSize contigLeftBegin = 0;
        if (contigBegin > seedBegin + extender.maxErrorsPerRead - errors)
            contigLeftBegin = contigBegin - (seedBegin + extender.maxErrorsPerRead - errors);

        TContigInfix contigLeft(contig, contigLeftBegin, contigBegin);
        TReadInfix readLeft(read, 0, seedBegin);

        if (!_extendLeft(extender, extender.patternStateRev, contigLeft, readLeft, errors, matchBegin))
            return false;
    }

    // This removes some duplicates.
    if (errors - seedErrors < extender.minErrorsPerRead)
        return false;

    // Extend right.
    TContigSeqSize matchEnd = contigBegin + extender.seedLength;

    if (seedBegin + extender.seedLength < readLength)
    {
        TContigSeqSize contigRightEnd = extender.contigSizes[contigId];
        if (contigRightEnd > contigBegin + readLength - seedBegin + extender.maxErrorsPerRead - errors)
            contigRightEnd = contigBegin + readLength - seedBegin + extender.maxErrorsPerRead - errors;

        if (contigBegin + extender.seedLength >= contigRightEnd)
            return false;

        TContigInfix contigRight(contig, contigBegin + extender.seedLength, contigRightEnd);
        TReadInfix readRight(read, seedBegin + extender.seedLength, readLength);

        if (!_extendRight(extender, extender.patternState, contigRight, readRight, errors, matchEnd))
            return false;
    }

    // This removes some duplicates.
    if (errors < extender.minErrorsPerRead)
        return false;

    bool reverseComplemented = _fixReverseComplemented(extender, readId);
    onMatch(extender.matchesDelegate, contigId, matchBegin, matchEnd, readId, errors, reverseComplemented);

    return true;
}

template <typename TMatchesDelegate, typename TSpec, typename TPatternState, typename TContigInfix, typename TReadInfix>
inline bool _extendLeft(Extender<TMatchesDelegate, EditDistance, TSpec> & extender,
                        TPatternState & patternState,
                        TContigInfix & contigInfix,
                        TReadInfix & readInfix,
                        TReadSeqSize & errors,
                        TContigSeqSize & matchBegin)
{
    typedef ModifiedString<TReadInfix, ModReverse>          TReadInfixRev;
    typedef ModifiedString<TContigInfix, ModReverse>        TContigInfixRev;
    typedef Finder<TContigInfixRev>                         TFinder;

    TContigInfixRev contigInfixRev(contigInfix);
    TReadInfixRev readInfixRev(readInfix);

    // Lcp trick.
    TContigSeqSize lcp = lcpLength(contigInfixRev, readInfixRev);
    if (lcp == length(readInfix))
    {
        matchBegin -= lcp;
        return true;
    }
    setEndPosition(contigInfix, endPosition(contigInfix) - lcp);
    setEndPosition(readInfix, endPosition(readInfix) - lcp);

    TReadSeqSize remainingErrors = extender.maxErrorsPerRead - errors;
    TReadSeqSize minErrors = remainingErrors + 1;
    TContigSeqSize endPos = 0;

    // Stop seed extension.
    if (!remainingErrors)
        return false;

    // Align.
    TFinder finder(contigInfixRev);
    patternState.leftClip = remainingErrors;

    // TODO(esiragusa): Use a generic type for errors.
    while (find(finder, readInfixRev, patternState, -static_cast<int>(remainingErrors)))
    {
        TReadSeqSize currentErrors = -getScore(patternState);

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = position(finder) + 1;
        }
    }

    errors += minErrors;
    matchBegin -= endPos + lcp;

    return errors <= extender.maxErrorsPerRead;
}

template <typename TMatchesDelegate, typename TSpec, typename TPatternState, typename TContigInfix, typename TReadInfix>
inline bool _extendRight(Extender<TMatchesDelegate, EditDistance, TSpec> & extender,
                         TPatternState & patternState,
                         TContigInfix & contigInfix,
                         TReadInfix & readInfix,
                         TReadSeqSize & errors,
                         TContigSeqSize & matchEnd)
{
    typedef Finder<TContigInfix>    TFinder;

    // Lcp trick.
    TContigSeqSize lcp = lcpLength(contigInfix, readInfix);
    if (lcp == length(readInfix))
    {
        matchEnd += lcp;
        return true;
    }
    else if (lcp == length(contigInfix))
    {
        errors += length(readInfix) - length(contigInfix);
        matchEnd += length(readInfix);
        return errors <= extender.maxErrorsPerRead;
    }
    setBeginPosition(contigInfix, beginPosition(contigInfix) + lcp);
    setBeginPosition(readInfix, beginPosition(readInfix) + lcp);

    // NOTE Uncomment this to disable lcp trick.
//    TContigSeqSize lcp = 0;

    TReadSeqSize remainingErrors = extender.maxErrorsPerRead - errors;
    TReadSeqSize minErrors = remainingErrors + 1;
    TContigSeqSize endPos = 0;

    // NOTE Comment this to disable lcp trick.
    // Stop seed extension.
    if (!remainingErrors)
        return false;

    // Remove last base.
    TContigInfix contigPrefix(contigInfix);
    TReadInfix readPrefix(readInfix);
    setEndPosition(contigPrefix, endPosition(contigPrefix) - 1);
    setEndPosition(readPrefix, endPosition(readPrefix) - 1);

    // Align.
    TFinder finder(contigPrefix);
    patternState.leftClip = remainingErrors;

    // TODO(esiragusa): Use a generic type for errors.
    while (find(finder, readPrefix, patternState, -static_cast<int>(remainingErrors)))
    {
        TContigSeqSize currentEnd = position(finder) + 1;
        TReadSeqSize currentErrors = -getScore(patternState);

        // Compare last base.
        if (contigInfix[currentEnd] != back(readInfix))
            if (++currentErrors > remainingErrors)
                continue;

        if (currentErrors <= minErrors)
        {
            minErrors = currentErrors;
            endPos = currentEnd;
        }
    }

    errors += minErrors;
    matchEnd += endPos + lcp + 1;

    return errors <= extender.maxErrorsPerRead;
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_EXTENDER_H_
