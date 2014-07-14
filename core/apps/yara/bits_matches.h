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
// This file contains classes for storing and manipulating matches.
// ==========================================================================

#ifndef APP_YARA_BITS_MATCHES_H_
#define APP_YARA_BITS_MATCHES_H_

// ============================================================================
// Extras
// ============================================================================

namespace seqan {

// ----------------------------------------------------------------------------
// Class Adder
// ----------------------------------------------------------------------------

template <typename TUnaryFunction, unsigned DELTA>
struct Adder
{
    TUnaryFunction const & f;

    Adder(TUnaryFunction const & f) : f(f) {}

    template <typename TValue>
    unsigned operator() (TValue const & val) const
    {
        return f(val) + DELTA;
    }
};

// ----------------------------------------------------------------------------
// Class KeyIndicator
// ----------------------------------------------------------------------------

template <typename TTarget, typename TKey, typename TSpec = void>
struct KeyIndicator
{
    TTarget &       target;
    TKey const &    key;

    KeyIndicator(TTarget & target, TKey const & key) :
        target(target),
        key(key)
    {}

    template <typename TValue>
    void operator() (TValue const & val) const
    {
        SEQAN_ASSERT_LT(key(val), length(target));
        target[key(val)] = true;
    }
};

// ----------------------------------------------------------------------------
// Class KeyCounter
// ----------------------------------------------------------------------------

template <typename TTarget, typename TKey, typename TThreading = Serial, typename TSpec = void>
struct KeyCounter
{
    TTarget &       target;
    TKey const &    key;

    KeyCounter(TTarget & target, TKey const & key) :
        target(target),
        key(key)
    {}

    template <typename TValue>
    void operator() (TValue const & val) const
    {
        SEQAN_ASSERT_LT(key(val), length(target));
        atomicInc(target[key(val)], TThreading());
    }
};

// ----------------------------------------------------------------------------
// Class KeySorter
// ----------------------------------------------------------------------------

template <typename TSource, typename TSpec = void>
struct KeySorter
{
    TSource const & source;

    KeySorter(TSource const & source) :
        source(source)
    {}

    template <typename TKey>
    inline bool operator()(TKey a, TKey b) const
    {
        return source[a] < source[b];
    }
};

// --------------------------------------------------------------------------
// Function bucket()
// --------------------------------------------------------------------------
// Bucket elements in the concat of a ConcatDirect StringSet.
// Remarks: the concat string must be already sorted by key.

template <typename TString, typename TSpec, typename TKeyGetter, typename TMaxKey, typename TThreading>
inline void bucket(StringSet<TString, Owner<ConcatDirect<TSpec > > > & me, TKeyGetter const & key, TMaxKey maxKey, TThreading const & tag)
{
    typedef StringSet<TString, Owner<ConcatDirect<TSpec > > >    TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type           TLimits;
    typedef Adder<TKeyGetter, 1u>                                TNextKey;
    typedef KeyCounter<TLimits, TNextKey, TThreading>            TCounter;

    if (empty(concat(me))) return;

    // Shift the counts by one.
    TNextKey nextKey(key);

    // Resize the limits string to count all keys.
    resize(me.limits, maxKey, 0, Exact());

    // Count the number of keys present in the concat string.
    forEach(concat(me), TCounter(me.limits, nextKey), tag);

    // Build the limits string from the key counts.
    partialSum(me.limits, tag);
}

// --------------------------------------------------------------------------
// Function bucket()
// --------------------------------------------------------------------------
// Bucket elements in the host of a Segment StringSet.
// Remarks: the host string must be already sorted by key.

template <typename THost, typename TSpec, typename TKeyGetter, typename TMaxKey, typename TThreading>
inline void bucket(StringSet<THost, Segment<TSpec> > & me, TKeyGetter const & key, TMaxKey maxKey, TThreading const & tag)
{
    typedef StringSet<THost, Segment<TSpec> >                    TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type           TLimits;
    typedef Adder<TKeyGetter, 1u>                                TNextKey;
    typedef KeyCounter<TLimits, TNextKey, TThreading>            TCounter;

    if (empty(host(me))) return;

    // Shift the key counts by one.
    TNextKey nextKey(key);

    // Resize the limits string to accomodate counts for all keys.
    resize(me.limits, maxKey + 1, 0, Exact());

    // Count the number of keys present in the host string.
    forEach(host(me), TCounter(me.limits, nextKey), tag);

    // Limits are the cumulated key counts.
    partialSum(me.limits, tag);

    // Positions are the shifted limits.
    assign(me.positions, prefix(me.limits, length(me.limits) - 1));
}

}

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

template <typename TObject, typename TMember>
struct Getter;

// ============================================================================
// Tags
// ============================================================================

struct SortErrors_;
typedef Tag<SortErrors_> const SortErrors;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Match
// ----------------------------------------------------------------------------

#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif

template <typename TSpec = void>
struct Match
{
    unsigned        readId       : YaraBits<TSpec>::READ_ID;
    unsigned char   contigId; // : YaraBits<TSpec>::CONTIG_ID;
    bool            isRev        : 1;
    unsigned        contigBegin  : YaraBits<TSpec>::CONTIG_SIZE;
    unsigned short  contigEnd    : YaraBits<TSpec>::READ_SIZE;
    unsigned        errors       : YaraBits<TSpec>::ERRORS;
}
#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
;

#ifdef PLATFORM_WINDOWS
      #pragma pack(pop)
#endif

// ----------------------------------------------------------------------------
// Class Match Getter
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Getter<Match<TSpec>, SortReadId>
{
    inline unsigned operator()(Match<TSpec> const & me) const
    {
        return getReadId(me);
    }
};

// ----------------------------------------------------------------------------
// Class MatchSorter
// ----------------------------------------------------------------------------

template <typename TMatch, typename TSpec>
struct MatchSorter;

template <typename TMatch>
struct MatchSorter<TMatch, SortReadId>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getReadId(a) < getReadId(b);
    }
};

template <typename TMatch>
struct MatchSorter<TMatch, SortBeginPos>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getSortKey(a, SortBeginPos()) < getSortKey(b, SortBeginPos());
    }
};

template <typename TMatch>
struct MatchSorter<TMatch, SortEndPos>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getSortKey(a, SortEndPos()) < getSortKey(b, SortEndPos());
    }
};

template <typename TMatch>
struct MatchSorter<TMatch, SortErrors>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getErrors(a) < getErrors(b);
    }
};

// ----------------------------------------------------------------------------
// Class MatchesCompactor
// ----------------------------------------------------------------------------

template <typename TCounts, typename TPosition>
struct MatchesCompactor
{
    TCounts &    unique;

    MatchesCompactor(TCounts & unique) :
        unique(unique)
    {}

    template <typename TIterator>
    void operator() (TIterator & it)
    {
        typedef typename Value<TIterator>::Type     TMatches;
        typedef typename Value<TMatches>::Type      TMatch;

        TMatches const & matches = value(it);

        sort(matches, MatchSorter<TMatch, TPosition>());
        unique[position(it) + 1] = compactUniqueMatches(matches, TPosition());
    }
};

// ----------------------------------------------------------------------------
// Class MatchesPicker
// ----------------------------------------------------------------------------

template <typename TMatches>
struct MatchesPicker
{
    typedef typename Value<TMatches>::Type TMatch;

    Rng<MersenneTwister> rng;
    TMatch invalid;

    MatchesPicker() :
        rng(0xABAD1DEA)
    {
        setInvalid(invalid);
    }

    TMatch operator() (TMatches const & matches)
    {
        return empty(matches) ? invalid : matches[pickRandomNumber(rng) % length(matches)];
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Match Setters
// ----------------------------------------------------------------------------

template <typename TSpec, typename TReadSeqs, typename TReadSeqId>
inline void setReadId(Match<TSpec> & me, TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    me.readId = getReadId(readSeqs, readSeqId);
    me.isRev = isRevReadSeq(readSeqs, readSeqId);
}

template <typename TSpec, typename TContigPos>
inline void setContigPosition(Match<TSpec> & me, TContigPos contigBegin, TContigPos contigEnd)
{
    SEQAN_ASSERT_EQ(getValueI1(contigBegin), getValueI1(contigEnd));
    SEQAN_ASSERT_LT(getValueI2(contigBegin), getValueI2(contigEnd));

    me.contigId = getValueI1(contigBegin);
    me.contigBegin = getValueI2(contigBegin);
    me.contigEnd = getValueI2(contigEnd) - getValueI2(contigBegin);
}

template <typename TSpec>
inline void setInvalid(Match<TSpec> & me)
{
    me.readId = 0;
    me.contigId = 0;
    me.isRev = 0;
    me.contigBegin = 0;
    me.contigEnd = 0;
    me.errors = YaraLimits<TSpec>::ERRORS;
}

// ----------------------------------------------------------------------------
// Match Getters
// ----------------------------------------------------------------------------
// NOTE(esiragusa): getReadId<TMatch>() could be getValue(Match<TSpec>, TMember())

template <typename TReadSeqs, typename TSpec>
inline typename Size<TReadSeqs>::Type
getReadSeqId(Match<TSpec> const & me, TReadSeqs const & readSeqs)
{
    return onForwardStrand(me) ? getFirstMateFwdSeqId(readSeqs, me.readId) : getFirstMateRevSeqId(readSeqs, me.readId);
}

template <typename TSpec>
inline unsigned getReadId(Match<TSpec> const & me)
{
    return me.readId;
}

template <typename TSpec>
inline unsigned char getContigId(Match<TSpec> const & me)
{
    return me.contigId;
}

template <typename TSpec>
inline unsigned getContigBegin(Match<TSpec> const & me)
{
    return me.contigBegin;
}

template <typename TSpec>
inline unsigned getContigEnd(Match<TSpec> const & me)
{
    return me.contigBegin + me.contigEnd;
}

template <typename TSpec>
inline bool onReverseStrand(Match<TSpec> const & me)
{
    return me.isRev;
}

template <typename TSpec>
inline bool onForwardStrand(Match<TSpec> const & me)
{
    return !onReverseStrand(me);
}

template <typename TSpec>
inline unsigned char getErrors(Match<TSpec> const & me)
{
    return me.errors;
}

template <typename TSpec>
inline bool isInvalid(Match<TSpec> const & me)
{
    return getContigBegin(me) == getContigEnd(me);
}

template <typename TSpec>
inline bool isValid(Match<TSpec> const & me)
{
    return !isInvalid(me);
}

// ----------------------------------------------------------------------------
// Function getErrors()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned getErrors(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return (unsigned)a.errors + b.errors;
}

// ----------------------------------------------------------------------------
// Function getTemplateLength()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned getTemplateLength(Match<TSpec> const & a, Match<TSpec> const & b)
{
    if (getContigBegin(a) < getContigBegin(b))
        return getContigEnd(b) - getContigBegin(a);
    else
        return getContigEnd(a) - getContigBegin(b);
}

// ----------------------------------------------------------------------------
// Function cigarLength()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned getCigarLength(Match<TSpec> const & me)
{
    return isInvalid(me) ? 0 : 2 * getErrors(me) + 1;
}

// ----------------------------------------------------------------------------
// Function getSortKey(SortBeginPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned long getSortKey(Match<TSpec> const & me, SortBeginPos)
{
    return ((unsigned long)getContigId(me)     << (1 + YaraBits<TSpec>::CONTIG_SIZE + YaraBits<TSpec>::ERRORS)) |
           ((unsigned long)onReverseStrand(me) << (YaraBits<TSpec>::CONTIG_SIZE + YaraBits<TSpec>::ERRORS))     |
           ((unsigned long)getContigBegin(me)  <<  YaraBits<TSpec>::ERRORS)                                     |
           ((unsigned long)getErrors(me));
}

// ----------------------------------------------------------------------------
// Function getSortKey(SortEndPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned long getSortKey(Match<TSpec> const & me, SortEndPos)
{
    return ((unsigned long)getContigId(me)     << (1 + YaraBits<TSpec>::CONTIG_SIZE + YaraBits<TSpec>::ERRORS)) |
           ((unsigned long)onReverseStrand(me) << (YaraBits<TSpec>::CONTIG_SIZE + YaraBits<TSpec>::ERRORS))     |
           ((unsigned long)getContigEnd(me)    <<  YaraBits<TSpec>::ERRORS)                                     |
           ((unsigned long)getErrors(me));
}

// ----------------------------------------------------------------------------
// Function strandEqual()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool strandEqual(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return !(onForwardStrand(a) ^ onForwardStrand(b));
}

// ----------------------------------------------------------------------------
// Function contigEqual()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool contigEqual(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return getContigId(a) == getContigId(b) && strandEqual(a, b);
}

// ----------------------------------------------------------------------------
// Function isDuplicate(SortBeginPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, SortBeginPos)
{
    return contigEqual(a, b) && getContigBegin(a) == getContigBegin(b);
}

// ----------------------------------------------------------------------------
// Function isDuplicate(SortEndPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, SortEndPos)
{
    return contigEqual(a, b) && getContigEnd(a) == getContigEnd(b);
}

// ----------------------------------------------------------------------------
// Function compactUniqueMatches()
// ----------------------------------------------------------------------------
// Compact unique matches at the beginning of their container.

template <typename TMatches, typename TPosition>
inline typename Size<TMatches>::Type
compactUniqueMatches(TMatches & matches, Tag<TPosition> const & posTag)
{
    typedef typename Iterator<TMatches, Standard>::Type         TMatchesIterator;

    TMatchesIterator matchesBegin = begin(matches, Standard());
    TMatchesIterator matchesEnd = end(matches, Standard());
    TMatchesIterator newIt = matchesBegin;
    TMatchesIterator oldIt = matchesBegin;

    while (oldIt != matchesEnd)
    {
        *newIt = *oldIt;

        ++oldIt;

        while (oldIt != matchesEnd && isDuplicate(*newIt, *oldIt, posTag)) ++oldIt;

        ++newIt;
    }

    return newIt - matchesBegin;
}

// ----------------------------------------------------------------------------
// Function removeDuplicates()
// ----------------------------------------------------------------------------
// Remove duplicate matches from a set of matches.

template <typename TMatchesSet, typename TThreading>
inline void removeDuplicates(TMatchesSet & matchesSet, TThreading const & threading)
{
    typedef typename StringSetLimits<TMatchesSet>::Type         TLimits;

    TLimits newLimits;
    resize(newLimits, length(stringSetLimits(matchesSet)), Exact());
    front(newLimits) = 0;

    // Sort matches by end position and move unique matches at the beginning.
    iterate(matchesSet, MatchesCompactor<TLimits, SortEndPos>(newLimits), Rooted(), threading);

    // Exclude duplicate matches at the end.
    assign(stringSetLimits(matchesSet), newLimits);
    _refreshStringSetLimits(matchesSet, threading);

    // Sort matches by begin position and move unique matches at the beginning.
    iterate(matchesSet, MatchesCompactor<TLimits, SortBeginPos>(newLimits), Rooted(), threading);

    // Exclude duplicate matches at the end.
    assign(stringSetLimits(matchesSet), newLimits);
    _refreshStringSetLimits(matchesSet, threading);
}

// ----------------------------------------------------------------------------
// Function countBestMatches()
// ----------------------------------------------------------------------------
// Count the number of cooptimal matches - ordering by errors is required.

template <typename TMatches>
inline typename Size<TMatches>::Type
countBestMatches(TMatches const & matches)
{
    typedef typename Iterator<TMatches const, Standard>::Type   TIter;
    typedef typename Size<TMatches>::Type                       TCount;

    TIter itBegin = begin(matches, Standard());
    TIter itEnd = end(matches, Standard());

    TCount count = 0;

    for (TIter it = itBegin; it != itEnd && getErrors(*it) <= getErrors(*itBegin); it++, count++) ;

    return count;
}

// ----------------------------------------------------------------------------
// Function removeSuboptimal()
// ----------------------------------------------------------------------------
// Remove suboptimal matches from a set of matches.

template <typename TMatchesSet, typename TThreading>
inline void removeSuboptimal(TMatchesSet & matchesSet, TThreading const & threading)
{
    typedef typename StringSetLimits<TMatchesSet>::Type         TLimits;
    typedef typename Suffix<TLimits>::Type                      TLimitsSuffix;
    typedef typename Value<TMatchesSet>::Type                   TMatches;

    TLimits newLimits;
    resize(newLimits, length(stringSetLimits(matchesSet)), Exact());
    SEQAN_ASSERT_GT(length(stringSetLimits(matchesSet)), 0);
    front(newLimits) = 0;

    // Count co-optimal matches.
    TLimitsSuffix counts = suffix(newLimits, 1);
    transform(counts, matchesSet, countBestMatches<TMatches>, threading);

    // Exclude suboptimal matches at the end.
    assign(stringSetLimits(matchesSet), newLimits);
    _refreshStringSetLimits(matchesSet, threading);
}

// ----------------------------------------------------------------------------
// Function findMatch()
// ----------------------------------------------------------------------------

template <typename TMatches, typename TMatch>
inline typename Iterator<TMatches const, Standard>::Type
findMatch(TMatches const & matches, TMatch const & match)
{
    typedef typename Iterator<TMatches const, Standard>::Type   TIter;

    TIter it = begin(matches, Standard());
    TIter itEnd = end(matches, Standard());

    for (; it != itEnd && !isDuplicate(*it, match, SortBeginPos()); it++) ;

    return it;
}

// ----------------------------------------------------------------------------
// Function sortMatches()
// ----------------------------------------------------------------------------

//template <typename TMatches, typename TKey>
//inline void sortMatches(TMatches & matches)
//{
//    typedef typename Value<TMatches>::Type  TMatch;
//
//    stableSort(matches, MatchSorter<TMatch, TKey>());
//}

template <typename TIterator, typename TKey>
inline void sortMatches(TIterator & it)
{
    typedef typename Value<TIterator>::Type TMatches;
    typedef typename Value<TMatches>::Type  TMatch;

    TMatches matches = value(it);
    stableSort(matches, MatchSorter<TMatch, TKey>());
}

// ----------------------------------------------------------------------------
// Function findSameContig()
// ----------------------------------------------------------------------------
// Find the first pair of matches on the same contig.

template <typename TMatchesIterator, typename TMatches>
inline bool findSameContig(TMatchesIterator & leftIt, TMatchesIterator & rightIt,
                           TMatches const & left, TMatches const & right)
{
    while (!atEnd(leftIt, left) && !atEnd(rightIt, right))
    {
        if (getContigId(*leftIt) < getContigId(*rightIt))
            findNextContig(leftIt, left, getContigId(*leftIt));
        else if (getContigId(*leftIt) > getContigId(*rightIt))
            findNextContig(rightIt, right, getContigId(*rightIt));
        else
            return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function findNextContig()
// ----------------------------------------------------------------------------
// Find the first match after given contigId.

template <typename TMatchesIterator, typename TMatches, typename TContigId>
inline void findNextContig(TMatchesIterator & it, TMatches const & matches, TContigId contigId)
{
    while (!atEnd(it, matches) && getContigId(*it) <= contigId) ++it;
}

// ----------------------------------------------------------------------------
// Function findReverseStrand()
// ----------------------------------------------------------------------------
// Find the first match on the reverse strand of the given contigId.

template <typename TMatchesIterator, typename TMatches, typename TContigId>
inline void findReverseStrand(TMatchesIterator & it, TMatches const & matches, TContigId contigId)
{
    while (!atEnd(it, matches) && (getContigId(*it) <= contigId) && onForwardStrand(*it)) ++it;
}

// ----------------------------------------------------------------------------
// Function bucketMatches()
// ----------------------------------------------------------------------------

template <typename TMatches, typename TDelegate>
inline void bucketMatches(TMatches const & left, TMatches const & right, TDelegate & delegate)
{
    typedef typename Iterator<TMatches const, Standard>::Type   TIterator;
    typedef typename Infix<TMatches const>::Type                TInfix;

    TIterator leftIt = begin(left, Standard());
    TIterator rightIt = begin(right, Standard());

    // Find matches on the same contig.
    while (findSameContig(leftIt, rightIt, left, right))
    {
        unsigned contigId = getContigId(*leftIt);

        TIterator leftBegin;
        TIterator rightBegin;

        // Find matches on forward strand.
        leftBegin = leftIt;
        rightBegin = rightIt;
        findReverseStrand(leftIt, left, contigId);
        findReverseStrand(rightIt, right, contigId);
        TInfix leftFwd = infix(left, position(leftBegin, left), position(leftIt, left));
        TInfix rightFwd = infix(right, position(rightBegin, right), position(rightIt, right));

        // Find matches on reverse strand.
        leftBegin = leftIt;
        rightBegin = rightIt;
        findNextContig(leftIt, left, contigId);
        findNextContig(rightIt, right, contigId);
        TInfix leftRev = infix(left, position(leftBegin, left), position(leftIt, left));
        TInfix rightRev = infix(right, position(rightBegin, right), position(rightIt, right));

        delegate(leftFwd, rightRev, FwdRev());
        delegate(leftFwd, rightFwd, FwdFwd());
        delegate(leftRev, rightFwd, RevFwd());
        delegate(leftRev, rightRev, RevRev());
    }
}

#endif  // #ifndef APP_YARA_BITS_MATCHES_H_
