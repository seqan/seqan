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

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ----------------------------------------------------------------------------
// Functor Getter
// ----------------------------------------------------------------------------

template <typename TObject, typename TTag>
struct Getter
{
    inline unsigned operator()(TObject const & me) const
    {
        return getValue(me, TTag());
    }
};

// ============================================================================
// Tags
// ============================================================================

struct MatchErrors_;
typedef Tag<MatchErrors_> const MatchErrors;

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
    typename Member<Match, SortReadId>::Type    readId       : MemberBits<Match, SortReadId>::VALUE;
    typename Member<Match, SortContigId>::Type  contigId; // : MemberBits<Match, SortContigId>::VALUE;
    bool                                        isRev        : 1;
    typename Member<Match, ContigSize>::Type    contigBegin  : MemberBits<Match, ContigSize>::VALUE;
    typename Member<Match, ReadSize>::Type      contigEnd    : MemberBits<Match, ReadSize>::VALUE;
    typename Member<Match, MatchErrors>::Type    errors       : MemberBits<Match, MatchErrors>::VALUE;
}
#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
;

#ifdef PLATFORM_WINDOWS
      #pragma pack(pop)
#endif

// ----------------------------------------------------------------------------
// Class MatchSorter
// ----------------------------------------------------------------------------

template <typename TMatch, typename TTag>
struct MatchSorter
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getValue(a, TTag()) < getValue(b, TTag());
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
    me.errors = MemberLimits<Match<TSpec>, MatchErrors>::VALUE;
}

// ----------------------------------------------------------------------------
// Match Getters
// ----------------------------------------------------------------------------

template <typename TSpec>
inline typename Member<Match<TSpec>, SortReadId>::Type
getValue(Match<TSpec> & me, SortReadId)
{
    return me.readId;
}
template <typename TSpec>
inline typename Member<Match<TSpec>, SortReadId>::Type
getValue(Match<TSpec> const & me, SortReadId)
{
    return me.readId;
}

template <typename TSpec>
inline typename Member<Match<TSpec>, SortContigId>::Type
getValue(Match<TSpec> const & me, SortContigId)
{
    return me.contigId;
}

template <typename TSpec>
inline typename Member<Match<TSpec>, SortBeginPos>::Type
getValue(Match<TSpec> const & me, SortBeginPos)
{
    return me.contigBegin;
}

template <typename TSpec>
inline typename Member<Match<TSpec>, SortEndPos>::Type
getValue(Match<TSpec> const & me, SortEndPos)
{
    return me.contigBegin + me.contigEnd;
}

template <typename TSpec>
inline typename Member<Match<TSpec>, MatchErrors>::Type
getValue(Match<TSpec> & me, MatchErrors)
{
    return me.errors;
}

template <typename TSpec>
inline typename Member<Match<TSpec>, MatchErrors>::Type
getValue(Match<TSpec> const & me, MatchErrors)
{
    return me.errors;
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

// ----------------------------------------------------------------------------
// Match Composite Getters
// ----------------------------------------------------------------------------

template <typename TReadSeqs, typename TSpec>
inline typename Size<TReadSeqs>::Type
getReadSeqId(Match<TSpec> const & me, TReadSeqs const & readSeqs)
{
    return onForwardStrand(me) ? getFirstMateFwdSeqId(readSeqs, me.readId) : getFirstMateRevSeqId(readSeqs, me.readId);
}

template <typename TSpec>
inline bool isInvalid(Match<TSpec> const & me)
{
    return getValue(me, SortBeginPos()) == getValue(me, SortEndPos());
}

template <typename TSpec>
inline bool isValid(Match<TSpec> const & me)
{
    return !isInvalid(me);
}

// ----------------------------------------------------------------------------
// Function getValue(, MatchErrors())
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
    if (getValue(a, SortBeginPos()) < getValue(b, SortBeginPos()))
        return getValue(b, SortEndPos()) - getValue(a, SortBeginPos());
    else
        return getValue(a, SortEndPos()) - getValue(b, SortBeginPos());
}

// ----------------------------------------------------------------------------
// Function cigarLength()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned getCigarLength(Match<TSpec> const & me)
{
    return isInvalid(me) ? 0 : 2 * getValue(me, MatchErrors()) + 1;
}

// ----------------------------------------------------------------------------
// Function getSortKey(SortBeginPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned long getSortKey(Match<TSpec> const & me, SortBeginPos)
{
    return ((unsigned long)getValue(me, SortContigId())     << (1 + MemberBits<Match<TSpec>, ContigSize>::VALUE + MemberBits<Match<TSpec>, MatchErrors>::VALUE)) |
           ((unsigned long)onReverseStrand(me) << (MemberBits<Match<TSpec>, ContigSize>::VALUE + MemberBits<Match<TSpec>, MatchErrors>::VALUE))     |
           ((unsigned long)getValue(me, SortBeginPos())  <<  MemberBits<Match<TSpec>, MatchErrors>::VALUE)                                                    |
           ((unsigned long)getValue(me, MatchErrors()));
}

// ----------------------------------------------------------------------------
// Function getSortKey(SortEndPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned long getSortKey(Match<TSpec> const & me, SortEndPos)
{
    return ((unsigned long)getValue(me, SortContigId())     << (1 + MemberBits<Match<TSpec>, ContigSize>::VALUE + MemberBits<Match<TSpec>, MatchErrors>::VALUE)) |
           ((unsigned long)onReverseStrand(me) << (MemberBits<Match<TSpec>, ContigSize>::VALUE + MemberBits<Match<TSpec>, MatchErrors>::VALUE))     |
           ((unsigned long)getValue(me, SortEndPos())    <<  MemberBits<Match<TSpec>, MatchErrors>::VALUE)                                                    |
           ((unsigned long)getValue(me, MatchErrors()));
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
    return getValue(a, SortContigId()) == getValue(b, SortContigId()) && strandEqual(a, b);
}

// ----------------------------------------------------------------------------
// Function isDuplicate(SortBeginPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, SortBeginPos)
{
    return contigEqual(a, b) && getValue(a, SortBeginPos()) == getValue(b, SortBeginPos());
}

// ----------------------------------------------------------------------------
// Function isDuplicate(SortEndPos)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, SortEndPos)
{
    return contigEqual(a, b) && getValue(a, SortEndPos()) == getValue(b, SortEndPos());
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

    for (TIter it = itBegin; it != itEnd && getValue(*it, MatchErrors()) <= getValue(*itBegin, MatchErrors()); it++, count++) ;

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
    SEQAN_ASSERT_GT(length(stringSetLimits(matchesSet)), 0u);
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
        if (getValue(*leftIt, SortContigId()) < getValue(*rightIt, SortContigId()))
            findNextContig(leftIt, left, getValue(*leftIt, SortContigId()));
        else if (getValue(*leftIt, SortContigId()) > getValue(*rightIt, SortContigId()))
            findNextContig(rightIt, right, getValue(*rightIt, SortContigId()));
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
    while (!atEnd(it, matches) && getValue(*it, SortContigId()) <= contigId) ++it;
}

// ----------------------------------------------------------------------------
// Function findReverseStrand()
// ----------------------------------------------------------------------------
// Find the first match on the reverse strand of the given contigId.

template <typename TMatchesIterator, typename TMatches, typename TContigId>
inline void findReverseStrand(TMatchesIterator & it, TMatches const & matches, TContigId contigId)
{
    while (!atEnd(it, matches) && (getValue(*it, SortContigId()) <= contigId) && onForwardStrand(*it)) ++it;
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
        unsigned contigId = getValue(*leftIt, SortContigId());

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
