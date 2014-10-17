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


namespace seqan {

// ============================================================================
// Functors
// ============================================================================

// ----------------------------------------------------------------------------
// Functor Getter
// ----------------------------------------------------------------------------

template <typename TObject, typename TTag>
struct Getter
{
    inline unsigned operator()(TObject const & me) const
    {
        return getMember(me, TTag());
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction MemberBits
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec>
struct MemberBits
{
    static const unsigned VALUE = BitsPerValue<typename Member<TObject, TSpec>::Type>::VALUE;
};

// ----------------------------------------------------------------------------
// Metafunction MemberLimits
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec>
struct MemberLimits
{
    static const unsigned VALUE = Power<2, MemberBits<TObject, TSpec>::VALUE>::VALUE - 1;
};

}


using namespace seqan;

// ============================================================================
// Tags
// ============================================================================

struct ReadId_;
typedef Tag<ReadId_> const ReadId;

struct ContigId_;
typedef Tag<ContigId_> const ContigId;

struct ReadSize_;
typedef Tag<ReadSize_> const ReadSize;

struct ContigSize_;
typedef Tag<ContigSize_> const ContigSize;

struct Errors_;
typedef Tag<Errors_> const Errors;

typedef ContigSize      ContigBegin;
typedef ReadSize        ContigEnd;

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
    typename Member<Match, ReadId>::Type        readId       : MemberBits<Match, ReadId>::VALUE;
    typename Member<Match, ContigId>::Type      contigId; // : MemberBits<Match, ContigId>::VALUE;
    bool                                        isRev        : 1;
    typename Member<Match, ContigBegin>::Type   contigBegin  : MemberBits<Match, ContigSize>::VALUE;
    typename Member<Match, ContigEnd>::Type     contigEnd    : MemberBits<Match, ReadSize>::VALUE;
    typename Member<Match, Errors>::Type        errors       : MemberBits<Match, Errors>::VALUE;
}
#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
;

#ifdef PLATFORM_WINDOWS
      #pragma pack(pop)
#endif

// ============================================================================
// Match Types
// ============================================================================

// ----------------------------------------------------------------------------
// Member Types
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec>
struct Member<Match<TSpec>, ReadId>
{
    typedef __uint32    Type;
};

template <typename TSpec>
struct Member<Match<TSpec>, ContigId>
{
    typedef __uint8     Type;
};

template <typename TSpec>
struct Member<Match<TSpec>, ContigSize>
{
    typedef __uint32    Type;
};

template <typename TSpec>
struct Member<Match<TSpec>, ReadSize>
{
    typedef __uint16    Type;
};

template <typename TSpec>
struct Member<Match<TSpec>, Errors>
{
    typedef __uint32    Type;
};
}

// ----------------------------------------------------------------------------
// Member Bits
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TSpec>
struct MemberBits<Match<TSpec>, ReadId>
{
    static const unsigned VALUE = 21;
};

template <typename TSpec>
struct MemberBits<Match<TSpec>, ContigId>
{
    static const unsigned VALUE = 8;
};

template <typename TSpec>
struct MemberBits<Match<TSpec>, ContigSize>
{
    static const unsigned VALUE = 30;
};

template <typename TSpec>
struct MemberBits<Match<TSpec>, ReadSize>
{
    static const unsigned VALUE = 14;
};

template <typename TSpec>
struct MemberBits<Match<TSpec>, Errors>
{
    static const unsigned VALUE = 6;
};
}

// ----------------------------------------------------------------------------
// Class MatchSorter
// ----------------------------------------------------------------------------

template <typename TMatch, typename TTag>
struct MatchSorter
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getMember(a, TTag()) < getMember(b, TTag());
    }
};

template <typename TMatch>
struct MatchSorter<TMatch, ContigBegin>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getSortKey(a, ContigBegin()) < getSortKey(b, ContigBegin());
    }
};

template <typename TMatch>
struct MatchSorter<TMatch, ContigEnd>
{
    inline bool operator()(TMatch const & a, TMatch const & b) const
    {
        return getSortKey(a, ContigEnd()) < getSortKey(b, ContigEnd());
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

template <typename TSpec, typename TContigBegin, typename TContigEnd>
inline void setContigPosition(Match<TSpec> & me, TContigBegin contigBegin, TContigEnd contigEnd)
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
    me.errors = MemberLimits<Match<TSpec>, Errors>::VALUE;
}

// ----------------------------------------------------------------------------
// Function getMember()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TTag>
inline typename Member<Match<TSpec>, TTag>::Type
getMember(Match<TSpec> & me, TTag const & tag)
{
    return getMember(static_cast<Match<TSpec> const &>(me), tag);
}

template <typename TSpec>
inline typename Member<Match<TSpec>, ReadId>::Type
getMember(Match<TSpec> const & me, ReadId)
{
    return me.readId;
}

template <typename TSpec>
inline typename Member<Match<TSpec>, ContigId>::Type
getMember(Match<TSpec> const & me, ContigId)
{
    return me.contigId;
}

template <typename TSpec>
inline typename Member<Match<TSpec>, ContigSize>::Type
getMember(Match<TSpec> const & me, ContigBegin)
{
    return me.contigBegin;
}

template <typename TSpec>
inline typename Member<Match<TSpec>, ContigSize>::Type
getMember(Match<TSpec> const & me, ContigEnd)
{
    return me.contigBegin + me.contigEnd;
}

template <typename TSpec>
inline typename Member<Match<TSpec>, Errors>::Type
getMember(Match<TSpec> const & me, Errors)
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
// Composite Getters
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
    return getMember(me, ContigBegin()) == getMember(me, ContigEnd());
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
inline typename Member<Match<TSpec>, Errors>::Type
getErrors(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return (typename Member<Match<TSpec>, Errors>::Type)getMember(a, Errors()) + getMember(b, Errors());
}

// ----------------------------------------------------------------------------
// Function getTemplateLength()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline typename Member<Match<TSpec>, ContigSize>::Type
getTemplateLength(Match<TSpec> const & a, Match<TSpec> const & b)
{
    typedef typename Member<Match<TSpec>, ContigSize>::Type TContigSize;

    if (getMember(a, ContigId()) != getMember(b, ContigId()))
        return MaxValue<TContigSize>::VALUE;

    if (getMember(a, ContigBegin()) < getMember(b, ContigBegin()))
        return getMember(b, ContigEnd()) - getMember(a, ContigBegin());
    else
        return getMember(a, ContigEnd()) - getMember(b, ContigBegin());
}

// ----------------------------------------------------------------------------
// Function getTemplateDeviation()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TSize>
inline typename Member<Match<TSpec>, ContigSize>::Type
getTemplateDeviation(Match<TSpec> const & a, Match<TSpec> const & b, TSize expectedLength)
{
    typedef typename Member<Match<TSpec>, ContigSize>::Type TContigSize;
    typedef typename MakeSigned<TContigSize>::Type          TSignedContigSize;

    if (isValid(a) && isValid(a))
        return _abs((TSignedContigSize)getTemplateLength(a, b) - (TSignedContigSize)expectedLength);
    else
        return MaxValue<TContigSize>::VALUE;
}

// ----------------------------------------------------------------------------
// Function cigarLength()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned getCigarLength(Match<TSpec> const & me)
{
    return isInvalid(me) ? 0 : 2 * getMember(me, Errors()) + 1;
}

// ----------------------------------------------------------------------------
// Function getSortKey(ContigBegin)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline __uint64 getSortKey(Match<TSpec> const & me, ContigBegin)
{
    typedef Match<TSpec>    TMatch;

    return ((__uint64)getMember(me, ContigId())      << (1 + MemberBits<TMatch, ContigSize>::VALUE + MemberBits<TMatch, Errors>::VALUE)) |
           ((__uint64)onReverseStrand(me)            << (MemberBits<TMatch, ContigSize>::VALUE + MemberBits<TMatch, Errors>::VALUE))     |
           ((__uint64)getMember(me, ContigBegin())   <<  MemberBits<TMatch, Errors>::VALUE)                                              |
           ((__uint64)getMember(me, Errors()));
}

// ----------------------------------------------------------------------------
// Function getSortKey(ContigEnd)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline __uint64 getSortKey(Match<TSpec> const & me, ContigEnd)
{
    typedef Match<TSpec>    TMatch;

    return ((__uint64)getMember(me, ContigId())     << (1 + MemberBits<TMatch, ContigSize>::VALUE + MemberBits<TMatch, Errors>::VALUE)) |
           ((__uint64)onReverseStrand(me)           << (MemberBits<TMatch, ContigSize>::VALUE + MemberBits<TMatch, Errors>::VALUE))     |
           ((__uint64)getMember(me, ContigEnd())    <<  MemberBits<TMatch, Errors>::VALUE)                                              |
           ((__uint64)getMember(me, Errors()));
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
    return getMember(a, ContigId()) == getMember(b, ContigId()) && strandEqual(a, b);
}

// ----------------------------------------------------------------------------
// Function isDuplicate(ContigBegin)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, ContigBegin)
{
    return contigEqual(a, b) && getMember(a, ContigBegin()) == getMember(b, ContigBegin());
}

// ----------------------------------------------------------------------------
// Function isDuplicate(ContigEnd)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, ContigEnd)
{
    return contigEqual(a, b) && getMember(a, ContigEnd()) == getMember(b, ContigEnd());
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
    iterate(matchesSet, MatchesCompactor<TLimits, ContigEnd>(newLimits), Rooted(), threading);

    // Exclude duplicate matches at the end.
    assign(stringSetLimits(matchesSet), newLimits);
    _refreshStringSetLimits(matchesSet, threading);

    // Sort matches by begin position and move unique matches at the beginning.
    iterate(matchesSet, MatchesCompactor<TLimits, ContigBegin>(newLimits), Rooted(), threading);

    // Exclude duplicate matches at the end.
    assign(stringSetLimits(matchesSet), newLimits);
    _refreshStringSetLimits(matchesSet, threading);
}

// ----------------------------------------------------------------------------
// Function countMatchesInStrata()
// ----------------------------------------------------------------------------
// Count the number of matches within the first strata - ordering by errors is required.

template <typename TMatches, typename TStrata>
inline typename Size<TMatches>::Type
countMatchesInStrata(TMatches const & matches, TStrata strata)
{
    typedef typename Iterator<TMatches const, Standard>::Type   TIter;
    typedef typename Size<TMatches>::Type                       TCount;

    TIter itBegin = begin(matches, Standard());
    TIter itEnd = end(matches, Standard());

    TCount count = 0;

    for (TIter it = itBegin; it != itEnd && getMember(*it, Errors()) <= getMember(*itBegin, Errors()) + strata; it++, count++) ;

    return count;
}

// ----------------------------------------------------------------------------
// Function countMatchesInBestStratum()
// ----------------------------------------------------------------------------

template <typename TMatches>
inline typename Size<TMatches>::Type
countMatchesInBestStratum(TMatches const & matches)
{
    return countMatchesInStrata(matches, 0u);
}

// ----------------------------------------------------------------------------
// Function clipMatches()
// ----------------------------------------------------------------------------
// Clip trailing matches from a set of matches.

template <typename TMatchesSet, typename TClipper, typename TThreading>
inline void clipMatches(TMatchesSet & matchesSet, TClipper clipper, TThreading const & threading)
{
    typedef typename StringSetLimits<TMatchesSet>::Type         TLimits;
    typedef typename Suffix<TLimits>::Type                      TLimitsSuffix;

    TLimits newLimits;
    resize(newLimits, length(stringSetLimits(matchesSet)), Exact());
    SEQAN_ASSERT_GT(length(stringSetLimits(matchesSet)), 0u);
    front(newLimits) = 0;

    // Count leading matches to preserve.
    TLimitsSuffix counts = suffix(newLimits, 1);
    transform(counts, matchesSet, clipper, threading);

    // Clip trailing matches.
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

    for (; it != itEnd && !isDuplicate(*it, match, ContigBegin()); it++) ;

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
//    sort(matches, MatchSorter<TMatch, TKey>());
//}

template <typename TIterator, typename TKey>
inline void sortMatches(TIterator & it)
{
    typedef typename Value<TIterator>::Type TMatches;
    typedef typename Value<TMatches>::Type  TMatch;

    TMatches matches = value(it);
    sort(matches, MatchSorter<TMatch, TKey>());
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
        if (getMember(*leftIt, ContigId()) < getMember(*rightIt, ContigId()))
            findNextContig(leftIt, left, getMember(*leftIt, ContigId()));
        else if (getMember(*leftIt, ContigId()) > getMember(*rightIt, ContigId()))
            findNextContig(rightIt, right, getMember(*rightIt, ContigId()));
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
    while (!atEnd(it, matches) && getMember(*it, ContigId()) <= contigId) ++it;
}

// ----------------------------------------------------------------------------
// Function findReverseStrand()
// ----------------------------------------------------------------------------
// Find the first match on the reverse strand of the given contigId.

template <typename TMatchesIterator, typename TMatches, typename TContigId>
inline void findReverseStrand(TMatchesIterator & it, TMatches const & matches, TContigId contigId)
{
    while (!atEnd(it, matches) && (getMember(*it, ContigId()) <= contigId) && onForwardStrand(*it)) ++it;
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
        unsigned contigId = getMember(*leftIt, ContigId());

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
