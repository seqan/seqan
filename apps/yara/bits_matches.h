// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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

// Forward declaration only.
template <typename TObject, typename TSpec>
struct MemberBits;

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

#pragma pack(push,1)
template <typename TSpec = void>
struct Match
{
    typename Member<Match, ReadId>::Type        readId       : MemberBits<Match, ReadId>::VALUE;
    typename Member<Match, ContigId>::Type      contigId     : MemberBits<Match, ContigId>::VALUE;
    bool                                        isRev        : 1;
    typename Member<Match, ContigBegin>::Type   contigBegin  : MemberBits<Match, ContigSize>::VALUE;
    typename Member<Match, ContigEnd>::Type     contigEnd    : MemberBits<Match, ReadSize>::VALUE;
    typename Member<Match, Errors>::Type        errors       : MemberBits<Match, Errors>::VALUE;
};
#pragma pack(pop)

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
    typedef uint32_t    Type;
};

template <typename TContigsSize, typename TContigsLen, typename TContigsSum>
struct Member<Match<Limits<TContigsSize, TContigsLen, TContigsSum> >, ContigId>
{
    // To remove GCC packed-bitfield-compat warning. See MemberBits below.
    typedef uint32_t  Type;
};

template <typename TContigsSize, typename TContigsLen, typename TContigsSum>
struct Member<Match<Limits<TContigsSize, TContigsLen, TContigsSum> >, ContigSize>
{
    typedef TContigsLen  Type;
};

template <typename TSpec>
struct Member<Match<TSpec>, ReadSize>
{
    typedef uint16_t    Type;
};

template <typename TSpec>
struct Member<Match<TSpec>, Errors>
{
    typedef uint32_t    Type;
};

// ----------------------------------------------------------------------------
// Member Bits
// ----------------------------------------------------------------------------

template <typename TObject, typename TSpec>
struct MemberBits
{
    static const unsigned VALUE = BitsPerValue<typename Member<TObject, TSpec>::Type>::VALUE;
};

template <typename TSpec>
struct MemberBits<Match<TSpec>, ReadId>
{
    static const unsigned VALUE = 21;
};

template <typename TContigsLen, typename TContigsSum>
struct MemberBits<Match<Limits<uint8_t, TContigsLen, TContigsSum> >, ContigId>
{
    // To remove GCC packed-bitfield-compat warning.
    static const unsigned VALUE = 8;
};

template <typename TContigsLen, typename TContigsSum>
struct MemberBits<Match<Limits<uint16_t, TContigsLen, TContigsSum> >, ContigId>
{
    // To remove GCC packed-bitfield-compat warning.
    static const unsigned VALUE = 16;
};

template <typename TContigsLen, typename TContigsSum>
struct MemberBits<Match<Limits<uint32_t, TContigsLen, TContigsSum> >, ContigId>
{
    static const unsigned VALUE = 30;
};

template <typename TContigsSize, typename TContigsSum>
struct MemberBits<Match<Limits<TContigsSize, uint64_t, TContigsSum> >, ContigSize>
{
    static const unsigned VALUE = 48;
};

template <typename TSpec>
struct MemberBits<Match<TSpec>, ReadSize>
{
    static const unsigned VALUE = 13;
};

template <typename TSpec>
struct MemberBits<Match<TSpec>, Errors>
{
    static const unsigned VALUE = 7;
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

// ============================================================================
// Match Setters
// ============================================================================

// ----------------------------------------------------------------------------
// Function setReadId()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TReadSeqs, typename TReadSeqId>
inline void setReadId(Match<TSpec> & me, TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    me.readId = getReadId(readSeqs, readSeqId);
    me.isRev = isRevReadSeq(readSeqs, readSeqId);
}

template <typename TSpec, typename TReadId>
inline void setReadId(Match<TSpec> & me, TReadId readId)
{
    me.readId = readId;
}

// ----------------------------------------------------------------------------
// Function setContigPosition()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TContigBegin, typename TContigEnd>
inline void setContigPosition(Match<TSpec> & me, TContigBegin contigBegin, TContigEnd contigEnd)
{
    SEQAN_ASSERT_EQ(getValueI1(contigBegin), getValueI1(contigEnd));
    SEQAN_ASSERT_LT(getValueI2(contigBegin), getValueI2(contigEnd));

    me.contigId = getValueI1(contigBegin);
    me.contigBegin = getValueI2(contigBegin);
    me.contigEnd = getValueI2(contigEnd) - getValueI2(contigBegin);
}

// ----------------------------------------------------------------------------
// Function addContigPosition()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TDelta, typename TContigSeqs>
inline void addContigPosition(Match<TSpec> & me, TDelta delta, TContigSeqs const & contigSeqs)
{
    typedef typename Member<Match<TSpec>, ContigSize>::Type TContigSize;

    SEQAN_ASSERT_GEQ(delta, 0u);
    TContigSize contigLength = length(contigSeqs[getMember(me, ContigId())]);
    me.contigBegin = (me.contigBegin + me.contigEnd + delta < contigLength) ?
                      me.contigBegin + delta : contigLength - me.contigEnd;
}

// ----------------------------------------------------------------------------
// Function subContigPosition()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TDelta>
inline void subContigPosition(Match<TSpec> & me, TDelta delta)
{
    SEQAN_ASSERT_GEQ(delta, 0u);
    me.contigBegin = (me.contigBegin > delta) ? me.contigBegin - delta : 0;
}

// ----------------------------------------------------------------------------
// Function flipStrand()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void flipStrand(Match<TSpec> & me)
{
    me.isRev = !me.isRev;
}

// ----------------------------------------------------------------------------
// Function setInvalid()
// ----------------------------------------------------------------------------

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

// ============================================================================
// Match Getters
// ============================================================================

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

// ----------------------------------------------------------------------------
// Function getCigarLength()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline unsigned getCigarLength(Match<TSpec> const & me)
{
    return isInvalid(me) ? 0 : 2 * getMember(me, Errors()) + 1;
}

// ----------------------------------------------------------------------------
// Function onReverseStrand()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool onReverseStrand(Match<TSpec> const & me)
{
    return me.isRev;
}

// ----------------------------------------------------------------------------
// Function onForwardStrand()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool onForwardStrand(Match<TSpec> const & me)
{
    return !onReverseStrand(me);
}

// ----------------------------------------------------------------------------
// Function isInvalid()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isInvalid(Match<TSpec> const & me)
{
    return getMember(me, ContigBegin()) == getMember(me, ContigEnd());
}

// ----------------------------------------------------------------------------
// Function isValid()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isValid(Match<TSpec> const & me)
{
    return !isInvalid(me);
}

// ----------------------------------------------------------------------------
// Function getReadSeqId()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getReadSeqId(Match<TSpec> const & me, TReadSeqs const & readSeqs)
{
    return onForwardStrand(me) ? getFirstMateFwdSeqId(readSeqs, me.readId) : getFirstMateRevSeqId(readSeqs, me.readId);
}

// ----------------------------------------------------------------------------
// Function getReadLength()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TReadSeqs>
inline typename Member<Match<TSpec>, ReadSize>::Type
getReadLength(Match<TSpec> const & me, TReadSeqs const & readSeqs)
{
    return length(readSeqs[getMember(me, ReadId())]);
}

// ----------------------------------------------------------------------------
// Function getErrorRate()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TReadSeqs>
inline float getErrorRate(Match<TSpec> const & me, TReadSeqs const & readSeqs)
{
    return (float)getMember(me, Errors()) / getReadLength(me, readSeqs);
}

// ----------------------------------------------------------------------------
// Function getSortKey(ContigBegin)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline uint64_t getSortKey(Match<TSpec> const & me, ContigBegin)
{
    typedef Match<TSpec>    TMatch;

    return ((uint64_t)getMember(me, ContigId())      << (1 + MemberBits<TMatch, ContigSize>::VALUE + MemberBits<TMatch, Errors>::VALUE)) |
           ((uint64_t)onReverseStrand(me)            << (MemberBits<TMatch, ContigSize>::VALUE + MemberBits<TMatch, Errors>::VALUE))     |
           ((uint64_t)getMember(me, ContigBegin())   <<  MemberBits<TMatch, Errors>::VALUE)                                              |
           ((uint64_t)getMember(me, Errors()));
}

// ----------------------------------------------------------------------------
// Function getSortKey(ContigEnd)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline uint64_t getSortKey(Match<TSpec> const & me, ContigEnd)
{
    typedef Match<TSpec>    TMatch;

    return ((uint64_t)getMember(me, ContigId())     << (1 + MemberBits<TMatch, ContigSize>::VALUE + MemberBits<TMatch, Errors>::VALUE)) |
           ((uint64_t)onReverseStrand(me)           << (MemberBits<TMatch, ContigSize>::VALUE + MemberBits<TMatch, Errors>::VALUE))     |
           ((uint64_t)getMember(me, ContigEnd())    <<  MemberBits<TMatch, Errors>::VALUE)                                              |
           ((uint64_t)getMember(me, Errors()));
}

// ============================================================================
// Match Pair Getters
// ============================================================================

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
    return getMember(a, ContigId()) == getMember(b, ContigId());
}

// ----------------------------------------------------------------------------
// Function isDuplicate(ContigBegin)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, ContigBegin)
{
    return contigEqual(a, b) && strandEqual(a, b) && getMember(a, ContigBegin()) == getMember(b, ContigBegin());
}

// ----------------------------------------------------------------------------
// Function isDuplicate(ContigEnd)
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isDuplicate(Match<TSpec> const & a, Match<TSpec> const & b, ContigEnd)
{
    return contigEqual(a, b) && strandEqual(a, b) && getMember(a, ContigEnd()) == getMember(b, ContigEnd());
}

// ----------------------------------------------------------------------------
// Function isEqual()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool isEqual(Match<TSpec> const & a, Match<TSpec> const & b)
{
    return getMember(a, ReadId()) == getMember(b, ReadId()) &&
           contigEqual(a, b) && strandEqual(a, b) &&
           getMember(a, ContigBegin()) == getMember(b, ContigBegin()) &&
           getMember(a, ContigEnd()) == getMember(b, ContigEnd()) &&
           getMember(a, Errors()) == getMember(b, Errors());
}

// ----------------------------------------------------------------------------
// Function getLibraryLength()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline typename Member<Match<TSpec>, ContigSize>::Type
getLibraryLength(Match<TSpec> const & a, Match<TSpec> const & b)
{
    typedef typename Member<Match<TSpec>, ContigSize>::Type TContigSize;

    if (isValid(a) && isValid(b) && contigEqual(a, b))
    {
        if (getMember(b, ContigEnd()) > getMember(a, ContigBegin()))
            return getMember(b, ContigEnd()) - getMember(a, ContigBegin());
        else
            return getMember(a, ContigEnd()) - getMember(b, ContigBegin());
    }
    else
    {
        return std::numeric_limits<TContigSize>::max();
    }
}

// ----------------------------------------------------------------------------
// Function getLibraryDeviation()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TSize>
inline typename Member<Match<TSpec>, ContigSize>::Type
getLibraryDeviation(Match<TSpec> const & a, Match<TSpec> const & b, TSize meanLength)
{
    typedef typename Member<Match<TSpec>, ContigSize>::Type TContigSize;
    typedef typename MakeSigned<TContigSize>::Type          TSignedContigSize;

    if (isValid(a) && isValid(b) && contigEqual(a, b))
        return _abs((TSignedContigSize)getLibraryLength(a, b) - (TSignedContigSize)meanLength);
    else
        return std::numeric_limits<TContigSize>::max();
}

// ----------------------------------------------------------------------------
// Function orientationProper()
// ----------------------------------------------------------------------------
// Check orientation --> ... <--

template <typename TSpec>
inline bool orientationProper(Match<TSpec> const & one, Match<TSpec> const & two)
{
    bool oneBeforeTwo = getMember(one, ContigBegin()) < getMember(two, ContigBegin());

    return ((onForwardStrand(one) && onReverseStrand(two) && oneBeforeTwo) ||
            (onForwardStrand(two) && onReverseStrand(one) && !oneBeforeTwo));
}

// ----------------------------------------------------------------------------
// Function isProper()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TMean, typename TStdDev>
inline bool isProper(Match<TSpec> const & one, Match<TSpec> const & two, TMean mean, TStdDev stdDev)
{
    return orientationProper(one, two) &&
           getLibraryDeviation(one, two, mean) < 6 * stdDev;
}

// ============================================================================
// Functions
// ============================================================================

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

template <typename TMatches, typename TKey>
inline void sortMatches(TMatches && matches)
{
    typedef typename Value<TMatches>::Type  TMatch;

    sort(matches, MatchSorter<TMatch, TKey>());
}

// ----------------------------------------------------------------------------
// Function findProperMates()
// ----------------------------------------------------------------------------

template <typename TMatches, typename TMatch, typename TReadSeqs, typename TContigSeqs, typename TNumber1, typename TNumber2>
inline typename Infix<TMatches const>::Type
findProperMates(TMatches const & mates, TMatch const & match,
                TReadSeqs const & readSeqs, TContigSeqs const & contigSeqs,
                TNumber1 mean, TNumber2 stdDev)
{
    typedef typename Iterator<TMatches const>::Type TIter;
    typedef typename Size<TReadSeqs>::Type          TReadId;
    typedef typename Value<TReadSeqs const>::Type   TReadSeq;
    typedef typename Size<TReadSeq>::Type           TReadSeqSize;
    typedef typename MakeSigned<TReadSeqSize>::Type TReadDelta;

    TReadId mateId = getMateId(readSeqs, getMember(match, ReadId()));
    TReadDelta mateLength = length(readSeqs[mateId]);

    // Create lower and upper bound for the mate.
    TMatch mateLeq = match;
    setReadId(mateLeq, mateId);
    flipStrand(mateLeq);
    TMatch mateGeq = mateLeq;
    mateLeq.errors = 0;
    mateGeq.errors = MemberLimits<TMatch, Errors>::VALUE;

    TReadSeqSize deltaMinus = std::max(static_cast<TReadDelta>(0),
                                       static_cast<TReadDelta>(mean) - static_cast<TReadDelta>(6 * stdDev) - mateLength);
    TReadSeqSize deltaPlus = std::max(static_cast<TReadDelta>(0),
                                      static_cast<TReadDelta>(mean) + static_cast<TReadDelta>(6 * stdDev) - mateLength);

    // --> ... mate
    if (onForwardStrand(match))
    {
        addContigPosition(mateLeq, deltaMinus, contigSeqs);
        addContigPosition(mateGeq, deltaPlus, contigSeqs);
    }
    // mate ... <--
    else
    {
        subContigPosition(mateLeq, deltaPlus);
        subContigPosition(mateGeq, deltaMinus);
    }

    TIter first = std::lower_bound(begin(mates, Standard()), end(mates, Standard()), mateLeq, MatchSorter<TMatch, ContigBegin>());
    TIter last = std::upper_bound(begin(mates, Standard()), end(mates, Standard()), mateGeq, MatchSorter<TMatch, ContigEnd>());

    // Return empty infix if no proper mates were found.
    if (first > last) return infix(mates, 0, 0);

    return infix(mates, position(first, mates), position(last, mates));
}

// ----------------------------------------------------------------------------
// Function forAllMatchesPairs()
// ----------------------------------------------------------------------------

template <typename TMatchesSet, typename TReadSeqs, typename TBinaryFunction, typename TThreading>
inline void forAllMatchesPairs(TMatchesSet const & matchesSet, TReadSeqs const & readSeqs, TBinaryFunction && func, TThreading)
{
    typedef typename Size<TReadSeqs>::Type             TReadId;
    typedef typename Value<TMatchesSet const>::Type    TMatches;

    forEach(seqan::Range<TReadId>(0, getPairsCount(readSeqs)), [&](TReadId pairId)
    {
        TReadId firstId = getFirstMateFwdSeqId(readSeqs, pairId);
        TReadId secondId = getSecondMateFwdSeqId(readSeqs, pairId);

        TMatches const & firstMatches = matchesSet[firstId];
        TMatches const & secondMatches = matchesSet[secondId];

        if (!empty(firstMatches) && !empty(secondMatches))
            func(firstMatches, secondMatches);
    },
    TThreading());
}

// ----------------------------------------------------------------------------
// Function getResidualWeight()
// ----------------------------------------------------------------------------
// Residual weight of unseen strata.

template <typename TErrorRate>
inline double getResidualWeight(TErrorRate errorRate)
{
    double p = std::pow(10.0, std::min(2.0, 100.0 * errorRate - 7.0));
    return p / (1.0 - p);
}

// ----------------------------------------------------------------------------
// Function getMatchWeight()
// ----------------------------------------------------------------------------
// Weight of one match.

template <typename TErrorRate>
inline double getMatchWeight(TErrorRate errorRate, TErrorRate optimalRate)
{
    return (1.0 - errorRate) / std::pow(10.0, 300.0 * (errorRate - optimalRate));
}

// ----------------------------------------------------------------------------
// Function getStratumWeight()
// ----------------------------------------------------------------------------
// Weight of one stratum.

template <typename TErrorRate, typename TCount>
inline double getStratumWeight(TErrorRate errorRate, TErrorRate optimalRate, TCount stratumCount)
{
    return stratumCount * getMatchWeight(errorRate, optimalRate);
}

// ----------------------------------------------------------------------------
// Function getFirstTwoStrataWeight()
// ----------------------------------------------------------------------------
// Weight of first two strata.

template <typename TErrorRate, typename TCount>
inline double getFirstTwoStrataWeight(TErrorRate optimalRate, TCount optimalCount, TCount subCount)
{
    return getStratumWeight(optimalRate, optimalRate, optimalCount) +
           getStratumWeight(optimalRate + 0.01, optimalRate, subCount) +
           getResidualWeight(optimalRate);
}

// ----------------------------------------------------------------------------
// Function getMatchProb()
// ----------------------------------------------------------------------------
// Single-end match.

template <typename TErrorRate, typename TCount>
inline double getMatchProb(TErrorRate errorRate, TErrorRate optimalRate, TCount optimalCount, TCount subCount)
{
    return getMatchWeight(errorRate, optimalRate) / getFirstTwoStrataWeight(optimalRate, optimalCount, subCount);
}

// ----------------------------------------------------------------------------
// Function getLibraryProb()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TMean, typename TStdDev>
inline double getLibraryProb(Match<TSpec> const & one, Match<TSpec> const & two, TMean mean, TStdDev stdDev)
{
    if (!isProper(one, two, mean, stdDev)) return 0.0009; // 0.001 * 0.90

    double libraryDev = getLibraryDeviation(one, two, mean);
    double libraryScore = (double)libraryDev / stdDev;

    // 1.0 - 2.0 * zScore(libraryScore)
    static const double SQRT_2 = 1.41421356237;
    return std::max(0.001, std::erfc((double)libraryScore / SQRT_2));
}

// ----------------------------------------------------------------------------
// Function findPrimaryMatch()
// ----------------------------------------------------------------------------

template <typename TMatches, typename TErrorRate, typename TCount, typename TReadSeqs, typename TContigSeqs,
          typename TMean, typename TStdDev>
inline Pair<typename Iterator<TMatches, Standard>::Type, double>
findPrimaryMatch(TMatches const & firstMatches,
                 TMatches const & secondMatches,
                 TErrorRate firstOptimalRate,
                 TErrorRate secondOptimalRate,
                 TCount secondOptimalCount,
                 TCount secondSubCount,
                 TReadSeqs const & readSeqs,
                 TContigSeqs const & contigSeqs,
                 TMean mean,
                 TStdDev stdDev)
{
    typedef typename Value<TMatches const>::Type                TMatch;
    typedef typename Iterator<TMatches const, Standard>::Type   TMatchesIt;
    typedef Pair<TMatchesIt, double>                            TPair;

    double firstMatchesWeightSum = 0;
    double firstMatchWeightMax = 0;
    TMatchesIt firstMatchBestIt = end(firstMatches, Standard());

    iterate(firstMatches, [&](TMatchesIt firstMatchesIt)
    {
        TMatch const & firstMatch = value(firstMatchesIt);
        double firstMatchWeight = 0;

        TCount secondOptimalImproperCount = secondOptimalCount;
        TCount secondSubImproperCount = secondSubCount;

        TMatches const & mates = findProperMates(secondMatches, firstMatch, readSeqs, contigSeqs, mean, stdDev);

        forEach(mates, [&](TMatch const & secondMatch)
        {
            // Sum weights of proper match.
            double secondErrorRate = getErrorRate(secondMatch, readSeqs);
            firstMatchWeight += getMatchWeight(secondErrorRate, secondOptimalRate) *
                                getLibraryProb(firstMatch, secondMatch, mean, stdDev);

            // Decrease number of improper matches.
            if (secondErrorRate == secondOptimalRate)
                --secondOptimalImproperCount;
            else if (secondErrorRate == secondOptimalRate + 0.01)
                --secondSubImproperCount;
        });

        // Sum weights of improper matches.
        firstMatchWeight += (getStratumWeight(secondOptimalRate, secondOptimalRate, secondOptimalImproperCount) +
                             getStratumWeight(secondOptimalRate + 0.01, secondOptimalRate, secondSubImproperCount) +
                             getResidualWeight(secondOptimalRate)) * 0.0009;

        double firstErrorRate = getErrorRate(firstMatch, readSeqs);
        firstMatchWeight *= getMatchWeight(firstErrorRate, firstOptimalRate);

        // Update best matches.
        if (firstMatchWeight > firstMatchWeightMax)
        {
            firstMatchWeightMax = firstMatchWeight;
            firstMatchBestIt = firstMatchesIt;
        }

        // Update sum of weights.
        firstMatchesWeightSum += firstMatchWeight;
    },
    Standard(), Serial());

    double firstMatchBestProb = firstMatchWeightMax / firstMatchesWeightSum;

    return TPair(firstMatchBestIt, firstMatchBestProb);
}

// ----------------------------------------------------------------------------
// Function write()
// ----------------------------------------------------------------------------
// Debug.

template <typename TStream, typename TSpec>
inline void write(TStream & stream, Match<TSpec> const & me)
{
    stream << getMember(me, ReadId()) << " @ " << (unsigned)getMember(me, ContigId()) << '/' << onReverseStrand(me)
           << " : " << Pair<unsigned>(getMember(me, ContigBegin()), getMember(me, ContigEnd()))
           << " + " << getMember(me, Errors()) << '\n';
}

#endif  // #ifndef APP_YARA_BITS_MATCHES_H_
