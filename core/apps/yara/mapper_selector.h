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

#ifndef APP_YARA_MAPPER_SELECTOR_H_
#define APP_YARA_MAPPER_SELECTOR_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class PairsSelector
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TSpec, typename Traits>
struct PairsSelector
{
    typedef typename Traits::TReadsContext     TReadsContext;
    typedef typename Traits::TReadSeqs         TReadSeqs;
    typedef typename Traits::TMatchesSet       TMatchesSet;
    typedef typename Traits::TMatches          TMatches;

    // Shared-memory read-write data.
    TMatches &          pairs;
    TReadsContext &     ctx;

    // Shared-memory read-only data.
    TReadSeqs const &   readSeqs;
    TMatchesSet const & firstMatchesSet;
    TMatchesSet const & secondMatchesSet;
    Options const &     options;

    PairsSelector(TMatches & pairs,
                  TReadsContext & ctx,
                  TReadSeqs const & readSeqs,
                  TMatchesSet const & firstMatchesSet,
                  TMatchesSet const & secondMatchesSet,
                  Options const & options) :
        pairs(pairs),
        ctx(ctx),
        readSeqs(readSeqs),
        firstMatchesSet(firstMatchesSet),
        secondMatchesSet(secondMatchesSet),
        options(options)
    {
        _selectPairsImpl(*this);
    }

    template <typename TIterator>
    void operator() (TIterator const & it)
    {
        _selectPairImpl(*this, it);
    }

    template <typename TMatch, typename TOrientation>
    void operator() (TMatch const & first, TMatch const & second, TOrientation const & tag)
    {
        _enumeratePairsImpl(*this, first, second, tag);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _selectPairsImpl()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
inline void _selectPairsImpl(PairsSelector<TSpec, Traits> & me)
{
    typedef typename Traits::TReadSeqs              TReadSeqs;
    typedef Segment<TReadSeqs const, PrefixSegment> TPrefix;

    TPrefix pairs(me.readSeqs, getPairsCount(me.readSeqs));

    // Iterate over all pairs.
    iterate(pairs, me, Rooted(), typename Traits::TThreading());
}

// ----------------------------------------------------------------------------
// Function _selectPairImpl()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TIterator>
inline void _selectPairImpl(PairsSelector<TSpec, Traits> & me, TIterator const & it)
{
    typedef typename Traits::TReadSeqs                  TReadSeqs;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    // Get pairId.
    TReadId pairId = position(it);

    TReadId firstId = getFirstMateFwdSeqId(me.readSeqs, pairId);
    TReadId secondId = getSecondMateFwdSeqId(me.readSeqs, pairId);

    bucketMatches(me.firstMatchesSet[firstId], me.secondMatchesSet[secondId], me);
}

// ----------------------------------------------------------------------------
// Function _enumeratePairsImpl()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatches>
inline void _enumeratePairsImpl(PairsSelector<TSpec, Traits> & me, TMatches const & first, TMatches const & second, FwdRev)
{
    if (me.options.libraryOrientation == FWD_REV || me.options.libraryOrientation == ANY)
        _enumeratePairs(me, first, second);
}

template <typename TSpec, typename Traits, typename TMatches>
inline void _enumeratePairsImpl(PairsSelector<TSpec, Traits> & me, TMatches const & first, TMatches const & second, RevFwd)
{
    if (me.options.libraryOrientation == FWD_REV || me.options.libraryOrientation == ANY)
        _enumeratePairs(me, second, first);
}

template <typename TSpec, typename Traits, typename TMatches>
inline void _enumeratePairsImpl(PairsSelector<TSpec, Traits> & me, TMatches const & first, TMatches const & second, FwdFwd)
{
    if (me.options.libraryOrientation == FWD_FWD || me.options.libraryOrientation == ANY)
    {
        _enumeratePairs(me, first, second);
        _enumeratePairs(me, second, first);
    }
}

template <typename TSpec, typename Traits, typename TMatches>
inline void _enumeratePairsImpl(PairsSelector<TSpec, Traits> & me, TMatches const & first, TMatches const & second, RevRev)
{
    if (me.options.libraryOrientation == REV_REV || me.options.libraryOrientation == ANY)
    {
        _enumeratePairs(me, first, second);
        _enumeratePairs(me, second, first);
    }
}

// ----------------------------------------------------------------------------
// Function _enumeratePairs()
// ----------------------------------------------------------------------------

//template <typename TSpec, typename Traits, typename TMatches>
//inline void _enumeratePairs(PairsSelector<TSpec, Traits> & me, TMatches const & left, TMatches const & right)
//{
//    typedef typename Iterator<TMatches const, Standard>::Type TIterator;
//
//    TIterator leftBegin = begin(left, Standard());
//    TIterator leftEnd = end(left, Standard());
//    TIterator rightBegin = begin(right, Standard());
//    TIterator rightEnd = end(right, Standard());
//
//    TIterator leftIt = leftBegin;
//
//    if (leftIt == leftEnd) return;
//
//    // Left queue C= right queue, i.e. leftTail >= rightTail && leftHead <= rightHead
//
//    // Get next right match.
//    for (TIterator rightIt = rightBegin; rightIt != rightEnd; ++rightIt)
//    {
//        // Compute the interval of feasible left matches from current right match.
//        unsigned rightHead = getMember(*rightIt, ContigEnd()) - me.options.libraryLength + me.options.libraryError;
//        unsigned rightTail = getMember(*rightIt, ContigEnd()) - me.options.libraryLength - me.options.libraryError;
//
//        // Seek first feasible left match - beyond the right tail.
//        while (leftIt != leftEnd && getMember(*leftIt, ContigBegin()) < rightTail)
//            ++leftIt;
//
//        // No left matches anymore.
//        TIterator leftTailIt = leftIt;
//        if (leftTailIt == leftEnd)
//            break;
//
//        // Continue with next right match if there are no feasible left matches anymore.
//        unsigned leftTail = getMember(*leftTailIt, ContigBegin());
//        if (leftTail >= rightHead)
//            continue;
//
//        // Seek first infeasible left match - beyond the right head.
//        while (leftIt != leftEnd && getMember(*leftIt, ContigBegin()) < rightHead)
//            ++leftIt;
//        TIterator leftHeadIt = leftIt;
//
//        // Couple all lefts matches in the queue with current right match.
//        for (TIterator leftQueueIt = leftTailIt; leftQueueIt != leftHeadIt; ++leftQueueIt)
//            _selectBestPair(me, *leftQueueIt, *rightIt);
//
//        // Empty left queue.
////        if (leftTailIt == leftHeadIt) continue;
//    }
//}

// ----------------------------------------------------------------------------
// Function _enumeratePairs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatches>
inline void _enumeratePairs(PairsSelector<TSpec, Traits> & me, TMatches const & left, TMatches const & right)
{
    typedef typename Iterator<TMatches const, Standard>::Type TIterator;

    TIterator leftIt = begin(left, Standard());
    TIterator leftEnd = end(left, Standard());
    TIterator rightBegin = begin(right, Standard());
    TIterator rightEnd = end(right, Standard());

    for (; leftIt != leftEnd; ++leftIt)
        for (TIterator rightIt = rightBegin; rightIt != rightEnd; ++rightIt)
             if (getTemplateDeviation(*leftIt, *rightIt, me.options.libraryLength) <= me.options.libraryError)
                _selectBestPair(me, *leftIt, *rightIt);
}

// ----------------------------------------------------------------------------
// Function _selectBestPair()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatch>
inline void _selectBestPair(PairsSelector<TSpec, Traits> & me, TMatch const & left, TMatch const & right)
{
    TMatch & bestLeft = me.pairs[getMember(left, ReadId())];
    TMatch & bestRight = me.pairs[getMember(right, ReadId())];

    if (isPaired(me.ctx, getMember(left, ReadId())) && isPaired(me.ctx, getMember(right, ReadId()))) return;

    if (getErrors(left, right) < getErrors(bestLeft, bestRight) ||
        (getErrors(left, right) == getErrors(bestLeft, bestRight) &&
            getTemplateDeviation(left, right, me.options.libraryLength) <=
            getTemplateDeviation(bestLeft, bestRight, me.options.libraryLength)))
    {
        bestLeft = left;
        bestRight = right;
    }
}

#endif  // #ifndef APP_YARA_MAPPER_SELECTOR_H_
