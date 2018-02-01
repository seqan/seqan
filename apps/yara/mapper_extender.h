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

#ifndef APP_YARA_MAPPER_EXTENDER_H_
#define APP_YARA_MAPPER_EXTENDER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class HitsExtender
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TSpec, typename Traits>
struct HitsExtender
{
    typedef typename Traits::TContigSeqs       TContigSeqs;
    typedef typename Traits::TContigsPos       TContigsPos;
    typedef typename Traits::TReadSeqs         TReadSeqs;
    typedef typename Traits::TReadSeq          TReadSeq;
    typedef typename Traits::TReadsContext     TReadsContext;
    typedef typename Traits::TMatchesAppender  TMatches;
    typedef typename Traits::TMatch            TMatch;
    typedef typename Traits::TSeeds            TSeeds;
    typedef typename Traits::THits             THits;
    typedef typename Traits::TRanks            TRanks;
    typedef typename Traits::TSA               TSA;

    typedef Extender<TContigSeqs, TReadSeq, EditDistance> TExtender;

    // Thread-private data.
    TExtender           extender;
    TMatch              prototype;

    // Shared-memory read-write data.
    TReadsContext &     ctx;
    TMatches &          matches;

    // Shared-memory read-only data.
    TContigSeqs const & contigSeqs;
    TReadSeqs &         readSeqs;
    TSeeds const &      seeds;
    THits const &       hits;
    TRanks const &      ranks;
    unsigned            seedErrors;
    TSA const &         sa;
    Options const &     options;

    HitsExtender(TReadsContext & ctx,
                 TMatches & matches,
                 TContigSeqs const & contigSeqs,
                 TSeeds const & seeds,
                 THits const & hits,
                 TRanks const & ranks,
                 unsigned seedErrors,
                 TSA const & sa,
                 Options const & options) :
        extender(contigSeqs),
        prototype(),
        ctx(ctx),
        matches(matches),
        contigSeqs(contigSeqs),
        readSeqs(host(seeds)),
        seeds(seeds),
        hits(hits),
        ranks(ranks),
        seedErrors(seedErrors),
        sa(sa),
        options(options)
    {
        _extendHitsImpl(*this, typename Traits::TStrategy());
    }

    template <typename THitsIterator>
    void operator() (THitsIterator const & hitsIt)
    {
        _extendHitImpl(*this, hitsIt, typename Traits::TStrategy());
    }

    template <typename TMatchPos, typename TErrors>
    void operator() (TMatchPos matchBegin, TMatchPos matchEnd, TErrors errors)
    {
        _addMatchImpl(*this, matchBegin, matchEnd, errors);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _extendHitsImpl()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TStrategy>
inline void _extendHitsImpl(HitsExtender<TSpec, Traits> & me, TStrategy const & /* tag */)
{
    // Iterate over all hits.
    iterate(me.hits, me, Standard(), typename Traits::TThreading());
}

template <typename TSpec, typename Traits>
inline void _extendHitsImpl(HitsExtender<TSpec, Traits> & me, Strata)
{
    typedef typename Traits::TReadSeqs              TReadSeqs;
    typedef Segment<TReadSeqs const, PrefixSegment> TPrefix;

    TPrefix reads(me.readSeqs, getReadsCount(me.readSeqs));

    // Iterate over all reads.
    iterate(reads, me, Rooted(), typename Traits::TThreading());
}

// ----------------------------------------------------------------------------
// Function _extendHitImpl()
// ----------------------------------------------------------------------------
// Extends one hit.

template <typename TSpec, typename Traits, typename THitsIterator, typename TStrategy>
inline void _extendHitImpl(HitsExtender<TSpec, Traits> & me, THitsIterator const & hitsIt, TStrategy const & /* tag */)
{
    typedef typename Traits::TContigsPos                TContigsPos;

    typedef typename Traits::TReadSeqs                  TReadSeqs;
    typedef typename Traits::TReadSeq                   TReadSeq;
    typedef typename Id<TReadSeqs>::Type                TReadId;
    typedef Pair<typename Position<TReadSeq>::Type>     TReadPos;
    typedef typename Size<TReadSeq>::Type               TReadSeqSize;
    typedef typename Size<TReadSeq>::Type               TErrors;

    typedef typename Traits::TSeeds                     TSeeds;
    typedef typename Id<TSeeds>::Type                   TSeedId;

    typedef typename Traits::THits                      THits;
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef typename Position<THit>::Type               THitRange;
    typedef TReadSeqSize                                THitErrors;

    typedef typename Traits::TMatch                     TMatch;

    typedef typename Traits::TSA                        TSA;
    typedef typename Size<TSA>::Type                    TSAPos;
    typedef typename Value<TSA>::Type                   TSAValue;


    // Get hit id.
    THitId hitId = position(hitsIt, me.hits);

    // Extract hit info.
    TSeedId seedId = getSeedId(me.hits, hitId);
    THitRange hitRange = getRange(me.hits, hitId);
    THitErrors hitErrors = getErrors(me.hits, hitId);

    // Get read.
    TReadId readSeqId = getReadSeqId(me.seeds, seedId);
    TReadSeq const & readSeq = me.readSeqs[readSeqId];

    // Fill readSeqId.
    setReadId(me.prototype, me.readSeqs, readSeqId);

    // Get position in read.
    TReadPos readPos = getPosInRead(me.seeds, seedId);
    TReadSeqSize seedLength = getValueI2(readPos) - getValueI1(readPos);

    for (TSAPos saPos = getValueI1(hitRange); saPos < getValueI2(hitRange); ++saPos)
    {
        // Invert SA value.
        TSAValue saValue = me.sa[saPos];
        SEQAN_ASSERT_GEQ(suffixLength(saValue, me.contigSeqs), seedLength);
        if (suffixLength(saValue, me.contigSeqs) < seedLength) continue;
        setSeqOffset(saValue, suffixLength(saValue, me.contigSeqs) - seedLength);

        // Compute position in contig.
        TContigsPos contigBegin = saValue;
        TContigsPos contigEnd = posAdd(contigBegin, seedLength);

        // Get absolute number of errors.
        TErrors maxErrors = getReadErrors<TMatch>(me.options, length(readSeq));

        extend(me.extender,
               readSeq,
               contigBegin, contigEnd,
               getValueI1(readPos), getValueI2(readPos),
               hitErrors, maxErrors,
               me);
    }
}

template <typename TSpec, typename Traits, typename TReadSeqsIterator>
inline void _extendHitImpl(HitsExtender<TSpec, Traits> & me, TReadSeqsIterator const & it, Strata)
{
    typedef typename Traits::TReadSeqs                  TReadSeqs;
    typedef typename Size<TReadSeqs>::Type              TReadId;
    typedef typename Traits::TReadSeq                   TReadSeq;
    typedef typename Size<TReadSeq>::Type               TReadSeqSize;

    typedef typename Traits::TSeeds                     TSeeds;
    typedef typename Id<TSeeds>::Type                   TSeedId;

    typedef typename Traits::THits                      THits;
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Iterator<THits const, Standard>::Type  THitsIt;

    typedef typename Traits::TRanks                     TRanks;
    typedef typename Reference<TRanks const>::Type      TRank;

    typedef typename Traits::TMatch                     TMatch;

    // Get readId.
    TReadId readId = position(it);
    TReadSeqSize readLength = length(me.readSeqs[readId]);
    TReadSeqSize readStrata = getReadStrata<TMatch>(me.options, readLength);

    TReadId fwdSeqId = getFirstMateFwdSeqId(me.readSeqs, readId);
    TReadId revSeqId = getFirstMateRevSeqId(me.readSeqs, readId);

    TRank fwdRank = me.ranks[fwdSeqId];
    TRank revRank = me.ranks[revSeqId];
    SEQAN_ASSERT_EQ(length(fwdRank), length(revRank));

    // TODO(esiragusa): Get hits of fwd and rev read seq.
//    THits fwdHits = getHits(me.hits, fwdSeqId);
//    THits revHits = getHits(me.hits, revSeqId);

    // Iterate fwd and rev seeds by rank.
    TSeedId seedRanks = length(fwdRank);
    for (TSeedId seedRank = 0; seedRank < seedRanks && !isMapped(me.ctx, readId); ++seedRank)
    {
        // Get seedIds by rank.
        TSeedId fwdSeedId = fwdRank[seedRank];
        TSeedId revSeedId = revRank[seedRank];

        // Get hits.
        THitIds fwdHitIds = getHitIds(me.hits, fwdSeedId);
        THitIds revHitIds = getHitIds(me.hits, revSeedId);

        // Verify seed hits.
        THitsIt hitsBegin = begin(me.hits, Standard());
        for (THitsIt it = hitsBegin + getValueI1(fwdHitIds); it != hitsBegin + getValueI2(fwdHitIds); ++it)
            _extendHitImpl(me, it, All());
        for (THitsIt it = hitsBegin + getValueI1(revHitIds); it != hitsBegin + getValueI2(revHitIds); ++it)
            _extendHitImpl(me, it, All());

        // Mark mapped reads.
        if (getMinErrors(me.ctx, readId) + readStrata <= seedRank * (me.seedErrors + 1))
            setMapped(me.ctx, readId);
    }
}

// ----------------------------------------------------------------------------
// Function _addMatchImpl()
// ----------------------------------------------------------------------------
// Adds one match.

template <typename TSpec, typename Traits, typename TMatchPos, typename TMatchErrors>
inline void _addMatchImpl(HitsExtender<TSpec, Traits> & me,
                          TMatchPos matchBegin,
                          TMatchPos matchEnd,
                          TMatchErrors matchErrors)
{
    typedef typename Traits::TReadSeqs       TReadSeqs;
    typedef typename Size<TReadSeqs>::Type   TReadSeqId;

    setContigPosition(me.prototype, matchBegin, matchEnd);
    me.prototype.errors = matchErrors;
    appendValue(me.matches, me.prototype, Generous(), typename Traits::TThreading());

    TReadSeqId readId = getMember(me.prototype, ReadId());
    setMinErrors(me.ctx, readId, matchErrors);
}

#endif  // #ifndef APP_YARA_MAPPER_EXTENDER_H_
