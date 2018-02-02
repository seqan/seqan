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

#ifndef APP_YARA_MAPPER_RANKER_H_
#define APP_YARA_MAPPER_RANKER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeedsRanker
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TSpec, typename Traits>
struct SeedsRanker
{
    typedef typename Traits::THitsCounts       THitsCounts;
    typedef typename Traits::TRanks            TRanks;
    typedef typename Traits::TSeeds            TSeeds;
    typedef typename Traits::THits             THits;
    typedef typename Traits::TReadSeqs         TReadSeqs;

    // Shared-memory read-write data.
    THitsCounts &       counts;
    TRanks &            ranks;

    // Shared-memory read-only data.
    TSeeds const &      seeds;
    THits const &       hits;
    TReadSeqs const &   readSeqs;
    Options const &     options;

    SeedsRanker(THitsCounts & counts,
                TRanks & ranks,
                TSeeds const & seeds,
                THits const & hits,
                Options const & options) :
        counts(counts),
        ranks(ranks),
        seeds(seeds),
        hits(hits),
        readSeqs(host(seeds)),
        options(options)
    {
        _rankAllSeedsImpl(*this);
    }

    template <typename TReadSeqsIterator>
    void operator() (TReadSeqsIterator const & it)
    {
        _countHitsImpl(*this, it);
    }
};

// ----------------------------------------------------------------------------
// Class SeedsSorter
// ----------------------------------------------------------------------------

template <typename THitsCounts, typename TSpec = void>
struct SeedsSorter
{
    THitsCounts const &    counts;

    SeedsSorter(THitsCounts const & counts) :
        counts(counts)
    {}

    template <typename TIterator>
    void operator() (TIterator & it)
    {
        typedef typename Value<TIterator>::Type     TSeedIds;

        TSeedIds const & seedIds = value(it);

        sort(seedIds, KeySorter<THitsCounts>(counts));
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _rankAllSeedsImpl()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
inline void _rankAllSeedsImpl(SeedsRanker<TSpec, Traits> & me)
{
    resize(me.counts, lengthSum(me.seeds), Exact());

    // One element per seed.
    resize(concat(me.ranks), lengthSum(me.seeds), Exact());
    // One bucket per read seq.
    resize(stringSetLimits(me.ranks), length(me.readSeqs) + 1, 0, Exact());

    // Fill counts and ranks.
    iterate(me.readSeqs, me, Standard(), typename Traits::TThreading());

    // Bucket the seeds by read seq.
    partialSum(stringSetLimits(me.ranks), typename Traits::TThreading());

    // Rank seeds of each read seq by number of hits.
    iterate(me.ranks, SeedsSorter<typename Traits::THitsCounts>(me.counts), Rooted(), typename Traits::TThreading());
}

// ----------------------------------------------------------------------------
// Function _countHitsImpl()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TReadSeqsIterator>
inline void _countHitsImpl(SeedsRanker<TSpec, Traits> & me, TReadSeqsIterator const & it)
{
    typedef typename Traits::THits                      THits;
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Size<THit>::Type                   THitSize;
    typedef typename Traits::TSeeds                     TSeeds;
    typedef typename Id<TSeeds>::Type                   TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename Traits::TReadSeqs                  TReadSeqs;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    TReadId readSeqId = position(it, me.readSeqs);

    // Count the number of seeds per read seq.
    TSeedIds readSeedIds = getSeedIds(me.seeds, readSeqId);
    assignValue(stringSetLimits(me.ranks), readSeqId + 1, getValueI2(readSeedIds) - getValueI1(readSeedIds));

    for (TSeedId seedId = getValueI1(readSeedIds); seedId < getValueI2(readSeedIds); ++seedId)
    {
        // Fill the seed ranks with seedIds to be sorted e.g. the identity permutation.
        assignValue(concat(me.ranks), seedId, seedId);

        // Count the number of hits per seed.
        THitIds seedHitIds = getHitIds(me.hits, seedId);
        THitSize seedHits = countHits<THitSize>(me.hits, seedHitIds);
        assignValue(concat(me.counts), seedId, seedHits);
    }

}

#endif  // #ifndef APP_YARA_MAPPER_RANKER_H_
