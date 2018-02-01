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

#ifndef APP_YARA_MAPPER_COLLECTOR_H_
#define APP_YARA_MAPPER_COLLECTOR_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class SeedsCollector
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TSpec, typename Traits>
struct SeedsCollector
{
    typedef typename Traits::TReadsContext      TReadsContext;
    typedef typename Traits::TSeeds             TSeeds;
    typedef typename Traits::TReadSeqs          TReadSeqs;
    typedef typename Traits::TSeedsCount        TSeedsCount;

    // Shared-memory read-write data.
    TReadsContext &     ctx;
    TSeeds &            seeds;
    TSeedsCount &       seedsCount;

    // Shared-memory read-only data.
    unsigned            seedErrors;
    TReadSeqs const &   readSeqs;
    Options const &     options;

    SeedsCollector(TReadsContext & ctx,
                   TSeeds & seeds,
                   TSeedsCount & seedsCount,
                   unsigned seedErrors,
                   TReadSeqs const & readSeqs,
                   Options const & options) :
        ctx(ctx),
        seeds(seeds),
        seedsCount(seedsCount),
        seedErrors(seedErrors),
        readSeqs(readSeqs),
        options(options)
    {
        _init(*this);

        // Iterate over all reads.
        iterate(readSeqs, *this, Rooted(), typename Traits::TThreading());

        _finalize(*this);
    }

    template <typename TReadSeqsIterator>
    void operator() (TReadSeqsIterator const & it)
    {
        _collectSeedsImpl(*this, it);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _init()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
inline void _init(SeedsCollector<TSpec, Traits> & me)
{
    partialSum(me.seedsCount, me.seedsCount, typename Traits::TThreading());

    // Resize space for seeds.
    clear(me.seeds);
    resize(me.seeds, back(me.seedsCount), Exact());
}

template <typename Traits>
inline void _init(SeedsCollector<Counter, Traits> & me)
{
    clear(me.seedsCount);
    resize(me.seedsCount, getReadSeqsCount(me.readSeqs), 0, Exact());
}

// ----------------------------------------------------------------------------
// Function _collectSeedsImpl()
// ----------------------------------------------------------------------------
// Collects seeds from one unseeded read.

template <typename TSpec, typename Traits, typename TReadSeqsIterator>
inline void _collectSeedsImpl(SeedsCollector<TSpec, Traits> & me, TReadSeqsIterator const & it)
{
    typedef typename Traits::TReadSeqs                  TReadSeqs;
    typedef typename Id<TReadSeqs>::Type                TReadSeqId;

    TReadSeqId readSeqId = position(it);
    TReadSeqId readId = getReadId(me.readSeqs, readSeqId);

    if (!isMapped(me.ctx, readId) && getSeedErrors(me.ctx, readSeqId) == me.seedErrors)
        _getSeeds(me, readSeqId);
}

// ----------------------------------------------------------------------------
// Function _getSeeds()
// ----------------------------------------------------------------------------
// Enumerates the seeds for a given read sequence.

template <typename TSpec, typename Traits, typename TReadSeqId>
inline void _getSeeds(SeedsCollector<TSpec, Traits> & me, TReadSeqId readSeqId)
{
    typedef typename Traits::TMatch                         TMatch;
    typedef typename Traits::TReadSeqs                      TReadSeqs;
    typedef typename StringSetPosition<TReadSeqs>::Type     TPos;
    typedef typename Value<TReadSeqs>::Type                 TReadSeq;
    typedef typename Size<TReadSeq>::Type                   TSize;

    TSize readLength = length(me.readSeqs[readSeqId]);
    TSize readErrors = getReadErrors<TMatch>(me.options, readLength);
    TSize seedErrors = getSeedErrors(me.ctx, readSeqId);
    TSize seedsCount = static_cast<TSize>(std::ceil((readErrors + 1) / (seedErrors + 1.0)));
    TSize seedsLength = readLength / seedsCount;

    for (TSize seedId = 0; seedId < seedsCount; ++seedId)
        _addSeed(me, TPos(readSeqId, seedId * seedsLength), seedsLength);
}

// ----------------------------------------------------------------------------
// Function _addSeed()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TPos, typename TSize>
inline void _addSeed(SeedsCollector<TSpec, Traits> & me, TPos seedPos, TSize seedLength)
{
    me.seedsCount[getSeqNo(seedPos)]--;
    assignInfixWithLength(me.seeds, me.seedsCount[getSeqNo(seedPos)], seedPos, seedLength);
}

template <typename Traits, typename TPos, typename TSize>
inline void _addSeed(SeedsCollector<Counter, Traits> & me, TPos seedPos, TSize /* seedLength */)
{
    me.seedsCount[getSeqNo(seedPos)]++;
}

// ----------------------------------------------------------------------------
// Function _finalize()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
inline void _finalize(SeedsCollector<TSpec, Traits> & me)
{
    _refreshStringSetLimits(me.seeds);
}

template <typename Traits>
inline void _finalize(SeedsCollector<Counter, Traits> & /* me */) {}

#endif  // #ifndef APP_YARA_MAPPER_COLLECTOR_H_
