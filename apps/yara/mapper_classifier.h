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

#ifndef APP_YARA_MAPPER_CLASSIFIER_H_
#define APP_YARA_MAPPER_CLASSIFIER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadsClassifier
// ----------------------------------------------------------------------------
// One instance per thread.

template <typename TSpec, typename TConfig>
struct ReadsClassifier
{
    typedef typename TConfig::TReadsContext     TReadsContext;
    typedef typename TConfig::THits             THits;
    typedef typename TConfig::TSeeds            TSeeds;
    typedef typename TConfig::TReadSeqs         TReadSeqs;

    // Shared-memory read-write data.
    TReadsContext &     ctx;
    THits &             hits;

    // Shared-memory read-only data.
    TSeeds const &      seeds;
    TReadSeqs const &   readSeqs;
    Options const &     options;

    ReadsClassifier(TReadsContext & ctx,
                    THits & hits,
                    TSeeds const & seeds,
                    Options const & options) :
        ctx(ctx),
        hits(hits),
        seeds(seeds),
        readSeqs(host(seeds)),
        options(options)
    {
        _classifyReadsImpl(*this, typename TConfig::TStrategy()); //, typename TConfig::TAnchoring());
    }

    template <typename TReadSeqsIterator>
    void operator() (TReadSeqsIterator const & it)
    {
        _classifyReadImpl(*this, it, typename TConfig::TStrategy()); //, typename TConfig::TAnchoring());
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _classifyReadsImpl(); Default
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void _classifyReadsImpl(ReadsClassifier<TSpec, TConfig> & me, All)
{
    // Iterate over all reads.
    iterate(me.readSeqs, me, Rooted(), typename TConfig::TThreading());
}

// ----------------------------------------------------------------------------
// Function _classifyReadsImpl(); AnchorOne
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void _classifyReadsImpl(ReadsClassifier<TSpec, TConfig> & me, Strata)//AnchorOne)
{
    typedef typename TConfig::TReadSeqs             TReadSeqs;
    typedef Segment<TReadSeqs const, PrefixSegment> TPrefix;

    TPrefix pairs(me.readSeqs, getReadsCount(me.readSeqs));

    // Iterate over all pairs.
    iterate(pairs, me, Rooted(), typename TConfig::TThreading());
}

// ----------------------------------------------------------------------------
// Function _classifyReadImpl(); Default
// ----------------------------------------------------------------------------
// Raises the seeds errors, mark for reseeding and clears the hits of hard reads.

template <typename TSpec, typename TConfig, typename TReadSeqsIterator>
inline void _classifyReadImpl(ReadsClassifier<TSpec, TConfig> & me, TReadSeqsIterator const & it, All)
{
    typedef typename TConfig::THits                     THits;
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Size<THit>::Type                   THitSize;
    typedef typename TConfig::TSeeds                    TSeeds;
    typedef typename Id<TSeeds>::Type                   TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename TConfig::TReadSeqs                 TReadSeqs;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    TReadId readSeqId = position(it);

    // Count the hits per read.
    TSeedIds readSeedIds = getSeedIds(me.seeds, readSeqId);
    THitIds readHitIds = getHitIds(me.hits, readSeedIds);
    THitSize readHits = countHits<THitSize>(me.hits, readHitIds);

    // Re-seed hard reads.
    if (readHits > me.options.hitsThreshold)
    {
        // Guess a good seeding stragegy.
        setSeedErrors(me.ctx, readSeqId, (readHits < 200 * me.options.hitsThreshold) ? 1 : 2);

        // Clear the hits of the read.
        clearHits(me.hits, readHitIds);
    }
}

// ----------------------------------------------------------------------------
// Function _classifyReadImpl(); Strata
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeqsIterator>
inline void _classifyReadImpl(ReadsClassifier<TSpec, TConfig> & me, TReadSeqsIterator const & it, Strata)
{
    typedef typename TConfig::THits                     THits;
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Size<THit>::Type                   THitSize;
    typedef typename TConfig::TSeeds                    TSeeds;
    typedef typename Id<TSeeds>::Type                   TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename TConfig::TReadSeqs                 TReadSeqs;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    // Get readSeqId.
    TReadId fwdSeqId = position(it);

    // Get mate id.
    TReadId revSeqId = getFirstMateRevSeqId(me.readSeqs, fwdSeqId);

    // Get seed ids.
    TSeedIds fwdSeedIds = getSeedIds(me.seeds, fwdSeqId);
    TSeedIds revSeedIds = getSeedIds(me.seeds, revSeqId);

    // Get hit ids.
    THitIds fwdHitIds = getHitIds(me.hits, fwdSeedIds);
    THitIds revHitIds = getHitIds(me.hits, revSeedIds);

    // Count the hits of each read.
    THitSize fwdHits = countHits<THitSize>(me.hits, fwdHitIds);
    THitSize revHits = countHits<THitSize>(me.hits, revHitIds);
    THitSize readHits = fwdHits + revHits;

    // Re-seed hard reads.
    if (readHits > me.options.hitsThreshold)
    {
        // Guess a good seeding stragegy.
        unsigned seedErrors = (readHits < 2 * 200 * me.options.hitsThreshold) ? 1 : 2;
        setSeedErrors(me.ctx, fwdSeqId, seedErrors);
        setSeedErrors(me.ctx, revSeqId, seedErrors);

        // Clear the hits of the read.
        clearHits(me.hits, fwdHitIds);
        clearHits(me.hits, revHitIds);
    }
}

// ----------------------------------------------------------------------------
// Function _classifyReadImpl(); AnchorOne
// ----------------------------------------------------------------------------
// Selects the mate to anchor; raises the seeds errors, mark for reseeding and clears the hits of hard anchors.

template <typename TSpec, typename TConfig, typename TReadSeqsIterator>
inline void _classifyReadImpl(ReadsClassifier<TSpec, TConfig> & me, TReadSeqsIterator const & it, All, AnchorOne)
{
    typedef typename TConfig::THits                     THits;
    typedef typename Value<THits>::Type                 THit;
    typedef typename Id<THit>::Type                     THitId;
    typedef Pair<THitId>                                THitIds;
    typedef typename Size<THit>::Type                   THitSize;
    typedef typename TConfig::TSeeds                    TSeeds;
    typedef typename Id<TSeeds>::Type                   TSeedId;
    typedef Pair<TSeedId>                               TSeedIds;
    typedef typename TConfig::TReadSeqs                 TReadSeqs;
    typedef typename Size<TReadSeqs>::Type              TReadId;

    // Get readSeqId.
    TReadId readSeqId = position(it);

    // Get mate id.
    TReadId mateSeqId = getMateSeqId(me.readSeqs, readSeqId);

    // Get seed ids.
    TSeedIds readSeedIds = getSeedIds(me.seeds, readSeqId);
    TSeedIds mateSeedIds = getSeedIds(me.seeds, mateSeqId);

    // Get hit ids.
    THitIds readHitIds = getHitIds(me.hits, readSeedIds);
    THitIds mateHitIds = getHitIds(me.hits, mateSeedIds);

    // Count the hits of each read.
    THitSize readHits = countHits<THitSize>(me.hits, readHitIds);
    THitSize mateHits = countHits<THitSize>(me.hits, mateHitIds);

    TReadId anchorSeqId;
    TReadId otherSeqId;
    THitIds anchorHitIds;
    THitIds otherHitIds;

    // Choose the easiest read as the anchor.
    THitSize anchorHits = std::min(readHits, mateHits);

    if (anchorHits == readHits)
    {
        anchorSeqId = readSeqId;
        anchorHitIds = readHitIds;
        otherSeqId = mateSeqId;
        otherHitIds = mateHitIds;
    }
    else
    {
        anchorSeqId = mateSeqId;
        anchorHitIds = mateHitIds;
        otherSeqId = readSeqId;
        otherHitIds = readHitIds;
    }

    // Clear the hits of the other read.
    clearHits(me.hits, otherHitIds);

    // Re-seed hard anchors.
    if (anchorHits > me.options.hitsThreshold)
    {
        // Guess a good seeding stragegy.
        setSeedErrors(me.ctx, anchorSeqId, (anchorHits < 200 * me.options.hitsThreshold) ? 1 : 2);

        // Clear the hits of the anchor.
        clearHits(me.hits, anchorHitIds);
    }
}

#endif  // #ifndef APP_YARA_MAPPER_CLASSIFIER_H_
