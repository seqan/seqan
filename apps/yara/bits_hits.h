// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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

#ifndef APP_YARA_BITS_HITS_H_
#define APP_YARA_BITS_HITS_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Hit<Exact>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec = Exact>
struct Hit
{
    typename Position<Hit>::Type    range;
};

// ----------------------------------------------------------------------------
// Class Hit<HammingDistance>
// ----------------------------------------------------------------------------

template <typename TSize>
struct Hit<TSize, HammingDistance>
{
    typename Position<Hit>::Type    range;
    typename Id<Hit>::Type          seedId;
    unsigned char                   errors;
};

namespace seqan
{

// ----------------------------------------------------------------------------
// Metafunction Size<Hit>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
struct Size<Hit<TSize, TSpec> >
{
    typedef TSize Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position<Hit>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
struct Position<Hit<TSize, TSpec> >
{
    typedef Pair<TSize> Type;
};

// ----------------------------------------------------------------------------
// Metafunction Id<Hit>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
struct Id<Hit<TSize, TSpec> >
{
    typedef unsigned int Type;
};

// ----------------------------------------------------------------------------
// Metafunction Spec<Hit>
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
struct Spec<Hit<TSize, TSpec> >
{
    typedef TSpec Type;
};
}

// ----------------------------------------------------------------------------
// Class HitsCounter
// ----------------------------------------------------------------------------

template <typename TSize, typename TThreading, typename TSpec = void>
struct HitsCounter
{
    TSize count;

    HitsCounter() :
        count(0)
    {}

    template <typename THit>
    void operator() (THit const & hit)
    {
        atomicAdd(count, getCount(hit), TThreading());
    }
};

// ----------------------------------------------------------------------------
// Class HitsSorterByXXX
// ----------------------------------------------------------------------------

template <typename THit>
struct HitsSorterByCount
{
    inline bool operator()(THit const & a, THit const & b) const
    {
        return getCount(a) < getCount(b);
    }
};

// ============================================================================
// Functions
// ============================================================================

template <typename TSize, typename TSpec>
inline bool operator< (Hit<TSize, TSpec> const & a, Hit<TSize, TSpec> const & b)
{
    return a.seedId < b.seedId;
}

// ----------------------------------------------------------------------------
// Function clearRange()
// ----------------------------------------------------------------------------

template <typename THit>
inline THit clearRange(THit hit)
{
    setValueI1(hit.range, 0);
    setValueI2(hit.range, 0);
    return hit;
}

// ----------------------------------------------------------------------------
// Function getCount()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
inline TSize
getCount(Hit<TSize, TSpec> const & hit)
{
    return getValueI2(hit.range) - getValueI1(hit.range);
}

// ----------------------------------------------------------------------------
// Function getHitErrors()
// ----------------------------------------------------------------------------

template <typename TSize, typename TSpec>
inline unsigned char
getErrors(Hit<TSize, TSpec> const & /* hit */)
{
    return 0;
}

template <typename TSize>
inline unsigned char
getErrors(Hit<TSize, HammingDistance> const & hit)
{
    return hit.errors;
}

// ----------------------------------------------------------------------------
// Function getErrors()
// ----------------------------------------------------------------------------

template <typename THits, typename THitId>
inline unsigned char
getErrors(THits const & hits, THitId hitId)
{
    return getErrors(hits[hitId]);
}

// ----------------------------------------------------------------------------
// Function getRange()
// ----------------------------------------------------------------------------

template <typename THits, typename THitId>
inline typename Position<typename Value<THits>::Type>::Type
getRange(THits const & hits, THitId hitId)
{
    return hits[hitId].range;
}

// ----------------------------------------------------------------------------
// Function getSeedId()
// ----------------------------------------------------------------------------

template <typename THits, typename THitId>
inline typename Id<typename Value<THits>::Type>::Type
getSeedId(THits const & hits, THitId hitId)
{
    typedef typename Value<THits>::Type THit;
    typedef typename Spec<THit>::Type   THitSpec;

    return _getSeedId(hits, hitId, THitSpec());
}

template <typename THits, typename THitId>
inline typename Id<typename Value<THits>::Type>::Type
_getSeedId(THits const & /* hits */, THitId hitId, Exact)
{
    return hitId;
}

template <typename THits, typename THitId>
inline typename Id<typename Value<THits>::Type>::Type
_getSeedId(THits const & hits, THitId hitId, HammingDistance)
{
    return hits[hitId].seedId;
}

// ----------------------------------------------------------------------------
// Function getHitIds()
// ----------------------------------------------------------------------------

template <typename THits, typename TSeedId>
inline Pair<typename Id<typename Value<THits>::Type>::Type>
getHitIds(THits const & hits, TSeedId seedId)
{
    return getHitIds(hits, Pair<TSeedId>(seedId, seedId + 1));
}

template <typename THits, typename TSeedId>
inline Pair<typename Id<typename Value<THits>::Type>::Type>
getHitIds(THits const & hits, Pair<TSeedId> seedIds)
{
    typedef typename Value<THits>::Type THit;
    typedef typename Spec<THit>::Type   THitSpec;

    return _getHitIds(hits, seedIds, THitSpec());
}

template <typename THits, typename TSeedId>
inline Pair<typename Id<typename Value<THits>::Type>::Type>
_getHitIds(THits const & /* hits */, Pair<TSeedId> seedIds, Exact)
{
    return seedIds;
}

template <typename THits, typename TSeedId>
inline Pair<typename Id<typename Value<THits>::Type>::Type>
_getHitIds(THits const & hits, Pair<TSeedId> seedIds, HammingDistance)
{
    typedef typename Value<THits>::Type                     THit;
    typedef typename Id<THit>::Type                         THitId;
    typedef Pair<THitId>                                    THitIds;
    typedef typename Iterator<THits const, Standard>::Type  THitsIterator;

    THit firstSeedHit;
    THit lastSeedHit;
    firstSeedHit.seedId = getValueI1(seedIds);
    lastSeedHit.seedId = getValueI2(seedIds) - 1;

    THitsIterator hitsBegin = begin(hits, Standard());
    THitsIterator hitsEnd = end(hits, Standard());

    THitsIterator firstHit = std::lower_bound(hitsBegin, hitsEnd, firstSeedHit);
    THitsIterator lastHit = std::upper_bound(hitsBegin, hitsEnd, lastSeedHit);

    return THitIds(position(firstHit, hits), position(lastHit, hits));
}

// ----------------------------------------------------------------------------
// Function sortHits()
// ----------------------------------------------------------------------------

template <typename THits, typename TThreading>
inline void sortHits(THits & hits, TThreading const & threading)
{
    typedef typename Value<THits>::Type THit;
    typedef typename Spec<THit>::Type   THitSpec;

    _sortHits(hits, THitSpec(), threading);
}

template <typename THits, typename TThreading>
inline void _sortHits(THits & /* hits */, Exact, TThreading const & /* threading */) {}

template <typename THits, typename TThreading>
inline void _sortHits(THits & hits, HammingDistance, TThreading const & threading)
{
    return sort(hits, threading);
}

// ----------------------------------------------------------------------------
// Function countHits()
// ----------------------------------------------------------------------------

template <typename TSize, typename THits, typename TThreading>
inline TSize countHits(THits const & hits, TThreading const & threading)
{
    return forEach(hits, HitsCounter<TSize, TThreading>(), threading).count;
}

// ----------------------------------------------------------------------------
// Function countHits()
// ----------------------------------------------------------------------------

template <typename TSize, typename THits, typename THitId>
inline TSize countHits(THits const & hits, Pair<THitId> hitIds)
{
    return forEach(infix(hits, getValueI1(hitIds), getValueI2(hitIds)),
                   HitsCounter<TSize, Serial>(), Serial()).count;
}

// ----------------------------------------------------------------------------
// Function clearHits()
// ----------------------------------------------------------------------------

template <typename THits, typename THitId>
inline void clearHits(THits & hits, Pair<THitId> hitIds)
{
    typedef typename Value<THits>::Type THit;

    std::transform(begin(hits, Standard()) + getValueI1(hitIds),
                   begin(hits, Standard()) + getValueI2(hitIds),
                   begin(hits, Standard()) + getValueI1(hitIds),
                   clearRange<THit>);
}

#endif  // #ifndef APP_YARA_BITS_HITS_H_
