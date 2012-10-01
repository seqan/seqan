// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Implementation of the Space Saving (list variant) hot list algorithm.
// ==========================================================================

#ifndef SEQAN_SYNOPSIS_HOT_LIST_SPACE_SAVING_H_
#define SEQAN_SYNOPSIS_HOT_LIST_SPACE_SAVING_H_

#include <algorithm>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Data structures for the lists maintaining the counts.

/**
.Spec.Space Saving HotList
..summary:Approximate hot list using the Space Saving algorithm (list variant).
..cat:Synopsis Data Structures
..general:Class.HotList
..signature:HotList<TValue, SpaceSaving>
..param.TValue:Type of the items to keep in hot list.
..remarks:See Ahmed Metwally,  Divyakant Agrawal, and  Amr El Abbadi. Efficient Computation of SpaceSaving and Top-k Elements in Data Streams. ICDT '05, pp. 555-566.
..include:seqan/synopsis.h

.Memfunc.Space Saving HotList#HotList
..class:Spec.Space Saving HotList
..summary:Constructor
..signature:HotList(k)
..param.k:The number of most frequent items to keep.
 */

struct SpaceSaving_;
typedef Tag<SpaceSaving_> SpaceSaving;

template <typename TValue, typename TCallback>
class HotList<TValue, SpaceSaving, TCallback>
{
public:
    typedef HotList<TValue, SpaceSaving> THotList_;
    typedef typename Size<THotList_>::Type TSize_;
    typedef CounterBuckets<TValue, TSize_> TCounterBuckets_;
    typedef typename Iterator<TCounterBuckets_, Entries>::Type TEntriesIter_;
    typedef std::tr1::unordered_map<TValue, TEntriesIter_> TMap_;

    TSize_ _n;
    unsigned _k;
    TSize_ _size;  // Number of elements in hot list.
    TCounterBuckets_ _buckets;
    TMap_ _map;

    Holder<TCallback> _callbackContext;
    
    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    HotList(unsigned k)
            : _n(0), _k(k), _size(0)
    {
        create(_callbackContext);
    }

    HotList(unsigned k, TCallback & callbackContext)
            : _n(0), _k(k), _size(0), _callbackContext(callbackContext)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCallback>
inline void
clear(HotList<TValue, SpaceSaving, TCallback> & hotList)
{
    hotList._n = 0;
    hotList._size = 0;
    clear(hotList._buckets);
    hotList._map.clear();
}

// ----------------------------------------------------------------------------
// Function registerItem()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCallback>
inline typename Size<HotList<TValue, SpaceSaving> >::Type
registerItem(HotList<TValue, SpaceSaving, TCallback> & hotList,
             TValue const & x)
{
    typedef HotList<TValue, SpaceSaving, TCallback> THotList;
    typedef typename THotList::TMap_ TMap;
    typedef typename THotList::TMap_::iterator TMapIterator;
    typedef typename THotList::TCounterBuckets_ TCounterBuckets;
    typedef typename Size<THotList>::Type TSize;
    typedef typename Iterator<TCounterBuckets, Entries>::Type TEntriesIterator;

    hotList._n += 1;

    // Case 1: Entry is in hot list.  We can simply increase the count.
    TMapIterator it = hotList._map.find(x);
    if (it != hotList._map.end()) {
        TSize prev = value(*it->second);
        (void) prev;
        increaseCounter(hotList._buckets, it->second);
        SEQAN_ASSERT_EQ(prev + 1, value(*it->second));
        return value(*it->second);
    }

    // Case 2: Entry is not yet in hot list and there still is space.
    if (hotList._size < hotList._k) {
        hotList._size += 1;
        TEntriesIterator itE = addCounter(hotList._buckets, x);
        hotList._map.insert(typename TMap::value_type(x, itE));
        return 1;
    }

    // Case 3: Entry is not yet in hot list and there is not enough space.  In
    // this case, we hijack one of the existing entries.
    //
    // Erase old.
    SEQAN_ASSERT_NOT(hotList._buckets.buckets.empty());
    TEntriesIterator itE = hotList._buckets.buckets.front().entries.begin();
    TMapIterator itM = hotList._map.find(itE->value);
    SEQAN_ASSERT(itM != hotList._map.end());
    hotList._map.erase(itM);
    // Insert new.
    itE->value = x;
    hotList._map.insert(typename TMap::value_type(x, itE));
    increaseCounter(hotList._buckets, itE);
    return value(*itE);
}

// ----------------------------------------------------------------------------
// Function removeItem()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCallback>
inline void
removeItem(HotList<TValue, SpaceSaving, TCallback> & hotList,
           TValue const & x)
{
    typedef HotList<TValue, SpaceSaving, TCallback> THotList;
    typedef typename THotList::TMap_::iterator TIterator;
    TIterator it = hotList._map.find(x);
    SEQAN_ASSERT_NOT(it == hotList.end());
    eraseCounter(hotList._buckets, it->second);
    hotList._map.erase(it);
}

// ----------------------------------------------------------------------------
// Function getItems()
// ----------------------------------------------------------------------------

template <typename TResult, typename TValue, typename TCallback>
inline void
getItems(TResult & items,
         HotList<TValue, SpaceSaving, TCallback> const & hotList)
{
    typedef HotList<TValue, Frequent, TCallback> THotList;
    typedef typename Size<THotList>::Type TSize;
    typedef typename std::list<Bucket_<TValue, TSize> >::const_reverse_iterator TBucketRIterator;
    typedef typename std::list<BucketEntry_<TValue, TSize> >::const_iterator TBucketEntryIterator;
    typedef Triple<TValue, TSize, TSize> TTriple;

    // Sum up buckets lengths up to k elements.
    unsigned sum = 0;
    for (TBucketRIterator it = hotList._buckets.buckets.rbegin(), itEnd = hotList._buckets.buckets.rend(); it != itEnd; ++it) {
        sum += it->entries.size();
        if (sum >= hotList._k)
            break;
    }

    // Allocate memory.
    clear(items);
    reserve(items, sum, Exact());

    // Collect results.
    for (TBucketRIterator it = hotList._buckets.buckets.rbegin(), itEnd = hotList._buckets.buckets.rend(); it != itEnd; ++it) {
        for (TBucketEntryIterator it2 = it->entries.begin(), it2End = it->entries.end(); it2 != it2End; ++it2)
            appendValue(items, TTriple(it2->value, it->count, MaxValue<TSize>::VALUE));
    }
}

}  // namespace seqan

#endif  // SEQAN_SYNOPSIS_HOT_LIST_SPACE_SAVING_H_
