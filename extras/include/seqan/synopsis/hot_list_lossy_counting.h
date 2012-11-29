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
// Implementation of the Lossy Counting (with deltas) hot list algorithm.
// ==========================================================================

#ifndef SEQAN_SYNOPSIS_HOT_LIST_LOSSY_COUNTING_H_
#define SEQAN_SYNOPSIS_HOT_LIST_LOSSY_COUNTING_H_

#include <algorithm>
#include <map>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Lossy Counting HotList
..cat:Synopsis Data Structures
..general:Class.HotList
..summary:Approximate hot list for iceberg queries using the Lossy Counting algorithm.
..signature:HotList<TValue, LossyCounting>
..param.TValue:Type of the items to keep in hot list.
..remarks:See Gurmeet Singh Manku and Rajeev Motwani.  Approximate Frequency Counts over Data Streams.  VLDB '02.
..include:seqan/synopsis.h

.Memfunc.Lossy Counting HotList#HotList
..class:Spec.Lossy Counting HotList
..summary:Constructor
..signature:HotList(k)
..param.k:Gives $1/epsilon$, will keep approximately the $epsilon * n$ most frequent items.
 */

struct LossyCounting_;
typedef Tag<LossyCounting_> LossyCounting;

template <typename TValue, typename TCallback>
class HotList<TValue, LossyCounting, TCallback>
{
public:
    typedef HotList<TValue, LossyCounting> THotList_;
    typedef typename Size<THotList_>::Type TSize_;
    typedef typename std::map<TValue, TSize_> TMap_;

    TSize_ _n;
    TSize_ _delta;  // One global delta.
    unsigned _k;
    TMap_ _map;  // dictionary

    Holder<TCallback> _callbackContext;
    
    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    HotList(unsigned k)
            : _n(0), _delta(0), _k(k)
    {
        create(_callbackContext);
    }

    HotList(unsigned k, TCallback & callbackContext)
            : _n(0), _delta(0), _k(k), _callbackContext(callbackContext)
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
clear(HotList<TValue, LossyCounting, TCallback> & hotList)
{
    hotList._n = 0;
    hotList._delta = 0;
    hotList._map.clear();
}

// ----------------------------------------------------------------------------
// Function registerItem()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCallback>
inline typename Size<HotList<TValue, LossyCounting, TCallback> >::Type
registerItem(HotList<TValue, LossyCounting, TCallback> & hotList, TValue const & x)
{
    typedef HotList<TValue, LossyCounting> THotList;
    typedef typename THotList::TMap_ TMap;
    typedef typename THotList::TMap_::iterator TIterator;
    typedef typename Size<THotList>::Type TSize;

    hotList._n += 1;

    // First, update / insert record in map.
    
    TIterator itX = hotList._map.find(x);
    if (itX != hotList._map.end()) {
        // Case 1: x can be found in hot list and we can simply increase
        // its counter.
        itX->second += 1;
    } else {
        // Case 2: Otherwise, register it with a count of 1.
        itX = hotList._map.insert(typename TMap::value_type(x, 1 + hotList._delta)).first;
        hotListItemAdded(_callbackContext(hotList), x, static_cast<TSize>(1), 1 + hotList._delta);
    }

    // Second, compress buckets if at window border.

    if (hotList._n / hotList._k == hotList._delta)
        return itX->second;  // Guard condition to reduce indentation.

    TSize result = itX->second;
    hotList._delta = hotList._n / hotList._k;
    for (TIterator it = hotList._map.begin(), itEnd = hotList._map.end(); it != itEnd;) {
        if (it->second < hotList._delta) {
            if (it == itX)  // Case: inserted item is immediately removed.
                result = 0;
            TIterator it2 = it;
            ++it;
            hotList._map.erase(it2);
        } else {
            ++it;
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
// Function removeItem()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCallback>
inline void
removeItem(HotList<TValue, LossyCounting, TCallback> & hotList, TValue const & x)
{
    // TODO(holtgrew): Duplicate from Frequent specialization.
    typedef HotList<TValue, LossyCounting> THotList;
    typedef typename THotList::TMap_ TMap;
    typedef typename THotList::TMap_::const_iterator TIterator;

    TIterator it = hotList._map.find(x);
    if (it == hotList._map.end())
        return;
    hotList.erase(it);
}

// ----------------------------------------------------------------------------
// Function getItems()
// ----------------------------------------------------------------------------

template <typename TResult, typename TValue, typename TCallback>
inline void
getItems(TResult & items,
         HotList<TValue, LossyCounting, TCallback> const & hotList)
{
    typedef HotList<TValue, Frequent> THotList;
    typedef typename THotList::TMap_ TMap;
    typedef typename THotList::TMap_::const_iterator TIterator;
    typedef typename Size<HotList<TValue, Frequent> >::Type TSize;
    typedef Triple<TValue, TSize, TSize> TTriple;

    clear(items);
    reserve(items, hotList._map.size(), Exact());

    for (TIterator it = hotList._map.begin(), itEnd = hotList._map.end(); it != itEnd; ++it)
        appendValue(items, TTriple(it->first, it->second, it->second + hotList._delta));
    typedef CmpTripleI2DescI3DescI1Asc_<TValue, TSize> TCmp;
    std::sort(begin(items, Standard()), end(items, Standard()), TCmp());
}

}  // namespace seqan

#endif  // SEQAN_SYNOPSIS_HOT_LIST_LOSSY_COUNTING_H_
