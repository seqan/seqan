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
// Implementation of the Frequent hot list algorithm.
// ==========================================================================

// TODO(holtgrew): We should include the slist+list trick by Karp et al.

#ifndef SEQAN_SYNOPSIS_HOT_LIST_FREQUENT_H_
#define SEQAN_SYNOPSIS_HOT_LIST_FREQUENT_H_

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
.Spec.Frequent HotList
..cat:Synopsis Data Structures
..general:Class.HotList
..summary:Approximate hot list for top-k queries using the Frequent algorithm.
..signature:HotList<TValue, Frequent>
..param.TValue:Type of the items to keep in hot list.
..remarks:Update time is $O(1)$, space usage is $O(k)$. Performance guarantee is $(-n/(k+1), 0)$
..remarks:See E. Demaine, A. Lopez-Ortiz, and J. I. Munro. Frequency estimation of internet packet streams with limited space. In European Symposium on Algorithms (ESA), 2002.
..remarks:See R. Karp, C. Papadimitriou, and S. Shenker. A simple algorithm for ﬁnding frequent elements in sets and bags.  ACM Transactions on Database Systems, 28:51–55, 2003. 
..include:seqan/synopsis.h

.Memfunc.Frequent HotList#HotList
..class:Spec.Frequent HotList
..summary:Constructor
..signature:HotList(k[, callback])
..signature:HotList(k)
..param.k:The number of most frequent items to keep.
..param.callback:Optional callback context object adhering to the HotList Callback Concept.
 */

struct Frequent_;
typedef Tag<Frequent_> Frequent;

template <typename TValue, typename TCallback>
class HotList<TValue, Frequent, TCallback>
{
public:
    typedef HotList<TValue, Frequent> THotList_;
    typedef typename Size<THotList_>::Type TSize_;
    typedef typename std::map<TValue, TSize_> TMap_;
    
    TSize_ _n;   // number of items seen so far, unnecessary?
    unsigned _k; // number of buckets
    TMap_ _map;  // dictionary

    Holder<TCallback> _callbackContext;
    
    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    HotList(unsigned k)
            : _n(0), _k(k)
    {
        create(_callbackContext);
    }

    HotList(unsigned k, TCallback & callbackContext)
            : _n(0), _k(k), _callbackContext(callbackContext)
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
clear(HotList<TValue, Frequent, TCallback> & hotList)
{
    hotList._n = 0;
    hotList._map.clear();
}

// ----------------------------------------------------------------------------
// Function registerItem()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCallback>
inline typename Size<HotList<TValue, Frequent> >::Type
registerItem(HotList<TValue, Frequent, TCallback> & hotList, TValue const & x)
{
    typedef HotList<TValue, Frequent, TCallback> THotList;
    typedef typename THotList::TMap_ TMap;
    typedef typename THotList::TMap_::iterator TIterator;
    typedef typename Size<THotList>::Type TSize;

    hotList._n += 1;

    TIterator it = hotList._map.find(x);
    if (it != hotList._map.end()) {
        // Case 1: x can be found in hot list and we can simply increase
        // its counter.
        it->second += 1;
        return it->second;
    }

    if (hotList._map.size() < hotList._k - 1) {
        // Case 2: x is not yet in hot list but there still is room.
        hotList._map.insert(typename TMap::value_type(x, 1));
        hotListItemAdded(_callbackContext(hotList), x, static_cast<TSize>(1), MaxValue<TSize>::VALUE);
        return 1;
    }

    // Case 3: value is not yet in hot list and we use it to trigger a round
    // of decreasing.
    for (TIterator it = hotList._map.begin(), itEnd = hotList._map.end(); it != itEnd;) {
        if (it->second == 1) {
            // Erasing only invalidates iterators that are equal to it.
            TIterator it2 = it;
            it++;
            hotList._map.erase(it2);
            hotListItemRemoved(_callbackContext(hotList), x, static_cast<TSize>(1), MaxValue<unsigned>::VALUE);
        } else {
            it->second -= 1;
            it++;
        }
    }

    // TODO(holtgrew): Can we insert x now?
    return 0;
}

// ----------------------------------------------------------------------------
// Function removeItem()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCallback>
inline void
removeItem(HotList<TValue, Frequent, TCallback> & hotList, TValue const & x)
{
    typedef HotList<TValue, Frequent, TCallback> THotList;
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
         HotList<TValue, Frequent, TCallback> const & hotList)
{
    typedef HotList<TValue, Frequent, TCallback> THotList;
    typedef typename THotList::TMap_ TMap;
    typedef typename THotList::TMap_::const_iterator TIterator;
    typedef typename Size<HotList<TValue, Frequent> >::Type TSize;
    typedef Triple<TValue, TSize, TSize> TTriple;

    clear(items);
    reserve(items, hotList._map.size(), Exact());

    for (TIterator it = hotList._map.begin(), itEnd = hotList._map.end(); it != itEnd; ++it)
        appendValue(items, TTriple(it->first, 0, it->second));
    typedef CmpTripleI2DescI3DescI1Asc_<TValue, TSize> TCmp;
    std::sort(begin(items, Standard()), end(items, Standard()), TCmp());
}

}  // namespace seqan

#endif  // SEQAN_SYNOPSIS_HOT_LIST_FREQUENT_H_
