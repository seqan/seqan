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
// Doubly linked list (bucket) of doubly linked lists.  Used for fast
// implementations of certain hot list specializations.
// ==========================================================================

// TODO(holtrew): The interface is semi-good. Enough as an abstraction for internal use, but needs polishing before users can be exposed to it.

#ifndef SEQAN_SYNOPSIS_COUNTER_BUCKETS_H_
#define SEQAN_SYNOPSIS_COUNTER_BUCKETS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Entries_;
typedef Tag<Entries_> Entries;

struct Buckets_;
typedef Tag<Buckets_> Buckets;

template <typename TValue, typename TCount>
struct Bucket_;

template <typename TValue, typename TCount>
struct BucketEntry_
{
    typedef Bucket_<TValue, TCount> TBucket;
    typedef typename std::list<TBucket>::iterator TIterator;

    TValue value;
    TIterator bucketIt;

    BucketEntry_(TValue v, TIterator it) : value(v), bucketIt(it) {}
};

template <typename TValue, typename TCount>
struct Bucket_
{
    typedef BucketEntry_<TValue, TCount> TBucketEntry;

    TCount count;
    std::list<TBucketEntry> entries;

    explicit
    Bucket_(unsigned c) : count(c) {}
};

template <typename TValue, typename TCount>
class CounterBuckets
{
public:
    typedef Bucket_<TValue, TCount> TBucket;
    typedef std::list<TBucket> TList;
    typedef typename TList::iterator TBucketIterator;
    
    TList buckets;
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TValue, typename TCount>
struct Iterator<CounterBuckets<TValue, TCount>, Entries>
{
    typedef typename std::list<BucketEntry_<TValue, TCount> >::iterator Type;
};

template <typename TValue, typename TCount>
struct Iterator<CounterBuckets<TValue, TCount>, Buckets>
{
    typedef typename std::list<Bucket_<TValue, TCount> >::iterator Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function value() for entries iterator.
// ----------------------------------------------------------------------------

template <typename TValue, typename TCount>
inline TCount
value(BucketEntry_<TValue, TCount> const & entry)
{
    return entry.bucketIt->count;
}

template <typename TValue, typename TCount>
inline TCount
value(BucketEntry_<TValue, TCount> & entry)
{
    return value(const_cast<BucketEntry_<TValue, TCount> const &>(entry));
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCount>
inline void
clear(CounterBuckets<TValue, TCount> & buckets)
{
    buckets.buckets.clear();
}

// ----------------------------------------------------------------------------
// Function addCounter()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCount, typename TValue2>
inline typename Iterator<CounterBuckets<TValue, TCount>, Entries>::Type
addCounter(CounterBuckets<TValue, TCount> & buckets, TValue2 const & v)
{
    typedef Bucket_<TValue, TCount> TBucket;
    typedef BucketEntry_<TValue, TCount> TBucketEntry;
    
    // Make sure the list head exists (list not empty) and its count is 1.
    // Insert it if necessary.
    if (buckets.buckets.empty() || buckets.buckets.front().count != static_cast<TCount>(1))
        buckets.buckets.push_front(TBucket(1));

    // Then, add new counter and return it.
    buckets.buckets.front().entries.push_front(TBucketEntry(v, buckets.buckets.begin()));
    return buckets.buckets.front().entries.begin();
}

// ----------------------------------------------------------------------------
// Function eraseCounter()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCount, typename TIterator>
inline void
eraseCounter(CounterBuckets<TValue, TCount> & buckets, TIterator it)
{
    typedef typename Iterator<CounterBuckets<TValue, TCount>, Buckets>::Type TIterator2;
    TIterator2 bucketIt = it->bucketIt;
    bucketIt->erase(it);
    if (bucketIt->empty())
        buckets.buckets.erase(bucketIt);
}

// ----------------------------------------------------------------------------
// Function increaseCounter()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCount, typename TIterator>
inline void
increaseCounter(CounterBuckets<TValue, TCount> & buckets,
                TIterator it)
{
    typedef typename Iterator<CounterBuckets<TValue, TCount>, Entries>::Type TEntryIterator;
    typedef typename Iterator<CounterBuckets<TValue, TCount>, Buckets>::Type TBucketIterator;
    typedef Bucket_<TValue, TCount> TBucket;

    // Set itBNext to point to the next bucket.  If this does not exist or the
    // counter is not current + 1 then insert new one or, if the list would be
    // empty after moving an item out of it, increase counter of current
    // bucket and reuse it.
    TBucketIterator itB = it->bucketIt;
    TBucketIterator itBNext = itB;
    ++itBNext;
    if (itBNext == buckets.buckets.end() || itBNext->count != itB->count + 1) {
        if (itB->entries.size() != 1u) {
            itBNext = buckets.buckets.insert(itBNext, TBucket(itB->count + 1));
        } else {
            itB->count += 1;
            return;
        }
    }

    // Now, move list entry to next bucket and possibly remove old bucket.
    itBNext->entries.splice(itBNext->entries.begin(), itB->entries, it);
    it->bucketIt = itBNext;
    if (itB->entries.empty())
        buckets.buckets.erase(itB);
}

// ----------------------------------------------------------------------------
// Function decreaseAllCounters()
// ----------------------------------------------------------------------------

template <typename TValue, typename TCount>
inline void
decreaseAllCounters(CounterBuckets<TValue, TCount> & buckets)
{
    typedef typename Iterator<CounterBuckets<TValue, TCount>, Buckets>::Type TIterator;
    for (TIterator it = buckets.buckets.begin(), itEnd = buckets.buckets.end(); it != itEnd; ++it)
        it->count -= 1;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SYNOPSIS_COUNTER_BUCKETS_H_
