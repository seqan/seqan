// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Thread-safe suspenable queue
// ==========================================================================
// This queue suspends the caller if it pops a value when the queue was empty
// or appends a value to a full fixed-size queue.

#ifndef SEQAN_PARALLEL_PARALLEL_QUEUE_SUSPENDABLE_H_
#define SEQAN_PARALLEL_PARALLEL_QUEUE_SUSPENDABLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConcurrentQueue
// ----------------------------------------------------------------------------
/*!
 * @class ConcurrentSuspendableQueue Concurrent Suspendable Queue
 * @extends ConcurrentQueue
 * @headerfile <seqan/parallel.h>
 * @brief Thread-safe suspendable queue for multiple producers and multiple consumers.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class ConcurrentQueue<TValue, Suspendable<TSpec> >;
 *
 * @tparam TValue Element type of the queue.
 * @tparam TSpec  Tag for further specializing the Concurrent Queue. Default is <tt>void</tt>.
 *
 * In contrast to the standard @Class.ConcurrentQueue@ this queue suspends the caller
 * if it pops a value when the queue was empty or appends a value to a full
 * fixed-size queue.
 *
 * The implementation uses Mutexes and Events to optionally suspend the calling
 * thread and uses a @Class.AllocString@ as ring buffer.
 *
 */

template <typename TSpec = void>
struct Suspendable;

template <typename TValue, typename TSpec>
class ConcurrentQueue<TValue, Suspendable<TSpec> >
{
public:
    typedef typename Host<ConcurrentQueue>::Type    TString;
    typedef typename Size<TString>::Type            TSize;

    size_t                  readerCount;
    size_t                  writerCount;

    TString                 data;
    TSize                   occupied;
    TSize                   back;
    TSize                   front;

    std::mutex              cs;
    std::condition_variable more;

    bool                    virgin;

    ConcurrentQueue():
        readerCount(0),
        writerCount(0),
        occupied(0),
        back(0),
        front(0),
        virgin(true)
    {}

    ~ConcurrentQueue()
    {
        SEQAN_ASSERT_EQ(writerCount, 0u);

        // wait for all pending readers to finish
        while (readerCount != 0u)
        {}

        typename Iterator<TString, Standard>::Type arrayBegin = begin(data, Standard());

        if (occupied != 0)
        {
            if (front < back)
            {
                arrayDestruct(arrayBegin + front, arrayBegin + back);
            }
            else
            {
                arrayDestruct(arrayBegin, arrayBegin + back);
                arrayDestruct(arrayBegin + front, arrayBegin + capacity(data));
            }
        }
        _setLength(data, 0);
    }
};

template <typename TValue>
class ConcurrentQueue<TValue, Suspendable<Limit> >:
    public ConcurrentQueue<TValue, Suspendable<> >
{
public:
    typedef ConcurrentQueue<TValue, Suspendable<> > TBase;
    typedef typename Host<ConcurrentQueue>::Type    TString;
    typedef typename Size<TString>::Type            TSize;

    std::condition_variable less;

    ConcurrentQueue(TSize maxSize):
        TBase()
    {
        reserve(this->data, maxSize, Exact());
        _setLength(this->data, maxSize);
    }

    ConcurrentQueue(ConcurrentQueue const & other):
        TBase((TBase const &)other)
    {}
};

template <typename TValue>
struct DefaultOverflowImplicit<ConcurrentQueue<TValue, Suspendable<Limit> > >
{
    typedef Limit Type;
};

template <typename TValue, typename TSpec>
inline void
lockReading(ConcurrentQueue<TValue, Suspendable<TSpec> > &)
{}

template <typename TValue, typename TSpec>
inline void
unlockReading(ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    std::unique_lock<std::mutex> lock(me.cs);
    if (--me.readerCount == 0u)
    {
        lock.unlock();
        me.less.notify_all();
    }
}

template <typename TValue, typename TSpec>
inline void
lockWriting(ConcurrentQueue<TValue, Suspendable<TSpec> > &)
{}

template <typename TValue, typename TSpec>
inline void
unlockWriting(ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    std::unique_lock<std::mutex> lk(me.cs);
    if (--me.writerCount == 0u)
    {
        lk.unlock();
        me.more.notify_all();
    }
}

template <typename TValue, typename TSize, typename TSpec>
inline void
setReaderCount(ConcurrentQueue<TValue, Suspendable<TSpec> > & me, TSize readerCount)
{
    std::unique_lock<std::mutex> lock(me.cs);
    me.readerCount = readerCount;
}

template <typename TValue, typename TSize, typename TSpec>
inline void
setWriterCount(ConcurrentQueue<TValue, Suspendable<TSpec> > & me, TSize writerCount)
{
    std::unique_lock<std::mutex> lock(me.cs);
    me.writerCount = writerCount;
}

template <typename TValue, typename TSize1, typename TSize2, typename TSpec>
inline void
setReaderWriterCount(ConcurrentQueue<TValue, Suspendable<TSpec> > & me, TSize1 readerCount, TSize2 writerCount)
{
    std::unique_lock<std::mutex> lock(me.cs);
    me.readerCount = readerCount;
    me.writerCount = writerCount;
}

template <typename TValue, typename TSize, typename TSpec>
inline bool
waitForMinSize(ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
               TSize minSize)
{
    std::unique_lock<std::mutex> lock(me.cs);
    while (me.occupied < minSize && me.writerCount > 0u)
        me.more.wait(lock);
    return me.occupied >= minSize;
}

template <typename TValue, typename TSpec>
inline bool
empty(ConcurrentQueue<TValue, Suspendable<TSpec> > const & me)
{
    return me.occupied == 0;
}

template <typename TValue, typename TSpec>
inline typename Size<ConcurrentQueue<TValue, Suspendable<TSpec> > >::Type
length(ConcurrentQueue<TValue, Suspendable<TSpec> > const & me)
{
    return me.occupied;
}


template <typename TValue, typename TSpec>
inline bool
_popFront(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
          std::unique_lock<std::mutex> & lk)
{
    typedef ConcurrentQueue<TValue, Suspendable<TSpec> >    TQueue;
    typedef typename Host<TQueue>::Type                     TString;
    typedef typename Size<TString>::Type                    TSize;
    typedef typename Iterator<TString, Standard>::Type      TIter;

    TSize cap = capacity(me.data);

    while (me.occupied == 0u && me.writerCount > 0u)
        me.more.wait(lk);

    if (me.occupied == 0u)
        return false;

    SEQAN_ASSERT_NEQ(me.occupied, 0u);

    // extract value and destruct it in the data string
    TIter it = begin(me.data, Standard()) + me.front;
    std::swap(result, *it);
    valueDestruct(it);

    me.front = (me.front + 1) % cap;
    me.occupied--;

    /* now: either me.occupied > 0 and me.nextout is the index
       of the next occupied slot in the buffer, or
       me.occupied == 0 and me.nextout is the index of the next
       (empty) slot that will be filled by a producer (such as
       me.nextout == me.nextin) */

    return true;
}

template <typename TValue, typename TSpec>
inline bool
_popBack(TValue & result,
         ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
         std::unique_lock<std::mutex> & lk)
{
    typedef ConcurrentQueue<TValue, Suspendable<TSpec> >    TQueue;
    typedef typename Host<TQueue>::Type                     TString;
    typedef typename Size<TString>::Type                    TSize;
    typedef typename Iterator<TString, Standard>::Type      TIter;

    TSize cap = capacity(me.data);

    while (me.occupied == 0u && me.writerCount > 0u)
        me.more.wait(lk);

    if (me.occupied == 0u)
        return false;

    SEQAN_ASSERT_NEQ(me.occupied, 0u);

    me.back = (me.back + cap - 1) % cap;

    // extract value and destruct it in the data string
    TIter it = begin(me.data, Standard()) + me.back;
    std::swap(result, *it);
    valueDestruct(it);

    me.occupied--;

    /* now: either me.occupied > 0 and me.nextout is the index
       of the next occupied slot in the buffer, or
       me.occupied == 0 and me.nextout is the index of the next
       (empty) slot that will be filled by a producer (such as
       me.nextout == me.nextin) */

    return true;
}

template <typename TValue, typename TSpec>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    std::unique_lock<std::mutex> lock(me.cs);
    return _popFront(result, me, lock);
}

template <typename TValue>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, Suspendable<Limit> > & me)
{
    {
        std::unique_lock<std::mutex> lk(me.cs);
        if (!_popFront(result, me, lk))
            return false;
    }
    me.less.notify_all();
    return true;
}

template <typename TValue, typename TSpec>
inline bool
popBack(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    std::unique_lock<std::mutex> lk(me.cs);
    return _popBack(result, me, lk);
}

template <typename TValue>
inline bool
popBack(TValue & result, ConcurrentQueue<TValue, Suspendable<Limit> > & me)
{
    {
        std::unique_lock<std::mutex> lk(me.cs);
        if (!_popBack(result, me, lk))
            return false;
    }
    me.less.notify_all();
    return true;
}


template <typename TValue, typename TValue2, typename TSpec, typename TExpand>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
            TValue2 && val,
            Tag<TExpand> expandTag)
{
    typedef ConcurrentQueue<TValue, Suspendable<TSpec> >    TQueue;
    typedef typename Host<TQueue>::Type                     TString;
    typedef typename Size<TString>::Type                    TSize;

    {
        std::lock_guard<std::mutex> lock(me.cs);
        TSize cap = capacity(me.data);

        if (me.occupied >= cap)
        {
            // increase capacity
            _setLength(me.data, cap);
            reserve(me.data, cap + 1, expandTag);
            TSize delta = capacity(me.data) - cap;

            // create a gap of delta many values between tail and head
            _clearSpace(me.data, delta, me.back, me.back, expandTag);
            if (me.occupied != 0 && me.back <= me.front)
                me.front += delta;

            cap += delta;
        }

        valueConstruct(begin(me.data, Standard()) + me.back, val);
        me.back = (me.back + 1) % cap;

        me.occupied++;
    }

    /* now: either me.occupied < BSIZE and me.nextin is the index
     of the next empty slot in the buffer, or
     me.occupied == BSIZE and me.nextin is the index of the
     next (occupied) slot that will be emptied by a consumer
     (such as me.nextin == me.nextout) */
    me.more.notify_all();
    return true;
}

template <typename TValue, typename TValue2, typename TSpec, typename TExpand>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<Limit> > & me,
            TValue2 && val,
            Tag<TExpand> expandTag);

template <typename TValue, typename TValue2>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<Limit> > & me,
            TValue2 && val,
            Limit)
{
    typedef ConcurrentQueue<TValue, Suspendable<Limit> >    TQueue;
    typedef typename Host<TQueue>::Type                     TString;
    typedef typename Size<TString>::Type                    TSize;

    {
        std::unique_lock<std::mutex> lock(me.cs);
        TSize cap = capacity(me.data);

        while (me.occupied >= cap && me.readerCount > 0u)
            me.less.wait(lock);

        if (me.occupied >= cap)
            return false;

        SEQAN_ASSERT_LT(me.occupied, cap);

        valueConstruct(begin(me.data, Standard()) + me.back, val);
        me.back = (me.back + 1) % cap;
        me.occupied++;
    }

    /* now: either me.occupied < BSIZE and me.nextin is the index
     of the next empty slot in the buffer, or
     me.occupied == BSIZE and me.nextin is the index of the
     next (occupied) slot that will be emptied by a consumer
     (such as me.nextin == me.nextout) */
    me.more.notify_all();
    return true;
}

template <typename TValue, typename TValue2, typename TSpec>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
            TValue2 && val)
{
    return appendValue(me, val, typename DefaultOverflowImplicit<ConcurrentQueue<TValue, Suspendable<TSpec> > >::Type());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_QUEUE_SUSPENDABLE_H_
