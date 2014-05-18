// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Thread-safe queue
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_QUEUE_H_
#define SEQAN_PARALLEL_PARALLEL_QUEUE_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConcurrentQueue
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec = void>
struct ConcurrentQueue
{
    typedef typename Host<ConcurrentQueue>::Type    TString;
    typedef typename Size<TString>::Type            TSize;

    TString                 data;
    mutable ReadWriteLock   lock;

    volatile TSize headPos;
    volatile TSize headReadPos;
    volatile TSize tailPos;
    volatile TSize tailWritePos;
    TSize roundSize;

    volatile unsigned readerCount;
    volatile unsigned writerCount;
    volatile unsigned pushed, popped;
    volatile bool firstValue;


    ConcurrentQueue() :
        headPos(0),
        headReadPos(0),
        tailPos(0),
        tailWritePos(0),
        roundSize(0),
        readerCount(0),
        writerCount(0),
        pushed(0), popped(0),
        firstValue(false)
    {}

    // you can set the initial capacity here
    ConcurrentQueue(TSize initCapacity) :
        headPos(0),
        headReadPos(0),
        tailPos(0),
        tailWritePos(0),
        readerCount(0),
        writerCount(0),
        pushed(0), popped(0),
        firstValue(false)
    {
        reserve(data, initCapacity + 1, Exact());
        roundSize = (TSize)1 << (log2(capacity(data) - 1) + 1);
    }

    ConcurrentQueue(TString & data):
        data(data),
        headPos(0),
        headReadPos(0),
        tailPos(length(data)),
        tailWritePos(length(data)),
        readerCount(0),
        writerCount(0),
        pushed(0), popped(0),
        firstValue(false)
    {
        roundSize = (TSize)1 << (log2(capacity(data) - 1) + 1);
    }

    ~ConcurrentQueue()
    {
        SEQAN_ASSERT_EQ(tailPos, tailWritePos);
        SEQAN_ASSERT_EQ(headPos, headReadPos);
        SEQAN_ASSERT(empty(lock));
        SEQAN_ASSERT_EQ(writerCount, 0u);

        // wait for all pending readers to finish
        while (readerCount != 0)
        {}

        TSize cap = capacity(data);
        if (tailPos < headPos)
        {
            _clearSpace(data, 0u, tailPos, headPos, Insist());
            _setLength(data, cap - (headPos - tailPos));
        }
        else
        {
            _setLength(data, tailPos);
            _clearSpace(data, 0u, (TSize)0, headPos, Insist());
            _setLength(data, tailPos - headPos);
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Value<ConcurrentQueue<TValue, TSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Host<ConcurrentQueue<TValue, TSpec> >
{
    typedef String<TValue> Type;
};

template <typename TValue, typename TSpec>
struct Size<ConcurrentQueue<TValue, TSpec> >:
    Size<Host<ConcurrentQueue<TValue, TSpec> > >
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function lockReading() / unlockReading()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void
lockReading(ConcurrentQueue<TValue, TSpec> & me)
{
    atomicInc(me.readerCount);
}

template <typename TValue, typename TSpec>
inline void
unlockReading(ConcurrentQueue<TValue, TSpec> & me)
{
    atomicDec(me.readerCount);
}

// ----------------------------------------------------------------------------
// Function lockWriting() / unlockWriting()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void
lockWriting(ConcurrentQueue<TValue, TSpec> & me)
{
    atomicInc(me.writerCount);
}

template <typename TValue, typename TSpec>
inline void
unlockWriting(ConcurrentQueue<TValue, TSpec> & me)
{
    atomicDec(me.writerCount);
}

// ----------------------------------------------------------------------------
// Function _cyclicInc()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue
_cyclicInc(TValue value, TValue modulo, TValue roundSize)
{
    TValue newVal = value + 1;
    if ((newVal & (roundSize - 1)) >= modulo)
        newVal += roundSize - modulo;
    return newVal;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool
empty(ConcurrentQueue<TValue, TSpec> const & me)
{
    ScopedWriteLock<> writeLock(me.lock);
    return me.headPos == me.tailPos;
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline typename Size<ConcurrentQueue<TValue, TSpec> >::Type
length(ConcurrentQueue<TValue, TSpec> const & me)
{
    typedef typename Size<ConcurrentQueue<TValue, TSpec> >::Type TSize;

    ScopedWriteLock<> writeLock(me.lock);
    TSize mask = me.roundSize - 1;
    if ((me.headPos & mask) <= (me.tailPos & mask))
        return me.tailPos - me.headPos;
    else
        return me.tailPos - me.headPos - (me.roundSize - capacity(me.data));
}

// ----------------------------------------------------------------------------
// Function capacity()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline typename Size<ConcurrentQueue<TValue, TSpec> >::Type
capacity(ConcurrentQueue<TValue, TSpec> const & me)
{
    ScopedReadLock<> writeLock(me.lock);
    return capacity(me.data) - 1;
}

// ----------------------------------------------------------------------------
// Function tryPopFront()
// ----------------------------------------------------------------------------
//
//  [  ?  ]  [  4  ]  [  3  ]  [  8  ]  [  0  ]  [  x  ]  [  ?  ]
//                       |                          ^
//                       v                          |
//             head            headRead   tail  tailWrite
//
// empty = (head == tail)
// full = (tail + 1 == head)
//
// valid data between  [headRead, tail)
// currently filled    [tail, tailWrite)
// currently removed   [head, headRead)

template <typename TValue, typename TSpec>
inline bool
tryPopFront(TValue & result, ConcurrentQueue<TValue, TSpec> & me)
{
    typedef ConcurrentQueue<TValue, TSpec>              TQueue;
    typedef typename Host<TQueue>::Type                 TString;
    typedef typename Size<TString>::Type                TSize;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    // try to extract a value
    ScopedReadLock<> readLock(me.lock);

    TSize cap = capacity(me.data);
    TSize roundSize = me.roundSize;
    TSize headReadPos;
    TSize newHeadReadPos;

    // wait for queue to become filled
    do {
        headReadPos = me.headReadPos;

        // return if queue is empty?
        if (headReadPos == me.tailPos)
            return false;

        newHeadReadPos = _cyclicInc(headReadPos, cap, roundSize);
    }
    while (!atomicCasBool(me.headReadPos, headReadPos, newHeadReadPos));

    // extract value and destruct it in the data string
    TIter it = begin(me.data, Standard()) + (headReadPos & (roundSize - 1));
    std::swap(result, *it);
    valueDestruct(it);

    atomicInc(me.popped);

    // wait for pending previous reads and synchronize headPos to headReadPos
    while (!atomicCasBool(me.headPos, headReadPos, newHeadReadPos))
    {}

    return true;
}

// ----------------------------------------------------------------------------
// Function waitForWriters()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void
waitForWriters(ConcurrentQueue<TValue, TSpec> & me, unsigned writerCount)
{
    while (me.writerCount < writerCount)
    {}
}

template <typename TValue, typename TSpec>
inline void
waitForFirstValue(ConcurrentQueue<TValue, TSpec> & me)
{
    while (!me.firstValue)
    {}
}

// ----------------------------------------------------------------------------
// Function popFront()
// ----------------------------------------------------------------------------

// returns if no writer is locked the queue and queue is empty
template <typename TValue, typename TSpec>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, TSpec> & me)
{
    while (me.writerCount != 0)
    {
        if (tryPopFront(result, me))
            return true;
    }
    // we have to give it another try if the queue was empty inside the loop
    // but after the check a writer pushes a value and zeroes the writerCount
    return (tryPopFront(result, me));
}

template <typename TValue, typename TSpec>
#ifdef SEQAN_CXX11_STANDARD
inline TValue &&
#else
inline TValue
#endif
popFront(ConcurrentQueue<TValue, TSpec> & me)
{
    TValue result;
    while (!tryPopFront(result, me))
    {}
#ifdef SEQAN_CXX11_STANDARD
    return std::move(result);
#else
    return result;
#endif
}

template <typename TValue, typename TSpec>
inline bool
#ifdef SEQAN_CXX11_STANDARD
_queueOverflow(ConcurrentQueue<TValue, TSpec> & me, TValue &&, Insist)
#else
_queueOverflow(ConcurrentQueue<TValue, TSpec> & me, TValue &, Insist)
#endif
{
    ignoreUnusedVariableWarning(me);
    SEQAN_ASSERT_GT(capacity(me.data), 1u);
    return false;
}

template <typename TValue, typename TSpec>
inline bool
#ifdef SEQAN_CXX11_STANDARD
_queueOverflow(ConcurrentQueue<TValue, TSpec> & me, TValue &&, Limit)
#else
_queueOverflow(ConcurrentQueue<TValue, TSpec> & me, TValue &, Limit)
#endif
{
    ignoreUnusedVariableWarning(me);
    SEQAN_ASSERT_GT(capacity(me.data), 1u);
    return false;
}

template <typename TValue, typename TSpec, typename TExpand>
inline bool
_queueOverflow(ConcurrentQueue<TValue, TSpec> & me,
#ifdef SEQAN_CXX11_STANDARD
               TValue && val,
#else
               TValue & val,
#endif
               Tag<TExpand> expandTag)
{
    typedef ConcurrentQueue<TValue, TSpec>              TQueue;
    typedef typename Host<TQueue>::Type                 TString;
    typedef typename Size<TString>::Type                TSize;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    // try to extend capacity

    ScopedWriteLock<> writeLock(me.lock);
    TSize cap = capacity(me.data);

    SEQAN_ASSERT_EQ(me.tailPos, me.tailWritePos);
    SEQAN_ASSERT_EQ(me.headPos, me.headReadPos);

    bool valueWasAppended = false;

    // did we reach the capacity limit (another thread could have done the upgrade already)?
    if (_cyclicInc(me.tailPos, cap, me.roundSize) >= me.headPos + me.roundSize)
    {
        if (cap != 0)
        {
            TIter it = begin(me.data, Standard()) + (me.tailPos & (me.roundSize - 1));
    //        valueConstruct(it, val, Move());
            valueConstruct(it);
            std::swap(*it, val);
            me.tailWritePos = me.tailPos = me.headPos + me.roundSize;
            me.firstValue = true;
            valueWasAppended = true;
            atomicInc(me.pushed);
        }

        SEQAN_ASSERT_EQ(me.tailPos, me.headPos + me.roundSize);

        // get positions of head/tail in current data sequence
        TSize headIdx = me.headPos & (me.roundSize - 1);
        TSize tailIdx = me.tailPos & (me.roundSize - 1);

        // increase capacity
        _setLength(me.data, cap);
        reserve(me.data, cap + 1, expandTag);
        TSize delta = capacity(me.data) - cap;
        me.roundSize = (TSize)1 << (log2(capacity(me.data) - 1) + 1);

        // create a gap of delta many values between tail and head
        _clearSpace(me.data, delta, headIdx, headIdx, expandTag);
        if (cap != 0)
        {
            me.headReadPos = me.headPos = headIdx + delta;
            me.tailWritePos = me.tailPos = tailIdx + me.roundSize;
        }
    }
    return valueWasAppended;
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------
// the queue is growing dynamically if expandTag is Generous,
// otherwise appendValue spinlocks until there is space to fill the value
template <typename TValue, typename TSpec, typename TExpand>
inline void
appendValue(ConcurrentQueue<TValue, TSpec> & me,
#ifdef SEQAN_CXX11_STANDARD
            TValue && val,
#else
            TValue val,
#endif
            Tag<TExpand> expandTag)
{
    typedef ConcurrentQueue<TValue, TSpec>              TQueue;
    typedef typename Host<TQueue>::Type                 TString;
    typedef typename Size<TString>::Type                TSize;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    while (true)
    {
        // try to append the value
        {
            ScopedReadLock<> readLock(me.lock);
            TSize cap = capacity(me.data);
            TSize roundSize = me.roundSize;

            while (true)
            {
                TSize tailWritePos = me.tailWritePos;
                TSize newTailWritePos = _cyclicInc(tailWritePos, cap, roundSize);
                if (newTailWritePos >= me.headPos + roundSize)
                    break;

                if (atomicCasBool(me.tailWritePos, tailWritePos, newTailWritePos))
                {
                    TIter it = begin(me.data, Standard()) + (tailWritePos & (roundSize - 1));
//                    valueConstruct(it, val, Move());
                    valueConstruct(it);
                    std::swap(*it, val);

                    // wait for pending previous writes and synchronize tailPos to tailWritePos
                    while (!atomicCasBool(me.tailPos, tailWritePos, newTailWritePos))
                    {}

                    me.firstValue = true;
                    atomicInc(me.pushed);
                    return;
                }
            }
        }

        if (_queueOverflow(me, val, expandTag))
            return;
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_QUEUE_H_
