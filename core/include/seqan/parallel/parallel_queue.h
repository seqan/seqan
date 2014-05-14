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

    ConcurrentQueue() :
        headPos(0),
        headReadPos(0),
        tailPos(0),
        tailWritePos(0)
    {}

    // you can set the initial capacity here
    ConcurrentQueue(TSize capacity) :
        headPos(0),
        headReadPos(0),
        tailPos(0),
        tailWritePos(0)
    {
        reserve(data, capacity + 1, Exact());
    }

    ConcurrentQueue(TString & data):
        data(data),
        headPos(0),
        headReadPos(0),
        tailPos(length(data)),
        tailWritePos(length(data))
    {}

    ~ConcurrentQueue()
    {
        SEQAN_ASSERT_EQ(tailPos, tailWritePos);
        SEQAN_ASSERT_EQ(headPos, headReadPos);
        SEQAN_ASSERT(empty(lock));

        TSize cap = capacity(data);
        if (tailPos < headPos)
        {
            _clearSpace(data, 0u, tailPos, headPos, Insist());
            _setLength(data, cap - (headPos - tailPos));
        }
        else
        {
            _clearSpace(data, 0u, (TSize)0, headPos, Insist());
            _setLength(data, cap - headPos);
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

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _cyclicInc()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue
_cyclicInc(TValue value, TValue modulo)
{
    TValue newVal = value + 1;
    return (newVal >= modulo)? 0 : newVal;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool
empty(ConcurrentQueue<TValue, TSpec> const & me)
{
    ScopedWriteLock(me.lock);
    return me.headPos == me.tailPos;
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline typename Size<ConcurrentQueue<TValue, TSpec> >::Type
length(ConcurrentQueue<TValue, TSpec> const & me)
{
    ScopedWriteLock(me.lock);
    if (me.headPos <= me.tailPos)
        return me.tailPos - me.headPos;
    else
        return capacity(me.data) - me.headPos + me.tailPos;
}

// ----------------------------------------------------------------------------
// Function capacity()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline typename Size<ConcurrentQueue<TValue, TSpec> >::Type
capacity(ConcurrentQueue<TValue, TSpec> const & me)
{
    return capacity(me.data);
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
tryPopFront(TValue & result,
            ConcurrentQueue<TValue, TSpec> & me)
{
    typedef ConcurrentQueue<TValue, TSpec>              TQueue;
    typedef typename Host<TQueue>::Type                 TString;
    typedef typename Size<TString>::Type                TSize;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    // try to extract a value
    ScopedReadLock(me.lock);

    TSize cap = capacity(me.data);
    TSize headReadPos;
    TSize newHeadReadPos;

    // wait for queue to become filled
    do {
        headReadPos = me.headReadPos;

        // return if queue is empty?
        if (headReadPos == me.tailPos)
            return false;

        newHeadReadPos = _cyclicInc(headReadPos, cap);
    }
    while (atomicCas(me.headReadPos, headReadPos, newHeadReadPos) != headReadPos);

    // extract value and destruct it in the data string
    TIter it = begin(me.data, Standard()) + headReadPos;
    std::swap(result, *it);
    valueDestruct(it);

    // wait for pending previous reads and synchronize headPos to headReadPos
    while (atomicCas(me.headPos, headReadPos, newHeadReadPos) != headReadPos);

    return true;
}

// ----------------------------------------------------------------------------
// Function popFront()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
#ifdef SEQAN_CXX11_STANDARD
inline TValue &&
#else
inline TValue
#endif
popFront(ConcurrentQueue<TValue, TSpec> & me)
{
    TValue val;
    while (!tryPopFront(val, me));
#ifdef SEQAN_CXX11_STANDARD
    return std::move(val);
#else
    return val;
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
    SEQAN_ASSERT_NEQ(capacity(me), 0u);
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
    SEQAN_ASSERT_NEQ(capacity(me), 0u);
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

    ScopedWriteLock(me.lock);
    TSize cap = capacity(me.data);

    SEQAN_ASSERT_EQ(me.tailPos, me.tailWritePos);
    SEQAN_ASSERT_EQ(me.headPos, me.headReadPos);

    bool valueWasAppended = false;

    // did we reach the capacity limit (another thread could have done the upgrade already)?
    if (_cyclicInc(me.tailPos, cap) == me.headPos)
    {
        if (cap != 0)
        {
            TIter it = begin(me.data, Standard()) + me.tailPos;
    //        valueConstruct(it, val, Move());
            valueConstruct(it);
            std::swap(*it, val);
            me.tailWritePos = me.tailPos = me.headPos;
            valueWasAppended = true;
        }

        SEQAN_ASSERT_EQ(me.tailPos, me.headPos);

        // increase capacity
        _setLength(me.data, cap);
        reserve(me.data, cap + 1, expandTag);
        TSize delta = capacity(me.data) - cap;

        // create a gap of delta many values between tail and head
        _clearSpace(me.data, delta, me.headPos, me.headPos, expandTag);
        if (cap != 0)
        {
            me.headPos += delta;
            me.headReadPos = me.headPos;
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
            ScopedReadLock(me.lock);

            TSize cap = capacity(me.data);

            while (true)
            {
                TSize tailWritePos = me.tailWritePos;
                TSize newTailWritePos = _cyclicInc(tailWritePos, cap);
                if (newTailWritePos == me.headPos)
                    break;

                if (atomicCas(me.tailWritePos, tailWritePos, newTailWritePos) == tailWritePos)
                {
                    TIter it = begin(me.data, Standard()) + tailWritePos;
//                    valueConstruct(it, val, Move());
                    valueConstruct(it);
                    std::swap(*it, val);

                    // wait for pending previous writes and synchronize tailPos to tailWritePos
                    while (atomicCas(me.tailPos, tailWritePos, newTailWritePos) != tailWritePos);
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
