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
// Thread-safe queue of either fixed/variable size
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

    TString         data;
    ReadWriteLock   lock;

    volatile TSize headPos;
    volatile TSize headReadPos;
    volatile TSize tailPos;
    volatile TSize tailWritePos;

    ConcurrentQueue() :
        headPos(0),
        tailPos(0),
        tailWritePos(0)
    {}

    ConcurrentQueue(TString & data):
        data(data),
        headPos(0),
        tailPos(0),
        tailWritePos(0)
    {}

    ~ConcurrentQueue()
    {
        SEQAN_ASSERT_EQ(tailPos, tailWritePos);
        SEQAN_ASSERT_EQ(headPos, headReadPos);
        SEQAN_ASSERT(empty(lock));

        TSize cap = capacity(data);
        if (tailPos < headPos)
        {
            _clearSpace(data, 0, tailPos, headPos);
            _setLength(data, cap - (headPos - tailPos));
        }
        else
        {
            _clearSpace(data, 0, 0, headPos);
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
    return (newVal == modulo)? 0 : newVal;
}

// ----------------------------------------------------------------------------
// Function _cyclicDec()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue
_cyclicDec(TValue value, TValue modulo)
{
    TValue newVal = (value == 0)? modulo : value;
    return newVal - 1;
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

template <typename TValue, typename TString, typename TSpec>
inline bool
tryPopFront(TValue & result,
            ConcurrentQueue<TString, TSpec> & me,
            Parallel)
{
    typedef typename Size<TString>::Type                TSize;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    // try to extract a value
    ScopedReadLock(me.lock);

    TString & string = host(me);

    TSize cap = capacity(string);
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

    // extract value and destruct it in the string
    TIter it = begin(string, Standard()) + headReadPos;
    swap(result, *it);
    valueDestruct(it);

    // wait for pending previous reads and synchronize headPos to headReadPos
    while (atomicCas(me.headPos, headReadPos, newHeadReadPos) != headReadPos);

    return true;
}

// ----------------------------------------------------------------------------
// Function popFront()
// ----------------------------------------------------------------------------

template <typename TValue, typename TString, typename TSpec, typename TParallel>
#ifdef SEQAN_CXX11_STANDARD
inline TValue &&
#else
inline TValue
#endif
popFront(ConcurrentQueue<TString, TSpec> & me,
         Tag<TParallel> parallelTag)
{
    TValue val;
    while (!tryPopFront(val, me, parallelTag));
#ifdef SEQAN_CXX11_STANDARD
    return std::move(val);
#else
    return val;
#endif
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TValue, typename TExpand, typename TParallel>
inline void
appendValue(ConcurrentQueue<TString, TSpec> & me,
            TValue const & val,
            Tag<TExpand> const & expandTag,
            Tag<TParallel> const & parallelTag)
{
    appendValue(host(me), val, expandTag, parallelTag);
}

template <typename TString, typename TSpec, typename TValue, typename TExpand>
inline void
appendValue(ConcurrentQueue<TString, TSpec> & me,
#ifdef SEQAN_CXX11_STANDARD
            TValue && val,
#else
            TValue val,
#endif
            Tag<TExpand> const & expandTag,
            Parallel)
{
    typedef typename Size<TString>::Type                TSize;
    typedef typename Iterator<TString, Standard>::Type  TIter;

    TString & string = host(me);

    while (true)
    {
        // try to append the value
        {
            ScopedReadLock(me.lock);

            TSize cap = capacity(string);
            SEQAN_ASSERT_NEQ(cap, 0u);

            while (true)
            {
                TValue tailWritePos = me.tailWritePos;
                TValue newTailWritePos = _cyclicInc(tailWritePos, cap);
                if (newTailWritePos == me.headPos)
                    break;

                if (atomicCas(me.tailWritePos, tailWritePos, newTailWritePos) == tailWritePos)
                {
                    TIter it = begin(string, Standard()) + tailWritePos;
                    valueConstruct(it);
                    swap(*it, val);

                    // wait for pending previous writes and synchronize tailPos to tailWritePos
                    while (atomicCas(me.tailPos, tailWritePos, newTailWritePos) != tailWritePos);
                    return;
                }
            }
        }

        // try to extend capacity
        {
            ScopedWriteLock(me.lock);
            TSize cap = capacity(string);

            SEQAN_ASSERT_EQ(me.tailPos, me.tailWritePos);
            SEQAN_ASSERT_EQ(me.headPos, me.headReadPos);

            // did we reach the capacity limit?
            if (_cyclicInc(me.tailPos) == me.headPos)
            {
                valueConstruct(begin(string, Standard()) + me.tailPos, val);
                me.tailWritePos = me.tailPos = me.headPos;

                // increase capacity
                _setLength(string, cap);
                reserve(string, cap + 1, expandTag);
                TSize delta = capacity(string) - cap;

                // create a gap of delta many values between tail and head
                _clearSpace(string, delta, me.headPos, me.headPos, expandTag);
                me.headPos += delta;
                me.headReadPos = me.headPos;
            }
        }
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_QUEUE_H_
