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

template <typename TString, typename TSpec = void>
struct ConcurrentQueue
{
    typedef typename Position<TString>::Type TPos;

    TString *       data_host;
    ReadWriteLock   lock;


    volatile TPos headPos;
    volatile TPos headReadPos;
    volatile TPos tailPos;
    volatile TPos tailWritePos;

    ConcurrentQueue() :
        data_host(),
        headPos(0),
        tailPos(0),
        tailWritePos(0)
    {}

    ConcurrentQueue(TString & string):
        headPos(0),
        tailPos(0),
        tailWritePos(0)
    {
        setHost(*this, string);
    }

    ~ConcurrentQueue()
    {
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Value<ConcurrentQueue<TString, TSpec> > : Value<TString> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline TString &
host(ConcurrentQueue<TString, TSpec> & me)
{
    return *me.data_host;
}

template <typename TString, typename TSpec>
inline TString const &
host(ConcurrentQueue<TString, TSpec> const & me)
{
    return *me.data_host;
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline void
setHost(ConcurrentQueue<TString, TSpec> & me, TString & string)
{
    me.data_host = &string;
    me.tailPos = length(string);
}

// ----------------------------------------------------------------------------
// Function _cyclicInc()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue
cyclicInc(TValue value, TValue modulo)
{
    TValue newVal = value + 1;
    return (newVal == modulo)? 0 : newVal;
}

template <typename TValue>
inline TValue
atomicCyclicInc(TValue volatile & value, TValue modulo)
{
    do {
        TValue curVal = value;
        TValue newVal = cyclicInc(curVal);
    }
    while (atomicCas(value, curVal, newVal) != curVal);
    return newVal;
}

// ----------------------------------------------------------------------------
// Function _cyclicDec()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue
cyclicDec(TValue value, TValue modulo)
{
    TValue newVal = (value == 0)? modulo : value;
    return newVal - 1;
}

template <typename TValue>
inline TValue
atomicCyclicDec(TValue volatile & value, TValue modulo)
{
    do {
        TValue curVal = value;
        TValue newVal = cyclicDec(curVal);
    }
    while (atomicCas(value, curVal, newVal) != curVal);
    return newVal;
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TValue, typename TExpand, typename TParallel>
inline void
popFront(ConcurrentQueue<TString, TSpec> & me,
            TValue const & val,
            Tag<TExpand> const & expandTag,
            Tag<TParallel> const & parallelTag)
{
    appendValue(host(me), val, expandTag, parallelTag);
}

template <typename TString, typename TSpec, typename TValue, typename TExpand, typename TParallel>
inline void
appendValue(ConcurrentQueue<TString, TSpec> & me,
            TValue const & val,
            Tag<TExpand> const & expandTag,
            Tag<TParallel> const & parallelTag)
{
    appendValue(host(me), val, expandTag, parallelTag);
}

//
//  [==0==] [==1==] [==2==] [==3==]
//             ^       ^       ^
//            head    tail  tailWrite
//
//
// empty = (head == tail)
// full = (tail + 1 == head)
//
// valid data between [head, tail)

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

    TSize headReadPos;
    TSize newHeadReadPos;

    // wait for queue to become filled
    do {
        headReadPos = me.headReadPos;

        // return if queue is empty?
        if (headReadPos == me.tailPos)
            return false;

        newHeadReadPos = cyclicInc(headReadPos, cap);
    }
    while (atomicCas(me.headReadPos, headReadPos, newHeadReadPos) != headReadPos);

    // extract value and destruct it in the string
    TString & string = host(me);
    TIter it = begin(string, Standard()) + headReadPos;
    swap(result, *it);
    valueDestruct(it);

    // wait for pending previous reads and synchronize headPos to headReadPos
    while (atomicCas(me.headPos, headReadPos, newHeadReadPos) != headReadPos);

    return true;
}

template <typename TValue, typename TString, typename TSpec, typename TParallel>
inline TValue &&
popFront(ConcurrentQueue<TString, TSpec> & me,
         Tag<TParallel> parallelTag)
{
    TValue val;
    while (!tryPopFront(val, me, parallelTag));
    return std::move(val);
}

template <typename TString, typename TSpec, typename TValue, typename TExpand>
inline void
appendValue(ConcurrentQueue<TString, TSpec> & me,
            TValue && val,
            Tag<TExpand> const & expandTag,
            Parallel)
{
    typedef typename Size<TString>::Type TSize;

    TString & string = host(me);

    while (true)
    {
        // try to append the value
        {
            ScopedReadLock(me.lock);

            TSize cap = capacity(string);
            SEQAN_ASSERT_NEQ(cap, 0u);

            TSize newTailWritePos = atomicCyclicInc(me.tailWritePos, cap);
            TSize tailWritePos = cyclicDec(newTailWritePos);

            if (newTailWritePos != me.headPos)
            {
                valueConstruct(begin(string, Standard()) + tailWritePos, val);
                // wait for pending previous writes and synchronize tailPos to tailWritePos
                while (atomicCas(me.tailPos, tailWritePos, newTailWritePos) != tailWritePos);
                return;
            }
            atomicCyclicDec(me.tailWritePos, cap);
        }

        // try to extend capacity
        {
            ScopedWriteLock(me.lock);
            TSize cap = capacity(string);

            // did we reach the capacity limit?
            if (cyclicInc(me.tailPos) == me.headPos)
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

template <typename TTargetValue, typename TTargetSpec, typename TValue>
inline void
appendValue(String<TTargetValue, TTargetSpec> & me, TValue const & _value, Insist, Parallel)
{
    valueConstruct(begin(me, Standard()) + _incLength(me, Parallel()) - 1, _value);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_QUEUE_H_
