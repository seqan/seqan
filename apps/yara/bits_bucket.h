// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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
// This file contains functions/ors to bucket (i.e. index) a container by key.
// ==========================================================================

#ifndef APP_YARA_BITS_BUCKET_H_
#define APP_YARA_BITS_BUCKET_H_

namespace seqan {

// ============================================================================
// Functors
// ============================================================================

// ----------------------------------------------------------------------------
// Class Adder
// ----------------------------------------------------------------------------

template <typename TUnaryFunction, unsigned DELTA>
struct Adder
{
    TUnaryFunction const & f;

    Adder(TUnaryFunction const & f) : f(f) {}

    template <typename TValue>
    unsigned operator() (TValue const & val) const
    {
        return f(val) + DELTA;
    }
};

// ----------------------------------------------------------------------------
// Class KeyIndicator
// ----------------------------------------------------------------------------

template <typename TTarget, typename TKey, typename TSpec = void>
struct KeyIndicator
{
    TTarget &       target;
    TKey const &    key;

    KeyIndicator(TTarget & target, TKey const & key) :
        target(target),
        key(key)
    {}

    template <typename TValue>
    void operator() (TValue const & val) const
    {
        SEQAN_ASSERT_LT(key(val), length(target));
        target[key(val)] = true;
    }
};

// ----------------------------------------------------------------------------
// Class KeyCounter
// ----------------------------------------------------------------------------

template <typename TTarget, typename TKey, typename TThreading = Serial, typename TSpec = void>
struct KeyCounter
{
    TTarget &       target;
    TKey const &    key;

    KeyCounter(TTarget & target, TKey const & key) :
        target(target),
        key(key)
    {}

    template <typename TValue>
    void operator() (TValue const & val) const
    {
        SEQAN_ASSERT_LT(key(val), length(target));
        atomicInc(target[key(val)], TThreading());
    }
};

// ----------------------------------------------------------------------------
// Class KeySorter
// ----------------------------------------------------------------------------

template <typename TSource, typename TSpec = void>
struct KeySorter
{
    TSource const & source;

    KeySorter(TSource const & source) :
        source(source)
    {}

    template <typename TKey>
    inline bool operator()(TKey a, TKey b) const
    {
        return source[a] < source[b];
    }
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function bucket()
// --------------------------------------------------------------------------
// Bucket elements in the concat of a ConcatDirect StringSet.
// Remarks: the concat string must be already sorted by key.

template <typename TString, typename TSpec, typename TKeyGetter, typename TMaxKey, typename TThreading>
inline void bucket(StringSet<TString, Owner<ConcatDirect<TSpec > > > & me, TKeyGetter const & key, TMaxKey maxKey, TThreading const & tag)
{
    typedef StringSet<TString, Owner<ConcatDirect<TSpec > > >    TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type           TLimits;
    typedef Adder<TKeyGetter, 1u>                                TNextKey;
    typedef KeyCounter<TLimits, TNextKey, TThreading>            TCounter;

    // Shift the counts by one.
    TNextKey nextKey(key);

    // Resize the limits string to count all keys.
    resize(stringSetLimits(me), maxKey, 0, Exact());

    // Count the number of keys present in the concat string.
    forEach(concat(me), TCounter(stringSetLimits(me), nextKey), tag);

    // Build the limits string from the key counts.
    partialSum(stringSetLimits(me), tag);
}

// --------------------------------------------------------------------------
// Function bucket()
// --------------------------------------------------------------------------
// Bucket elements in the host of a Segment StringSet.
// Remarks: the host string must be already sorted by key.

template <typename THost, typename TSpec, typename TKeyGetter, typename TMaxKey, typename TThreading>
inline void bucket(StringSet<THost, Segment<TSpec> > & me, TKeyGetter const & key, TMaxKey maxKey, TThreading const & tag)
{
    typedef StringSet<THost, Segment<TSpec> >                    TStringSet;
    typedef typename StringSetLimits<TStringSet>::Type           TLimits;
    typedef Adder<TKeyGetter, 1u>                                TNextKey;
    typedef KeyCounter<TLimits, TNextKey, TThreading>            TCounter;

    // Shift the key counts by one.
    TNextKey nextKey(key);

    // Resize the limits string to accomodate counts for all keys.
    resize(stringSetLimits(me), maxKey + 1, 0, Exact());

    // Count the number of keys present in the host string.
    forEach(host(me), TCounter(stringSetLimits(me), nextKey), tag);

    // Limits are the cumulated key counts.
    partialSum(stringSetLimits(me), tag);

    // Positions are the shifted limits.
    assign(stringSetPositions(me), prefix(stringSetLimits(me), length(stringSetLimits(me)) - 1), Exact());
}

}

#endif  // #ifndef APP_YARA_BITS_BUCKET_H_
