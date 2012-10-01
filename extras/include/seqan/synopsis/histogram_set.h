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

#ifndef SEQAN_SYNOPSIS_HISTOGRAM_SET_H_
#define SEQAN_SYNOPSIS_HISTOGRAM_SET_H_

#include <algorithm>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.HistogramSet
..cat:Synopsis Data Structures
..summary:Data structure for keeping a compact set of histograms.
..signature:HistogramSet<TValue>
..param.TValue:Type of the histogram counters.
..include:seqan/synopsis.h
 */

template <typename TValue>
class HistogramSet
{
public:
    typedef StringSet<String<TValue>, Owner<ConcatDirect<> > > TStringSet_;

    unsigned _count;          // Number of histograms.
    unsigned _histogramSize;  // Number of entries in each histogram.
    TStringSet_ _stringSet;   // Hosting string set.

    HistogramSet() {}

    HistogramSet(unsigned count, unsigned size)
            : _count(count), _histogramSize(size)
    {
        resize(_stringSet.concat, _histogramSize * _count, 0);
        resize(_stringSet.limits, _count);
        for (unsigned i = 0, sum = 0; i < _count; ++i) {
            _stringSet.limits[i] = sum;
            sum += _histogramSize;
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T:Class.HistogramSet
///.Metafunction.Value.class:Class.HistogramSet

template <typename TValue>
struct Value<HistogramSet<TValue> >
{
    typedef TValue Type;
};

template <typename TValue>
struct Value<HistogramSet<TValue> const>
{
    typedef TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

///.Metafunction.Size.param.T:Class.HotList
///.Metafunction.Size.class:Class.HotList

template <typename TValue>
struct Size<HistogramSet<TValue> >
{
    typedef unsigned Type;
};

template <typename TValue>
struct Size<HistogramSet<TValue> const>
{
    typedef unsigned Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

///.Metafunction.Position.param.T:Class.HotList
///.Metafunction.Position.class:Class.HotList

template <typename TValue>
struct Position<HistogramSet<TValue> >
{
    typedef unsigned Type;
};

template <typename TValue>
struct Position<HistogramSet<TValue> const>
{
    typedef unsigned Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
clear(HistogramSet<TValue> & hSet)
{
    ::std::fill(begin(hSet._stringSet.concat, Standard()), end(hSet._stringSet.concat, Standard()), 0);
}

template <typename TValue>
inline void
clear(HistogramSet<TValue> & hSet, unsigned idx)
{
    ::std::fill(begin(hSet._stringSet[idx], Standard()), end(hSet._stringSet[idx], Standard()), 0);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TValue>
inline typename Reference<HistogramSet<TValue> const>::Type
value(HistogramSet<TValue> const & hSet, unsigned hIdx, unsigned entryIdx)
{
    SEQAN_ASSERT_LEQ(hIdx, hSet._count);
    SEQAN_ASSERT_LEQ(entryIdx, hSet._histogramSize);
    return hSet._stringSet[hIdx][entryIdx];
}

template <typename TValue>
inline typename Reference<HistogramSet<TValue> >::Type
value(HistogramSet<TValue> & hSet, unsigned hIdx, unsigned entryIdx)
{
    SEQAN_ASSERT_LEQ(hIdx, hSet._count);
    SEQAN_ASSERT_LEQ(entryIdx, hSet._histogramSize);
    return hSet._stringSet[hIdx][entryIdx];
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue
getValue(HistogramSet<TValue> const & hSet, unsigned hIdx, unsigned entryIdx)
{
    SEQAN_ASSERT_LEQ(hIdx, hSet._count);
    SEQAN_ASSERT_LEQ(entryIdx, hSet._histogramSize);
    return hSet._stringSet[hIdx][entryIdx];
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TValue2>
inline void
setValue(HistogramSet<TValue> & hSet, unsigned hIdx, unsigned entryIdx, TValue2 const & x)
{
    SEQAN_ASSERT_LEQ(hIdx, hSet._count);
    SEQAN_ASSERT_LEQ(entryIdx, hSet._histogramSize);
    hSet._stringSet[hIdx][entryIdx] = x;
}

// ----------------------------------------------------------------------------
// Function incrementValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TValue2>
inline void
incrementValue(HistogramSet<TValue> & hSet, unsigned hIdx, unsigned entryIdx, TValue2 const & x)
{
    SEQAN_ASSERT_LEQ(hIdx, hSet._count);
    SEQAN_ASSERT_LEQ(entryIdx, hSet._histogramSize);
    hSet._stringSet[hIdx][entryIdx] += x;
}

}  // namespace seqan

#endif  // SEQAN_SYNOPSIS_HISTOGRAM_SET_H_
