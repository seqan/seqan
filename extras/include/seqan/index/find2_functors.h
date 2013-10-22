// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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
// Some useful functors to be used along with the Finder class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_INDEX_FIND_FUNCTORS_H_
#define SEQAN_EXTRAS_INDEX_FIND_FUNCTORS_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Counts_
// ----------------------------------------------------------------------------

struct Counts_;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class OccurrencesCounter
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec = typename ExecSpec<TIndex>::Type>
struct OccurrencesCounter
{
    typename Member<OccurrencesCounter, Counts_>::Type    counts;

    OccurrencesCounter() {}

    template <typename TPattern>
    OccurrencesCounter(TPattern const & pattern)
    {
        init(*this, pattern);
    }

    template <typename TFinder>
    inline SEQAN_HOST_DEVICE void
    operator() (TFinder const & finder)
    {
        counts[getThreadId()] += countOccurrences(textIterator(finder));
    }
};

// ----------------------------------------------------------------------------
// Member Counts_; Count
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct Member<OccurrencesCounter<TIndex, TSpec>, Counts_>
{
    typedef String<typename Size<TIndex>::Type> Type;
};

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
struct Member<OccurrencesCounter<TIndex, Device<TSpec> >, Counts_>
{
    typedef thrust::device_vector<typename Size<TIndex>::Type>  Type;
};
#endif

template <typename TIndex, typename TSpec>
struct Member<OccurrencesCounter<TIndex, View<Device<TSpec> > >, Counts_>
{
    typedef typename View<typename Member<OccurrencesCounter<TIndex, Device<TSpec> >, Counts_>::Type>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct View<OccurrencesCounter<TIndex, TSpec> >
{
    typedef OccurrencesCounter<TIndex, View<TSpec> >  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Device
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct Device<OccurrencesCounter<TIndex, TSpec> >
{
    typedef OccurrencesCounter<TIndex, Device<TSpec> >  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function init()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TPattern>
inline void
init(OccurrencesCounter<TIndex, TSpec> & counter, TPattern const & /* pattern */)
{
    resize(counter.counts, omp_get_max_threads(), 0, Exact());
}

template <typename TIndex, typename TSpec, typename TPattern>
inline void
init(OccurrencesCounter<TIndex, Device<TSpec> > & counter, TPattern const & pattern)
{
    resize(counter.counts, length(needle(pattern)), Exact());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
typename View<OccurrencesCounter<TIndex, TSpec> >::Type
view(OccurrencesCounter<TIndex, TSpec> & counter)
{
    typename View<OccurrencesCounter<TIndex, TSpec> >::Type counterView;

    counterView.counts = view(counter.counts);

    return counterView;
}

// ----------------------------------------------------------------------------
// Function getCount()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
inline typename Size<TIndex>::Type
getCount(OccurrencesCounter<TIndex, TSpec> & counter)
{
    typedef typename Size<TIndex>::Type TSize;

    // TODO(esiragusa): Add function reduce().
    TSize count = 0;
    for (TSize i = 0; i < length(counter.counts); ++i)
        count += counter.counts[i];

    return count;
}

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
inline typename Size<TIndex>::Type
getCount(OccurrencesCounter<TIndex, Device<TSpec> > & counter)
{
    return thrust::reduce(begin(counter.counts, Standard()), end(counter.counts, Standard()));
}
#endif

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_FIND_FUNCTORS_H_
