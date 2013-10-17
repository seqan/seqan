// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
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

#ifndef INDEX_FM_RANK_DICTIONARY_BASE_H_
#define INDEX_FM_RANK_DICTIONARY_BASE_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag FibreRanks
// ----------------------------------------------------------------------------

struct FibreRanks_;

typedef Tag<FibreRanks_>
const FibreRanks;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class RankDictionary
// ----------------------------------------------------------------------------

/**
.Class.RankDictionary:
..cat:Index
..summary:A rank dictionary is a data structure to store the rank of an element of a sequence at every position of the 
sequence.
..signature:RankDictionary<TValue, TSpec>
..param.TSpec:The rank dictionary specialisation.
...type:Spec.WaveletTree
...type:Spec.SequenceBitMask
...default:@Spec.WaveletTree@
..include:seqan/index.h
*/
template <typename TValue, typename TSpec>
struct RankDictionary;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------
 
template <typename TValue, typename TSpec>
struct Size<RankDictionary<TValue, TSpec> >
{
    typedef typename Size<String<TValue, TSpec> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Value<RankDictionary<TValue, TSpec> >
{
    typedef TValue  Type;
};

template <typename TValue, typename TSpec>
struct Value<RankDictionary<TValue, TSpec> const> :
    Value<RankDictionary<TValue, TSpec> > {};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryFibreSpec
// ----------------------------------------------------------------------------

template <typename TRankDictionary>
struct RankDictionaryFibreSpec
{
    typedef Alloc<> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
SEQAN_HOST_DEVICE inline typename Fibre<RankDictionary<TValue, TSpec>, FibreRanks>::Type &
getFibre(RankDictionary<TValue, TSpec> & dict, FibreRanks)
{
    return dict.ranks;
}

template <typename TValue, typename TSpec>
SEQAN_HOST_DEVICE inline typename Fibre<RankDictionary<TValue, TSpec>, FibreRanks>::Type const &
getFibre(RankDictionary<TValue, TSpec> const & dict, FibreRanks)
{
    return dict.ranks;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void clear(RankDictionary<TValue, TSpec> & dict)
{
    clear(getFibre(dict, FibreRanks()));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
SEQAN_HOST_DEVICE inline bool empty(RankDictionary<TValue, TSpec> const & dict)
{
    return empty(getFibre(dict, FibreRanks()));
}

// ----------------------------------------------------------------------------
// Function createRankDictionary()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TText>
inline void
createRankDictionary(RankDictionary<TValue, TSpec> & dict, TText const & text)
{
    typedef typename Iterator<TText const, Standard>::Type      TTextIterator;

    // Resize the RankDictionary.
    resize(dict, length(text), Exact());

    // Assign the text value by value.
    TTextIterator textBegin = begin(text, Standard());
    TTextIterator textEnd = end(text, Standard());
    for (TTextIterator textIt = textBegin; textIt != textEnd; ++textIt)
        setValue(dict, textIt - textBegin, value(textIt));

    // Update all ranks.
    updateRanks(dict);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool open(RankDictionary<TValue, TSpec> & dict, const char * fileName, int openMode)
{
    return open(getFibre(dict, FibreRanks()), fileName, openMode);
}

template <typename TValue, typename TSpec>
inline bool open(RankDictionary<TValue, TSpec> & dict, const char * fileName)
{
    return open(dict, fileName, DefaultOpenMode<RankDictionary<TValue, TSpec> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool save(RankDictionary<TValue, TSpec> const & dict, const char * fileName, int openMode)
{
    return save(getFibre(dict, FibreRanks()), fileName, openMode);
}

template <typename TValue, typename TSpec>
inline bool save(RankDictionary<TValue, TSpec> const & dict, const char * fileName)
{
    return save(dict, fileName, DefaultOpenMode<RankDictionary<TValue, TSpec> >::VALUE);
}

}

#endif  // INDEX_FM_RANK_DICTIONARY_BASE_H_
