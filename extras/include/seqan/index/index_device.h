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
// This file contains Index specializations for thrust::device_vector.
// ==========================================================================

#ifndef SEQAN_EXTRAS_INDEX_DEVICE_H_
#define SEQAN_EXTRAS_INDEX_DEVICE_H_

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Device                                                  [Index]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct Device<Index<TText, TSpec> >
{
    typedef Index<typename Device<TText>::Type, TSpec>          Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsDevice                                                [Index]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
struct IsDevice<Index<TText, TSpec> > : IsDevice<TText> {};

// ----------------------------------------------------------------------------
// Metafunction IsDevice                                                 [Iter]
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct IsDevice<Iter<TIndex, VSTree<TSpec> > > : IsDevice<TIndex> {};

// ----------------------------------------------------------------------------
// Metafunction FibreSA                                          [Device Index]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreSA>
{
    typedef thrust::device_vector<typename SAValue<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type> Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreSA>
{
    typedef thrust::device_vector<typename SAValue<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec> >::Type>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreLcp                                         [Device Index]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreLcp>
{
    typedef thrust::device_vector<typename Size<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>    Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreLcp>
{
    typedef thrust::device_vector<typename Size<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec> >::Type>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreChildtab                                    [Device Index]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreChildtab>
{
    typedef thrust::device_vector<typename Size<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>    Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreChildtab>
{
    typedef thrust::device_vector<typename Size<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec> >::Type>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreBwt                                         [Device Index]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreBwt>
{
    typedef thrust::device_vector<typename Value<Index<thrust::device_vector<TValue, TAlloc>, TSpec> >::Type>   Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreBwt>
{
    typedef thrust::device_vector<typename Value<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec> >::Type> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreSA                                        [Device FMIndex]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TOccSpec, typename TSpec>
struct Fibre<Index<thrust::device_vector<TValue, TAlloc>, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
    typedef CompressedSA<thrust::device_vector<TValue, TAlloc>, TSpec>      Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TOccSpec, typename TSpec>
struct Fibre<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
    typedef CompressedSA<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibrePrefixSum                                 [Device LfTable]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<LfTable<thrust::device_vector<TValue, TAlloc>, TSpec>, FibrePrefixSum>
{
    typedef thrust::device_vector<TValue, TAlloc>               TText_;
    typedef typename Value<LfTable<TText_, TSpec> >::Type       TValue_;
    typedef typename MakeUnsigned<TValue_>::Type                TUValue_;

    typedef PrefixSumTable<TUValue_, Device<TSpec> >            Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<LfTable<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibrePrefixSum>
{
    typedef StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec> TText_;
    typedef typename Value<LfTable<TText_, TSpec> >::Type               TValue_;
    typedef typename MakeUnsigned<TValue_>::Type                        TUValue_;

    typedef PrefixSumTable<TUValue_, Device<TSpec> >                    Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreValues                                    [Device LfTable]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<LfTable<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreValues>
{
    typedef thrust::device_vector<TValue, TAlloc>               TText_;
    typedef typename Value<LfTable<TText_, TSpec> >::Type       TValue_;

    typedef RankDictionary<TwoLevels<TValue_, Device<TSpec> > > Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<LfTable<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreValues>
{
    typedef StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec> TText_;
    typedef typename Value<LfTable<TText_, TSpec> >::Type               TValue_;

    typedef RankDictionary<TwoLevels<TValue_, Device<TSpec> > >         Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreSentinels                                 [Device LfTable]
// ----------------------------------------------------------------------------

// NOTE(esiragusa): Single text sentinel rank dictionary is an integer.
//template <typename TText, typename TViewSpec, typename TSpec>
//struct Fibre<LfTable<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreSentinels>
//{
//    typedef typename Fibre<LfTable<TText, TSpec>, FibreSentinels>::Type   Type;
//};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<LfTable<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreSentinels>
{
//    typedef RankDictionary<TwoLevels<bool, Device<TSpec> > >    Type;
    typedef RankDictionary<Naive<bool, Device<TSpec> > >    Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreEntries                            [Device PrefixSumTable]
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec>
struct Fibre<PrefixSumTable<TChar, Device<TSpec> >, FibreEntries>
{
    typedef thrust::device_vector<unsigned>         Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreRanks                              [Device RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TwoLevels<TValue, Device<TSpec> > >, FibreRanks>
{
    typedef TwoLevels<TValue, TSpec>        TRankDictionarySpec_;

    typedef thrust::device_vector<RankDictionaryEntry_<TRankDictionarySpec_> >  Type;
};

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<Naive<TValue, Device<TSpec> > >, FibreRanks>
{
    typedef Naive<TValue, TSpec>                                TRankDictionarySpec_;
    typedef RankDictionary<TRankDictionarySpec_>                TRankDictionary_;
    typedef typename Size<TRankDictionary_>::Type               TSize_;

    typedef thrust::device_vector<TSize_>                       Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreSparseString                         [Device CompressedSA]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<CompressedSA<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreSparseString>
{
    typedef thrust::device_vector<TValue/*, TAlloc*/>       TText_;
    typedef typename SAValue<TText_>::Type                  TSAValue_;
    typedef thrust::device_vector<TSAValue_/*, TAlloc*/>    TString_;

    typedef SparseString<TString_, TSpec>                   Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Fibre<CompressedSA<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TSpec>, FibreSparseString>
{
    typedef StringSet<thrust::device_vector<TValue/*, TAlloc*/>, TSSetSpec> TText_;
    typedef typename SAValue<TText_>::Type                                  TSAValue_;
    typedef thrust::device_vector<TSAValue_/*, TAlloc*/>                    TString_;

    typedef SparseString<TString_, TSpec>                                   Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreIndicators                           [Device SparseString]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<SparseString<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreIndicators>
{
    typedef RankDictionary<TwoLevels<bool, Device<TSpec> > >        Type;
};

// ----------------------------------------------------------------------------
// Metafunction FibreValues                               [Device SparseString]
// ----------------------------------------------------------------------------
// NOTE(esiragusa): FibreValues is already defined as below.

//template <typename TValue, typename TAlloc, typename TSpec>
//struct Fibre<SparseString<thrust::device_vector<TValue, TAlloc>, TSpec>, FibreValues>
//{
//    typedef thrust::device_vector<TValue, TAlloc>       Type;
//};

// ----------------------------------------------------------------------------
// Metafunction HistoryStack_                              [Device VSTree Iter]
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TIndexSpec, typename TSpec>
struct HistoryStack_<Iter<Index<thrust::device_vector<TValue, TAlloc>, TIndexSpec>, VSTree<TSpec> > >
{
    typedef thrust::device_vector<TValue, TAlloc>       TText_;
    typedef Index<TText_, TIndexSpec>                   TIndex_;
    typedef Iter<TIndex_, VSTree<TSpec> >               TIter_;
    typedef typename HistoryStackEntry_<TIter_>::Type   TEntry_;
    typedef thrust::device_vector<TEntry_>              Type;
};

template <typename TValue, typename TAlloc, typename TSSetSpec, typename TIndexSpec, typename TSpec>
struct HistoryStack_<Iter<Index<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, TIndexSpec>, VSTree<TSpec> > >
{
    typedef StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>     TText_;
    typedef Index<TText_, TIndexSpec>                                       TIndex_;
    typedef Iter<TIndex_, VSTree<TSpec> >                                   TIter_;
    typedef typename HistoryStackEntry_<TIter_>::Type                       TEntry_;
    typedef thrust::device_vector<TEntry_>                                  Type;
};

// ============================================================================
// Functions
// ============================================================================
// NOTE(esiragusa): Functions assign() are not specific to Device Index.

// ----------------------------------------------------------------------------
// Function assign()                                                  [FMIndex]
// ----------------------------------------------------------------------------

template <typename TText, typename TOccSpec, typename TSpec, typename TText2, typename TOccSpec2, typename TSpec2>
inline void
assign(Index<TText, FMIndex<TOccSpec, TSpec> > & index, Index<TText2, FMIndex<TOccSpec2, TSpec2> > const & source)
{
    assign(indexText(index), indexText(source));
    assign(indexSA(index), indexSA(source));
    assign(indexLF(index), indexLF(source));
}

template <typename TText, typename TOccSpec, typename TSpec, typename TText2, typename TOccSpec2, typename TSpec2>
inline void
assign(Index<TText, FMIndex<TOccSpec, TSpec> > & index, Index<TText2, FMIndex<TOccSpec2, TSpec2> > & source)
{
    assign(index, reinterpret_cast<Index<TText2, FMIndex<TOccSpec2, TSpec2> > const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                                  [LfTable]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TText2, typename TSpec2>
inline void
assign(LfTable<TText, TSpec> & lfTable, LfTable<TText2, TSpec2> const & source)
{
    assign(getFibre(lfTable, FibrePrefixSum()), getFibre(source, FibrePrefixSum()));
    assign(getFibre(lfTable, FibreValues()), getFibre(source, FibreValues()));
    assign(getFibre(lfTable, FibreSentinels()), getFibre(source, FibreSentinels()));
    assign(lfTable.sentinelSubstitute, source.sentinelSubstitute);
}

template <typename TText, typename TSpec, typename TText2, typename TSpec2>
inline void
assign(LfTable<TText, TSpec> & lfTable, LfTable<TText2, TSpec2> & source)
{
    assign(lfTable, reinterpret_cast<LfTable<TText2, TSpec2> const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                           [PrefixSumTable]
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec, typename TChar2, typename TSpec2>
inline void
assign(PrefixSumTable<TChar, TSpec> & pst, PrefixSumTable<TChar2, TSpec2> const & source)
{
    assign(getFibre(pst, FibreEntries()), getFibre(source, FibreEntries()));
}

template <typename TChar, typename TSpec, typename TChar2, typename TSpec2>
inline void
assign(PrefixSumTable<TChar, TSpec> & pst, PrefixSumTable<TChar2, TSpec2> & source)
{
    assign(pst, reinterpret_cast<PrefixSumTable<TChar2, TSpec2> const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline void
assign(RankDictionary<TwoLevels<TValue, TSpec> > & dict, RankDictionary<TwoLevels<TValue2, TSpec2> > const & source)
{
    assign(getFibre(dict, FibreRanks()), getFibre(source, FibreRanks()));
    assign(dict._length, source._length);
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline void
assign(RankDictionary<TwoLevels<TValue, TSpec> > & dict, RankDictionary<TwoLevels<TValue2, TSpec2> > & source)
{
    assign(dict, reinterpret_cast<RankDictionary<TwoLevels<TValue2, TSpec2> > const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                           [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline void
assign(RankDictionary<Naive<TValue, TSpec> > & dict, RankDictionary<Naive<TValue2, TSpec2> > const & source)
{
    assign(getFibre(dict, FibreRanks()), getFibre(source, FibreRanks()));
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline void
assign(RankDictionary<Naive<TValue, TSpec> > & dict, RankDictionary<Naive<TValue2, TSpec2> > & source)
{
    assign(dict, reinterpret_cast<RankDictionary<Naive<TValue2, TSpec2> > const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                             [CompressedSA]
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TText2, typename TSpec2>
inline void
assign(CompressedSA<TText, TSpec> & sa, CompressedSA<TText2, TSpec2> const & source)
{
    assign(getFibre(sa, FibreSparseString()), getFibre(source, FibreSparseString()));
    assign(getFibre(sa, FibreLF()), getFibre(source, FibreLF()));
}

template <typename TText, typename TSpec, typename TText2, typename TSpec2>
inline void
assign(CompressedSA<TText, TSpec> & sa, CompressedSA<TText2, TSpec2> & source)
{
    assign(sa, reinterpret_cast<CompressedSA<TText2, TSpec2> const &>(source));
}

// ----------------------------------------------------------------------------
// Function assign()                                             [SparseString]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TString2, typename TSpec2>
inline void
assign(SparseString<TString, TSpec> & sparseString, SparseString<TString2, TSpec2> const & source)
{
    assign(getFibre(sparseString, FibreValues()), getFibre(source, FibreValues()));
    assign(getFibre(sparseString, FibreIndicators()), getFibre(source, FibreIndicators()));
    assign(sparseString._length, source._length);
}

template <typename TString, typename TSpec, typename TString2, typename TSpec2>
inline void
assign(SparseString<TString, TSpec> & sparseString, SparseString<TString2, TSpec2> & source)
{
    assign(sparseString, reinterpret_cast<SparseString<TString2, TSpec2> const &>(source));
}

}

#endif  // #ifndef SEQAN_EXTRAS_INDEX_DEVICE_H_
