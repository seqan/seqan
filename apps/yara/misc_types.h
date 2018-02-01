// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2018, Enrico Siragusa, FU Berlin
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
// This file contains type overloadings.
// ==========================================================================

#ifndef APP_YARA_MISC_TYPES_H_
#define APP_YARA_MISC_TYPES_H_

using namespace seqan;

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Limits
// ----------------------------------------------------------------------------

template <typename T1 = void, typename T2 = void, typename T3 = void>
struct Limits {};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

namespace seqan
{
template <typename T1, typename T2, typename T3>
struct Value<Limits<T1, T2, T3>, 1>
{
    typedef T1 Type;
};

template <typename T1, typename T2, typename T3>
struct Value<Limits<T1, T2, T3>, 2>
{
    typedef T2 Type;
};

template <typename T1, typename T2, typename T3>
struct Value<Limits<T1, T2, T3>, 3>
{
    typedef T3 Type;
};

// ============================================================================
// StringSet
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction LengthSum
// ----------------------------------------------------------------------------

template <typename TString, typename TSum>
struct LengthSum<StringSet<TString, Owner<ConcatDirect<Limits<TSum> > > > >
{
    typedef TSum Type;
};
}

// ============================================================================
// SeqStore Config
// ============================================================================

// ----------------------------------------------------------------------------
// Reads SeqsStore Config
// ----------------------------------------------------------------------------

typedef SeqConfig<void>         YaraReadsConfig;

// ----------------------------------------------------------------------------
// Contigs SeqsStore Config
// ----------------------------------------------------------------------------

template <typename TSpec = Alloc<> >
struct YaraContigsConfig
{
    typedef Dna5                                    TAlphabet;
    typedef Packed<TSpec>                           TSeqSpec;
    typedef Owner<ConcatDirect<Limits<uint64_t> > > TSeqsSpec;
    typedef Owner<ConcatDirect<Limits<uint32_t> > > TSeqNamesSpec;
};

// ----------------------------------------------------------------------------
// FormattedFile Context
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TString, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Bam, Output, StringSet<TString, Owner<ConcatDirect<TSpec> > > >, TStorageSpec>
{
    typedef StringSet<CharString, Owner<ConcatDirect<TSpec> > >     TNameStore;
    typedef NameStoreCache<TNameStore>                              TNameStoreCache;
    typedef BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> Type;
};

// ----------------------------------------------------------------------------
// Metafunction SeqFileBuffer_
// ----------------------------------------------------------------------------
// Overloaded to load FastQ files with Iupac symbols.

#ifdef YARA_INDEXER
template <typename TString, typename TSSetSpec, typename TSpec>
struct SeqFileBuffer_<StringSet<TString, TSSetSpec>, TSpec>
{
    typedef String<Iupac>   Type;
};
#endif
}

// ============================================================================
// FMIndex
// ============================================================================

// ----------------------------------------------------------------------------
// FMIndex Config
// ----------------------------------------------------------------------------

template <typename TSize, typename TLen, typename TSum, typename TAlloc = Alloc<> >
struct YaraFMConfig
{
    typedef YaraFMConfig<TSize, TLen, TSum, TAlloc>     TMe;

    // Text.
    typedef Owner<ConcatDirect<TMe> >                   TSSetSpec_;
    typedef StringSet<String<Dna, TAlloc>, TSSetSpec_>  Text;

    // LF LengthSum.
    typedef TSum                                        LengthSum;

    // LF's RankDictionary Config.
    typedef Levels<void, TMe>                           Bwt;
    typedef typename If<IsSameType<TSize, uint8_t>,
                        Naive<void, TMe>,
                        Levels<void, TMe> >::Type       Sentinels;

    // RankDictionary Config.
    typedef TAlloc                                      Fibre;
    typedef TSum                                        Size;

    // Sparse SA sampling rate.
    static const unsigned SAMPLING =                    10;
    static const unsigned WORDS_PER_BLOCK =             0;
    static const unsigned LEVELS =                      1;
};

// ----------------------------------------------------------------------------
// FMIndex SAValue
// ----------------------------------------------------------------------------

namespace seqan {
template <typename TValue, typename TSpec, typename TSize, typename TLen, typename TSum, typename TAlloc>
struct SAValue<StringSet<String<TValue, TSpec>, Owner<ConcatDirect<YaraFMConfig<TSize, TLen, TSum, TAlloc> > > > >
{
    typedef Pair<TSize, TLen, Pack>   Type;
};

// ----------------------------------------------------------------------------
// SparseString FibreIndicators
// ----------------------------------------------------------------------------
// NOTE(esiragusa): CompressedSA propagates down TConfig as TSpec.
// TODO(esiragusa): remove this crap once the CSA gets refactored.

template <typename TValue, typename TAlloc, typename TSpec>
struct Fibre<SparseString<String<TValue, TAlloc>, TSpec>, FibreIndicators>
{
    typedef RankDictionary<bool, Levels<void, TSpec> > Type;
};

// ----------------------------------------------------------------------------
// SparseString Size
// ----------------------------------------------------------------------------
// TODO(esiragusa): remove this crap once the CSA gets refactored.

template <typename TString, typename TSize, typename TLen, typename TSum, typename TAlloc>
struct Size<SparseString<TString, YaraFMConfig<TSize, TLen, TSum, TAlloc> > >
{
    typedef TSum Type;
};
}

#endif  // #ifndef APP_YARA_MISC_TYPES_H_
