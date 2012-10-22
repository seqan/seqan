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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains basic type definitions.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_STORE_H_
#define SEQAN_EXTRAS_MASAI_STORE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include "store/store_io.h"

using namespace seqan;

// ============================================================================
// Fragment Store Types
// ============================================================================

// ----------------------------------------------------------------------------
// Fragment Store Configuration
// ----------------------------------------------------------------------------

struct FragStoreConfig
{
    typedef String<Dna5Q>           TReadSeq;
    typedef String<Dna5>            TContigSeq;

    typedef double                  TMean;
    typedef double                  TStd;
    typedef signed char             TMappingQuality;

    typedef void                    TReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >  TReadSeqStoreSpec;
    typedef MMap<>                  TReadNameSpec;
    typedef Owner<ConcatDirect<> >  TReadNameStoreSpec;
    typedef void                    TMatePairStoreElementSpec;
    typedef void                    TLibraryStoreElementSpec;
    typedef void                    TContigStoreElementSpec;
    typedef void                    TContigFileSpec;
    typedef void                    TAlignedReadStoreElementSpec;
    typedef Owner<ConcatDirect<> >  TAlignedReadTagStoreSpec;
    typedef void                    TAnnotationStoreElementSpec;
};

// ----------------------------------------------------------------------------
// Genome Type
// ----------------------------------------------------------------------------

typedef StringSet<FragStoreConfig::TContigSeq, Dependent<> >                TGenome;

// TODO(esiragusa):Overload Size/Limits metafunctions and test
//namespace seqan
//{
//    template <>
//    struct Size<FragStoreConfig::TContigSeq>
//    {
//        typedef unsigned int            Type;
//    };
//
//    template <>
//    struct StringSetLimits<TGenome>
//    {
//        typedef String<unsigned char>   Type;
//    };
//}

// ----------------------------------------------------------------------------
// Fragment Store Type
// ----------------------------------------------------------------------------

typedef FragmentStore<void, FragStoreConfig>            TFragmentStore;

// ----------------------------------------------------------------------------
// Fragment Store Genome Types
// ----------------------------------------------------------------------------

typedef TFragmentStore::TContigStore                    TContigStore;
typedef Size<TContigStore>::Type                        TContigStoreSize;
typedef Value<TContigStore>::Type                       TContigStoreElement;
typedef TFragmentStore::TContigSeq                      TContigSeq;
typedef Size<TContigSeq>::Type                          TContigSeqSize;
typedef Segment<TContigSeq, InfixSegment>               TContigInfix;

// ----------------------------------------------------------------------------
// Fragment Store Reads Types
// ----------------------------------------------------------------------------

typedef TFragmentStore::TReadSeqStore                   TReadSeqStore;
typedef Size<TReadSeqStore>::Type                       TReadSeqStoreSize;
typedef Value<TReadSeqStore>::Type const                TReadSeq;
typedef Size<TReadSeq>::Type                            TReadSeqSize;

// ----------------------------------------------------------------------------
// Fragment Store Mapped Reads Types
// ----------------------------------------------------------------------------

typedef TFragmentStore::TAlignedReadStore               TAlignedReadStore;
typedef Value<TAlignedReadStore>::Type                  TAlignedReadStoreElement;
typedef TFragmentStore::TAlignQualityStore              TAlignQualityStore;
typedef Value<TAlignQualityStore>::Type                 TAlignQualityStoreElement;
typedef TFragmentStore::TAlignedReadTagStore            TAlignedReadTagStore;
typedef Value<TAlignedReadTagStore>::Type               TAlignedReadTagStoreElement;


// ============================================================================
// Dna5 specializations to deal with uncalled bases
// ============================================================================

namespace seqan {
static unsigned char __MASK_DNA5_EQ[]  = {1, 2, 4, 8, 0};
static unsigned char __MASK_DNA5_LT[]  = {0, 1, 2, 3, 4};
static unsigned char __MASK_DNA5Q_LT[] = {0, 1, 2, 3, 5};

// ----------------------------------------------------------------------------
// Comparison operators for (Dna5 vs Dna5Q) and (Dna5Q vs Dna5)
// ----------------------------------------------------------------------------

template <>
inline bool operator==(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)];
}

template <>
inline bool operator==(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)];
}

template <>
inline bool operator!=(Dna5 const & left_, Dna5Q const & right_)
{
    return !(__MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)]);
}

template <>
inline bool operator!=(Dna5Q const & left_, Dna5 const & right_)
{
    return !(__MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)]);
}

template <>
inline bool operator<(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] < __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool operator<(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5Q_LT[ordValue(left_)] < __MASK_DNA5_LT[ordValue(right_)];
}

template <>
inline bool operator>(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] > __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool operator>(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5Q_LT[ordValue(left_)] > __MASK_DNA5_LT[ordValue(right_)];
}

template <>
inline bool operator<=(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] <= __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool operator<=(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5Q_LT[ordValue(left_)] <= __MASK_DNA5_LT[ordValue(right_)];
}

template <>
inline bool operator>=(Dna5 const & left_, Dna5Q const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] >= __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool operator>=(Dna5Q const & left_, Dna5 const & right_)
{
    return __MASK_DNA5Q_LT[ordValue(left_)] >= __MASK_DNA5_LT[ordValue(right_)];
}

// ----------------------------------------------------------------------------
// ordLess/Equal/Greater() functions for (Dna5 vs Dna5)
// ----------------------------------------------------------------------------
    
template <>
inline bool ordLess(Dna5 const & left_, Dna5 const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] < __MASK_DNA5Q_LT[ordValue(right_)];
}

template <>
inline bool ordEqual(Dna5 const & left_, Dna5 const & right_)
{
    return __MASK_DNA5_EQ[ordValue(left_)] & __MASK_DNA5_EQ[ordValue(right_)];
}

template <>
inline bool ordGreater(Dna5 const & left_, Dna5 const & right_)
{
    return __MASK_DNA5_LT[ordValue(left_)] > __MASK_DNA5Q_LT[ordValue(right_)];
}
}


#endif  // #ifndef SEQAN_EXTRAS_MASAI_STORE_H_
