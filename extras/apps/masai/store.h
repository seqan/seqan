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
#include <seqan/seq_io.h>
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
    typedef Alloc<>                 TReadNameSpec;
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

namespace seqan
{
    template <>
    struct Size<FragStoreConfig::TContigSeq>
    {
        typedef unsigned int            Type;
    };

    // NOTE(esiragusa): Genome can be at most 2^32 bp in total
    template <>
    struct StringSetLimits<TGenome>
    {
        typedef String<unsigned int>    Type;
    };
}

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

typedef TFragmentStore::TReadStore                      TReadStore;
typedef Value<TReadStore>::Type                         TReadStoreElement;
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

template <typename TReadSeq>
void _storeReadSeq(TFragmentStore & store, TReadSeq const & seq)
{
    appendValue(store.readSeqStore, seq, Generous());
}

template <typename TSize>
void _reserveReadStore(TFragmentStore & store, TSize space, True const & /* tag */)
{
    reserve(store.readStore, space, Exact());
}

template <typename TSize>
void _reserveReadStore(TFragmentStore &, TSize, False const & /* tag */)
{}

// ============================================================================

template <typename TSize, typename TLength>
void _reserveReadNameStore(TFragmentStore & store, TSize count, TLength length, True const & /* tag */)
{
    reserve(store.readNameStore.concat, count * length, Exact());
    reserve(store.readNameStore, count, Exact());
}

template <typename TSize, typename TLength>
void _reserveReadNameStore(TFragmentStore &, TSize, TLength, False const & /* tag */)
{}

// ============================================================================

template <typename TReadSeqName>
void _storeReadName(TFragmentStore & store, TReadSeqName const & seqName, True const & /* tag */)
{
    appendValue(store.readNameStore, seqName, Generous());
}

template <typename TReadSeqName>
void _storeReadName(TFragmentStore &, TReadSeqName const &, False const & /* tag */)
{}

// ============================================================================

template <typename TReadId>
void _storeReadId(TFragmentStore & store, TReadId const & matePairId, True const & /* tag */)
{
	typename Value<TFragmentStore::TReadStore>::Type r;
	r.matePairId = matePairId;

	appendValue(store.readStore, r, Generous());
}

template <typename TReadId>
void _storeReadId(TFragmentStore &, TReadId const &, False const & /* tag */)
{}

// ============================================================================

template <typename TSeqName, typename TSeq>
unsigned long _estimateRecordLength(TSeqName const & seqName, TSeq const & seq, Fastq const & /* tag */)
{
    // 6 stands for: @, +, and four \n.
    return length(seqName) + 2 * length(seq) + 6;
}

template <typename TSeqName, typename TSeq>
unsigned long _estimateRecordLength(TSeqName const & seqName, TSeq const & seq, Fasta const & /* tag */)
{
    // 3 stands for: >, and two \n.
    return length(seqName) + length(seq) + 3;
}

// ============================================================================

template <typename TRecordReader, typename TFormat, typename TUseReadStore, typename TUseReadNameStore>
bool _loadReads(TFragmentStore & store,
                TRecordReader & reader,
                unsigned long fileSize,
                TFormat const & /* tag */,
                TUseReadStore const & /* tag */,
                TUseReadNameStore const & /* tag */)
{
    CharString seqName;
    FragStoreConfig::TReadSeq seq;

    // Read first record.
    if (readRecord(seqName, seq, reader, TFormat()) != 0)
        return false;

    // Estimate record size.
    unsigned long recordLength = _estimateRecordLength(seqName, seq, TFormat());

    // Estimate number of records.
    unsigned long numberOfRecords = fileSize / recordLength;

    // Reserve space in the readSeqStore, also considering reverse complemented reads.
    reserve(store.readSeqStore.concat, 2 * numberOfRecords * length(seq), Exact());
    reserve(store.readSeqStore, 2 * numberOfRecords, Exact());

    // Reserve space in the readStore.
    _reserveReadStore(store, numberOfRecords, TUseReadStore());

    // Reserve space in the readNameStore.
    _reserveReadNameStore(store, numberOfRecords, length(seqName), TUseReadNameStore());

    // Store first record.
    _storeReadSeq(store, seq);
    _storeReadName(store, seqName, TUseReadNameStore());
    _storeReadId(store, TReadStoreElement::INVALID_ID, TUseReadStore());

    // Read whole file.
    while (!atEnd(reader))
    {
        if (readRecord(seqName, seq, reader, TFormat()) != 0)
            return false;

        _storeReadSeq(store, seq);
        _storeReadName(store, seqName, TUseReadNameStore());
        _storeReadId(store, TReadStoreElement::INVALID_ID, TUseReadStore());
    }

    return true;
}

// ----------------------------------------------------------------------------
// Function loadReads()                                         [FragmentStore]
// ----------------------------------------------------------------------------

// TODO(esiragusa): Implement paired-end loadReads()
template <typename TFileName, typename TUseReadStore, typename TUseReadNameStore>
bool loadReads(TFragmentStore & store,
               TFileName & fileName,
               TUseReadStore const & /* tag */,
               TUseReadNameStore const & /* tag */)
{
    typedef std::fstream                            TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;

    TStream file(toCString(fileName), std::ios::binary | std::ios::in);

    // Compute file size.
    file.seekg(0, std::ios::end);
    unsigned long fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    TRecordReader reader(file);

    // Autodetect file format.
    AutoSeqStreamFormat tagSelector;
    checkStreamFormat(reader, tagSelector);

    switch (tagSelector.tagId)
    {
    case 1:
        return _loadReads(store, reader, fileSize, Fasta(), TUseReadStore(), TUseReadNameStore());

    case 2:
        return _loadReads(store, reader, fileSize, Fastq(), TUseReadStore(), TUseReadNameStore());

    default:
        return false;
    }
}


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
