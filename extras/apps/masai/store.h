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
#include "tags.h"

using namespace seqan;

// ============================================================================
// Types
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
// Contigs Type
// ----------------------------------------------------------------------------

typedef StringSet<FragStoreConfig::TContigSeq, Dependent<> >                TContigs;

namespace seqan
{
    template <>
    struct Size<FragStoreConfig::TContigSeq>
    {
        typedef unsigned int            Type;
    };

    // NOTE(esiragusa): Genome can be at most 2^32 bp in total
    template <>
    struct StringSetLimits<TContigs>
    {
        typedef String<unsigned int>    Type;
    };
}

// ----------------------------------------------------------------------------
// Fragment Store Type
// ----------------------------------------------------------------------------

typedef FragmentStore<void, FragStoreConfig>            TFragmentStore;

// ----------------------------------------------------------------------------
// Fragment Store Contig Types
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
typedef TFragmentStore::TReadNameStore                  TReadNameStore;
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
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Genome
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Genome
{
    TFragmentStore          & _store;
    TContigs                contigs;

    Genome(TFragmentStore & store) :
        _store(store)
    {}
};

// ----------------------------------------------------------------------------
// Class ReadsConfig
// ----------------------------------------------------------------------------

template <typename TUseReadStore_       = True,
          typename TUseReadNameStore_   = True,
          typename TForward_            = True,
          typename TReverse_            = True>
struct ReadsConfig
{
    typedef TUseReadStore_      TUseReadStore;
    typedef TUseReadNameStore_  TUseReadNameStore;
    typedef TForward_           TForward;
    typedef TReverse_           TReverse;
};

// ----------------------------------------------------------------------------
// Class Reads
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = ReadsConfig<> >
struct Reads
{
    Holder<TFragmentStore>          _store;
    unsigned                        _avgSeqLengthEstimate;
    unsigned                        _avgNameLengthEstimate;
    unsigned                        _countEstimate;
    unsigned                        readsCount;

    Reads() :
        _store(),
        _avgSeqLengthEstimate(0),
        _avgNameLengthEstimate(0),
        _countEstimate(0),
        readsCount(0)
    {}

    Reads(TFragmentStore & store) :
        _store(store),
        _avgSeqLengthEstimate(0),
        _avgNameLengthEstimate(0),
        _countEstimate(0),
        readsCount(0)
    {}
};

// ----------------------------------------------------------------------------
// Class ReadsLoader
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = ReadsConfig<> >
struct ReadsLoader
{
    // TODO(esiragusa): Use MMap FileReader. Support compressed streams.
    typedef std::fstream                            TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;
    typedef Reads<TSpec, TConfig>                   TReads;

    TStream                         _file;
    unsigned long                   _fileSize;
    AutoSeqStreamFormat             _fileFormat;
    std::auto_ptr<TRecordReader>    _reader;
    Holder<TReads>                  reads;

    ReadsLoader(TReads & reads) :
        _fileSize(0),
        reads(reads)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function load()                                                     [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TString>
bool load(Genome<TSpec> & genome, TString const & genomeFile)
{
    // TODO(esiragusa): Use record reader instead of loadContigs() from FragmentStore.
    if (!loadContigs(genome._store, genomeFile))
        return false;

    reserve(genome.contigs, length(genome._store.contigStore));
    for (unsigned contigId = 0; contigId < length(genome._store.contigStore); ++contigId)
    {
        shrinkToFit(genome._store.contigStore[contigId].seq);
        appendValue(genome.contigs, genome._store.contigStore[contigId].seq);
    }

    return true;
}

// ----------------------------------------------------------------------------
// Function reverse()                                                  [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec>
void reverse(Genome<TSpec> & genome)
{
    for (unsigned contigId = 0; contigId < length(genome._store.contigStore); ++contigId)
        reverse(genome._store.contigStore[contigId].seq);
}

// ----------------------------------------------------------------------------
// Function clear()                                                    [Genome]
// ----------------------------------------------------------------------------

template <typename TSpec>
void clear(Genome<TSpec> & genome)
{
    clear(genome.contigs);
    clearContigs(genome._store);
}

// ----------------------------------------------------------------------------
// Function reserve()                                                   [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void reserve(Reads<TSpec, TConfig> & reads)
{
    reserve(reads, reads._countEstimate);
}

template <typename TConfig>
inline void reserve(Reads<PairedEnd, TConfig> & /* reads */)
{
}

template <typename TSpec, typename TConfig, typename TSize>
inline void reserve(Reads<TSpec, TConfig> & reads, TSize count)
{
    // Reserve space in the readSeqStore, also considering reverse complemented reads.
    reserve(value(reads._store).readSeqStore.concat, 2 * count * reads._avgSeqLengthEstimate, Exact());
    reserve(value(reads._store).readSeqStore, 2 * count, Exact());

    // Reserve space for ids.
    reserveIds(reads, count, typename TConfig::TUseReadStore());

    // Reserve space for names.
    reserveNames(reads, count, reads._avgNameLengthEstimate, typename TConfig::TUseReadNameStore());
}

// ----------------------------------------------------------------------------
// Function reserveIds()                                                [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
inline void reserveIds(Reads<TSpec, TConfig> & /* reads */, TSize /* space */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TSize>
inline void reserveIds(Reads<TSpec, TConfig> & reads, TSize space, True const & /* tag */)
{
    reserve(getSeqs(reads), space, Exact());
}

// ----------------------------------------------------------------------------
// Function reserveNames()                                              [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize, typename TLength>
inline void reserveNames(Reads<TSpec, TConfig> & /* reads */, TSize /* count */, TLength /* length */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TSize, typename TLength>
inline void reserveNames(Reads<TSpec, TConfig> & reads, TSize count, TLength length, True const & /* tag */)
{
    reserve(getNames(reads).concat, count * length, Exact());
    reserve(getNames(reads), count, Exact());
}

// ----------------------------------------------------------------------------
// Function clear()                                                     [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void clear(Reads<TSpec, TConfig> & reads)
{
    clearReads(value(reads._store));
}

// ----------------------------------------------------------------------------
// Function appendSeq()                                                 [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadSeq>
inline void appendSeq(Reads<TSpec, TConfig> & reads, TReadSeq const & seq)
{
    appendValue(getSeqs(reads), seq, Generous());
}

// ----------------------------------------------------------------------------
// Function appendName()                                                [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadName>
inline void appendName(Reads<TSpec, TConfig> & /* reads */, TReadName const & /* seqName */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TReadName>
inline void appendName(Reads<TSpec, TConfig> & reads, TReadName const & seqName, True const & /* tag */)
{
    appendValue(getNames(reads), seqName, Generous());
}

// ----------------------------------------------------------------------------
// Function appendId()                                                  [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadId>
inline void appendId(Reads<TSpec, TConfig> & /* reads */, TReadId const & /* matePairId */, False const & /* tag */)
{}

template <typename TSpec, typename TConfig, typename TReadId>
inline void appendId(Reads<TSpec, TConfig> & reads, TReadId const & matePairId, True const & /* tag */)
{
	typename Value<TFragmentStore::TReadStore>::Type r;
	r.matePairId = matePairId;

	appendValue(getIds(reads), r, Generous());
}

// ----------------------------------------------------------------------------
// Function getSeqs()                                                   [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline TReadSeqStore &
getSeqs(Reads<TSpec, TConfig> const & reads)
{
    return value(reads._store).readSeqStore;
}

// ----------------------------------------------------------------------------
// Function getNames()                                                  [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline TReadNameStore &
getNames(Reads<TSpec, TConfig> const & reads)
{
    return value(reads._store).readNameStore;
}

// ----------------------------------------------------------------------------
// Function getIds()                                                    [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline TReadStore &
getIds(Reads<TSpec, TConfig> const & reads)
{
    return value(reads._store).readStore;
}

// ----------------------------------------------------------------------------
// Function getReadId()                                                 [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TReadId>
inline TReadId getReadId(Reads<TSpec, TConfig> const & reads, TReadId readId)
{
    // Deal with reverse complemented reads.
    if (readId >= reads.readsCount)
        return readId - reads.readsCount;

    return readId;
}

// ----------------------------------------------------------------------------
// Function avgSeqLength()                                              [Reads]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
TReadSeqSize avgSeqLength(Reads<TSpec, TConfig> const & reads)
{
    return reads._avgSeqLengthEstimate;
}

// ----------------------------------------------------------------------------
// Function _estimateRecordSize()                                 [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
unsigned long _estimateRecordSize(ReadsLoader<TSpec, TConfig> const & loader, Fastq const & /* tag */)
{
    // 6 stands for: @, +, and four \n.
    return value(loader.reads)._avgNameLengthEstimate + 2 * value(loader.reads)._avgSeqLengthEstimate + 6;
}

template <typename TSpec, typename TConfig>
unsigned long _estimateRecordSize(ReadsLoader<TSpec, TConfig> const & loader, Fasta const & /* tag */)
{
    // 3 stands for: >, and two \n.
    return value(loader.reads)._avgNameLengthEstimate + value(loader.reads)._avgSeqLengthEstimate + 3;
}

// ----------------------------------------------------------------------------
// Function _estimateReadsStatistics()                            [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void _estimateReadsStatistics(ReadsLoader<TSpec, TConfig> & loader)
{
    CharString                  seqName;
    FragStoreConfig::TReadSeq   seq;

    // Read first record.
    if (readRecord(seqName, seq, *(loader._reader), loader._fileFormat) != 0)
        return;

    // Estimate read seqs and names length.
    value(loader.reads)._avgSeqLengthEstimate = length(seq);
    value(loader.reads)._avgNameLengthEstimate = length(seqName);

    // Estimate record size.
    unsigned long recordSize;
    switch (loader._fileFormat.tagId)
    {
    case Find<AutoSeqStreamFormat, Fasta>::VALUE:
        recordSize = _estimateRecordSize(loader, Fasta());
        break;
    case Find<AutoSeqStreamFormat, Fastq>::VALUE:
        recordSize = _estimateRecordSize(loader, Fastq());
        break;
    default:
        recordSize = 0;
        break;
    }

    // Estimate number of reads in file.
    if (recordSize > 0)
        value(loader.reads)._countEstimate = loader._fileSize / recordSize;
    else
        value(loader.reads)._countEstimate = 0;
}

// ----------------------------------------------------------------------------
// Function open()                                                [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TString>
bool open(ReadsLoader<TSpec, TConfig> & loader, TString const & readsFile)
{
    typedef ReadsLoader<TSpec, TConfig>             TReadsLoader;
    typedef typename TReadsLoader::TRecordReader    TRecordReader;

    // Open file.
    loader._file.open(toCString(readsFile), std::ios::binary | std::ios::in);

    if (!loader._file.is_open())
        return false;

    // Compute file size.
    loader._file.seekg(0, std::ios::end);
    loader._fileSize = loader._file.tellg();
    loader._file.seekg(0, std::ios::beg);

    // Initialize record reader.
    loader._reader.reset(new TRecordReader(loader._file));

    // Autodetect file format.
    if (!guessStreamFormat(*(loader._reader), loader._fileFormat))
        return false;

    // Estimate statistics for reads in file.
    _estimateReadsStatistics(loader);

    // Reinitialize record reader.
    loader._file.seekg(0, std::ios::beg);
    loader._reader.reset(new TRecordReader(loader._file));

    return true;
}

template <typename TConfig, typename TString>
bool open(ReadsLoader<PairedEnd, TConfig> & loader, TString const & readsLeftFile, TString const & readsRightFile)
{
    if (!loadReads(value(value(loader.reads)._store), readsLeftFile, readsRightFile))
        return false;

    value(loader.reads).readsCount = length(getSeqs(value(loader.reads)));

    _loadReverseComplement(loader);

    return true;
}

// ----------------------------------------------------------------------------
// Function close()                                               [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
bool close(ReadsLoader<TSpec, TConfig> & loader)
{
    loader._file.close();

    return !loader._file.is_open();
}

template <typename TConfig>
bool close(ReadsLoader<PairedEnd, TConfig> & /* loader */)
{
    return true;
}

// ----------------------------------------------------------------------------
// Function load()                                                [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
bool load(ReadsLoader<TSpec, TConfig> & loader)
{
    return load(loader, MaxValue<unsigned long>::VALUE);
}

template <typename TConfig>
bool load(ReadsLoader<PairedEnd, TConfig> & /* loader */)
{
    return true;
}

template <typename TSpec, typename TConfig, typename TSize>
bool load(ReadsLoader<TSpec, TConfig> & loader, TSize count)
{
    switch (loader._fileFormat.tagId)
    {
    case Find<AutoSeqStreamFormat, Fasta>::VALUE:
        return load(loader, count, Fasta());
    case Find<AutoSeqStreamFormat, Fastq>::VALUE:
        return load(loader, count, Fastq());
    default:
        return false;
    }
}

template <typename TSpec, typename TConfig, typename TSize, typename TFormat>
bool load(ReadsLoader<TSpec, TConfig> & loader, TSize count, TFormat const & /* tag */)
{
    CharString                  seqName;
    FragStoreConfig::TReadSeq   seq;

    // Read records.
    for (; count > 0 && !atEnd(loader); count--)
    {
        // NOTE(esiragusa): AutoFormat seems to thrash memory allocation.
//        if (readRecord(seqName, seq, *(loader._reader), loader._fileFormat) != 0)

        if (readRecord(seqName, seq, *(loader._reader), TFormat()) != 0)
            return false;

        appendSeq(value(loader.reads), seq);
        appendName(value(loader.reads), seqName, typename TConfig::TUseReadNameStore());
        appendId(value(loader.reads), TReadStoreElement::INVALID_ID, typename TConfig::TUseReadStore());
    }

    value(loader.reads).readsCount = length(getSeqs(value(loader.reads)));

    // Append reverse complemented reads.
    _loadReverseComplement(loader);

    return true;
}

// ----------------------------------------------------------------------------
// Function _loadReverseComplement()                              [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void _loadReverseComplement(ReadsLoader<TSpec, TConfig> & loader)
{
    for (TReadSeqStoreSize readId = 0; readId < value(loader.reads).readsCount; ++readId)
    {
        TReadSeq & read = getSeqs(value(loader.reads))[readId];
        appendSeq(value(loader.reads), read);
        reverseComplement(back(getSeqs(value(loader.reads))));
    }
}

// ----------------------------------------------------------------------------
// Function atEnd()                                               [ReadsLoader]
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline bool atEnd(ReadsLoader<TSpec, TConfig> & reads)
{
    return atEnd(*(reads._reader));
}

// ============================================================================
// Operators
// ============================================================================

// Dna5 specializations to deal with uncalled bases

namespace seqan {
static unsigned char __MASK_DNA5_EQ[]  = {1, 2, 4, 8, 0};
static unsigned char __MASK_DNA5_LT[]  = {0, 1, 2, 3, 4};
static unsigned char __MASK_DNA5Q_LT[] = {0, 1, 2, 3, 5};

// ----------------------------------------------------------------------------
// Operators ==, !=, <, >, <=, >=                [Dna5 vs Dna5Q, Dna5Q vs Dna5]
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
// Functions ordLess/Equal/Greater()                             [Dna5 vs Dna5]
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
