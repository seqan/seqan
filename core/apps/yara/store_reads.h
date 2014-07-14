// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
// This file contains the Reads and ReadsLoader classes.
// ==========================================================================
// NOTE(esiragusa): the FragmentStore should provide these functionalities.

#ifndef APP_YARA_STORE_READS_H_
#define APP_YARA_STORE_READS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ReadsConfig
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct ReadsConfig
{
    typedef String<Dna5Q>           TReadSeq;
    typedef Owner<ConcatDirect<> >  TReadSpec;
    typedef Owner<ConcatDirect<> >  TReadNameSpec;
    typedef Nothing                 TInputType;
};

// ----------------------------------------------------------------------------
// Class Reads
// ----------------------------------------------------------------------------

template <typename TSpec = SingleEnd, typename TConfig = ReadsConfig<TSpec> >
struct Reads
{
    typedef typename TConfig::TReadSeq                  TReadSeq;
    typedef typename TConfig::TReadSpec                 TReadSpec;
    typedef typename TConfig::TReadNameSpec             TReadNameSpec;

    typedef StringSet<TReadSeq, TReadSpec>              TReadSeqs;
	typedef StringSet<CharString, TReadNameSpec>        TReadNames;
	typedef NameStoreCache<TReadNames, CharString>      TReadNamesCache;

    TReadSeqs           seqs;
    TReadNames          names;
    TReadNamesCache     namesCache;

    Reads() :
        seqs(),
        names(),
		namesCache(names)
    {}
};

// ----------------------------------------------------------------------------
// Metafunction InputStream
// ----------------------------------------------------------------------------

template <typename TSpec>
struct InputStream
{
    typedef std::fstream    Type;
};

template <>
struct InputStream<GZFile>
{
    typedef Stream<GZFile>  Type;
};

#if SEQAN_HAS_BZIP2
template <>
struct InputStream<BZ2File>
{
    typedef Stream<BZ2File> Type;
};
#endif

// ----------------------------------------------------------------------------
// Class ReadsLoader
// ----------------------------------------------------------------------------

template <typename TSpec = SingleEnd, typename TConfig = ReadsConfig<> >
struct ReadsLoader
{
    typedef typename TConfig::TInputType            TInputType;
    typedef typename InputStream<TInputType>::Type  TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;

    TStream                         _file;
    AutoSeqStreamFormat             _fileFormat;
    std::auto_ptr<TRecordReader>    _reader;
};

// ----------------------------------------------------------------------------
// Class ReadsLoader; PairedEnd
// ----------------------------------------------------------------------------

template <typename TConfig>
struct ReadsLoader<PairedEnd, TConfig>
{
    typedef typename TConfig::TInputType            TInputType;
    typedef typename InputStream<TInputType>::Type  TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;

    TStream                             _file1;
    TStream                             _file2;
    Pair<AutoSeqStreamFormat>           _fileFormat;
    Pair<std::auto_ptr<TRecordReader> > _reader;
};

// ----------------------------------------------------------------------------
// Class LoadReadsWorker
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
struct LoadReadsWorker
{
    Reads<TSpec, TConfig> *         reads;
    ReadsLoader<TSpec, TConfig> &   readsLoader;
    unsigned                        readsCount;

    LoadReadsWorker(Reads<TSpec, TConfig> * reads, ReadsLoader<TSpec, TConfig> & readsLoader, unsigned readsCount) :
        reads(reads),
        readsLoader(readsLoader),
        readsCount(readsCount)
    {}

    void operator() ()
    {
        load(*reads, readsLoader, readsCount);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------
// TODO(esiragusa): remove this sh*t when the new i/o gets merged.

namespace seqan {
template <typename TSpec>
inline bool open(Stream<TSpec> & stream, const char *fileName, int /* mode */)
{
    return open(stream, fileName, "r");
}
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TString>
void open(ReadsLoader<TSpec, TConfig> & me, TString const & readsFile)
{
    typedef ReadsLoader<TSpec, TConfig>             TReadsLoader;
    typedef typename TReadsLoader::TRecordReader    TRecordReader;

    // Open file.
    if (!open(me._file, toCString(readsFile), OPEN_RDONLY))
        throw RuntimeError("Error while opening reads file.");

    // Initialize record reader.
    me._reader.reset(new TRecordReader(me._file));

    // Autodetect file format.
    if (!guessStreamFormat(*(me._reader), me._fileFormat))
        throw RuntimeError("Error while guessing reads file format.");
}

template <typename TConfig, typename TString>
void open(ReadsLoader<PairedEnd, TConfig> & me, Pair<TString> const & readsFile)
{
    typedef ReadsLoader<PairedEnd, TConfig>         TReadsLoader;
    typedef typename TReadsLoader::TRecordReader    TRecordReader;

    // Open files.
    if (!open(me._file1, toCString(readsFile.i1), OPEN_RDONLY))
        throw RuntimeError("Error while opening reads file.");

    if (!open(me._file2, toCString(readsFile.i2), OPEN_RDONLY))
        throw RuntimeError("Error while opening reads file.");

    // Initialize record reader.
    me._reader.i1.reset(new TRecordReader(me._file1));
    me._reader.i2.reset(new TRecordReader(me._file2));

    // Autodetect file format.
    if (!guessStreamFormat(*(me._reader.i1), me._fileFormat.i1) ||
        !guessStreamFormat(*(me._reader.i2), me._fileFormat.i2))
        throw RuntimeError("Error while guessing reads file format.");
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void close(ReadsLoader<TSpec, TConfig> & me)
{
    close(me._file);
}

template <typename TConfig>
void close(ReadsLoader<PairedEnd, TConfig> & me)
{
    close(me._file1);
    close(me._file2);
}

// ----------------------------------------------------------------------------
// Function load()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
void load(Reads<TSpec, TConfig> & reads, ReadsLoader<TSpec, TConfig> & me, TSize count)
{
    _load(reads, count, *(me._reader), me._fileFormat);
}

template <typename TConfig, typename TSize>
void load(Reads<PairedEnd, TConfig> & reads, ReadsLoader<PairedEnd, TConfig> & me, TSize count)
{
    _load(reads, count, *(me._reader.i1), me._fileFormat.i1);
    _load(reads, count, *(me._reader.i2), me._fileFormat.i2);
}

template <typename TSpec, typename TConfig, typename TSize, typename TReader, typename TFormat>
void _load(Reads<TSpec, TConfig> & reads, TSize count, TReader & reader, TFormat & format)
{
    typedef Reads<TSpec, TConfig>           TReads;
    typedef typename TReads::TReadSeq       TReadSeq;

    CharString  seqName;
    TReadSeq    seq;

    // Read records.
    for (; count > 0 && !atEnd(reader); count--)
    {
        if (readRecord(seqName, seq, reader, format) != 0)
            throw RuntimeError("Error while reading read record.");

        appendValue(reads.seqs, seq, Generous());
        appendValue(reads.names, prefix(seqName, lastOf(seqName, IsSpace())), Generous());
    }
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline bool atEnd(ReadsLoader<TSpec, TConfig> & reads)
{
    return atEnd(*(reads._reader));
}

template <typename TConfig>
inline bool atEnd(ReadsLoader<PairedEnd, TConfig> & reads)
{
    return atEnd(*(reads._reader.i1)) && atEnd(*(reads._reader.i2));
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
void clear(Reads<TSpec, TConfig> & me)
{
    clear(me.seqs);
    clear(me.names);
//    clear(me.namesCache);
}

// ----------------------------------------------------------------------------
// Function appendReverseComplement()
// ----------------------------------------------------------------------------
// Append reverse complemented reads.

template <typename TSpec, typename TConfig>
void appendReverseComplement(Reads<TSpec, TConfig> & me)
{
    typedef Reads<TSpec, TConfig>           TReads;
    typedef typename TReads::TReadSeqs      TReadSeqs;
    typedef typename Value<TReadSeqs>::Type TReadSeq;
    typedef typename Size<TReadSeqs>::Type  TReadSeqId;

    TReadSeqId readSeqsCount = length(me.seqs);

    reserve(me.seqs, 2 * readSeqsCount, Exact());
    reserve(concat(me.seqs), 2 * lengthSum(me.seqs), Exact());

    for (TReadSeqId readSeqId = 0; readSeqId < readSeqsCount; ++readSeqId)
    {
        TReadSeq const & read = me.seqs[readSeqId];
        appendValue(me.seqs, read);
        reverseComplement(back(me.seqs));
    }
}

// ----------------------------------------------------------------------------
// Functions on ReadSeqsStore
// ----------------------------------------------------------------------------

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getReadSeqsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs);
}

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getReadsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs) / 2;
}

template <typename TReadSeqs>
inline typename Size<TReadSeqs>::Type
getPairsCount(TReadSeqs const & readSeqs)
{
    return length(readSeqs) / 4;
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isFwdReadSeq(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return readSeqId < getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isRevReadSeq(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return !isFwdReadSeq(readSeqs, readSeqId);
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isFirstMate(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return isFwdReadSeq(readSeqs, readSeqId) ? readSeqId < getPairsCount(readSeqs) : readSeqId < getPairsCount(readSeqs) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline bool isSecondMate(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));
    return !isFirstMate(readSeqs, readSeqId);
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getFirstMateFwdSeqId(TReadSeqs const & /* readSeqs */, TPairId pairId)
{
    return pairId;
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getSecondMateFwdSeqId(TReadSeqs const & readSeqs, TPairId pairId)
{
    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    return pairId + getPairsCount(readSeqs);
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getFirstMateRevSeqId(TReadSeqs const & readSeqs, TPairId pairId)
{
//    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    return getFirstMateFwdSeqId(readSeqs, pairId) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TPairId>
inline typename Size<TReadSeqs>::Type
getSecondMateRevSeqId(TReadSeqs const & readSeqs, TPairId pairId)
{
    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    return getSecondMateFwdSeqId(readSeqs, pairId) + getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline typename Size<TReadSeqs>::Type
getReadId(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    return isFwdReadSeq(readSeqs, readSeqId) ? readSeqId : readSeqId - getReadsCount(readSeqs);
}

template <typename TReadSeqs, typename TReadSeqId>
inline typename Size<TReadSeqs>::Type
getPairId(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));

    typename Size<TReadSeqs>::Type pairId = readSeqId;

    if (isRevReadSeq(readSeqs, readSeqId))
        pairId -= getReadsCount(readSeqs);

    if (isSecondMate(readSeqs, readSeqId))
        pairId -= getPairsCount(readSeqs);

    SEQAN_ASSERT_LT(pairId, getPairsCount(readSeqs));
    
    return pairId;
}

template <typename TReadSeqs, typename TReadId>
inline typename Size<TReadSeqs>::Type
getMateId(TReadSeqs const & readSeqs, TReadId readId)
{
    typename Size<TReadSeqs>::Type pairId = getPairId(readSeqs, readId);

    return isFirstMate(readSeqs, readId) ? getSecondMateFwdSeqId(readSeqs, pairId) : getFirstMateFwdSeqId(readSeqs, pairId);
}

template <typename TReadSeqs, typename TReadSeqId>
inline typename Size<TReadSeqs>::Type
getMateSeqId(TReadSeqs const & readSeqs, TReadSeqId readSeqId)
{
    SEQAN_ASSERT_LT(readSeqId, getReadSeqsCount(readSeqs));

    typename Size<TReadSeqs>::Type pairId = getPairId(readSeqs, readSeqId);

    if (isFirstMate(readSeqs, readSeqId))
    {
        if (isFwdReadSeq(readSeqs, readSeqId))
            return getSecondMateRevSeqId(readSeqs, pairId);
        else
            return getSecondMateFwdSeqId(readSeqs, pairId);
    }
    else
    {
        if (isFwdReadSeq(readSeqs, readSeqId))
            return getFirstMateRevSeqId(readSeqs, pairId);
        else
            return getFirstMateFwdSeqId(readSeqs, pairId);
    }
}

#endif  // #ifndef APP_YARA_STORE_READS_H_
