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
// This file contains the Contigs and ContigsLoader classes.
// ==========================================================================
// NOTE(esiragusa): the FragmentStore should provide these functionalities.

#ifndef APP_YARA_STORE_CONTIGS_H_
#define APP_YARA_STORE_CONTIGS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/random.h>

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ContigsConfig
// ----------------------------------------------------------------------------

template <typename TSpec = Alloc<> >
struct ContigsConfig
{
    typedef Dna5                        TContigValue;
    typedef Packed<TSpec>               TContigSpec;
    typedef Owner<ConcatDirect<> >      TContigsSpec;
    typedef Owner<ConcatDirect<> >      TContigNameSpec;
};

// ----------------------------------------------------------------------------
// Class Contigs
// ----------------------------------------------------------------------------

template <typename TSpec = Alloc<>, typename TConfig = ContigsConfig<TSpec> >
struct Contigs
{
    typedef typename TConfig::TContigValue              TContigValue;
    typedef typename TConfig::TContigSpec               TContigSpec;
    typedef typename TConfig::TContigsSpec              TContigsSpec;
    typedef typename TConfig::TContigNameSpec           TContigNameSpec;

    typedef String<TContigValue, TContigSpec>           TContigSeq;
    typedef StringSet<TContigSeq, TContigsSpec>         TContigSeqs;
	typedef StringSet<CharString, TContigNameSpec>      TContigNames;
	typedef NameStoreCache<TContigNames, CharString>    TContigNamesCache;

    TContigSeqs             seqs;
    TContigNames            names;
    TContigNamesCache       namesCache;

    Contigs() :
        seqs(),
        names(),
		namesCache(names)
    {}
};

// ----------------------------------------------------------------------------
// Class ContigsLoader
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = ContigsConfig<> >
struct ContigsLoader
{
    typedef std::fstream                            TStream;
    typedef RecordReader<TStream, SinglePass<> >    TRecordReader;

    TStream                         _file;
    unsigned long                   _fileSize;
    AutoSeqStreamFormat             _fileFormat;
    std::auto_ptr<TRecordReader>    _reader;

    ContigsLoader() :
        _fileSize(0)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void clear(Contigs<TSpec, TConfig> & me)
{
    clear(me.seqs);
    clear(me.names);
    clear(me.namesCache);
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TSize>
inline void reserve(Contigs<TSpec, TConfig> & me, TSize newCapacity)
{
    reserve(me.seqs, newCapacity, Exact());
}

// ----------------------------------------------------------------------------
// Function reverse()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void reverse(Contigs<TSpec, TConfig> & me)
{
    for (unsigned contigId = 0; contigId < length(me.seqs); ++contigId)
        reverse(me.seqs[contigId]);
}

// ----------------------------------------------------------------------------
// Function removeNs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TContigId, typename TRng>
inline void _removeNs(Contigs<TSpec, TConfig> & me, TContigId contigId, TRng & rng)
{
    typedef Contigs<TSpec, TConfig>                         TContigs;
    typedef typename TContigs::TContigSeqs                  TContigSeqs;
    typedef typename Value<TContigSeqs>::Type               TContigSeq;
    typedef typename Value<TContigSeq>::Type                TAlphabet;
    typedef typename Iterator<TContigSeq, Standard>::Type   TContigIt;

    TContigIt cIt = begin(me.seqs[contigId], Standard());
    TContigIt cEnd = end(me.seqs[contigId], Standard());

    while (cIt != cEnd)
    {
        for (; cIt != cEnd && value(cIt) != TAlphabet('N'); ++cIt) ;

        if (cIt == cEnd) break;

        for (; cIt != cEnd && value(cIt) == TAlphabet('N'); ++cIt)
            value(cIt) = pickRandomNumber(rng) % ValueSize<Dna>::VALUE;
    }
}

template <typename TSpec, typename TConfig>
inline void removeNs(Contigs<TSpec, TConfig> & me)
{
    Rng<MersenneTwister> rng(0xDEADBEEF);

    for (unsigned contigId = 0; contigId < length(me.seqs); ++contigId)
        _removeNs(me, contigId, rng);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileName>
inline bool open(Contigs<TSpec, TConfig> & me, TFileName const & fileName)
{
    CharString name;

    name = fileName;    append(name, ".txt");
    if (!open(me.seqs, toCString(name))) return false;

    name = fileName;    append(name, ".rid");
    if (!open(me.names, toCString(name))) return false;

    refresh(me.namesCache);

    return true;
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileName>
inline bool save(Contigs<TSpec, TConfig> const & me, TFileName const & fileName)
{
    CharString name;

    name = fileName;    append(name, ".txt");
    if (!save(me.seqs, toCString(name))) return false;

    name = fileName;    append(name, ".rid");
    if (!save(me.names, toCString(name))) return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function load()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void load(Contigs<TSpec, TConfig> & me, ContigsLoader<TSpec, TConfig> & loader)
{
    // Reserve space for contigs.
    reserve(me, loader._fileSize);

    CharString contigName;
    Dna5String contigSeq;

    // Read records.
    while (!atEnd(*(loader._reader)))
    {
        if (readRecord(contigName, contigSeq, *(loader._reader), loader._fileFormat) != 0)
            throw RuntimeError("Error while reading contig.");

        appendValue(me.seqs, contigSeq);
        appendValue(me.names, contigName);
    }

    refresh(me.namesCache);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TFileName>
inline void open(ContigsLoader<TSpec, TConfig> & loader, TFileName const & contigsFile)
{
    typedef ContigsLoader<TSpec, TConfig>            TContigsLoader;
    typedef typename TContigsLoader::TRecordReader   TRecordReader;

    // Open file.
    if (!open(loader._file, toCString(contigsFile), OPEN_RDONLY))
        throw RuntimeError("Error while opening contigs file.");

    // Compute file size.
    loader._file.seekg(0, std::ios::end);
    loader._fileSize = loader._file.tellg();
    loader._file.seekg(0, std::ios::beg);

    // Initialize record reader.
    loader._reader.reset(new TRecordReader(loader._file));

    // Autodetect file format.
    if (!guessStreamFormat(*(loader._reader), loader._fileFormat))
        throw RuntimeError("Error while guessing contigs file format.");
}

#endif  // #ifndef APP_YARA_STORE_CONTIGS_H_
