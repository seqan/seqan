// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
//         David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Basic definitions for the stream module.
// ==========================================================================

#ifndef SEQAN_STREAM_STREAM_BASE_H_
#define SEQAN_STREAM_STREAM_BASE_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction BasicStream
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection>
struct BasicStream :
    If<
        IsSameType<TDirection, Input>,
        std::basic_istream<TValue>,
        typename If<
            IsSameType<TDirection, Output>,
            std::basic_ostream<TValue>,
            std::basic_iostream<TValue>
            >::Type
        >
{};

template <typename TDirection, typename TDummy = void>
struct IosOpenMode;


template <typename TDummy>
struct IosOpenMode<Input, TDummy>
{
    static const int VALUE;
};

template <typename TDummy>
struct IosOpenMode<Output, TDummy>
{
    static const int VALUE;
};

template <typename TDummy>
struct IosOpenMode<Bidirectional, TDummy>
{
    static const int VALUE;
};


template <typename TDummy>
const int IosOpenMode<Input, TDummy>::VALUE = std::ios::in;

template <typename TDummy>
const int IosOpenMode<Output, TDummy>::VALUE = std::ios::out;

template <typename TDummy>
const int IosOpenMode<Bidirectional, TDummy>::VALUE = std::ios::in | std::ios::out;

// --------------------------------------------------------------------------
// Sequence Format Tags
// --------------------------------------------------------------------------

struct TagFasta_;
typedef Tag<TagFasta_> Fasta;

struct TagFastq_;
typedef Tag<TagFastq_> Fastq;

// --------------------------------------------------------------------------
// Compression Type Tags
// --------------------------------------------------------------------------

struct GZFile_;
typedef Tag<GZFile_> GZFile;

struct BgzfFile_;
typedef Tag<BgzfFile_> BgzfFile;

struct BZ2File_;
typedef Tag<BZ2File_> BZ2File;

// --------------------------------------------------------------------------
// Metafunction MagicHeader
// --------------------------------------------------------------------------

template <typename TTag, typename T = void>
struct MagicHeader;

template <typename T>
struct MagicHeader<Nothing, T>
{
    static unsigned char const * VALUE;
};

template <typename T>
unsigned char const * MagicHeader<Nothing, T>::VALUE = NULL;


template <typename T>
struct MagicHeader<GZFile, T>
{
    static unsigned char const VALUE[3];
};

template <typename T>
unsigned char const MagicHeader<GZFile, T>::VALUE[3] = { 0x1f, 0x8b, 0x08 };  // gzip's magic number


template <typename T>
struct MagicHeader<BgzfFile, T>
{
    static unsigned char const VALUE[3];
};

template <typename T>
unsigned char const MagicHeader<BgzfFile, T>::VALUE[3] = { 0x1f, 0x8b, 0x08 };  // gzip's magic number


template <typename T>
struct MagicHeader<BZ2File, T>
{
    static unsigned char const VALUE[3];
};

template <typename T>
unsigned char const MagicHeader<BZ2File, T>::VALUE[3] = { 0x42, 0x5a, 0x68 };  // bzip2's magic number



// TODO(weese:) The following defines makes the old guessFormat functions in file_format_mmap.h obsolete. Disable them!
template <typename T>
struct MagicHeader<Fasta, T>
{
    static unsigned char const VALUE[1];
};

template <typename T>
unsigned char const MagicHeader<Fasta, T>::VALUE[1] = { '>' };  // Fasta's first character


template <typename T>
struct MagicHeader<Fastq, T>
{
    static unsigned char const VALUE[1];
};

template <typename T>
unsigned char const MagicHeader<Fastq, T>::VALUE[1] = { '@' };  // Fastq's first character

// --------------------------------------------------------------------------
// Metafunction FileFormatExtensions
// --------------------------------------------------------------------------

// TODO(weese:) rename FileFormatExtensions to FileTypeExtensions or FileExtensions
template <typename TFormat, typename T = void>
struct FileFormatExtensions;

template <typename T>
struct FileFormatExtensions<Nothing, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileFormatExtensions<Nothing, T>::VALUE[1] =
{
    ""
};              // default output extension


template <typename T>
struct FileFormatExtensions<GZFile, T>
{
    static char const * VALUE[3];
};

template <typename T>
char const * FileFormatExtensions<GZFile, T>::VALUE[3] =
{
    ".gz",      // default output extension
    ".Z",
    ".zip"
};


template <typename T>
struct FileFormatExtensions<BgzfFile, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileFormatExtensions<BgzfFile, T>::VALUE[1] =
{
    ".bgzf"       // default output extension
};


template <typename T>
struct FileFormatExtensions<BZ2File, T>
{
    static char const * VALUE[2];
};

template <typename T>
char const * FileFormatExtensions<BZ2File, T>::VALUE[2] =
{
    ".bz2",      // default output extension
    ".bz"
};


// --------------------------------------------------------------------------
// Metafunction Chunk
// --------------------------------------------------------------------------

// Chunking is not support for any object (default fallback).
template <typename TObject>
struct Chunk
{
    typedef Nothing Type;
};

// --------------------------------------------------------------------------
// Concept InputStreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT(InputStreamConcept, (TStream))
{
    typedef typename Value<TStream>::Type       TValue;
    typedef typename Size<TStream>::Type        TSize;
    typedef typename Position<TStream>::Type    TPosition;

    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TPosition>));

    SEQAN_CONCEPT_USAGE(InputStreamConcept)
    {}
};

// --------------------------------------------------------------------------
// Concept OutputStreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT(OutputStreamConcept, (TStream))
{
    typedef typename Value<TStream>::Type       TValue;
    typedef typename Size<TStream>::Type        TSize;
    typedef typename Position<TStream>::Type    TPosition;

    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TPosition>));

    SEQAN_CONCEPT_USAGE(OutputStreamConcept)
    {}
};

// --------------------------------------------------------------------------
// Concept BidirectionalStreamConcept
// --------------------------------------------------------------------------

SEQAN_CONCEPT_REFINE(BidirectionalStreamConcept, (TStream), (InputStreamConcept)(OutputStreamConcept))
{};

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_STREAM_BASE_H_
