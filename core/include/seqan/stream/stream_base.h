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
// Direction Tags
// --------------------------------------------------------------------------

struct Input_;
typedef Tag<Input_> Input;

struct Output_;
typedef Tag<Output_> Output;

// --------------------------------------------------------------------------
// Metafunction BasicStream
// --------------------------------------------------------------------------

template <typename TValue, typename TDirection>
struct BasicStream:
    If<
        IsSameType<TDirection, Input>,
        std::basic_istream<TValue>,
        std::basic_ostream<TValue> >
{};

template <typename TDirection>
struct IosOpenMode;

template <>
struct IosOpenMode<Input>
{
    enum { VALUE = std::ios_base::in };
};

template <>
struct IosOpenMode<Output>
{
    enum { VALUE = std::ios_base::out };
};

// --------------------------------------------------------------------------
// Compression Type Tags
// --------------------------------------------------------------------------

#if SEQAN_HAS_ZLIB
struct GZFile_;
typedef Tag<GZFile_> GZFile;
#endif

#if SEQAN_HAS_BZIP2
struct BZ2File_;
typedef Tag<BZ2File_> BZ2File;
#endif

// --------------------------------------------------------------------------
// Metafunction Chunk
// --------------------------------------------------------------------------

// Chunking is not support for any object (default fallback).
template <typename TObject>
struct Chunk
{
    typedef Nothing Type;
};

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_STREAM_BASE_H_
