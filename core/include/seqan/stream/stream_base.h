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
// ==========================================================================
// Base class for streams.
//
// See header stream_put.h for the generic implementation of streamPut().
// ==========================================================================

#ifndef SEQAN_STREAM_STREAM_BASE_H_
#define SEQAN_STREAM_STREAM_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class Stream
 * @implements StreamConcept
 * @headerfile <seqan/stream.h>
 * @brief Abstract base class to fulfill the @link StreamConcept stream concept @endlink.
 *
 * @signature template <typename TSpec>
 *            class Stream;
 *
 * @tparam TSpec The specializing type.
 */

/**
.Class.Stream
..cat:Input/Output
..implements:Concept.StreamConcept
..signature:Stream<TSpec>
..summary:Abstract base class to fulfill the @Concept.StreamConcept@ concept.
..concept:Concept.StreamConcept
..include:seqan/stream.h
 */

template <typename TPointer = char *>
struct CharArray;

#if SEQAN_HAS_ZLIB  // Enable Stream<GZFile> if available.
struct GZFile_;
typedef Tag<GZFile_> GZFile;
#endif  // #if SEQAN_HAS_ZLIB

#if SEQAN_HAS_BZIP2  // Enable Stream<BZ2File> if available.
struct BZ2File_;
typedef Tag<BZ2File_> BZ2File;
#endif  // #if SEQAN_HAS_ZLIB

template <typename TSpec>
class Stream;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/*!
 * @fn Stream#open
 * @brief Open a stream.
 *
 * @signature bool open(stream, fileName, mode);
 *
 * @param[in,out] stream   The stream to open.
 * @param[in]     fileName The path to the file to open.  Type: <tt>char const *</tt>.
 * @param[in]     mode     The mode for opening.
 */

/*!
 * @fn Stream#close
 * @brief Close a stream.
 *
 * @signature void close(stream);
 *
 * @param[in,out] stream The Stream to close.
 */

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

///.Function.atEnd.param.iterator.type:Class.Stream
///.Function.atEnd.class:Class.Stream

template <typename TSpec>
inline bool
atEnd(Stream<TSpec> & stream)
{
    return streamEof(stream);
}

template <typename TSpec>
inline bool
atEnd(Stream<TSpec> const & stream)
{
    return streamEof(stream);
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_STREAM_BASE_H_
