// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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

#ifndef SEQAN_STREAM_RECORD_READER_BASE_H_
#define SEQAN_STREAM_RECORD_READER_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Mapped_;
typedef Tag<Mapped_> Mapped;

template <typename TSpec = void>
struct SinglePass {};

template <typename TSpec = void>
struct DoublePass {};

/**
.Class.RecordReader
..cat:Input/Output
..summary:Buffer management for streams.
..signature:RecordReader<TStream, TSpec>
..param.TStream:The @Concept.Stream@ type to work on.
..param.TSpec:The record reader specialization to chose.
..see:Concept.Stream
..include:seqan/stream.h

.Memfunc.RecordReader#RecordReader
..class:Class.RecordReader
..summary:Constructor
..signature:RecordReader()
..signature:RecordReader(file[, bufferSize])
..param.file:The file/String to use for reading.
...type:nolink:The $TFile$ type of $RecordReader$.
..param.bufferSize:The size of the buffer to use.
...type:nolink:$unsigned$
 */

template <typename TStream, typename TSpec>
class RecordReader;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafuction Position
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
struct Position<RecordReader<TStream, TSpec> >
{
    typedef typename Position<TStream>::Type Type;
};

template <typename TStream, typename TSpec>
struct Position<RecordReader<TStream, TSpec> const> :
            Position<RecordReader<TStream, TSpec> >
{};

// ============================================================================
// Functions
// ============================================================================

/**
.Function.resultCode
..class:Class.RecordReader
..cat:Input/Output
..summary:Returns $int$ current status code for reader (0 on success).
..signature:resultCode(recordReader)
..param.recordReader:The @Class.RecordReader@ to query the state of.
...type:Class.RecordReader
..returns:$int$, zero if there was no error reading, non-zero value on errors. Note that zero is also returned on EOF.
..include:seqan/stream.h

.Function.value
..class:Class.RecordReader
..cat:Input/Output
..param.object.type:Class.RecordReader
..include:seqan/stream.h

.Function.goNext
..cat:Input/Output
..signature:goNext(recordReader)
..param.recordReader:The @Class.RecordReader@ to advance the position in.
...type:Class.RecordReader
..include:seqan/stream.h

.Function.nextIs
..class:Class.RecordReader
..cat:Input/Output
..signature:nextIs(recordReader, tag)
..summary:Query whether the next record is of a given type.
..param.recordReader:The @Class.RecordReader@ to peek into.
...type:Class.RecordReader
...remarks:Stays unchanged.
..param.tag:Tag to select the given record type.
..returns:$bool$ indicating whether the next record in $recordReader$ is of type given by tag.
..remarks:The checks are mostly heuristic, mostly looking at one or few characters from recordReader.
..include:seqan/stream.h

.Function.atEnd
..class:Class.RecordReader
..cat:Input/Output
..summary:Returns $true$ if there is no more data to be read.
..signature:atEnd(recordReader)
..param.recordReader:The @Class.RecordReader@ to query the state of.
...type:Class.RecordReader
..returns:This function returns $true$ if the file is at end or there was an error reading. It returns $false$ if there is more data to read. In parsing functions, you can use @Function.resultCode@ to get the result to return from your parsing function.
..see:Function.value
..see:Function.resultCode
..include:seqan/stream.h
 */

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_RECORD_READER_BASE_H_
