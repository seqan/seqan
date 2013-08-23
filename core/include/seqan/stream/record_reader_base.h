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

#ifndef SEQAN_STREAM_RECORD_READER_BASE_H_
#define SEQAN_STREAM_RECORD_READER_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct StringReader_;
typedef Tag<StringReader_> StringReader;

template <typename TSpec = void>
struct SinglePass {};

template <typename TSpec = void>
struct DoublePass {};

/*!
 * @class RecordReader
 * @headerfile <seqan/stream.h>
 * @brief Buffer management when reading from @link StreamConcept streams @endlink.
 *
 * @signature template <typename TStream, typename TSpec>
 *            class RecordReader;
 *
 * @tparam TStream The @link ConceptStream @endlink to read from.
 * @tparam TSpec   The record reader specialization type.
 *
 * @section Examples
 *
 * @include demos/input_output/record_reader.cpp
 *
 * The output is as follows:
 *
 * @include demos/input_output/record_reader.cpp.stdout
 */

/*!
 * @fn RecordReader::RecordReader
 * @brief Constructor.
 *
 * @signature RecordReader::RecordReader();
 * @signature RecordReader::RecordReader(stream[, bufferSize]);
 *
 * @param[in] stream     The @link StreamConcept @endlink to read from.
 * @param[in] bufferSize The size of the buffer to use.  Type: <tt>unsigned</tt>.
 */

/**
.Class.RecordReader
..cat:Input/Output
..summary:Buffer management for streams.
..signature:RecordReader<TStream, TSpec>
..param.TStream:The @Concept.StreamConcept@ type to work on.
..param.TSpec:The record reader specialization to chose.
..see:Concept.StreamConcept
..example.file:demos/input_output/record_reader.cpp
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

/*!
 * @mfn RecordReader#Position
 * @brief Returns the position tpye to use in @link RecordReader#position @endlink and @link RecordReader#setPosition @endlink.
 *
 * @signature Position<TReader>::Type;
 *
 * @tparam TReader The RecordReader to query for its position type.
 *
 * @return Type The resulting position type.
 *
 * @see RecordReader#position
 * @see RecordReader#setPosition
 */

/**
.Metafunction.RecordReader#Position
..class:Class.RecordReader
..cat:Input/Output
..summary:Returns the position type to use in @Function.RecordReader#position@ and @Function.RecordReader#setPosition@.
..signature:Position<TReader>::Type
..param.TReader:The @Class.RecordReader@ type to query for its position type.
..returns:The position type related to the record reader type.
..include:seqan/stream.h
..see:Function.RecordReader#position
..see:Function.RecordReader#setPosition
*/

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

/*!
 * @fn RecordReader#resultCode
 * @brief Return current status code of the reader/underlying stream.
 *
 * @signature int resultCode(reader);
 *
 * @param[in] reader The RecordReader to query for its status code.
 *
 * @return int The status code, 0 for no errors, non-0 value for errors.
 */

/*!
 * @fn RecordReader#value
 * @brief Return the current value (character) of the reader.
 *
 * This is undefined if the reader is at the end or an error occured.
 *
 * @signature char value(reader);
 *
 * @param[in] reader The RecordReader to read from.
 *
 * @return char The resulting character.
 */

/*!
 * @fn RecordReader#goNext
 * @brief Advance to the next position in the stream.
 *
 * @signature bool goNext(reader);
 *
 * @param[in,out] reader The RecordReader to advance.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

/*!
 * @fn RecordReader#nextIs
 * @brief Query whether the next record is of a given type.
 * 
 * @signature bool nextIs(reader, tag)
 * 
 * @param[in,out] reader The @link RecordReader @endlink to peek into.  Remains unchanged.
 * @param[in]     tag    Tag to select the given record type.
 * 
 * @return bool Indicating whether the next record in <tt>reader</tt> is of type given by tag.
 * 
 * @section Remarks
 * 
 * The checks are mostly heuristic, mostly looking at one or few characters from reader.
 */

/*!
 * @fn RecordReader#atEnd
 * @brief Returns <tt>true</tt> if there is no more data to be read.
 * 
 * @signature bool atEnd(reader);
 * 
 * @param reader The @link RecordReader @endlink to query the state of.
 * 
 * @return bool This function returns <tt>true</tt> if the file is at end or there was an error reading.  It returns
 *              <tt>false</tt> if there is more data to read.  In parsing functions, you can use @link
 *              RecordReader#resultCode @endlink to get the result to return from your parsing function.
 * 
 * @see RecordReader#value
 * @see RecordReader#resultCode
 */

/*!
 * @fn RecordReader#position
 * @brief Returns the current position of the reader.
 *
 * @signature TPosition position(reader);
 *
 * @param[in] reader The RecordReader to query.
 *
 * @return TPosition The resulting position.  Use @link RecordReader#Position @endlink to retrieve.
 */

/*!
 * @fn RecordReader#setPosition
 * @brief Sets the current position of the reader.
 *
 * @signature void setPosition(reader, pos);
 *
 * @param[in,out] reader The RecordReader to set the position of.
 * @param[in]     pos    The position to set.
 *
 * @section Remarks
 *
 * The underlying data source has to support setting the position.  For RecordReader objects reading from strings this
 * always works.  When reading from streams, setting the position must be supported by the underlying stream
 *
 * @see RecordReader#position
 */

/**
.Function.RecordReader#resultCode
..cat:Input/Output
..class:Class.RecordReader
..signature:int resultCode(reader)
..summary:Returns $int$ current status code for reader (0 on success).
..param.reader:The @Class.RecordReader@ to query the state of.
...type:Class.RecordReader
..returns:$int$, zero if there was no error reading, non-zero value on errors. Note that zero is also returned on EOF.
..include:seqan/stream.h

.Function.RecordReader#value
..summary:Returns the current value of the reader.
..cat:Input/Output
..class:Class.RecordReader
..signature:char value(reader)
..param.reader:The @Class.RecordReader@ to get the current value from.
...type:Class.RecordReader
..include:seqan/stream.h

.Function.RecordReader#goNext
..cat:Input/Output
..class:Class.RecordReader
..signature:bool goNext(recordReader)
..summary:Advance record reader to next position.
..param.recordReader:The @Class.RecordReader@ to advance the position in.
...type:Class.RecordReader
..include:seqan/stream.h

.Function.RecordReader#nextIs
..class:Class.RecordReader
..cat:Input/Output
..signature:bool nextIs(recordReader, tag)
..summary:Query whether the next record is of a given type.
..param.recordReader:The @Class.RecordReader@ to peek into.
...type:Class.RecordReader
...remarks:Stays unchanged.
..param.tag:Tag to select the given record type.
..returns:$bool$ indicating whether the next record in $recordReader$ is of type given by tag.
..remarks:The checks are mostly heuristic, mostly looking at one or few characters from recordReader.
..include:seqan/stream.h

.Function.RecordReader#atEnd
..cat:Input/Output
..class:Class.RecordReader
..summary:Returns $true$ if there is no more data to be read.
..signature:bool atEnd(recordReader)
..param.recordReader:The @Class.RecordReader@ to query the state of.
...type:Class.RecordReader
..returns:This function returns $true$ if the file is at end or there was an error reading. It returns $false$ if there is more data to read. In parsing functions, you can use @Function.RecordReader#resultCode@ to get the result to return from your parsing function.
..see:Function.RecordReader#value
..see:Function.RecordReader#resultCode
..include:seqan/stream.h

.Function.RecordReader#position
..cat:Input/Output
..class:Class.RecordReader
..summary:Returns the current position of the reader.
..signature:TPosition position(reader)
..param.reader:The @Class.RecordReader@ to query for its position.
...type:Class.RecordReader
..param.TPosition:The @Metafunction.RecordReader#Position@ of the reader.
..returns:The current position of the record reader.
..see:Metafunction.RecordReader#Position
..include:seqan/stream.h

.Function.RecordReader#setPosition
..cat:Input/Output
..class:Class.RecordReader
..summary:Returns the current position of the reader.
..remarks:
The underlying data source has to support setting the position.
For String RecordReader objects, this works nicely, when reading from streams, setting the position of the record reader is supported if the underlying string supports it.
..signature:void setPosition(reader, pos)
..param.reader:The @Class.RecordReader@ to query for its position.
...type:Class.RecordReader
..param.pos:The position to set the reader to.
...type:Metafunction.RecordReader#Position
..returns:$void$
..see:Metafunction.RecordReader#Position
..see:Function.RecordReader#position
..include:seqan/stream.h
 */

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_RECORD_READER_BASE_H_
