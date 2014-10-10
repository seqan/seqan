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
// SmartFile for .intervals format
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_INTERVALS_IO_INTERVALS_FILE_H_
#define CORE_INCLUDE_SEQAN_INTERVALS_IO_INTERVALS_FILE_H_

#include <seqan/intervals_io/intervals_io.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef IntervalsFileIn
// ----------------------------------------------------------------------------

/*!
 * @class IntervalsFileIn
 * @extends SmartFile
 * @headerfile <seqan/ucsc_io.h>
 * @brief @link SmartFile @endlink for reading .intervals files.
 *
 * @signature typedef SmartFile<Intervals, Input> IntervalsFileIn;
 */
typedef SmartFile<Intervals, Input>   IntervalsFileIn;

// ----------------------------------------------------------------------------
// Typedef IntervalsFileOut
// ----------------------------------------------------------------------------

/*!
 * @class IntervalsFileInOut
 * @extends SmartFile
 * @headerfile <seqan/ucsc_io.h>
 * @brief @link SmartFile @endlink for reading .intervals files.
 *
 * @signature typedef SmartFile<Intervals, Output> IntervalsFileOut;
 */
typedef SmartFile<Intervals, Output> IntervalsFileOut;


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Intervals, T> : public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Metafunction FileFormatExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileFormatExtensions<Intervals, T>
{
    static char const * VALUE[1];  // default is one extension
};

template <typename T>
char const * FileFormatExtensions<Intervals, T>::VALUE[1] =
{
    ".intervals"      // default output extension
};

// ----------------------------------------------------------------------------
// Metafunction SmartFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct SmartFileContext<SmartFile<Intervals, TDirection, TSpec>, TStorageSpec>
{
    typedef IntervalsIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<SmartFile<Intervals, TDirection, TSpec> >
{
    typedef Intervals Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn IntervalsFileIn#readRecord
 * @brief Reading records from a IntervalsFileIn.
 *
 * @signature void readRecord(record, file);
 *
 * @param[out] record The resulting @link GenomicRegion @endlink.
 * @param[out] file   The IntervalsFileIn to read from.
 *
 * @throw IOError in case of problems.
 */
template <typename TSpec>
void readRecord(GenomicRegion & record, SmartFile<Intervals, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function write()
// ----------------------------------------------------------------------------

/*!
 * @fn IntervalsFileIn#writeRecord
 * @brief Writing records to a IntervalsFileIn.
 *
 * @signature void readRecord(record, file);
 *
 * @param[out] file   The IntervalsFileIn to write to.
 * @param[out] record The @link GenomicReegion @endlink to write out.
 *
 * @throw IOError in case of problems.
 */
template <typename TSpec>
void writeRecord(SmartFile<Intervals, Output, TSpec> & file, GenomicRegion & record)
{
    writeRecord(file.iter, context(file), record, file.format);
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_INTERVALS_IO_INTERVALS_FILE_H_
