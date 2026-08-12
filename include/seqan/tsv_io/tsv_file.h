// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: David Weese <dave.weese@gmail.com>
// ==========================================================================
// Class for reading/writing files in TSV (tab-separated values) format.
// ==========================================================================

#ifndef SEQAN_TSV_IO_TSV_FILE_H_
#define SEQAN_TSV_IO_TSV_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Tsv_;
typedef Tag<Tsv_> Tsv;

// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef TsvFileIn
// ----------------------------------------------------------------------------

/*!
 * @class TsvFileIn
 * @signature typedef FormattedFile<Tsv, Input> TsvFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/tsv_io.h>
 * @brief Class for reading TSV files.
 *
 * @see TsvHeader
 * @see TsvRecord
 */

typedef FormattedFile<Tsv, Input>   TsvFileIn;

// ----------------------------------------------------------------------------
// Typedef TsvFileOut
// ----------------------------------------------------------------------------

/*!
 * @class TsvFileOut
 * @signature typedef FormattedFile<Tsv, Output> TsvFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/tsv_io.h>
 * @brief Class for writing TSV files.
 *
 * @see TsvHeader
 * @see TsvRecord
 */

typedef FormattedFile<Tsv, Output>  TsvFileOut;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Tsv, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Tsv, T>::VALUE[1] =
{
    ".tsv"     // default output extension
};

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Tsv, T> :
    public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Tsv, TDirection, TSpec>, TStorageSpec>
{
    typedef TsvIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Tsv, TDirection, TSpec> >
{
    typedef Tsv Type;
};

// ----------------------------------------------------------------------------
// Function readHeader(); TsvHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readHeader(TsvHeader & header, FormattedFile<Tsv, Input, TSpec> & file)
{
    readHeader(header, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecord(); TsvRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(TsvRecord & record, FormattedFile<Tsv, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function writeHeader(); TsvHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeHeader(FormattedFile<Tsv, Output, TSpec> & file, TsvHeader & header)
{
    writeHeader(file.iter, header, context(file), file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); TsvRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<Tsv, Output, TSpec> & file, TsvRecord & record)
{
    writeRecord(file.iter, record, context(file), file.format);
}

}  // namespace seqan

#endif // SEQAN_TSV_IO_TSV_FILE_H_
