// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Lily Shellhammer <lily.shellhammer@gmail.com>
// ==========================================================================
// Class for reading/writing files in Connect (.ct) files
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_FILE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_FILE_H_


namespace seqan {

// ============================================================================
// Classes, Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------
// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef RnaFileIn
// ----------------------------------------------------------------------------

/*!
 * @class ConnectFileIn
 * @signature typedef FormattedFile<Rna, Input> ConenctFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/rna_format_io.h>
 * @brief Class for reading Rna files.
 *
 * @see RnaRecord
 */
typedef FormattedFile<Connect, Input> ConnectFileIn;

// ----------------------------------------------------------------------------
// Typedef ConnectFileOut
// ----------------------------------------------------------------------------

/*!
 * @class ConnectFileOut
 * @signature typedef FormattedFile<Rna, Output> ConenctFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/Rna_format_io.h>
 * @brief Class for writing Rna files.
 *
 * @see RnaRecord
 */

typedef FormattedFile<Connect, Output> ConnectFileOut;

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Connect, TDirection, TSpec>, TStorageSpec>
{
    typedef RnaIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Connect, TDirection, TSpec> >
{
    typedef Connect Type;
};

// ----------------------------------------------------------------------------
// Function readRecord(); RnaHeader, Connect File
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(RnaRecord & record, FormattedFile<Connect, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}



// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, Connect File
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<Connect, Output, TSpec> & file, RnaRecord & record)
{
    writeRecord(file.iter, record, file.format);
}

} //seqan namespace

#endif  //SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_FILE_H_
