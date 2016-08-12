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
// Author: Lily Shellhammer
// ==========================================================================
// Class for reading/writing files in Connect (.ct) files
// ==========================================================================

#ifndef SEQAN_CONNECT_IO_FILE_H_
#define SEQAN_CONNECT_IO_FILE_H_


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
// Typedef RNAFileIn
// ----------------------------------------------------------------------------

/*!
 * @class RNAFileIn
 * @signature typedef FormattedFile<RNA, Input> ConenctFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/RNA_format_io.h>
 * @brief Class for reading RNA files.
 *
 * @see RNAHeader
 * @see RNARecord
 */
typedef FormattedFile<Connect, Input>   ConnectFileIn;

// ----------------------------------------------------------------------------
// Typedef VcfFileOut
// ----------------------------------------------------------------------------

/*!
 * @class RNAFileOut
 * @signature typedef FormattedFile<RNA, Output> ConenctFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/RNA_format_io.h>
 * @brief Class for writing RNA files.
 *
 * @see RNAHeader
 * @see RNARecord
 */

typedef FormattedFile<Connect, Output>  ConnectFileOut;

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Connect, TDirection, TSpec>, TStorageSpec>
{
    typedef RNAIOContext Type;
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
// Function readHeader(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readHeader(RNAHeader & header, FormattedFile<Connect, Input, TSpec> & file)
{
    readHeader(header, context(file), file.iter, file.format);
}


// ----------------------------------------------------------------------------
// Function readRecord(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(RNARecord & record, FormattedFile<Connect, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}



// ----------------------------------------------------------------------------
// Function writeRecord(); RNARecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<Connect, Output, TSpec> & file, RNARecord & record)
{
    writeRecord(file.iter, record, file.format);
}

//--------------------------------------------------------------------------
// Function writeHeader(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeHeader(FormattedFile<Connect, Output, TSpec> & file, RNAHeader & header)
{
    writeHeader(file.iter, header, file.format);
}


} //seqan namespace

#endif  //SEQAN_CONNECT_IO_FILE_H