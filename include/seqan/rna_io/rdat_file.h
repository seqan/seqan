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
// Class for reading/writing files in Rdat (.rdat) files
// ==========================================================================

#ifndef SEQAN_RDAT_IO_FILE_H_
#define SEQAN_RDAT_IO_FILE_H_


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
// Typedef RdatFileIn
// ----------------------------------------------------------------------------

/*!
 * @class RdatFileIn
 * @signature typedef FormattedFile<Rdat, Input> RdatFileIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/Rdat_io.h>
 * @brief Class for reading Rdat files.
 *
 * @see RdatHeader
 * @see RdatRecord
 */
typedef FormattedFile<Rdat, Input>   RdatFileIn;

// ----------------------------------------------------------------------------
// Typedef VcfFileOut
// ----------------------------------------------------------------------------

/*!
 * @class RdatFileOut
 * @signature typedef FormattedFile<Rdat, Output> RdatFileOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/Rdat_io.h>
 * @brief Class for writing Rdat files.
 *
 * @see RdatHeader
 * @see RdatRecord
 */

typedef FormattedFile<Rdat, Output>  RdatFileOut;

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Rdat, TDirection, TSpec>, TStorageSpec>
{
    typedef RNAIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Rdat, TDirection, TSpec> >
{
    typedef Rdat Type;
};

// ----------------------------------------------------------------------------
// Function readHeader(); RdatHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readHeader(RNAHeader & header, FormattedFile<Rdat, Input, TSpec> & file)
{
    readHeader(header, context(file), file.iter, file.format);
}


// ----------------------------------------------------------------------------
// Function readRecord(); RdatHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(RNARecord & record, FormattedFile<Rdat, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}



// ----------------------------------------------------------------------------
// Function writeRecord(); RdatRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<Rdat, Output, TSpec> & file, RNARecord & record)
{
    writeRecord(file.iter, record, file.format);
}

//--------------------------------------------------------------------------
// Function writeHeader(); RdatHeader
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeHeader(FormattedFile<Rdat, Output, TSpec> & file, RNAHeader & header)
{
    writeHeader(file.iter, header, file.format);
}


} //seqan namespace

#endif  //SEQAN_RDAT_IO_FILE_H