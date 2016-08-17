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
// Class for writing files in dot-bracket files
// ==========================================================================

#ifndef SEQAN_RNA_DOT_BRACKET_IO_FILE_H_
#define SEQAN_RNA_DOT_BRACKET_IO_FILE_H_


namespace seqan {

// ============================================================================
// Classes, Tags
// ============================================================================
// ============================================================================
// Typedefs
// ============================================================================

// ----------------------------------------------------------------------------
// Typedef DotBracketFileIn
// ----------------------------------------------------------------------------

/*!
 * @class DotBracketFileIn
 * @signature typedef FormattedFile<DotBracket, Input> DotBracketFileIn;
 * @extends FormattedFileIn
 * @recordfile <seqan/Rna_io.h>
 * @brief Class for reading DotBracket files.
 *
 * @see DotBracketrecord
 * @see DotBracketRecord
 */

typedef FormattedFile<DotBracket, Input>   DotBracketFileIn;

// ----------------------------------------------------------------------------
// Typedef DotBracketFileOut
// ----------------------------------------------------------------------------

/*!
 * @class DotBracketFileOut
 * @signature typedef FormattedFile<DotBracket, Output> DotBracketFileOut;
 * @extends FormattedFileOut
 * @recordfile <seqan/Rna_io.h>
 * @brief Class for writing DotBracket files.
 *
 * @see DotBracketrecord
 * @see DotBracketRecord
 */

typedef FormattedFile<DotBracket, Output>  DotBracketFileOut;

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<DotBracket, TDirection, TSpec>, TStorageSpec>
{
    typedef RnaIOContext Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<DotBracket, TDirection, TSpec> >
{
    typedef DotBracket Type;
};

// ----------------------------------------------------------------------------
// Function writeRecord(); DotBracketRecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
writeRecord(FormattedFile<DotBracket, Output, TSpec> & file, RnaRecord & record)
{
    writeRecord(file.iter, record, file.format);
}

// ----------------------------------------------------------------------------
// Function readrecord(); Rnarecord
// ----------------------------------------------------------------------------

template <typename TSpec>
inline void
readRecord(RnaRecord & record, FormattedFile<DotBracket, Input, TSpec> & file)
{
    readRecord(record, context(file), file.iter, file.format);
}



} //seqan namespace

#endif	//SEQAN_RNA_DOT_BRACKET_IO_FILE_H_