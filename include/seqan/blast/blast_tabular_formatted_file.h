// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// This file contains routines to read BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_BLAST_BLAST_TABULAR_FORMATTED_FILE_H_
#define SEQAN_BLAST_BLAST_TABULAR_FORMATTED_FILE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct BlastTabular_;
typedef Tag<BlastTabular_> BlastTabular;

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

// TODO adapt magicheader?
template <typename T>
struct MagicHeader<BlastTabular, T> :
    public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

//TODO change this to static constexpr
template <typename T>
struct FileExtensions<BlastTabular, T>
{
    static char const * VALUE[5];    // default is one extension
};

template <typename T>
char const * FileExtensions<BlastTabular, T>::VALUE[1] =
{
    ".blast",
    ".m8",
    ".bm8",
    ".m9",
    ".bm9"
};


// ----------------------------------------------------------------------------
// Type BlastFileIn
// ----------------------------------------------------------------------------

/*!
 * @class BlastTabularFileIn
 * @signature typedef FormattedFile<BlastTabular, Input> BlastTabularIn;
 * @extends FormattedFileIn
 * @headerfile <seqan/blast.h>
 * @brief FormattedFileIn abstraction for a subset of BlastFormats
 *
 * Only @link BlastFormatFile @endlink::TABULAR or
 * @link BlastFormatFile @endlink::TABULAR_WITH_HEADER and
 * @link BlastFormatGeneration @endlink == ::BLAST_PLUS
 * are supported with this interface. For more options, see
 * @link BlastFormat @endlink.
 *
 * @see BlastRecord
 */

typedef FormattedFile<BlastTabular, Input> BlastTabularIn;

// ----------------------------------------------------------------------------
// Type BlastFileOut
// ----------------------------------------------------------------------------

/*!
 * @class BlastTabularFileOut
 * @signature typedef FormattedFile<BlastTabular, Output> BlastTabularOut;
 * @extends FormattedFileOut
 * @headerfile <seqan/blast.h>
 * @brief FormattedFileOut abstraction for a subset of BlastFormats
 *
 * Only @link BlastFormatFile @endlink::TABULAR or
 * @link BlastFormatFile @endlink::TABULAR_WITH_HEADER and
 * @link BlastFormatGeneration @endlink == ::BLAST_PLUS
 * are supported with this interface. For more options, see
 * @link BlastFormat @endlink.
 * @see BlastRecord
 */

typedef FormattedFile<BlastTabular, Output> BlastTabularOut;

// ============================================================================
// Typedefs
// ============================================================================


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FormattedFileContext
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec, typename TStorageSpec>
struct FormattedFileContext<FormattedFile<Blast, TDirection, TSpec>, TStorageSpec>
{
    typedef BlastDbSpecs<> Type;
};

// ----------------------------------------------------------------------------
// Metafunction FileFormats
// ----------------------------------------------------------------------------

template <typename TDirection, typename TSpec>
struct FileFormat<FormattedFile<Blast, TDirection, TSpec> >
{
    typedef std::pair<BlastFormatFile, BlastFormatProgram> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(); BlastRecord
// ----------------------------------------------------------------------------

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec,
          BlastFormatProgram p>
inline void
__readRecord<BlastFormatFile::TABULAR_WITH_HEADER>(
    BlastRecord<TQId, TSId, TPos, TAlign> & record,
    FormattedFile<BlastTabular, Input, TSpec> & file)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    readRecord(record, context(file).dbName, file.iter, TFormat();
}

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec,
          BlastFormatProgram p>
inline void
__readRecord<BlastFormatFile::TABULAR>(
    BlastRecord<TQId, TSId, TPos, TAlign> & record,
    FormattedFile<BlastTabular, Input, TSpec> & file)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    readRecord(record, file.iter, TFormat();
}

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec,
          BlastFormatFile f>
inline void
_readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & record,
            FormattedFile<BlastTabular, Input, TSpec> & file)
{
    switch(std::get<1>(file.format))
    {
        case BlastFormatProgram::BLASTN:
            return __readRecord<f, BlastFormatProgram::BLASTN>(record, file);
        case BlastFormatProgram::BLASTP:
            return __readRecord<f, BlastFormatProgram::BLASTP>(record, file);
        case BlastFormatProgram::BLASTX:
            return __readRecord<f, BlastFormatProgram::BLASTX>(record, file);
        case BlastFormatProgram::TBLASTN:
            return __readRecord<f, BlastFormatProgram::TBLASTN>(record, file);
        case BlastFormatProgram::TBLASTX:
            return __readRecord<f, BlastFormatProgram::TBLASTX>(record, file);
        default:
            SEQAN_FAIL("Not a valid BlastFormatProgram.");
    }
}

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & record,
           FormattedFile<BlastTabular, Input, TSpec> & file)
{
    switch(std::get<0>(file.format))
    {
        case BlastFormatFile::TABULAR:
            return _readRecord<BlastFormatFile::TABULAR>(record, file);
        case BlastFormatFile::TABULAR_WITH_HEADER:
            return _readRecord<BlastFormatFile::TABULAR_WITH_HEADER>(record, file);
        default:
            SEQAN_FAIL("Only TABULAR and TABULAR_WITH_HEADER supported in BlastTabularIn");
    }
}

// ----------------------------------------------------------------------------
// Function write(); BlastRecord
// ----------------------------------------------------------------------------

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec,
          BlastFormatProgram p>
inline void
__writeRecord<BlastFormatFile::TABULAR_WITH_HEADER>(
    FormattedFile<BlastTabular, Input, TSpec> & file,
    BlastRecord<TQId, TSId, TPos, TAlign> & record)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    writeRecord(file.iter, record, context(file).dbName, TFormat();
}

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec,
          BlastFormatProgram p>
inline void
__writeRecord<BlastFormatFile::TABULAR>(
    FormattedFile<BlastTabular, Input, TSpec> & file,
    BlastRecord<TQId, TSId, TPos, TAlign> & record)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    writeRecord(file.iter, record, TFormat();
}

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec,
          BlastFormatFile f>
inline void
_writeRecord(FormattedFile<BlastTabular, Output, TSpec> & file,
             BlastRecord<TQId, TSId, TPos, TAlign> & record)
{
    switch(std::get<1>(file.format))
    {
        case BlastFormatProgram::BLASTN:
            return __writeRecord<f, BlastFormatProgram::BLASTN> (file, record);
        case BlastFormatProgram::BLASTP:
            return __writeRecord<f, BlastFormatProgram::BLASTP> (file, record);
        case BlastFormatProgram::BLASTX:
            return __writeRecord<f, BlastFormatProgram::BLASTX> (file, record);
        case BlastFormatProgram::TBLASTN:
            return __writeRecord<f, BlastFormatProgram::TBLASTN>(file, record);
        case BlastFormatProgram::TBLASTX:
            return __writeRecord<f, BlastFormatProgram::TBLASTX>(file, record);
        default:
            SEQAN_FAIL("Not a valid BlastFormatProgram.");
    }
}

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec>
inline void
writeRecord(FormattedFile<BlastTabular, Output, TSpec> & file,
            BlastRecord<TQId, TSId, TPos, TAlign> & record)
{
    switch(std::get<0>(file.format))
    {
        case BlastFormatFile::TABULAR:
            return _writeRecord<BlastFormatFile::TABULAR>(file, record);
        case BlastFormatFile::TABULAR_WITH_HEADER:
            return _writeRecord<BlastFormatFile::TABULAR_WITH_HEADER>(file, record);
        default:
            SEQAN_FAIL("Only TABULAR and TABULAR_WITH_HEADER supported in BlastTabularIn");
    }
}


} // namespace seqan

#endif // SEQAN_BLAST_READ_BLAST_TABULAR_H_
