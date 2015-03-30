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
 * @remarks
 *
 * This is a FormattedFile abstraction of the BlastIO module. It is the
 * interface most basic (and easy to use), but it is slightly slower than other
 * interfaces and supports less features. See @link BlastFormat @endlink
 * for possibilities to write Blast compatible files with more fine-grained
 * control.
 *
 * The subset of @link BlastFormat @endlinks supported are those that
 * have @link BlastFormatFile @endlink == ::TABULAR or
 * ::TABULAR_WITH_HEADER and @link BlastFormatGeneration @endlink ==
 * ::BLAST_PLUS.
 *
 * It only contains the @link FormattedFileOut#writeRecord @endlink
 * function (@link FormattedFileOut#writeHeader @endlink is a NOOP).
 * Before calling it, make sure that the database name inside the context
 * is set (see below).
 *
 * @example
 * @code{.cpp}
 * BlastTabularOut out("/tmp/example.blast");
 *
 * context(out).dbSpecs.dbName = "Legendary Nucleotide Database";
 *
 * BlastRecord<> r;
 * r.qId = "FIRSTREAD abcdefg";
 *
 * for (...)
 * {
 *     BlastMatch<> m;
 *
 *     // "fill" the match object
 *
 *     appendValue(r.matches, m);
 * }
 *
 * writeRecord(out, r);
 * @endcode
 *
 * @see BlastRecord
 */

typedef FormattedFile<BlastTabular, Output> BlastTabularOut;

// ============================================================================
// Typedefs
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function guessFormat()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool guessFormat(FormattedFile<BlastTabular, Input, TSpec> & file)
{
    if (value(file.iter) == '#')
        setFormat(file,
                  BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                              BlastFormatProgram::BLASTN,
                              BlastFormatGeneration::BLAST_PLUS>());
    else
        setFormat(file,
                  BlastFormat<BlastFormatFile::TABULAR,
                              BlastFormatProgram::BLASTN,
                              BlastFormatGeneration::BLAST_PLUS>());
    // Always BLASTN, but doesn't matter for reading
    return true;
}

template <typename TSpec>
inline bool guessFormat(FormattedFile<BlastTabular, Output, TSpec> &)
{
    return true;
}

// ----------------------------------------------------------------------------
// Function guessFormatFromFilename()
// ----------------------------------------------------------------------------

template <typename TString>
inline bool
guessFormatFromFilename(TString const &,
                        TagSelector<BlastTabularFormats> &)
{
    return true;
}

// ----------------------------------------------------------------------------
// Function setFormat()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularFileIn#setFormat
 * @brief Convenience function in addition to FormattedFile#setFormat
 *
 * @signature void setFormat(file, tag);
 *
 * @param[in,out] file The BlastTabularFileIn to change.
 * @param[in]     tag  The @link BlastFormat @endlink to set.
 */

/*!
 * @fn BlastTabularFileOut#setFormat
 * @brief Convenience function in addition to FormattedFile#setFormat
 *
 * @signature void setFormat(file, tag);
 *
 * @param[in,out] file The BlastTabularFileOut to change.
 * @param[in]     tag  The @link BlastFormat @endlink to set.
 */

// template <typename TDirection,
//           typename TSpec,
//           BlastFormatFile f,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// setFormat(FormattedFile<BlastTabular, TDirection, TSpec> & file,
//           BlastFormat<f, p, g> const &,
//           TagList<BlastFormat<f, p, g>> const &)
// {
//     file.format.tagId = LENGTH<BlastTabularFormats>::VALUE - 1;
// }
//
// template <typename TDirection,
//           typename TSpec,
//           typename TSubList,
//           BlastFormatFile f,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// setFormat(FormattedFile<BlastTabular, TDirection, TSpec> & file,
//           BlastFormat<f, p, g> const &,
//           TagList<BlastFormat<f, p, g>, TSubList> const &)
// {
//     file.format.tagId = LENGTH<BlastTabularFormats>::VALUE -
//                         LENGTH<TSubList>::VALUE + 1;
// }
//
// template <typename TDirection,
//           typename TSpec,
//           typename TTag,
//           BlastFormatFile f,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// setFormat(FormattedFile<BlastTabular, TDirection, TSpec> &,
//           BlastFormat<f, p, g> const &,
//           TagList<TTag> const &)
// {
//     SEQAN_FAIL("Given BlastFormat cannot be used as tag for FormattedFile.");
// }
//
// template <typename TDirection,
//           typename TSpec,
//           typename TTag,
//           typename TSubList,
//           BlastFormatFile f,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// setFormat(FormattedFile<BlastTabular, TDirection, TSpec> & file,
//           BlastFormat<f, p, g> const &,
//           TagList<TTag, TSubList> const &)
// {
//     typedef BlastFormat<f, p, g> TFormat;
//     setFormat(file, TFormat(), TSubList());
// }

template <typename TDirection,
          typename TSpec,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
setFormat(FormattedFile<BlastTabular, TDirection, TSpec> & file,
          BlastFormat<f, p, g> const &)
{
    typedef BlastFormat<f, p, g> TFormat;
    assign(file.format, TFormat());
}

// ----------------------------------------------------------------------------
// Function readRecord(); BlastRecord
// ----------------------------------------------------------------------------

// template <typename TQId,
//           typename TSId,
//           typename TPos,
//           typename TAlign,
//           typename TSpec,
//           BlastFormatProgram p>
// inline void
// readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & record,
//            FormattedFile<BlastTabular, Input, TSpec> & file,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
//                        p,
//                        BlastFormatGeneration::BLAST_PLUS> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
//                         p,
//                         BlastFormatGeneration::BLAST_PLUS> TFormat;
//     readRecord(record, context(file).dbSpecs.dbName, file.iter,
//                context(file), TFormat());
// }
//
// template <typename TQId,
//           typename TSId,
//           typename TPos,
//           typename TAlign,
//           typename TSpec,
//           BlastFormatProgram p>
// inline void
// readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & record,
//            FormattedFile<BlastTabular, Input, TSpec> & file,
//            BlastFormat<BlastFormatFile::TABULAR,
//                        p,
//                        BlastFormatGeneration::BLAST_PLUS> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR,
//                         p,
//                         BlastFormatGeneration::BLAST_PLUS> TFormat;
//     readRecord(record, file.iter, context(file), TFormat());
// }

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> &,
           FormattedFile<BlastTabular, Input, TSpec> &,
           TagSelector<> const &)
{
    SEQAN_FAIL("Invalid TagId");
}

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec,
          typename TTagList>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & record,
           FormattedFile<BlastTabular, Input, TSpec> & file,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        readRecord(record,  file.iter, context(file), TFormat());
    else
        readRecord(record,
                   file,
                   static_cast<typename TagSelector<TTagList>::Base const &>(format));
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
    readRecord(record, file, file.format);
}

// ----------------------------------------------------------------------------
// Function writeRecord(); Blast
// ----------------------------------------------------------------------------

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec>
inline void
writeRecord(FormattedFile<BlastTabular, Output, TSpec> & ,
            BlastRecord<TQId, TSId, TPos, TAlign> const & ,
            TagSelector<> const &)
{
    SEQAN_FAIL("Invalid TagId");
}

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec,
          typename TTagList>
inline void
writeRecord(FormattedFile<BlastTabular, Output, TSpec> & file,
            BlastRecord<TQId, TSId, TPos, TAlign> const & record,
            TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;

    if (isEqual(format, TFormat()))
        writeRecord(file.iter, record, context(file).dbSpecs, TFormat());
    else
        writeRecord(file,
                    record,
                    static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// convient BlastFile variant
template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TSpec>
inline void
writeRecord(FormattedFile<BlastTabular, Output, TSpec> & file,
            BlastRecord<TQId, TSId, TPos, TAlign> const & record)
{
    writeRecord(file, record, file.format);
}

} // namespace seqan

#endif // SEQAN_BLAST_READ_BLAST_TABULAR_H_
