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
// This file contains routines to generate BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_BLAST_BLAST_TABULAR_WRITE_H_
#define SEQAN_BLAST_BLAST_TABULAR_WRITE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions and global const-expressions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// _firstOcc
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Size<TString>::Type
_firstOcc(TString const & str, typename Value<TString>::Type const & val)
{
    typedef typename Size<TString>::Type S;
    for (S s = 0; s < length(str); ++s)
        if (value(str, s) == val)
            return s;
    return length(str);
}

// ----------------------------------------------------------------------------
// Function _writeField() [match object given]
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeField(TStream & s,
            BlastMatch<TQId, TSId, TPos, TAlign> const & match,
            typename BlastMatchField<g>::Enum const fieldId,
            BlastFormat<f, p, g> const &)
{
    typedef BlastFormat<f, p, g> TFormat;
    typedef typename BlastMatchField<g>::Enum ENUM;
    switch (fieldId)
    {
        case ENUM::STD:
             _writeFields(s,
                          match,
                          BlastMatchField<g>::defaults,
                          TFormat(),
                          false /*appendNewline*/);
            break;
        case ENUM::Q_SEQ_ID:
            write(s, prefix(match.qId, _firstOcc(match.qId, ' ')));
            break;
//         case ENUM::Q_GI: write(s,  * ); break;
//         case ENUM::Q_ACC: write(s,  * ); break;
//         case ENUM::Q_ACCVER: write(s,  * ); break;
        case ENUM::Q_LEN:
            write(s, match.qLength);
            break;
        case ENUM::S_SEQ_ID:
            write(s, prefix(match.sId, _firstOcc(match.sId, ' ')));
            break;
//         case ENUM::S_ALL_SEQ_ID: write(s,  * ); break;
//         case ENUM::S_GI: write(s,  * ); break;
//         case ENUM::S_ALL_GI: write(s,  * ); break;
//         case ENUM::S_ACC: write(s,  * ); break;
//         case ENUM::S_ACCVER: write(s,  * ); break;
//         case ENUM::S_ALLACC: write(s,  * ); break;
        case ENUM::S_LEN:
            write(s, match.sLength);
            break;
        case ENUM::Q_START:
        {
            TPos effectiveQStart    = match.qStart;
            TPos effectiveQEnd      = match.qEnd;
            _untranslatePositions(effectiveQStart, effectiveQEnd,
                                  match.qFrameShift, match.qLength,
                                  QHasRevComp<TFormat>(),
                                  QIsTranslated<TFormat>());
            write(s, effectiveQStart);
        } break;
        case ENUM::Q_END:
        {
            TPos effectiveQStart    = match.qStart;
            TPos effectiveQEnd      = match.qEnd;
            _untranslatePositions(effectiveQStart, effectiveQEnd,
                                  match.qFrameShift, match.qLength,
                                  QHasRevComp<TFormat>(),
                                  QIsTranslated<TFormat>());
            write(s, effectiveQEnd);
        } break;
        case ENUM::S_START:
        {
            TPos effectiveSStart    = match.sStart;
            TPos effectiveSEnd      = match.sEnd;

            _untranslatePositions(effectiveSStart, effectiveSEnd,
                                  match.sFrameShift, match.sLength,
                                  SHasRevComp<TFormat>(),
                                  SIsTranslated<TFormat>());
            write(s, effectiveSStart);
        } break;
        case ENUM::S_END:
        {
            TPos effectiveSStart    = match.sStart;
            TPos effectiveSEnd      = match.sEnd;

            _untranslatePositions(effectiveSStart, effectiveSEnd,
                                  match.sFrameShift, match.sLength,
                                  SHasRevComp<TFormat>(),
                                  SIsTranslated<TFormat>());
            write(s, effectiveSEnd);
        } break;
//         case ENUM::Q_SEQ: write(s,  * ); break;
//         case ENUM::S_SEQ: write(s,  * ); break;
        case ENUM::E_VALUE:
        {
            std::string formatString;
            // imported from NCBI code
            if (match.eValue < 1.0e-180)
                formatString = "%3.1lf";
            else if (match.eValue < 1.0e-99)
                formatString = "%2.0le";
            else if (match.eValue < 0.0009)
                formatString = "%3.0le";
            else if (match.eValue < 0.1)
                formatString = "%4.3lf";
            else if (match.eValue < 1.0)
                formatString = "%3.2lf";
            else if (match.eValue < 10.0)
                formatString = "%2.1lf";
            else
                formatString = "%5.0lf";

            write(s, FormattedNumber<double>(formatString.c_str(),
                                             match.eValue));
        } break;
        case ENUM::BIT_SCORE:
        {
            std::string formatString;
            // imported from NCBI code
            if (match.bitScore > 9999)
                formatString = "%4.3le";
            else if (match.bitScore > 99.9)
                formatString = "%4.0ld";
            else
                formatString = "%4.1lf";

            write(s, FormattedNumber<double>(formatString.c_str(),
                                             match.bitScore));
        } break;
        case ENUM::SCORE:
            write(s, match.score);
            break;
        case ENUM::LENGTH:
            write(s, match.aliLength);
            break;
        case ENUM::P_IDENT:
            write(s, FormattedNumber<double>("%.2f",
                                             double(match.identities) * 100
                                             / match.aliLength));
            break;
        case ENUM::N_IDENT:
            write(s, match.identities);
            break;
        case ENUM::MISMATCH:
            write(s, match.mismatches);
            break;
        case ENUM::POSITIVE:
            write(s, match.positives);
            break;
        case ENUM::GAP_OPEN:
            write(s, match.gapOpenings);
            break;
        case ENUM::GAPS:
            write(s, match.gaps);
            break;
        case ENUM::P_POS:
            write(s, FormattedNumber<double>("%.2f",
                                             double(match.positives) * 100
                                             / match.aliLength));
            break;
        case ENUM::FRAMES:
            // for formats that don't have frames, blast says 0 instead of +1
            if (qNumFrames(TFormat()) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.qFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            write(s, '/');
            if (sNumFrames(TFormat()) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.qFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            break;
        case ENUM::Q_FRAME:
            if (qNumFrames(TFormat()) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.qFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            break;
        case ENUM::S_FRAME:
            if (sNumFrames(TFormat()) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.qFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            break;
//         case ENUM::BTOP: write( * ); break;
//         case ENUM::S_TAX_IDS: write( * ); break;
//         case ENUM::S_SCI_NAMES: write( * ); break;
//         case ENUM::S_COM_NAMES: write( * ); break;
//         case ENUM::S_BLAST_NAMES: write( * ); break;
//         case ENUM::S_S_KINGDOMS: write( * ); break;
//         case ENUM::S_TITLE: write( * ); break;
//         case ENUM::S_ALL_TITLES: write( * ); break;
//         case ENUM::S_STRAND: write( * ); break;
//         case ENUM::Q_COV_S: write( * ); break;
//         case ENUM::Q_COV_HSP:
        default:
            write(s, "n/i"); // not implemented
    };
}

// ----------------------------------------------------------------------------
// Function _writeFields() [match object given]
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFieldList,
          BlastFormatFile f,
          BlastFormatProgram p>
inline void
_writeFields(TStream & stream,
             BlastMatch<TQId, TSId, TPos, TAlign> const & match,
             TFieldList const & fields,
             BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> const &,
             bool const appendNewline = true)
{
    typedef BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> TFormat;

    for (auto it = seqan::begin(fields), itB = it, itEnd = seqan::end(fields);
         it != itEnd;
         ++it)
    {
        if (it != itB)
            write(stream, _seperatorString(TFormat()));

        _writeField(stream, match, *it, TFormat());
    }
    if (appendNewline)
        write(stream, '\n');
}

template <typename TStream,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFieldList,
          BlastFormatFile f,
          BlastFormatProgram p>
inline void
_writeFields(TStream & stream,
             BlastMatch<TQId, TSId, TPos, TAlign> const & match,
             TFieldList const &,
             BlastFormat<f, p, BlastFormatGeneration::BLAST_LEGACY> const &)
{
    // we override the Generation here and discard the list, since LEGACY
    // only supports defaults anyway
    typedef BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> TFormat;
    _writeFields(stream,
                 match,
                 BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::defaults,
                 TFormat());
}

// ----------------------------------------------------------------------------
// Function _writeFieldLabels()
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TFieldList,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeFieldLabels(TStream & stream,
                  TFieldList const & fields,
                  BlastFormat<f, p, g> const &)
{
    typedef BlastFormat<f, p, g> TFormat;
    for (auto it = seqan::begin(fields), itB = it, itEnd = seqan::end(fields);
         it != itEnd;
         ++it)
    {
        if (it != itB)
            write(stream, _seperatorString(TFormat()));

        write(stream, BlastMatchField<g>::columnLabels[static_cast<uint8_t>(*it)]);
    }

    writeValue(stream, '\n');
}

// ----------------------------------------------------------------------------
// Function _writeFields() or labels [no match object given]
// ----------------------------------------------------------------------------

template <typename TStream, BlastFormatFile f, BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeFields(TStream & /**/,
             BlastFormat<f, p, g> const & /*tag*/)
{
}

template <typename TStream, typename TField, typename... TFields,
          BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
inline void
_writeFields(TStream & stream,
             BlastFormat<f, p, g> const & /*tag*/,
             TField const & field1, TFields const & ... fields)
{
    typedef BlastFormat<f, p, g> TFormat;

    write(stream, _seperatorString(TFormat()));
    write(stream, field1);
    _writeFields(stream, TFormat(), fields... );
}

// ----------------------------------------------------------------------------
// Function writeHeader()
// ----------------------------------------------------------------------------

template <typename TStream, typename TString1, typename TString2,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
_writeHeaderWithoutColumnLabels(TStream & stream,
                                TString1 const & qryId,
                                TString2 const & dbName,
                                BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                            p,
                                            g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,p,g> TFormat;

    write(stream, "# ");
    write(stream, _programTagToString(TFormat()));
    write(stream, " I/O Module of SeqAn-");
    write(stream, SEQAN_VERSION_MAJOR);
    write(stream, '.');
    write(stream, SEQAN_VERSION_MINOR);
    write(stream, '.');
    write(stream, SEQAN_VERSION_PATCH);
    write(stream, " (http://www.seqan.de)\n# Query: ");
    write(stream, qryId);
    write(stream, "\n# Database: ");
    write(stream, dbName);
    write(stream, '\n');
}

/*!
 * @fn BlastFormat#writeHeader
 * @headerfile seqan/blast.h
 * @brief write the header for @link BlastFormatFile
 * @endlink::TABULAR_WITH_HEADER  without @link BlastRecord @endlink-object
 * @signature writeHeader(stream, qryId, dbname, [matchCount,] blastFormatTag[, labels...])
 *
 * This function writes the header of a record if @link BlastFormatFile @endlink
 * is TABULAR_WITH_HEADER (for TABULAR this is a no-op). In contrast to
 * @link BlastRecord#writeHeader @endlink no @link BlastRecord @endlink-object
 * is required. Custom column labels can be specified as variadic
 * columnlabel arguments (just pass a printable argument for each label).
 *
 * @param[in,out] stream    The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] qryId     The full ID of the query sequence (e.g. FASTA
 * one-line description)
 * @param[in] dbName    The name of the database (NOT the ID of a subject sequence,
 * but the name of the database, e.g. "NCBI NR" or path to a file)
 * @param[in] matchCount The amount of matches in the record. Mandatory parameter
 * for @link BlastFormatGeneration::BLAST_PLUS @endlink, optional for
 * @link BlastFormatGeneration::BLAST_LEGACY @endlink
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 * @param[in] labels...   Optional custom column labels
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 */

// default fields
template <typename TStream,
          typename TString,
          typename TString2,
          typename TNumeric,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TString const & qryId,
            TString2 const & dbName,
            TNumeric const hitCount,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

    _writeHeaderWithoutColumnLabels(stream, qryId, dbName, TFormat());

    write(stream, "# Fields: ");
    write(stream,
          BlastMatchField<g>::columnLabels
          [static_cast<uint8_t>(BlastMatchField<g>::Enum::STD)/* == 0 */]);
    write(stream, '\n');

    if (g == BlastFormatGeneration::BLAST_PLUS)
    {
        write(stream, "# ");
        write(stream, hitCount);
        write(stream, " hits found\n");
    }
}

// custom fields
// this is not exposed by the documentation above, but only 
template <typename TStream,
          typename TString,
          typename TString2,
          typename TNumeric,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TString const & qryId,
            TString2 const & dbName,
            TNumeric const hitCount,
            TFieldList const & fieldList,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "writeMatch() with custom fields only supported with "
                  "BlastFormatGeneration::BLAST_PLUS; use default fields or "
                  "the third signature (arbitrary columns) instead.");
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    _writeHeaderWithoutColumnLabels(stream, qryId, dbName, TFormat());

    write(stream, "# Fields: ");
    _writeFieldLabels(stream, fieldList, TFormat());

    write(stream, "# ");
    write(stream, hitCount);
    write(stream, " hits found\n");
}

// custom columns
template <typename TStream,
          typename TString,
          typename TString2,
          typename TNumeric,
          typename TField,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TString const & qryId,
            TString2 const & dbName,
            TNumeric const hitCount,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /**/,
            TField const & field1,
            TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

    _writeHeaderWithoutColumnLabels(stream, qryId, dbName, TFormat());

    write(stream, "# Fields: ");
    write(stream, field1);
    _writeFields(stream, TFormat(), fields...);
    write(stream, '\n');

    if (g == BlastFormatGeneration::BLAST_PLUS)
    {
        write(stream, "# ");
        write(stream, hitCount);
        write(stream, " hits found\n");
    }
}

// BlastFormatGeneration::BLAST_LEGACY works without hitCount aswell
template <typename TStream,
          typename TString1,
          typename TString2,
          BlastFormatProgram p,
          typename std::enable_if<IsSequence<TString1>::VALUE>::type = 0>
inline void
writeHeader(TStream & stream,
            TString1 const & qryId,
            TString2 const & dbName,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;

    writeHeader(stream, qryId, dbName, 0, TFormat());
}

// BlastFormatGeneration::BLAST_LEGACY works without hitCount aswell, custom fields
template <typename TStream,
          typename TString,
          typename TString2,
          typename TField,
          typename... TFields,
          BlastFormatProgram p>
inline void
writeHeader(TStream & stream,
            TString const & qryId,
            TString2 const & dbName,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> const & /**/,
            TField const & field1,
            TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;

    writeHeader(stream, qryId, dbName, 0, TFormat(), field1, fields...);
}

/*!
 * @fn BlastRecord#writeHeader
 * @headerfile seqan/blast.h
 * @brief write the header of a @link BlastRecord @endlink to file
 * @signature writeHeader(stream, blastRecord, dpSpecs, [fieldList,] blastFormatTag)
 *
 * This function writes the header of a record if @link BlastFormatFile @endlink
 * is TABULAR_WITH_HEADER (for TABULAR this is a no-op). Custom column labels
 * can be specified by passing a sequence of
 * @link BlastMatchField::Enum @endlink. Use this in conjunction with the same 
 * options of @link BlastMatch#writeMatch @endlink .
 *
 * @param[in,out] stream     The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] blastRecord    The @link BlastRecord @endlink whose header you want to print.
 * @param[in] dbSpecs        A @link BlastDbSpecs @endlink object with at least .dbname set.
 * @param[in] fieldList      Sequence of @link BlastMatchField::Enum @endlink
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 */

// Record parameter
template <typename TStream,
          typename TRecord,
          typename TSpec,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TRecord const & record,
            BlastDbSpecs<TSpec> const & dbSpecs,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g>  TFormat;
    writeHeader(stream,
                record.qId,
                dbSpecs.dbName,
                record.matches.size(),
                TFormat());
}

// Record parameter, custom fields
template <typename TStream,
          typename TRecord,
          typename TSpec,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & stream,
            TRecord const & record,
            BlastDbSpecs<TSpec> const & dbSpecs,
            TFieldList const & fieldList,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /**/)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "writeMatch() with custom fields only supported with "
                  "BlastFormatGeneration::BLAST_PLUS; use default fields or "
                  "the second signature (arbitrary columns) instead.");
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    writeHeader(stream, record.qId, dbSpecs.dbName,
                record.matches.size(), fieldList, TFormat());
}

// NO-OPS for TABULAR
template <typename TStream, typename TString1, typename TString2,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeHeader(TStream &, TString1 const &, TString2 const &,
            BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
}

template <typename TStream,
          typename T,
          typename T2,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & /**/,
            T const & /**/,
            T2 const & /**/,
            BlastFormat<BlastFormatFile::TABULAR,p,g> const & /*tag*/,
            TFields const & ... /**/)
{
}

template <typename TStream,
          typename T,
          typename T2,
          typename T3,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & /**/,
            T const & /**/,
            T2 const & /**/,
            T3 const & /**/,
            BlastFormat<BlastFormatFile::TABULAR,p,g> const & /*tag*/,
            TFields const & ... /**/)
{
}

template <typename TStream,
          typename T,
          typename T2,
          typename T3,
          typename T4,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TStream & /**/,
            T const & /**/,
            T2 const & /**/,
            T3 const & /**/,
            T4 const & /**/,
            BlastFormat<BlastFormatFile::TABULAR,p,g> const & /*tag*/,
            TFields const & ... /**/)
{
}

// ----------------------------------------------------------------------------
// Function writeMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastMatch#writeMatch
 * @headerfile seqan/blast.h
 * @brief write a @link BlastMatch @endlink to file
 * @signature writeMatch(stream, blastMatch, [fields,] blastFormatTag)
 *
 * This function writes a single match if @link BlastFormatFile @endlink
 * is BlastFormatFile::TABULAR_WITH_HEADER or
 * BlastFormatFile::TABULAR. Not specifying <tt>fields</tt> causes default
 * columns to be printed, which is described <a href="https://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/blastall/blastall_node93.html">
 * here</a>. Please note that BLAST is 1-indexed and considers the last position
 * to be the back, not the end, i.e. last one included in a match/sequence/...,
 * not the one behind it (as SeqAn does); this functions corrects for both of
 * these bahaviours, so you don't have to. Additionally, based on your
 * @link BlastFormatProgram @endlink, positions are transformed back to DNA space, if
 * translation has taken place.
 * Please note also that query and subject IDs are truncated at the first space
 * character in NCBI BLAST, this is also done by default here.
 *
 * By passing the <tt>fields</tt> variable you manually specify which columns you want to
 * print; the same conversions mentioned above will me made. See
 * @link BlastMatchField::Enum @endlink for a list of fields available. <tt>fields</tt>
 * is only available with @link BlastFormatGeneration @endlink ==
 * @link BlastFormatGeneration::BLAST_PLUS @endlink.
 *
 * Many guides recommend always printing the default 12 columns and using only
 * additional columns with additional (custom) data.
 * 
 * Please see @link BlastFormat#writeMatch @endlink for an implementation that
 * does not require a @link BlastMatch @endlink object.
 *
 * @param[in,out] stream    The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] blastMatch    The @link BlastMatch @endlink you wish to print.
 * @param[in] fields        A Sequence of @link BlastMatchField::Enum @endlink
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 * @see BlastRecord#writeHeader
 */

// Function for match object and manually specified fields
template <typename TStream,
          typename TBlastMatch,
          typename TBlastMatchFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           TBlastMatch const & match,
           TBlastMatchFields const & fields,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "writeMatch() with custom fields only supported with "
                  "BlastFormatGeneration::BLAST_PLUS; use default fields or "
                  "the second signature (arbitrary columns) instead.");
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    _writeFields(stream, match, fields, TFormat());
}

// Function for match object and default colummns
template <typename TStream, typename TBlastMatch,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           TBlastMatch const & match,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _writeFields(stream, match, BlastMatchField<g>::defaults, TFormat());
}

// writeMatch TABULAR_WITH_HEADER equal to TABULAR
template <typename TStream,
          typename TBlastMatch,
          typename TBlastMatchFields,
          BlastFormatProgram p>
inline void
writeMatch(TStream & stream,
           TBlastMatch const & match,
           TBlastMatchFields const & fields,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_PLUS> const &/**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    writeMatch(stream, match, fields, TFormat());
}

template <typename TStream, typename TBlastMatch,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           TBlastMatch const & match,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &/**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    writeMatch(stream, match, TFormat());
}

/*!
 * @fn BlastFormat#writeMatch
 * @headerfile seqan/blast.h
 * @brief write blast tabular output without a @link BlastMatch @endlink object
 * @signature writeMatch(stream, blastFormatTag, columns...)
 *
 * This function writes a single match if @link BlastFormatFile @endlink
 * is TABULAR_WITH_HEADER or TABULAR. In contrast to 
 * @link BlastMatch#writeMatch @endlink (which works on @link BlastMatch @endlink)
 * this signature allows an arbitrary amount of and
 * arbitrary typed columns to be printed. If you do this and you use the
 * TABULAR_WITH_HEADER format, you should also print custom column labels
 * with @link BlastFormat#writeHeader @endlink. Please note that this function
 * does none of the adjustments mentioned in the other signature,
 * so they have to be done by yourself (if you use
 * the mentioned fields and you want compatability).
 *
 * Many guides recommend always printing the default 12 columns and using only
 * additional columns with additional (custom) data.
 *
 * @param[in,out] stream    The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in] blastFormatTag The @link BlastFormat @endlink specifier.
 * @param[in] columns...   Custom columns
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 * @see BlastRecord#writeHeader
 * @see BlastFormat#writeHeader
 */

// Function for arbitrary number and typed fields
template <typename TStream, typename TField, typename... TFields,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           BlastFormat<BlastFormatFile::TABULAR,
                       p,
                       g> const & /*tag*/,
            TField const & field1,
            TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    write(stream, field1);

    _writeFields(stream, TFormat(), fields...);
    write(stream, '\n');
}

// writeMatch TABULAR_WITH_HEADER equal to TABULAR
template <typename TStream,
          typename TField,
          typename... TFields,
          BlastFormatProgram p, BlastFormatGeneration g>
inline void
writeMatch(TStream & stream,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const & /**/,
           TField const & field1,
           TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    writeMatch(stream, TFormat(), field1, fields... );
}

// ----------------------------------------------------------------------------
// Function writeTop()
// ----------------------------------------------------------------------------

// TABULAR formats have no special "top"

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

// dox in blast_base.h

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeRecordImplTab(TStream                    & stream,
                    TRecord              const & record,
                    TDbSpecs             const & dbSpecs,
                    BlastFormat<f, p, g> const & /*tag*/)
{
    typedef BlastFormat<f, p, g> TFormat;

    //TODO if debug, do lots of sanity checks on record

    //NOOP for TABULAR
    writeHeader(stream, record, dbSpecs, TFormat());
    for (auto it = record.matches.begin(); it != record.matches.end(); ++it)
    {
        //SOME SANITY CHECKS
        SEQAN_ASSERT_EQ(CharString(record.qId), CharString(it->qId));

        writeMatch(stream, *it, TFormat());
    }
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          typename TFieldList,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeRecordImplTab(TStream                    & stream,
                    TRecord              const & record,
                    TDbSpecs             const & dbSpecs,
                    TFieldList           const & fields,
                    BlastFormat<f, p, g> const & /*tag*/)
{
    typedef BlastFormat<f, p, g> TFormat;

    //TODO if debug, do lots of sanity checks on record

    //NOOP for TABULAR
    writeHeader(stream, record, dbSpecs, fields, TFormat());
    for (auto it = record.matches.begin(); it != record.matches.end(); ++it)
    {
        //SOME SANITY CHECKS
        SEQAN_ASSERT_EQ(CharString(record.qId), CharString(it->qId));

        writeMatch(stream, *it, fields, TFormat());
    }
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _writeRecordImplTab(stream, record, dbSpecs, TFormat());
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    _writeRecordImplTab(stream, record, dbSpecs, TFormat());
}


template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            TFieldList    const & fields,
            BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _writeRecordImplTab(stream, record, dbSpecs, fields, TFormat());
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            TFieldList    const & fields,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    _writeRecordImplTab(stream, record, dbSpecs, fields, TFormat());
}

// ----------------------------------------------------------------------------
// Function writeBottom()
// ----------------------------------------------------------------------------


// template <typename TStream,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// writeBottom(TStream                           & /**/,
//             TDbSpecs                    const & /**/,
//             TBlastScoringAdapater       const & /**/,
//             BlastFormat<BlastFormatFile::TABULAR, p,g> const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
// // }
// 
// template <typename TStream,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// writeBottom(TStream                           & /**/,
//             TDbSpecs                    const & /**/,
//             TBlastScoringAdapater       const & /**/,
//             BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p,g> const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
// // }

} // namespace seqan
#endif // header guard
