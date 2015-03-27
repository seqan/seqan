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
// Function _writeFieldLabels()
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeFieldLabels(TFwdIterator & stream,
                  BlastIOContext & context,
                  BlastFormat<f, p, g> const &)
{
    typedef BlastFormat<f, p, g> TFormat;
    if (!isEmpty(context.fieldsAsStrings)) // give preference to string labels
    {
         write(stream, concat(context.fieldsAsStrings, ", "));
    }
    else
    {
        for (auto it = seqan::begin(context.fields),
                  itB = it,
                  itEnd = seqan::end(context.fields);
            it != itEnd;
            ++it)
        {
            if (it != itB)
                write(stream, ", ");//_seperatorString(TFormat()));

            write(stream, BlastMatchField<g>::columnLabels[static_cast<uint8_t>(*it)]);
        }
    }

    writeValue(stream, '\n');
}

// ----------------------------------------------------------------------------
// Function writeHeader()
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeHeaderWithoutColumnLabels(TFwdIterator & stream,
                                BlastIOContext & context,
                                BlastRecord<TQId, TSId, TPos, TAlign> const & r,
                                BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                            p,
                                            g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,p,g> TFormat;

    write(stream, "# ");
    if (!isEmpty(context.versionString))
    {
        write(stream, context.versionString);
    }
    else
    {
        write(stream, _programTagToString(TFormat()));
        write(stream, " I/O Module of SeqAn-");
        write(stream, SEQAN_VERSION_MAJOR);
        write(stream, '.');
        write(stream, SEQAN_VERSION_MINOR);
        write(stream, '.');
        write(stream, SEQAN_VERSION_PATCH);
        write(stream, " (http://www.seqan.de)");
    }
    write(stream, "\n# Query: ");
    write(stream, r.qId);
    write(stream, "\n# Database: ");
    write(stream, context.dbSpecs.dbName);
    write(stream, '\n');
}

/*!
 * @fn BlastRecord#writeHeader
 * @headerfile seqan/blast.h
 * @brief write the header of a @link BlastRecord @endlink to file
 * @signature writeHeader(stream, context, blastRecord, tag)
 *
 * @param[in,out] stream      The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in,out] context     A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastRecord The @link BlastRecord @endlink whose header you want to print.
 * @param[in]     tag         The @link BlastFormat @endlink specifier.
 *
 * This function writes the header of a record if @link BlastFormatFile @endlink
 * is TABULAR_WITH_HEADER (for TABULAR this is a no-op).
 *
 * If context.@link BlastIOContext::versionString @endlink is set, this will be written,
 * otherwise one is generated. If either context.@link BlastIOContext::fields @endlink or
 * context.@link BlastIOContext::fieldsAsStrings @endlink is specified these will be printed
 * as column labels. If both are specified than @link BlastIOContext::fieldsAsStrings @endlink
 * are given preference. Please note that it is recommended to use
 * @link BlastIOContext::fields @endlink and not @link BlastIOContext::fieldsAsStrings @endlink to stay
 * "standards"-compliant. Also only @link BlastIOContext::fields @endlink has an influence on
 * the values printed by @link BlastMatch#writeMatch @endlink.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 */

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p>
inline void
writeHeader(TFwdIterator & stream,
            BlastIOContext & context,
            BlastRecord<TQId, TSId, TPos, TAlign> const & r,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    _writeHeaderWithoutColumnLabels(stream, context, r, TFormat());

    if (length(r.matches) > 0)
    {
        write(stream, "# Fields: ");
        _writeFieldLabels(stream, context, TFormat());
    }

    write(stream, "# ");
    write(stream, length(r.matches));
    write(stream, " hits found\n");
}

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p>
inline void
writeHeader(TFwdIterator & stream,
            BlastIOContext & context,
            BlastRecord<TQId, TSId, TPos, TAlign> const & r,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_LEGACY> TFormat;

    if ((length(context.fields) != 1) || (context.fields[0] != TEnum::STD))
    {
        if (SEQAN_ENABLE_DEBUG)
            std::cerr << "custom fields set, but will be reset for "
                         "BlastFormatGeneration::BLAST_LEGACY\n";
        clear(context.fields);
        appendValue(context.fields,
                    BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::STD);
    }

    _writeHeaderWithoutColumnLabels(stream, context, r, TFormat());

    write(stream, "# Fields: ");
    _writeFieldLabels(stream, context, TFormat());
}

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeHeader(TFwdIterator &,
            BlastIOContext &,
            BlastRecord<TQId, TSId, TPos, TAlign> const &,
            BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
    // NOOP
}

// ----------------------------------------------------------------------------
// Function _writeField() [match object given]
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeField(TFwdIterator & s,
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

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatFile f,
          BlastFormatProgram p>
inline void
_writeFields(TFwdIterator & stream,
             BlastIOContext & context,
             BlastMatch<TQId, TSId, TPos, TAlign> const & match,
             BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> const &,
             bool const appendNewline = true)
{
    typedef BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> TFormat;

    for (auto it = seqan::begin(context.fields),
              itB = it,
              itEnd = seqan::end(context.fields);
         it != itEnd;
         ++it)
    {
        if (it != itB)
            write(stream, '\t');//_seperatorString(TFormat()));

        _writeField(stream, match, *it, TFormat());
    }
    if (appendNewline)
        write(stream, '\n');
}

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatFile f,
          BlastFormatProgram p>
inline void
_writeFields(TFwdIterator & stream,
             BlastIOContext & context,
             BlastMatch<TQId, TSId, TPos, TAlign> const & match,
             BlastFormat<f, p, BlastFormatGeneration::BLAST_LEGACY> const &)
{
    // we override the Generation here
    typedef BlastFormat<f, p, BlastFormatGeneration::BLAST_PLUS> TFormat;
    _writeFields(stream, context, match, TFormat());
}

// ----------------------------------------------------------------------------
// Function writeMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastMatch#writeMatch
 * @headerfile seqan/blast.h
 * @brief write a @link BlastMatch @endlink to file
 * @signature writeMatch(stream, context, blastMatch, tag)
 *
 * @param[in,out] stream     The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in,out] context    A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastMatch The @link BlastMatch @endlink you wish to print.
 * @param[in]     tag        The @link BlastFormat @endlink specifier.
 *
 * This function writes a single match if @link BlastFormatFile @endlink
 * is BlastFormatFile::TABULAR_WITH_HEADER or
 * BlastFormatFile::TABULAR.
 * Please note that BLAST is 1-indexed and considers the last position
 * to be the back, not the end, i.e. last one included in a match/sequence/...,
 * not the one behind it (as SeqAn does); this functions corrects for both of
 * these bahaviours, so you don't have to. Additionally, based on your
 * @link BlastFormatProgram @endlink, positions are transformed back to DNA space, if
 * translation has taken place.
 * Please note also that query and subject IDs are truncated at the first space
 * character in NCBI BLAST, this is also done by default here.
 *
 * By setting context.@link BlastIOContext::fields @endlink you can specify which columns you
 * wish to print; the same conversions mentioned above will me made. See
 * @link BlastMatchField::Enum @endlink for a list of fields available. For
 * @link BlastFormatGeneration @endlink ==
 * @link BlastFormatGeneration::BLAST_LEGACY @endlink the @link BlastIOContext::fields @endlink
 * variable is ignored.
 *
 * Many guides recommend always printing the default 12 columns and using only
 * additional columns with additional (custom) data.
 *
 * Please see @link BlastFormat#writeMatch @endlink for an implementation that
 * does not require a @link BlastMatch @endlink object.
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 * @see BlastRecord#writeHeader
 */

template <typename TQId,
          typename TSId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeMatch(TFwdIterator & stream,
           BlastIOContext & context,
           BlastMatch<TQId, TSId, TPos, TAlign> const & match,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    if ((g == BlastFormatGeneration::BLAST_LEGACY) &&
        ((length(context.fields) != 1) || (context.fields[0] != TEnum::STD)))
    {
        if (SEQAN_ENABLE_DEBUG)
            std::cerr << "custom fields set, but will be reset for "
                         "BlastFormatGeneration::BLAST_LEGACY\n";
        clear(context.fields);
        appendValue(context.fields,
                    BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::STD);
    }
    _writeFields(stream, context, match, TFormat());
}

template <typename TQId,
          typename TSId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeMatch(TFwdIterator & stream,
           BlastIOContext & context,
           BlastMatch<TQId, TSId, TPos, TAlign> const & match,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;

    _writeFields(stream, context, match, TFormat());
}

// ----------------------------------------------------------------------------
// Function _writeFields() or labels [no match object given]
// ----------------------------------------------------------------------------

// template <typename TFwdIterator, BlastFormatFile f, BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// _writeFields(TFwdIterator & /**/,
//              BlastFormat<f, p, g> const & /*tag*/)
// {
// }
//
// template <typename TFwdIterator, typename TField, typename... TFields,
//           BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
// inline void
// _writeFields(TFwdIterator & stream,
//              BlastFormat<f, p, g> const & /*tag*/,
//              TField const & field1, TFields const & ... fields)
// {
//     typedef BlastFormat<f, p, g> TFormat;
//
//     write(stream, '\t');//_seperatorString(TFormat()));
//     write(stream, field1);
//     _writeFields(stream, TFormat(), fields... );
// }

/*
 * @fn BlastFormat#writeMatch
 * @headerfile seqan/blast.h
 * @brief write blast tabular output without a @link BlastMatch @endlink object
 * @signature writeMatch(stream, tag, columns...)
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
 * @param[in] tag The @link BlastFormat @endlink specifier.
 * @param[in] columns...   Custom columns
 *
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 * @see BlastRecord#writeHeader
 * @see BlastFormat#writeHeader
 */

// Function for arbitrary number and typed fields
// template <typename TFwdIterator, typename TField, typename... TFields,
//           BlastFormatProgram p, BlastFormatGeneration g>
// inline void
// writeMatch(TFwdIterator & stream,
//            BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/,
//             TField const & field1,
//             TFields const & ... fields)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
//     write(stream, field1);
//
//     _writeFields(stream, TFormat(), fields...);
//     write(stream, '\n');
// }
//
// // writeMatch TABULAR_WITH_HEADER equal to TABULAR
// template <typename TFwdIterator,
//           typename TField,
//           typename... TFields,
//           BlastFormatProgram p, BlastFormatGeneration g>
// inline void
// writeMatch(TFwdIterator & stream,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const & /**/,
//            TField const & field1,
//            TFields const & ... fields)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
//     writeMatch(stream, TFormat(), field1, fields... );
// }

// ----------------------------------------------------------------------------
// Function writeTop()
// ----------------------------------------------------------------------------

// TABULAR formats have no special "top"

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

// dox in blast_base.h

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_writeRecordTab(TFwdIterator & stream,
                BlastIOContext & context,
                BlastRecord<TQId, TSId, TPos, TAlign> const & r,
                BlastFormat<f, p, g> const & /*tag*/)
{
    typedef BlastFormat<f, p, g> TFormat;

    //TODO if debug, do lots of sanity checks on record

    //NOOP for TABULAR
    writeHeader(stream, context, r, TFormat());
    for (auto it = r.matches.begin(); it != r.matches.end(); ++it)
    {
        //SOME SANITY CHECKS
        SEQAN_ASSERT(startsWith(r.qId, it->qId));

        writeMatch(stream, context, *it, TFormat());
    }
}

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TFwdIterator & stream,
            BlastIOContext & context,
            BlastRecord<TQId, TSId, TPos, TAlign> const & r,
            BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _writeRecordTab(stream, context, r, TFormat());
}

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
writeRecord(TFwdIterator & stream,
            BlastIOContext & context,
            BlastRecord<TQId, TSId, TPos, TAlign> const & r,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    _writeRecordTab(stream, context, r, TFormat());
}

// ----------------------------------------------------------------------------
// Function writeBottom()
// ----------------------------------------------------------------------------


// template <typename TFwdIterator,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// writeBottom(TFwdIterator & /**/,
//             TDbSpecs const & /**/,
//             TBlastScoringAdapater const & /**/,
//             BlastFormat<BlastFormatFile::TABULAR, p,g> const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
// // }
// 
// template <typename TFwdIterator,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// writeBottom(TFwdIterator & /**/,
//             TDbSpecs const & /**/,
//             TBlastScoringAdapater const & /**/,
//             BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p,g> const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
// // }

} // namespace seqan
#endif // header guard
