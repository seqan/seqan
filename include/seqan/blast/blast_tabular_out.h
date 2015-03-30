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

// ----------------------------------------------------------------------------
// Type BlasabularOut
// ----------------------------------------------------------------------------

/*!
 * @class BlastTabularOut
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
 * function (@link FormattedFileOut#writeRecordHeader @endlink is a NOOP).
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
// Metafunctions and global const-expressions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function guessFormat()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool guessFormat(FormattedFile<BlastTabular, Output, TSpec> &)
{
    return true;
}

// ----------------------------------------------------------------------------
// Function _firstOcc
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
          typename TScore,
          typename TConString,
          BlastFormatProgram p,
          BlastTabularSpec h>
inline void
_writeFieldLabels(TFwdIterator & stream,
                  BlastIOContext<TScore, TConString, p, h> & context,
                  BlastTabular const &)
{
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
                write(stream, ", ");//_seperatorString(BlastTabular()));

            write(stream, BlastMatchField<>::columnLabels[static_cast<uint8_t>(*it)]);
        }
    }

    writeValue(stream, '\n');
}

// ----------------------------------------------------------------------------
// Function writeRecordHeader()
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastTabularSpec h>
inline void
_writeRecordHeaderWithoutColumnLabels(TFwdIterator & stream,
                                BlastIOContext<TScore, TConString, p, h> & context,
                                BlastRecord<TQId, TSId, TPos, TAlign> const & r,
                                BlastTabular const & /*tag*/)
{
    write(stream, "# ");
    if (isEmpty(context.versionString))
        context._setDefaultVersionString();
    write(stream, context.versionString);

    write(stream, "\n# Query: ");
    write(stream, r.qId);
    write(stream, "\n# Database: ");
    write(stream, context.dbSpecs.dbName);
    write(stream, '\n');
}

/*!
 * @fn BlastRecord#writeRecordHeader
 * @headerfile seqan/blast.h
 * @brief write the header of a @link BlastRecord @endlink to file
 * @signature writeRecordHeader(stream, context, blastRecord, tag)
 *
 * @param[in,out] stream      The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in,out] context     A @link BlastIOContext<TScore, TConString, p, h> @endlink with parameters and buffers.
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
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastTabularSpec h>
inline void
writeRecordHeader(TFwdIterator & stream,
                  BlastIOContext<TScore, TConString, p, h> & context,
                  BlastRecord<TQId, TSId, TPos, TAlign> const & r,
                  BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> BlastTabular();

    if (getBlastTabularSpec(context) == BlastTabularSpec::NO_HEADER)
        return;

    _writeRecordHeaderWithoutColumnLabels(stream, context, r, BlastTabular());

    if (length(r.matches) > 0)
    {
        write(stream, "# Fields: ");
        if (SEQAN_LIKELY(!context.legacyFormat))
        {
            _writeFieldLabels(stream, context, BlastTabular());
        }
        else
        {
            write(stream, BlastMatchField<>::legacyColumnLabels);

            if (SEQAN_ENABLE_DEBUG && ((length(context.fields) != 1) || (context.fields[0] != TEnum::STD)))
                std::cerr << "Attention: custom fields set, but will be ignored, because legacyFormat is also set.\n";

        }

    }
    if (SEQAN_LIKELY(!context.legacyFormat))
    {
        write(stream, "# ");
        write(stream, length(r.matches));
        write(stream, " hits found\n");
    }
}

// ----------------------------------------------------------------------------
// Function _writeField() [match object given]
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastTabularSpec h>
inline void
_writeField(TFwdIterator & s,
            BlastIOContext<TScore, TConString, p, h> & context,
            BlastMatch<TQId, TSId, TPos, TAlign> const & match,
            typename BlastMatchField<>::Enum const fieldId,
            BlastTabular const &)
{
    switch (fieldId)
    {
        case BlastMatchField<>::Enum::STD:
             _writeFields(s,
                          match,
                          BlastMatchField<>::defaults,
                          BlastTabular(),
                          false /*appendNewline*/);
            break;
        case BlastMatchField<>::Enum::Q_SEQ_ID:
            write(s, prefix(match.qId, _firstOcc(match.qId, ' ')));
            break;
//         case ENUM::Q_GI: write(s,  * ); break;
//         case ENUM::Q_ACC: write(s,  * ); break;
//         case ENUM::Q_ACCVER: write(s,  * ); break;
        case BlastMatchField<>::Enum::Q_LEN:
            write(s, match.qLength);
            break;
        case BlastMatchField<>::Enum::S_SEQ_ID:
            write(s, prefix(match.sId, _firstOcc(match.sId, ' ')));
            break;
//         case ENUM::S_ALL_SEQ_ID: write(s,  * ); break;
//         case ENUM::S_GI: write(s,  * ); break;
//         case ENUM::S_ALL_GI: write(s,  * ); break;
//         case ENUM::S_ACC: write(s,  * ); break;
//         case ENUM::S_ACCVER: write(s,  * ); break;
//         case ENUM::S_ALLACC: write(s,  * ); break;
        case BlastMatchField<>::Enum::S_LEN:
            write(s, match.sLength);
            break;
        case BlastMatchField<>::Enum::Q_START:
        {
            TPos effectiveQStart    = match.qStart;
            TPos effectiveQEnd      = match.qEnd;
            _untranslateQPositions(effectiveQStart, effectiveQEnd, match.qFrameShift, match.qLength,
                                   context.blastProgram, BlastProgramTag<p>());
            write(s, effectiveQStart);
        } break;
        case BlastMatchField<>::Enum::Q_END:
        {
            TPos effectiveQStart    = match.qStart;
            TPos effectiveQEnd      = match.qEnd;
            _untranslateQPositions(effectiveQStart, effectiveQEnd, match.qFrameShift, match.qLength,
                                   context.blastProgram, BlastProgramTag<p>());
            write(s, effectiveQEnd);
        } break;
        case BlastMatchField<>::Enum::S_START:
        {
            TPos effectiveSStart    = match.sStart;
            TPos effectiveSEnd      = match.sEnd;
            _untranslateSPositions(effectiveSStart, effectiveSEnd, match.sFrameShift, match.sLength,
                                   context.blastProgram, BlastProgramTag<p>());
            write(s, effectiveSStart);
        } break;
        case BlastMatchField<>::Enum::S_END:
        {
            TPos effectiveSStart    = match.sStart;
            TPos effectiveSEnd      = match.sEnd;
            _untranslateSPositions(effectiveSStart, effectiveSEnd, match.sFrameShift, match.sLength,
                                   context.blastProgram, BlastProgramTag<p>());
            write(s, effectiveSEnd);
        } break;
//         case ENUM::Q_SEQ: write(s,  * ); break;
//         case ENUM::S_SEQ: write(s,  * ); break;
        case BlastMatchField<>::Enum::E_VALUE:
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
        case BlastMatchField<>::Enum::BIT_SCORE:
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
        case BlastMatchField<>::Enum::SCORE:
            write(s, match.score);
            break;
        case BlastMatchField<>::Enum::LENGTH:
            write(s, match.aliLength);
            break;
        case BlastMatchField<>::Enum::P_IDENT:
            write(s, FormattedNumber<double>("%.2f",
                                             double(match.identities) * 100
                                             / match.aliLength));
            break;
        case BlastMatchField<>::Enum::N_IDENT:
            write(s, match.identities);
            break;
        case BlastMatchField<>::Enum::MISMATCH:
            write(s, match.mismatches);
            break;
        case BlastMatchField<>::Enum::POSITIVE:
            write(s, match.positives);
            break;
        case BlastMatchField<>::Enum::GAP_OPEN:
            write(s, match.gapOpenings);
            break;
        case BlastMatchField<>::Enum::GAPS:
            write(s, match.gaps);
            break;
        case BlastMatchField<>::Enum::P_POS:
            write(s, FormattedNumber<double>("%.2f",
                                             double(match.positives) * 100
                                             / match.aliLength));
            break;
        case BlastMatchField<>::Enum::FRAMES:
            // for formats that don't have frames, blast says 0 instead of +1
            if (qNumFrames(BlastTabular()) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.qFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            write(s, '/');
            if (sNumFrames(BlastTabular()) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.qFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            break;
        case BlastMatchField<>::Enum::Q_FRAME:
            if (qNumFrames(BlastTabular()) > 1)
                write(s, FormattedNumber<int8_t>("%i", match.qFrameShift));
            else
                write(s, FormattedNumber<int8_t>("%i", 0));
            break;
        case BlastMatchField<>::Enum::S_FRAME:
            if (sNumFrames(BlastTabular()) > 1)
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
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastTabularSpec h>
inline void
_writeFields(TFwdIterator & stream,
             BlastIOContext<TScore, TConString, p, h> & context,
             BlastMatch<TQId, TSId, TPos, TAlign> const & match,
             BlastTabular const &,
             bool const appendNewline = true)
{
    if (SEQAN_LIKELY(!context.legacyFormat))
    {
        for (auto it = seqan::begin(context.fields), itB = it, itEnd = seqan::end(context.fields); it != itEnd; ++it)
        {
            if (it != itB)
                write(stream, '\t');

            _writeField(stream, match, *it, BlastTabular());
        }
    }
    else
    {
        for (auto it = seqan::begin(BlastMatchField<>::defaults), itB = it,
             itEnd = seqan::end(BlastMatchField<>::defaults); it != itEnd; ++it)
        {
            if (it != itB)
                write(stream, '\t');

            _writeField(stream, match, *it, BlastTabular());
        }

        if (SEQAN_ENABLE_DEBUG &&
            ((length(context.fields) != 1) || (context.fields[0] != BlastMatchField<>::Enum::STD)))
            std::cerr << "Attention: custom fields set, but will be ignored, because legacyFormat is also set.\n";

    }

    if (appendNewline)
        write(stream, '\n');
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
 * @param[in,out] context    A @link BlastIOContext<TScore, TConString, p, h> @endlink with parameters and buffers.
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
 * @see BlastRecord#writeRecordHeader
 */

template <typename TQId,
          typename TSId,
          typename TFwdIterator,
          typename TScore,
          typename TConString,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastTabularSpec h>
inline void
writeMatch(TFwdIterator & stream,
           BlastIOContext<TScore, TConString, p, h> & context,
           BlastMatch<TQId, TSId, TPos, TAlign> const & match,
           BlastTabular const & /*tag*/)
{
    _writeFields(stream, context, match, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function _writeFields() or labels [no match object given]
// ----------------------------------------------------------------------------

// template <typename TFwdIterator, BlastFormatFile f, BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// _writeFields(TFwdIterator & /**/,
//              BlastTabular const & /*tag*/)
// {
// }
//
// template <typename TFwdIterator,
//           typename TScore,
//           typename TConString, typename TField, typename... TFields,
//           BlastFormatFile f, BlastFormatProgram p, BlastFormatGeneration g>
// inline void
// _writeFields(TFwdIterator & stream,
//              BlastTabular const & /*tag*/,
//              TField const & field1, TFields const & ... fields)
// {
// //
//     write(stream, '\t');//_seperatorString(BlastTabular()));
//     write(stream, field1);
//     _writeFields(stream, BlastTabular(), fields... );
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
 * with @link BlastFormat#writeRecordHeader @endlink. Please note that this function
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
 * @see BlastRecord#writeRecordHeader
 * @see BlastFormat#writeRecordHeader
 */

// Function for arbitrary number and typed fields
// template <typename TFwdIterator,
//           typename TScore,
//           typename TConString, typename TField, typename... TFields,
//           BlastFormatProgram p, BlastFormatGeneration g>
// inline void
// writeMatch(TFwdIterator & stream,
//            BlastTabular const & /*tag*/,
//             TField const & field1,
//             TFields const & ... fields)
// {
// //     write(stream, field1);
//
//     _writeFields(stream, BlastTabular(), fields...);
//     write(stream, '\n');
// }
//
// // writeMatch TABULAR_WITH_HEADER equal to TABULAR
// template <typename TFwdIterator,
//           typename TScore,
//           typename TConString,
//           typename TField,
//           typename... TFields,
//           BlastFormatProgram p, BlastFormatGeneration g>
// inline void
// writeMatch(TFwdIterator & stream,
//            BlastTabular const & /**/,
//            TField const & field1,
//            TFields const & ... fields)
// {
// //     writeMatch(stream, BlastTabular(), field1, fields... );
// }

// ----------------------------------------------------------------------------
// Function writeTop()
// ----------------------------------------------------------------------------

// TABULAR formats have no special "top"

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabular#writeRecord
 * @headerfile seqan/blast.h
 * @brief write a @link BlastRecord @endlink including it's @link BlastMatch @endlinkes and possible headers to a file.
 * @signature void writeRecord(stream, context, blastRecord, blastTabular);
 *
 * @param[in,out] stream       The file to write to (FILE, fstream, @link OutputStreamConcept @endlink ...)
 * @param[in,out] context      A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastRecord  The @link BlastRecord @endlink you wish to print.
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 *
 * See @link BlastIOContext @endlink for ways to influence the output.
 */

template <typename TFwdIterator,
          typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastTabularSpec h>
inline void
writeRecord(TFwdIterator & stream,
            BlastIOContext<TScore, TConString, p, h> & context,
            BlastRecord<TQId, TSId, TPos, TAlign> const & r,
            BlastTabular const & /*tag*/)
{
    //TODO if debug, do lots of sanity checks on record

    //NOOP for TABULAR
    writeRecordHeader(stream, context, r, BlastTabular());
    for (auto it = r.matches.begin(); it != r.matches.end(); ++it)
    {
        //SOME SANITY CHECKS
        SEQAN_ASSERT(startsWith(r.qId, it->qId));

        writeMatch(stream, context, *it, BlastTabular());
    }
}

/*!
 * @fn BlastTabularOut#writeRecord
 * @headerfile seqan/blast.h
 * @brief write a @link BlastRecord @endlink including it's @link BlastMatch @endlinkes and possible headers to a file.
 * @signature void writeRecord(blastTabularOut, blastRecord);
 *
 * @param[in,out] blastTabularOut A @link BlastTabularOut @endlink formattedFile.
 * @param[in]     blastRecord     The @link BlastRecord @endlink you wish to print.
 *
 * See @link BlastIOContext @endlink for ways to influence the output.
 */

template <typename TScore,
          typename TConString,
          typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastTabularSpec h>
inline void
writeRecord(BlastTabularOut<TScore, TConString, p, h> & formattedFile,
            BlastRecord<TQId, TSId, TPos, TAlign> const & r)
{
    writeRecord(formattedFile.iter, context(formattedFile), r, BlastTabular());
}


// ----------------------------------------------------------------------------
// Function writeBottom()
// ----------------------------------------------------------------------------


// template <typename TFwdIterator,
//           typename TScore,
//           typename TConString,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// writeBottom(TFwdIterator & /**/,
//             TDbSpecs const & /**/,
//             TBlastScoringAdapater const & /**/,
//             BlastTabular const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
// // }
// 
// template <typename TFwdIterator,
//           typename TScore,
//           typename TConString,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// writeBottom(TFwdIterator & /**/,
//             TDbSpecs const & /**/,
//             TBlastScoringAdapater const & /**/,
//             BlastTabular const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
// // }

} // namespace seqan
#endif // header guard
