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

#ifndef SEQAN_BLAST_BLAST_TABULAR_READ_H_
#define SEQAN_BLAST_BLAST_TABULAR_READ_H_

#include <regex>

/* IMPLEMENTATION NOTES

BLAST TABULAR example:

The format of a blast tabular output file is less simple than it looks, here's
the general form

HEADER
 MATCH
 MATCH
 MATCH
 MATCH
HEADER
 MATCH
HEADER
HEADER
...

=> Header for each sequence, 0-n Matchs for each sequence
=> Each record is one-line, each Header is multiline

A Header usually consists of:

# Program Tag [e.g BLASTX 2.2.27+ ]
# Query Id [ID of the query *sequence*]
# Database Id [note that this is not the name of the sequence in the db, but of
  the database itself, e.g. "nr" -> usually the same for each file]
# Fields: [Headers of columns]
# n "hits found"

The first three lines are always written.
The Fields line is always writen by NCBI Blast, but only when hits > 0 by NCBI Blast+.
The "number of hits"-line is always printed by NCBI Blast+, and never by NCBI Blast.

Possibly other lines can be written as comments.

Because 0 matches are allowed, multiple Headers can succeed each other, the
criterium for seperation employed by this implementation is that an NCBI Blast
headers always ends after the "Fields" line and NCBI Blast+ records end after
the "number of hits"-line.
*/

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Type BlastFileIn
// ----------------------------------------------------------------------------

/*!
 * @class BlastTabularIn
 * @signature template <typename TBlastIOContext>
 * using BlastTabularIn = FormattedFile<BlastTabular, Input, TBlastIOContext>;
 * @extends FormattedFileIn
 * @headerfile <seqan/blast.h>
 * @brief FormattedFileIn abstraction for @link BlastTabular @endlink
 *
 * This is a @link FormattedFile @endlink specialization for reading @link BlastTabular @endlink formats. For details
 * on how to influence the reading of files and how to differentiate between the tabular format without headers and the
 * one with headers, see @link BlastIOContext @endlink. TODO change
 * Please note that you have specify the type of the context as a template parameter to BlastTabularIn, see the example
 * below.
 *
 * TODO example
 *
 * @see BlastRecord
 */

template <typename TBlastIOContext = BlastIOContext<>>
using BlastTabularIn = FormattedFile<BlastTabular, Input, TBlastIOContext>;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

//TODO document?
template <typename TScore,
          typename TString,
          typename TFwdIterator,
          BlastProgram p,
          BlastTabularSpec h>
inline bool
atEnd(BlastIOContext<TScore, TString, p, h> const & context,
      TFwdIterator const & iter,
      BlastTabular const &)
{
    return (atEnd(iter) && empty(context._lineBuffer));
}

template <typename TContext>
inline bool
atEnd(BlastTabularIn<TContext> & formattedFile)
{
    return atEnd(context(formattedFile), formattedFile.iter, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function guessFormat()
// ----------------------------------------------------------------------------

template <typename TSpec>
inline bool guessFormat(FormattedFile<BlastTabular, Input, TSpec> & file)
{
    if (value(file.iter) == '#')
        setBlastTabularSpec(context(file), BlastTabularSpec::HEADER);
    else
        setBlastTabularSpec(context(file), BlastTabularSpec::NO_HEADER);
    return true;
}

// ----------------------------------------------------------------------------
// Function strSplit()
// ----------------------------------------------------------------------------

template <typename TString,
          typename TSpec,
          typename TSequence>
inline void
_strSplit(StringSet<TString, TSpec> & result,
         TSequence const & sequence,
         std::regex const & regex)
{
    static thread_local std::string tmp;
    tmp = std::regex_replace(sequence, regex, "\x7F");
    strSplit(result, tmp, EqualsChar<'\x7F'>(), true);
}

// ----------------------------------------------------------------------------
// Function onMatch()
// ----------------------------------------------------------------------------

//NOTE(h-2): dox disabled to clean-up interface
/*
 * @fn BlastTabular#onMatch
 * @brief Returns whether the iterator is on the beginning of a match line.
 *
 * @signature bool onMatch(stream, blastTabular)
 *
 * @param[in] iter    An input iterator over a stream or any fwd-iterator over a string
 * @param[in] blastTabular     @link BlastTabular @endlink specialization
 *
 * @throw IOError On low-level I/O errors.
 *
 * @return    true or false
 */

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline bool
onMatch(BlastIOContext<TScore, TString, p, h> & context,
        BlastTabular const &)
{
    return (!startsWith(context._lineBuffer, "#") && !empty(context._lineBuffer));
}

// ----------------------------------------------------------------------------
// Function onRecord()
// ----------------------------------------------------------------------------

//TODO document publicly
/*!
 * @fn BlastTabular#onRecord
 * @brief Returns whether the iterator is on the beginning of a match line.
 *
 * @signature bool onRecord(context, blastTabular)
 *
 * @param[in] iter    An input iterator over a stream or any fwd-iterator over a string
 * @param[in] blastTabular     @link BlastTabular @endlink specialization
 *
 * @throw IOError On low-level I/O errors.
 *
 * @return    true or false
 */

template <typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline bool
onRecord(BlastIOContext<TScore, TString, p, h> & context,
         BlastTabular const &)
{
    if (getBlastTabularSpec(context) == BlastTabularSpec::NO_HEADER)
        return onMatch(context, BlastTabular());

    //      RecordHeader                                  Footer
    return (startsWith(context._lineBuffer, "#") && (!startsWith(context._lineBuffer, "# BLAST processed ")));
}

// ----------------------------------------------------------------------------
// Function goNextLine()
// ----------------------------------------------------------------------------

template <typename TScore,
          typename TString,
          typename TFwdIterator,
          BlastProgram p,
          BlastTabularSpec h>
inline void
goNextLine(BlastIOContext<TScore, TString, p, h> & context,
           TFwdIterator & iter,
           BlastTabular const &)
{
    clear(context._lineBuffer);
    if (!atEnd(iter))
        readLine(context._lineBuffer, iter);
}

// ----------------------------------------------------------------------------
// Function readRecordHeader()
// ----------------------------------------------------------------------------

//NOTE(h-2): dox disabled to clean-up interface
/*
 * @fn BlastTabular#readRecordHeader
 * @brief read a the header of a record from a Blast tabular file
 * @headerfile seqan/blast.h
 *
 * @signature void readRecordHeader(blastRecord, stream, context, blastTabular);
 *
 * @param[out]    blastRecord  A @link BlastRecord @endlink.
 * @param[in,out] stream       An input iterator over a stream or any fwd-iterator over a string.
 * @param[in,out] context      An @link BlastIOContext @endlink with buffers.
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 *
 * @section Remarks
 *
 * Call this function on every line beginning that is not "onMatch". Will
 * set the @link BlastRecord::qId @endlink member of the record and resize the @link BlastRecord::matches @endlink
 * member to the expected number of matches succeeding the header. [In case the format is detected as being legacy
 * (see below), the matches member will not be resized and you need to append matches instead of assigning them if
 * using @link BlastTabular#readMatch @endlink ].
 *
 * This function also sets many properties of the @link BlastIOContext @endlink, including these members:
 * <li> @link BlastIOContext::versionString @endlink: version string of the header</li>
 * <li> @link BlastIOContext::dbName @endlink: name of the database</li>
 * <li> @link BlastIOContext::fields @endlink: descriptors for the columns.</li>
 * <li> @link BlastIOContext::fieldsAsStrings @endlink: labels of the columns
 * as they appear in the file.</li>
 * <li> @link BlastIOContext::conformancyErrors @endlink: if this StringSet is not empty, then there are issues in the
 * header.</li>
 * <li> @link BlastIOContext::otherLines @endlink: any lines that cannot be interpreted; these always also imply
 * conformancyErrors.</li>
 * <li>@link BlastIOContext::legacyFormat @endlink: whether the record header (and likely the entire file) is in
 * legacyFormat.</li>
 *
 * It also sets the blast program run-time parameter of the context depending on the information found in the header. If
 * the compile time parameter was set on the context and they are different this will result in a conformancyError.
 * Otherwise you can read the value with @link BlastIOContext#getBlastProgram @endlink.
 *
 * Please note that for @link BlastIOContext::legacyFormat @endlink the @link BlastIOContext::fields @endlink member
 * is always ignored, however @link BlastIOContext::fieldsAsStrings @endlink is still read from the header, in case
 * you want to process it.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFwdIterator,
          typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readRecordHeaderImpl(BlastRecord<TQId, TSId, TPos, TAlign> & r,
                      TFwdIterator & iter,
                      BlastIOContext<TScore, TString, p, h> & context,
                      BlastTabular const &)
{
    // this is a match instead of a header
    if (onMatch(context, BlastTabular()))
        SEQAN_THROW(ParseError("Header expected, but no header found."));
    else
        setBlastTabularSpec(context, BlastTabularSpec::HEADER);

    int queryLinePresent = 0;
    int dbLinePresent = 0;
    int fieldsLinePresent = 0;
    int hitsLinePresent = 0;

    // TODO replace the startsWith and suffix() with std::regex
    do
    {
        assign(context._lineBuffer, suffix(context._lineBuffer, 2)); // strip "# "

        if (startsWith(context._lineBuffer, "BLAST") || startsWith(context._lineBuffer, "TBLAST"))
        {

            // last line of file
            if (SEQAN_UNLIKELY(startsWith(context._lineBuffer, "# BLAST processed ") && !context.legacyFormat))
            {
                //TODO we landed in the footer by mistake
            }
            else // first line of record header
            {
                std::swap(context.versionString, context._lineBuffer);
                context.blastProgram = _programStringToTag(prefix(context.versionString,
                                                                  _firstOcc(context.versionString, ' ')));

    // TODO(h4nn3s): regex test
    //             std::regex versionRE("[0-9]+\\.[0-9]+\\.[0-9]+\\+", std::regex::awk);
    //             context.legacyFormat = std::regex_search(seqan::begin(key), seqan::end(key), versionRE);

                context.legacyFormat = (context.versionString.find('+') == std::string::npos);
            }
        }
        else if (startsWith(context._lineBuffer, "Query:"))
        {
            ++queryLinePresent;
            assign(r.qId, suffix(context._lineBuffer, 6));

        }
        else if (startsWith(context._lineBuffer, "Database:"))
        {
            ++dbLinePresent;
            assign(context.dbName, suffix(context._lineBuffer, 9));
        }
        else if (startsWith(context._lineBuffer, "Fields:"))
        {
            ++fieldsLinePresent;
            _strSplit(context.fieldsAsStrings,
                      static_cast<TString>(suffix(context._lineBuffer, 7)),
                      std::regex(", "));

            if (context.legacyFormat)
            {
                // assume defaults for LEGACY
                if (!context.ignoreFieldsInHeader)
                    appendValue(context.fields, BlastMatchField<>::Enum::STD);

//                 break; // header is finished
            }
            else if (!context.ignoreFieldsInHeader)
            {
                bool defaults = true;
                resize(context.fields, length(context.fieldsAsStrings));
                for (uint8_t j = 0; j < length(context.fieldsAsStrings); ++j)
                {
                    for (uint8_t i = 0; i < length(BlastMatchField<>::columnLabels); ++i)
                    {
                        if (context.fieldsAsStrings[j] ==
                            BlastMatchField<>::columnLabels[i])
                        {
                            context.fields[j] =
                                static_cast<BlastMatchField<>::Enum>(i);

                            if ((j >= length(BlastMatchField<>::defaults)) &&
                                (static_cast<BlastMatchField<>::Enum>(i) !=
                                BlastMatchField<>::defaults[j]))
                                defaults = false;
                            break;
                        }
                    }
                }
                if (defaults) // replace multiple fields with the default meta field
                {
                    clear(context.fields);
                    appendValue(context.fields, BlastMatchField<>::Enum::STD);
                }
            }
        }
        else
        {
            if (!context.legacyFormat)
            {
                // is hits counter?
                if (endsWith(context._lineBuffer, "hits found"))
                {
                    clear(context._stringBuffer);
                    for (unsigned i = 0; (i < length(context._lineBuffer) && isdigit(context._lineBuffer[i])); ++i)
                        appendValue(context._stringBuffer, context._lineBuffer[i], Generous());

                    uint64_t hits = lexicalCast<uint64_t>(context._stringBuffer);

                    if (hits)
                    {
                        resize(r.matches, hits);
                    }
                    else  // hits = 0 means no fieldList, restore default
                    {
                        appendValue(context.fields, BlastMatchField<>::Enum::STD);
                        _strSplit(context.fieldsAsStrings, BlastMatchField<>::columnLabels[0], std::regex(", "));
                    }

                    ++hitsLinePresent;
//                     break; // header is finished
                }
                else
                {
                    appendValue(context.otherLines, context._lineBuffer, Generous());
                }
            }
            else
            {
                appendValue(context.otherLines, context._lineBuffer, Generous());
            }
        }

        goNextLine(context, iter, BlastTabular());

    } while (startsWith(context._lineBuffer, "#") && !startsWith(context._lineBuffer, "# Blast")); // next header

    if (getBlastProgram(context) == BlastProgram::UNKNOWN)
        appendValue(context.conformancyErrors,
                    "Type of BlastProgram could not be determined from header, you are advised to look "
                    "at context.versionString and context.otherLines.");
    else if ((p != BlastProgram::UNKNOWN) && (context.blastProgram != p))
        appendValue(context.conformancyErrors,
                    std::string("You fixed the BlastProgramType to ") +
                    std::string(_programTagToString(p)) +
                    std::string (" at compile-time, but the type ") +
                    std::string(_programTagToString(getBlastProgram(context))) +
                    std::string(" was detected in the file!"));

    if (queryLinePresent != 1)
        appendValue(context.conformancyErrors, "No or multiple query lines present.");

    if (dbLinePresent != 1)
        appendValue(context.conformancyErrors, "No or multiple database lines present.");

    if (context.legacyFormat)
    {
        if (fieldsLinePresent != 1)
            appendValue(context.conformancyErrors, "No or multiple fields lines present.");
    }
    else
    {
        // is ommitted in BLAST_PLUS when there are no hits
        if ((fieldsLinePresent != 1) && (length(r.matches) > 0))
            appendValue(context.conformancyErrors, "No or multiple fields lines present.");
    }

    if (length(context.otherLines) != 0)
        appendValue(context.conformancyErrors, "Unexpected lines present, see context.otherLines.");

}

template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFwdIterator,
          typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
readRecordHeader(BlastRecord<TQId, TSId, TPos, TAlign> & r,
                 TFwdIterator & iter,
                 BlastIOContext<TScore, TString, p, h> & context,
                 BlastTabular const &)
{
    ++context.numberOfRecords;

    if (getBlastTabularSpec(context) == BlastTabularSpec::NO_HEADER)
        return;

    clear(r);
    clear(context.versionString);
    clear(context.dbName);
    clear(context.otherLines);
    clear(context.conformancyErrors);
    clear(context.fieldsAsStrings);

    if (!context.ignoreFieldsInHeader)
        clear(context.fields);

    _readRecordHeaderImpl(r, iter, context, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function skipHeader()
// ----------------------------------------------------------------------------

//NOTE(h-2): dox disabled to clean-up interface
/*
 * @fn BlastTabular#skipRecordHeader
 * @brief skip a header from a Blast tabular input file.
 *
 * @signature void skipRecordHeader(stream, context, blastTabular);
 *
 * @param[in,out] stream       An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context      A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 *
 * @section Remarks
 *
 * Call this function whenever you want to skip exactly one header (and possibly examine the members of the
 * context). If you want
 * to go directly to the beginning of the next match (possibly skipping multiple
 * headers that have no succeeding matches) use @link BlastTabular#skipUntilMatch @endlink instead.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @see BlastTabular#skipUntilMatch
 * @headerfile seqan/blast.h
 */

template <typename TFwdIterator,
          typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
skipRecordHeader(TFwdIterator & iter,
                 BlastIOContext<TScore, TString, p, h> & context,
                 BlastTabular const & /*tag*/)
{
    readRecordHeader(context.bufRecord, iter, context, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function readMatch()
// ----------------------------------------------------------------------------

template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readField(BlastMatch<TQId, TSId, TPos, TAlign> & match,
           BlastIOContext<TScore, TString, p, h> & context,
           typename BlastMatchField<>::Enum const fieldId)
{
    switch (fieldId)
    {
        case BlastMatchField<>::Enum::STD: // this is cought in the calling function
            break;
        case BlastMatchField<>::Enum::Q_SEQ_ID:
            match.qId = context._stringBuffer;
            break;
//         case ENUM::Q_GI: write(s,  * ); break;
//         case ENUM::Q_ACC: write(s,  * ); break;
//         case ENUM::Q_ACCVER: write(s,  * ); break;
        case BlastMatchField<>::Enum::Q_LEN:
            match.qLength = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::S_SEQ_ID:
            match.sId = context._stringBuffer;
            break;
//         case ENUM::S_ALL_SEQ_ID: write(s,  * ); break;
//         case ENUM::S_GI: write(s,  * ); break;
//         case ENUM::S_ALL_GI: write(s,  * ); break;
//         case ENUM::S_ACC: write(s,  * ); break;
//         case ENUM::S_ACCVER: write(s,  * ); break;
//         case ENUM::S_ALLACC: write(s,  * ); break;
        case BlastMatchField<>::Enum::S_LEN:
            match.sLength = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::Q_START:
            match.qStart = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::Q_END:
            match.qEnd = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::S_START:
            match.sStart = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::S_END:
            match.sEnd = lexicalCast<TPos>(context._stringBuffer);
            break;
//         case ENUM::Q_SEQ: write(s,  * ); break;
//         case ENUM::S_SEQ: write(s,  * ); break;
        case BlastMatchField<>::Enum::E_VALUE:
            match.eValue = lexicalCast<double>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::BIT_SCORE:
            match.bitScore = lexicalCast<double>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::SCORE:
            match.alignStats.alignmentScore = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::LENGTH:
            match.alignStats.alignmentLength = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::P_IDENT:
            match.alignStats.alignmentIdentity = lexicalCast<double>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::N_IDENT:
            match.alignStats.numMatches = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::MISMATCH:
            match.alignStats.numMismatches = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::POSITIVE:
            match.alignStats.numPositiveScores = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::GAP_OPEN:
            match.alignStats.numGapOpens = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::GAPS:
            match.alignStats.numGaps = lexicalCast<TPos>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::P_POS:
            match.alignStats.alignmentSimilarity = lexicalCast<double>(context._stringBuffer);
        case BlastMatchField<>::Enum::FRAMES:
        {
            clear(context._setBuffer2);
            strSplit(context._setBuffer2, context._stringBuffer, EqualsChar<'/'>());
            if (length(context._setBuffer2) != 2)
                SEQAN_THROW(ParseError("Could not process frame string."));
            match.qFrameShift = lexicalCast<int8_t>(context._setBuffer2[0]);
            match.sFrameShift = lexicalCast<int8_t>(context._setBuffer2[1]);
        } break;
        case BlastMatchField<>::Enum::Q_FRAME:
            match.qFrameShift = lexicalCast<int8_t>(context._stringBuffer);
            break;
        case BlastMatchField<>::Enum::S_FRAME:
            match.sFrameShift = lexicalCast<int8_t>(context._stringBuffer);
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
            SEQAN_THROW(ParseError("The requested column type is not yet "
                                   "implemented."));
    };
}

/*
 * @fn BlastTabular#readMatch
 * @brief read a match from a file in BlastTabular format
 *
 * @signature void readMatch(blastMatch, stream, blastTabular);
 *
 * @param[out]    blastMatch   A @link BlastMatch @endlink object to hold all relevant info
 * @param[in,out] stream       An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context      A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 *
 * @section Remarks
 *
 * The @link BlastIOContext::fields @endlink member of the context
 * can be specified if you expect a custom column composition. Specifying this
 * parameter implies that you know the columns are not
 * default and that they instead represent the given order and types.
 * You may specify less
 * fields than are actually present, in this case the additional fields will be
 * discarded. The parameter is ignored if @link BlastIOContext::legacyFormat @endlink is set.
 *
 * To differentiate between members of a @link BlastMatch @endlink that were read from the file and those that have
 * not been set, the latter are initialized to their respective max-values.
 *
 * Please note that the only transformations made to the data are the following:
 *
 *  <li> computation of the number of identities (from the percentage) [default]</li>
 *  <li> computation of the number of positives (from the percentage) [if given]</li>
 *  <li> number of gaps computed from other values [default]</li>
 *
 * In contrast to @link BlastTabular#writeMatch @endlink no other transformations
 * are made, e.g. the positions are still one-indexed and
 * flipped for reverse strand matches. This is due to the required fields for
 * retransformation (sequence lengths, frames) not being available in the
 * default columns.
 *
 * Instead of using this signature you may also use @link BlastTabular#readMatch0 @endlink which works without a
 * @link BlastMatch @endlink parameter and supports arbitrary columns.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @headerfile seqan/blast.h
 */

template <typename TQId,
          typename TSId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
readMatch(BlastMatch<TQId, TSId, TPos, TAlign> & match,
          TFwdIterator & iter,
          BlastIOContext<TScore, TString, p, h> & context,
          BlastTabular const &)
{
    if (context.legacyFormat)
    {
        if (SEQAN_ENABLE_DEBUG)
        {
            if ((length(context.fields) != 1) || (context.fields[0] != BlastMatchField<>::Enum::STD))
                std::cerr << "Warning: custom fields set, but will be reset, because legacyFormat is also set.\n";
        }

        // set defaults
        clear(context.fields);
        appendValue(context.fields,  BlastMatchField<>::Enum::STD);
    }

    // header should have been read or skipped
    if (SEQAN_UNLIKELY(!onMatch(context, BlastTabular())))
        SEQAN_THROW(ParseError("Not on beginning of Match (you should have skipped comments)."));

    match._maxInitialize(); // mark all members as "not set"

    clear(context._setBuffer1);

    auto & fields = context._setBuffer1;
    strSplit(fields, context._lineBuffer, IsTab());

    // line in buffer was processed so new line is read from iter into buffer
    goNextLine(context, iter, BlastTabular());

    unsigned n = 0u;
    for (BlastMatchField<>::Enum const f : context.fields)
    {
        // this field represents multiple fields
        if (f == BlastMatchField<>::Enum::STD)
        {
            for (typename BlastMatchField<>::Enum const f2 : BlastMatchField<>::defaults)
            {
                if (SEQAN_UNLIKELY(n >= length(fields)))
                    SEQAN_THROW(ParseError("More columns expected than were present in file."));

                context._stringBuffer = static_cast<decltype(context._stringBuffer)>(fields[n++]);
                _readField(match, context, f2);
            }
        } else
        {
            if (SEQAN_UNLIKELY(n >= length(fields)))
                SEQAN_THROW(ParseError("More columns expected than were present in file."));

            context._stringBuffer = static_cast<decltype(context._stringBuffer)>(fields[n++]);
            _readField(match, context, f);
        }
    }

    // retransform the percentages to real numbers
    if (!_memberIsSet(match.alignStats.numMatches) &&
        _memberIsSet(match.alignStats.alignmentLength) &&
        _memberIsSet(match.alignStats.alignmentIdentity))
        match.alignStats.numMatches = std::lround(double(match.alignStats.alignmentLength) *
                                            match.alignStats.alignmentIdentity / 100.0);
    if (!_memberIsSet(match.alignStats.numPositiveScores) &&
        _memberIsSet(match.alignStats.alignmentLength) &&
        _memberIsSet(match.alignStats.alignmentSimilarity))
        match.alignStats.numPositiveScores = std::lround(double(match.alignStats.alignmentLength) *
                                            match.alignStats.alignmentSimilarity / 100.0);

    //TODO also transform numbers to percentages; possibly do other calculations

    // and compute gaps from other values
    // since gaps are included in the mismatch count in BLAST_LEGACY they cannot be computed here
    if (!context.legacyFormat &&
        !_memberIsSet(match.alignStats.numGaps) &&
        _memberIsSet(match.alignStats.alignmentLength) &&
        _memberIsSet(match.alignStats.numMatches) &&
        _memberIsSet(match.alignStats.numMismatches))
        match.alignStats.numGaps = match.alignStats.alignmentLength -
                                   match.alignStats.numMatches -
                                   match.alignStats.numMismatches;
}

// ----------------------------------------------------------------------------
// Function skipMatch()
// ----------------------------------------------------------------------------

//NOTE(h-2): dox disabled to clean-up interface
/*
 * @fn BlastTabular#skipMatch
 * @brief skip a line that contains a match
 * @signature void skipMatch(stream, context, blastTabular);
 * @headerfile seqan/blast.h
 *
 * @param[in,out] stream  An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 *
 * @section Remarks
 *
 * This function always checks whether it is on a line that is not a comment and
 * it does verify the line read as looking like a match. If you do not want this
 * behaviour you can instead <tt>skipLine(iter)</tt> which is also
 * faster.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TFwdIterator,
          typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
skipMatch(TFwdIterator & iter,
          BlastIOContext<TScore, TString, p, h> & context,
          BlastTabular const &)
{
    readMatch(context.bufMatch, iter, context, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TAlign,
          typename TPos,
          typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readRecordWithHeader(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
                      TFwdIterator & iter,
                      BlastIOContext<TScore, TString, p, h> & context,
                      BlastTabular const &)
{
    readRecordHeader(blastRecord, iter, context, BlastTabular());

    if (!context.legacyFormat) // this is detected from the header
    {
        // .matches already resized for us
        for (auto & m : blastRecord.matches)
        {
            if (!onMatch(context, BlastTabular()))
            {
                appendValue(context.conformancyErrors,
                            "Less matches than promised by header");
                break;
            }

            readMatch(m, iter, context, BlastTabular());
        }

        if ((!atEnd(context, iter, BlastTabular())) && onMatch(context, BlastTabular()))
            appendValue(context.conformancyErrors,
                        "More matches than promised by header");

        while ((!atEnd(context, iter, BlastTabular())) && onMatch(context, BlastTabular()))
        {
            blastRecord.matches.emplace_back();
            readMatch(back(blastRecord.matches), iter, context, BlastTabular());
        }
    } else
    {
        while ((!atEnd(context, iter, BlastTabular())) && onMatch(context, BlastTabular()))
        {
            blastRecord.matches.emplace_back();
            readMatch(back(blastRecord.matches), iter, context, BlastTabular());
        }
    }
}

// TABULAR WITHOUT HEADER
template <typename TQId,
          typename TSId,
          typename TAlign,
          typename TPos,
          typename TFwdIterator,
          typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readRecordNoHeader(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
                    TFwdIterator & iter,
                    BlastIOContext<TScore, TString, p, h> & context,
                    BlastTabular const &)
{
    ++context.numberOfRecords;

    clear(blastRecord);

    auto it = seqan::begin(context._lineBuffer); // move into line below when seqan supports && properly
    readUntil(blastRecord.qId, it, IsTab());

    TString curIdPlusTab = blastRecord.qId;
    appendValue(curIdPlusTab, '\t');

    std::vector<typename BlastMatchField<>::Enum> fieldListMinusFirst;
    if (context.fields[0] == BlastMatchField<>::Enum::STD)
    {
        fieldListMinusFirst = BlastMatchField<>::_defaultsMinusFirst;
        if (length(context.fields) > 1)
            fieldListMinusFirst.insert(std::end(fieldListMinusFirst),
                                       std::begin(context.fields) + 1,
                                       std::end(context.fields));
    }
    else if (context.fields[0] == BlastMatchField<>::Enum::Q_SEQ_ID)
    {
        assign(fieldListMinusFirst, suffix(context.fields, 1));
    } else
    {
        SEQAN_FAIL("readRecord interface on header-less format with custom "
                   "fields not supported, unless first custom field is "
                   "Q_SEQ_ID. Use the readMatch interface instead.");
    }

    // use fieldListMinusFirst
    context.fields.swap(fieldListMinusFirst);

    while ((!atEnd(context, iter, BlastTabular())) && onMatch(context, BlastTabular()))
    {
        blastRecord.matches.emplace_back();
        // read remainder of line
        readMatch(back(blastRecord.matches), iter, context, BlastTabular());

        if (!startsWith(context._lineBuffer, curIdPlusTab)) // next record reached
            break;
    }

    // revert to original behaviour
    context.fields.swap(fieldListMinusFirst);

    if (length(blastRecord.matches) == 0)
        SEQAN_THROW(ParseError("No Matches could be read."));
}

/*!
 * @fn BlastTabular#readRecord
 * @brief Read a record from a file in BlastTabular format.
 * @headerfile seqan/blast.h
 * @signature void readRecord(blastRecord, stream, context, blastTabular);
 *
 * @param[out]    blastRecord  A @link BlastRecord @endlink to hold all information related to one query sequence.
 * @param[in,out] stream       An input iterator over a stream or any fwd-iterator over a string.
 * @param[in,out] context      A @link BlastIOContext @endlink with further parameters and buffers.
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 *
 * @section Remarks
 *
 * This function will read an entire record from a blast tabular file, i.e. it will read the record header
 * (if it exists) and 0-n @link BlastMatch @endlinkes belonging to one query.
 *
 * Please note that if there are no headers in the file the boundary
 * between records is inferred from the indentity of the first field, i.e.
 * non-standard field configurations must also have Q_SEQ_ID as their first @link BlastMatchField @endlink.
 *
 * @subsection Record header
 *
 * The @link BlastRecord::qId @endlink member of the record is read from the header and  the
 * @link BlastRecord::matches @endlink are resized to the expected number of matches succeeding the header.
 *
 * This function also sets many properties of the @link BlastIOContext @endlink, including these members:
 * <li> @link BlastIOContext::versionString @endlink: version string of the header</li>
 * <li> @link BlastIOContext::dbName @endlink: name of the database</li>
 * <li> @link BlastIOContext::fields @endlink: descriptors for the columns.</li>
 * <li> @link BlastIOContext::fieldsAsStrings @endlink: labels of the columns
 * as they appear in the file.</li>
 * <li> @link BlastIOContext::conformancyErrors @endlink: if this StringSet is not empty, then there are issues in the
 * header.</li>
 * <li> @link BlastIOContext::otherLines @endlink: any lines that cannot be interpreted; these always also imply
 * conformancyErrors.</li>
 * <li>@link BlastIOContext::legacyFormat @endlink: whether the record header (and likely the entire file) is in
 * legacyFormat.</li>
 *
 * It also sets the blast program run-time parameter of the context depending on the information found in the header. If
 * the compile time parameter was set on the context and they are different this will result in a conformancyError.
 * Otherwise you can read the value with @link BlastIOContext#getBlastProgram @endlink.
 *
 * Please note that for @link BlastIOContext::legacyFormat @endlink the @link BlastIOContext::fields @endlink member
 * is always ignored, however @link BlastIOContext::fieldsAsStrings @endlink is still read from the header, in case
 * you want to process it.
 *
 * In case you do not wish the @link BlastIOContext::fields @endlink to be read from the header, you can set
 * context.@link BlastIOContext::ignoreFieldsInHeader @endlink to true. This will be prevent it from being read and will
 * allow you to specify it manually which might be relevant for reading the match lines.
 *
 * @subsection Matches
 *
 * A match line contains 1 - n columns or fields, 12 by default.
 * The @link BlastIOContext::fields @endlink member of the context is considered when reading these fields. It is
 * usually extracted from the header but can also be set by yourself if there is no header or if you want to overwrite
 * the headers information (see above).
 * You may specify less
 * fields than are actually present, in this case the additional fields will be
 * discarded. The parameter is ignored if @link BlastIOContext::legacyFormat @endlink is set.
 *
 * To differentiate between members of a @link BlastMatch @endlink that were read from the file and those that have
 * not been set (e.g. both could be 0), the latter are initialized to their respective max-values.
 *
 * Please note that the only transformations made to the data are the following:
 *
 *  <li> computation of the number of identities (from the percentage) [default]</li>
 *  <li> computation of the number of positives (from the percentage) [if given]</li>
 *  <li> number of gaps computed from other values [default]</li>
 *
 * In contrast to @link BlastTabular#writeRecord @endlink no other transformations
 * are made, e.g. the positions are still one-indexed and
 * flipped for reverse strand matches. This is due to the required fields for
 * retransformation (sequence lengths, frames) not being available in the
 * default columns.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TQId,
          typename TSId,
          typename TAlign,
          typename TPos,
          typename TFwdIterator,
          typename TScore,
          typename TString,
          BlastProgram p>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
          TFwdIterator & iter,
          BlastIOContext<TScore, TString, p, BlastTabularSpec::NO_HEADER> & context,
          BlastTabular const &)
{
    _readRecordNoHeader(blastRecord, iter, context, BlastTabular());
}

template <typename TQId,
          typename TSId,
          typename TAlign,
          typename TPos,
          typename TFwdIterator,
          typename TScore,
          typename TString,
          BlastProgram p>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
          TFwdIterator & iter,
          BlastIOContext<TScore, TString, p, BlastTabularSpec::HEADER> & context,
          BlastTabular const &)
{
    _readRecordWithHeader(blastRecord, iter, context, BlastTabular());
}

template <typename TQId,
          typename TSId,
          typename TAlign,
          typename TPos,
          typename TFwdIterator,
          typename TScore,
          typename TString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
          TFwdIterator & iter,
          BlastIOContext<TScore, TString, p, h> & context,
          BlastTabular const &)
{
    if (getBlastTabularSpec(context) == BlastTabularSpec::HEADER)
        _readRecordNoHeader(blastRecord, iter, context, BlastTabular());
    else
        _readRecordWithHeader(blastRecord, iter, context, BlastTabular());
}

/*!
 * @fn BlastTabularIn#readRecord
 * @brief read a record from a blast tabular formattedFile
 * @headerfile seqan/blast.h
 * @signature void readRecord(blastRecord, formattedFile);
 *
 * @param[out]    blastRecord   A @link BlastRecord @endlink to hold all information related to one query sequence.
 * @param[in,out] formattedFile The @link BlastTabularOut @endlink file.
 *
 * See @link BlastTabular#readRecord @endlink for details.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TQId,
          typename TSId,
          typename TAlign,
          typename TPos,
          typename TContext>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
           BlastTabularIn<TContext> & formattedFile)
{
    readRecord(blastRecord, formattedFile.iter, context(formattedFile), BlastTabular());
}

// ----------------------------------------------------------------------------
// Function readHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabular#readHeader
 * @headerfile seqan/blast.h
 * @brief Read the header (top-most section) of a BlastTabular file (this is a NOOP).
 * @signature void readHeader(stream, context, blastTabular);
 *
 * @param[in,out] stream         The file to read to (FILE, fstream, @link InputStreamConcept @endlink ...)
 * @param[in,out] context        A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastTabular   The @link BlastTabular @endlink tag.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TFwdIterator,
          typename TScore,
          typename TConString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
readHeader(BlastIOContext<TScore, TConString, p, h> & context,
           TFwdIterator & iter,
           BlastTabular const & /*tag*/)
{
    readLine(context._lineBuffer, iter); // fill the line-buffer the first time

    if (getBlastTabularSpec(context) == BlastTabularSpec::UNKNOWN)
    {
        if (onMatch(context, BlastTabular()))
            setBlastTabularSpec(context, BlastTabularSpec::NO_HEADER);
        else
            setBlastTabularSpec(context, BlastTabularSpec::HEADER);
    }
}

/*!
 * @fn BlastTabularIn#readHeader
 * @headerfile seqan/blast.h
 * @brief Read the header (top-most section) of a BlastTabular file (this is a NOOP).
 * @signature void readHeader(blastTabularIn);
 *
 * @param[in,out] blastTabularIn A @link BlastTabularIn @endlink formattedFile.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TContext>
inline void
readHeader(BlastTabularIn<TContext> & formattedFile)
{
    readHeader(context(formattedFile), formattedFile.iter, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function readFooter()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabular#readFooter
 * @headerfile seqan/blast.h
 * @brief Read the footer of a BlastTabular file (this is a NOOP).
 * @signature void readFooter(stream, context, blastTabular);
 *
 * @param[in,out] stream         The file to read to (FILE, fstream, @link InputStreamConcept @endlink ...)
 * @param[in,out] context        A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     blastTabular   The @link BlastTabular @endlink tag.
 *
 * @section Remarks
 *
 * While there is a footer in the tabular format, it is not necessary to call this function (it will be read by the
 * last call to readRecord).
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TFwdIterator,
          typename TScore,
          typename TConString,
          BlastProgram p,
          BlastTabularSpec h>
inline void
readFooter(BlastIOContext<TScore, TConString, p, h> & context,
           TFwdIterator & iter,
           BlastTabular const & /*tag*/)
{
    clear(context.otherLines);
    clear(context.conformancyErrors);
    clear(context.fieldsAsStrings);
    clear(context.fields);

    if (getBlastTabularSpec(context) == BlastTabularSpec::HEADER)
    {
        if (SEQAN_UNLIKELY(!startsWith(context._lineBuffer, "# BLAST processed")))
            SEQAN_FAIL("ERROR: Tried to read footer, but was not on footer.");

        clear(context._stringBuffer);
        auto it = seqan::begin(context._lineBuffer);
        it += 16; // skip "BLAST processed "
        readUntil(context._stringBuffer, it,  IsBlank());

        uint64_t numRecords = lexicalCast<uint64_t>(context._stringBuffer);

        clear(context.conformancyErrors);
        if (context.numberOfRecords < numRecords)
            appendValue(context.conformancyErrors, "The file claims to contain more records than you read.");
        else if (context.numberOfRecords > numRecords)
            appendValue(context.conformancyErrors, "The file claims to contain less records than you read.");

        if (!atEnd(iter))
            appendValue(context.conformancyErrors, "Footer line read, but not at end of file.");
    }

    goNextLine(context, iter, BlastTabular());
}

/*!
 * @fn BlastTabularIn#readFooter
 * @headerfile seqan/blast.h
 * @brief Read the footer of a BlastTabular file (this is a NOOP).
 * @signature void readFooter(blastTabularIn);
 *
 * @param[in,out] blastTabularIn A @link BlastTabularIn @endlink formattedFile.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TContext>
inline void
readFooter(BlastTabularIn<TContext> & formattedFile)
{
    readFooter(context(formattedFile), formattedFile.iter, BlastTabular());
}



} // namespace seqan

#endif // SEQAN_BLAST_READ_BLAST_TABULAR_H_
