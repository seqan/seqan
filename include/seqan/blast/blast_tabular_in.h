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

/*!
 * @fn BlastTabular#onMatch
 * @brief returns whether the iterator is on beginning of match
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

template <typename TFwdIterator>
inline bool
onMatch(TFwdIterator & iter,
        BlastTabular const &)
{
    return (value(iter) != '#');
}

// ----------------------------------------------------------------------------
// Function readRecordHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabular#readRecordHeader
 * @brief read a the header of a record from a Blast tabular file
 * @headerfile seqan/blast.h
 *
 * @signature void readRecordHeader(blastRecord, stream, context, blastTabular);
 *
 * @param[out]    blastRecord  A @link BlastRecord @endlink whose
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
 * Please note that for or @link BlastIOContext::legacyFormat @endlink the @link BlastIOContext::fields @endlink member
 * is always ignored, however @link BlastIOContext::fieldsAsStrings @endlink is still read from the header, in case
 * you want to process it.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @see BlastTabular#onMatch
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

    // this is a record instead of a header
    if (onMatch(iter, BlastTabular()))
        SEQAN_THROW(ParseError("Header expected, but no header found."));

    int queryLinePresent = 0;
    int dbLinePresent = 0;
    int fieldsLinePresent = 0;
    int hitsLinePresent = 0;
    int lastLine = 0;

//     clear(context.buffer1);
//     clear(context.buffer2);
    auto & buffer = context.buffer1;
    auto & buffer2 = context.buffer2;

    while ((!atEnd(iter)) && (!onMatch(iter, BlastTabular())))// in Header
    {
        // skip '#'
        goNext(iter);
        //skip blanks following the '#'
        skipUntil(iter, IsGraph());

        clear(buffer);
        readUntil(buffer, iter, IsBlank());

        if (startsWith(buffer, "BLAST") || startsWith(buffer, "TBLAST"))
        {
            context.blastProgram = _programStringToTag(buffer);
            readLine(buffer, iter);
            context.versionString = buffer;

// TODO(h4nn3s): regex test
//             std::regex versionRE("[0-9]+\\.[0-9]+\\.[0-9]+\\+", std::regex::awk);
//             context.legacyFormat = std::regex_search(seqan::begin(key), seqan::end(key), versionRE);

            context.legacyFormat = (buffer.find('+') == std::string::npos);
        }
        else if (buffer == "Query:")
        {
            skipUntil(iter, IsGraph());
            readLine(r.qId, iter);
            ++queryLinePresent;
        }
        else if (buffer == "Database:")
        {
            skipUntil(iter, IsGraph());
            readLine(context.dbName, iter);
            ++dbLinePresent;
        }
        else if (buffer == "Fields:")
        {
            skipUntil(iter, IsGraph());

            clear(buffer);
            readLine(buffer, iter);
            _strSplit(context.fieldsAsStrings, buffer, std::regex(", "));

            ++fieldsLinePresent;
            if (context.legacyFormat)
            {
                // assume defaults for LEGACY
                if (!context.ignoreFieldsInHeader)
                    appendValue(context.fields, BlastMatchField<>::Enum::STD);

                break; // header is finished
            }

            if (context.ignoreFieldsInHeader)
                continue;

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
        else
        {
            readLine(buffer, iter);

            if (!context.legacyFormat)
            {
                // last line of BlastPlus Format
                if (startsWith(buffer, "BLAST processed"))
                {
                    ++lastLine;
                }
                // is hits counter?
                else if (endsWith(buffer, "hits found"))
                {
                    clear(buffer2);
                    for (unsigned i = 0; (i < length(buffer) && isdigit(buffer[i])); ++i)
                        appendValue(buffer2, buffer[i], Generous());

                    uint64_t hits = lexicalCast<uint64_t>(buffer2);

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
                    break; // header is finished
                }
                else
                {
                    appendValue(context.otherLines, buffer, Generous());
                }
            }
            else
            {
                appendValue(context.otherLines, buffer, Generous());
            }
        }
    }



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

    setBlastTabularSpec(context, BlastTabularSpec::HEADER);
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

template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFwdIterator,
          typename TScore,
          typename TString,
          BlastProgram p>
inline void
readRecordHeader(BlastRecord<TQId, TSId, TPos, TAlign> &,
                 TFwdIterator &,
                 BlastIOContext<TScore, TString, p, BlastTabularSpec::NO_HEADER> &,
                 BlastTabular const &)
{
    // NOOP for TABULAR without header
}

// ----------------------------------------------------------------------------
// Function skipHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabular#skipHeader
 * @brief skip a header from a Blast tabular input file.
 *
 * @signature void skipHeader(stream, context, blastTabular);
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
skipHeader(TFwdIterator & iter,
           BlastIOContext<TScore, TString, p, h> & context,
           BlastTabular const & /*tag*/)
{
    readRecordHeader(context.bufRecord, iter, context, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function skipUntilMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabular#skipUntilMatch
 * @brief skip arbitrary number of headers and/or comment lines until the beginning of a match is reached.
 *
 * @signature void skipUntilMatch(stream, blastTabular);
 *
 * @param[in,out] stream       An input iterator over a stream or any fwd-iterator over a string
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 *
 * @section Remarks
 *
 * Call this function whenever you are on a comment character ('#') in the file
 * and want to jump to the beginning of the next match. If you want to skip only
 * a single header (to count skipped headers or to verify its conformance
 * to standards), use BlastTabular#skipHeader instead.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @see BlastTabular#skipHeader
 * @headerfile seqan/blast.h
 */

template <typename TFwdIterator>
inline void
skipUntilMatch(TFwdIterator & iter,
               BlastTabular const & /*tag*/)
{
    while ((!atEnd(iter)) && value(iter) == '#') // skip comments
        skipLine(iter);
    if (atEnd(iter))
        SEQAN_THROW(ParseError("EOF reached without finding Match."));
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
            match.qId = context.buffer2;
            break;
//         case ENUM::Q_GI: write(s,  * ); break;
//         case ENUM::Q_ACC: write(s,  * ); break;
//         case ENUM::Q_ACCVER: write(s,  * ); break;
        case BlastMatchField<>::Enum::Q_LEN:
            match.qLength = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::S_SEQ_ID:
            match.sId = context.buffer2;
            break;
//         case ENUM::S_ALL_SEQ_ID: write(s,  * ); break;
//         case ENUM::S_GI: write(s,  * ); break;
//         case ENUM::S_ALL_GI: write(s,  * ); break;
//         case ENUM::S_ACC: write(s,  * ); break;
//         case ENUM::S_ACCVER: write(s,  * ); break;
//         case ENUM::S_ALLACC: write(s,  * ); break;
        case BlastMatchField<>::Enum::S_LEN:
            match.sLength = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::Q_START:
            match.qStart = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::Q_END:
            match.qEnd = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::S_START:
            match.sStart = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::S_END:
            match.sEnd = lexicalCast<TPos>(context.buffer2);
            break;
//         case ENUM::Q_SEQ: write(s,  * ); break;
//         case ENUM::S_SEQ: write(s,  * ); break;
        case BlastMatchField<>::Enum::E_VALUE:
            match.eValue = lexicalCast<double>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::BIT_SCORE:
            match.bitScore = lexicalCast<double>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::SCORE:
            match.alignStats.alignmentScore = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::LENGTH:
            match.alignStats.alignmentLength = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::P_IDENT:
            match.alignStats.alignmentIdentity = lexicalCast<double>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::N_IDENT:
            match.alignStats.numMatches = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::MISMATCH:
            match.alignStats.numMismatches = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::POSITIVE:
            match.alignStats.numPositiveScores = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::GAP_OPEN:
            match.alignStats.numGapOpens = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::GAPS:
            match.alignStats.numGaps = lexicalCast<TPos>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::P_POS:
            match.alignStats.alignmentSimilarity = lexicalCast<double>(context.buffer2);
        case BlastMatchField<>::Enum::FRAMES:
        {
            clear(context.buffers2);
            strSplit(context.buffers2, context.buffer2, EqualsChar<'/'>());
            if (length(context.buffers2) != 2)
                SEQAN_THROW(ParseError("Could not process frame string."));
            match.qFrameShift = lexicalCast<int8_t>(context.buffers2[0]);
            match.sFrameShift = lexicalCast<int8_t>(context.buffers2[1]);
        } break;
        case BlastMatchField<>::Enum::Q_FRAME:
            match.qFrameShift = lexicalCast<int8_t>(context.buffer2);
            break;
        case BlastMatchField<>::Enum::S_FRAME:
            match.sFrameShift = lexicalCast<int8_t>(context.buffer2);
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

/*!
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
    if (SEQAN_UNLIKELY(!onMatch(iter, BlastTabular())))
        SEQAN_THROW(ParseError("Not on beginning of Match (you should have skipped comments)."));

    match._maxInitialize(); // mark all members as not set

    clear(context.buffer1);
    clear(context.buffers1);

    auto & line = context.buffer1;
    readLine(line, iter);

    auto & fields = context.buffers1;
    strSplit(fields, line, EqualsChar<'\t'>());

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

                context.buffer2 = static_cast<decltype(context.buffer2)>(fields[n++]);
                _readField(match, context, f2);
            }
        } else
        {
            if (SEQAN_UNLIKELY(n >= length(fields)))
                SEQAN_THROW(ParseError("More columns expected than were present in file."));

            context.buffer2 = static_cast<decltype(context.buffer2)>(fields[n++]);
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

/*!
 * @fn BlastTabular#readMatch0
 * @brief read arbitrary columns from a file in BlastTabular
 *
 * @signature void readMatch0(stream, tag, args ...);
 *
 * @param[in,out] stream        An input iterator over a stream or any fwd-iterator over a string
 * @param[in]     blastTabular The @link BlastTabular @endlink tag.
 * @param[out]    args          Arbitrary typed variables
 *
 * @section Remarks
 *
 * Use this signature only if you do not or cannot use @link BlastMatch
 * @endlinkes. You can specify any number of arguments that are expected
 * to be able to hold the values in the columns read, i.e. if you pass a
 * double as argument and the value in the column cannot be successfully cast
 * to double, an exception will be thrown. If you want to be on the safe side,
 * you can pass CharStrings and evaluate them in another way.
 *
 * You may specify less columns than are available in the file, all but the first
 * n will be discarded.
 *
 * No transformations are made on the data, e.g. the positions are still
 * one-indexed and flipped for reverse strand matches.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @headerfile seqan/blast.h
 */

// arbitrary columns
template <typename TTarget>
inline SEQAN_FUNC_ENABLE_IF(IsSequence<TTarget>)
_assignOrCast(TTarget & target, std::string const & source)
{
    assign(target, source);
}

template <typename TTarget>
inline SEQAN_FUNC_ENABLE_IF(Is<NumberConcept<TTarget>>)
_assignOrCast(TTarget & target, std::string const & source)
{
    target = lexicalCast<TTarget>(source);
}

template <typename TFwdIterator,
          typename TArg>
inline void
_readMatchImplBlastTab(TFwdIterator & iter,
                       TArg & arg)
{
    static thread_local std::string buffer;
    clear(buffer);

    readUntil(buffer, iter, OrFunctor<IsTab,IsNewline>());
    _assignOrCast(arg, buffer);

    // as this is the last requested field, go to beginning of next line
    skipLine(iter);
}

template <typename TFwdIterator,
          typename TArg,
          typename... TArgs>
inline void
_readMatchImplBlastTab(TFwdIterator & iter,
                       TArg & arg,
                       TArgs & ... args)
{
    static thread_local std::string buffer;
    clear(buffer);

    readUntil(buffer, iter, IsTab());
    skipOne(iter, IsTab());
    _assignOrCast(arg, buffer);

    // recurse to next argument
    _readMatchImplBlastTab(iter, args...);
}

// custom arguments
template <typename TFwdIterator,
          typename... TArgs>
inline void
readMatch0(TFwdIterator & iter,
           BlastTabular const &,
           TArgs & ... args)
{
    // header should have been read or skipped
    if (SEQAN_UNLIKELY(!onMatch(iter, BlastTabular())))
        SEQAN_THROW(ParseError("ERROR: Not on beginning of Match (you should have skipped comments)."));

    _readMatchImplBlastTab(iter, args...);
}

// ----------------------------------------------------------------------------
// Function skipMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabular#skipMatch
 * @brief skip a line that contains a match
 *
 * @signature void skipMatch(stream, context, blastTabular);
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
 *
 * @see BlastTabular#skipHeader
 * @headerfile seqan/blast.h
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
_readRecordHeader(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
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
            if (!onMatch(iter, BlastTabular()))
            {
                appendValue(context.conformancyErrors,
                            "Less matches than promised by header");
                break;
            }

            readMatch(m, iter, context, BlastTabular());
        }

        if ((!atEnd(iter)) && onMatch(iter, BlastTabular()))
            appendValue(context.conformancyErrors,
                        "More matches than promised by header");

        while ((!atEnd(iter)) && onMatch(iter, BlastTabular()))
        {
            blastRecord.matches.emplace_back();
            readMatch(back(blastRecord.matches), iter, context, BlastTabular());
        }
    } else
    {
        while ((!atEnd(iter)) && onMatch(iter, BlastTabular()))
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
    clear(blastRecord);

    auto curId = context.lastId;

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

    while ((!atEnd(iter)) && onMatch(iter, BlastTabular()))
    {
        // no ID in buffer, yet
        if (SEQAN_LIKELY(empty(context.lastId) || empty(curId)))
        {
            // read current read ID
            readUntil(curId, iter, OrFunctor<IsTab,IsNewline>());
            skipOne(iter, IsTab());
            if (SEQAN_UNLIKELY(empty(context.lastId)))
                context.lastId = curId;
        }

        if (curId != context.lastId)
        {
            context.lastId = curId;
            break; // new Record reached
        }

        blastRecord.matches.emplace_back();
        // read remainder of line
        readMatch(back(blastRecord.matches), iter, context, BlastTabular());
        back(blastRecord.matches).qId = curId;
        std::swap(context.lastId, curId);
        clear(curId);
    }

    // revert to original behaviour
    context.fields.swap(fieldListMinusFirst);

    if (length(blastRecord.matches) == 0)
        SEQAN_THROW(ParseError("No Matches could be read."));

    blastRecord.qId = blastRecord.matches.front().qId;
}

/*!
 * @fn BlastTabular#readRecord
 * @brief read a record from a file in BlastTabular
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
 * (if it exists) and 0-n @link BlastMatch @endlinkes belonging to one query. For more details on this see
 * @link BlastTabular#readRecordHeader @endlink and @link BlastTabular#readMatch @endlink.
 *
 * Please note that if a record header exists the @link BlastIOContext::fields @endlink member is read from
 * it first and then used as in-parameter for reading the matches, i.e.
 * @link BlastTabular#readRecordHeader @endlink reads the
 * column labels and after that @link BlastTabular#readMatch @endlink expects the values in the columns to correspond
 * to their labels. If you do not want this behaviour, because you expect the header
 * to be non-standard -- as is the case with many tools -- you can set
 * context.@link BlastIOContext::ignoreFieldsInHeader @endlink to true.
 * In this case the you can also set
 * context.@link BlastIOContext::fields @endlink manually and these columns
 * will be expected for the matches, independent of the header. The latter
 * behaviour is also default when there are no headers.
 *
 * Please note also that if there are no headers the boundary
 * between records is inferred from the indentity of the first field, i.e.
 * non-standard field configurations must also have Q_SEQ_ID as their first @link BlastMatchField @endlink.
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
    _readRecordHeader(blastRecord, iter, context, BlastTabular());
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
    if (onMatch(iter, BlastTabular()))
    {
        setBlastTabularSpec(context, BlastTabularSpec::NO_HEADER);
        _readRecordNoHeader(blastRecord, iter, context, BlastTabular());
    }
    else
    {
        setBlastTabularSpec(context, BlastTabularSpec::HEADER);
        _readRecordHeader(blastRecord, iter, context, BlastTabular());
    }
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
 * @brief read the header (top-most section) of a BlastTabular file (this is a NOOP)
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
readHeader(BlastIOContext<TScore, TConString, p, h> &,
           TFwdIterator &,
           BlastTabular const & /*tag*/)
{
}

/*!
 * @fn BlastTabularIn#readHeader
 * @headerfile seqan/blast.h
 * @brief read the header (top-most section) of a BlastTabular file (this is a NOOP)
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
 * @brief read the footer of a BlastTabular file (currently NOOP)
 * @signature void readFooter(stream, context, blastTabular);
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
readFooter(BlastIOContext<TScore, TConString, p, h> &,
           TFwdIterator &,
           BlastTabular const & /*tag*/)
{
    //TODO check if this is really NO-OP for this format
}

/*!
 * @fn BlastTabularIn#readFooter
 * @headerfile seqan/blast.h
 * @brief read the footer of a BlastTabular file (currently NOOP)
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
