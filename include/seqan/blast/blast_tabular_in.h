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

/* IMPLEMENTATION NOTES

BLAST TABULAR example:

The format of a blast tabular output file is less simple than it looks, here's
the general form

COMMENTLINES
 MATCH
 MATCH
 MATCH
 MATCH
COMMENTLINES
 MATCH
COMMENTLINES
COMMENTLINES
...

=> COMMENTLINES for each sequence, 0-n Matchs for each sequence
=> Each record is one-line, each COMMENTLINES is multiline

COMMENTLINES usually consists of:

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

Because 0 matches are allowed, multiple COMMENTLINES can succeed each other, the
criterium for separation employed by this implementation is that an NCBI Blast
COMMENTLINES always ends after the "Fields" line and NCBI Blast+ COMMENTLINES end after
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
 * @class BlastTabularFileIn
 * @signature template <typename TBlastIOContext>
 * using BlastTabularFileIn = FormattedFile<BlastTabular, Input, TBlastIOContext>;
 * @extends FormattedFileIn
 * @headerfile <seqan/blast.h>
 * @brief FormattedFileIn abstraction for @link BlastTabular @endlink
 *
 * This is a @link FormattedFile @endlink specialization for reading @link BlastTabular @endlink formats. For details
 * on how to influence the reading of files and how to differentiate between the tabular format without comment lines
 * and the one with comment lines, see @link BlastIOContext @endlink.
 * Please note that you have specify the type of the context as a template parameter to BlastTabularFileIn.
 *
 * @section Overview
 *
 * <ul>
 * <li> open @link BlastTabularFileIn @endlink</li>
 * <li> @link BlastTabularFileIn#readHeader @endlink </li>
 * <li> while @link BlastTabularFileIn#onRecord @endlink </li>
 * <ul><li> @link BlastTabularFileIn#readRecord @endlink</li></ul>
 * <li> @link BlastTabularFileIn#readFooter @endlink </li>
 * </ul>
 *
 * For a detailed example have a look at the
 * <a href="http://seqan.readthedocs.io/en/develop/Tutorial/InputOutput/BlastIO.html">Blast IO tutorial</a>.
 *
 * @see BlastRecord
 */

template <typename TBlastIOContext = BlastIOContext<>>
using BlastTabularFileIn = FormattedFile<BlastTabular, Input, TBlastIOContext>;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

template <typename TScore,
          typename TFwdIterator,
          BlastProgram p,
          BlastTabularSpec h>
inline bool
atEnd(BlastIOContext<TScore, p, h> const & context,
      TFwdIterator const & iter,
      BlastTabular const &)
{
    return (atEnd(iter) && empty(context._lineBuffer));
}

template <typename TContext>
inline bool
atEnd(BlastTabularFileIn<TContext> & formattedFile)
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
        context(file).tabularSpec = BlastTabularSpec::COMMENTS;
    else
        context(file).tabularSpec = BlastTabularSpec::NO_COMMENTS;
    return true;
}

// ----------------------------------------------------------------------------
// Function _onMatch()
// ----------------------------------------------------------------------------

template <typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline bool
_onMatch(BlastIOContext<TScore, p, h> & context,
        BlastTabular const &)
{
    return (!startsWith(context._lineBuffer, "#") && !empty(context._lineBuffer));
}

// ----------------------------------------------------------------------------
// Function onRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularFileIn#onRecord
 * @brief Returns whether the currently buffered line looks like the start of a record.
 * @signature bool onRecord(blastTabularIn);
 * @headerfile seqan/blast.h
 *
 * @param[in,out] blastTabularIn A @link BlastTabularFileIn @endlink formattedFile.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @return bool true or false
 */

template <typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline bool
onRecord(BlastIOContext<TScore, p, h> & context,
         BlastTabular const &)
{
    if (context.tabularSpec == BlastTabularSpec::NO_COMMENTS)
        return _onMatch(context, BlastTabular());

    //      comment lines                                  Footer
    return (startsWith(context._lineBuffer, "#") && (!startsWith(context._lineBuffer, "# BLAST processed ")));
}

template <typename TContext>
inline bool
onRecord(BlastTabularFileIn<TContext> & formattedFile)
{
    return onRecord(context(formattedFile), BlastTabular());
}

// ----------------------------------------------------------------------------
// Function _goNextLine()
// ----------------------------------------------------------------------------

template <typename TScore,
          typename TFwdIterator,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_goNextLine(BlastIOContext<TScore, p, h> & context,
           TFwdIterator & iter,
           BlastTabular const &)
{
    clear(context._lineBuffer);
    if (!atEnd(iter))
        readLine(context._lineBuffer, iter);
}

// ----------------------------------------------------------------------------
// Function _readCommentLines()
// ----------------------------------------------------------------------------

template <typename ... TSpecs,
          typename TFwdIterator,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readCommentLinesImpl(BlastRecord<TSpecs...> & r,
                      TFwdIterator & iter,
                      BlastIOContext<TScore, p, h> & context,
                      BlastTabular const &)
{
    // this is a match instead of comment lines
    if (_onMatch(context, BlastTabular()))
        SEQAN_THROW(ParseError("Commen lines expected, but not found."));
    else
        context.tabularSpec = BlastTabularSpec::COMMENTS;

    int queryLinePresent = 0;
    int dbLinePresent = 0;
    int fieldsLinePresent = 0;
    int hitsLinePresent = 0;

    do
    {
        if (std::regex_search(std::begin(context._lineBuffer),
                              std::end(context._lineBuffer),
                              std::regex("^# T?BLAST")))
        {
            // last line of file
            if (SEQAN_UNLIKELY(startsWith(context._lineBuffer, "# BLAST processed ") && !context.legacyFormat))
            {
                SEQAN_FAIL("ERROR: You called readRecord() when you should have called readFooter()."
                           "Always check onRecord() before calling readRecord().");
            }
            else // first line of the comments
            {
                assign(context.versionString, suffix(context._lineBuffer, 2));
                context.blastProgram = _programStringToTag(prefix(context.versionString,
                                                                  std::find(begin(context.versionString, Standard()),
                                                                            end(context.versionString, Standard()),
                                                                            ' ')
                                                                  - begin(context.versionString, Standard())));

                context.legacyFormat = !std::regex_search(std::begin(context.versionString),
                                                          std::end(context.versionString),
                                                          std::regex("\\d\\.\\d\\.\\d{1,2}\\+"));
            }
        }
        else if (startsWith(context._lineBuffer, "# Query:"))
        {
            ++queryLinePresent;
            assign(r.qId, suffix(context._lineBuffer, 9));

        }
        else if (startsWith(context._lineBuffer, "# Database:"))
        {
            ++dbLinePresent;
            assign(context.dbName, suffix(context._lineBuffer, 12));
        }
        else if (startsWith(context._lineBuffer, "# Fields:"))
        {
            ++fieldsLinePresent;
            clear(context._stringBuffer);
            // make sure target is big enough
            resize(context._stringBuffer, length(context._lineBuffer) - 10);
            // use escape character as placeholder for replacement since strSplit only handles single characters
            std::regex_replace(std::begin(context._stringBuffer),
                               std::begin(context._lineBuffer) + 10, // skip "# Fields:"
                               std::end(context._lineBuffer),
                               std::regex(", "),
                               "\x7F");
            // shrink back down to actual size (replacing two letters with one makes string shorter!)
            resize(context._stringBuffer, length(context._stringBuffer.c_str()));
            strSplit(context.fieldsAsStrings, context._stringBuffer, EqualsChar<'\x7F'>(), true);

            if (context.legacyFormat)
            {
                // assume defaults for LEGACY
                if (!context.ignoreFieldsInComments)
                    appendValue(context.fields, BlastMatchField<>::Enum::STD);
            }
            else if (!context.ignoreFieldsInComments)
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
                    for (unsigned i = 2; (i < length(context._lineBuffer) && isdigit(context._lineBuffer[i])); ++i)
                        appendValue(context._stringBuffer, context._lineBuffer[i], Generous());

                    uint64_t hits = lexicalCast<uint64_t>(context._stringBuffer);

                    if (hits)
                    {
                        resize(r.matches, hits);
                    }
                    else  // hits = 0 means no fieldList, restore default
                    {
                        appendValue(context.fields, BlastMatchField<>::Enum::STD);

                        context._stringBuffer = std::regex_replace(BlastMatchField<>::columnLabels[0],
                                                                   std::regex(", "),
                                                                   "\x7F");
                        strSplit(context.fieldsAsStrings, context._stringBuffer, EqualsChar<'\x7F'>(), true);
                    }

                    ++hitsLinePresent;
//                     break; // comments are finished
                }
                else
                {
                    appendValue(context.otherLines, suffix(context._lineBuffer, 1), Generous());
                }
            }
            else
            {
                appendValue(context.otherLines, suffix(context._lineBuffer, 1), Generous());
            }
        }

        _goNextLine(context, iter, BlastTabular());

    } while (startsWith(context._lineBuffer, "#") &&                   // still on comments
             !std::regex_search(std::begin(context._lineBuffer),// start of next record
                                std::end(context._lineBuffer),
                                std::regex("^# T?BLAST")));

    if (context.blastProgram == BlastProgram::UNKNOWN)
        appendValue(context.conformancyErrors,
                    "Type of BlastProgram could not be determined from comments, you are advised to look "
                    "at context.versionString and context.otherLines.");

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
        // is omitted in BLAST_PLUS when there are no hits
        if ((fieldsLinePresent != 1) && (length(r.matches) > 0))
            appendValue(context.conformancyErrors, "No or multiple fields lines present.");
    }

    if (!empty(context.otherLines))
        appendValue(context.conformancyErrors, "Unexpected lines present, see context.otherLines.");
}

template <typename ... TSpecs,
          typename TFwdIterator,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readCommentLines(BlastRecord<TSpecs...> & r,
                  TFwdIterator & iter,
                  BlastIOContext<TScore, p, h> & context,
                  BlastTabular const &)
{
    ++context._numberOfRecords;

    if (context.tabularSpec == BlastTabularSpec::NO_COMMENTS)
        return;

    clear(r);
    clear(context.versionString);
    clear(context.dbName);
    clear(context.otherLines);
    clear(context.conformancyErrors);
    clear(context.fieldsAsStrings);

    if (!context.ignoreFieldsInComments)
        clear(context.fields);

    _readCommentLinesImpl(r, iter, context, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function _readField()
// ----------------------------------------------------------------------------

template <typename TAlignRow0,
          typename TAlignRow1,
          typename TPos,
          typename TQId,
          typename TSId,
          typename TScore,
          typename ... TSpecs,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readField(BlastMatch<TAlignRow0, TAlignRow1, TPos, TQId, TSId> & match,
           BlastRecord<TSpecs...> & record,
           BlastIOContext<TScore, p, h> & context,
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
        case BlastMatchField<>::Enum::Q_ACC:
            appendValue(record.qAccs, context._stringBuffer);
            break;
//         case ENUM::Q_ACCVER: write(s,  * ); break;
        case BlastMatchField<>::Enum::Q_LEN:
            match.qLength = lexicalCast<TPos>(context._stringBuffer);
            record.qLength = match.qLength;
            break;
        case BlastMatchField<>::Enum::S_SEQ_ID:
            match.sId = context._stringBuffer;
            break;
//         case ENUM::S_ALL_SEQ_ID: write(s,  * ); break;
//         case ENUM::S_GI: write(s,  * ); break;
//         case ENUM::S_ALL_GI: write(s,  * ); break;
        case BlastMatchField<>::Enum::S_ACC:
            appendValue(match.sAccs, context._stringBuffer);
            break;
//         case ENUM::S_ACCVER: write(s,  * ); break;
        case BlastMatchField<>::Enum::S_ALLACC:
            strSplit(match.sAccs, context._stringBuffer, EqualsChar<';'>());
            break;
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
            break;
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
        case BlastMatchField<>::Enum::S_TAX_IDS:
        {
            StringSet<CharString> temp;
            strSplit(temp, context._stringBuffer, EqualsChar<';'>());
            for (auto const & s : temp)
                appendValue(match.sTaxIds, lexicalCast<uint64_t>(s));
        } break;
//         case ENUM::S_SCI_NAMES: write( * ); break;
//         case ENUM::S_COM_NAMES: write( * ); break;
//         case ENUM::S_BLAST_NAMES: write( * ); break;
//         case ENUM::S_S_KINGDOMS: write( * ); break;
//         case ENUM::S_TITLE: write( * ); break;
//         case ENUM::S_ALL_TITLES: write( * ); break;
//         case ENUM::S_STRAND: write( * ); break;
//         case ENUM::Q_COV_S: write( * ); break;
//         case ENUM::Q_COV_HSP:
        case BlastMatchField<>::Enum::LCA_ID:
            record.lcaId = context._stringBuffer;
            break;
        case BlastMatchField<>::Enum::LCA_TAX_ID:
            record.lcaTaxId = lexicalCast<uint32_t>(context._stringBuffer);
            break;
        default:
            SEQAN_THROW(ParseError("The requested column type is not yet "
                                   "implemented."));
    };
}

// ----------------------------------------------------------------------------
// Function _readMatch()
// ----------------------------------------------------------------------------

template <typename TAlignRow0,
          typename TAlignRow1,
          typename TPos,
          typename TQId,
          typename TSId,
          typename TFwdIterator,
          typename TScore,
          typename ... TSpecs,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readMatch(BlastMatch<TAlignRow0, TAlignRow1, TPos, TQId, TSId> & match,
           BlastRecord<TSpecs...> & record,
           TFwdIterator & iter,
           BlastIOContext<TScore, p, h> & context,
           BlastTabular const &)
{
    if (context.legacyFormat)
    {
        #if SEQAN_ENABLE_DEBUG
        if ((length(context.fields) != 1) || (context.fields[0] != BlastMatchField<>::Enum::STD))
            std::cerr << "Warning: custom fields set, but will be reset, because legacyFormat is also set.\n";
        #endif
        // set defaults
        clear(context.fields);
        appendValue(context.fields,  BlastMatchField<>::Enum::STD);
    }

    // comments should have been read or skipped
    if (SEQAN_UNLIKELY(!_onMatch(context, BlastTabular())))
        SEQAN_THROW(ParseError("Not on beginning of Match (you should have skipped comments)."));

    match._maxInitialize(); // mark all members as "not set"

    clear(context._setBuffer1);

    auto & fields = context._setBuffer1;
    strSplit(fields, context._lineBuffer, IsTab());

    // line in buffer was processed so new line is read from iter into buffer
    _goNextLine(context, iter, BlastTabular());

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
                _readField(match, record, context, f2);
            }
        } else
        {
            if (SEQAN_UNLIKELY(n >= length(fields)))
                SEQAN_THROW(ParseError("More columns expected than were present in file."));

            context._stringBuffer = static_cast<decltype(context._stringBuffer)>(fields[n++]);
            _readField(match, record, context, f);
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
// Function readRecord()
// ----------------------------------------------------------------------------

template <typename TFwdIterator,
          typename ... TSpecs,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readRecordWithCommentLines(BlastRecord<TSpecs...> & blastRecord,
                            TFwdIterator & iter,
                            BlastIOContext<TScore, p, h> & context,
                            BlastTabular const &)
{
    _readCommentLines(blastRecord, iter, context, BlastTabular());

    if (!context.legacyFormat) // this is detected from the comments
    {
        // .matches already resized for us
        for (auto & m : blastRecord.matches)
        {
            if (!_onMatch(context, BlastTabular()))
            {
                appendValue(context.conformancyErrors,
                            "Less matches than promised by the comments");
                break;
            }

            _readMatch(m, blastRecord, iter, context, BlastTabular());
        }

        if ((!atEnd(context, iter, BlastTabular())) && _onMatch(context, BlastTabular()))
            appendValue(context.conformancyErrors,
                        "More matches than promised by the comments");

        while ((!atEnd(context, iter, BlastTabular())) && _onMatch(context, BlastTabular()))
        {
            blastRecord.matches.emplace_back();
            _readMatch(back(blastRecord.matches), blastRecord, iter, context, BlastTabular());
        }
    } else
    {
        while ((!atEnd(context, iter, BlastTabular())) && _onMatch(context, BlastTabular()))
        {
            blastRecord.matches.emplace_back();
            _readMatch(back(blastRecord.matches), blastRecord, iter, context, BlastTabular());
        }
    }
}

template <typename ... TSpecs,
          typename TFwdIterator,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
_readRecordWithoutCommentLines(BlastRecord<TSpecs...> & blastRecord,
                               TFwdIterator & iter,
                               BlastIOContext<TScore, p, h> & context,
                               BlastTabular const &)
{
    ++context._numberOfRecords;

    clear(blastRecord);

    auto it = begin(context._lineBuffer, Rooted()); // move into line below when seqan supports && properly
    readUntil(blastRecord.qId, it, IsTab());

    auto curIdPlusTab = blastRecord.qId;
    appendValue(curIdPlusTab, '\t');

    if ((context.fields[0] != BlastMatchField<>::Enum::STD) &&
        (context.fields[0] != BlastMatchField<>::Enum::Q_SEQ_ID))
    {
        SEQAN_FAIL("ERROR: readRecord interface on comment-less format with custom "
                   "fields not supported, unless first custom field is "
                   "Q_SEQ_ID. Use the lowlevel readMatch interface instead.");
    }

    while ((!atEnd(context, iter, BlastTabular())) && _onMatch(context, BlastTabular()))
    {
        blastRecord.matches.emplace_back();
        // read remainder of line
        _readMatch(back(blastRecord.matches), blastRecord, iter, context, BlastTabular());

        if (!startsWith(context._lineBuffer, curIdPlusTab)) // next record reached
            break;
    }

    if (empty(blastRecord.matches))
        SEQAN_THROW(ParseError("No Matches could be read."));
}

/*!
 * @fn BlastTabularFileIn#readRecord
 * @brief Read a record from a file in BlastTabular format.
 * @headerfile seqan/blast.h
 * @signature void readRecord(blastRecord, blastTabularIn);
 *
 * @param[out]    blastRecord  A @link BlastRecord @endlink to hold all information related to one query sequence.
 * @param[in,out] blastTabularIn A @link BlastTabularFileIn @endlink formattedFile.
 *
 * @section Remarks
 *
 * This function will read an entire record from a blast tabular file, i.e. it will read the comment lines
 * (if the format is @link BlastTabularSpec::COMMENTS @endlink) and 0-n @link BlastMatch @endlinkes belonging
 * to one query.
 *
 * Please note that if there are no comment lines in the file the boundary
 * between records is inferred from the indentity of the first field, i.e.
 * non-standard field configurations must also have Q_SEQ_ID as their first @link BlastMatchField @endlink.
 *
 * @subsection Comment lines
 *
 * The @link BlastRecord::qId @endlink member of the record is read from the comment lines and the
 * @link BlastRecord::matches @endlink are resized to the expected number of matches succeeding the comments.
 *
 * This function also sets many properties of blastTabularIn's @link BlastIOContext @endlink, including these members:
 * <li> @link BlastIOContext::versionString @endlink: version string of the program.</li>
 * <li> @link BlastIOContext::dbName @endlink: name of the database.</li>
 * <li> @link BlastIOContext::fields @endlink: descriptors for the columns.</li>
 * <li> @link BlastIOContext::fieldsAsStrings @endlink: labels of the columns
 * as they appear in the file.</li>
 * <li> @link BlastIOContext::conformancyErrors @endlink: if this StringSet is not empty, then there are issues in the
 * comments.</li>
 * <li> @link BlastIOContext::otherLines @endlink: any lines that cannot be interpreted; these always also imply
 * conformancyErrors.</li>
 * <li>@link BlastIOContext::legacyFormat @endlink: whether the record (and likely the entire file) is in
 * legacyFormat.</li>
 *
 * It also sets the blast program run-time parameter of the context depending on the information found in the comments. If
 * the compile time parameter was set on the context and they are different this will result in a critical error.
 *
 * Please note that for @link BlastIOContext::legacyFormat @endlink the @link BlastIOContext::fields @endlink member
 * is always ignored, however @link BlastIOContext::fieldsAsStrings @endlink is still read from the comments, in case
 * you want to process it.
 *
 * In case you do not wish the @link BlastIOContext::fields @endlink to be read from the comments, you can set
 * context.@link BlastIOContext::ignoreFieldsInComments @endlink to true. This will be prevent it from being read and will
 * allow you to specify it manually which might be relevant for reading the match lines.
 *
 * If the format is @link BlastTabularSpec::NO_COMMENTS @endlink none of the above happens and
 * @link BlastRecord::qId @endlink is derived from the first match.
 *
 * @subsection Matches
 *
 * A match line contains 1 - n columns or fields, 12 by default.
 * The @link BlastIOContext::fields @endlink member of the context is considered when reading these fields. It is
 * usually extracted from the comment lines but can also be set by yourself if there are no comments or if you want
 * to overwrite the comments' information (see above).
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
 * In contrast to @link BlastTabularFileOut#writeRecord @endlink no other transformations
 * are made, e.g. the positions are still one-indexed and
 * flipped for reverse strand matches. This is due to the required fields for
 * retransformation (sequence lengths, frames) not being available in the
 * default columns.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename ... TSpecs,
          typename TFwdIterator,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
readRecord(BlastRecord<TSpecs...> & blastRecord,
           TFwdIterator & iter,
           BlastIOContext<TScore, p, h> & context,
           BlastTabular const &)
{
    if (context.tabularSpec == BlastTabularSpec::NO_COMMENTS)
        _readRecordWithoutCommentLines(blastRecord, iter, context, BlastTabular());
    else
        _readRecordWithCommentLines(blastRecord, iter, context, BlastTabular());
}

template <typename ... TSpecs,
          typename TContext>
inline void
readRecord(BlastRecord<TSpecs...> & blastRecord,
           BlastTabularFileIn<TContext> & formattedFile)
{
    readRecord(blastRecord, formattedFile.iter, context(formattedFile), BlastTabular());
}

// ----------------------------------------------------------------------------
// Function readHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularFileIn#readHeader
 * @headerfile seqan/blast.h
 * @brief Read the header (top-most section) of a BlastTabular file.
 * @signature void readHeader(blastTabularIn);
 *
 * @param[in,out] blastTabularIn A @link BlastTabularFileIn @endlink formattedFile.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TFwdIterator,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
readHeader(BlastIOContext<TScore, p, h> & context,
           TFwdIterator & iter,
           BlastTabular const & /*tag*/)
{
    readLine(context._lineBuffer, iter); // fill the line-buffer the first time

    if (context.tabularSpec == BlastTabularSpec::UNKNOWN)
    {
        if (_onMatch(context, BlastTabular()))
            context.tabularSpec = BlastTabularSpec::NO_COMMENTS;
        else
            context.tabularSpec = BlastTabularSpec::COMMENTS;
    }
}

template <typename TContext>
inline void
readHeader(BlastTabularFileIn<TContext> & formattedFile)
{
    readHeader(context(formattedFile), formattedFile.iter, BlastTabular());
}

// ----------------------------------------------------------------------------
// Function readFooter()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastTabularFileIn#readFooter
 * @headerfile seqan/blast.h
 * @brief Read the footer (bottom-most section) of a BlastTabular file.
 * @signature void readFooter(blastTabularIn);
 *
 * @param[in,out] blastTabularIn A @link BlastTabularFileIn @endlink formattedFile.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 */

template <typename TFwdIterator,
          typename TScore,
          BlastProgram p,
          BlastTabularSpec h>
inline void
readFooter(BlastIOContext<TScore, p, h> & context,
           TFwdIterator & iter,
           BlastTabular const & /*tag*/)
{
    clear(context.otherLines);
    clear(context.conformancyErrors);
    clear(context.fieldsAsStrings);
    clear(context.fields);

    if ((context.tabularSpec == BlastTabularSpec::COMMENTS) && !context.legacyFormat)
    {
        if (SEQAN_UNLIKELY(!startsWith(context._lineBuffer, "# BLAST processed")))
        {
            std::cout << "\"" << context._lineBuffer << "\"" << std::endl;
            SEQAN_FAIL("ERROR: Tried to read footer, but was not on footer.");
        }
        clear(context._stringBuffer);
        auto it = seqan::begin(context._lineBuffer);
        it += 18; // skip "BLAST processed "
        readUntil(context._stringBuffer, it,  IsBlank());

        uint64_t numRecords = lexicalCast<uint64_t>(context._stringBuffer);

        clear(context.conformancyErrors);
        if (context._numberOfRecords < numRecords)
            appendValue(context.conformancyErrors, "The file claims to contain more records than you read.");
        else if (context._numberOfRecords > numRecords)
            appendValue(context.conformancyErrors, "The file claims to contain less records than you read.");
    }

    if (!atEnd(iter))
        appendValue(context.conformancyErrors, "Expected to be at end of file.");

    _goNextLine(context, iter, BlastTabular());
}

template <typename TContext>
inline void
readFooter(BlastTabularFileIn<TContext> & formattedFile)
{
    readFooter(context(formattedFile), formattedFile.iter, BlastTabular());
}

} // namespace seqan

#endif // SEQAN_BLAST_READ_BLAST_TABULAR_H_
