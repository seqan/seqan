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

template <typename TSpec = void>
struct PrivateSpace_
{
    static std::vector<BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum>
    const defaultsMinusQId;

};

template <typename TSpec>
std::vector<BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum>
const PrivateSpace_<TSpec>::defaultsMinusQId =
{
    {
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::S_SEQ_ID,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::P_IDENT,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::LENGTH,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::MISMATCH,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::GAP_OPEN,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::Q_START,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::Q_END,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::S_START,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::S_END,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::E_VALUE,
        BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum::BIT_SCORE
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function strSplit()
// ----------------------------------------------------------------------------

template <typename TString,
          typename TSpec,
          typename TSequence>
inline void
strSplit(StringSet<TString, TSpec> & result,
         TSequence const & sequence,
         std::regex const & regex,
         bool allowEmptyStrings = true)
{
    if (!allowEmptyStrings)
        SEQAN_FAIL("This strSplit overload works only if emptyStrings allowed");

    std::string tmp;
    tmp = std::regex_replace(sequence, regex, "\x7F");
    strSplit(result, tmp, EqualsChar<'\x7F'>(), true);
}

// ----------------------------------------------------------------------------
// Function onMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#onMatch
 * @brief returns whether the iterator is on beginning of match
 *
 * @signature bool onMatch(iter, tag)
 *
 * @param[in] iter    An input iterator over a stream or any fwd-iterator over a string
 * @param[in] tag     @link BlastFormat @endlink specialization
 *
 * @throw IOError On low-level I/O errors.
 *
 * @return    true or false
 */

template <typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline bool
onMatch(TFwdIterator & iter,
        BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    return (value(iter) != '#');
}

template <typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline bool
onMatch(TFwdIterator & iter,
        BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    return onMatch(iter, BlastFormat<BlastFormatFile::TABULAR, p, g>());
}

// ----------------------------------------------------------------------------
// Function _verifyFields()
// ----------------------------------------------------------------------------

// // fields as string of blastmatchfield::enum
// template <typename TFieldList,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline SEQAN_FUNC_ENABLE_IF(And<IsSequence<TFieldList>,
//                                 IsSameType<typename Value<TFieldList>::Type,
//                                            typename BlastMatchField<g>::Enum>>)
// _verifyFields(BlastIOContext &,
//               TFieldList const & fields,
//               unsigned const hits, // irrelevant for traditional header
//               BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     // In BLAST_PLUS, iff hits are zero no fields line is printed
//     if (g == BlastFormatGeneration::BLAST_PLUS)
//         if ((hits == 0) && length(fields) == 0 )
//             return;
//
//     if (length(fields) != BlastMatchField<g>::defaults)
//         SEQAN_THROW(RecoverableParseError("Default fields and header "
//                                           "verification requested, but "
//                                           "fields were not default."));
//
//     for (unsigned i = 0; i < length(fields); ++i)
//         if (fields[i] != BlastMatchField<g>::defaults[i])
//             SEQAN_THROW(RecoverableParseError("Default fields and header "
//                                               "verification requested, but "
//                                               "fields were not default."));
// }
//
// // fields as strings
// template <typename TString,
//           typename TSpec,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// _verifyFields(BlastIOContext & context,
//               StringSet<TString, TSpec> const & fields,
//               unsigned const hits, // irrelevant for traditional header
//               BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g>  TFormat;
//
//     // In BLAST_PLUS, iff hits are zero no fields line is printed
//     if (g == BlastFormatGeneration::BLAST_PLUS)
//         if ((hits == 0) && length(fields) == 0 )
//             return;
//
//     context.buffer1 = concat(fields, _seperatorString(TFormat()));
//
//     // will always be zero but this is cleaner
//     const uint8_t std = static_cast<uint8_t>(BlastMatchField<g>::Enum::STD);
//
//     // compare the first twelve fields
//     if (prefix(context.buffer1, length(BlastMatchField<g>::columnLabels[std])) !=
//         BlastMatchField<g>::columnLabels[std])
//         SEQAN_THROW(RecoverableParseError("Default fields and header "
//                                           "verification requested, but "
//                                           "fields were not default."));
// }

// ----------------------------------------------------------------------------
// Function readHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastRecord#readHeader
 * @brief read a Header from a Blast output file
 * @headerfile seqan/blast.h
 *
 * @signature readHeader(blastRecord, iter, context, tag);
 *
 * @param[out]  blastRecord A @link BlastRecord @endlink whose
 * @param[in,out] iter      An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context   An @link BlastIOContext @endlink with buffers.
 * @param[in]   tag         @link BlastFormat @endlink tag, only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported.
 *
 * Call this function on every line beginning that is not "onMatch". Will
 * set the @link BlastRecord::qId @endlink member of the record. For
 * @link BlastFormatGeneration::BLAST_PLUS @endlink it will also resize
 * the @link BlastRecord::matches @endlink member to the expected number of
 * matches that will succeed the header.
 *
 * This function sets the following members of the @link BlastIOContext @endlink
 * which may or may not intereset you:
 * <li> @link BlastDbSpecs::dbName @endlink member of
 * @link BlastIOContext::dbSpecs @endlink (name of the database)</li>
 * <li> @link BlastIOContext::headerConformsToStandards @endlink a bool
 * signafying whether the header contained no anomalies.</li>
 * <li> @link BlastIOContext::fields @endlink: descriptors for the columns</li>
 * <li> @link BlastIOContext::fieldsAsStrings @endlink: labels for the columns
 * as they appear in the file</li>
 * <li> @link BlastIOContext::otherLines @endlink: any unknown lines</li>
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @see BlastFormat#onMatch
 */

template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_readHeaderImpl(BlastRecord<TQId, TSId, TPos, TAlign> & r,
                TFwdIterator & iter,
                BlastIOContext & context,
                BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    typedef BlastMatchField<BlastFormatGeneration::BLAST_PLUS> TMatchField;

    // this is a record instead of a header
    if (onMatch(iter, TFormat()))
        SEQAN_THROW(ParseError("Header expected, but no header found."));

    bool keyCorrect = false;
    int queryLinePresent = 0;
    int dbLinePresent = 0;
    int fieldsLinePresent = 0;
    int hitsLinePresent = 0;
    int lastLine = 0;

    clear(context.buffer1);
    clear(context.buffer2);
    auto & key = context.buffer1;
    auto & buf = context.buffer2;

    while ((!atEnd(iter)) && (!onMatch(iter, TFormat())))// in Header
    {
        // skip '#'
        goNext(iter);
        //skip blanks following the '#'
        skipUntil(iter, IsGraph());

        clear(key);
        readUntil(key, iter, IsBlank());

        if (startsWith(key, "BLAST"))
        {
            if (key == _programTagToString(TFormat()))
                keyCorrect = true;

            readLine(key, iter);
            context.versionString = key;
        }
        else if (key == "Query:")
        {
            skipUntil(iter, IsGraph());
            readLine(r.qId, iter);
            ++queryLinePresent;
        }
        else if (key == "Database:")
        {
            skipUntil(iter, IsGraph());
            readLine(context.dbSpecs.dbName, iter);
            ++dbLinePresent;
        }
        else if (key == "Fields:")
        {
            skipUntil(iter, IsGraph());

            clear(buf);
            readLine(buf, iter);
            strSplit(context.fieldsAsStrings, buf, std::regex(", "));

            ++fieldsLinePresent;
            if (g == BlastFormatGeneration::BLAST_LEGACY)
            {
                // assume defaults for LEGACY
                if (!context.ignoreFieldsInHeader)
                    appendValue(context.fields, TMatchField::Enum::STD);

                break; // header is finished
            }

            if (context.ignoreFieldsInHeader)
                continue;

            bool defaults = true;
            resize(context.fields, length(context.fieldsAsStrings));
            for (uint8_t j = 0; j < length(context.fieldsAsStrings); ++j)
            {
                for (uint8_t i = 0; i < length(TMatchField::columnLabels); ++i)
                {
                    if (context.fieldsAsStrings[j] ==
                        TMatchField::columnLabels[i])
                    {
                        context.fields[j] =
                            static_cast<TMatchField::Enum>(i);

                        if ((j >= length(TMatchField::defaults)) &&
                            (static_cast<TMatchField::Enum>(i) !=
                             TMatchField::defaults[j]))
                            defaults = false;
                        break;
                    }
                }
            }
            if (defaults) // replace multiple fields with the default meta field
            {
                clear(context.fields);
                appendValue(context.fields, TMatchField::Enum::STD);
            }
        }
        else
        {
            readLine(key, iter);

            if (g == BlastFormatGeneration::BLAST_PLUS)
            {
                // last line of BlastPlus Format
                if (startsWith(key, "BLAST processed"))
                {
                    ++lastLine;
                }
                // is hits counter?
                else if (endsWith(key, "hits found"))
                {
                    clear(buf);
                    for (unsigned i = 0;
                         (i < length(key) && isdigit(key[i]));
                         ++i)
                        append(buf, key[i], Generous());

                    uint64_t hits = lexicalCast<uint64_t>(buf);
//                     if (ret && strict)
//                         throw std::ios_base::failure(BAD_FORMAT);

                    if (hits)
                    {
                        //resize(r.matches, hits);
                        r.matches.resize(hits); // TODO fix this
                    }
                    else  // hits = 0 means no fieldList, restore default
                    {
                        appendValue(context.fields, TMatchField::Enum::STD);
                        strSplit(context.fieldsAsStrings,
                                 TMatchField::columnLabels[0],
                                 std::regex(", "));
                    }

                    ++hitsLinePresent;
                    break; // header is finished
                }
                else
                {
                    appendValue(context.otherLines, key, Generous());
                }
            }
            else
            {
                appendValue(context.otherLines, key, Generous());
            }
        }
    }


    if (g == BlastFormatGeneration::BLAST_LEGACY)
    {
        if (!keyCorrect)
            appendValue(context.conformancyErrors,
                        std::string("Key incorrect, expected ") +
                        std::string(_programTagToString(TFormat())));

        if (queryLinePresent != 1)
            appendValue(context.conformancyErrors,
                        "No or multiple query lines present.");

        if (dbLinePresent != 1)
            appendValue(context.conformancyErrors,
                        "No or multiple database lines present.");

        if (fieldsLinePresent != 1)
            appendValue(context.conformancyErrors,
                        "No or multiple fields lines present.");

        if (length(context.otherLines) != 0)
            appendValue(context.conformancyErrors,
                        "Unexpected lines present, see context.otherLines.");
    }

    if ((g == BlastFormatGeneration::BLAST_PLUS) &&
        !((lastLine == 1) && atEnd(iter)))
    {
        if (!keyCorrect)
            appendValue(context.conformancyErrors,
                        std::string("Key incorrect, expected ") +
                        std::string(_programTagToString(TFormat())));

        if (queryLinePresent != 1)
            appendValue(context.conformancyErrors,
                        "No or multiple query lines present.");

        if (dbLinePresent != 1)
            appendValue(context.conformancyErrors,
                        "No or multiple database lines present.");
        // is ommitted in BLAST_PLUS when there are no hits
        if ((fieldsLinePresent != 1) && (length(r.matches) > 0))
            appendValue(context.conformancyErrors,
                        "No or multiple fields lines present.");

        if (length(context.otherLines) != 0)
            appendValue(context.conformancyErrors,
                        "Unexpected lines present, see context.otherLines.");
    }
}

template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readHeader(BlastRecord<TQId, TSId, TPos, TAlign> & r,
           TFwdIterator & iter,
           BlastIOContext & context,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    clear(r);
    clear(context.versionString);
    clear(context.dbSpecs);
    clear(context.otherLines);
    clear(context.conformancyErrors);
    clear(context.fieldsAsStrings);

    if (!context.ignoreFieldsInHeader)
        clear(context.fields);

    _readHeaderImpl(r, iter, context, TFormat());
}

template <typename TQId,
          typename TSId,
          typename TPos,
          typename TAlign,
          typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readHeader(BlastRecord<TQId, TSId, TPos, TAlign> &,
           TFwdIterator &,
           BlastIOContext &,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    // NOOP for TABULAR without header
}

// ----------------------------------------------------------------------------
// Function skipHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#skipHeader
 * @brief skip a header from a Blast tabular output file, optionally verifying it for format compliance.
 *
 * @signature int skipHeader(iter, [strict,] tag);
 *
 * @param[in,out] iter    An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     tag     @link BlastFormat @endlink tag, only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported.
 *
 * Call this function whenever you want to skip exactly one header. If you want
 * to go directly to the beginning of the next match (possibly skipping multiple
 * headers that have no succeeding matches) use skipUntilMatch instead.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @see BlastFormat#skipUntilMatch
 * @headerfile seqan/blast.h
 */

template <typename TFwdIterator,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
skipHeader(TFwdIterator & iter,
           BlastIOContext & context,
           BlastFormat<m, p, g> const & /*tag*/)
{
    BlastRecord<> r;
    readHeader(r, iter, context, BlastFormat<m, p, g>());
}

/*
 * @fn BlastFormat#readHeader
 * @brief read a Header from a Blast output file
 * @headerfile seqan/blast.h
 *
 * @signature readHeader(qId, dbName, versionString, [hits,] [fields, [otherLines,]] iter, strict, tag);
 *
 * @param[out]  qId         String to hold the query ID from the header
 * @param[out]  dbName      String to hold the database name from the header
 * @param[out]  versionString  String to hold the Blast program Tag and version
 * @param[out]  hits        Numerical to hold the number of hits that will follow the header (only available in BlastPlus spec of BlastFormat)
 * @param[out]  fields      StringSet to hold column identifiers, useful if non-defaults are expected
 * @param[out]  otherLines  StringSet to hold any comment or header lines that are not identified otherwise
 * @param[in,out] iter      An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context   An @link BlastIOContext @endlink with buffers.
 * @param[in]   strict      bool to signify whether the function should error on a non-conforming header or just "get whatever possible". If not using strict, it is recommended to pass fields and otherLines and verify these manually.
 * @param[in]   tag         @link BlastFormat @endlink tag, only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported.
 *
 * Call this function on every line beginning that is not "onMatch". Please note
 * that the strictness behaviour is slightly different, depending on whether you
 * specify a fieldList parameter or not. If you don't, strictness will check
 * whether the fields are default, if you do, strictness will not (because
 * you can check yourself!).
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 * @throw RecoverableParseError if the strictness condition is violated;
 * it is safe to catch this (and e.g. print some debug message) and then just
 * continue program operation.
 *
 * @see BlastFormat#onMatch
 */

// // default traditional Blast or BlastPlus
// template <typename TQId,
//           typename TDBName,
//           typename TVersionString,
//           typename TFwdIterator,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// readHeader(TQId & qId,
//            TDBName & dbName,
//            TVersionString & versionString,
//            TFwdIterator & iter,
//            BlastIOContext & context,
//            bool const strict,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
//
//     clear(context.buffers1);
//     clear(context.buffers2);
//     auto & otherLines   = context.buffers1;
//     auto & fields       = context.buffers2;
//
//     unsigned long hits = 0;
//
//     _readHeaderImplBlastTab(qId,
//                             dbName,
//                             versionString,
//                             fields,
//                             hits,
//                             otherLines,
//                             iter,
//                             context,
//                             strict,
//                             TFormat());
//
//     if (strict)
//         _verifyFields(context, fields, hits, TFormat());
// }
//
// // BlastPlus with hit count
// template <typename TQId,
//           typename TDBName,
//           typename TVersionString,
//           typename TFwdIterator,
//           BlastFormatProgram p>
// inline void
// readHeader(TQId & qId,
//            TDBName & dbName,
//            TVersionString & versionString,
//            unsigned long & hits,
//            TFwdIterator & iter,
//            BlastIOContext & context,
//            bool const strict,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
//                        p,
//                        BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
//                         p,
//                         BlastFormatGeneration::BLAST_PLUS> TFormat;
//
//     clear(context.buffers1);
//     clear(context.buffers2);
//     auto & otherLines   = context.buffers1;
//     auto & fields       = context.buffers2;
//
//     _readHeaderImplBlastTab(qId,
//                             dbName,
//                             versionString,
//                             fields,
//                             hits,
//                             otherLines,
//                             iter,
//                             context,
//                             strict,
//                             TFormat());
//
//     if (strict)
//          _verifyFields(context, fields, hits, TFormat());
// }
//
// // with fields
// template <typename TQId,
//           typename TDBName,
//           typename TVersionString,
//           typename TFieldList,
//           typename TFwdIterator,
//           BlastFormatProgram p,
//           BlastFormatGeneration g,
//           typename std::enable_if<
//             Is<ContainerConcept<TFieldList>>::VALUE, int>::type = 0>
// inline void
// readHeader(TQId & qId,
//            TDBName & dbName,
//            TVersionString & versionString,
//            TFieldList & fields,
//            TFwdIterator & iter,
//            BlastIOContext & context,
//            bool const strict,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
//
//     clear(context.buffers1);
//     auto & otherLines   = context.buffers1;
//
//     unsigned long hits = 0;
//
//     _readHeaderImplBlastTab(qId,
//                             dbName,
//                             versionString,
//                             fields,
//                             hits,
//                             otherLines,
//                             iter,
//                             context,
//                             strict,
//                             TFormat());
//
//     // don't verify fields, if user specified that he wants the list
//     // of fields, because that implies that he expects non-defaults
// }
//
// // with fields and hits count
// template <typename TQId,
//           typename TDBName,
//           typename TVersionString,
//           typename TFieldList,
//           typename TFwdIterator,
//           BlastFormatProgram p,
//           typename std::enable_if<
//             Is<ContainerConcept<TFieldList>>::VALUE, int>::type = 0>
// inline void
// readHeader(TQId & qId,
//            TDBName & dbName,
//            TVersionString & versionString,
//            unsigned long & hits,
//            TFieldList & fields,
//            TFwdIterator & iter,
//            BlastIOContext & context,
//            bool const  strict,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
//                        p,
//                        BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
//                         p,
//                         BlastFormatGeneration::BLAST_PLUS> TFormat;
//
//     clear(context.buffers1);
//     auto & otherLines   = context.buffers1;
//
//     _readHeaderImplBlastTab(qId,
//                             dbName,
//                             versionString,
//                             fields,
//                             hits,
//                             otherLines,
//                             iter,
//                             context,
//                             strict,
//                             TFormat());
//
//     // don't verify fields, if user specified that he wants the list
//     // of fields, because that implies that he expects non-defaults
// }
//
//
// // with fields and otherLines
// template <typename TQId,
//           typename TDBName,
//           typename TVersionString,
//           typename TFieldList,
//           typename TOtherString,
//           typename TFwdIterator,
//           BlastFormatProgram p,
//           BlastFormatGeneration g,
//           typename std::enable_if<
//             Is<ContainerConcept<TFieldList>>::VALUE, int>::type = 0>
// inline void
// readHeader(TQId & qId,
//            TDBName & dbName,
//            TVersionString & versionString,
//            TFieldList & fields,
//            StringSet<TOtherString> & otherLines,
//            TFwdIterator & iter,
//            BlastIOContext & context,
//            bool const strict,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
//
//     unsigned long hits = 0;
//
//     _readHeaderImplBlastTab(qId,
//                             dbName,
//                             versionString,
//                             fields,
//                             hits,
//                             otherLines,
//                             iter,
//                             context,
//                             strict,
//                             TFormat());
//
//     // don't verify fields, if user specified that he wants the list
//     // of fields, because that implies that he expects non-defaults
// }
//
// // with fields and otherLines and hits count
// template <typename TQId,
//           typename TDBName,
//           typename TVersionString,
//           typename TFieldList,
//           typename TOtherString,
//           typename TFwdIterator,
//           BlastFormatProgram p,
//           typename std::enable_if<
//             Is<ContainerConcept<TFieldList>>::VALUE, int>::type = 0>
// inline void
// readHeader(TQId & qId,
//            TDBName & dbName,
//            TVersionString & versionString,
//            unsigned long & hits,
//            TFieldList & fields,
//            StringSet<TOtherString> & otherLines,
//            TFwdIterator & iter,
//            BlastIOContext & context,
//            bool const strict,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
//                        p,
//                        BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
//                         p,
//                         BlastFormatGeneration::BLAST_PLUS> TFormat;
//
//     _readHeaderImplBlastTab(qId,
//                             dbName,
//                             versionString,
//                             fields,
//                             hits,
//                             otherLines,
//                             iter,
//                             context,
//                             strict,
//                             TFormat());
//
//     // don't verify fields, if user specified that he wants the list
//     // of fields, because that implies that he expects non-defaults
// }

// ----------------------------------------------------------------------------
// Function skipHeader()
// ----------------------------------------------------------------------------

/*
 * @fn BlastFormat#skipHeader
 * @brief skip a header from a Blast tabular output file, optionally verifying it for format compliance.
 *
 * @signature int skipHeader(iter, [strict,] tag);
 *
 * @param[in,out] iter    An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     strict  bool to signify whether the function should error on a non-conforming header or just "skip whatever possible".
 * @param[in]     tag     @link BlastFormat @endlink tag, only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported.
 *
 * Call this function whenever you want to skip exactly one header. If you want
 * to go directly to the beginning of the next match (possibly skipping multiple
 * headers that have no succeeding matches) use skipUntilMatch instead.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 * @throw RecoverableParseError if the strictness condition is violated;
 * it is safe to catch this (and e.g. print some debug message) and then just
 * continue program operation.
 *
 * @see BlastFormat#skipUntilMatch
 * @headerfile seqan/blast.h
 */

// template <typename TFwdIterator,
//           BlastFormatFile m,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// skipHeader(TFwdIterator & iter,
//            BlastIOContext & context,
//            bool const strict,
//            BlastFormat<m, p, g> const & /*tag*/)
// {
//     clear(context.buffer1);
//     clear(context.buffers1);
//     clear(context.buffers2);
//     auto & qId              = context.buffer1;
//     auto & dbName           = context.buffer1;
//     auto & versionString    = context.buffer1;
//     auto & fields           = context.buffers1;
//     auto & otherLines       = context.buffers2;
//
//     unsigned long hits = 0;
//
//     _readHeaderImplBlastTab(qId,
//                             dbName,
//                             versionString,
//                             fields,
//                             hits,
//                             otherLines,
//                             iter,
//                             context,
//                             strict,
//                             BlastFormat<m, p, g>());
// }
//
// template <typename TFwdIterator,
//           BlastFormatFile m,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// skipHeader(TFwdIterator & iter,
//            BlastIOContext & context,
//            BlastFormat<m, p, g> const & /*tag*/)
// {
//     skipHeader(iter, context, false, BlastFormat<m, p, g>());
// }

// ----------------------------------------------------------------------------
// Function skipUntilMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#skipUntilMatch
 * @brief skip arbitrary number of headers and/or comment lines until the beginning of a match is reached
 *
 * @signature skipUntilMatch(iter, tag);
 *
 * @param[in,out]   iter  An input iterator over a stream or any fwd-iterator over a string
 * @param[in]       tag   @link BlastFormat @endlink specialization,
 * with BlastFormatFile == BlastFormatFile::TABULAR || BlastFormatFile::TABULAR_WITH_HEADER
 *
 * Call this function whenever you are on a comment character ('#') in the file
 * and want to jump to the beginning of the next match. If you want to skip only
 * a single header (to count skipped headers or to verify its conformance
 * to standards), use BlastFormat#skipHeader instead.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @see BlastFormat#skipHeader
 * @headerfile seqan/blast.h
 */

template <typename TFwdIterator,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
skipUntilMatch(TFwdIterator & iter,
               BlastFormat<m, p, g> const & /*tag*/)
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
          typename TAlign>
inline void
_readField(BlastMatch<TQId, TSId, TPos, TAlign> & match,
           BlastIOContext & context,
           typename BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum
             const fieldId)
{
    typedef decltype(fieldId) ENUM;
    switch (fieldId)
    {
        case ENUM::STD: // this is cought in the calling function
            break;
        case ENUM::Q_SEQ_ID:
            match.qId = context.buffer2;
            break;
//         case ENUM::Q_GI: write(s,  * ); break;
//         case ENUM::Q_ACC: write(s,  * ); break;
//         case ENUM::Q_ACCVER: write(s,  * ); break;
        case ENUM::Q_LEN:
            match.qLength = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::S_SEQ_ID:
            match.sId = context.buffer2;
            break;
//         case ENUM::S_ALL_SEQ_ID: write(s,  * ); break;
//         case ENUM::S_GI: write(s,  * ); break;
//         case ENUM::S_ALL_GI: write(s,  * ); break;
//         case ENUM::S_ACC: write(s,  * ); break;
//         case ENUM::S_ACCVER: write(s,  * ); break;
//         case ENUM::S_ALLACC: write(s,  * ); break;
        case ENUM::S_LEN:
            match.sLength = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::Q_START:
            match.qStart = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::Q_END:
            match.qEnd = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::S_START:
            match.sStart = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::S_END:
            match.sEnd = lexicalCast<TPos>(context.buffer2);
            break;
//         case ENUM::Q_SEQ: write(s,  * ); break;
//         case ENUM::S_SEQ: write(s,  * ); break;
        case ENUM::E_VALUE:
            match.eValue = lexicalCast<double>(context.buffer2);
            break;
        case ENUM::BIT_SCORE:
            match.bitScore = lexicalCast<double>(context.buffer2);
            break;
        case ENUM::SCORE:
            match.score = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::LENGTH:
            match.aliLength = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::P_IDENT:
            // we don't have pIdent in the object, so instead we check if
            // N_IDENT is set already. If not we set N_IDENT to P_IDENT*100
            // and fix it later
//             std::cout << __LINE__ << "  " << match.identities << "  "
//                         << std::numeric_limits<TPos>::max()
//                         << "  " << _memberIsSet(match.identities) << "\n";
            if (!_memberIsSet(match.identities))
            {
                match.identities = lexicalCast<double>(context.buffer2) * 100;
//                 std::cout << __LINE__ << "  " << toCString(context.buffer2)
//                           << lexicalCast<double>(context.buffer2) << "  "
//                           << TPos(lexicalCast<double>(context.buffer2)) << "\n";
            }
            break;
        case ENUM::N_IDENT:
            match.identities = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::MISMATCH:
            match.mismatches = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::POSITIVE:
            match.positives = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::GAP_OPEN:
            match.gapOpenings = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::GAPS:
            match.gaps = lexicalCast<TPos>(context.buffer2);
            break;
        case ENUM::P_POS:
            // we don't have P_POS in the object, so instead we check if
            // POSITIVE is set already. If not we set POSITIVE to P_POS*100
            // and fix it later
            if (!_memberIsSet(match.positives))
                match.positives = lexicalCast<double>(context.buffer2) * 100;
        case ENUM::FRAMES:
        {
            clear(context.buffers2);
            strSplit(context.buffers2, context.buffer2, EqualsChar<'/'>());
            if (length(context.buffers2) != 2)
                SEQAN_THROW(ParseError("Could not process frame string."));
            match.qFrameShift = lexicalCast<int8_t>(context.buffers2[0]);
            match.sFrameShift = lexicalCast<int8_t>(context.buffers2[1]);
        } break;
        case ENUM::Q_FRAME:
            match.qFrameShift = lexicalCast<int8_t>(context.buffer2);
            break;
        case ENUM::S_FRAME:
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

template <typename TQId,
          typename TSId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          BlastFormatProgram    p>
inline void
_readMatchImpl(BlastMatch<TQId, TSId, TPos, TAlign> & match,
               TFwdIterator & iter,
               BlastIOContext & context,
               BlastFormat<BlastFormatFile::TABULAR,
                           p,
                           BlastFormatGeneration::BLAST_PLUS> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    typedef BlastMatchField<BlastFormatGeneration::BLAST_PLUS> TField;
    // header should have been read or skipped
    if (SEQAN_UNLIKELY(!onMatch(iter, TFormat())))
        SEQAN_THROW(ParseError("Not on beginning of Match (you should have "
                               "skipped comments)."));

    match._maxInitialize(); // mark all members as not set

    clear(context.buffer1);
    clear(context.buffers1);

    auto & line = context.buffer1;
    readLine(line, iter);

    auto & fields = context.buffers1;
    strSplit(fields, line, EqualsChar<'\t'>());

    int identIsPercent = -1;
    int posIsPercent = -1;

    unsigned n = 0u;
    for (typename TField::Enum const f : context.fields)
    {
        // this field represents multiple fields
        if (f == TField::Enum::STD)
        {
            for (typename TField::Enum const f2 : TField::defaults)
            {
                if (SEQAN_UNLIKELY(n >= length(fields)))
                    SEQAN_THROW(ParseError("More columns expected than were "
                                           "present in file."));

                context.buffer2 = static_cast<decltype(context.buffer2)>(fields[n++]);
                _readField(match, context, f2);
            }
            if (identIsPercent == -1)
                identIsPercent = 1;
        } else
        {
            if (SEQAN_UNLIKELY(n >= length(fields)))
                SEQAN_THROW(ParseError("More columns expected than were "
                                        "present in file."));

            context.buffer2 = static_cast<decltype(context.buffer2)>(fields[n++]);
            _readField(match, context, f);

            if ((f == TField::Enum::P_IDENT) && (identIsPercent == -1))
                identIsPercent = 1;
            else if (f == TField::Enum::N_IDENT)
                identIsPercent = 0;
            else if ((f == TField::Enum::P_POS) && (posIsPercent == -1))
                posIsPercent = 1;
            else if (f == TField::Enum::POSITIVE)
                posIsPercent = 0;
        }
    }

    // retransform the percentages to real numbers
    if (_memberIsSet(match.aliLength) && (identIsPercent == 1))
        match.identities = ROUND(double(match.aliLength * match.identities) / 10000);
    if (_memberIsSet(match.aliLength) && (posIsPercent == 1))
        match.positives = ROUND(double(match.aliLength * match.positives) / 10000);
    // and compute gaps from other values
    if (_memberIsSet(match.aliLength) && _memberIsSet(match.identities) &&
         _memberIsSet(match.mismatches) && !_memberIsSet(match.gaps))
        match.gaps = match.aliLength - match.mismatches - match.identities;
}

template <typename TQId,
          typename TSId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          BlastFormatProgram    p>
inline void
_readMatchImpl(BlastMatch<TQId, TSId, TPos, TAlign> & match,
               TFwdIterator & iter,
               BlastIOContext & context,
               BlastFormat<BlastFormatFile::TABULAR,
                           p,
                           BlastFormatGeneration::BLAST_LEGACY> const &)
{
    // BlastFormatGeneration::BLAST_LEGACY doesn't have individual fields
    // because only defaults are supported.
    // Since defaults are the same, we can overload the type here as BLAST_PLUS
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    typedef typename BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum TEnum;

    if (SEQAN_ENABLE_DEBUG)
    {
        if ((length(context.fields) != 1) || (context.fields[0] != TEnum::STD))
            std::cerr << "custom fields set, but will be ignored for "
                         "BlastFormatGeneration::BLAST_LEGACY\n";
    }

    // set defaults
    clear(context.fields);
    appendValue(context.fields, TEnum::STD);

    _readMatchImpl(match, iter, context, TFormat());

    // since gaps are included in the mismatch count in BLAST_LEGACY the gaps
    // computed here will always be zero, so we instead reset them to signify
    // that they are unknown
    match.gaps = std::numeric_limits<TPos>::max();
}

/*!
 * @fn BlastMatch#readMatch
 * @brief read a match from a file in BlastFormat
 *
 * @signature readMatch(blastMatch, iter, tag);
 *
 * @param[out]    blastMatch  A @link BlastMatch @endlink object to hold all relevant info
 * @param[in,out] iter        An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context     A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     tag         @link BlastFormat @endlink tag, only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported.
 *
 * This signature works with @link BlastMatch @endlinkes which is the recommended
 * way.
 *
 * The @link BlastIOContext::fields @endlink member of the context
 * can be specified if you expect a custom column composition. Specifying this
 * parameter implies that you know the columns are not
 * default and that they instead represent the given order and types.
 * You may specify less
 * fields than are actually present, in this case the additional fields will be
 * discarded. The parameter is not available when BlastFormatGeneration ==
 * BLAST_LEGACY.
 *
 * All members of the @link BlastMatch @endlink that are not set after the read
 * operation will be initialized to their respective max-values to signify this.
 *
 * Please note that the only transformations made to the data are the following:
 *
 *  <li> computation of the number of identities (from the percentage) [default]</li>
 *  <li> computation of the number of positives (from the percentage) [if given]</li>
 *  <li> number of gaps computed from other values [default]</li>
 *
 * In contrast to @link BlastMatch#writeMatch @endlink no other transformations
 * are made, e.g. the positions are still one-indexed and
 * flipped for reverse strand matches. This is due to the required fields for
 * retransformation (sequence lengths, frames) not being available in the
 * default columns.
 *
 * Instead of using this signature you may also use @link BlastFormat#readMatch
 * @endlink which works without a @link BlastMatch @endlink parameter and
 * supports arbitrary columns.
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
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readMatch(BlastMatch<TQId, TSId, TPos, TAlign> & match,
          TFwdIterator & iter,
          BlastIOContext & context,
          BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    clear(match);
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _readMatchImpl(match, iter, context, TFormat());
}

template <typename TQId,
          typename TSId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readMatch(BlastMatch<TQId, TSId, TPos, TAlign> & match,
          TFwdIterator & iter,
          BlastIOContext & context,
          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    readMatch(match, iter, context, TFormat());
}

// template <typename TQId,
//           typename TSId,
//           typename TFwdIterator,
//           typename TPos,
//           typename TAlign,
//           typename TFieldList,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// readMatch(BlastMatch<TQId, TSId, TPos, TAlign> & match,
//           TFwdIterator & iter,
//           BlastIOContext & context,
//           TFieldList const & fieldList,
//           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     static_assert(g == BlastFormatGeneration::BLAST_PLUS,
//                   "fieldList parameter may only be specified for "
//                   "BlastFormatGeneration::BLAST_PLUS");
//
//     typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
//     _readMatchImpl(match, iter, context, fieldList, TFormat());
// }
//
//
// // default fields
// template <typename TQId,
//           typename TSId,
//           typename TFwdIterator,
//           typename TPos,
//           typename TAlign,
//           BlastFormatProgram    p,
//           BlastFormatGeneration g>
// inline void
// readMatch(BlastMatch<TQId, TSId, TPos, TAlign> & match,
//           TFwdIterator & iter,
//           BlastIOContext & context,
//           BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
//     _readMatchImpl(match, iter, context, BlastMatchField<g>::defaults, TFormat());
// }
//
// // default fields
// template <typename TQId,
//           typename TSId,
//           typename TFwdIterator,
//           typename TPos,
//           typename TAlign,
//           BlastFormatProgram    p,
//           BlastFormatGeneration g>
// inline void
// readMatch(BlastMatch<TQId, TSId, TPos, TAlign> & match,
//           TFwdIterator & iter,
//           BlastIOContext & context,
//           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat; // same as above
//     _readMatchImpl(match, iter, context, BlastMatchField<g>::defaults, TFormat());
// }

/*!
 * @fn BlastFormat#readMatch
 * @brief read arbitrary columns from a file in BlastFormat
 *
 * @signature readMatch(iter, tag, args ...);
 *
 * @param[in,out] iter    An input iterator over a stream or any fwd-iterator over a string
 * @param[in]     tag     Only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported
 * @param[out]    args    Arbitrary typed variables
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
    std::string buffer;
    SEQAN_TRY
    {
        readUntil(buffer, iter, OrFunctor<IsTab,IsNewline>());
        _assignOrCast(arg, buffer);
    } SEQAN_CATCH (UnexpectedEnd const &)
    {
        return;
    }
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
    std::string buffer;
    readUntil(buffer, iter, IsTab());
    skipOne(iter, IsTab());
    _assignOrCast(arg, buffer);

    // recurse to next argument
    _readMatchImplBlastTab(iter, args...);
}

// custom arguments
template <typename TFwdIterator,
          typename... TArgs,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline void
readMatch(TFwdIterator & iter,
          BlastFormat<BlastFormatFile::TABULAR, p, g> const &,
          TArgs & ... args)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    // header should have been read or skipped
    if (SEQAN_UNLIKELY(!onMatch(iter, TFormat())))
        SEQAN_THROW(ParseError("Not on beginning of Match (you should have "
                               "skipped comments)."));

    _readMatchImplBlastTab(iter, args...);
}

template <typename TFwdIterator,
          typename ... TArgs,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline void
readMatch(TFwdIterator & iter,
          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &,
          TArgs & ... args)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat; // same as above
    readMatch(iter, TFormat(), args...);
}

// ----------------------------------------------------------------------------
// Function skipMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastMatch#skipMatch
 * @brief skip a line that contains a match
 *
 * @signature skipMatch(iter, context, tag);
 *
 * @param[in,out] iter    An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     tag     @link BlastFormat @endlink specialization,
 * with BlastFormatFile == BlastFormatFile::TABULAR || BlastFormatFile::TABULAR_WITH_HEADER
 *
 * This function always checks whether it is on a line that is not a comment and
 * it does verify the line read as looking like a match. If you do not want this
 * behaviour you can instead <texttt>skipLine(iter)</texttt> which is also
 * faster.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @see BlastFormat#skipHeader
 * @headerfile seqan/blast.h
 */

// template <typename TFwdIterator,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// skipMatch(TFwdIterator & iter,
//           BlastIOContext &,
//           BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
// {
//     typedef  BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
//     // header should have been read or skipped
//     if (SEQAN_UNLIKELY(!onMatch(iter, TFormat())))
//         SEQAN_THROW(ParseError("Not on beginning of Match (you should have "
//                                "skipped comments)."));
//     skipLine(iter);
// }
//
// template <typename TFwdIterator,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// skipMatch(TFwdIterator & iter,
//           BlastIOContext & context,
//           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     typedef  BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
//     skipMatch(iter, context, TFormat());
// }

template <typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
skipMatch(TFwdIterator & iter,
          BlastIOContext & context,
          BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    typedef  BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;

    BlastMatch<> m;
    readMatch(m, iter, context, TFormat());
}

template <typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
skipMatch(TFwdIterator & iter,
          BlastIOContext & context,
          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef  BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    skipMatch(iter, context, TFormat());
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastRecord#readRecord
 * @brief read a record from a file in BlastFormat
 *
 * @signature readRecord(blastRecord, iter, context, tag);
 *
 * @param[out]    blastRecord A @link BlastRecord @endlink to hold all relevant info
 * @param[in,out] iter        An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context     A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     tag         @link BlastFormat @endlink specialization
 *
 * This function will read an entire record from a blast (tabular) file, i.e.
 * a @link BlastRecord @endlink object including 0-n @link BlastMatch @endlink
 * objects.
 *
 * In the case of the TABULAR_WITH_HEADER format, this function also sets the
 * following members of the @link BlastIOContext @endlink
 * which may or may not intereset you:
 * <li> @link BlastDbSpecs::dbName @endlink member of
 * @link BlastIOContext::dbSpecs @endlink (name of the database)</li>
 * <li> @link BlastIOContext::headerConformsToStandards @endlink a bool
 * signafying whether the header contained no anomalies.</li>
 * <li> @link BlastIOContext::fields @endlink: descriptors for the columns</li>
 * <li> @link BlastIOContext::fieldsAsStrings @endlink: labels for the columns
 * as they appear in the file</li>
 * <li> @link BlastIOContext::otherLines @endlink: any unknown lines</li>
 *
 * In @link BlastIOContext::fields @endlink member is read from the header and
 * used as in-parameter for reading the matches, i.e. it reads the column
 * labels and after that expects the values in the columns to correspond
 * to this. If you do not want this behaviour, because you expect the header
 * to be non-standard -- as is the case with many tools -- you can set
 * context.@link BlastIOContext::ignoreFieldsInHeader @endlink to true.
 * In this case the you can also set
 * context.@link BlastIOContext::fields @endlink manually and these columns
 * will be expected for the matches, independent of the header. The latter
 * behaviour is also default for the TABULAR format without headers.
 *
 * For @link BlastFormatGeneration @endlink::BLAST_LEGACY the fields member is
 * always ignored, however @link BlastIOContext::fieldsAsStrings @endlink is
 * still read from the header.
 *
 * Please note that for the TABULAR format (without headers) the boundary
 * between records is inferred from the indentity of the first field, i.e.
 * custom columns must also have Q_SEQ_ID as first field. Also you must pass
 * an std::string buffer as in-out parameter to cache this field.
 * @link BlastMatch#readMatch @endlink.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @headerfile seqan/blast.h
 */

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TAlign,
          typename TPos,
          BlastFormatProgram p>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
           TFwdIterator & iter,
           BlastIOContext & context,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_PLUS> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_PLUS> TFormat;

    readHeader(blastRecord, iter, context, TFormat());

    // .matches already resized for us
    for (auto & m : blastRecord.matches)
    {
        if (!onMatch(iter, TFormat()))
        {
            appendValue(context.conformancyErrors,
                        "Less matches than promised by header");
            break;
        }

        readMatch(m, iter, context, TFormat());
    }

    if ((!atEnd(iter)) && onMatch(iter, TFormat()))
        appendValue(context.conformancyErrors,
                    "More matches than promised by header");

    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
    {
        blastRecord.matches.emplace_back();
        readMatch(back(blastRecord.matches), iter, context, TFormat());
    }
}

template <typename TFwdIterator,
          typename TQId,
          typename TSId,
          typename TAlign,
          typename TPos,
          BlastFormatProgram p>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
           TFwdIterator & iter,
           BlastIOContext & context,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_LEGACY> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_LEGACY> TFormat;

    readHeader(blastRecord, iter, context, TFormat());

    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
    {
        blastRecord.matches.emplace_back();
        readMatch(back(blastRecord.matches), iter, context, TFormat());
    }
}

// TABULAR WITHOUT HEADER
template <typename TQId,
          typename TSId,
          typename TAlign,
          typename TPos,
          typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
           TFwdIterator & iter,
           BlastIOContext & context,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        BlastFormatGeneration::BLAST_PLUS>  TFormat;
    typedef BlastMatchField<BlastFormatGeneration::BLAST_PLUS> TField;
    clear(blastRecord);

    auto curId = context.lastId;

    std::vector<typename TField::Enum> fieldListMinusFirst;
    if (context.fields[0] == TField::Enum::STD)
    {
        fieldListMinusFirst = PrivateSpace_<>::defaultsMinusQId;
        if (length(context.fields) > 1)
            fieldListMinusFirst.insert(std::end(fieldListMinusFirst),
                                       std::begin(context.fields) + 1,
                                       std::end(context.fields));
    }
    else if (context.fields[0] == TField::Enum::Q_SEQ_ID)
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

    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
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
        readMatch(back(blastRecord.matches), iter, context, TFormat());
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

/*
 * @fn BlastRecord#readRecord
 * @brief read a record from a file in BlastFormat
 *
 * @signature readRecord(blastRecord, dbName, [fieldList,] iter, [fieldList,] tag); [WITH_HEADER]
 * readRecord(blastRecord, buffer, iter, [fieldList,] tag); [WITHOUT_HEADER]
 *
 * @param[out]    blastRecord A @link BlastRecord @endlink to hold all relevant info
 * @param[out]    dbName      String to hold the database name from the header
 * @param[out]    fieldList   Read the @link BlastMatchField @endlinks from the header
 * (only TABULAR_WITH_HEADER)
 * @param[in,out] iter        An input iterator over a stream or any fwd-iterator over a string
 * @param[in,out] context     A @link BlastIOContext @endlink with parameters and buffers.
 * @param[in]     fieldList   Expect this fieldList in the matches.
 * @param[in]     tag         @link BlastFormat @endlink specialization
 *
 * The fieldList parameter can indeed be specified as either an in-parameter
 * (both TABULAR and TABULAR_WITH_HEADER) or an out-parameter (only
 * TABULAR_WITH_HEADER). Using it as an in-parameter indicates that you know
 * that the columns are the specified composition and you don't care whatever
 * is written in the header. Using it as an out-parameter will result in the
 * fieldList being read from the header and then used as in-parameter for reading
 * the matches (and being returned to you via the parameter).
 *
 * Both fieldList parameters can only be used with
 * @link BlastFormatGeneration @endlink::BLAST_PLUS and the out-parameter is
 * only available with @link BlastFormatFile @endlink::TABULAR_WITH_HEADER.
 *
 * Please note that for the TABULAR format (without headers) the boundary
 * between records is inferred from the indentity of the first field, i.e.
 * custom columns must also have Q_SEQ_ID as first field. Also you must pass
 * an std::string buffer as in-out parameter to cache this field.
 * @link BlastMatch#readMatch @endlink.
 *
 * @throw IOError On low-level I/O errors.
 * @throw ParseError On high-level file format errors.
 *
 * @headerfile seqan/blast.h
 */

// template <typename TFwdIterator,
//           typename TDbName,
//           typename TQId,
//           typename TSId,
//           typename TAlign,
//           typename TPos,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
//            TDbName & dbName,
//            TFwdIterator & iter,
//            BlastIOContext & context,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
//
//     clear(blastRecord);
//     clear(dbName);
//
//     clear(context.buffer1);
//     auto & versionString = context.buffer1;
//
//     readHeader(blastRecord.qId,
//                dbName,
//                versionString,
//                iter,
//                context,
//                false,
//                TFormat());
//
//     while ((!atEnd(iter)) && onMatch(iter, TFormat()))
//     {
//         blastRecord.matches.emplace_back();
//         readMatch(back(blastRecord.matches), iter, context, TFormat());
//     }
// }
//
// // custom fields
// template <typename TFwdIterator,
//           typename TDbName,
//           typename TQId,
//           typename TSId,
//           typename TAlign,
//           typename TPos,
//           typename TFieldList,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
//            TDbName & dbName,
//            TFwdIterator & iter,
//            BlastIOContext & context,
//            TFieldList const & fieldList,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     static_assert(g == BlastFormatGeneration::BLAST_PLUS,
//                   "fieldList parameter may only be specified for "
//                   "BlastFormatGeneration::BLAST_PLUS");
//
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
//
//     clear(blastRecord);
//     clear(dbName);
//
//     clear(context.buffer1);
//     auto & versionString = context.buffer1;
//
//     readHeader(blastRecord.qId,
//                dbName,
//                versionString,
//                iter,
//                context,
//                false,
//                TFormat());
//
//     while ((!atEnd(iter)) && onMatch(iter, TFormat()))
//     {
//         blastRecord.matches.emplace_back();
//         //                              only difference to above ↓
//         readMatch(back(blastRecord.matches), iter, context, fieldList,
//                   TFormat());
//     }
// }
//
// // detect fields
// template <typename TFwdIterator,
//           typename TDbName,
//           typename TQId,
//           typename TSId,
//           typename TAlign,
//           typename TPos,
//           typename TFieldList,
//           BlastFormatProgram p,
//           BlastFormatGeneration g,
//           typename std::enable_if<IsSequence<TFieldList>::VALUE, int>::type = 0>
// inline void
// readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
//            TDbName & dbName,
//            TFieldList & fieldList,
//            TFwdIterator & iter,
//            BlastIOContext & context,
//            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
// {
//     static_assert(g == BlastFormatGeneration::BLAST_PLUS,
//                   "fieldList parameter may only be specified for "
//                   "BlastFormatGeneration::BLAST_PLUS");
//
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
//
//     clear(blastRecord);
//     clear(dbName);
//
//     clear(context.buffer1);
//     auto & versionString = context.buffer1;
//
//     readHeader(blastRecord.qId,
//                dbName,
//                versionString,
//                fieldList, // out-parameter to readHeader()
//                iter,
//                context,
//                false,
//                TFormat());
//
//     while ((!atEnd(iter)) && onMatch(iter, TFormat()))
//     {
//         blastRecord.matches.emplace_back();
//         //                             in-parameter to readMatch ↓
//         readMatch(back(blastRecord.matches), iter, context, fieldList,
//                   TFormat());
//     }
// }

// TABULAR FORMAT has no headers, so we compare first query id
// template <typename TFwdIterator,
//           typename TQId = CharString,
//           typename TSId = CharString,
//           typename TAlign = Align<CharString, ArrayGaps>,
//           typename TPos = unsigned,
//           BlastFormatProgram p,
//           BlastFormatGeneration g,
//           typename std::enable_if<
//             Is<ReversibleContainerConcept<TFwdIterator> >::VALUE,
//             int> = 0>
// inline void
// readRecord(BlastRecord<TQId, TSId, TPos, TAlign> & blastRecord,
//            TFwdIterator & iter,
//            BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
// {
//     typedef BlastFormat<BlastFormatFile::TABULAR, p, g>  TFormat;
//
//     clear(blastRecord);
//
//     std::string curId;
//     std::string lastId;
//
//     //TODO for custom fields we have to change this, i.e. read entire line
//     //and see if query id is contained
//     while ((!atEnd(iter)) && onMatch(iter, TFormat()))
//     {
//         clear(curId);
//         // read current read ID and then rewind stream to beginning of line
//         readUntil(curId, iter, OrFunctor<IsTab,IsNewline>());
//         iter -= length(curId);
//
//         if ((curId != lastId) && (!empty(lastId)))
//             break; // new Record reached
//
//         blastRecord.matches.emplace_back();
//         readMatch(back(blastRecord.matches), iter, TFormat());
//         lastId = curId;
//     }
//
//     if (length(blastRecord.matches) == 0)
//         SEQAN_THROW(ParseError("No Matches could be read."));
//
//     blastRecord.qId = blastRecord.matches.front().qId;
// }

// same as above but lastId comes from outside


} // namespace seqan

#endif // SEQAN_BLAST_READ_BLAST_TABULAR_H_
