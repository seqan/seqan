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

HEADER
 RECORD
 RECORD
 RECORD
 RECORD
HEADER
 RECORD
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

Because 0 records are allowed, multiple Headers can succeed each other, the criterium for seperation employed by this implementation is that an NCBI Blast
records always ends after the "Fields" line and NCBI Blast+ records end after
the "number of hits"-line.
A file is considered NCBI Blast+ format, when it has comment line that starts
with "# BLAST" and ends with "+". In all other cases it is considered
traditional Blast format.
*/

#define BAD_FORMAT "File is non-standard and could not be read."

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct DetectFields_;
typedef Tag<DetectFields_> DetectFields;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/*!
 * @fn BlastFormat#onMatch
 * @brief returns whether RecordReader is on beginning of match
 *
 * @signature bool onMatch(iter, tag)
 *
 * @param[in] iter    FwdIterator or Stream
 * @param[in] tag     BlastFormat specialization
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

// fields as string of blastmatchfield::enum
template <typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename std::enable_if<IsSequence<TFieldList>::VALUE>::type = 0,
          typename std::enable_if<
            std::is_same<typename Value<TFieldList>::Type,
                         typename BlastMatchField<g>::Enum>::value>::type = 0>
inline void
_verifyFields(TFieldList const & fields,
              unsigned const hits, // irrelevant for traditional header
              BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g>  TFormat;

    if (length(fields) != BlastMatchField<g>::defaults)
        throw std::ios_base::failure(BAD_FORMAT);

    for (unsigned i = 0; i < length(fields); ++i)
        if (fields[i] != BlastMatchField<g>::defaults[i])
            throw std::ios_base::failure(BAD_FORMAT);
}


// fields as strings
template <typename TString,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_verifyFields(StringSet<TString> const & fields,
              unsigned const   hits, // irrelevant for traditional header
              BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g>  TFormat;

//     // TODO(h4nn3s): what was this supposed to do?
//     if (g == BlastFormatGeneration::BLAST_PLUS)
//         if ((hits == 0) && length(fields) )
//             0;

    CharString fieldStr;
    joinStringSet(fieldStr, fields, _seperatorString(TFormat()));

    // will always be zero but this is cleaner
    const uint8_t std = static_cast<uint8_t>(BlastMatchField<g>::Enum::STD);

    // compare the first twelve fields
    if (prefix(fieldStr, length(BlastMatchField<g>::columnLabels[std])) !=
        BlastMatchField<g>::columnLabels[std])
        throw std::ios_base::failure(BAD_FORMAT);
}

// ----------------------------------------------------------------------------
// Function readHeader()
// ----------------------------------------------------------------------------

template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFwdIterator,
          typename TString,
          typename TString2,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_readHeaderImplBlastTab(TqId                                    & qId,
                        TDBName                                 & dbName,
                        TVersionString                          & versionString,
                        StringSet<TString>                      & fields,
                        unsigned long                           & hits,
                        StringSet<TString2>                     & otherLines,
                        // any other lines
                        TFwdIterator & iter,
                        bool                              const   strict,
                        BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;

    // this is a record instead of a header
    if (onMatch(iter, TFormat()))
        throw std::ios_base::failure(BAD_FORMAT);

    bool keyCorrect = false;
    int queryLinePresent = 0;
    int dbLinePresent = 0;
    int fieldsLinePresent = 0;
    int hitsLinePresent = 0;
    int lastLine = 0;

    while ((!atEnd(iter)) && (!onMatch(iter, TFormat())))// in Header
    {
        // skip '#'
        goNext(iter);
        //skip blanks following the '#'
        skipBlanks(iter);

        CharString key = "";
        readUntilWhitespace(key, iter);

        if (key == _programTagToString(TFormat()))
        {
            // get whole line
            CharString buf;
            readLine(buf, iter);
            CharString &fullLine = key;
            append(fullLine, buf);

            versionString = fullLine;
            keyCorrect = true;
        }
        else if (key == "Query:")
        {
            skipBlanks(iter);
            readLine(qId, iter);
            ++queryLinePresent;
        }
        else if (key == "Database:")
        {
            skipBlanks(iter);
            readLine(dbName, iter);
            ++dbLinePresent;
        }
        else if (key == "Fields:")
        {
            skipBlanks(iter);

            CharString buf;
            readLine(buf, iter);
            strSplit(fields, buf, ", ");

            ++fieldsLinePresent;
            if (g == BlastFormatGeneration::BLAST_LEGACY)
                break; // header is finished
        }
        else
        {
            // get whole line
            CharString buf;
            readLine(buf, iter);
            CharString &fullLine = key;
            append(fullLine, buf);

            if (g == BlastFormatGeneration::BLAST_PLUS)
            {
                // last line of BlastPlus Format
                if (hasPrefix(fullLine, "BLAST processed"))
                {
                    ++lastLine;
                }
                // is hits counter?
                else if (endsWith(fullLine, "hits found"))
                {
                    CharString hitsString = "";
                    for (unsigned i = 0;
                         (i < length(fullLine) && isdigit(fullLine[i]));
                         ++i)
                        append(hitsString, fullLine[i], Generous());

                    hits = lexicalCast<decltype(hits)>(hitsString);
//                     if (ret && strict)
//                         throw std::ios_base::failure(BAD_FORMAT);

                    ++hitsLinePresent;
                    break; // header is finished
                }
            }
            else
            {
                appendValue(otherLines, fullLine, Generous());
            }
        }
    }

    if (!strict)
        return;

    if (g == BlastFormatGeneration::BLAST_LEGACY)
        if (  keyCorrect                &&
             (queryLinePresent   == 1)  &&
             (dbLinePresent      == 1)  &&
             (fieldsLinePresent  == 1)  &&
             (length(otherLines) == 0)   )
            return;

    if (g == BlastFormatGeneration::BLAST_PLUS)
        if (( keyCorrect                                   &&
             (queryLinePresent    == 1)                    &&
             (dbLinePresent       == 1)                    &&
             (hitsLinePresent     == 1)                    &&
             ( (fieldsLinePresent == 1) || (hits==0) )     &&
             (length(otherLines)  == 0)                     ) ||
             ( (lastLine          == 1) && atEnd(iter))   )
            return;

    throw std::ios_base::failure(BAD_FORMAT);
}

// fieldList is actually a list of fields and not a string
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFwdIterator,
          typename TString,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename std::enable_if<IsSequence<TFieldList>::VALUE>::type = 0,
          typename std::enable_if<
            std::is_same<typename Value<TFieldList>::Type,
                         typename BlastMatchField<g>::Enum>::value>::type = 0>
inline void
_readHeaderImplBlastTab(TqId                                    & qId,
                        TDBName                                 & dbName,
                        TVersionString                          & versionString,
                        TFieldList                              & fields,
                        unsigned long                           & hits,
                        StringSet<TString>                      & otherLines,
                        // any other lines
                        TFwdIterator & iter,
                        bool                              const   strict,
                        BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                    p,
                                    BlastFormatGeneration::BLAST_PLUS> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;
    typedef BlastMatchField<BlastFormatGeneration::BLAST_PLUS> TMatchField;

    StringSet<CharString> fieldStrings;

    _readHeaderImplBlastTab(qId,
                            dbName,
                            versionString,
                            fieldStrings,
                            hits,
                            otherLines,
                            iter,
                            strict,
                            TFormat());

    for (CharString const & s : fieldStrings)
    {
        for (uint8_t i = 0; i < length(TMatchField::columnLabels); ++i)
        {
            if (s == TMatchField::columnLabels[i])
            {
                appendValue(fields, static_cast<TMatchField::Enum>(i));
                break;
            }
        }
    }
}

/*!
 * @fn BlastFormat#readHeader
 * @brief read a Header from a Blast output file
 *
 * @signature int readHeader(qId, dbName, versionString, [hits,] [fields, [otherLines,]] iterator/stream, strict, tag);
 *
 * @param[out]  qId     String to hold the query ID from the header
 * @param[out]  dbName  String to hold the database name from the header
 * @param[out]  versionString  String to hold the Blast program Tag and version
 * @param[out]  hits    Numerical to hold the number of hits that will follow the header (only available in BlastPlus spec of BlastFormat)
 * @param[out]  fields  StringSet to hold column identifiers, useful if non-defaults are expected
 * @param[out]  otherLines  StringSet to hold any comment or header lines that are not identified otherwise
 * @param[in,out]   iter  RecordReader
 * @param[in]   strict  bool to signify whether the function should error on a non-conforming header or just "get whatever possible". If not using strict, it is recommended to pass fields and otherLines and verify these manually.
 * @param[in]   tag     @link BlastFormat @endlink tag, only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported.
 *
 * @see BlastFormat#onMatch
 * @headerfile seqan/blast.h
 * @section Remarks
 *
 * call this function on every line beginning that is not "onMatch"
 */

// default traditional Blast or BlastPlus
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readHeader(TqId                                                 & qId,
           TDBName                                              & dbName,
           TVersionString                                       & versionString,
           TFwdIterator & iter,
           bool                                           const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;

    StringSet<CharString> otherLines;
    StringSet<CharString> fields;
    unsigned long hits = 0;

    _readHeaderImplBlastTab(qId,
                            dbName,
                            versionString,
                            fields,
                            hits,
                            otherLines,
                            iter,
                            strict,
                            TFormat());

    if (strict)
        _verifyFields(fields, hits, TFormat());
}

// BlastPlus with hit count
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFwdIterator,
          BlastFormatProgram p>
inline void
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           unsigned long                                    & hits,
           TFwdIterator & iter,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    StringSet<CharString> otherLines;
    StringSet<CharString> fields;

    _readHeaderImplBlastTab(qId,
                            dbName,
                            versionString,
                            fields,
                            hits,
                            otherLines,
                            iter,
                            strict,
                            TFormat());

    if (strict)
         _verifyFields(fields, hits, TFormat());
}


// with fields
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldList,
          typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           TFieldList                                       & fields,
           TFwdIterator & iter,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;

    StringSet<CharString> otherLines;
    unsigned long hits = 0;

    _readHeaderImplBlastTab(qId,
                            dbName,
                            versionString,
                            fields,
                            hits,
                            otherLines,
                            iter,
                            strict,
                            TFormat());

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
}

// with fields and hits count
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldList,
          typename TFwdIterator,
          BlastFormatProgram p>
inline void
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           unsigned long                                    & hits,
           TFieldList                                       & fields,
           TFwdIterator & iter,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    StringSet<CharString> otherLines;

    _readHeaderImplBlastTab(qId,
                            dbName,
                            versionString,
                            fields,
                            hits,
                            otherLines,
                            iter,
                            strict,
                            TFormat());

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
}


// with fields and otherLines
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldList,
          typename TOtherString,
          typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           TFieldList                                       & fields,
           StringSet<TOtherString>                          & otherLines,
           TFwdIterator & iter,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;

    unsigned long hits = 0;

    _readHeaderImplBlastTab(qId,
                            dbName,
                            versionString,
                            fields,
                            hits,
                            otherLines,
                            iter,
                            strict,
                            TFormat());

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
}

// with fields and otherLines and hits count
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldList,
          typename TOtherString,
          typename TFwdIterator,
          BlastFormatProgram p>
inline void
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           unsigned long                                    & hits,
           TFieldList                                       & fields,
           StringSet<TOtherString>                          & otherLines,
           TFwdIterator & iter,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    _readHeaderImplBlastTab(qId,
                            dbName,
                            versionString,
                            fields,
                            hits,
                            otherLines,
                            iter,
                            strict,
                            TFormat());

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
}

// ----------------------------------------------------------------------------
// Function skipHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#skipHeader
 * @brief skip a header from a Blast tabular output file, optionally verifying it for format compliance.
 *
 * @signature int skipHeader(iterator/stream, [strict,] tag);
 *
 * @param[in,out]   iter  RecordReader
 * @param[in]   strict  bool to signify whether the function should error on a non-conforming header or just "skip whatever possible".
 * @param[in]   tag     @link BlastFormat @endlink tag, only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported.
 *
 * @see BlastFormat#skipUntilMatch
 * @    0 on success, and non-zero otherwise
 * @headerfile seqan/blast.h
 * @section Remarks
 *
 * call this function whenever you want to skip exactly one header. If you want
 * to go directly to the beginning of the next match (possibly skipping multiple
 * headers that have no succeeding matches) use skipUntilMatch instead.
 */

template <typename TFwdIterator,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
skipHeader(TFwdIterator & iter,
           bool                         const   strict,
           BlastFormat<m, p, g>           const & /*tag*/)
{
    CharString qId;
    CharString dbName;
    CharString versionString;
    StringSet<CharString> fields;
    unsigned long hits = 0;
    StringSet<CharString> otherLines;

    _readHeaderImplBlastTab(qId,
                            dbName,
                            versionString,
                            fields,
                            hits,
                            otherLines,
                            iter,
                            strict,
                            BlastFormat<m, p, g>());
}

template <typename TFwdIterator,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
skipHeader(TFwdIterator & iter,
           BlastFormat<m, p, g> const & /*tag*/)
{
    skipHeader(iter, false, BlastFormat<m, p, g>());
}

// ----------------------------------------------------------------------------
// Function skipUntilMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#skipUntilMatch
 * @brief skip arbitrary number of headers and/or comment lines until the beginning of a match is reached
 *
 * @signature skipUntilMatch(iterator/stream, tag);
 *
 * @param[in,out]   iter  RecordReader
 * @param[in]       tag     BlastFormat specialization,
 * with BlastFormat::_m == BlastFormatFile::TABULAR || BlastFormatFile::TABULAR_WITH_HEADER
 *
 * @    0 on success, and non-zero otherwise
 * @see BlastFormat#skipHeader
 * @headerfile seqan/blast.h
 * @section Remarks
 *
 * call this function whenever you are on a comment character ('#') in the file and want to jump to the beginning of the next record. If you want skip only a single header (to count skipped headers or verify its conformance
to standards), use BlastFormat#skipHeader instead.
 */

template <typename TFwdIterator,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
skipUntilMatch(TFwdIterator & iter,
               BlastFormat<m, p, g>      const & /*tag*/)
{
    while ((!atEnd(iter)) && value(iter) == '#') // skip comments
        skipLine(iter);
}

// ----------------------------------------------------------------------------
// Function readMatch()
// ----------------------------------------------------------------------------

template <typename TqId,
          typename TsId,
          typename TPos,
          typename TAlign,
          BlastFormatGeneration g>
inline void
_readField(BlastMatch<TqId, TsId, TAlign, TPos> & match,
           CharString const & buffer,
           typename BlastMatchField<g>::Enum const fieldId)
{
    typedef typename BlastMatchField<g>::Enum ENUM;
    switch (fieldId)
    {
        case ENUM::STD: // this is cought in the calling function
            break;
        case ENUM::Q_SEQ_ID:
            match.qId = buffer;
            break;
//         case ENUM::Q_GI: write(s,  * ); break;
//         case ENUM::Q_ACC: write(s,  * ); break;
//         case ENUM::Q_ACCVER: write(s,  * ); break;
        case ENUM::Q_LEN:
            match.qLength = lexicalCast<TPos>(buffer);
            break;
        case ENUM::S_SEQ_ID:
            match.sId = buffer;
            break;
//         case ENUM::S_ALL_SEQ_ID: write(s,  * ); break;
//         case ENUM::S_GI: write(s,  * ); break;
//         case ENUM::S_ALL_GI: write(s,  * ); break;
//         case ENUM::S_ACC: write(s,  * ); break;
//         case ENUM::S_ACCVER: write(s,  * ); break;
//         case ENUM::S_ALLACC: write(s,  * ); break;
        case ENUM::S_LEN:
            match.sLength = lexicalCast<TPos>(buffer);
            break;
        case ENUM::Q_START:
            match.qStart = lexicalCast<TPos>(buffer);
            break;
        case ENUM::Q_END:
            match.qEnd = lexicalCast<TPos>(buffer);
            break;
        case ENUM::S_START:
            match.sStart = lexicalCast<TPos>(buffer);
            break;
        case ENUM::S_END:
            match.sEnd = lexicalCast<TPos>(buffer);
            break;
//         case ENUM::Q_SEQ: write(s,  * ); break;
//         case ENUM::S_SEQ: write(s,  * ); break;
        case ENUM::E_VALUE:
            match.eValue = lexicalCast<double>(buffer);
            break;
        case ENUM::BIT_SCORE:
            match.bitScore = lexicalCast<double>(buffer);
            break;
        case ENUM::SCORE:
            match.score = lexicalCast<TPos>(buffer);
            break;
        case ENUM::LENGTH:
            match.aliLength = lexicalCast<TPos>(buffer);
            break;
        case ENUM::P_IDENT:
            // we don't have pIdent in the object, so instead we check if
            // N_IDENT is set already. If not we set N_IDENT to P_IDENT*100
            // and fix it later
            if (!_memberIsSet(match.identities))
                match.identities = lexicalCast<TPos>(buffer) * 100;
            break;
        case ENUM::N_IDENT:
            match.identities = lexicalCast<TPos>(buffer);
            break;
        case ENUM::MISMATCH:
            match.mismatches = lexicalCast<TPos>(buffer);
            break;
        case ENUM::POSITIVE:
            match.positives = lexicalCast<TPos>(buffer);
            break;
        case ENUM::GAP_OPEN:
            match.gapOpenings = lexicalCast<TPos>(buffer);
            break;
        case ENUM::GAPS:
            match.gaps = lexicalCast<TPos>(buffer);
            break;
        case ENUM::P_POS:
            // we don't have P_POS in the object, so instead we check if
            // POSITIVE is set already. If not we set POSITIVE to P_POS*100
            // and fix it later
            if (!_memberIsSet(match.positives))
                match.positives = lexicalCast<TPos>(buffer) * 100;
        case ENUM::FRAMES:
        {
            StringSet<CharString> buffers;
            strSplit(buffers, buffer, '/');
            if (length(buffers) != 2)
                throw std::ios_base::failure(BAD_FORMAT);
            match.qFrameShift = lexicalCast<int8_t>(buffers[0]);
            match.sFrameShift = lexicalCast<int8_t>(buffers[1]);
        } break;
        case ENUM::Q_FRAME:
            match.qFrameShift = lexicalCast<int8_t>(buffer);
            break;
        case ENUM::S_FRAME:
            match.sFrameShift = lexicalCast<int8_t>(buffer);
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
    };
}

template <typename TqId,
          typename TsId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          typename TFieldList,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline void
_readMatchImpl(BlastMatch<TqId, TsId, TAlign, TPos>      & match,
               TFwdIterator & iter,
               TFieldList const & fieldList,
               BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    // header should have been read or skipped
    if (SEQAN_UNLIKELY(!onMatch(iter, TFormat())))
        throw std::ios_base::failure(BAD_FORMAT);

    match._maxInitialize(); // mark all members as not set

    CharString line;
    readLine(line, iter);

    StringSet<CharString> fields;
    strSplit(fields, line, _seperatorString(TFormat()) /* == '\t' */);

    bool hasPIDENT = false;
    bool hasPPOS = false;

    auto it = begin(fields);
    auto itEnd = end(fields);
    for (typename BlastMatchField<g>::Enum const f : fieldList)
    {
        // this field represents multiple fields
        if (SEQAN_UNLIKELY(f == BlastMatchField<g>::Enum::STD))
        {
            for (typename BlastMatchField<g>::Enum const f2 :
                 BlastMatchField<g>::defaults)
            {
                if (SEQAN_UNLIKELY(it == itEnd))
                    throw std::ios_base::failure(BAD_FORMAT);

                _readField(match, *it, f2);
                ++it;
            }
            hasPIDENT = true;
        } else
        {
            if (SEQAN_UNLIKELY(it == itEnd))
                throw std::ios_base::failure(BAD_FORMAT);

            _readField(match, *it, f);
            ++it;
            if (f == BlastMatchField<g>::Enum::P_IDENT)
                hasPIDENT = true;
            if (f == BlastMatchField<g>::Enum::P_POS)
                hasPPOS = true;
        }
        skipLine(iter); // skip possibly remaining fields and goto next line
    }

    // retransform the percentages to real numbers and compute gaps
    if (_memberIsSet(match.aliLength) && (hasPIDENT))
        match.identities = ROUND((match.aliLength * match.identities) / 10000);
    if (_memberIsSet(match.aliLength) && (hasPPOS))
        match.positives = ROUND((match.aliLength * match.positives) / 10000);
     if (_memberIsSet(match.aliLength) && _memberIsSet(match.identities) &&
         _memberIsSet(match.mismatches) && !_memberIsSet(match.gaps))
        match.gaps = match.aliLength - match.mismatches - match.identities;
}

/*!
 * @fn BlastMatch#readMatch
 * @brief read a match from a file in BlastFormat
 *
 * @signature int readMatch(blastMatch, iterator/stream, [fieldList,] BlastFormat);
 *
 * @param[out]      blastMatch  A @link BlastMatch @endlink object to hold all relevant info
 * @param[in,out]   iter        An interator or stream
 * @param[in]       fieldList   A Sequence of @link BlastMatchField @endlink
 * @param[in]       tag         @link BlastFormat @endlink tag, only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported.
 *
 * This signature works with @link BlastMatch @endlinkes which is the recommended
 * way. The fieldList parameter can be specified if the file if you expect a
 * custom column composition (this only works for BlastFormatGeneration ==
 * BLAST_PLUS).
 *
 * Specifying the fieldList parameter implies that you know the columns are not
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
 *  * computation of the number of identities (from the percentage) [default]
 *  * computation of the number of positives (from the percentage) [if given]
 *  * number of gaps computed from other values [default]
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
 * @headerfile seqan/blast.h
 */

template <typename TqId,
          typename TsId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readMatch(BlastMatch<TqId, TsId, TAlign, TPos> & match,
          TFwdIterator & iter,
          TFieldList const & fieldList,
          BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "fieldList parameter may only be specified for "
                  "BlastFormatGeneration::BLAST_PLUS");

    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _readMatchImpl(match, iter, fieldList, TFormat());
}

template <typename TqId,
          typename TsId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readMatch(BlastMatch<TqId, TsId, TAlign, TPos> & match,
          TFwdIterator & iter,
          TFieldList const & fieldList,
          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "fieldList parameter may only be specified for "
                  "BlastFormatGeneration::BLAST_PLUS");

    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _readMatchImpl(match, iter, fieldList, TFormat());
}


// default fields
template <typename TqId,
          typename TsId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline void
readMatch(BlastMatch<TqId, TsId, TAlign, TPos>      & match,
          TFwdIterator & iter,
          BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    _readMatchImpl(match, iter, BlastMatchField<g>::defaults, TFormat());
}

// default fields
template <typename TqId,
          typename TsId,
          typename TFwdIterator,
          typename TPos,
          typename TAlign,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline void
readMatch(BlastMatch<TqId, TsId, TAlign, TPos>      & match,
          TFwdIterator & iter,
          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat; // same as above
    _readMatchImpl(match, iter, BlastMatchField<g>::defaults, TFormat());
}

/*!
 * @fn BlastFormat#readMatch
 * @brief read arbitrary columsn from a file in BlastFormat
 *
 * @signature int readMatch(iterator/stream, BlastFormat, args ...);
 *
 * @param[in,out]   iter        An interator or stream
 * @param[in]       tag         Only BlastFormatFile == TABULAR || TABULAR_WITH_HEADER supported
 * @param[out]      args        Arbitrary typed variables
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
 * @headerfile seqan/blast.h
 */

// arbitrary columns
template <typename TTarget,
          typename std::enable_if<IsSequence<TTarget>::VALUE>::type = 0>
inline void
_assignOrCast(TTarget & target, CharString const & source)
{
    target = source;
}

template <typename TTarget,
          typename std::enable_if<Is<NumberConcept<TTarget>>::VALUE>::type = 0>
inline void
_assignOrCast(TTarget & target, CharString const & source)
{
    target = lexicalCast<TTarget>(source);
}

template <typename TFwdIterator,
          typename TPos,
          typename TArg>
inline void
_readMatchImplBlastTab(TFwdIterator & iter, TArg & arg)
{
    CharString buf;
    try
    {
        readUntilChar(buf, iter, OrFunctor<IsTab,IsNewline>());
        _assignOrCast(arg, buf);
    } catch (UnexpectedEnd const &)
    {
        return;
    }
    // as this is the last requested field, go to beginning of next line
    skipLine(iter);
}

template <typename TFwdIterator,
          typename TPos,
          typename TArg,
          typename... TArgs>
inline void
_readMatchImplBlastTab(TFwdIterator & iter, TArg & arg, TArgs & ... args)
{
    CharString buf;
    readUntilChar(buf, iter, IsTab());
    goNext(iter); //skip '\t', go to begin of next field
    _assignOrCast(arg, buf);

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
        throw std::ios_base::failure(BAD_FORMAT);

    _readMatchImplBlastTab(iter, args...);
}

template <typename TFwdIterator,
          typename ... TArgs,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline void
readMatch(TFwdIterator & iter,
          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &,
          TArgs ... args)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat; // same as above
    readMatch(iter, TFormat(), args...);
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastRecord#readRecord
 * @brief read a record from a file in BlastFormat
 *
 * @signature int readRecord(blastRecord, iterator/stream, BlastFormat);
 *
 * @param[out]      blastRecord A BlastMatch object to hold all relevant info
 * @param[in,out]   iter      RecordReader
 * @param[in]       tag         BlastFormat specialization
 *
 *
 * @headerfile seqan/blast.h
 */

template <typename TFwdIterator,
          typename TDbName = CharString,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign>   & blastRecord,
           TDbName & dbName,
           TFwdIterator & iter,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    typedef BlastMatch<TQId, TSId, TAlign, TPos> TBlastMatch;

    CharString versionString; // -> /dev/null

    readHeader(blastRecord.qId,
               dbName,
               versionString,
               iter,
               false,
               TFormat());

    TBlastMatch match;
    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
    {
        readMatch(match, iter, TFormat());
        appendValue(blastRecord.matches, match);
    }
}

// custom fields
template <typename TFwdIterator,
          typename TDbName = CharString,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign>   & blastRecord,
           TDbName & dbName,
           TFwdIterator & iter,
           TFieldList const & fieldList,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "fieldList parameter may only be specified for "
                  "BlastFormatGeneration::BLAST_PLUS");

    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    typedef BlastMatch<TQId, TSId, TAlign, TPos> TBlastMatch;

    CharString versionString; // -> /dev/null

    readHeader(blastRecord.qId,
               dbName,
               versionString,
               iter,
               false,
               TFormat());

    TBlastMatch match;
    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
    {
        // only difference to above ↓
        readMatch(match, iter, fieldList, TFormat());
        appendValue(blastRecord.matches, match);
    }
}

// detect fields
template <typename TFwdIterator,
          typename TDbName = CharString,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign>   & blastRecord,
           TDbName & dbName,
           TFwdIterator & iter,
           DetectFields const &,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &)
{
    static_assert(g == BlastFormatGeneration::BLAST_PLUS,
                  "DetectFields parameter may only be specified for "
                  "BlastFormatGeneration::BLAST_PLUS");

    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    typedef BlastMatch<TQId, TSId, TAlign, TPos> TBlastMatch;
    typedef std::vector<BlastMatchField<BlastFormatGeneration::BLAST_PLUS>::Enum>
            TFieldList;

    CharString versionString; // -> /dev/null

    TFieldList fieldList;
    readHeader(blastRecord.qId,
               dbName,
               versionString,
               fieldList, // out-parameter to readHeader()
               iter,
               false,
               TFormat());

    TBlastMatch match;
    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
    {
        // in-parameter to readMatch ↓
        readMatch(match, iter, fieldList, TFormat());
        appendValue(blastRecord.matches, match);
    }
}

// TABULAR FORMAT had no headers, so we compare first query id
template <typename TFwdIterator,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign>   & blastRecord,
           TFwdIterator & iter,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g>  TFormat;
    typedef BlastMatch<TQId, TSId, TAlign, TPos> TBlastMatch;

    CharString curId, lastId;
    TBlastMatch match;

    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
    {
        // read current read ID and then rewind stream to beginning of line
        readUntilChar(curId, iter, OrFunctor<IsTab,IsNewline>());
        seek(iter, -length(curId)-1, SEEK_CURRENT);

        if ((curId != lastId) && (lastId != ""))
            return; // new Record reached

        readMatch(match, iter);

        appendValue(blastRecord.matches, match);
    }

    if (length(blastRecord.matches) == 0)
        throw std::ios_base::failure(BAD_FORMAT);

    blastRecord.qId = blastRecord.matches[0].qId;
}

// custom fields
template <typename TFwdIterator,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned,
          typename TFieldList,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readRecord(BlastRecord<TQId, TSId, TPos, TAlign>   & blastRecord,
           TFwdIterator & iter,
           TFieldList & fieldList,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const &)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g>  TFormat;
    typedef BlastMatch<TQId, TSId, TAlign, TPos> TBlastMatch;

    CharString curId, lastId;
    TBlastMatch match;

    //TODO for custom fields we have to change this, i.e. read entire line
    //and see if query id is contained
    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
    {
        // read current read ID and then rewind stream to beginning of line
        readUntilChar(curId, iter, OrFunctor<IsTab,IsNewline>());
        seek(iter, -length(curId)-1, SEEK_CURRENT);

        if ((curId != lastId) && (lastId != ""))
            return; // new Record reached

        readMatch(match, iter, fieldList);

        appendValue(blastRecord.matches, match);
    }

    if (length(blastRecord.matches) == 0)
        throw std::ios_base::failure(BAD_FORMAT);

    blastRecord.qId = blastRecord.matches[0].qId;
}

} // namespace seqan

#endif // SEQAN_BLAST_READ_BLAST_TABULAR_H_
