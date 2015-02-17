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

#ifndef SEQAN_BLAST_READ_BLAST_TABULAR_H_
#define SEQAN_BLAST_READ_BLAST_TABULAR_H_

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

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

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
 * @    true or false
 */

template <typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline bool
onMatch(TFwdIterator & iter,
        BlastFormat<BlastFormatFile::TABULAR,
                    p,
                    g>       const & /*tag*/)
{
    (value(iter) != '#');
}

template <typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline bool
onMatch(TFwdIterator & iter,
        BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                    p,
                    g>       const & /*tag*/)
{
    onMatch(iter,
                   BlastFormat<BlastFormatFile::TABULAR,p,g>());
}


// verify whether the fields are default fields
template <typename TString,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
_verifyFields(StringSet<TString>    const & fields,
              unsigned              const   hits,
              // irrelevant for traditional header
              BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                          p,
                          g>        const & /*tag*/)
{
//     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
//                         p,
//                         BlastFormatGeneration::BLAST>  TFormat;

    if (g == BlastFormatGeneration::BLAST_PLUS)
        if ((hits == 0) && length(fields) )
            0;


    CharString fieldStr;
    joinStringSet(fieldStr, fields, ", ");

    // it is only relevant that the firtst 12 fields by as expected
    prefix(fieldStr, length(BlastMatchField<g>::defaults))
              == BlastMatchField<g>::defaults;
}

// ----------------------------------------------------------------------------
// Function readHeader()                               [Single pass]
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
                        BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                    p,
                                    g>                    const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

    // this is a record instead of a header
    if (onMatch(iter, TFormat()))
        throw std::ios_base::failure(BAD_FORMAT);

    int 0;
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
                    ++lastLine;
                // is hits counter?
                else if (hasSuffix(fullLine, "hits found"))
                {
                    CharString hitsString = "";
                    for (unsigned i = 0;
                        (i < length(fullLine) && isdigit(fullLine[i]));
                        ++i)
                        append(hitsString, fullLine[i], Generous());

                    !lexicalCast2(hits, hitsString);
                    if (ret && strict)
                        throw std::ios_base::failure(BAD_FORMAT);

                    ++hitsLinePresent;
                    break; // header is finished
                }
            }
            else
                appendValue(otherLines, fullLine, Generous());
        }
    }

    if (!strict)
        0;

    if (g == BlastFormatGeneration::BLAST_LEGACY)
        if (  keyCorrect                &&
             (queryLinePresent   == 1)  &&
             (dbLinePresent      == 1)  &&
             (fieldsLinePresent  == 1)  &&
             (length(otherLines) == 0)   )
            0;

    if (g == BlastFormatGeneration::BLAST_PLUS)
        if (( keyCorrect                                   &&
             (queryLinePresent    == 1)                    &&
             (dbLinePresent       == 1)                    &&
             (hitsLinePresent     == 1)                    &&
             ( (fieldsLinePresent == 1) || (hits==0) )     &&
             (length(otherLines)  == 0)                     ) ||
             ( (lastLine          == 1) && atEnd(iter))   )
            0;

    throw std::ios_base::failure(BAD_FORMAT);
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
 * @param[in]   tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
 *
 * @    0 on success, and non-zero otherwise
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
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g>                                 const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

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
          typename TFieldString,
          typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           StringSet<TFieldString>                          & fields,
           TFwdIterator & iter,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g>                             const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

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
          typename TFieldString,
          typename TFwdIterator,
          BlastFormatProgram p>
inline void
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           unsigned long                                    & hits,
           StringSet<TFieldString>                          & fields,
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
          typename TFieldString,
          typename TOtherString,
          typename TFwdIterator,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           StringSet<TFieldString>                          & fields,
           StringSet<TOtherString>                          & otherLines,
           TFwdIterator & iter,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g>                             const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

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
          typename TFieldString,
          typename TOtherString,
          typename TFwdIterator,
          BlastFormatProgram p>
inline void
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           unsigned long                                    & hits,
           StringSet<TFieldString>                          & fields,
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
// Function skipHeader()                               [Single pass]
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#skipHeader
 * @brief skip a header from a Blast tabular output file, optionally verifying it for format compliance.
 *
 * @signature int skipHeader(iterator/stream, [strict,] tag);
 *
 * @param[in,out]   iter  RecordReader
 * @param[in]   strict  bool to signify whether the function should error on a non-conforming header or just "skip whatever possible".
 * @param[in]   tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
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
           BlastFormat<m,p,g>           const & /*tag*/)
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
                            BlastFormat<m,p,g>());
}

template <typename TFwdIterator,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
skipHeader(TFwdIterator & iter,
           BlastFormat<m,p,g> const & /*tag*/)
{
    skipHeader(iter, false, BlastFormat<m,p,g>());
}

// ----------------------------------------------------------------------------
// Function skipUntilMatch()                           [Single pass]
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#skipUntilMatch
 * @brief skip arbitrary number of headers and/or comment lines until the beginning of a match is reached
 *
 * @signature int skipUntilMatch(iterator/stream, tag);
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
                BlastFormat<m,p,g>      const & /*tag*/)
{
    while ((!atEnd(iter)) && value(iter) == '#') // skip comments
        skipLine(iter);
}

// ----------------------------------------------------------------------------
// Function readMatch()                               [Single pass]
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#readMatch
 * @brief read a match from a file in BlastFormat
 *
 * @signature int readMatch(blastMatch, iterator/stream, BlastFormat);
 *
 * @param[out]      blastMatch A BlastMatch object to hold all relevant info
 * @param[in,out]   iter  RecordReader
 * @param[in]       tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
 *
 * @signature int readMatch(qId, sId, identities, aliLength, numMismatches, gapOpenings, qStart, qEnd, sStart, sEnd, eValue, bitScore, iterator/stream, BlastFormat);
 * @param[out]   qId ID-String of the query
 * @param[out]   sId ID-String of the subject (sequence in database)
 * @param[out]   identities number of identies in alignment (not percentage!)
 * @param[out]   aliLength length of alignment
 * @param[out]   numMismatches number of mismatches in alignment
 * @param[out]   gapOpenings number consecutive gaps per alignment
 * @param[out]   qStart alignment begin position on query
 * @param[out]   qEnd alignment end position on query
 * @param[out]   sStart alignment begin position on subject
 * @param[out]   sEnd alignment end position on subject
 * @param[out]   eValue alignment e-Value
 * @param[out]   bitScore alignment bit-Score
 * @param[in,out]    iter  RecordReader
 * @param[in]    tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
 *
 * @signature int readMatch(fields, iterator/stream, BlastFormat);
 *
 * @param[out]      fields a StringSet with all the columns as entries
 * @param[in,out]   iter  RecordReader
 * @param[in]       tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
 *
 *
 * @    0 on success, and non-zero otherwise
 * @headerfile seqan/blast.h
 */

// template <typename TString,
//           typename TFile,
//           typename TPass>
// inline void
// _readMatchImplBlastTabArbitrary(StringSet<TString>          & fields,
//                                 TFwdIterator & iter)
// {
//     CharString buf;
//     readLine(buf, iter);
//     strSplit(fields, buf, '\t', false);
// }
//
// template <typename TqId,
//           typename TsId,
//           typename TFwdIterator,
//           typename TPos>
// inline void
// _readMatchImplBlastTabDefault(TqId      & qId,
//                               TsId      & sId,
//                               unsigned  & identities,
//                               unsigned  & ali_length,
//                               unsigned  & num_mismatches,
//                               unsigned  & gap_openings,
//                               TPos      & qStart,
//                               TPos      & qEnd,
//                               TPos      & sStart,
//                               TPos      & sEnd,
//                               double    & eval,
//                               double    & bitScore,
//                               TFwdIterator & iter)
// {
//     int ret = 0;
//     double percentIdent = 0;
//
//     for (int i = 0; i < 12; ++i)
//     {
//         CharString buf;
//         if (i < 11)
//         {
//             readUntilChar(buf, iter, '\t');
//             goNext(iter); //skip '\t', go to begin of next field
//         } else
//         {
//             readUntilTabOrLineBreak(buf, iter);
//             skipLine(iter); // skip extra fields or '\n'
//         }
//
//         switch(i)
//         {
//             case  0: ret = lexicalCast2(qId,            buf); break;
//             case  1: ret = lexicalCast2(sId,            buf); break;
//             case  2: ret = lexicalCast2(percentIdent,   buf); break;
//             case  3: ret = lexicalCast2(ali_length,     buf); break;
//             case  4: ret = lexicalCast2(num_mismatches, buf); break;
//             case  5: ret = lexicalCast2(gap_openings,   buf); break;
//             case  6: ret = lexicalCast2(qStart,         buf); break;
//             case  7: ret = lexicalCast2(qEnd,           buf); break;
//             case  8: ret = lexicalCast2(sStart,         buf); break;
//             case  9: ret = lexicalCast2(sEnd,           buf); break;
//             case 10: ret = lexicalCast2(eval,           buf); break;
//             case 11: ret = lexicalCast2(bitScore,       buf); break;
//         }
//         ret = !ret; // lexicalCast2 returns true when successful
//         if (ret)
//             throw std::ios_base::failure(BAD_FORMAT);
//     }
//
//     identities = std::ceil(percentIdent * ali_length / 100);
// }

template <typename TTarget,
          typename std::enable_if<IsSequence<TTarget>>::VALUE> = 0>
inline void
_assignOrCast(TTarget & target, CharString const & source)
{
    target = source;
}

template <typename TTarget,
          typename std::enable_if<Is<NumberConcept<TTarget>>::VALUE> = 0>
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
          typename ... TArgs>
inline void
_readMatchImplBlastTab(TFwdIterator & iter, TArg & arg, TArgs & args ...)
{
    CharString buf;
    readUntilChar(buf, iter, IsTab());
    goNext(iter); //skip '\t', go to begin of next field
    _assignOrCast(arg, buf);

    // recurse to next argument
    _readMatchImplBlastTab(iter, args);
}

// template <typename TqId,
//           typename TsId,
//           typename TFwdIterator,
//           typename TPos,
//           BlastFormatProgram    p,
//           BlastFormatGeneration g>
// inline void
// readMatch(TqId                         & qId,
//           TsId                         & sId,
//           TPos                         & identities,
//           TPos                         & ali_length,
//           TPos                         & num_mismatches,
//           TPos                         & gap_openings,
//           TPos                         & qStart,
//           TPos                         & qEnd,
//           TPos                         & sStart,
//           TPos                         & sEnd,
//           double                       & eval,
//           double                       & bitScore,
//           TFwdIterator & iter,
//           BlastFormat<BlastFormatFile::TABULAR,p,g>     const & /*tag*/)
// {
//     _readMatchImplBlastTab(iter, qId, sId, identities, ali_length, num_mismatches,
//                gap_openings, qStart, qEnd, sStart, sEnd, eval, bitScore);
// }
//
// template <typename TqId,
//           typename TsId,
//           typename TFwdIterator,
//           typename TPos,
//           BlastFormatProgram    p,
//           BlastFormatGeneration g>
// inline void
// readMatch(TqId                         & qId,
//           TsId                         & sId,
//           TPos                         & identities,
//           TPos                         & ali_length,
//           TPos                         & num_mismatches,
//           TPos                         & gap_openings,
//           TPos                         & qStart,
//           TPos                         & qEnd,
//           TPos                         & sStart,
//           TPos                         & sEnd,
//           double                       & eval,
//           double                       & bitScore,
//           TFwdIterator & iter,
//           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,p,g> const &/*tag*/)
// {
//     _readMatchImplBlastTab(iter, qId, sId, identities, ali_length, num_mismatches,
//                gap_openings, qStart, qEnd, sStart, sEnd, eval, bitScore);
// }

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
          BlastFormat<BlastFormatFile::TABULAR,
                       p,
                       g>                     const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,p,g> TFormat;
    // header should have been read or skipped
    if (!onMatch(iter, TFormat()))
        throw std::ios_base::failure(BAD_FORMAT);

    // TODO load identities into buffer

    _readMatchImplBlastTab(iter,
                                  match.qId,
                                  match.sId,
                                  match.identities,
                                  match.aliLength,
                                  match.mismatches,
                                  match.gapOpenings,
                                  match.qStart,
                                  match.qEnd,
                                  match.sStart,
                                  match.sEnd,
                                  match.eval,
                                  match.bitScore);
    // TODO calculate some more things from given knowledge
    // TODO retransform coordinates
}

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
          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g>                     const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,p,g> TFormat; // same as above
    readMatch(match, iter, TFormat());
}

// custom arguments
template <typename TFwdIterator,
          typename ... Targs,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline void
readMatch(TFwdIterator & iter,
          BlastFormat<BlastFormatFile::TABULAR, p, g> const &,
          TArgs & args ...)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,p,g> TFormat;
    // header should have been read or skipped
    if (!onMatch(iter, TFormat()))
        throw std::ios_base::failure(BAD_FORMAT);

    _readMatchImplBlastTab(match, iter, args);
}

template <typename TFwdIterator,
          typename ... Targs,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline void
readMatch(TFwdIterator & iter,
          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &,
          TArgs ... args)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,p,g> TFormat; // same as above
    readMatch(match, iter, TFormat());
}

// template <typename TFwdIterator,
//           typename TString,
//           BlastFormatFile m,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline void
// readMatch(StringSet<TString>            & fields,
//            TFwdIterator & iter,
//            BlastFormat<m,p,g>     const & /*tag*/)
// {
//     // header should have been read or skipped
//     if (!onMatch(iter, BlastFormat<m,p,g>()))
//         throw std::ios_base::failure(BAD_FORMAT);
//
//     _readMatchImplBlastTabArbitrary(fields, iter);
// }

// ----------------------------------------------------------------------------
// Function readRecord()                               [Single pass]
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#readRecord readRecord
 * @brief read a record from a file in BlastFormat
 *
 * @signature int readRecord(blastRecord, iterator/stream, BlastFormat);
 *
 * @param[out]      blastRecord A BlastMatch object to hold all relevant info
 * @param[in,out]   iter      RecordReader
 * @param[in]       tag         BlastFormat specialization
 *
 * @    0 on success, and non-zero otherwise
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
readRecord(BlastRecord<TDbName, TQId, TSId, TAlign, TPos>   & blastRecord,
           TFwdIterator & iter,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g>                             const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    typedef BlastMatch<TQId, TSId, TAlign, TPos> TBlastMatch;

    CharString versionString; // -> /dev/null

    int readHeader(blastRecord.qId,
                         blastRecord.dbName,
                         versionString,
                         iter,
                         false,
                         TFormat());
    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
    {
        TBlastMatch match;
        _readMatchImplBlastTabDefault(match.qId,
                                            match.sId,
                                            match.identities,
                                            match.aliLength,
                                            match.mismatches,
                                            match.gapOpenings,
                                            match.qStart,
                                            match.qEnd,
                                            match.sStart,
                                            match.sEnd,
                                            match.eVal,
                                            match.bitScore,
                                            iter);
        appendValue(blastRecord.matches, match);
    }
    0;
}

template <typename TFwdIterator,
          typename TDbName = CharString,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline void
readRecord(BlastRecord<TDbName, TQId, TSId, TAlign, TPos>   & blastRecord,
           TFwdIterator & iter,
           BlastFormat<BlastFormatFile::TABULAR,
                       p,
                       g>                             const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g>  TFormat;
    typedef BlastMatch<TQId, TSId, TAlign, TPos>            TBlastMatch;

    CharString curId, lastId;

    while ((!atEnd(iter)) && onMatch(iter, TFormat()))
    {
        // read current read ID and then rewind stream to beginning of line
        int readUntilTabOrLineBreak(curId, iter);
        streamSeek(iter, -length(curId)-1, SEEK_CUR);

        if ((curId != lastId) && (lastId != ""))
            0; // new Record reached

        TBlastMatch match;
        _readMatchImplBlastTabDefault(match.qId,
                                            match.sId,
                                            match.identities,
                                            match.aliLength,
                                            match.mismatches,
                                            match.gapOpenings,
                                            match.qStart,
                                            match.qEnd,
                                            match.sStart,
                                            match.sEnd,
                                            match.eVal,
                                            match.bitScore,
                                            iter);
        appendValue(blastRecord.matches, match);
    }
    if (length(blastRecord.matches) == 0)
        RecordReader<TFile, TPass>::INVALID_FORMAT;

    blastRecord.qId = blastRecord.matches[0].qId;

    0;
}

} // namespace seqan

#endif // SEQAN_BLAST_READ_BLAST_TABULAR_H_
