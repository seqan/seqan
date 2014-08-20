// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2014, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// This file contains routines to read BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_READ_BLAST_TABULAR_H_
#define SEQAN_EXTRAS_BLAST_READ_BLAST_TABULAR_H_

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
// ============================================================================v //

/*!
 * @fn BlastFormat#onMatch
 * @brief returns whether RecordReader is on beginning of match
 *
 * @signature bool onMatch(reader, tag)
 *
 * @param[in] reader    RecordReader
 * @param[in] tag       BlastFormat specialization
 *
 * @return     true or false
 */

template <typename TFile,
          typename TPass,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline bool
onMatch(RecordReader<TFile, TPass > & reader,
         BlastFormat<BlastFormatFile::TABULAR,
                     p,
                     g>       const & /*tag*/)
{
    if (value(reader) == '#')
        return false;
    return true;
}

template <typename TFile,
          typename TPass,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline bool
onMatch(RecordReader<TFile, TPass > & reader,
         BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                     p,
                     g>       const & /*tag*/)
{
    return onMatch(reader,
                    BlastFormat<BlastFormatFile::TABULAR,p,g>());
}


// verify whether the fields are default fields
template <typename TString,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
_verifyFields(StringSet<TString>    const & fields,
              unsigned              const   hits,
              // irrelevant for traditional header
              BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                          p,
                          g>        const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST>  TFormat;

    if (g == BlastFormatGeneration::BLAST_PLUS)
        if ((hits == 0) && length(fields) )
            return 0;


    CharString fieldStr;
    joinStringSet(fieldStr, fields, ", ");

    // it is only relevant that the firtst 12 fields by as expected
    return prefix(fieldStr, length(_defaultFields(TFormat())))
              == _defaultFields(TFormat());
}

// ----------------------------------------------------------------------------
// Function readHeader()                               [Single pass]
// ----------------------------------------------------------------------------

template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFile,
          typename TPass,
          typename TString,
          typename TString2,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
_readHeaderImplBlastTab(TqId                                    & qId,
                        TDBName                                 & dbName,
                        TVersionString                          & versionString,
                        StringSet<TString>                      & fields,
                        unsigned long                           & hits,
                        StringSet<TString2>                     & otherLines,
                        // any other lines
                        RecordReader<TFile, TPass >             & reader,
                        bool                              const   strict,
                        BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                    p,
                                    g>                    const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

    // this is a record instead of a header
    if (onMatch(reader, TFormat()))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    int ret = 0;
    bool keyCorrect = false;
    int queryLinePresent = 0;
    int dbLinePresent = 0;
    int fieldsLinePresent = 0;
    int hitsLinePresent = 0;
    int lastLine = 0;

    while ((!atEnd(reader)) && (!onMatch(reader, TFormat())))// in Header
    {
        // skip '#'
        ret = goNext(reader);
        if (ret)
            return ret;
        //skip blanks following the '#'
        ret = skipBlanks(reader);
        if (ret)
            return ret;

        CharString key = "";
        ret = readUntilWhitespace(key, reader);
        if (ret)
            return ret;

        if (key == _programTagToString(TFormat()))
        {
            // get whole line
            CharString buf;
            ret = readLine(buf, reader);
            if (ret)
                return ret;
            CharString &fullLine = key;
            append(fullLine, buf);

            versionString = fullLine;
            keyCorrect = true;
        }
        else if (key == "Query:")
        {
            ret = skipBlanks(reader);
            if (ret)
                return ret;
            ret = readLine(qId, reader);
            if (ret)
                return ret;
            ++queryLinePresent;
        }
        else if (key == "Database:")
        {
            ret = skipBlanks(reader);
            if (ret)
                return ret;
            ret = readLine(dbName, reader);
            if (ret)
                return ret;
            ++dbLinePresent;
        }
        else if (key == "Fields:")
        {
            ret = skipBlanks(reader);
            if (ret)
                return ret;

            CharString buf;
            ret = readLine(buf, reader);
            if (ret)
                return ret;
            strSplit(fields, buf, ", ");

            ++fieldsLinePresent;
            if (g == BlastFormatGeneration::BLAST)
                break; // header is finished
        }
        else
        {
            // get whole line
            CharString buf;
            ret = readLine(buf, reader);
            if (ret)
                return ret;
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

                    ret = !lexicalCast2(hits, hitsString);
                    if (ret && strict)
                        return RecordReader<TFile, TPass >::INVALID_FORMAT;

                    ++hitsLinePresent;
                    break; // header is finished
                }
            }
            else
                appendValue(otherLines, fullLine, Generous());
        }
    }

    if (!strict)
        return 0;

    if (g == BlastFormatGeneration::BLAST)
        if (  keyCorrect                &&
             (queryLinePresent   == 1)  &&
             (dbLinePresent      == 1)  &&
             (fieldsLinePresent  == 1)  &&
             (length(otherLines) == 0)   )
            return 0;

    if (g == BlastFormatGeneration::BLAST_PLUS)
        if (( keyCorrect                                   &&
             (queryLinePresent    == 1)                    &&
             (dbLinePresent       == 1)                    &&
             (hitsLinePresent     == 1)                    &&
             ( (fieldsLinePresent == 1) || (hits==0) )     &&
             (length(otherLines)  == 0)                     ) ||
             ( (lastLine          == 1) && atEnd(reader))   )
            return 0;

    return RecordReader<TFile, TPass >::INVALID_FORMAT;
}

/*!
 * @fn BlastFormat#readHeader
 * @brief read a Header from a Blast output file
 *
 * @signature int readHeader(qId, dbName, versionString, [hits,] [fields, [otherLines,]] recordReader, strict, tag);
 *
 * @param[out]  qId     String to hold the query ID from the header
 * @param[out]  dbName  String to hold the database name from the header
 * @param[out]  versionString  String to hold the Blast program Tag and version
 * @param[out]  hits    Numerical to hold the number of hits that will follow the header (only available in BlastPlus spec of BlastFormat)
 * @param[out]  fields  StringSet to hold column identifiers, useful if non-defaults are expected
 * @param[out]  otherLines  StringSet to hold any comment or header lines that are not identified otherwise
 * @param[in,out]   reader  RecordReader
 * @param[in]   strict  bool to signify whether the function should return error on a non-conforming header or just "get whatever possible". If not using strict, it is recommended to pass fields and otherLines and verify these manually.
 * @param[in]   tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
 *
 * @return     0 on success, and non-zero otherwise
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
          typename TFile,
          typename TPass,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
readHeader(TqId                                                 & qId,
           TDBName                                              & dbName,
           TVersionString                                       & versionString,
           RecordReader<TFile, TPass >                          & reader,
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

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    if (strict)
    {
        ret =  _verifyFields(fields, hits, TFormat());
        if (ret)
            return RecordReader<TFile, TPass >::INVALID_FORMAT;
    }

    return 0;
}

// BlastPlus with hit count
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFile,
          typename TPass,
          BlastFormatProgram p>
inline int
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           unsigned long                                    & hits,
           RecordReader<TFile, TPass >                      & reader,
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

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    if (strict)
    {
        ret =  _verifyFields(fields, hits, TFormat());
        if (ret)
            RecordReader<TFile, TPass >::INVALID_FORMAT;
    }

    return 0;
}


// with fields
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldString,
          typename TFile,
          typename TPass,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           StringSet<TFieldString>                          & fields,
           RecordReader<TFile, TPass >                      & reader,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g>                             const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

    StringSet<CharString> otherLines;
    unsigned long hits =0;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
    return 0;
}

// with fields and hits count
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldString,
          typename TFile,
          typename TPass,
          BlastFormatProgram p>
inline int
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           unsigned long                                    & hits,
           StringSet<TFieldString>                          & fields,
           RecordReader<TFile, TPass >                      & reader,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    StringSet<CharString> otherLines;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
    return 0;
}


// with fields and otherLines
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldString,
          typename TOtherString,
          typename TFile,
          typename TPass,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           StringSet<TFieldString>                          & fields,
           StringSet<TOtherString>                          & otherLines,
           RecordReader<TFile, TPass >                      & reader,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g>                             const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> TFormat;

    unsigned long hits =0;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
    return 0;
}

// with fields and otherLines and hits count
template <typename TqId,
          typename TDBName,
          typename TVersionString,
          typename TFieldString,
          typename TOtherString,
          typename TFile,
          typename TPass,
          BlastFormatProgram p>
inline int
readHeader(TqId                                             & qId,
           TDBName                                          & dbName,
           TVersionString                                   & versionString,
           unsigned long                                    & hits,
           StringSet<TFieldString>                          & fields,
           StringSet<TOtherString>                          & otherLines,
           RecordReader<TFile, TPass >                      & reader,
           bool                                       const   strict,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       BlastFormatGeneration::BLAST_PLUS> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST_PLUS> TFormat;

    int ret =  _readHeaderImplBlastTab(qId,
                                       dbName,
                                       versionString,
                                       fields,
                                       hits,
                                       otherLines,
                                       reader,
                                       strict,
                                       TFormat());
    if (ret)
        return ret;

    // don't verify fields, if user specified that he wants the list
    // of fields, because that implies that he expects non-defaults
    return 0;
}

// ----------------------------------------------------------------------------
// Function skipHeader()                               [Single pass]
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#skipHeader
 * @brief skip a header from a Blast tabular output file, optionally verifying it for format compliance.
 *
 * @signature int skipHeader(recordReader, [strict,] tag);
 *
 * @param[in,out]   reader  RecordReader
 * @param[in]   strict  bool to signify whether the function should return error on a non-conforming header or just "skip whatever possible".
 * @param[in]   tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
 *
 * @see BlastFormat#skipUntilMatch
 * @return     0 on success, and non-zero otherwise
 * @headerfile seqan/blast.h
 * @section Remarks
 *
 * call this function whenever you want to skip exactly one header. If you want
 * to go directly to the beginning of the next match (possibly skipping multiple
 * headers that have no succeeding matches) use skipUntilMatch instead.
 */

template <typename TFile,
          typename TPass,
          BlastFormatOptions::M m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
skipHeader(RecordReader<TFile, TPass >        & reader,
           bool                         const   strict,
           BlastFormat<m,p,g>           const & /*tag*/)
{
    CharString qId;
    CharString dbName;
    CharString versionString;
    StringSet<CharString> fields;
    unsigned long hits = 0;
    StringSet<CharString> otherLines;

    return  _readHeaderImplBlastTab(qId,
                                    dbName,
                                    versionString,
                                    fields,
                                    hits,
                                    otherLines,
                                    reader,
                                    strict,
                                    BlastFormat<m,p,g>());
}

template <typename TFile,
          typename TPass,
          BlastFormatOptions::M m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
skipHeader(RecordReader<TFile, TPass > & reader,
           BlastFormat<m,p,g> const & /*tag*/)
{
    return skipHeader(reader, false, BlastFormat<m,p,g>());
}

// ----------------------------------------------------------------------------
// Function skipUntilMatch()                           [Single pass]
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#skipUntilMatch
 * @brief skip arbitrary number of headers and/or comment lines until the beginning of a match is reached
 *
 * @signature int skipUntilMatch(recordReader, tag);
 *
 * @param[in,out]   reader  RecordReader
 * @param[in]       tag     BlastFormat specialization,
 * with BlastFormat::_m == BlastFormatFile::TABULAR || BlastFormatFile::TABULAR_WITH_HEADER
 *
 * @return     0 on success, and non-zero otherwise
 * @see BlastFormat#skipHeader
 * @headerfile seqan/blast.h
 * @section Remarks
 *
 * call this function whenever you are on a comment character ('#') in the file and want to jump to the beginning of the next record. If you want skip only a single header (to count skipped headers or verify its conformance
to standards), use BlastFormat#skipHeader instead.
 */

template <typename TFile,
          typename TPass,
          BlastFormatOptions::M m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
skipUntilMatch(RecordReader<TFile, TPass >    & reader,
                BlastFormat<m,p,g>      const & /*tag*/)
{
    int ret = 0;
    while ((!atEnd(reader)) && value(reader) == '#') // skip comments
    {
        ret = skipLine(reader);
        if (ret)
            return ret;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function readMatch()                               [Single pass]
// ----------------------------------------------------------------------------

template <typename TString,
          typename TFile,
          typename TPass>
inline int
_readMatchImplBlastTabArbitrary(StringSet<TString>          & fields,
                                RecordReader<TFile, TPass > & reader)
{
    CharString buf;
    int ret = readLine(buf, reader);
    if ( (ret) && (!atEnd(reader)) ) // file could end w/o empty newline
        return ret;

    strSplit(fields, buf, '\t', false);
    return 0;
}

template <typename TqId,
          typename TsId,
          typename TFile,
          typename TPass,
          typename TPos>
inline int
_readMatchImplBlastTabDefault(TqId      & qId,
                              TsId      & sId,
                              unsigned  & identities,
                              unsigned  & ali_length,
                              unsigned  & num_mismatches,
                              unsigned  & gap_openings,
                              TPos      & qStart,
                              TPos      & qEnd,
                              TPos      & sStart,
                              TPos      & sEnd,
                              double    & eval,
                              double    & bitScore,
                              RecordReader<TFile, TPass > & reader)
{
    int ret = 0;
    double percentIdent = 0;

    for (int i = 0; i < 12; ++i)
    {
        CharString buf;
        if (i < 11)
        {
            ret = readUntilChar(buf, reader, '\t');
            if (ret)
                return ret;
            ret = goNext(reader); //skip '\t', go to begin of next field
            if (ret)
                return ret;
        } else
        {
            ret = readUntilTabOrLineBreak(buf, reader);
            if (ret)
                return ret;
            ret = skipLine(reader); // skip extra fields or '\n'
            if ( (ret) && (!atEnd(reader)) ) // file may end w/o empty newline
                return ret;
        }

        switch(i)
        {
            case  0: ret = lexicalCast2(qId,            buf); break;
            case  1: ret = lexicalCast2(sId,            buf); break;
            case  2: ret = lexicalCast2(percentIdent,   buf); break;
            case  3: ret = lexicalCast2(ali_length,     buf); break;
            case  4: ret = lexicalCast2(num_mismatches, buf); break;
            case  5: ret = lexicalCast2(gap_openings,   buf); break;
            case  6: ret = lexicalCast2(qStart,         buf); break;
            case  7: ret = lexicalCast2(qEnd,           buf); break;
            case  8: ret = lexicalCast2(sStart,         buf); break;
            case  9: ret = lexicalCast2(sEnd,           buf); break;
            case 10: ret = lexicalCast2(eval,           buf); break;
            case 11: ret = lexicalCast2(bitScore,       buf); break;
        }
        ret = !ret; // lexicalCast2 returns true when successful
        if (ret)
            return RecordReader<TFile, TPass >::INVALID_FORMAT;
    }

    identities = std::ceil(percentIdent * ali_length / 100);

    return 0;
}

/*!
 * @fn BlastFormat#readMatch
 * @brief read a match from a file in BlastFormat
 *
 * @signature int readMatch(blastMatch, recordReader, BlastFormat);
 *
 * @param[out]      blastMatch A BlastMatch object to hold all relevant info
 * @param[in,out]   reader  RecordReader
 * @param[in]       tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
 *
 * @signature int readMatch(qId, sId, identities, aliLength, numMismatches, gapOpenings, qStart, qEnd, sStart, sEnd, eValue, bitScore, recordReader, BlastFormat);
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
 * @param[in,out]    reader  RecordReader
 * @param[in]    tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
 *
 * @signature int readMatch(fields, recordReader, BlastFormat);
 *
 * @param[out]      fields a StringSet with all the columns as entries
 * @param[in,out]   reader  RecordReader
 * @param[in]       tag     BlastFormat specialization, with BlastFormat::_m == TABULAR || TABULAR_WITH_HEADER
 *
 *
 * @return     0 on success, and non-zero otherwise
 * @headerfile seqan/blast.h
 */

template <typename TqId,
          typename TsId,
          typename TFile,
          typename TPass,
          typename TPos,
          BlastFormatOptions::M          m,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline int
_readMatch(TqId                         & qId,
           TsId                         & sId,
           TPos                         & identities,
           TPos                         & ali_length,
           TPos                         & num_mismatches,
           TPos                         & gap_openings,
           TPos                         & qStart,
           TPos                         & qEnd,
           TPos                         & sStart,
           TPos                         & sEnd,
           double                       & eval,
           double                       & bitScore,
           RecordReader<TFile, TPass >  & reader,
           BlastFormat<m,p,g>     const & /*tag*/)
{
    // header should have been read or skipped
    if (!onMatch(reader, BlastFormat<m,p,g>()))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    return _readMatchImplBlastTabDefault(qId,
                                          sId,
                                          identities,
                                          ali_length,
                                          num_mismatches,
                                          gap_openings,
                                          qStart,
                                          qEnd,
                                          sStart,
                                          sEnd,
                                          eval,
                                          bitScore,
                                          reader);
}

template <typename TqId,
          typename TsId,
          typename TFile,
          typename TPass,
          typename TPos,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline int
readMatch(TqId                         & qId,
          TsId                         & sId,
          TPos                         & identities,
          TPos                         & ali_length,
          TPos                         & num_mismatches,
          TPos                         & gap_openings,
          TPos                         & qStart,
          TPos                         & qEnd,
          TPos                         & sStart,
          TPos                         & sEnd,
          double                       & eval,
          double                       & bitScore,
          RecordReader<TFile, TPass >  & reader,
          BlastFormat<BlastFormatFile::TABULAR,p,g>     const & /*tag*/)
{
    return _readMatch(qId, sId, identities, ali_length, num_mismatches,
                      gap_openings, qStart, qEnd, sStart, sEnd, eval, bitScore,
                      reader, BlastFormat<BlastFormatFile::TABULAR,p,g>());
}

template <typename TqId,
          typename TsId,
          typename TFile,
          typename TPass,
          typename TPos,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline int
readMatch(TqId                         & qId,
          TsId                         & sId,
          TPos                         & identities,
          TPos                         & ali_length,
          TPos                         & num_mismatches,
          TPos                         & gap_openings,
          TPos                         & qStart,
          TPos                         & qEnd,
          TPos                         & sStart,
          TPos                         & sEnd,
          double                       & eval,
          double                       & bitScore,
          RecordReader<TFile, TPass >  & reader,
          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,p,g> const &/*tag*/)
{
    return _readMatch(qId, sId, identities, ali_length, num_mismatches,
                      gap_openings, qStart, qEnd, sStart, sEnd, eval, bitScore,
                      reader,
                      BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,p,g>());
}

template <typename TqId,
          typename TsId,
          typename TFile,
          typename TPass,
          typename TPos,
          typename TAlign,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline int
readMatch(BlastMatch<TqId, TsId, TAlign, TPos>      & match,
           RecordReader<TFile, TPass >              & reader,
           BlastFormat<BlastFormatFile::TABULAR,
                       p,
                       g>                     const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR,p,g> TFormat;
    // header should have been read or skipped
    if (!onMatch(reader, TFormat()))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    return _readMatchImplBlastTabDefault(match.qId,
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
                                         match.bitScore,
                                         reader);
}

template <typename TqId,
          typename TsId,
          typename TFile,
          typename TPass,
          typename TPos,
          typename TAlign,
          BlastFormatProgram    p,
          BlastFormatGeneration g>
inline int
readMatch(BlastMatch<TqId, TsId, TAlign, TPos>      & match,
           RecordReader<TFile, TPass >              & reader,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g>                     const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,p,g> TFormat;
    // header should have been read or skipped
    if (!onMatch(reader, TFormat()))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    return _readMatchImplBlastTabDefault(match.qId,
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
                                         match.bitScore,
                                         reader);
}

template <typename TFile,
          typename TPass,
          typename TString,
          BlastFormatOptions::M m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
readMatch(StringSet<TString>            & fields,
           RecordReader<TFile, TPass >  & reader,
           BlastFormat<m,p,g>     const & /*tag*/)
{
    // header should have been read or skipped
    if (!onMatch(reader, BlastFormat<m,p,g>()))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    return _readMatchImplBlastTabArbitrary(fields, reader);
}

// ----------------------------------------------------------------------------
// Function readRecord()                               [Single pass]
// ----------------------------------------------------------------------------

/*!
 * @fn BlastFormat#readRecord readRecord
 * @brief read a record from a file in BlastFormat
 *
 * @signature int readRecord(blastRecord, recordReader, BlastFormat);
 *
 * @param[out]      blastRecord A BlastMatch object to hold all relevant info
 * @param[in,out]   reader      RecordReader
 * @param[in]       tag         BlastFormat specialization
 *
 * @return     0 on success, and non-zero otherwise
 *
 * @headerfile seqan/blast.h
 */


template <typename TFile,
          typename TPass,
          typename TDbName = CharString,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
readRecord(BlastRecord<TDbName, TQId, TSId, TAlign, TPos>   & blastRecord,
           RecordReader<TFile, TPass >                      & reader,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g>                             const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    typedef BlastMatch<TQId, TSId, TAlign, TPos> TBlastMatch;

    CharString versionString; // -> /dev/null

    int ret = readHeader(blastRecord.qId,
                         blastRecord.dbName,
                         versionString,
                         reader,
                         false,
                         TFormat());
    if (ret)
        return ret;
    while ((!atEnd(reader)) && onMatch(reader, TFormat()))
    {
        TBlastMatch match;
        ret = _readMatchImplBlastTabDefault(match.qId,
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
                                            reader);
        if (ret)
            return ret;
        appendValue(blastRecord.matches, match);
    }
    return 0;
}

template <typename TFile,
          typename TPass,
          typename TDbName = CharString,
          typename TQId = CharString,
          typename TSId = CharString,
          typename TAlign = Align<CharString, ArrayGaps>,
          typename TPos = unsigned,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
readRecord(BlastRecord<TDbName, TQId, TSId, TAlign, TPos>   & blastRecord,
           RecordReader<TFile, TPass >                      & reader,
           BlastFormat<BlastFormatFile::TABULAR,
                       p,
                       g>                             const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g>  TFormat;
    typedef BlastMatch<TQId, TSId, TAlign, TPos>            TBlastMatch;

    CharString curId, lastId;

    while ((!atEnd(reader)) && onMatch(reader, TFormat()))
    {
        // read current read ID and then rewind stream to beginning of line
        int ret = readUntilTabOrLineBreak(curId, reader);
        if (ret)
            return ret;
        ret = streamSeek(reader, -length(curId)-1, SEEK_CUR);
        if (ret)
            return ret;

        if ((curId != lastId) && (lastId != ""))
            return 0; // new Record reached

        TBlastMatch match;
        ret = _readMatchImplBlastTabDefault(match.qId,
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
                                            reader);
        if (ret)
            return ret;
        appendValue(blastRecord.matches, match);
    }
    if (length(blastRecord.matches) == 0)
        return RecordReader<TFile, TPass>::INVALID_FORMAT;

    blastRecord.qId = blastRecord.matches[0].qId;

    return 0;
}



} // namespace seqan
#endif // header guard
