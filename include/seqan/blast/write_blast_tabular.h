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
// This file contains routines to generate BLAST tab-seperated output
// ==========================================================================

#ifndef SEQAN_EXTRAS_BLAST_WRITE_BLAST_TABULAR_H_
#define SEQAN_EXTRAS_BLAST_WRITE_BLAST_TABULAR_H_

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

template <BlastFormatProgram p, BlastFormatGeneration g>
constexpr
const char *
_seperatorString(BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                              p,
                              g> const & /*tag*/)
{
    return ", ";
}

template <BlastFormatProgram p, BlastFormatGeneration g>
constexpr
const char *
_seperatorString(BlastFormat<BlastFormatFile::TABULAR,
                              p,
                              g> const & /*tag*/)
{
    return "\t";
}

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Helper functions for printing n columns (requires C++11 variadic templates)
// ----------------------------------------------------------------------------

template <typename TStream, typename TTag>
inline int
_writeFields(TStream & /**/, TTag const & /*tag*/)
{
    return 0;
}

template <typename TStream, typename TField, typename... TFields, typename TTag>
inline int
_writeFields(TStream & stream,
             TTag const & /*tag*/,
             TField const & field1, TFields const & ... fields)
{
    int ret = streamPut(stream, _seperatorString(TTag()));
    if (ret)
        return ret;

    ret = streamPut(stream, field1);
    if (ret)
        return ret;

    return _writeFields(stream, TTag(), fields... );
}

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
// Function writeHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastRecord#writeHeader
 * @headerfile seqan/blast.h
 * @brief write the header of a @link BlastRecord @endlink to file
 * @signature int writeHeader(stream, blastRecord, blastFormatTag)
 * @signature int writeHeader(stream, qryId, dbname, [matchCount,] blastFormatTag[, labels...])
 *
 * This function writes the header of a record if @link BlastFormatFile @endlink
 * is TABULAR_WITH_HEADER (for TABULAR this is a no-op). Custom column labels
 * can be specified with the variadic fields argument (just pass a printable
 * argument for each label); use this if you use custom columns.
 *
 * If you use this function, you also need to print your matches manually with
 * @link BlastMatch#writeMatch @endlink .
 *
 * @param stream    The file to write to (FILE, fstream, @link Stream @endlink ...)
 * @param blastRecord    The @link BlastRecord @endlink whose header you want to print.
 * @param qryId     Alternatively the full ID of the query sequence (e.g. FASTA
 * one-line description)
 * @param dbName    The name of the database (NOT the ID of a subject sequence,
 * but the name of the database, e.g. "NCBI NR")
 * @param matchCount The amount of matches in the record. Mandatory parameter
 * for BlastFormatGeneration::BLAST_PLUS, optional for ::BLAST
 * @param blastFormatTag The @link BlastFormat @endlink specifier.
 * @param labels...   Optional custom column labels
 *
 * @return non-zero on I/O-error
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 */

template <typename TStream, typename TString1, typename TString2,
          BlastFormatProgram p, BlastFormatGeneration g>
inline int
_writeHeaderWithoutFields(TStream & stream,
                          TString1 const & qryId,
                          TString2 const & dbName,
                          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                      p,
                                      g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,p,g> TFormat;

    int ret = streamPut(stream, "# ");
    if (ret)
        return ret;
    ret = streamPut(stream, _programTagToString(TFormat()));
    if (ret)
        return ret;
    ret = streamPut(stream, " I/O Module of SeqAn-");
    if (ret)
        return ret;
    ret = streamPut(stream, SEQAN_VERSION_MAJOR);
    if (ret)
        return ret;
    ret = streamPut(stream, '.');
    if (ret)
        return ret;
    ret = streamPut(stream, SEQAN_VERSION_MINOR);
    if (ret)
        return ret;
    ret = streamPut(stream, '.');
    if (ret)
        return ret;
    ret = streamPut(stream, SEQAN_VERSION_PATCH);
    if (ret)
        return ret;
    ret = streamPut(stream, " (http://www.seqan.de)\n# Query: ");
    if (ret)
        return ret;
    ret = streamPut(stream, qryId);
    if (ret)
        return ret;
    ret = streamPut(stream, "\n# Database: ");
    if (ret)
        return ret;
    ret = streamPut(stream, dbName);
    if (ret)
        return ret;
    ret = streamPut(stream, '\n');
    return ret;
}

// default fields
template <typename TStream,
          typename TString,
          typename TString2,
          typename TNumeric,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
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

    int ret = _writeHeaderWithoutFields(stream, qryId, dbName, TFormat());
    if (ret)
        return ret;

    ret = streamPut(stream, "# Fields: ");
    if (ret)
        return ret;

    ret = streamPut(stream, _defaultFields(TFormat()));
    if (ret)
        return ret;

    ret = streamPut(stream, '\n');
    if (ret)
        return ret;

    if (g == BlastFormatGeneration::BLAST_PLUS)
    {
        ret = streamPut(stream, "# ");
        if (ret)
            return ret;
        ret = streamPut(stream, hitCount);

        if (ret)
            return ret;

        ret = streamPut(stream, " hits found\n");
    }
    return ret;
}

// custom fields
template <typename TStream,
          typename TString,
          typename TString2,
          typename TNumeric,
          typename TField,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
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

    int ret = _writeHeaderWithoutFields(stream, qryId, dbName, TFormat());
    if (ret)
        return ret;

    ret = streamPut(stream, "# Fields: ");
    if (ret)
        return ret;

    ret = streamPut(stream, field1);
    if (ret)
        return ret;

    ret = _writeFields(stream, TFormat(), fields...);
    if (ret)
        return ret;

    ret = streamPut(stream, '\n');
    if (ret)
        return ret;

    if (g == BlastFormatGeneration::BLAST_PLUS)
    {
        ret = streamPut(stream, "# ");
        if (ret)
            return ret;
        ret = streamPut(stream, hitCount);

        if (ret)
            return ret;

        ret = streamPut(stream, " hits found\n");
    }
    return ret;
}

// BlastFormatGeneration::BLAST works without hitCount aswell
template <typename TStream,
          typename TString1,
          typename TString2,
          BlastFormatProgram p,
          typename std::enable_if<IsSequence<TString1>::VALUE>::type = 0>
inline int
writeHeader(TStream & stream,
            TString1 const & qryId,
            TString2 const & dbName,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST> TFormat;

    return writeHeader(stream,
                       qryId,
                       dbName,
                       0,
                       TFormat());
}

// BlastFormatGeneration::BLAST works without hitCount aswell, custom fields
template <typename TStream,
          typename TString,
          typename TString2,
          typename TField,
          typename... TFields,
          BlastFormatProgram p>
inline int
writeHeader(TStream & stream,
            TString const & qryId,
            TString2 const & dbName,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST> const & /**/,
            TField const & field1,
            TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        BlastFormatGeneration::BLAST> TFormat;

    return writeHeader(stream,
                       qryId,
                       dbName,
                       0,
                       TFormat(),
                       field1,
                       fields...);
}

// Record parameter BLAST_PLUS
template <typename TStream,
          typename TRecord,
          typename TSpec,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
writeHeader(TStream & stream,
            TRecord const & record,
            BlastDbSpecs<TSpec> const & dbSpecs,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g>  TFormat;
    return writeHeader(stream, record.qId, dbSpecs.dbName,
                       record.matches.size(),TFormat());
}

// NO-OPS for TABULAR
template <typename TStream, typename TString1, typename TString2,
          BlastFormatProgram p, BlastFormatGeneration g>
constexpr int
writeHeader(TStream &, TString1 const &, TString2 const &,
            BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
    return 0;
}

template <typename TStream,
          typename T,
          typename T2,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
constexpr int
writeHeader(TStream & /**/,
            T const & /**/,
            T2 const & /**/,
            BlastFormat<BlastFormatFile::TABULAR,p,g> const & /*tag*/,
            TFields const & ... /**/)
{
    return 0;
}

template <typename TStream,
          typename T,
          typename T2,
          typename TNumeric,
          typename... TFields,
          BlastFormatProgram p,
          BlastFormatGeneration g>
constexpr int
writeHeader(TStream & /**/,
            T const & /**/,
            T2 const & /**/,
            TNumeric const /**/,
            BlastFormat<BlastFormatFile::TABULAR,p,g> const & /*tag*/,
            TFields const & ... /**/)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function writeMatch()
// ----------------------------------------------------------------------------

/*!
 * @fn BlastMatch#writeMatch
 * @headerfile seqan/blast.h
 * @brief write a @link BlastMatch @endlink to file
 * @signature int writeMatch(stream, blastMatch, blastFormatTag)
 * @signature int writeMatch(stream, blastFormatTag, columns...)
 *
 * This function writes a single match if @link BlastFormatFile @endlink
 * is TABULAR_WITH_HEADER or TABULAR. The first signature causes default
 * columns to be printed, which is described <a href="https://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/blastall/blastall_node93.html">
 * here</a>. Please note that BLAST is 1-indexed and considers the last position
 * to be the back, not the end, i.e. last one included in a match/sequence/...,
 * not the one behind it (as SeqAn does); this functions corrects for both
 * these bahaviours, so you don't have to. Additionally, based on your
 * @link BlastFormat @endlink, positions are transformed back to DNA space, if
 * translation has taken place.
 * Please note also that query and subject IDs are truncated at the first space
 * character in NCBI BLAST, this is also done by default here.
 *
 * The second signature allows an arbitrary amount of and
 * arbitrary typed columns to be printed. If you do this and you use the
 * TABULAR_WITH_HEADER format, you should also print custom column labels
 * with @link BlastRecord#writeHeader @endlink. Please note that in this case
 * the adjustments mentioned above have to be done by yourself (if you use
 * the mentioned fields and you want compatability).
 *
 * Many guides recommend always printing the default columns and using only
 * additional columns with custom data. In this case you still have to use the
 * second signature and correct for BLAST's conventions.
 *
 * @param stream    The file to write to (FILE, fstream, @link Stream @endlink ...)
 * @param blastMatch    The @link BlastMatch @endlink you wish to print.
 * @param blastFormatTag The @link BlastFormat @endlink specifier.
 * @param columns...   Custom columns
 *
 * @return non-zero on I/O-error
 * @see BlastFormat
 * @see BlastRecord
 * @see BlastRecord#writeRecord
 * @see BlastRecord#writeHeader
 */


// Functions for arbitrary number and typed fields
template <typename TStream, typename TField, typename... TFields,
          BlastFormatProgram p, BlastFormatGeneration g>
inline int
writeMatch(TStream & stream,
           BlastFormat<BlastFormatFile::TABULAR,
                       p,
                       g> const & /*tag*/,
            TField const & field1,
            TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    int ret = streamPut(stream, field1);
    if (ret)
        return ret;

    ret = _writeFields(stream, TFormat(), fields...);
    if (ret)
        return ret;
    return streamPut(stream, '\n');
}

// BlastTabHdr Record equal to BlastTab Record
template <typename TStream, typename TField, typename... TFields,
          BlastFormatProgram p, BlastFormatGeneration g>
inline int
writeMatch(TStream & stream,
           BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                       p,
                       g> const & /*tag*/,
            TField const & field1,
            TFields const & ... fields)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    return writeMatch(stream, TFormat(), field1, fields... );
}

template <typename TStream, typename TBlastMatch,
          BlastFormatProgram p, BlastFormatGeneration g>
inline int
writeMatch(TStream & stream, TBlastMatch const & match,
           BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    typedef decltype(match.qStart) TPos;

    TPos effectiveQStart    = match.qStart;
    TPos effectiveQEnd      = match.qEnd;
    TPos effectiveSStart    = match.sStart;
    TPos effectiveSEnd      = match.sEnd;

    _untranslatePositions(effectiveQStart, effectiveQEnd, match.qFrameShift,
                          match.qLength, QHasRevComp<TFormat>(),
                          QHasFrames<TFormat>());
    _untranslatePositions(effectiveSStart, effectiveSEnd, match.sFrameShift,
                          match.sLength, SHasRevComp<TFormat>(),
                          SHasFrames<TFormat>());

    return writeMatch(stream,
                      TFormat(),
                      prefix(match.qId, _firstOcc(match.qId, ' ')),
                      prefix(match.sId, _firstOcc(match.sId, ' ')),
                      double(match.identities) * 100 / match.aliLength,
                      match.aliLength,
                      match.mismatches,
                      match.gapOpenings,
                      effectiveQStart,
                      effectiveQEnd,
                      effectiveSStart,
                      effectiveSEnd,
                      match.eValue,
                      match.bitScore);
}

template <typename TStream, typename TBlastMatch,
          BlastFormatProgram p, BlastFormatGeneration g>
inline int
writeMatch(TStream & stream, TBlastMatch const & match,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const &/**/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    return writeMatch(stream, match, TFormat());
}

// ----------------------------------------------------------------------------
// Function writeTop()
// ----------------------------------------------------------------------------

// template <typename TStream,
//           typename TDbSpecs,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline int
// writeTop(TStream                    & /**/,
//          TDbSpecs             const & /**/,
//          BlastFormat<BlastFormatFile::TABULAR, p, g> const & /*tag*/)
// {
//     return 0;
// }
// 
// template <typename TStream,
//           typename TDbSpecs,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline int
// writeTop(TStream                    & /**/,
//          TDbSpecs             const & /**/,
//          BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> const & /*tag*/)
// {
//     return 0;
// }

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatFile f,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
_writeRecordImplTab(TStream                    & stream,
                    TRecord              const & record,
                    TDbSpecs             const & dbSpecs,
                    BlastFormat<f, p, g> const & /*tag*/)
{
    typedef BlastFormat<f, p, g> TFormat;

    //TODO if debug, do lots of sanity checks on record

    //NOOP for TABULAR
    int ret = writeHeader(stream, record, dbSpecs, TFormat());
    if (ret)
        return ret;
    for (auto it = record.matches.begin(); it != record.matches.end(); ++it)
    {
        //SOME SANITY CHECKS
        SEQAN_ASSERT_EQ(CharString(record.qId), CharString(it->qId));

        ret = writeMatch(stream, *it, TFormat());
        if (ret)
            return ret;
    }
    return 0;
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            BlastFormat<BlastFormatFile::TABULAR,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR, p, g> TFormat;
    return _writeRecordImplTab(stream, record, dbSpecs, TFormat());
}

template <typename TStream,
          typename TRecord,
          typename TDbSpecs,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
writeRecord(TStream             & stream,
            TRecord       const & record,
            TDbSpecs      const & dbSpecs,
            BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p, g> TFormat;
    return _writeRecordImplTab(stream, record, dbSpecs, TFormat());
}

// ----------------------------------------------------------------------------
// Function writeBottom()
// ----------------------------------------------------------------------------


// template <typename TStream,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline int
// writeBottom(TStream                           & /**/,
//             TDbSpecs                    const & /**/,
//             TBlastScoringAdapater       const & /**/,
//             BlastFormat<BlastFormatFile::TABULAR, p,g> const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
//     return 0;
// }
// 
// template <typename TStream,
//           typename TDbSpecs,
//           typename TBlastScoringAdapater,
//           BlastFormatProgram p,
//           BlastFormatGeneration g>
// inline int
// writeBottom(TStream                           & /**/,
//             TDbSpecs                    const & /**/,
//             TBlastScoringAdapater       const & /**/,
//             BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER, p,g> const & /*tag*/)
// {
//     //TODO check if this is really NO-OP forr this format
//     return 0;
// }

} // namespace seqan
#endif // header guard
