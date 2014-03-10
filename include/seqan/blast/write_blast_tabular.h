// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
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

#include <sstream>

#include <seqan/basic.h>
#include <seqan/version.h>

#include <seqan/blast/blast_base.h>

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

#ifdef SEQAN_CPP11
template <BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
constexpr
const char *
_seperatorString(BlastFormat<BlastFormatOptions::TabularWithHeader,
                              p,
                              g> const & /*tag*/)
{
    return ", ";
}

template <BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
constexpr
const char *
_seperatorString(BlastFormat<BlastFormatOptions::Tabular,
                              p,
                              g> const & /*tag*/)
{
    return "\t";
}


#endif //SEQAN_CPP11

// ============================================================================
// Functions
// ============================================================================

#ifdef SEQAN_CPP11
// ----------------------------------------------------------------------------
// Helper functions for printing n columns (requires C++11 variadic templates)
// ----------------------------------------------------------------------------


template <typename TStream, typename TTag>
inline int
_writeFields(TStream & stream, TTag const & /*tag*/)
{
    return streamPut(stream, '\n');
}

template <typename TStream, typename TField, typename... TFields, typename TTag>
inline int
_writeFields(TStream & stream,
             TTag const & /*tag*/,
             TField const & field1, const TFields&... fields)
{
    int ret = streamPut(stream,  _seperatorString(TTag()));
    if (ret)
        return ret;

    ret = streamPut(stream, field1);
    if (ret)
        return ret;

    return _writeFields(stream, TTag(), fields... );
}
#endif //SEQAN_CPP11

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

/**
.Function.BLAST I/O#writeHeader
..signature:int writeHeader(stream, query_id, db_name, [fields,] BlastFormat)
..signature:int writeHeader(stream, query_id, db_name, BlastFormat[, field1, ... fieldN])
..param.stream:The stream to write to.
...type:Concept.Stream
..param.query_id:ID of the query sequence.
...type:nolink:String-type
..param.db_name:Name of the database / file name.
...type:nolink:String-type
..param.fields:A StringSet with column identifiers to print
...type;StringSet
..param.BlastFormat: The format tag, note that BlastFormat must be further specified
...type:Class.BlastFormat
..param.fieldN:Manually supply headers of columns, defaults to Blast's defaul 12 columns (only supported with C++11)
...type:nolink:String-type
..remarks:For BlastFormat-types that have no header, this is a NOOP.
..include:seqan/blast.h
*/

template <typename TStream, typename TString1, typename TString2,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeHeader(TStream &, TString1 const &, TString2 const &,
            BlastFormat<BlastFormatOptions::Tabular, p, g> const & /*tag*/)
{
    return 0;
}


template <typename TStream, typename TString1, typename TString2,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
_writeHeaderWithoutFields(TStream & stream,
                          TString1 const & qryId, TString2 const & dbName,
                          BlastFormat<BlastFormatOptions::TabularWithHeader,
                                      p,
                                      g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,p,g> TFormat;

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

// default case
template <typename TStream, typename TString1, typename TString2,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeHeader(TStream & stream,
            TString1 const & qryId, TString2 const & dbName,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,p,g> TFormat;

    int ret = _writeHeaderWithoutFields(stream, qryId, dbName, TFormat());
    if (ret)
        return ret;

    streamPut(stream, "# ");
    if (ret)
        return ret;

    ret = streamPut(stream, _defaultFields(TFormat()));
    if (ret)
        return ret;

    return streamPut(stream, '\n');
}

// 
template <typename TStream, typename TqId, typename TdbName,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeHeader(TStream & stream,
            TqId const & qryId, TdbName const & dbName,
            StringSet<CharString> const & fields,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader,p,g> TFormat;

    int ret = _writeHeaderWithoutFields(stream, qryId, dbName, TFormat());
    if (ret)
        return ret;

    ret = streamPut(stream, "# Fields: ");
    if (ret)
        return ret;

    for (unsigned i = 0; i < length(fields); ++i)
    {
        ret = streamPut(stream, fields[i]);
        if (ret)
            return ret;

        if (i == length(fields) -1)
            ret = streamPut(stream, '\n');
        else
            ret = streamPut(stream, _seperatorString(TFormat()));

        if (ret)
            return ret;
    }

    return 0;
}

#ifdef SEQAN_CPP11
// Functions for arbitrary number and typed fields

template <typename TStream, typename TString,
          typename TField, typename... TFields,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeHeader(TStream & stream,
            TString const & qryId, TString const & dbName,
            BlastFormat<m,p,g> const & /*tag*/,
            TField const & field1,
            const TFields&... fields)
{
    return 0;
}

template <typename TStream, typename TString,
          typename TField, typename... TFields
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeHeader(TStream & stream,
            TString const & qryId, TString const & dbName,
            BlastFormat<m,p,g> const & /*tag*/,
            TField const & field1,
            const TFields&... fields)
{
    typedef BlastFormat<m,p,g> Format;

    int ret = _writeHeaderWithoutFields(stream, qryId, dbName, Format());
    if (ret)
        return ret;

    ret = streamPut(stream, "# Fields: ");
    if (ret)
        return ret;

    ret = streamPut(stream, field1);

    return _writeFields(stream, Format(), fields...);
}

#endif //SEQAN_CPP11

// ----------------------------------------------------------------------------
// Function writeMatch()
// ----------------------------------------------------------------------------

/**
.Function.BLAST I/O#writeMatch
..signature:int writeMatch(stream, query_id, subject_id, percentIdent, ali_length, num_mismatches, gap_openings, qStart, qEnd, sStart, sEnd, eval, bitScore, BlastFormat)
..signature:int writeMatch(stream, fields, BlastFormat)
..signature:int writeMatch(stream, query_id, subject_id, fields, BlastFormat)
..signature:int writeMatch(stream, BlastFormat, field1,[ ... fieldN])
..param.stream:The stream to write to.
...type:Concept.Stream
..param.query_id:ID of the query sequence.
...type:nolink:String-type
..param.subject_id:ID of the database sequence.
...type:nolink:String-type
..param.num_identies:number of identical position in the alignment
...type:nolink:double
..param.ali_length:length of alignment
...type:nolink:unsigned
..param.num_mismatches: number of mismatched positions in alignment
...type:nolink:unsigned
..param.gap_openings: number of gap_openings in alignment
...type:nolink:unsigned
..param.qStart: begin position of the alignment on the query
...type:nolink:unsigned
..param.qEnd: end position of the alignment on the query
...type:nolink:unsigned
..param.sStart: begin position of the alignment on the subject
...type:nolink:unsigned
..param.sEnd: end position of the alignment on the subject
...type:nolink:unsigned
..param.fields:A String with fields to print (note that these have to be of
the same type, e.g. String(Set) of Strings, or String of Double)
..type:String
..param.BlastFormat: The format tag, note that BlastFormat must be further specified
...type:Class.BlastFormat
..param.fieldN:parameters of differently typed values to print (only supported with C++11)
...type:nolink:anything printable by @streamPut@
..include:seqan/blast.h
*/

template <typename TStream, typename TqId, typename TsId,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream,
            TqId const & qId, TsId const & sId, unsigned const & num_identies,
            unsigned const & ali_length, unsigned const & num_mismatches,
            unsigned const & gap_openings,
            unsigned long const & qStart, unsigned long const & qEnd,
            unsigned long const & sStart, unsigned long const & sEnd,
            double const & eval, double const & bitScore,
            BlastFormat<BlastFormatOptions::Tabular, p, g> const & /*tag*/)
{
    int ret = streamPut(stream, prefix(qId, _firstOcc(qId, ' ')));
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, prefix(sId, _firstOcc(sId, ' ')));
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, double(num_identies) * 100 / ali_length );
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, ali_length);
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, num_mismatches);
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, gap_openings);
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, qStart);
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, qEnd);
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, sStart);
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, sEnd);
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, eval);
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, bitScore);
    if (ret)
        return ret;
    ret = streamPut(stream, '\n');
    return ret;

//     std::stringstream ss (std::stringstream::in | std::stringstream::out);
//     ss << _untilFirstWhitespace(qId) << '\t' <<_untilFirstWhitespace(sId)<< '\t'
//        << double(num_identies) * 100 / ali_length << '\t'
//        << ali_length << '\t' << num_mismatches << '\t' << gap_openings << '\t'
//        << qStart << '\t' << qEnd << '\t' << sStart << '\t' << sEnd << '\t'
//        << eval << '\t' << bitScore << '\n';
// 
//     ss.seekg(0, std::ios::beg);
//     return streamPut(stream, ss);
}

template <typename TStream, typename TqId, typename TsId, typename TField,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream,
            TqId const & qId, TsId const & sId,
            String<TField> const & fields,
            BlastFormat<BlastFormatOptions::Tabular, p, g> const & /*tag*/)
{
    int ret = streamPut(stream, prefix(qId, _firstOcc(qId, ' ')));
    if (ret)
        return ret;
    ret = streamPut(stream, '\t');
    if (ret)
        return ret;
    ret = streamPut(stream, prefix(sId, _firstOcc(sId, ' ')));
    if (ret)
        return ret;

    for (int i = 0; i < length(fields); ++i)
    {
        ret = streamPut(stream, '\t');
        if (ret)
            return ret;

        ret = streamPut(stream, fields[i]);
        if (ret)
            return ret;
    }
    ret = streamPut(stream, '\n');
    return ret;
}

template <typename TStream, typename TField,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream,
            String<TField> const & fields,
            BlastFormat<BlastFormatOptions::Tabular, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeMatch(stream, "", "", fields, TFormat());
}


// BlastTabHdr Record equal to BlastTab Record
template <typename TStream, typename TqId, typename TsId,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream,
            TqId const & qId, TsId const & sId, double const & percentIdent,
            unsigned const & ali_length, unsigned const & num_mismatches,
            unsigned const & gap_openings,
            unsigned long const & qStart, unsigned long const & qEnd,
            unsigned long const & sStart, unsigned long const & sEnd,
            double const & eval, double const & bitScore,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeMatch(stream, qId, sId, percentIdent, ali_length,
                       num_mismatches, gap_openings, qStart, qEnd, sStart,
                       sEnd, eval, bitScore, TFormat());
}

// BlastTabHdr Record equal to BlastTab Record
template <typename TStream, typename TqId, typename TsId, typename TField,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream,
            TqId const & qId, TsId const & sId,
            String<TField> const & fields,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeMatch(stream, qId, sId, fields, TFormat());
}

// BlastTabHdr Record equal to BlastTab Record
template <typename TStream, typename TField,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream,
            String<TField> const & fields,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeMatch(stream, "", "", fields, TFormat());
}


#ifdef SEQAN_CPP11
// Functions for arbitrary number and typed fields

template <typename TStream, typename TField, typename... TFields,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream,
            BlastFormat<BlastFormatOptions::Tabular,
                        p,
                        g> const & /*tag*/)
            TField const & field1,
            const TFields&... fields)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    int ret = streamPut(stream, field1);
    if (ret)
        return ret;

    return _writeFields(stream, TFormat(), fields...);
}

// BlastTabHdr Record equal to BlastTab Record
template <typename TStream, typename TField, typename... TFields,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/,
            TField const & field1,
            const TFields&... fields)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeMatch(stream, TFormat(), field1, fields... );
}
#endif //SEQAN_CPP11


template <typename TStream, typename TBlastMatch,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream, TBlastMatch const & match,
            BlastFormat<BlastFormatOptions::Tabular, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeMatch(stream,
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
                      match.eVal,
                      match.bitScore,
                      TFormat());
}

template <typename TStream, typename TBlastMatch,
          BlastFormatOptions::Program p, BlastFormatOptions::Generation g>
inline int
writeMatch(TStream & stream, TBlastMatch const & match,
            BlastFormat<BlastFormatOptions::TabularWithHeader, p, g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return writeMatch(stream,
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
                      match.eVal,
                      match.bitScore,
                      TFormat());
}


// ----------------------------------------------------------------------------
// Function writeTop()
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TString,
          typename TNum1,
          typename TNum2,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeTop(TStream             & /**/,
         TString       const & /**/,
         TNum1         const & /**/,
         TNum2         const & /**/,
            BlastFormat<BlastFormatOptions::Tabular,
                        p,
                        g> const & /*tag*/)
{
    //TODO check if this is really empty
    return 0;
}

template <typename TStream,
          typename TString,
          typename TNum1,
          typename TNum2,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeTop(TStream             & /**/,
         TString       const & /**/,
         TNum1         const & /**/,
         TNum2         const & /**/,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    //TODO check if this is really empty
    return 0;
}



// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

template <typename TStream,
          typename TRecord,
          BlastFormatOptions::M m,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
_writeRecordImplTab(TStream                    & stream,
                    TRecord              const & record,
                    BlastFormat<m, p, g> const & /*tag*/)
{
    typedef BlastFormat<m, p, g> TFormat;

    //TODO if debug, do lots of sanity checks on record

    //NOOP for Tabular
    int ret = writeHeader(stream, record.qId, record.dbName, TFormat());
    if (ret)
        return ret;

    for (auto const & bm : record.matches)
    {
        ret = writeMatch(stream, bm, TFormat());
        if (ret)
            return ret;
    }
    return 0;
}

template <typename TStream,
          typename TRecord,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeRecord(TStream             & stream,
            TRecord       const & record,
            BlastFormat<BlastFormatOptions::Tabular,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::Tabular, p, g> TFormat;
    return _writeRecordImplTab(stream, record, TFormat());
}

template <typename TStream,
          typename TRecord,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeRecord(TStream             & stream,
            TRecord       const & record,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    typedef BlastFormat<BlastFormatOptions::TabularWithHeader, p, g> TFormat;
    return _writeRecordImplTab(stream, record, TFormat());
}


// ----------------------------------------------------------------------------
// Function writeBottom()
// ----------------------------------------------------------------------------



template <typename TStream,
          typename TValue, typename TSpec,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeBottom(TStream                      & /**/,
            Score<TValue, TSpec>   const & /**/,
            BlastFormat<BlastFormatOptions::Tabular,
                        p,
                        g> const & /*tag*/)
{
    //TODO check if this is really empty
    return 0;
}

template <typename TStream,
          typename TValue, typename TSpec,
          BlastFormatOptions::Program p,
          BlastFormatOptions::Generation g>
inline int
writeBottom(TStream                      & /**/,
            Score<TValue, TSpec>   const & /**/,
            BlastFormat<BlastFormatOptions::TabularWithHeader,
                        p,
                        g> const & /*tag*/)
{
    //TODO check if this is really empty
    return 0;
}

} // namespace seqan
#endif // header guard
