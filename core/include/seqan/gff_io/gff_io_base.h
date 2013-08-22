// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_
#define CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Gff
// ----------------------------------------------------------------------------

/*!
 * @defgroup GffIO GFF I/O
 * @brief I/O functionality for the GFF and GTF file formats.
 */

/*!
 * @tag GffIO#Gff
 * @brief Tag for selecting the GFF format.
 *
 * @signature typedef Tag<TagGff_> Gff;
 */

/**
.Tag.File Format.tag.Gff:
    Gff annotation file.
..include:seqan/gff_io.h
*/

// TODO(singer): const should be non const, but is const elsewhere
struct TagGff_;
typedef Tag<TagGff_> Gff;

// ----------------------------------------------------------------------------
// Tag Gtf
// ----------------------------------------------------------------------------

/*!
 * @tag GffIO#Gtf
 * @brief Tag for selecting the GTF format.
 *
 * @signature typedef Tag<TagGtf_> Gtf;
 */

/**
.Tag.File Format.tag.Gtf:
    Gtf annotation file.
..include:seqan/gff_io.h
*/

// TODO(singer): const should be non const, but is const elsewhere
struct TagGtf_;
typedef Tag<TagGtf_> Gtf;

// ----------------------------------------------------------------------------
// Class GffRecord
// ----------------------------------------------------------------------------

/*!
 * @class GffRecord
 * @headerfile <seqan/gff_io.h>
 * @brief Represent a record from a Gff file.
 *
 * @signature class GffRecord;
 *
 * @var __int32 GffRecord::INVALID_POS
 * @brief Static member with invalid/sentinel position value.
 *
 * @var __int32 GffRecord::INVALID_IDX
 * @brief Static member with invalid/sentinel rID value.
 *
 * @fn GffRecord::INVALID_SCORE
 * @brief Returns invalid score (NaN float value).
 *
 * @signature float GffRecord::INVALID_SCORE();
 * 
 * @var CharString GffRecord::ref
 * @brief The sequence name of the record.
 *
 * @var __int32 GffRecord::rID
 * @brief Integer representing ref, defaults to INVALID_IDX.
 * 
 * @var CharString GffRecord::source
 * @brief The source of the record.
 * 
 * @var CharString GffRecord::type
 * @brief The type of the record.
 * 
 * @var __int32 GffRecord::beginPos
 * @brief The begin position of the record.
 * 
 * @var __int32 GffRecord::endPos
 * @brief The end position of the record.
 * 
 * @var float GffRecord::score
 * @brief The score of the record.
 * 
 * @var char GffRecord::strand
 * @brief The strand the record belongs to.
 * 
 * @var char GffRecord::phase
 * @brief The phase of the record.
 * 
 * @section Remarks
 * 
 * For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.  The
 * phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of
 * this feature to reach the first base of the next codon
 * 
 * @var TCharStringSet GffRecord::tagName
 * @brief The names of the attributes of the record, StringSet of CharString.
 * 
 * @var TCharStringSet GffRecord::tagValue
 * @brief The values of the attributes of the record, StringSet of CharString.
 * 
 * @section Remarks
 * 
 * For each value there is a name associated in GffRecord::tagName.
 *
 * @section Remarks
 * 
 * For each name there is a value associated in GffRecord::tagValue.
 */

/**
.Class.GffRecord
..cat:BAM I/O
..summary:Represent a record from a Gff file.
..include:seqan/gff_io.h

.Memvar.GffRecord#INVALID_POS
..class:Class.GffRecord
..summary:Static member with invalid/sentinel position value.
..type:nolink:$__uint32$

.Memvar.GffRecord#INVALID_SCORE
..class:Class.GffRecord
..summary:Static member with invalid score value.
..type:nolink:$float$

.Memvar.GffRecord#ref
..class:Class.GffRecord
..summary:The sequence id of the record.
..type:Shortcut.CharString

.Memvar.GffRecord#source
..class:Class.GffRecord
..summary:The source of the record.
..type:Shortcut.CharString

.Memvar.GffRecord#type
..class:Class.GffRecord
..summary:The type of the record.
..type:Shortcut.CharString

.Memvar.GffRecord#beginPos
..class:Class.GffRecord
..summary:The begin position of the record.
..type:nolink:$__uint32$

.Memvar.GffRecord#endPos
..class:Class.GffRecord
..summary:The end position of the record.
..type:nolink:$__uint32$

.Memvar.GffRecord#score
..class:Class.GffRecord
..summary:The score of the record.
..type:nolink:$float$

.Memvar.GffRecord#strand
..class:Class.GffRecord
..summary:The strand the record belongs to.
..type:nolink:$char$

.Memvar.GffRecord#phase
..class:Class.GffRecord
..summary:The phase of the record.
..remarks:For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon
..type:nolink:$char$

.Memvar.GffRecord#tagName
..class:Class.GffRecord
..summary:The names of the attributes of the record.
..type:Class.StringSet
..remarks:For each name there is a value associated in $Memvar.GffRecord#tagValue$

.Memvar.GffRecord#tagValue
..class:Class.GffRecord
..summary:The values of the attributes of the record.
..type:Class.StringSet
..remarks:For each value there is a name associated in $Memvar.GffRecord#tagName$
*/

struct GffRecord
{
    static __int32 const INVALID_POS = 2147483647;  // TODO(singer): Should be MaxValue<__int32>::VALUE, but that is not a constant expression :(
    static __int32 const INVALID_IDX = -1;

    // The member descriptions are taken from: http://gmod.org/wiki/GFF

    // TODO(singer): Maybe use a I/O context object and store ids as integers
    // The ID of the landmark used to establish the coordinate system for the current feature.
    String<char> ref;
    int rID;

    // The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature.
    String<char> source;

    // The type of the feature
    String<char> type;

    // A list of feature attributes in the format tag=value.
    StringSet<String<char> > tagName;
    StringSet<String<char> > tagValue;

    // The start and end of the feature, in 1-based integer coordinates, relative to the landmark given in column 1
    __uint32 beginPos;
    __uint32 endPos;

    // The score of the feature
    float score;

    // The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for features that are not stranded.
    char strand;

    // For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.
    // The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon.
    char phase;

    static float INVALID_SCORE()
    {
        union
        {
            __uint32 u;
            float f;
        } tmp;
        tmp.u = 0x7F800001;
        return tmp.f;
    }

    GffRecord() :
        rID(INVALID_IDX), beginPos(-1), endPos(-1), score(INVALID_SCORE()),
        strand('.'), phase('.')
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _parseReadGffKeyValue
// ----------------------------------------------------------------------------

template <typename TReader, typename TKeyString, typename TValueString>
inline int
_parseReadGffKeyValue(TValueString & outValue, TKeyString & key, TReader & reader)
{
    char c = value(reader);
    if (c == ' ' || c == '\t' || c == '\n' || c == '=')
        return 1;  // Key cannot be empty.

    for (; !atEnd(reader); goNext(reader))
    {
        c = value(reader);
        if (c == ' ' || c == '\t' || c == '\n' || c == '=' || c == ';')
            break;
        appendValue(key, c);
    }
    if (!atEnd(reader) && value(reader) == ';')
    {
        goNext(reader);
        return 0;
    }
    if (!atEnd(reader) && (value(reader) == '\r' || value(reader) == '\n'))
        return 0;

    if (skipWhitespaces(reader) != 0)
        return 1;

    if (value(reader) == '=')
    {
        goNext(reader);
        if (skipWhitespaces(reader) != 0)
            return 1;

        if (atEnd(reader))
            return 1;
    }

    if (value(reader) == '"')
    {
        // Handle the case of a string literal.

        goNext(reader);
        // Append all characters in the literal to outValue until the first i
        // line break or the closing '"'.
        for (; !atEnd(reader); goNext(reader))
        {
            if (value(reader) == '\n')
                return 1;

            if (value(reader) == '"')
            {
                goNext(reader);
                break;
            }
            appendValue(outValue, value(reader));
        }
        // Go over the trailing semicolon and any trailing space.
        while (!atEnd(reader) && (value(reader) == ';' || value(reader) == ' '))
            goNext(reader);
    }
    else
    {
        // Handle the non literal case.

        // Read until the first semicolon, return at whitespace.
        for (; !atEnd(reader); goNext(reader))
        {
            if (value(reader) == ';' || value(reader) == '\n' || value(reader) == '\r')
                break;
            appendValue(outValue, value(reader));
        }
        // Skip semicolon and spaces if any.
        while (!atEnd(reader) && (value(reader) == ';' || value(reader) == ' '))
            goNext(reader);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

/*!
 * @fn GffRecord#clear
 * @brief Reset a @link GffRecord @endlink object.
 * 
 * @signature void clear(record);
 * 
 * @param[in,out] record The GffRecord to reset.
 */

/**
.Function.GffRecord#clear
..class:Class.GffRecord
..cat:Input/Output
..signature:clear(record)
..param.record:The @Class.GffRecord@ to reset.
...type:Class.GffRecord
..summary:Reset a @Class.GffRecord@ object.
..include:seqan/gff_io.h
*/

inline void clear(GffRecord & record)
{
    clear(record.ref);
    clear(record.source);
    clear(record.type);
    clear(record.tagName);
    clear(record.tagValue);
}

// ----------------------------------------------------------------------------
// Function readRecord
// ----------------------------------------------------------------------------

/*!
 * @fn GffIO#readRecord
 * @brief Read one GFF/GTF record from a SinglePassRecordReader.
 *
 * @signature int readRecord(record, reader[, context, [, tag]]);
 *
 * @param[out]    record  The GffRecord to write the results to.
 * @param[in,out] reader  The SinglePassRecordReader to use for reading.
 * @param[in,out] context The GffIOContext to use for reading.  If present then ref will be translated to rID using the
 *                        reference name store from context.
 * @param[in]     tag     The format to read from, one of Gtf and Gff.  Note that the parser transparently parses both
 *                        GFF and GTF.
 *
 * @return int A status code, 0 on success, a different value on failures.
 */

/**
.Function.GffRecord#readRecord
..class:Class.GffRecord
..cat:Input/Output
..summary:Read one gff record.
..signature:readRecord(record, reader)
..param.record:The gff record.
...type:Class.GffRecord
..param.reader: The record reader object.
...type:Class.RecordReader
..include:seqan/gff_io.h
*/

template <typename TStream, typename TRecordReaderSpec>
inline int
_readGffRecord(GffRecord & record, RecordReader<TStream, TRecordReaderSpec> & reader)
{
    clear(record);

    // read column 1: seqid
    // The letters until the first whitespace will be read.
    // Then, we skip until we hit the first tab character.
    if (readUntilTabOrLineBreak(record.ref, reader))
        return 1;
    record.rID = GffRecord::INVALID_IDX;

    if (!empty(record.ref) && record.ref[0] == '#')
    {
        if (skipLine(reader))
            return 1;

        return 1;
    }
    if (skipWhitespaces(reader))
        return 1;

    // read column 2: source
    if (readUntilTabOrLineBreak(record.source, reader))
        return 1;

    if (record.source == ".")
        clear(record.source);

    if (skipWhitespaces(reader))
        return 1;

    // read column 3: type
    if (readUntilTabOrLineBreak(record.type, reader))
        return 1;

    if (skipWhitespaces(reader))
        return 1;

    // read column 4: begin position
    String<char> temp;
    if (readDigits(temp, reader))
        return 1;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(record.beginPos, temp))
            return 1;

        --record.beginPos;  // Translate from 1-based to 0-based.
    }
    else
    {
        record.beginPos = GffRecord::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return 1;
    }
    if (skipBlanks(reader))
        return 1;

    // read column 5: end position
    clear(temp);
    if (readDigits(temp, reader))
        return 1;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(record.endPos, temp))
            return 1;
    }
    else
    {
        record.endPos = GffRecord::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return 1;
    }
    if (skipBlanks(reader))
        return 1;


    // read column 6: score
    clear(temp);
    readFloat(temp, reader);

    if (length(temp) > 0u)
    {
        if (temp != ".")
        {
            if (!lexicalCast2(record.score, temp))
                return 1;
        }
        else
        {
            record.score = GffRecord::INVALID_SCORE();
            if (skipUntilWhitespace(reader))
                return 1;
        }
    }
    else
    {
        return 1;
    }

    if (skipBlanks(reader))
        return 1;

    // read column 7: strand
    clear(temp);
    if (readUntilTabOrLineBreak(temp, reader))
        return 1;

    if (temp[0] != '-' && temp[0] != '+')
    {
        record.strand = '.';
    }
    else
    {
        record.strand = temp[0];
    }

    if (skipBlanks(reader))
        return 1;

    // read column 8: phase
    clear(temp);
    if (readUntilTabOrLineBreak(temp, reader))
        return 1;

    if (temp != "0" && temp != "1" && temp != "2")
    {
        record.phase = '.';
    }
    else
    {
        record.phase = temp[0];
    }

    if (skipBlanks(reader))
        return 1;

    // It's fine if there are no attributes and the line ends here.
    if (atEnd(reader))
        return 0;
    if (value(reader) == '\n' || value(reader) == '\r')
        return skipLine(reader);

    // read column 9: attributes
    while (!atEnd(reader))
    {

        String<char> _key;
        String<char> _value;
        // Read next key/value pair.
        if (_parseReadGffKeyValue(_value, _key, reader) != 0)
            return 1;

        appendValue(record.tagName, _key);
        appendValue(record.tagValue, _value);

        clear(_key);
        clear(_value);

        // At end of line:  Skip EOL and break.
        if (!atEnd(reader) && (value(reader) == '\r' || value(reader) == '\n'))
        {
            if (skipLine(reader) != 0)
                return 1;

            break;
        }
    }
    return 0;
}

template <typename TRecordReader>
inline int
readRecord(GffRecord & record, TRecordReader & reader, Gff /*tag*/)
{
    return _readGffRecord(record, reader);
}

template <typename TRecordReader>
inline int
readRecord(GffRecord & record, TRecordReader & reader, Gtf /*tag*/)
{
    return _readGffRecord(record, reader);
}

// TODO(singer): Needs proper documentation!!! Check the length of the stores!!!
template <typename TRecordReader, typename TContextSpec, typename TContextSpec2>
inline int
_readGffRecord(GffRecord & record, TRecordReader & reader, GffIOContext<TContextSpec, TContextSpec2> & context)
{
    // Read record with string ref from GFF file.
    int res = _readGffRecord(record, reader);
    if (res != 0)
        return res;

    // Translate ref to rID using the context.  If there is no such sequence name in the context yet then we add it.
    unsigned idx = 0;
    if (!getIdByName(nameStore(context), record.ref, idx, nameStoreCache(context)))
    {
        idx = length(nameStore(context));
        appendName(nameStore(context), record.ref, nameStoreCache(context));
    }
    record.rID = idx;

    return 0;
}

template <typename TRecordReader, typename TContextSpec, typename TContextSpec2>
inline int
readRecord(GffRecord & record, TRecordReader & reader, GffIOContext<TContextSpec, TContextSpec2> & context, Gff /*tag*/)
{
    return _readGffRecord(record, reader, context);
}

template <typename TRecordReader, typename TContextSpec, typename TContextSpec2>
inline int
readRecord(GffRecord & record, TRecordReader & reader, GffIOContext<TContextSpec, TContextSpec2> & context, Gtf /*tag*/)
{
    return _readGffRecord(record, reader, context);
}

// ----------------------------------------------------------------------------
// Function _writeSemicolonSensitive()
// ----------------------------------------------------------------------------

// This function checks if the string to be written contains a semicolon. If
// this is the case then quotes are written around the string.
// Returns false on success.

template <typename TTargetStream, typename TString>
inline bool
_writeInQuotes(TTargetStream & target, TString & temp)
{
    // TODO(jsinger): What about escaping quote chars '"'?
    if (streamWriteChar(target, '"') || streamWriteBlock(target, begin(temp, Standard()), length(temp)) < length(temp) ||
        streamWriteChar(target, '"'))
        return true;

    return false;
}

template <typename TTargetStream, typename TString, typename TMustBeQuotedFunctor>
inline bool
_writePossiblyInQuotes(TTargetStream & target, TString & source, TMustBeQuotedFunctor const &func)
{
    // TODO(jsinger): What about escaping quote chars '"'?
    typedef typename Iterator<TString, Standard>::Type TIter;

    TIter itEnd = end(source, Standard());
    for (TIter it = begin(source, Standard()); it != itEnd; ++it)
    {
        // we have a problem if the string contains a '"' or a line break
        if (*it == '\n' || *it == '"')
            return 1;

        if (func(*it))
            return _writeInQuotes(target, source);
    }
    return (streamWriteBlock(target, begin(source, Standard()), length(source)) < length(source));
}

// ----------------------------------------------------------------------------
// Function writeRecord
// ----------------------------------------------------------------------------

/*!
 * @fn GffIO#writeRecord
 * @brief Writes on GFF/GTF record to a stream.
 *
 * @signature int writeRecord(stream, record);
 *
 * @param[in,out] stream The StreamConcept to write to.
 * @param[in]     record The GffRecord to write.
 *
 * @return int A status code, 0 on success, a different value on errors.
 */

/**
.Function.GffRecord#writeRecord
..class:Class.GffRecord
..cat:Input/Output
..summary:Writes one gff record to a stream.
..signature:writeRecord(TSreamm stream, GffRecord record)
..param.stream:The output stream.
...type:Concept.StreamConcept
..param.record:The gff record.
...type:Class.GffRecord
..include:seqan/gff_io.h
*/

template <typename TFormatTag>
struct GffRecordKeyMustBeQuoted_;

template <typename TFormatTag>
struct GffRecordValueMustBeQuoted_;

// GFF quotation rules

template <>
struct GffRecordKeyMustBeQuoted_<Gff>
{
    bool operator() (char c) const
    {
        return c == ';' || c == '=';
    }
};

template <>
struct GffRecordValueMustBeQuoted_<Gff> :
    GffRecordKeyMustBeQuoted_<Gff> {};

// GTF quotation rules

template <>
struct GffRecordKeyMustBeQuoted_<Gtf>
{
    bool operator() (char c) const
    {
        return c == ';' || c == ' ';
    }
};

template <>
struct GffRecordValueMustBeQuoted_<Gtf>
{
    bool operator() (char c) const
    {
//        return c == ';' || c == ' ' || !isdigit(c);
        return !isdigit(c);     // is equivalent to the above, quote everything except integral values
    }
};



template <typename TStream, typename TTag>
inline int
_writeAttributes(TStream & stream, GffRecord const & record, TTag)
{
    const char separatorBetweenTagAndValue = (IsSameType<TTag, Gff>::VALUE)? '=' : ' ';
    for (unsigned i = 0; i < length(record.tagName); ++i)
    {
        if (i != 0)
        {
            if (streamWriteChar(stream, ';'))
                return 1;

            // In GTF files a space follows the semicolon
            if (IsSameType<TTag, Gtf>::VALUE && streamWriteChar(stream, ' '))
                return 1;
        }

        if (_writePossiblyInQuotes(stream, record.tagName[i], GffRecordKeyMustBeQuoted_<TTag>()))
            return 1;

        if (!empty(record.tagValue[i]))
        {
            if (streamWriteChar(stream, separatorBetweenTagAndValue))
                return 1;

            if (_writePossiblyInQuotes(stream, record.tagValue[i], GffRecordValueMustBeQuoted_<TTag>()))
                return 1;
        }
    }
    
    // In GTF files each (especially the last) attribute must end with a semi-colon
    if (IsSameType<TTag, Gtf>::VALUE && !empty(record.tagName) && streamWriteChar(stream, ';'))
        return 1;

    return 0;
}

template <typename TStream, typename TSeqId, typename TTag>
inline int
_writeRecordImpl(TStream & stream, GffRecord const & record, TSeqId const & ref, TTag tag)
{
    // ignore empty annotations, i.e. annotations that are 'guessed' by implicit information from their children (in GFF)
    if (empty(ref))
        return 0;

    // write column 1: seqid
    if (streamWriteBlock(stream, begin(ref, Standard()), length(ref)) != length(ref))
        return 1;

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 2: source
    if (empty(record.source))
    {
        if (streamWriteChar(stream, '.'))
            return 1;
    }
    else
    {
        if (streamWriteBlock(stream, begin(record.source, Standard()), length(record.source)) != length(record.source))
            return 1;
    }

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 3: type
    if (streamWriteBlock(stream, begin(record.type, Standard()), length(record.type)) != length(record.type))
        return 1;

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 4: begin position
    if (record.beginPos != (unsigned)-1)
    {
        if (streamPut(stream, record.beginPos + 1))
            return 1;
    }
    else
    {
        if (streamWriteChar(stream, '.'))
            return 1;
    }

    if (streamWriteChar(stream, '\t'))
        return 1;


    // write column 5: end position
    if (record.endPos != (unsigned)-1)
    {
        if (streamPut(stream, record.endPos))
            return 1;
    }
    else
    {
        if (streamWriteChar(stream, '.'))
            return 1;
    }

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 6: score
    if (record.score != record.score)
    {
        if (streamWriteChar(stream, '.'))
            return 1;
    }
    else
    {
        if (streamPut(stream, record.score))
            return 1;
    }

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 7: strand
    if (streamWriteChar(stream, record.strand))
        return 1;

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 8: phase
    if (streamWriteChar(stream, record.phase))
        return 1;

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 9: attributes
    // only until length - 1, because there is no semicolon at the end of the line

    _writeAttributes(stream, record, tag);

    if (streamWriteChar(stream, '\n'))
        return 1;

    return 0;
}

template <typename TStream, typename TTag>
inline int
writeRecord(TStream & stream, GffRecord const & record, TTag const tag)
{
    return _writeRecordImpl(stream, record, record.ref, tag);
}

template <typename TStream, typename TContextSpec, typename TContextSpec2, typename TTag>
inline int
writeRecord(TStream & stream, GffRecord const & record, GffIOContext<TContextSpec, TContextSpec2> & context, TTag const tag)
{
    if (record.rID != GffRecord::INVALID_IDX)
    {
        String<char> tempSeqId = nameStore(context)[record.rID];
        return _writeRecordImpl(stream, record, tempSeqId, tag);
    }
    return _writeRecordImpl(stream, record, record.ref, tag);
}

}  // namespace seqan

#endif  // CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_
