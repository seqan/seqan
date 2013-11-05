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

struct GffContext
{
    String<char> buffer;
};

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

template <typename TForwardIter, typename TKeyString, typename TValueString>
inline void
_parseReadGffKeyValue(TValueString & outValue, TKeyString & key, TForwardIter & iter)
{
    IsWhitespace isWhitespace;

    //TODO(singer): AssertList functor would be need
    char c = value(iter);
    if (c == ' ' || c == '\t' || c == '\n' || c == '=')
    {
        throw std::runtime_error("The key field of an attribute is empty!");
        return;  // Key cannot be empty.
    }

    for (; !atEnd(iter); goNext(iter))
    {
        c = value(iter);
        //if (IsWhitespace(c) || c == '=' || c == ';')
        if (c == '\n' || c == '\r' || c == ' ' || c == '=' || c == ';')
            break;
        appendValue(key, c);
    }
    if (!atEnd(iter) && value(iter) == ';')
    {
        skipOne(iter);
        return;
    }

    if(value(iter) == '\r' || value(iter) == '\n')
        return;

    skipUntil(iter, NotFunctor<IsWhitespace>());

    if (value(iter) == '=')
    {
        skipOne(iter);
    }

    if (value(iter) == '"')
    {
        // Handle the case of a string literal.
        skipOne(iter);
        readUntil(outValue, iter, OrFunctor<EqualsChar<'"'>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
        skipOne(iter);

        // Go over the trailing semicolon and any trailing space.
        while (!atEnd(iter) && (value(iter) == ';' || value(iter) == ' '))
            goNext(iter);
    }
    else
    {
        // Read until the first semicolon, return at whitespace.
        readUntil(outValue, iter, OrFunctor<EqualsChar<';'>, IsNewline>());

        // Skip semicolon and spaces if any.
        while (!atEnd(iter) && (value(iter) == ';' || value(iter) == ' '))
            goNext(iter);
    }
    return;
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

//TODO(singer): dont we need to reset beginPos, score ...
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

//TODO(singer): no checking if record is complete
//TODO(singer): no checking whether lexicalCast is working
template <typename TFwdIterator>
inline void
_readGffRecord(GffRecord & record, TFwdIterator & iter, GffContext & context)
{
    IsNewline isNewline;

    clear(record);
    record.rID = GffRecord::INVALID_IDX;

    skipUntil(iter, NotFunctor<OrFunctor<EqualsChar<'#'>, IsWhitespace> >());  //skip commments and empty lines

    // read column 1: seqid
    readUntil(record.ref, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    skipOne(iter);

    // read column 2: source
    readUntil(record.source, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());

    if (record.source == ".")
        clear(record.source);

    skipOne(iter);

    // read column 3: type
    readUntil(record.type, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    skipOne(iter);

    // read column 4: begin position
    readUntil(context.buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    record.beginPos = lexicalCast<__uint32>(context.buffer);
    --record.beginPos;  // Translate from 1-based to 0-based.
    skipOne(iter);

    // read column 5: end position
    clear(context.buffer);
    readUntil(context.buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    record.endPos = lexicalCast<__uint32>(context.buffer);
    skipOne(iter);

    //check if end < begin
    if (record.endPos < record.beginPos)
        throw std::runtime_error("The begin position of the record is larger than the end position!");

    // read column 6: score
    clear(context.buffer);
    readUntil(context.buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    if (context.buffer != ".")
        record.score = lexicalCast<float>(context.buffer);
    skipOne(iter, IsTab());

    // read column 7: strand
    //TODO(singer): readUntil taking a char would be good!
    //TODO(singer): readOne
    record.strand = value(iter);
    if (record.strand != '-' &&record.strand != '+')
    {
        record.strand = '.';
    }
    skipOne(iter);
    skipOne(iter, IsTab());

    // read column 8: phase
    record.phase = value(iter);
    if (record.phase != '0' && record.phase != '1' && record.phase != '2')
    {
        record.phase = '.';
    }
    skipOne(iter);
    skipOne(iter, IsTab());

    // It's fine if there are no attributes and the line ends here.
    if (atEnd(iter) || isNewline(value(iter)))
    {
        skipLine(iter);
        return;
    }

    // read column 9: attributes
    while (!atEnd(iter))
    {

        String<char> _key;
        String<char> _value;
        // Read next key/value pair.
        _parseReadGffKeyValue(_value, _key, iter);

        appendValue(record.tagName, _key);
        appendValue(record.tagValue, _value);

        clear(_key);
        clear(_value);

        // At end of line:  Skip EOL and break.
        if (!atEnd(iter) && isNewline(value(iter)))
        {
            skipOne(iter);
            break;
        }
    }
    return;
}

template <typename TFwdIterator, typename TTag>
inline void 
readRecord(GffRecord & record, TFwdIterator & iter, Tag<TTag> const & /*tag*/, GffContext & context)
{
    _readGffRecord(record, iter, context);
}

template <typename TFwdIterator, typename TTag>
inline void 
readRecord(GffRecord & record, TFwdIterator & iter, Tag<TTag> const & tag)
{
    GffContext context;
    readRecord(record, iter, tag, context);
}

// TODO(singer): Needs proper documentation!!! Check the length of the stores!!!
template <typename TFwdIterator, typename TContextSpec, typename TContextSpec2>
inline void
_readGffRecord(GffRecord & record, TFwdIterator & iter, GffIOContext<TContextSpec, TContextSpec2> & ioContext, GffContext & gffContext)
{
    // Read record with string ref from GFF file.
    _readGffRecord(record, iter, gffContext);

    // Translate ref to rID using the context.  If there is no such sequence name in the context yet then we add it.
    unsigned idx = 0;
    if (!getIdByName(nameStore(ioContext), record.ref, idx, nameStoreCache(ioContext)))
    {
        idx = length(nameStore(ioContext));
        appendName(nameStore(ioContext), record.ref, nameStoreCache(ioContext));
    }
    record.rID = idx;
}

template <typename TFwdIterator, typename TContextSpec, typename TContextSpec2>
inline void
_readGffRecord(GffRecord & record, TFwdIterator & iter, GffIOContext<TContextSpec, TContextSpec2> & ioCcontext)
{
    GffContext gffContext;
    _readGffRecord(record, iter, ioCcontext, gffContext);
}


template <typename TRecordReader, typename TContextSpec, typename TContextSpec2, typename TTag>
inline void
readRecord(GffRecord & record, TRecordReader & reader, GffIOContext<TContextSpec, TContextSpec2> & context, Tag<TTag> const & /*tag*/)
{
    _readGffRecord(record, reader, context);
}

// ----------------------------------------------------------------------------
// Function _writeSemicolonSensitive()
// ----------------------------------------------------------------------------

// This function checks if the string to be written contains a semicolon. If
// this is the case then quotes are written around the string.
// Returns false on success.

template <typename TTargetStream, typename TString>
inline void
_writeInQuotes(TTargetStream & target, TString & temp)
{
    // TODO(jsinger): What about escaping quote chars '"'?
    writeValue(target, '"');
    write(target, temp);
    writeValue(target, '"');
}

template <typename TTarget, typename TString, typename TMustBeQuotedFunctor>
inline void
_writePossiblyInQuotes(TTarget& target, TString & source, TMustBeQuotedFunctor const &func)
{
    // TODO(jsinger): What about escaping quote chars '"'?
    typedef typename Iterator<TString>::Type TIter;
    TIter itEnd = end(source, Standard());
    for (TIter it = begin(source, Standard()); it != itEnd; ++it)
    {
        // we have a problem if the string contains a '"' or a line break
        if (*it == '\n' || *it == '"')
            throw std::runtime_error("Attribute contains illegal character!");

        if (func(*it))
        {
            _writeInQuotes(target, source);
            return;
        }
    }
    write(target, source);
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

template <typename TTarget>
inline void
_writeAdditionalSeperator(TTarget const & /*target*/, Gff)
{
    return;
}

template <typename TTarget>
inline void
_writeAdditionalSeperator(TTarget & target, Gtf)
{
    writeValue(target, ' ');
    return;
}


template <typename TTarget, typename TTag>
inline void
_writeAttributes(TTarget & target, GffRecord const & record, TTag const & tag)
{
    const char separatorBetweenTagAndValue = (IsSameType<TTag, Gff>::VALUE)? '=' : ' ';
    for (unsigned i = 0; i < length(record.tagName); ++i)
    {
        if (i != 0)
        {
            writeValue(target, ';');

            // In GTF files a space follows the semicolon
            _writeAdditionalSeperator(target, tag);
       }

        _writePossiblyInQuotes(target, record.tagName[i], GffRecordKeyMustBeQuoted_<TTag>());

        if (!empty(record.tagValue[i]))
        {
            writeValue(target, separatorBetweenTagAndValue);
            _writePossiblyInQuotes(target, record.tagValue[i], GffRecordValueMustBeQuoted_<TTag>());
        }
    }

    // In GTF files each (especially the last) attribute must end with a semi-colon
    if (IsSameType<TTag, Gtf>::VALUE && !empty(record.tagName))
        writeValue(target, ';');

    return;
}

//TODO(singer): No check whether the record is complete!
template <typename TTarget, typename TSeqId, typename TTag>
inline void
_writeRecordImpl(TTarget & target, GffRecord const & record, TSeqId const & ref, TTag tag)
{
    // ignore empty annotations, i.e. annotations that are 'guessed' by implicit information from their children (in GFF)
    if (empty(ref))
        return;

    // write column 1: seqid
    //typename Iterator<TSeqId const, Rooted>::Type itRef = begin(ref);
    write(target, ref);
    writeValue(target, '\t');

    // write column 2: source
    if (empty(record.source))
    {
        writeValue(target, '.');
    }
    else
    {
        write(target, record.source);
    }
    writeValue(target, '\t');

    // write column 3: type
    write(target, record.type);
    writeValue(target, '\t');

    // write column 4: begin position
    if (record.beginPos != (unsigned)-1)
        appendNumber(target, record.beginPos + 1);
    else
        throw std::runtime_error("No start position!");
    writeValue(target, '\t');

    // write column 5: end position
    if (record.endPos != (unsigned)-1 && record.beginPos <= record.endPos)
        appendNumber(target, record.endPos);
    else
        throw std::runtime_error("No end position!");
    writeValue(target, '\t');

    // write column 6: score
    if (record.score != record.score)
        writeValue(target, '.');
    else
        writeValue(target, record.score);
    writeValue(target, '\t');

    // write column 7: strand
    writeValue(target, record.strand);
    writeValue(target, '\t');

    // write column 8: phase
    writeValue(target, record.phase);
    writeValue(target, '\t');

    // write column 9: attributes
    // only until length - 1, because there is no semicolon at the end of the line

    _writeAttributes(target, record, tag);

    writeValue(target, '\n');
    return;
}

template <typename TTarget, typename TTag>
inline void
writeRecord(TTarget & target, GffRecord const & record, TTag const tag)
{
    _writeRecordImpl(target, record, record.ref, tag);
}

template <typename TTarget, typename TContextSpec, typename TContextSpec2, typename TTag>
inline void
writeRecord(TTarget & target, GffRecord const & record, GffIOContext<TContextSpec, TContextSpec2> & context, TTag const tag)
{
    if (record.rID != GffRecord::INVALID_IDX)
    {
        String<char> tempSeqId = nameStore(context)[record.rID];
        return _writeRecordImpl(target, record, tempSeqId, tag);
    }
    return _writeRecordImpl(target, record, record.ref, tag);
}

}  // namespace seqan

#endif  // CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_

