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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_READ_VCF_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_READ_VCF_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Vcf
// ----------------------------------------------------------------------------

struct Vcf_;
typedef Tag<Vcf_> Vcf;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function splitString
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document.

// Split at splitter, split at whitespace/non-whitespace boundary if splitter == '\xff'.

inline
void splitString(StringSet<CharString> & result,
                 CharString const & str,
                 char splitter = '\xff',
                 char quoteChar = '"')
{
    typedef Iterator<CharString const>::Type TIter;

    CharString buffer;

    if (splitter == '\xff')
    {
        // Algorithm for splitting at whitespace/non-whitespace transition.
        bool inWhitespace = false;
        for (TIter it = begin(str, Standard()); it != end(str, Standard()); ++it)
        {
            if (inWhitespace && !isspace(*it))
            {
                appendValue(result, buffer);
                clear(buffer);
            }

            if (isspace(*it))
            {
                inWhitespace = true;
            }
            else
            {
                inWhitespace = false;
                appendValue(buffer, *it);
            }
        }
        if (!empty(buffer))
            appendValue(result, buffer);
    }
    else
    {
        // Algorithm for splitting at splitter character.
        bool inQuote = false;
        bool isEmpty = true;
        for (TIter it = begin(str, Standard()); it != end(str, Standard()); ++it)
        {
            if (*it == quoteChar)
            {
                inQuote = !inQuote;
                continue;
            }

            if (*it == splitter && !inQuote)
            {
                appendValue(result, buffer);
                clear(buffer);
            }
            else
            {
                isEmpty = false;
                appendValue(buffer, *it);
            }
        }
        if (!isEmpty)
            appendValue(result, buffer);
    }
}

// ----------------------------------------------------------------------------
// Function read()                                                  [VcfHeader]
// ----------------------------------------------------------------------------

/**
.Function.VCF I/O#read
..cat:VCF I/O
..summary:Read a @Class.VcfHeader@.
..signature:int read(header, reader, context, Vcf())
..param.header:The @Class.VcfHeader@ to read into.
...type:Class.VcfHeader
..param.reader:The @Spec.Single-Pass RecordReader@ to read from.
...type:Spec.Single-Pass RecordReader
..param.context:The @Class.VcfIOContext@ to use for reading.
...class:Class.VcfIOContext
..return:$0$ on success, $1$ on failure.
..include:seqan/vcf_io.h
*/

inline void _parseVcfContig(CharString & chromName, CharString const & headerValue)
{
    if (length(headerValue) <= 2u)
        return;

    CharString tmp = infix(headerValue, 1, length(headerValue) - 2);
    StringSet<CharString> tmp2;
    splitString(tmp2, tmp, ',');
    for (unsigned i = 0; i < length(tmp2); ++i)
    {
        if (!startsWith(tmp2[i], "ID="))
            continue;
        chromName = suffix(tmp2[i], 3);
    }
}

template <typename TStream>
int read(VcfHeader & header,
         RecordReader<TStream, SinglePass<> > & reader,
         VcfIOContext & context,
         Vcf const & /*tag*/)
{
    clear(header);
    CharString buffer;
    int res = 0;

    while (!atEnd(reader) && value(reader) == '#')
    {
        goNext(reader);
        if (atEnd(reader))
            return EOF_BEFORE_SUCCESS;
        if (value(reader) == '#')
        {
            // Is header line.
            goNext(reader);
            if (atEnd(reader))
                return EOF_BEFORE_SUCCESS;

            // Read header key.
            VcfHeaderRecord record;
            if ((res = readUntilChar(record.key, reader, '=')) != 0)
                return res;  // Error or EOF.

            // Skip '='.
            goNext(reader);
            if (atEnd(reader))
                return EOF_BEFORE_SUCCESS;

            // Read header value.
            if ((res = readLine(record.value, reader)) != 0)
                return 1;  // Error.

            appendValue(header.headerRecords, record);

            // Parse out name if headerRecord is a contig field.
            if (record.key == "contig")
            {
                CharString chromName;
                _parseVcfContig(chromName, record.value);
                appendValue(header.sequenceNames, chromName);
            }
            refresh(context.sequenceNamesCache);
        }
        else
        {
            clear(buffer);

            // Is line "#CHROM\t...".
            res = readLine(buffer, reader);
            if (!startsWith(buffer, "CHROM"))
                return 1;  // Invalid line with samples.

            // Split line, get sample names.
            StringSet<CharString> fields;
            splitString(fields, buffer);
            if (length(fields) < 9u)
                return 1;  // Not enough fields.

            // Get sample names.
            for (unsigned i = 9; i < length(fields); ++i)
                appendValue(header.sampleNames, fields[i]);
            refresh(context.sampleNamesCache);
        }
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [VcfHeader]
// ----------------------------------------------------------------------------

/**
.Function.VCF I/O#readRecord
..cat:VCF I/O
..summary:Read a @Class.VcfRecord@ from a @Spec.Single-Pass RecordReader@.
..signature:int readRecord(record, reader, context, Vcf())
..param.record:The @Class.VcfRecord@ to read into.
...type:Class.VcfRecord
..param.reader:The @Spec.Single-Pass RecordReader@ to read from.
...type:Spec.Single-Pass RecordReader
..param.context:The @Class.VcfIOContext@ to use for reading.
...class:Class.VcfIOContext
..return:$0$ on success, $1$ on failure.
..include:seqan/vcf_io.h
*/

// Read record, updating list of known sequences if new one occurs.

template <typename TStream>
int readRecord(VcfRecord & record,
               RecordReader<TStream, SinglePass<> > & reader,
               VcfIOContext & context,
               Vcf const & /*tag*/)
{
    clear(record);
    CharString buffer;
    int res = 0;

    // CHROM
    res = readUntilWhitespace(buffer, reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.
    if (!getIdByName(*context.sequenceNames, buffer, record.rID, context.sequenceNamesCache))
    {
        record.rID = length(*context.sequenceNames);
        appendName(*context.sequenceNames, buffer, context.sequenceNamesCache);
    }

    res = skipWhitespaces(reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    // POS
    clear(buffer);
    res = readUntilWhitespace(buffer, reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.
    if (!lexicalCast2(record.beginPos, buffer))
        return 1;  // Could not cast number.
    record.beginPos -= 1;  // Translate from 1-based to 0-based.

    res = skipWhitespaces(reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    // ID
    res = readUntilWhitespace(record.id, reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    res = skipWhitespaces(reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    // REF
    res = readUntilWhitespace(record.ref, reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    res = skipWhitespaces(reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    // ALT
    res = readUntilWhitespace(record.alt, reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    res = skipWhitespaces(reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    // QUAL
    clear(buffer);
    res = readUntilWhitespace(buffer, reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    if (buffer == ".")
        record.qual = VcfRecord::MISSING_QUAL();
    else if (!lexicalCast2(record.qual, buffer))
        return 1;  // Could not cast.

    res = skipWhitespaces(reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    // FILTER
    res = readUntilWhitespace(record.filter, reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    res = skipWhitespaces(reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    // INFO
    res = readUntilWhitespace(record.info, reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    res = skipWhitespaces(reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    // FORMAT
    res = readUntilWhitespace(record.format, reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    res = skipWhitespaces(reader);
    if (res != 0)
        return res;  // Could be EOF_BEFORE_SUCCESS.

    // The samples.
    for (unsigned i = 0; i < length(*context.sampleNames); ++i)
    {
        clear(buffer);
        res = readUntilWhitespace(buffer, reader);
        if (res != 0 && res != EOF_BEFORE_SUCCESS)
            return 1;  // Could not read sample information.
        appendValue(record.genotypeInfos, buffer);
        if (atEnd(reader))
        {
            if ((i + 1) != length(*context.sampleNames))
                return 1;  // Not enough fields.
            else
                break;  // Done
        }
        if (value(reader) == '\r' || value(reader) == '\n')
        {
            res = skipLine(reader);
            if (res != 0 && res != EOF_BEFORE_SUCCESS)
                return res;  // Error skipping beyond the end of the line.
        }
        else
        {
            res = skipWhitespaces(reader);
            if (res != 0)
                return res;  // Could be EOF_BEFORE_SUCCESS.
        }
    }

    // Skip empty lines, necessary for getting to EOF if there is an empty line at the ned of the file.
    while (!atEnd(reader) && (value(reader) == '\r' || value(reader) == '\n'))
        if ((res = skipLine(reader)) != 0)
            return res;

    return 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_READ_VCF_H_
