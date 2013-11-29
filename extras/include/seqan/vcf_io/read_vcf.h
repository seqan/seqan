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

// TODO(holtgrew): Move?

/*!
 * @fn splitString
 * @headerfile <seqan/vcf_io.h>
 * @brief Split at splitter or whitespace/non-whitespace boundary.
 * 
 * The string is split at whitespace/non-whitespace boundary if splitter == '\xff'.
 *
 * @signature void splitString(result, str, splitter, quoteChar);
 *
 * @param[out] result    A StringSet of CharString for the result.
 * @param[in]  str       A CharString to split.
 * @param[in]  splitter  A to use as the splitter.  If this is '\xff' then <tt>str</tt> will be split at
 *                       whitespace/non-whitespace boundaries.  Defaults to '\xff'.
 * @param[in]  quoteChar The char to split at, defaults to '"'.
 */

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

/*!
 * @defgroup VcfIO VCF I/O
 * @brief Routines for VCF I/O.
 */

/*!
 * @fn VcfIO#read
 * @headerfile <seqan/vcf_io.h>
 * @brief Read a VcfHeader.
 *
 * @signature int read(header, reader, context, Vcf());
 *
 * @param[out]    header  The VcfHeader to read into.
 * @param[in,out] reader  The SinglePassRecordReader to use for reading.
 * @param[in,out] context VcfIOContext to use.
 *
 * @return int A status code, 0 on success, a different value otherwise.
 */

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

template <typename TForwardIter>
void 
read(VcfHeader & header,
     TForwardIter & iter,
     VcfIOContext & context,
     Vcf const & /*tag*/)
{
    clear(header);
    CharString buffer;
    int res = 0;

    while (!atEnd(iter) && value(iter) == '#')
    {
        skipOne(iter);
        if (value(iter) == '#')
        {
            // Is header line.
            skipOne(iter);

            // Read header key.
            VcfHeaderRecord record;
            readUntil(record.key, iter, OrFunctor<EqualsChar<'='>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Vcf> >());

            // Skip '='.
            skipOne(iter);

            // Read header value.
            readUntil(record.value, iter, IsNewline());
            skipOne(iter); // skip Newlin
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
            readUntil(buffer, iter, IsNewline());
            skipOne(iter); // skip Newlin
            if (!startsWith(buffer, "CHROM"))
                std::runtime_error("Invalid line with samples.");

            // Split line, get sample names.
            StringSet<CharString> fields;
            splitString(fields, buffer);
            if (length(fields) < 9u)
                std::runtime_error("Not enough fields.");

            // Get sample names.
            for (unsigned i = 9; i < length(fields); ++i)
                appendValue(header.sampleNames, fields[i]);
            refresh(context.sampleNamesCache);
        }
    }
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [VcfHeader]
// ----------------------------------------------------------------------------

/*!
 * @fn VcfIO#readRecord
 * @headerfile <seqan/vcf_io.h>
 * @brief Read a VcfRecord.
 *
 * @signature int readRecord(record, reader, context, Vcf());
 *
 * @param[out]    record  The VcfRecord to read into.
 * @param[in,out] reader  The SinglePassRecordReader to use for reading.
 * @param[in,out] context VcfIOContext to use.
 *
 * @return int A status code, 0 on success, a different value otherwise.
 */

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

struct VcfContext
{
    String<char> buffer;
};

// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter>
void
readRecord(VcfRecord & record,
           TForwardIter & iter,
           VcfIOContext & context,
           VcfContext & vcfContext,
           Vcf const & /*tag*/)
{
    typedef OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Vcf> > NextEntry;

    clear(record);

    // CHROM
    readUntil(vcfContext.buffer, iter, NextEntry());
    if (!getIdByName(*context.sequenceNames, vcfContext.buffer, record.rID, context.sequenceNamesCache))
    {
        record.rID = length(*context.sequenceNames);
        appendName(*context.sequenceNames, vcfContext.buffer, context.sequenceNamesCache);
    }
    skipOne(iter);

    // POS
    // TODO(signer): Do we need this clear?
    clear(vcfContext.buffer);
    readUntil(vcfContext.buffer, iter, NextEntry());
    record.beginPos = lexicalCast<__int32>(vcfContext.buffer) - 1; // Translate from 1-based to 0-based.
    skipOne(iter);

    // ID
    readUntil(record.id, iter, NextEntry());
    skipOne(iter);

    // REF
    readUntil(record.ref, iter, NextEntry());
    skipOne(iter);

    // ALT
    readUntil(record.alt, iter, NextEntry());
    skipOne(iter);

    // QUAL
    clear(vcfContext.buffer);
    readUntil(vcfContext.buffer, iter, NextEntry());

    if (vcfContext.buffer == ".")
        record.qual = VcfRecord::MISSING_QUAL();
    else
        record.qual = lexicalCast<float>(vcfContext.buffer);

    skipOne(iter);

    // FILTER
    readUntil(record.filter, iter, NextEntry());
    skipOne(iter);

    // INFO
    readUntil(record.info, iter, IsWhitespace());
    if (value(iter) == '\n')
    {
        skipOne(iter);
        return;
    }
    else
        skipOne(iter);

    // FORMAT
    readUntil(record.format, iter, NextEntry());
    skipOne(iter);

    // The samples.
    for (unsigned i = 0; i < length(*context.sampleNames); ++i)
    {
        clear(vcfContext.buffer);
        readUntil(vcfContext.buffer, iter, IsWhitespace());

        appendValue(record.genotypeInfos, vcfContext.buffer);
        if (atEnd(iter))
        {
            if ((i + 1) != length(*context.sampleNames))
                throw std::runtime_error("Not enough fields");
            else
                break;  // Done
        }
        if (value(iter) == '\r' || value(iter) == '\n')
        {
            skipOne(iter);
            return;
        }
        else
            skipOne(iter);
    }

    // Skip empty lines, necessary for getting to EOF if there is an empty line at the ned of the file.
    while (!atEnd(iter) && (value(iter) == '\r' || value(iter) == '\n'))
        skipOne(iter);
}

template <typename TForwardIter>
void
readRecord(VcfRecord & record,
           TForwardIter & iter,
           VcfIOContext & context,
           Vcf const & tag)
{
    VcfContext vcfContext;
    readRecord(record, iter, context, vcfContext, tag);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_READ_VCF_H_
