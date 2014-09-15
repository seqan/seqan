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
// Author: David Weese <david.weese@fu-berlin.de>
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

inline void
_parseVcfContig(CharString & chromName, CharString const & headerValue)
{
    if (length(headerValue) < 3u)
        return;

    CharString tmp = infix(headerValue, 1, length(headerValue) - 2);
    StringSet<CharString> tmp2;
    strSplit(tmp2, tmp, EqualsChar<','>());
    for (unsigned i = 0; i < length(tmp2); ++i)
    {
        if (!startsWith(tmp2[i], "ID="))
            continue;
        chromName = suffix(tmp2[i], 3);
    }
}

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(VcfHeader & header,
           VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Vcf const & /*tag*/)
{
    clear(header);
    CharString buffer;
    VcfHeaderRecord record;

    while (!atEnd(iter) && value(iter) == '#')
    {
        skipOne(iter);
        clear(buffer);

        if (value(iter) == '#')
        {
            // Is header line.
            skipOne(iter);
            clear(record);

            // Read header key.
            readUntil(record.key, iter, OrFunctor<EqualsChar<'='>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Vcf> >());

            // Skip '='.
            skipOne(iter);

            // Read header value.
            readLine(record.value, iter);
            appendValue(header, record);

            // Parse out name if headerRecord is a contig field.
            if (record.key == "contig")
            {
                _parseVcfContig(buffer, record.value);
                appendName(contigNamesCache(context), buffer);
            }
        }
        else
        {
            // Is line "#CHROM\t...".
            readLine(buffer, iter);
            if (!startsWith(buffer, "CHROM"))
                ParseError("Invalid line with samples.");

            // Split line, get sample names.
            StringSet<CharString> fields;
            strSplit(fields, buffer, IsTab());
            if (length(fields) < 9u)
                ParseError("Not enough fields.");

            // Get sample names.
            for (unsigned i = 9; i < length(fields); ++i)
                appendName(sampleNamesCache(context), fields[i]);
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

// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(VcfRecord & record,
           VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Vcf const & /*tag*/)
{
    typedef OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Vcf> > NextEntry;

    clear(record);
    CharString &buffer = context.buffer;

    // CHROM
    clear(buffer);
    readUntil(buffer, iter, NextEntry());
    if (empty(buffer))
        SEQAN_THROW(EmptyFieldError("CHROM"));
    record.rID = nameToId(contigNamesCache(context), buffer);
    skipOne(iter);

    // POS
    clear(buffer);
    readUntil(buffer, iter, NextEntry());
    if (empty(buffer))
        SEQAN_THROW(EmptyFieldError("POS"));
    record.beginPos = lexicalCast<__int32>(buffer) - 1; // Translate from 1-based to 0-based.
    skipOne(iter);

    // ID
    readUntil(record.id, iter, NextEntry());
    if (empty(record.id))
        SEQAN_THROW(EmptyFieldError("ID"));
    skipOne(iter);

    // REF
    readUntil(record.ref, iter, NextEntry());
    if (empty(record.id))
        SEQAN_THROW(EmptyFieldError("REF"));
    skipOne(iter);

    // ALT
    readUntil(record.alt, iter, NextEntry());
    if (empty(record.id))
        SEQAN_THROW(EmptyFieldError("ALT"));
    skipOne(iter);

    // QUAL
    clear(buffer);
    readUntil(buffer, iter, NextEntry());
    if (empty(buffer))
        SEQAN_THROW(EmptyFieldError("QUAL"));

    if (buffer == ".")
        record.qual = VcfRecord::MISSING_QUAL();
    else
        lexicalCastWithException(record.qual, buffer);

    skipOne(iter);

    // FILTER
    readUntil(record.filter, iter, NextEntry());
    if (empty(record.filter))
        SEQAN_THROW(EmptyFieldError("FILTER"));
    skipOne(iter);

    // INFO
    readUntil(record.info, iter, OrFunctor<IsTab, IsNewline>());
    if (empty(record.info))
        SEQAN_THROW(EmptyFieldError("INFO"));

    // the following columns are optional
    if (atEnd(iter) || IsNewline()(value(iter)))
    {
        skipLine(iter);
        return;
    }
    skipOne(iter);

    // FORMAT
    readUntil(record.format, iter, NextEntry());
    if (empty(record.format))
        SEQAN_THROW(EmptyFieldError("FORMAT"));
    skipOne(iter);

    // The samples.
    unsigned numSamples = length(sampleNames(context));
    for (unsigned i = 0; i < numSamples; ++i)
    {
        clear(buffer);
        if (i + 1 != numSamples)
        {
            readUntil(buffer, iter, NextEntry());
            skipOne(iter);
        }
        else
        {
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
        }

        if (empty(buffer))
        {
            char buffer[30];    // == 9 (GENOTYPE_) + 20 (#digits in MIN_INT64) + 1 (trailing zero)
            sprintf(buffer, "GENOTYPE_%u", i + 1);
            SEQAN_THROW(EmptyFieldError(buffer));
        }
        appendValue(record.genotypeInfos, buffer);
    }

    // skip line break and optional additional columns
    skipLine(iter);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_READ_VCF_H_
