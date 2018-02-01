// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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

#ifndef SEQAN_INCLUDE_SEQAN_VCF_READ_VCF_H_
#define SEQAN_INCLUDE_SEQAN_VCF_READ_VCF_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Vcf
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Vcf
 * @headerfile <seqan/vcf_io.h>
 * @brief Variant callinf format file.
 *
 * @signature typedef Tag<Vcf_> Vcf;
 */
struct Vcf_;
typedef Tag<Vcf_> Vcf;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                            [VcfHeader]
// ----------------------------------------------------------------------------

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TString>
inline void
_readVcfContig(VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context, TString const & headerValue)
{
    typedef OrFunctor<EqualsChar<','>, EqualsChar<'>'> >            IsCommaOrGt;
    typedef typename DirectionIterator<TString const, Input>::Type  TIter;

    TIter headerIter = directionIterator(headerValue, Input());
    CharString &buffer = context.buffer;

    skipOne(headerIter, EqualsChar<'<'>());

    // Seek contig ID key.
    while (!atEnd(headerIter))
    {
        clear(buffer);
        readUntil(buffer, headerIter, EqualsChar<'='>());
        if (buffer == "ID") break;
        skipUntil(headerIter, IsCommaOrGt());
        skipOne(headerIter);
    }

    if (atEnd(headerIter))
        SEQAN_THROW(ParseError("Contig ID key not found in header."));

    // Read contig ID value.
    clear(buffer);
    skipOne(headerIter, EqualsChar<'='>());
    readUntil(buffer, headerIter, IsCommaOrGt());
    if (empty(buffer))
        SEQAN_THROW(ParseError("Contig ID value not found in header."));
    appendName(contigNamesCache(context), buffer);
}

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readHeader(VcfHeader & header,
           VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Vcf const & /*tag*/)
{
    clear(header);
    CharString &buffer = context.buffer;
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
                _readVcfContig(context, record.value);
        }
        else
        {
            // Is line "#CHROM\t...".
            readLine(buffer, iter);
            if (!startsWith(buffer, "CHROM"))
                SEQAN_THROW(ParseError("Invalid line with samples."));

            // Split line, get sample names.
            StringSet<CharString> fields;
            strSplit(fields, buffer, IsTab());
            if (length(fields) < 8u)
                SEQAN_THROW(ParseError("Not enough fields."));

            // Get sample names.
            for (unsigned i = 8; i < length(fields); ++i)
            {
                if(i == 8 && fields[i] == "FORMAT")
                    continue;
                appendName(sampleNamesCache(context), fields[i]);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [VcfRecord]
// ----------------------------------------------------------------------------
// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(VcfRecord & record,
           VcfIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Vcf const & /*tag*/)
{
    clear(record);
    CharString &buffer = context.buffer;

    // get the next line on the buffer.
    clear(buffer);

    readLine(buffer, iter);
    // Split line, get field and sample values.
    // The first 8(9) columns are fields and the rest are values for samples
    //"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    StringSet<CharString> field_values;
    strSplit(field_values, buffer, IsTab(), false);

    unsigned numSamples = length(sampleNames(context));

    if (length(field_values) < 8u + numSamples)
        SEQAN_THROW(ParseError("Not enough values in a line."));

    record.rID      = nameToId(contigNamesCache(context), field_values[0]);
    record.beginPos = lexicalCast<int32_t>(field_values[1]) - 1; // Translate from 1-based to 0-based.
    record.id       = field_values[2];
    record.ref      = field_values[3];
    record.alt      = field_values[4];

    if (field_values[5] == ".")
        record.qual = VcfRecord::MISSING_QUAL();
    else
        lexicalCastWithException(record.qual, field_values[5]);

    record.filter   = field_values[6];
    record.info     = field_values[7];

    //check if we have a spare column for FORMAT
    unsigned samplesColStart = 8;
    if (length(field_values) > 8u + numSamples) // we have extara column for FORMAT
    {
        record.format = field_values[8];
        samplesColStart = 9;
    }

    // Get sample name values .
    for (unsigned i = samplesColStart; i < length(field_values); ++i)
    {
        appendValue(record.genotypeInfos, field_values[i]);
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_VCF_READ_VCF_H_
