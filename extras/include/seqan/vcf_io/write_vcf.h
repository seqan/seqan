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

#ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_WRITE_VCF_H_
#define SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_WRITE_VCF_H_

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

// ----------------------------------------------------------------------------
// Function write()                                                 [VcfHeader]
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should become writeRecord?

/*!
 * @fn VcfIO#write
 * @headerfile <seqan/vcf_io.h>
 * @brief Write a VcfHeader.
 *
 * @signature int write(target, header, context, Vcf());
 *
 * @param[in,out] target  The StreamConcept to write to.
 * @param[out]    header  The VcfHeader to write.
 * @param[in,out] context VcfIOContext to use.
 *
 * @return int A status code, 0 on success, a different value otherwise.
 */

/**
.Function.VCF I/O#write
..cat:VCF I/O
..summary:Write a @Class.VcfHeader@.
..signature:int write(stream, header, context, Vcf())
..param.stream:The @Concept.StreamConcept@ to write to.
...type:Concept.StreamConcept
..param.header:The @Class.VcfHeader@ to write.
...type:Class.VcfHeader
..param.context:The @Class.VcfIOContext@ to use for writing.
...class:Class.VcfIOContext
..return:$0$ on success, $1$ on failure.
..include:seqan/vcf_io.h
*/

template <typename TTarget>
void
write(TTarget & target,
      VcfHeader const & header,
      VcfIOContext const & /*vcfIOContext*/,
      Vcf const & /*tag*/)
{
    for (unsigned i = 0; i < length(header.headerRecords); ++i)
    {
        write(target, "##");
        write(target, header.headerRecords[i].key);
        writeValue(target, '=');
        write(target, header.headerRecords[i].value);
        writeValue(target, '\n');
    }

    write(target, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (unsigned i = 0; i < length(header.sampleNames); ++i)
    {
        writeValue(target, '\t');
        write(target, header.sampleNames[i]);
    }
    writeValue(target, '\n');
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                           [VcfRecord]
// ----------------------------------------------------------------------------


/*!
 * @fn VcfIO#writeRecord
 * @headerfile <seqan/vcf_io.h>
 * @brief Write a VcfRecord.
 *
 * @signature int writeRecord(target, record, context, Vcf());
 *
 * @param[in,out] target  The StreamConcept to write to.
 * @param[out]    record  The VcfRecord to write.
 * @param[in,out] context VcfIOContext to use.
 *
 * @return int A status code, 0 on success, a different value otherwise.
 */

/**
.Function.VCF I/O#writeRecord
..cat:VCF I/O
..summary:Write a @Class.VcfRecord@.
..signature:int writeRecord(stream, record, context, Vcf())
..param.stream:The @Concept.StreamConcept@ to write to.
...type:Concept.StreamConcept
..param.record:The @Class.VcfRecord@ to write.
...type:Class.VcfRecord
..param.context:The @Class.VcfIOContext@ to use for writing.
...class:Class.VcfIOContext
..return:$0$ on success, $1$ on failure.
..include:seqan/vcf_io.h
*/

template <typename TTarget>
void 
writeRecord(TTarget & target,
           VcfRecord const & record,
           VcfIOContext const & vcfIOContext,
           Vcf const & /*tag*/)
{
    write(target, (*vcfIOContext.sequenceNames)[record.rID]);
    writeValue(target, '\t');
    appendNumber(target, record.beginPos + 1);
    writeValue(target, '\t');
    if (empty(record.id))
        writeValue(target, '.');
    else
        write(target, record.id);
    writeValue(target, '\t');
    if (empty(record.ref))
        writeValue(target, '.');
    else
        write(target, record.ref);
    writeValue(target, '\t');
    if (empty(record.alt))
        writeValue(target, '.');
    else
        write(target, record.alt);
    writeValue(target, '\t');
    if (record.qual != record.qual)  // only way to test for nan
        writeValue(target, '.');
    else
        appendNumber(target, record.qual);
    writeValue(target, '\t');
    if (empty(record.filter))
        writeValue(target, '.');
    else
        write(target, record.filter);
    writeValue(target, '\t');
    if (empty(record.info))
        writeValue(target, '.');
    else
        write(target, record.info);
    writeValue(target, '\t');
    if (empty(record.format))
        writeValue(target, '.');
    else
        write(target, record.format);
    for (unsigned i = 0; i < length(record.genotypeInfos); ++i)
    {
        writeValue(target, '\t');
        if (empty(record.genotypeInfos[i]))
            writeValue(target, '.');
        else
            write(target, record.genotypeInfos[i]);
    }
    writeValue(target, '\n');
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INCLUDE_SEQAN_VCF_WRITE_VCF_H_
