// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_UCSC_UCSC_IO_H_
#define SEQAN_CORE_INCLUDE_SEQAN_UCSC_UCSC_IO_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Ucsc
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Ucsc_;

/*!
 * @defgroup UcscFileIO
 * @brief Tags for UCSC I/O.
 */

/*!
 * @tag UcscFileIO#Ucsc
 * @headerfile <seqan/ucsc_io.h>
 * @brief UCSC Genome browser annotation file (aka knownGene format).
 *
 * @signature typedef Tag<Ucsc_<UcscKnownGene_> > const Ucsc;
 */

struct UcscKnownGene_;
typedef Tag<Ucsc_<UcscKnownGene_> > Ucsc;

// ----------------------------------------------------------------------------
// Tag Ucsc
// ----------------------------------------------------------------------------

/*!
 * @tag UcscFileIO#UcscIsoforms
 * @headerfile <seqan/ucsc_io.h>
 * @brief UCSC Genome browser annotation file (aka knownGene format).
 *
 * @signature typedef Tag<Ucsc_<UcscKnownGene_> > const Ucsc;
 */

struct UcscIsoforms_;
typedef Tag<Ucsc_<UcscIsoforms_> > UcscIsoforms;

// ----------------------------------------------------------------------------
// Class UcscIOContext
// ----------------------------------------------------------------------------

/*!
 * @class UcscIOContext
 * @headerfile <sean/ucsc_io.h>
 * @brief Context to use for UCSC file I/O.
 *
 * @signature struct UcscIOContext;
 */
struct UcscIOContext
{
    /*!
     * @var CharString UcscIOContext::buffer
     * @brief Buffer used during UCSC I/O.
     */
    CharString buffer;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Ucsc, T> : public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Ucsc, T>
{
    static char const * VALUE[2];  // default is one extension
};

template <typename T>
char const * FileExtensions<Ucsc, T>::VALUE[2] =
{
    ".txt",     // default output extension
    ".tsv"
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecod
// ----------------------------------------------------------------------------

/*!
 * @fn UcscFileIO#readRecord
 * @headerfile <seqan/ucsc_io.h>
 * @brief Low-level reading for @link UcscRecord @endlink.
 *
 * @signature void recordRecord(record, ucscIOContext, iter, tag);
 *
 * @param[out]    record       The @link UcscRecord @endlink to store results in.
 * @param[in,out] uscscContext The @link UcscIOContext @endlink to use during the rading.
 * @param[in,out] iter         The @link ForwadIteratorConcept forward iterator @endlink to read from.
 * @param[in]     tag          Fixed to @link Ucsc @endlink.
 *
 * @throw IOError    in the case of I/O errors
 * @throw ParseError in the case of problems with parsing
 */
template <typename TForwardIter>
void readRecord(UcscRecord & record,
                UcscIOContext & ucscIOContext,
                TForwardIter & iter,
                Ucsc const & /*tag*/)
{
    OrFunctor<IsWhitespace, AssertFunctor<NotFunctor<IsNewline>, ParseError, Ucsc> > nextRecord;

    clear(record);

    // read column 1: transcript name
    // The letters until the first whitespace will be read.
    // Then, we skip until we hit the first tab character.

    readUntil(record.transName, iter, IsWhitespace());
    if (!empty(record.transName) && record.transName[0] == '#')
    {
        skipLine(iter);
        return;
    }
    skipOne(iter, IsTab());

    // read column 2: contig name
    readUntil(record.contigName, iter, IsWhitespace());
    if (empty(record.contigName))
        throw ParseError("unexpected end of record.");

    // read column 3: orientation
    char orientation;
    readOne(orientation, iter);
    if (IsNewline()(orientation))
    {
        record.format = record.KNOWN_ISOFORMS;
        insert(record.transName, 0, "GENE");
        return;
    }
    readOne(orientation, iter, OrFunctor<EqualsChar<'+'>, EqualsChar<'-'> >());
    skipOne(iter, IsTab());

    record.format = record.KNOWN_GENE;

    // read column 4: transcript begin position
    clear(ucscIOContext.buffer);
    readUntil(ucscIOContext.buffer, iter, nextRecord);
    record.annotationBeginPos = lexicalCast<__uint32>(ucscIOContext.buffer);
    skipOne(iter, IsTab());

    // read column 5: transcript end position
    clear(ucscIOContext.buffer);
    readUntil(ucscIOContext.buffer, iter, nextRecord);
    record.annotationEndPos = lexicalCast<__uint32>(ucscIOContext.buffer);
    skipOne(iter, IsTab());

    // read column 6: CDS begin position
    clear(ucscIOContext.buffer);
    readUntil(ucscIOContext.buffer, iter, nextRecord);
    record.cdsBegin = lexicalCast<__uint32>(ucscIOContext.buffer);
    skipOne(iter, IsTab());

    // read column 7: CDS end position
    clear(ucscIOContext.buffer);
    readUntil(ucscIOContext.buffer, iter, nextRecord);
    record.cdsEnd = lexicalCast<__uint32>(ucscIOContext.buffer);
    skipOne(iter, IsTab());

    // read column 8: exon count
    unsigned int exons;
    clear(ucscIOContext.buffer);
    readUntil(ucscIOContext.buffer, iter, nextRecord);
    exons = lexicalCast<unsigned int>(ucscIOContext.buffer);
    skipOne(iter, IsTab());

    // read column 9: exon begin positions
    for (unsigned int i = 0; i < exons; ++i)
    {
        clear(ucscIOContext.buffer);
        readUntil(ucscIOContext.buffer, iter, OrFunctor<OrFunctor<EqualsChar<';'>, EqualsChar<','> >, AssertFunctor<NotFunctor<IsNewline>, ParseError, Ucsc> >());

        unsigned long long tempBegin;
        tempBegin = lexicalCast<__uint32>(ucscIOContext.buffer);
        appendValue(record.exonBegin, tempBegin);
        skipOne(iter);
    }
    skipOne(iter, IsTab());

    // read column 10: exon end positions
    for (unsigned int i = 0; i < exons; ++i)
    {
        clear(ucscIOContext.buffer);
        readUntil(ucscIOContext.buffer, iter, OrFunctor<OrFunctor<EqualsChar<';'>, EqualsChar<','> >, AssertFunctor<NotFunctor<IsNewline>, ParseError, Ucsc> >());

        unsigned long long tempEnd;
        tempEnd =  lexicalCast<__uint32>(ucscIOContext.buffer);
        appendValue(record.exonEnds, tempEnd);
        skipOne(iter);
    }
    skipOne(iter, IsTab());

    // read column 11: protein name
    readUntil(record.proteinName, iter, IsWhitespace());
    skipOne(iter, IsTab());

    // skip column 12
    skipLine(iter);

    // adapt positions
    if (orientation == '-')
    {
        __uint32 tmp = record.annotationBeginPos;
        record.annotationBeginPos = record.annotationEndPos;
        record.annotationEndPos = tmp;
        tmp = record.cdsBegin;
        record.cdsBegin = record.cdsEnd;
        record.cdsEnd = tmp;
        for (unsigned int i = 0; i < exons; ++i)
        {
            tmp = record.exonBegin[i];
            record.exonBegin[i] = record.exonEnds[i];
            record.exonEnds[i] = tmp;
        }
    }
}

// ----------------------------------------------------------------------------
// Function writeRecord
// ----------------------------------------------------------------------------

// TODO(holtgrew): Dave/Enrico: please verify the type of target here

/*!
 * @fn UcscFileIO#writeRecord
 * @headerfile <seqan/ucsc_io.h>
 * @brief Low-level writing of @link UcscRecord @endlink.
 *
 * @signature void writeRecord(target, record, tag);
 *
 * @param[in,out] target @link OutputIteratorConcept Output iterator @endlink or @link ContainerConcept container
 *                       @endlink to write to.
 * @param[in]     record @link UcscRecord @endlink to write.
 * @param[in]     tag    Fixed to @link Ucsc @endlink.
 *
 * @throw IOError in case of I/O problems
 */
template <typename TTarget, typename TInnerTag>
void writeRecord(TTarget & target,
                 UcscRecord const & record,
                 Tag<Ucsc_<TInnerTag> > const & /*tag*/)
{
    unsigned suf = 0;
    if (record.format == record.KNOWN_ISOFORMS && length(record.transName) >= 4 && prefix(record.transName, 4) == "GENE")
        suf = 4;

    // write column 1: transcript name
    // The letters until the first whitespace will be write.
    // Then, we skip until we hit the first tab character.
    write(target, suffix(record.transName, suf));
    writeValue(target, '\t');

    // write column 2: contig name
    write(target, suffix(record.contigName, suf));

    if (record.format == record.KNOWN_ISOFORMS)
    {
        write(target, '\n');
        return;
    }
    writeValue(target, '\t');

    // write column 3: orientation
    __uint32 transBeginPos, transEndPos;
    __uint32 cdsBeginPos, cdsEndPos;
    if (record.annotationBeginPos < record.annotationEndPos)
    {
        writeValue(target, '+');

        transBeginPos = record.annotationBeginPos;
        transEndPos = record.annotationEndPos;
        cdsBeginPos = record.cdsBegin;
        cdsEndPos = record.cdsEnd;
    }
    else
    {
        writeValue(target, '-');

        transEndPos = record.annotationBeginPos;
        transBeginPos = record.annotationEndPos;
        cdsEndPos = record.cdsBegin;
        cdsBeginPos = record.cdsEnd;
    }
    writeValue(target, '\t');

    // write column 4: transcript begin position
    appendNumber(target, transBeginPos);
    writeValue(target, '\t');

    // write column 5: transcript end position
    appendNumber(target, transEndPos);
    writeValue(target, '\t');

    // write column 6: CDS begin position
    appendNumber(target, cdsBeginPos);
    writeValue(target, '\t');

    // write column 7: CDS end position
    appendNumber(target, cdsEndPos);
    writeValue(target, '\t');

    // write column 8: exon count
    appendNumber(target, length(record.exonBegin));
    writeValue(target, '\t');

    // write column 9: exon begin positions
    for (unsigned i = 0; i < length(record.exonBegin); ++i)
    {
        appendNumber(target, _min(record.exonBegin[i], record.exonEnds[i]));
        writeValue(target, ',');
    }
    writeValue(target, '\t');

    // write column 10: exon end positions
    for (unsigned i = 0; i < length(record.exonBegin); ++i)
    {
        appendNumber(target, _max(record.exonBegin[i], record.exonEnds[i]));
        writeValue(target, ',');
    }
    writeValue(target, '\t');

    // write column 11: protein name
    if (length(record.proteinName) > 0u)
    {
        write(target, record.proteinName);
    }

    writeValue(target, '\t');

    // skip column 12
    if (length(record.transName) > 0u)
    {
        write(target, record.transName);
    }

    writeValue(target, '\n');
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_UCSC_UCSC_IO_H_
