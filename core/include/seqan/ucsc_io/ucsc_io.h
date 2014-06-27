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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_IO_UCSC_H
#define SEQAN_HEADER_STORE_IO_UCSC_H

/* IOREV
 *
 * _doc_
 *
 *
 * maybe move this to file/ because its a file format
 *
 */


namespace SEQAN_NAMESPACE_MAIN {

template <typename TSpec>
struct Ucsc_;

/*!
 * @tag Ucsc
 * @headerfile <seqan/store.h>
 * @brief UCSC Genome browser annotation file (aka knownGene format).
 *
 * @signature typedef Tag<Ucsc<UcscKnownGene_> > const Ucsc;
 */

/**
.Tag.File Format.tag.Ucsc:
    Ucsc Genome Browser annotation file (a.k.a. knownGene format).
..include:seqan/store.h
*/

struct UcscKnownGene_;
typedef Tag<Ucsc_<UcscKnownGene_> > const Ucsc;

/*!
 * @tag UcscIsoforms
 * @headerfile <seqan/store.h>
 * @brief UCSC Genome browser annotation file (aka knownIsoforms format).
 *
 * @signature typedef Tag<Ucsc<UcscKnownIsoforms_> > const UcscIsoforms;
 */

/**
.Tag.File Format.tag.UcscIsoforms:
    Ucsc Genome Browser isoforms file (a.k.a. knownIsoforms format).
..include:seqan/store.h
*/
struct UcscKnownIsoforms_;
typedef Tag<Ucsc_<UcscKnownIsoforms_> > const UcscIsoforms;

//////////////////////////////////////////////////////////////////////////////
// Read Ucsc
//////////////////////////////////////////////////////////////////////////////

struct UcscContext
{
    String<char> buffer;
};

template <typename TForwardIter>
inline void
readRecord(
    UcscRecord & record,
    TForwardIter & iter,
    UcscContext & ucscContext)
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
    readUntil(record.contigName, iter, nextRecord);
    skipOne(iter);

    // read column 3: orientation
    clear(ucscContext.buffer);
    readUntil(ucscContext.buffer, iter, nextRecord);

    if (length(ucscContext.buffer) > 0)
    {
        if (length(ucscContext.buffer) == 1u && ucscContext.buffer[0] != '+' && ucscContext.buffer[0] != '-')
        {
            record.format = record.KNOWN_ISOFORMS;
            insert(record.transName, 0, "GENE");
            skipUntil(iter, IsNewline());
            if(!atEnd(iter))
                skipOne(iter);

            return;
        }
        record.format = record.KNOWN_GENE;
        char orientation = ucscContext.buffer[0];
        skipOne(iter);
    }

    // read column 4: transcript begin position
    clear(ucscContext.buffer);
    readUntil(ucscContext.buffer, iter, nextRecord);
    record.annotationBeginPos = lexicalCast<__uint64>(ucscContext.buffer);
    skipOne(iter);

    // read column 5: transcript end position
    clear(ucscContext.buffer);
    readUntil(ucscContext.buffer, iter, nextRecord);
    record.annotationEndPos = lexicalCast<__uint64>(ucscContext.buffer);
    skipOne(iter);

    // read column 6: CDS begin position
    clear(ucscContext.buffer);
    readUntil(ucscContext.buffer, iter, nextRecord);
    record.cdsBegin = lexicalCast<__uint64>(ucscContext.buffer);
    skipOne(iter);

    // read column 7: CDS end position
    clear(ucscContext.buffer);
    readUntil(ucscContext.buffer, iter, nextRecord);
    record.cdsEnd = lexicalCast<__uint64>(ucscContext.buffer);
    skipOne(iter);

    // read column 8: exon count
    int exons;
    clear(ucscContext.buffer);
    readUntil(ucscContext.buffer, iter, nextRecord);
    exons = lexicalCast<int>(ucscContext.buffer);
    skipOne(iter);

    // read column 9: exon begin positions
    for (int i = 0; i < exons; ++i)
    {
        clear(ucscContext.buffer);
        readUntil(ucscContext.buffer, iter, OrFunctor<OrFunctor<EqualsChar<';'>, EqualsChar<','> >, AssertFunctor<NotFunctor<IsNewline>, ParseError, Ucsc> >());

        unsigned long long tempBegin;
        tempBegin = lexicalCast<__uint64>(ucscContext.buffer);
        appendValue(record.exonBegin, tempBegin);
        skipOne(iter);
    }
    skipOne(iter, IsTab());

    // read column 10: exon end positions
    for (int i = 0; i < exons; ++i)
    {
        clear(ucscContext.buffer);
        readUntil(ucscContext.buffer, iter, OrFunctor<OrFunctor<EqualsChar<';'>, EqualsChar<','> >, AssertFunctor<NotFunctor<IsNewline>, ParseError, Ucsc> >());

        unsigned long long tempEnd;
        tempEnd =  lexicalCast<__uint64>(ucscContext.buffer);
        appendValue(record.exonEnd, tempEnd);
        skipOne(iter);
    }
    skipOne(iter, IsTab());

    // read column 11: protein name
    readUntil(record.proteinName, iter, IsWhitespace());

    skipOne(iter);

    // skip column 12
    skipLine(iter);

    // adapt positions
    if (orientation == '-')
    {
        __uint64 tmp = record.annotationBeginPos;
        record.annotationBeginPos = record.annotationEndPos;
        record.annotationEndPos = tmp;
        tmp = record.cdsBegin;
        record.cdsBegin = record.cdsEnd;
        record.cdsEnd = tmp;
        for (int i = 0; i < exons; ++i)
        {
            tmp = record.exonBegin[i];
            record.exonBegin[i] = record.exonEnd[i];
            record.exonEnd[i] = tmp;
        }
    }
}

template <typename TForwardIter, typename TSpec>
inline void
read(
    String<UcscRecord, TSpec> & records,
    TForwardIter & iter)
{
    typedef typename Iterator<String<UcscRecord, TSpec> >::Type TIter;

    clear(records);
    UcscContext context;
    unsigned i = 0;
    while (!atEnd(iter))
    {
        resize(records, length(records) + 1);
        try
        {
            readRecord(records[i], iter, context);
            ++i;
        }
        catch(std::runtime_error())
        {
            
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
// Write Ucsc
//////////////////////////////////////////////////////////////////////////////

template <typename TTarget>
inline void
writeRecord(
    TTarget & target,
    UcscRecord & record)
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
    __uint64 transBeginPos, transEndPos;
    __uint64 cdsBeginPos, cdsEndPos;
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
        appendNumber(target, _min(record.exonBegin[i], record.exonEnd[i]));
        writeValue(target, ',');
    }
    writeValue(target, '\t');

    // write column 10: exon end positions
    for (unsigned i = 0; i < length(record.exonBegin); ++i)
    {
        appendNumber(target, _max(record.exonBegin[i], record.exonEnd[i]));
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

template <typename TTarget, typename TSpec>
inline void
write(
    TTarget & target,
    String<UcscRecord, TSpec> & records)
{
    typedef typename Iterator<String<UcscRecord, TSpec>, Rooted>::Type TIter;

    TIter recordsIter = begin(records);

    while (!atEnd(recordsIter))
    {
        try
        {
            writeRecord(target, value(recordsIter));
        }
        catch(std::runtime_error())
        {

        }
        ++recordsIter;
    }
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
