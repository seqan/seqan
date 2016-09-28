// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_BPSEQ_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_BPSEQ_READ_WRITE_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Bpseq
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Bpseq
 * @headerfile <seqan/bpseq_io.h>
 * @brief Variant callinf format file.
 *
 * @signature typedef Tag<Bpseq_> Bpseq;
 */
struct Bpseq_;
typedef Tag<Bpseq_> Bpseq;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BpseqHeader]
// ----------------------------------------------------------------------------


template <typename TForwardIter>
inline void
readHeader(RnaHeader & header, RnaIOContext & context, TForwardIter & iter, Bpseq const & /*tag*/)
{
    // bpseq has no header
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BpseqRecord]
// ----------------------------------------------------------------------------
// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, RnaIOContext & context, TForwardIter & iter, Bpseq const & /*tag*/)
{
    typedef OrFunctor<IsSpace, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bpseq> > NextEntry;
    std::string buffer;
    clear(record);
    clear(context);

    unsigned currPos{0};
    TRnaRecordGraph graph;

    while (!atEnd(iter) && value(iter) == '#') // All the information stored in the # lines are saved in a single line
    {
        skipOne(iter);
        readLine(buffer, iter);
        append(record.comment, buffer);
        clear(buffer);
    }

    while (!atEnd(iter))
    {
        if (value(iter) == '#')
        {
            skipLine(iter);
            continue;
        }

        // read index position
        if (currPos == 0)
        {
            readUntil(buffer, iter, NextEntry());
            if (empty(buffer))
                SEQAN_THROW(EmptyFieldError("BEGPOS"));
            if (!lexicalCast(record.offset, buffer))
                throw BadLexicalCast(record.offset, buffer);
            currPos = record.offset;
            clear(buffer);
        }
        else
        {
            skipUntil(iter, NextEntry());
            ++currPos;
        }

        // read nucleotide
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readUntil(buffer, iter, IsWhitespace());
        appendValue(record.sequence, buffer[0]);
        addVertex(graph);                 // add base to graph
        if (empty(record.sequence))
            SEQAN_THROW(EmptyFieldError("SEQUENCE"));
        clear(buffer);

        // read paired index
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readUntil(buffer, iter, IsWhitespace());
        if (empty(buffer))
            SEQAN_THROW(EmptyFieldError("PAIR"));
        unsigned pairPos;
        if (!lexicalCast(pairPos, buffer))
            throw BadLexicalCast(pairPos, buffer);
        if (pairPos != 0 && currPos > pairPos)    // add edge if base is connected
            addEdge(graph, pairPos - record.offset, currPos - record.offset, 1.);
        clear(buffer);

        skipLine(iter);
    }
    append(record.fixedGraphs, RnaInterGraph(graph));
    record.seqLen = currPos;  //set amount of records

    return;
}

// ----------------------------------------------------------------------------
// Function writeHeader()                                           [BpseqHeader]
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeHeader(TTarget & target,
            RnaHeader const & header,
            RnaIOContext & context,
            Bpseq const & /*tag*/)
{
    // bpseq has no header
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                           [BpseqRecord]
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target,
            RnaRecord const & record,
            RnaIOContext & context,
            Bpseq const & /*tag*/)
{
    if (empty(record.sequence) && length(rows(record.align)) != 1)
        throw std::runtime_error("ERROR: Connect formatted file cannot contain an alignment.");
    if (length(record.fixedGraphs) != 1)
        throw std::runtime_error("ERROR: Connect formatted file cannot contain multiple structure graphs.");

    clear(context);
    if (!empty(record.name))
    {
        write(target, "# ");
        write(target, record.name);
        writeValue(target, '\n');
    }
    if (!empty(record.comment))
    {
        write(target, "# ");
        write(target, record.comment);
        writeValue(target, '\n');
    }

    unsigned offset = record.offset > 0 ? record.offset : 1;
    for (unsigned i = 0; i < record.seqLen; ++i)
    {
        write(target, offset + i);
        writeValue(target, ' ');
        write(target, record.sequence[i]);
        writeValue(target, ' ');
        if (degree(record.fixedGraphs[0].inter, i) != 0)
        {
            TRnaAdjacencyIterator adj_it(record.fixedGraphs[0].inter, i);
            write(target, value(adj_it) + offset);
        }
        else
        {
            writeValue(target, '0');
        }
        writeValue(target, '\n');
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_BPSEQ_READ_WRITE_H_
