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
// Authors: Gianvito Urgese <gianvito.urgese@polito.it>
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_BPSEQ_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_BPSEQ_READ_WRITE_H_

#include <seqan/rna_io.h>

namespace seqan {

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag Bpseq
// --------------------------------------------------------------------------

/*!
 * @tag FileFormats#Bpseq
 * @headerfile <seqan/rna_io.h>
 * @brief Bpseq format for RNA structures (*.bpseq).
 * @signature typedef Tag<Bpseq_> Bpseq;
 * @see FileFormats#RnaStruct
 */
struct Bpseq_;
typedef Tag<Bpseq_> Bpseq;

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Bpseq, T> :
    public MagicHeader<Nothing, T> {};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Bpseq, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Bpseq, T>::VALUE[1] =
{
    ".bpseq"     // default output extension
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                          [BpseqRecord]
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, TForwardIter & iter, Bpseq const & /*tag*/)
{
    typedef OrFunctor<IsSpace, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bpseq> > NextEntry;
    std::string buffer;
    clear(record);

    skipUntil(iter, NotFunctor<IsWhitespace>());
    while (!atEnd(iter) && value(iter) == '#')
    {  // All the information stored in the # lines are saved in a single line
        skipOne(iter);
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readLine(buffer, iter);
        if (empty(record.name))
        {
            record.name = buffer;
        }
        else
        {
            appendValue(record.comment, ' ');
            append(record.comment, buffer);
        }
        clear(buffer);
        skipUntil(iter, NotFunctor<IsWhitespace>());
    }

    RnaStructureGraph graph;
    unsigned currPos{};
    while (!atEnd(iter) && value(iter) != '#')
    {
        // read index position
        if (currPos == 0)
        {
            readUntil(buffer, iter, NextEntry());
            if (empty(buffer))
                SEQAN_THROW(EmptyFieldError("BEGPOS"));
            if (!lexicalCast(record.offset, buffer))
                SEQAN_THROW(BadLexicalCast(record.offset, buffer));

            currPos = record.offset;
            clear(buffer);
        }
        else
        {
            skipUntil(iter, NextEntry());
            ++currPos;
        }

        // read nucleotide
        skipUntil(iter, NotFunctor<IsBlank>());
        readUntil(buffer, iter, NextEntry());
        appendValue(record.sequence, buffer[0]);
        addVertex(graph.inter);                 // add base to graph
        if (empty(buffer))
            SEQAN_THROW(EmptyFieldError("SEQUENCE"));

        clear(buffer);

        // read paired index
        skipUntil(iter, NotFunctor<IsBlank>());
        readUntil(buffer, iter, IsWhitespace());
        if (empty(buffer))
            SEQAN_THROW(EmptyFieldError("PAIR"));

        unsigned pairPos;
        if (!lexicalCast(pairPos, buffer))
            SEQAN_THROW(BadLexicalCast(pairPos, buffer));

        if (pairPos != 0 && currPos > pairPos)    // add edge if base is connected
            addEdge(graph.inter, pairPos - record.offset, currPos - record.offset, 1.0);

        clear(buffer);

        skipLine(iter);
        skipUntil(iter, NotFunctor<IsWhitespace>());
    }
    append(record.fixedGraphs, graph);
    record.seqLen = currPos;  //set amount of records
}

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, RnaIOContext & /*context*/, TForwardIter & iter, Bpseq const & /*tag*/)
{
    readRecord(record, iter, Bpseq());
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                         [BpseqRecord]
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, Bpseq const & /*tag*/)
{
    if (empty(record.sequence) && length(rows(record.align)) != 1)
        SEQAN_THROW(ParseError("ERROR: Bpseq formatted file cannot contain an alignment."));
    if (length(record.fixedGraphs) != 1)
        SEQAN_THROW(ParseError("ERROR: Bpseq formatted file cannot contain multiple structure graphs."));

    if (!empty(record.name))
    {
        write(target, "# ");
        write(target, record.name);
        writeValue(target, '\n');
    }
    if (!empty(record.comment))
    {
        writeValue(target, '#');
        write(target, record.comment);
        writeValue(target, '\n');
    }

    unsigned offset = record.offset > 0 ? record.offset : 1;
    for (unsigned i = 0; i < record.seqLen; ++i)
    {
        appendNumber(target, offset + i);
        writeValue(target, '\t');
        write(target, record.sequence[i]);
        writeValue(target, '\t');
        if (degree(record.fixedGraphs[0].inter, i) != 0)
        {
            RnaAdjacencyIterator adj_it(record.fixedGraphs[0].inter, i);
            appendNumber(target, value(adj_it) + offset);
        }
        else
        {
            writeValue(target, '0');
        }
        writeValue(target, '\n');
    }
}

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, RnaIOContext & /*context*/, Bpseq const & /*tag*/)
{
    writeRecord(target, record, Bpseq());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_BPSEQ_READ_WRITE_H_
