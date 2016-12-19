// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Authors: Lily Shellhammer <lily.shellhammer@gmail.com>
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// This file contains routines to read and write to connect format files (.ct)
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_READ_WRITE_H_

/* IMPLEMENTATION NOTES

Rna FORMAT example:

=> record START : number of bases in the sequence
=> record END: title of the structure
=> Each line has information about a base pair in the sequence
    Each line is a base, with this order of information:
        - Base number: index n
        - Base (A, C, G, T, U, X)
        - Index n-1
        - Index n+1
        - Number of the base to which n is paired. No pairing is indicated by 0 (zero).
        - Natural numbering. RnaStructure ignores the actual value given in natural numbering,
            so it is easiest to repeat n here.

CT Files can hold multiple structures of a single sequence.
This is done by repeating the format for each structure without any blank lines between structures.

record
 N  SEQUENCE    N-1     N+1     J POSITION      N
 1  G           0       2       72              1
 2  C           1       3       71              2
 3  G           2       4       70              3

*/

namespace seqan{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag Connect
// --------------------------------------------------------------------------

/*!
 * @tag FileFormats#Connect
 * @headerfile <seqan/rna_io.h>
 * @brief Connect format for RNA structures (*.ct).
 * @signature typedef Tag<Connect_> Connect;
 * @see FileFormats#RnaStruct
 */
struct Connect_;
typedef Tag<Connect_> Connect;

// --------------------------------------------------------------------------
// Class MagicHeader
// --------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Connect, T> :
    public MagicHeader<Nothing, T> {};

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Connect, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileExtensions<Connect, T>::VALUE[1] =
{
    ".ct"      // default output extension
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function readRecord(); RnaRecord, Connect
// --------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, TForwardIter & iter, Connect const & /*tag*/)
{
    RnaStructureGraph graph;
    std::string buffer;
    clear(record);

    // read number of entries (sequence length)
    skipUntil(iter, NotFunctor<IsWhitespace>());
    readUntil(buffer, iter, IsWhitespace());
    if (!lexicalCast(record.seqLen, buffer))
        SEQAN_THROW(BadLexicalCast(record.seqLen, buffer));

    clear(buffer);

    // read energy
    skipUntil(iter, NotFunctor<IsWhitespace>());
    readUntil(buffer, iter, OrFunctor<IsWhitespace, EqualsChar<'='> >());
    if (startsWith(buffer, "ENERGY"))
    {
        skipUntil(iter, EqualsChar<'='>());
        skipOne(iter);
        skipUntil(iter, NotFunctor<IsWhitespace>());
        clear(buffer);
        readUntil(buffer, iter, IsWhitespace());
        if (!lexicalCast(graph.energy, buffer))
            SEQAN_THROW(BadLexicalCast(graph.energy, buffer));

        skipUntil(iter, NotFunctor<IsWhitespace>());
    }
    else
    {
        record.name = buffer;
    }
    clear(buffer);

    // read name
    readUntil(buffer, iter, IsNewline());
    append(record.name, buffer);
    clear(buffer);

    /*
    Example records:

     3  G           2       4       70          3
     N  SEQUENCE   N-1     N+1    J POSITION  N
    */

    // read nucleotides with pairs
    unsigned currPos{};
    record.offset = 1;
    while (!atEnd(iter) && currPos < record.seqLen + record.offset - 1)
    {
        // offset
        skipUntil(iter, NotFunctor<IsWhitespace>());
        if (currPos > 0)
        {
            skipUntil(iter, IsWhitespace());
            ++currPos;
        }
        else
        {
            readUntil(buffer, iter, IsWhitespace());
            if (!lexicalCast(record.offset, buffer))
                SEQAN_THROW(BadLexicalCast(record.offset, buffer));

            currPos = record.offset;
            clear(buffer);
        }
        skipUntil(iter, NotFunctor<IsWhitespace>());

        // nucleotide
        readUntil(buffer, iter, IsWhitespace());
        append(record.sequence, buffer);
        clear(buffer);
        addVertex(graph.inter);

        // skip redundant indices
        skipUntil(iter, NotFunctor<IsWhitespace>());
        skipUntil(iter, IsWhitespace());
        skipUntil(iter, NotFunctor<IsWhitespace>());
        skipUntil(iter, IsWhitespace());
        skipUntil(iter, NotFunctor<IsWhitespace>());

        // paired position: add undirected edge (weight=1.0) if connected
        unsigned pairPos;
        readUntil(buffer, iter, IsWhitespace());
        if (!lexicalCast(pairPos, buffer))
            SEQAN_THROW(BadLexicalCast(pairPos, buffer));

        if (pairPos != 0 && currPos > pairPos)
        {
            if (pairPos >= record.offset)
                addEdge(graph.inter, pairPos - record.offset, currPos - record.offset, 1.0);
            else
                SEQAN_THROW(ParseError("ERROR: Incompatible pairing position in input file."));
        }

        clear(buffer);
        skipLine(iter);
    }
    append(record.fixedGraphs, graph);
    SEQAN_ASSERT_EQ(record.seqLen, length(record.sequence));
}

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, RnaIOContext & /*context*/, TForwardIter & iter, Connect const & /*tag*/)
{
    readRecord(record, iter, Connect());
}

// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, Connect
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, Connect const & /*tag*/)
{
    if (empty(record.sequence) && length(rows(record.align)) != 1)
        SEQAN_THROW(ParseError("ERROR: Connect formatted file cannot contain an alignment."));
    if (length(record.fixedGraphs) != 1)
        SEQAN_THROW(ParseError("ERROR: Connect formatted file cannot contain multiple structure graphs."));

    Rna5String const sequence = empty(record.sequence) ? source(row(record.align, 0)) : record.sequence;
    RnaStructureGraph const & graph = record.fixedGraphs[0];

    // write "header"
    appendNumber(target, record.seqLen);
    if (graph.energy != 0.0f)
    {
        writeValue(target, '\t');
        write(target, "ENERGY = ");
        appendNumber(target, graph.energy);
    }
    writeValue(target, '\t');
    write(target, record.name);
    writeValue(target, '\n');

    // write "body"
    unsigned offset = record.offset > 0 ? record.offset : 1;
    for (TSizeRna5String idx = 0; idx < length(sequence); ++idx)
    {
        writeValue(target, ' ');  // All records start with a space
        appendNumber(target, idx + offset);
        writeValue(target, '\t');
        write(target, sequence[idx]);
        writeValue(target, '\t');
        appendNumber(target, idx + offset - 1);
        writeValue(target, '\t');
        appendNumber(target, idx + offset + 1);
        writeValue(target, '\t');
        if (degree(graph.inter, idx) != 0)
        {
            RnaAdjacencyIterator adjIter(graph.inter, idx);
            appendNumber(target, value(adjIter) + offset);
        }
        else
        {
            writeValue(target, '0');
        }
        writeValue(target, '\t');
        appendNumber(target, idx + offset);
        writeValue(target, '\n');
    }
}

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, RnaIOContext & /*context*/, Connect const & /*tag*/)
{
    writeRecord(target, record, Connect());
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_READ_WRITE_H_
