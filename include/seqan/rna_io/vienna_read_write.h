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
// Author: Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// This file contains routines to read/write Vienna format files
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_VIENNA_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_VIENNA_READ_WRITE_H_

#include <seqan/stream.h>
#include <seqan/rna_io/dot_bracket_read_write.h>  // for bracket-graph transformation
#include <stack>
#include <array>

namespace seqan{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag Vienna
// --------------------------------------------------------------------------

/*!
 * @tag FileFormats#Vienna
 * @headerfile <seqan/rna_io.h>
 * @brief Vienna format for RNA structures without pseudoknots (*.dbv).
 * @signature typedef Tag<Vienna_> Vienna;
 * @see FileFormats#RnaStruct
 */
struct Vienna_;
typedef Tag<Vienna_> Vienna;

// --------------------------------------------------------------------------
// Class Magicheader
// --------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Vienna, T> :
    public MagicHeader<Nothing, T> {};

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Vienna, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileExtensions<Vienna, T>::VALUE[1] =
{
    ".dbv"      // default output extension
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(); RnaRecord, Vienna
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, TForwardIter & iter, Vienna const & /*tag*/)
{
    std::string buffer;
    clear(record);

    // read name (and offset)
    skipOne(iter);  // ">" symbol
    readLine(buffer, iter);
    std::string::size_type pos = buffer.find_last_of('/');
    if (pos == std::string::npos)
    {
        record.name = buffer;
    }
    else
    {
        record.name = buffer.substr(0, pos);
        std::string posStr{buffer.substr(pos + 1)};
        pos = posStr.find('-');
        if (pos == std::string::npos || !lexicalCast(record.offset, posStr.substr(0, pos)))
        {
            record.name = buffer;
            record.offset = 1;
        }
    }
    clear(buffer);

    // read sequence
    readLine(record.sequence, iter);
    record.seqLen = length(record.sequence);

    // read bracket string and build graph
    readUntil(buffer, iter, IsWhitespace());
    if (length(buffer) != record.seqLen)
        SEQAN_THROW(ParseError("ERROR: Bracket string must be as long as sequence."));

    RnaStructureGraph graph;
    typedef typename Size<std::string>::Type TStdStringSize;
    for (TStdStringSize idx = 0; idx < length(buffer); ++idx)
        addVertex(graph.inter);

    std::stack<TStdStringSize> stack;
    for (TStdStringSize idx = 0; idx < length(buffer); ++idx)
    {
        if (buffer[idx] == '(')
        {
            stack.push(idx);
        }
        else if (buffer[idx] == ')')
        {
            if (!stack.empty())
            {
                addEdge(graph.inter, idx, stack.top(), 1.0);
                stack.pop();
            }
            else
            {
                SEQAN_THROW(ParseError("Invalid bracket notation: unpaired closing bracket"));
            }
        }
    }
    if(!stack.empty())
        SEQAN_THROW(ParseError("Invalid bracket notation: unpaired opening bracket"));

    append(record.fixedGraphs, graph);
    clear(buffer);

    // read energy if present
    skipUntil(iter, OrFunctor<IsNewline, EqualsChar<'('> >());
    if (!atEnd(iter) && value(iter) == '(')
    {
        skipOne(iter);
        readUntil(buffer, iter, EqualsChar<')'>());
        if (!lexicalCast(record.fixedGraphs[0].energy, buffer))
        {
            SEQAN_THROW(BadLexicalCast(record.fixedGraphs[0].energy, buffer));
        }
        clear(buffer);
    }
    if (!atEnd(iter))
        skipLine(iter);
}

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, RnaIOContext & /*context*/, TForwardIter & iter, Vienna const & /*tag*/)
{
    readRecord(record, iter, Vienna());
}

// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, Vienna
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, Vienna const & /*tag*/)
{
    if (empty(record.sequence) && length(rows(record.align)) != 1u)
        SEQAN_THROW(ParseError("ERROR: Vienna formatted file cannot contain an alignment."));
    if (length(record.fixedGraphs) != 1u)
        SEQAN_THROW(ParseError("ERROR: Vienna formatted file cannot contain multiple structure graphs."));

    Rna5String const sequence = empty(record.sequence) ? source(row(record.align, 0)) : record.sequence;
    RnaStructureGraph const & graph = record.fixedGraphs[0];

    // write opening character for new record entry
    writeValue(target, '>');
    // write name
    write(target, record.name);
    // write index beg-end
    writeValue(target, '/');
    appendNumber(target, record.offset);
    writeValue(target, '-');
    appendNumber(target, record.offset + record.seqLen - 1);
    writeValue(target, '\n');

    // write sequence
    write(target, sequence);
    writeValue(target, '\n');

    // write bracket string
    std::string bracketStr;
    resize(bracketStr, numVertices(graph.inter), ' ');
    std::stack<unsigned> stack;

    for (typename Size<std::string>::Type idx = 0; idx < length(bracketStr); ++idx)  // write pairs in bracket notation
    {
        if (degree(graph.inter, idx) == 0)                  // unpaired
        {
            bracketStr[idx] = '.';
            continue;
        }

        RnaAdjacencyIterator adj_it(graph.inter, idx);
        if (idx < value(adj_it))                            // open bracket
        {
            bracketStr[idx] = '(';
            stack.push(value(adj_it));
        }
        else                                                // close bracket
        {
            bracketStr[idx] = ')';
            if (stack.empty())
            {
                SEQAN_FAIL("Cannot reach here.");
            }
            if (stack.top() == idx)
            {
                stack.pop();
            }
            else
            {
                SEQAN_THROW(ParseError("ERROR: Vienna format does not allow pseudoknots."));
            }
        }
    }
    write(target, bracketStr);

    // write energy
    if (graph.energy != 0.0f)
    {
        write(target, " (");
        appendNumber(target, graph.energy);
        writeValue(target, ')');
    }
    writeValue(target, '\n');
}

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, RnaIOContext & /*context*/, Vienna const & /*tag*/)
{
    writeRecord(target, record, Vienna());
}

}  // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_VIENNA_READ_WRITE_H_
