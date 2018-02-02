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
// Authors: Lily Shellhammer <lily.shellhammer@gmail.com>
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// This file contains routines to write to DotBracket format files (.ct)
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_DOT_BRACKET_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_DOT_BRACKET_READ_WRITE_H_

#include <seqan/stream.h>
#include <stack>
#include <array>
#include <cstddef>

/* IMPLEMENTATION NOTES

DotBracket FORMAT example:

>S.cerevisiae_tRna-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)

-Header with amount of records
-Strand
-Dot-Bracket notation for knots or pseudoknots | Energy of rna strand
*/

namespace seqan{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag DotBracket
// --------------------------------------------------------------------------

/*!
 * @tag FileFormats#DotBracket
 * @headerfile <seqan/rna_io.h>
 * @brief Dot Bracket format for RNA structures (*.dbn).
 * @signature typedef Tag<DotBracket_> DotBracket;
 * @see FileFormats#RnaStruct
 */
struct DotBracket_;
typedef Tag<DotBracket_> DotBracket;

// --------------------------------------------------------------------------
// Class Magicheader
// --------------------------------------------------------------------------

template <typename T>
struct MagicHeader<DotBracket, T> :
    public MagicHeader<Nothing, T> {};

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<DotBracket, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileExtensions<DotBracket, T>::VALUE[1] =
{
    ".dbn"      // default output extension
};

// --------------------------------------------------------------------------
// Metafunction DotBracketArgs
// --------------------------------------------------------------------------

template <typename T = void>
struct DotBracketArgs
{
    constexpr static std::array<char, 30> const open = {{
        '(', '{', '<', '[', 'A', 'B', 'C', 'D', 'E', 'F',
        'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
        'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'}};
    constexpr static std::array<char, 30> const close = {{
        ')', '}', '>', ']', 'a', 'b', 'c', 'd', 'e', 'f',
        'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
        'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'}};
    constexpr static std::array<char, 6> const unpaired = {{
        '.', ',', ':', '_', '-', '='}};
};

template <typename T>
constexpr std::array<char, 30> const DotBracketArgs<T>::open;

template <typename T>
constexpr std::array<char, 30> const DotBracketArgs<T>::close;

template <typename T>
constexpr std::array<char, 6> const DotBracketArgs<T>::unpaired;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Function for converting a bracket string to an undirected graph
// ----------------------------------------------------------------------------

inline void bracket2graph(String<RnaStructureGraph> & graphSet, CharString const & bracketStr)
{
    typedef typename Size<CharString>::Type TCharStringSize;
    RnaStructureGraph graph;
    for (TCharStringSize idx = 0; idx < length(bracketStr); ++idx)
        addVertex(graph.inter);

    // declare stacks for different bracket pairs
    std::stack<TCharStringSize> stack[length(DotBracketArgs<>::open)];

    for (TCharStringSize idx = 0; idx < length(bracketStr); ++idx)
    {
        // skip unpaired
        auto unp_idx = std::find(DotBracketArgs<>::unpaired.begin(), DotBracketArgs<>::unpaired.end(), bracketStr[idx]);
        if (unp_idx != std::end(DotBracketArgs<>::unpaired))
            continue;

        // search opening bracket
        auto open_idx = std::find(DotBracketArgs<>::open.begin(), DotBracketArgs<>::open.end(), bracketStr[idx]);
        if (open_idx != std::end(DotBracketArgs<>::open))
        {
            std::ptrdiff_t brIndex = open_idx - DotBracketArgs<>::open.begin();
            stack[brIndex].push(idx);
            continue;
        }

        // search closing bracket
        auto clos_idx = std::find(DotBracketArgs<>::close.begin(), DotBracketArgs<>::close.end(), bracketStr[idx]);
        if (clos_idx != std::end(DotBracketArgs<>::close))
        {
            std::ptrdiff_t brIndex = clos_idx - DotBracketArgs<>::close.begin();
            if (!stack[brIndex].empty())
            {
                addEdge(graph.inter, idx, stack[brIndex].top(), 1.0);
                stack[brIndex].pop();
            }
            else
            {
                SEQAN_THROW(ParseError("Invalid bracket notation: unpaired closing bracket"));
            }
        }
        else
        {
            SEQAN_THROW(ParseError("Invalid bracket notation: unknown symbol"));
        }
    }

    for (decltype(length(DotBracketArgs<>::open)) idx = 0; idx < length(DotBracketArgs<>::open); ++idx)
    {
        if (!stack[idx].empty())
            SEQAN_THROW(ParseError("Invalid bracket notation: unpaired opening bracket"));
    }

    append(graphSet, graph);
}

// ----------------------------------------------------------------------------
// Helper Function for converting a bracket string to an undirected graph
// ----------------------------------------------------------------------------

inline CharString const graph2bracket(RnaStructureGraph const & graph)
{
    CharString bracketStr{};
    std::stack<unsigned> endpos_stack;                  // stack stores endpos of outer bracket
    String<unsigned> colors;                            // colors for bracket pairs
    resize(colors, numVertices(graph.inter), 0);        // color 0 means 'not set'
    resize(bracketStr, numVertices(graph.inter), ' ');

    unsigned unprocessed = 1; /* any value > 0 */
    for (unsigned col = 1; unprocessed > 0; ++col)
    {
        unsigned bracket_end = 0;                       // end position of current bracket
        unprocessed = 0;

        for (unsigned idx = 0; idx < numVertices(graph.inter); ++idx)
        {
            if (degree(graph.inter, idx) == 0 || (colors[idx] > 0 && colors[idx] < col))
                continue;                               // skip processed and unpaired entries

            RnaAdjacencyIterator adj_it(graph.inter, idx);
            unsigned const p_end = value(adj_it);       // paired end bracket

            if (p_end < bracket_end && idx < p_end)     // open bracket inside previous bracket
            {
                endpos_stack.push(bracket_end);
                bracket_end = p_end;
                colors[idx] = colors[p_end] = col;
            }
            else if (idx >= bracket_end)                // bracket behind previous bracket
            {
                if (endpos_stack.empty())
                {
                    if (idx < p_end)                    // open bracket on base level
                    {
                        bracket_end = p_end;
                        colors[idx] = colors[p_end] = col;
                    }
                    else                                // close bracket on base level
                    {
                        bracket_end = idx;
                    }
                }
                else                                    // close bracket, recover endpos from stack
                {
                    bracket_end = endpos_stack.top();
                    endpos_stack.pop();
                }
            }
            else                                        // bracket will get different color
            {
                ++unprocessed;
            }
        }

        while (!endpos_stack.empty())                   // reset stack for next color
        {
            endpos_stack.pop();
        }
    }

    // write pairs in bracket notation
    for (typename Size<String<unsigned> >::Type idx = 0; idx < length(colors); ++idx)
    {
        if (degree(graph.inter, idx) == 0)              // unpaired
        {
            SEQAN_ASSERT(colors[idx] == 0);
            bracketStr[idx] = '.';
            continue;
        }

        RnaAdjacencyIterator adj_it(graph.inter, idx);
        if (idx < value(adj_it))                        // open bracket
        {
            SEQAN_ASSERT(colors[idx] > 0);
            bracketStr[idx] = DotBracketArgs<>::open[colors[idx]-1];
        }
        else                                            // close bracket
        {
            SEQAN_ASSERT(colors[idx] > 0);
            bracketStr[idx] = DotBracketArgs<>::close[colors[idx]-1];
        }
    }
    return bracketStr;
}

// ----------------------------------------------------------------------------
// Function readRecord(); RnaRecord, DotBracket
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, TForwardIter & iter, DotBracket const & /*tag*/)
{
    std::string buffer;
    clear(record);

    // read name (and offset)
    skipOne(iter);                                                      // ">" symbol
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
    {
        SEQAN_THROW(ParseError("ERROR: Bracket string (" + std::to_string(length(buffer)) +
                               ") must be as long as sequence (" + std::to_string(record.seqLen) + ")."));
    }

    bracket2graph(record.fixedGraphs, buffer);
    clear(buffer);

    // read energy if present
    skipUntil(iter, OrFunctor<IsNewline, EqualsChar<'('> >());
    if (!atEnd(iter) && value(iter) == '(')
    {
        skipOne(iter);
        readUntil(buffer, iter, EqualsChar<')'>());
        if (!lexicalCast(record.fixedGraphs[0].energy, buffer))
            SEQAN_THROW(BadLexicalCast(record.fixedGraphs[0].energy, buffer));
        clear(buffer);
    }
    if (!atEnd(iter))
        skipLine(iter);
}

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, RnaIOContext & /*context*/, TForwardIter & iter, DotBracket const & /*tag*/)
{
    readRecord(record, iter, DotBracket());
}

// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, DotBracket
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, DotBracket const & /*tag*/)
{
    if (empty(record.sequence) && length(rows(record.align)) != 1u)
        SEQAN_THROW(ParseError("ERROR: DotBracket formatted file cannot contain an alignment."));
    if (length(record.fixedGraphs) != 1u)
        SEQAN_THROW(ParseError("ERROR: DotBracket formatted file cannot contain multiple structure graphs."));

    Rna5String const sequence = empty(record.sequence) ? source(row(record.align, 0)) : record.sequence;

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
    RnaStructureGraph const & graph = record.fixedGraphs[0];
    write(target, graph2bracket(graph));

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
writeRecord(TTarget & target, RnaRecord const & record, RnaIOContext & /*context*/,
            DotBracket const & /*tag*/)
{
    writeRecord(target, record, DotBracket());
}

}
#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_DOT_BRACKET_READ_WRITE_H_
