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

// ============================================================================
// Classes, Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag DotBracket
// ----------------------------------------------------------------------------
struct DotBracket_;
typedef Tag<DotBracket_> DotBracket;

// ----------------------------------------------------------------------------
// Class Magicheader
// ----------------------------------------------------------------------------
template <typename T>
struct MagicHeader<DotBracket, T> :
    public MagicHeader<Nothing, T> {};

// ============================================================================
// Metafunctions
// ============================================================================

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
    ".dt"      // default output extension
};

// ----------------------------------------------------------------------------
// Metafunction DotBracketArgs
// ----------------------------------------------------------------------------
template <typename T = void>
struct DotBracketArgs
{
    constexpr static std::array<char, 30> const OPEN = {{
        '(', '{', '<', '[', 'A', 'B', 'C', 'D', 'E', 'F',
        'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
        'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'}};
    constexpr static std::array<char, 30> const CLOSE = {{
        ')', '}', '>', ']', 'a', 'b', 'c', 'd', 'e', 'f',
        'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
        'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'}};
    constexpr static std::array<char, 6> const UNPAIRED = {{
        '.', ',', ':', '_', '-', '='}};
};

template <typename T>
constexpr std::array<char, 30> const DotBracketArgs<T>::OPEN;

template <typename T>
constexpr std::array<char, 30> const DotBracketArgs<T>::CLOSE;

template <typename T>
constexpr std::array<char, 6> const DotBracketArgs<T>::UNPAIRED;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Function for converting a bracket string to an undirected graph
// ----------------------------------------------------------------------------
static inline void bracket2graph(TRnaRecordGraph & graph, CharString const & bracket_str)
{
    for (unsigned idx = 0; idx < length(bracket_str); ++idx)
        addVertex(graph);

    // declare stacks for different bracket pairs
    std::stack<unsigned> stack[length(DotBracketArgs<>::OPEN)];

    for(unsigned idx = 0; idx < length(bracket_str); ++idx)
    {
        // skip unpaired
        auto elem = std::find(DotBracketArgs<>::UNPAIRED.begin(), DotBracketArgs<>::UNPAIRED.end(), bracket_str[idx]);
        if (elem != std::end(DotBracketArgs<>::UNPAIRED))
            continue;

        // search opening bracket
        elem = std::find(DotBracketArgs<>::OPEN.begin(), DotBracketArgs<>::OPEN.end(), bracket_str[idx]);
        if (elem != std::end(DotBracketArgs<>::OPEN))
        {
            std::size_t br_index = elem - DotBracketArgs<>::OPEN.begin();
            stack[br_index].push(idx);
            continue;
        }

        // search closing bracket
        elem = std::find(DotBracketArgs<>::CLOSE.begin(), DotBracketArgs<>::CLOSE.end(), bracket_str[idx]);
        if (elem != std::end(DotBracketArgs<>::CLOSE))
        {
            std::size_t br_index = elem - DotBracketArgs<>::CLOSE.begin();
            if (!stack[br_index].empty())
            {
                addEdge(graph, idx, stack[br_index].top(), 1.);
                stack[br_index].pop();
            }
            else
            {
                throw ParseError("Invalid bracket notation: unpaired closing bracket");
            }
        }
        else
        {
            throw ParseError("Invalid bracket notation: unknown symbol");
        }
    }

    for(unsigned idx = 0; idx < length(DotBracketArgs<>::OPEN); ++idx)
        if(!stack[idx].empty())
            throw ParseError("Invalid bracket notation: unpaired opening bracket");
}

// ----------------------------------------------------------------------------
// Helper Function for converting a bracket string to an undirected graph
// ----------------------------------------------------------------------------
static inline void graph2bracket(CharString & bracket_str, TRnaRecordGraph const & graph)
{
    std::stack<unsigned> endpos_stack;                  // stack stores endpos of outer bracket
    String<unsigned> colors;                            // colors for bracket pairs
    resize(colors, numVertices(graph), 0);              // color 0 means 'not set'
    clear(bracket_str);
    resize(bracket_str, numVertices(graph), ' ');

    unsigned unprocessed = 1; /* any value > 0 */
    for (unsigned col = 1; unprocessed > 0; ++col)
    {
        unsigned bracket_end = 0;                       // end position of current bracket
        unprocessed = 0;

        for (unsigned idx = 0; idx < numVertices(graph); ++idx)
        {
            if (degree(graph, idx) == 0 || (colors[idx] > 0 && colors[idx] < col))
                continue;                               // skip processed and unpaired entries

            TAdjacencyIterator adj_it(graph, idx);
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

    for (unsigned idx = 0; idx < length(colors); ++idx) // write pairs in bracket notation
    {
        if (degree(graph, idx) == 0)                    // unpaired
        {
            SEQAN_ASSERT(colors[idx] == 0);
            bracket_str[idx] = '.';
            continue;
        }

        TAdjacencyIterator adj_it(graph, idx);
        if (idx < value(adj_it))                        // open bracket
        {
            SEQAN_ASSERT(colors[idx] > 0);
            bracket_str[idx] = DotBracketArgs<>::OPEN[colors[idx]-1];
        }
        else                                            // close bracket
        {
            SEQAN_ASSERT(colors[idx] > 0);
            bracket_str[idx] = DotBracketArgs<>::CLOSE[colors[idx]-1];
        }
    }
}

// ----------------------------------------------------------------------------
// Function readRecord(); RnaRecord, DotBracket
// ----------------------------------------------------------------------------
template <typename TForwardIter>
inline void 
readRecord(RnaRecord & record, RnaIOContext & context, TForwardIter & iter, DotBracket const & /*tag*/)
{
    clear(context);
    clear(record);
    //read beginning
    //>S.cerevisiae_tRna-PHE M10740/1-73
    skipOne(iter);

    readUntil(record.name, iter, IsWhitespace());
    //Check that this line doesn't end after name
    if(*iter != '\n')
    {
        skipUntil(iter, EqualsChar<'/'>());
        skipOne(iter);
        if(*iter != '\n')
        {
            readUntil(context.buffer, iter, EqualsChar<'-'>());
            if (!lexicalCast(record.begPos, context.buffer))
                throw BadLexicalCast(record.begPos, context.buffer);    //Read in beginning index position
            clear(context.buffer);
            skipOne(iter);
            readUntil(context.buffer, iter, IsNewline());
            if (!lexicalCast(record.endPos, context.buffer))
                throw BadLexicalCast(record.endPos, context.buffer);    //Read in ending index position
            clear(context.buffer);
            skipOne(iter);

            if (record.endPos >= record.begPos)
            {
                record.amount = record.endPos - record.begPos + 1;     //Set record.amount
            }
            else
            {
                std::cerr << "ERROR: End position is greater than beginning position";
            }
        }
    }

    /*
    End part of record:
    
    GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
    (((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
    */

    //declare opening and closing brackets. the reason I use two diff strings is becasue I can see if we have a pair by comparing open[i]==close[i]
    Rna5String rec_base;
    readUntil(rec_base, iter, IsNewline());
    appendValue(record.sequence, rec_base);

    skipOne(iter);
    readUntil(context.buffer, iter, IsWhitespace());

    bracket2graph(record.graph, context.buffer);
    clear(context.buffer);
    if(!atEnd(iter))
    {
        skipUntil(iter, NotFunctor<IsWhitespace>());
        if(*iter == '(')
        {
            skipOne(iter);    
            readUntil(context.buffer, iter, EqualsChar<')'>());
            if (!lexicalCast(record.energy, context.buffer))
                throw BadLexicalCast(record.energy, context.buffer);
            clear(context.buffer);
        }
    }
}


// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, DotBracket 
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, DotBracket const & /*tag*/)     
{
    //write opening character for new record entry
    writeValue(target, '>');
    //write name
    write(target, record.name);
    writeValue(target, ' ');
    //write index beg/end
    writeValue(target, '/');
    if(record.begPos == -1)
    {
        write(target, "1-");
        if(record.amount != 0)  
            write(target, record.amount);
        else
        {
            std::cerr << "ERROR: There is no amount provided." << std::endl;
            return;
        }
    }
    else if(record.endPos >= record.begPos)
    {
        write(target, record.begPos);
        writeValue(target, '-');
        write(target, record.endPos);
    }

    writeValue(target, '\n');
    //write base
    write(target, record.sequence[0]);
    writeValue(target, '\n');

    CharString bracket_str("");
    graph2bracket(bracket_str, record.graph);
    write(target, bracket_str);
    writeValue(target, ' ');
    writeValue(target, '(');
    write(target, record.energy);
    writeValue(target, ')');
    writeValue(target, '\n');
}

}
#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_DOT_BRACKET_READ_WRITE_H_
