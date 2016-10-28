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
// Author: Lily Shellhammer <>
// ==========================================================================
// This file contains routines to write to DotBracket format files (.ct)
// ==========================================================================

#ifndef SEQAN_DOTBRACKET_FORMAT_IO_H
#define SEQAN_DOTBRACKET_FORMAT_IO_H

#include <seqan/stream.h>
#include <stack>
#include <string>

/* IMPLEMENTATION NOTES

DotBracket FORMAT example:

>S.cerevisiae_tRna-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)

*/
namespace seqan{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// ============================================================================
// Classes, Tags
// ============================================================================
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


// ==========================================================================
// Functions
// ==========================================================================

struct group {
    int position;
    char bracket;
};

// ----------------------------------------------------------------------------
// Function readRecord(); RnaRecord, DotBracket
// ----------------------------------------------------------------------------
template <typename TForwardIter>
inline void 
readRecord(RnaRecord & record, RnaIOContext & context, TForwardIter & iter, DotBracket const & /*tag*/)
{
    //read beginning
        //>S.cerevisiae_tRna-PHE M10740/1-73
    //if(iter == '>')       HOW?
    skipOne(iter);

    readUntil(record.name, iter, IsWhitespace());
    //CHECK IF BLANK AFTER NAME
    skipUntil(iter, EqualsChar<'/'>());
    skipOne(iter);
    readUntil(context.buffer, iter, EqualsChar<'-'>());
    if (!lexicalCast(record.begPos, context.buffer))
        throw BadLexicalCast(record.begPos, context.buffer);
    clear(context.buffer);
    skipOne(iter);
    readUntil(context.buffer, iter, IsNewline());
    if (!lexicalCast(record.endPos, context.buffer))
        throw BadLexicalCast(record.endPos, context.buffer);
    clear(context.buffer);
    skipOne(iter);
    //read body
    //GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
    //(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)

    //declare opening and closing brackets. the reason I use two diff strings is becasue I can see if we have a pair by comparing open[i]==close[i]
    readUntil(record.base, iter, IsNewline());
    skipOne(iter);

    resize(record.pair, length(record.base));

    std::string const SEQAN_DOTBRACKET_OPEN  = "({<[abcdefghijklmnopqrstuvwxyz";
    std::string const SEQAN_DOTBRACKET_CLOSE = ")}>]ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    group holder;
    //declare vector as stack
    std::stack<group> v[30];

    readUntil(context.buffer, iter, IsWhitespace());
    for(unsigned i = 0; i < length(record.base); ++i)
    {
        if (context.buffer[i] == '.')
        {
            //take position i and add 0 to pair
            record.pair[i] = 0;
            continue;
        }
        
        std::size_t br_index = SEQAN_DOTBRACKET_OPEN.find(context.buffer[i]);   // search opening bracket
        if (br_index != std::string::npos)
        {
            holder.position = i+1;              //holder gets position
            holder.bracket = context.buffer[i]; //holder gets opening or closing
            v[br_index].push(holder);
            continue;
        }
        
        br_index = SEQAN_DOTBRACKET_CLOSE.find(context.buffer[i]);              // search closing bracket
        if (br_index != std::string::npos && !v[br_index].empty())
        {
            holder = v[br_index].top();
            record.pair[i] = holder.position;
            std::cout << "record.pair[" << i << "] = " << holder.position-1 << std::endl;
            record.pair[holder.position-1] = i+1;
            //string at i gets holder.position for pair
            //position of holder gets i for pair
            v[br_index].pop();
        }
        else
        {
            throw ParseError("Invalid bracket notation");
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
        if(record.amount != -1)  
            write(target, record.amount);
    }
    else if(record.endPos >= record.begPos)
    {
        write(target, record.begPos);
        writeValue(target, '-');
        write(target, record.endPos);
    }
    writeValue(target, '\n');
    //write base
    write(target, record.base);
    writeValue(target, '\n');
    
    //write pairs
    unsigned const NUM = length(record.pair);
    std::stack<unsigned> endpos_stack;    // store endpos of outer bracket
    String<unsigned> colors;     // contains colors for bracket pairs
    resize(colors, NUM, 0);
    
    unsigned cnt = NUM;                 // count number of unprocessed entries
    for (unsigned col = 1; cnt > 0; ++col)
    {
        int bracket_end = -1;           // end position of current bracket
        cnt = 0;
        
        for (unsigned idx = 0; idx < NUM; ++idx)
        {
            unsigned const p_end = record.pair[idx];
            
            if (colors[idx] != 0 || p_end == 0 || idx > p_end)
            {
                continue;               // skip processed and unpaired entries
            }
            else if (p_end < bracket_end)
            {   // bracket inside another bracket
                endpos_stack.push(bracket_end);
                bracket_end = p_end;
                colors[idx] = colors[p_end-1] = col;
            }
            else if (idx > bracket_end)
            {     // bracket behind another bracket
                if (endpos_stack.empty())
                {
                    bracket_end = p_end;
                    colors[idx] = colors[p_end-1] = col;
                }
                else
                {
                    bracket_end = endpos_stack.top();
                    endpos_stack.pop();
                }
            }
            else
            {
                ++cnt;      // idx receives different color
            }
        }
        // reset stack for next color
        while (!endpos_stack.empty())
        {
            endpos_stack.pop();
        }
    }
    
    std::string const SEQAN_DOTBRACKET_OPEN  = "({<[abcdefghijklmnopqrstuvwxyz";
    std::string const SEQAN_DOTBRACKET_CLOSE = ")}>]ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    
    for (unsigned i = 0; i < NUM; ++i)
    {
        if (record.pair[i] == 0)
            writeValue(target, '.');
        else if (i+1 < record.pair[i])
            write(target, SEQAN_DOTBRACKET_OPEN[colors[i]-1]);
        else
            write(target, SEQAN_DOTBRACKET_CLOSE[colors[i]-1]);
    }
    writeValue(target, '\n');
}


}
#endif // SEQAN_DOTBRACKET_FORMAT_IO_H

