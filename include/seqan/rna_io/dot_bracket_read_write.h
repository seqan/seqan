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


/* IMPLEMENTATION NOTES

DotBracket FORMAT example:

>S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)

*/
namespace seqan{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// ============================================================================
// Forwards
// ============================================================================
// --------------------------------------------------------------------------
// Tag RNA
// --------------------------------------------------------------------------

struct DotBracket_;
typedef Tag<DotBracket_> DotBracket;

template <typename T>
struct MagicHeader<DotBracket, T>
{
    static char const VALUE[1];
};
template <typename T>
char const MagicHeader<DotBracket, T>::VALUE[1] = { '>' };  // DotBracket's first character


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
// ----------------------------------------------------------------------------
// Function readHeader(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readHeader(RNAHeader & header, RNAIOContext & context, TForwardIter & iter, DotBracket const & /*tag*/)
{
    //>S.cerevisiae_tRNA-PHE M10740/1-73
    skipOne(iter);
    readUntil(header.name, iter, IsWhitespace());
    skipUntil(iter, EqualsChar<'/'>());
    skipOne(iter);
    readUntil(context.buffer, iter, EqualsChar<'-'>());
    if (!lexicalCast(header.begPos, context.buffer))
        throw BadLexicalCast(header.begPos, context.buffer);
    clear(context.buffer);
    skipOne(iter);
    readUntil(context.buffer, iter, IsNewline());
    if (!lexicalCast(header.endPos, context.buffer))
        throw BadLexicalCast(header.endPos, context.buffer);
    clear(context.buffer);

    
}


// ----------------------------------------------------------------------------
// Function readRecord(); RNARecord, DotBracket
// ----------------------------------------------------------------------------
template <typename TForwardIter>
inline void 
readRecord(RNARecord & record, RNAIOContext & context, TForwardIter & iter, DotBracket const & /*tag*/)
{
    //GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
    //(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
    CharString brackets = "()<>[]aaAA";
    readUntil(record.base, iter, IsNewline());
    for(unsigned i = 0; i < length(record.base); ++i)
    {
        readOne(context.buffer, iter);
       //s if(context.buffer == )
            
    }


}
// ----------------------------------------------------------------------------
// Function writeHeader(); RNAHeader, DotBracket 
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeHeader(TTarget & target, RNAHeader const & header, DotBracket const & /*tag*/)
{   
    write(target, header.name);
    writeValue(target, ' ');
    writeValue(target, '/');
    write(target, "1-");
    write(target, header.amount);
    writeValue(target, '\n');

}


// ----------------------------------------------------------------------------
// Function writeRecord(); RNARecord, DotBracket 
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RNARecord const & record, DotBracket const & /*tag*/)     
{
    write(target, record.base);
  
    writeValue(target, '\n');
    for(unsigned i =0; i < length(record.pair); ++i)
    {
        write(target, record.pair[i]);
    }
    writeValue(target, '\n');
}


}
#endif // SEQAN_DOTBRACKET_FORMAT_IO_H