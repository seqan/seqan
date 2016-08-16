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
// This file contains routines to read and write to DotBracket format files (.ct)
// ==========================================================================

#ifndef SEQAN_DOTBRACKET_FORMAT_IO_H
#define SEQAN_DOTBRACKET_FORMAT_IO_H

#include <seqan/stream.h>


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
// Forwards
// ============================================================================
// --------------------------------------------------------------------------
// Tag Rna
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
// Function readHeader(); RnaHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readHeader(RnaHeader & header, RnaIOContext & context, TForwardIter & iter, DotBracket const & /*tag*/)
{
    //>S.cerevisiae_tRna-PHE M10740/1-73

    typedef OrFunctor<IsTab, IsNewline> TNextEntry;

    clear(header);

    skipOne(iter, EqualsChar<'>'>());   //skip first header character

    readUntil(header.name, iter, TNextEntry());    //read name 

    skipUntil(iter, NotFunctor<IsWhitespace>());    //skip whitespace
    skipUntil(iter, EqualsChar<'/'>());
    skipOne(iter);
    unsigned front, end;
    readUntil(context.buffer, iter, EqualsChar<'-'>());
    if (!lexicalCast(front, context.buffer))
        throw BadLexicalCast(front, context.buffer);
    clear(context.buffer);
    skipOne(iter, EqualsChar<'-'>());
    readUntil(context.buffer, iter, TNextEntry());
    if (!lexicalCast(end, context.buffer))
        throw BadLexicalCast(front, context.buffer);
    clear(context.buffer);
    header.amount = end-front;

    //SKIP LINES TO GET ENERGY AT END
    
}

template <typename TForwardIter>
inline void 
readRecord(RnaRecord & record, RnaIOContext & context, TForwardIter & iter, DotBracket const & /*tag*/)
{

    readUntil(record.base, iter, IsNewline());
    readUntil(record.pair, iter, IsWhitespace());   


}


// ----------------------------------------------------------------------------
// Function writeRecord(DotBracket);
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, DotBracket const & /*tag*/)     
{
    write(target, record.base);
  
    //Resize charString thing to length(record.base), fill in positions with the positions declared in __
    writeValue(target, '\n');
    write(target, record.pair);
}

// ----------------------------------------------------------------------------
// Function writeHeader(); RnaHeader
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeHeader(TTarget & target, RnaHeader const & header, DotBracket const & /*tag*/)
{   
    write(target, header.name);
    writeValue(target, '\n');
    write(target, )

}

} //namespace seqan

#endif // SEQAN_DOTBRACKET_FORMAT_IO_H