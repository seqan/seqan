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
// This file contains routines to read and write to connect format files (.ct)
// ==========================================================================

#ifndef SEQAN_RNA_FORMAT_READ_H
#define SEQAN_RNA_FORMAT_READ_H

#include <seqan/stream.h>


/* IMPLEMENTATION NOTES

RNA FORMAT example:

=> HEADER START : number of bases in the sequence
=> HEADER END: title of the structure
=> Each line has information about a base pair in the sequence
	Each line is a base, with this order of information:
		- Base number: index n
		- Base (A, C, G, T, U, X)
		- Index n-1
		- Index n+1
		- Number of the base to which n is paired. No pairing is indicated by 0 (zero).
		- Natural numbering. RNAstructure ignores the actual value given in natural numbering, 
			so it is easiest to repeat n here.

CT Files can hold multiple structures of a single sequence.
This is done by repeating the format for each structure without any blank lines between structures. 

HEADER
 N  SEQUENCE   N-1  	 N+1	J POSITION  N  
 1 	G       	0    	2   	72    		1
 2 	C       	1    	3   	71    		2
 3 	G       	2    	4   	70    		3


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

struct Connect_;
typedef Tag<Connect_> Connect;

template <typename T>
struct MagicHeader<Connect, T> :
    public MagicHeader<Nothing, T> {};

// ============================================================================
// Metafunctions
// ============================================================================

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

// ----------------------------------------------------------------------------
// Function readHeader(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readHeader(RNAHeader & header, RNAIOContext & context, TForwardIter & iter, Connect const & /*tag*/)
{
    typedef OrFunctor<IsTab, IsNewline> TNextEntry;

    //  73 ENERGY =     -17.50    S.cerevisiae_tRNA-PHE

    clear(header);
    clear(context.buffer); 
    skipUntil(iter, NotFunctor<IsWhitespace>()); 
    readUntil(context.buffer, iter, IsWhitespace());
    if (!lexicalCast(header.amount, context.buffer))
        throw BadLexicalCast(header.amount, context.buffer);
    clear(context.buffer);

    skipUntil(iter, NotFunctor<IsWhitespace>());
    skipUntil(iter, IsWhitespace());
    skipUntil(iter, EqualsChar<'='>());
    skipOne(iter);
    skipUntil(iter, NotFunctor<IsWhitespace>());
    skipOne(iter);
    
    readUntil(context.buffer, iter, IsWhitespace());
    if (!lexicalCast(header.energy, context.buffer))
        throw BadLexicalCast(header.energy, context.buffer);
    clear(context.buffer);

    readUntil(context.buffer, iter, NotFunctor<IsWhitespace>());

    readUntil(header.name,  iter, IsNewline());   
    
}


template <typename TForwardIter>
inline void 
readRecord(RNARecord & record, RNAIOContext & context, TForwardIter & iter, Connect const & /*tag*/)
{
// 3  G           2       4       70          3
// N  SEQUENCE   N-1     N+1    J POSITION  N  

    clear(context);
    unsigned counter = 0;
    skipUntil(iter, NotFunctor<IsWhitespace>()); 
    while (!atEnd(iter))
    {
        skipUntil(iter, NotFunctor<IsWhitespace>());          //skip beginning positiona and whitespace
        

        readUntil(context.buffer, iter, IsWhitespace());
        if (!lexicalCast(context.number, context.buffer))
            throw BadLexicalCast(context.number, context.buffer);
        clear(context.buffer);
        append(record.index, context.number);

        skipUntil(iter, NotFunctor<IsWhitespace>());    

        readUntil(context.base, iter, IsWhitespace());               //read base 
        append(record.base, context.base);
        clear(context.base);

        skipUntil(iter, NotFunctor<IsWhitespace>());    //skip whitespace
        skipUntil(iter, IsWhitespace());
        skipUntil(iter, NotFunctor<IsWhitespace>());
        skipUntil(iter, IsWhitespace());
        skipUntil(iter, NotFunctor<IsWhitespace>());
        
        readUntil(context.buffer, iter, IsWhitespace());  //read pair position
        if (!lexicalCast(context.number, context.buffer))
        throw BadLexicalCast(context.number, context.buffer);
        clear(context.buffer);
        append(record.pair, context.number);
        clear(context.pair);

        skipUntil(iter, IsNewline());      //skip until newline
        counter++;
    }
}


// ----------------------------------------------------------------------------
// Function writeRecord(RNA);
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RNARecord const & record, Connect const & /*tag*/)     
{
    for (unsigned i = 0; i < length(record.base); i++)
    {
        writeValue(target, ' ');    //All records start with a space
        appendNumber(target, record.index[i]);
        writeValue(target, ' ');
        write(target, record.base[i]);
        writeValue(target, '\t');
        appendNumber(target, record.index[i]-1);
        writeValue(target, '\t');
        appendNumber(target, record.index[i]+1);
        writeValue(target, '\t');
        write(target, record.pair[i]);
        writeValue(target, '\t');
        appendNumber(target, record.index[i]);
        writeValue(target, '\n');
    }
}

// ----------------------------------------------------------------------------
// Function writeHeader(); RNAHeader
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeHeader(TTarget & target, RNAHeader const & header, Connect const & /*tag*/)
{
    appendNumber(target, header.amount);
    writeValue(target, ' ');
    write(target, "ENERGY = ");
    writeValue(target, '\t');
    appendNumber(target, header.energy);
    writeValue(target, '\t');
    write(target, header.name);
    writeValue(target, '\n');
}


} //namespace seqan

#endif // SEQAN_RNA_FORMAT_READ_H