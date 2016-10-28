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
// Author: Lily Shellhammer <lily.shellhammer@gmail.com>
// ==========================================================================
// This file contains routines to read and write to connect format files (.ct)
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_READ_WRITE_H_

#include <seqan/stream.h>


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
		- Natural numbering. Rnastructure ignores the actual value given in natural numbering, 
			so it is easiest to repeat n here.

CT Files can hold multiple structures of a single sequence.
This is done by repeating the format for each structure without any blank lines between structures. 

record
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
// Tag Connect
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
// Function readRecord(); RnaRecord, Connect
// ----------------------------------------------------------------------------
template <typename TForwardIter>
inline void 
readRecord(RnaRecord & record, RnaIOContext & context, TForwardIter & iter, Connect const & /*tag*/)
{
    //First line record example:
    //  73 ENERGY =     -17.50    S.cerevisiae_tRna-PHE

    clear(record);
    clear(context.buffer); 
    skipUntil(iter, NotFunctor<IsWhitespace>()); 
    readUntil(context.buffer, iter, IsWhitespace());
    if (!lexicalCast(record.amount, context.buffer))
        throw BadLexicalCast(record.amount, context.buffer);
    clear(context.buffer);

    skipUntil(iter, NotFunctor<IsWhitespace>());
    skipUntil(iter, IsWhitespace());
    skipUntil(iter, EqualsChar<'='>());
    skipOne(iter);
    skipUntil(iter, NotFunctor<IsWhitespace>());
    
    readUntil(context.buffer, iter, IsWhitespace());
    if (!lexicalCast(record.energy, context.buffer))
        throw BadLexicalCast(record.energy, context.buffer);
    clear(context.buffer);

    readUntil(context.buffer, iter, NotFunctor<IsWhitespace>());

    readUntil(record.name,  iter, IsNewline());  

    /* 
    Example records:

     3  G           2       4       70          3
     N  SEQUENCE   N-1     N+1    J POSITION  N  
    */
    clear(context);
    Rna5String rec_base;
    unsigned counter = 0;
    skipUntil(iter, NotFunctor<IsWhitespace>()); 
    while (!atEnd(iter))
    {
        skipUntil(iter, NotFunctor<IsWhitespace>());          //skip beginning positiona and whitespace
        

        readUntil(context.buffer, iter, IsWhitespace());
        if (!lexicalCast(context.number, context.buffer))
            throw BadLexicalCast(context.number, context.buffer);
        clear(context.buffer);
        if (counter == 0)
        {
            record.begPos = context.number;
            if(record.amount != 0)  //in case amount is zero, leave it blank until counter is finshed
                record.endPos = record.amount - record.begPos + 1;
        }



        skipUntil(iter, NotFunctor<IsWhitespace>());    

        readUntil(context.base, iter, IsWhitespace());               //read base 
        append(rec_base, context.base);
        clear(context.base);
        addVertex(record.graph);

        skipUntil(iter, NotFunctor<IsWhitespace>());    //skip whitespace
        skipUntil(iter, IsWhitespace());
        skipUntil(iter, NotFunctor<IsWhitespace>());
        skipUntil(iter, IsWhitespace());
        skipUntil(iter, NotFunctor<IsWhitespace>());
        
        readUntil(context.buffer, iter, IsWhitespace());  //read pair position
        if (!lexicalCast(context.number, context.buffer))
            throw BadLexicalCast(context.number, context.buffer);
        clear(context.buffer);
        if (context.number != 0 && counter > context.number - 1)        // add edge if base is connected
            addEdge(record.graph, context.number - 1, counter, 1.);
        clear(context.pair);

        skipUntil(iter, IsNewline());      //skip until newline
        counter++;
    }
    appendValue(record.sequence, rec_base);

    if (record.amount == 0)
    {
        record.amount = counter;
        if (record.begPos != 1)
            record.endPos = record.amount - record.begPos + 1;
    }

}


// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, Connect
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, Connect const & /*tag*/)     
{
    //write old "header"
    if(record.amount != 0)
        appendNumber(target, record.amount);
    else
    {
        std::cerr << "ERROR. No amount of records specified";
        return;
    }
    writeValue(target, ' ');
    write(target, "ENERGY = ");
    writeValue(target, '\t');
    appendNumber(target, record.energy);
    writeValue(target, '\t');
    write(target, record.name);
    writeValue(target, '\n');
    //write "body"
    int begPos = 1;
    if(record.begPos != 1)
        begPos = record.begPos;
    for (unsigned i = 0; i < length(record.sequence[0]); i++)
    {
        writeValue(target, ' ');    //All records start with a space
        appendNumber(target, i+begPos);
        writeValue(target, ' ');
        write(target, record.sequence[0][i]);
        writeValue(target, '\t');
        appendNumber(target, i+begPos-1);
        writeValue(target, '\t');
        appendNumber(target, i+begPos+1);
        writeValue(target, '\t');
        if (degree(record.graph, i) != 0)
        {
            TAdjacencyIterator adj_it(record.graph, i);
            write(target, value(adj_it) + 1);
        }
        else
        {
            writeValue(target, '0');
        }
        writeValue(target, '\t');
        appendNumber(target, i+begPos);
        writeValue(target, '\n');
    }
}
 


} //namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_READ_WRITE_H_
