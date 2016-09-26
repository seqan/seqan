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
    std::string buffer;
    clear(record);
    clear(context);

    // read number of entries (sequence length)
    skipUntil(iter, NotFunctor<IsWhitespace>());
    readUntil(buffer, iter, IsWhitespace());
    if (!lexicalCast(record.seqLen, buffer))
        throw BadLexicalCast(record.seqLen, buffer);
    clear(buffer);

    //read energy
    skipUntil(iter, NotFunctor<IsWhitespace>());
    skipUntil(iter, IsWhitespace());
    skipUntil(iter, EqualsChar<'='>());
    skipOne(iter);
    skipUntil(iter, NotFunctor<IsWhitespace>());
    readUntil(buffer, iter, IsWhitespace());
    if (!lexicalCast(record.energy, buffer))
        throw BadLexicalCast(record.energy, buffer);
    clear(buffer);

    // read name
    readUntil(buffer, iter, NotFunctor<IsWhitespace>());
    readUntil(record.name,  iter, IsNewline());
    clear(buffer);

    /* 
    Example records:

     3  G           2       4       70          3
     N  SEQUENCE   N-1     N+1    J POSITION  N  
    */

    // read nucleotides with pairs
    TRnaRecordGraph graph;
    unsigned currPos {0};
    while (!atEnd(iter))
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
                throw BadLexicalCast(record.offset, buffer);
            currPos = record.offset;
            clear(buffer);
        }
        skipUntil(iter, NotFunctor<IsWhitespace>());

        // nucleotide
        readUntil(buffer, iter, IsWhitespace());
        append(record.sequence, buffer);
        clear(buffer);
        addVertex(graph);

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
            throw BadLexicalCast(pairPos, buffer);
        if (pairPos != 0 && currPos > pairPos)
        {
            if (pairPos >= record.offset)
                addEdge(graph, pairPos - record.offset, currPos - record.offset, 1.0);
            else
                throw std::runtime_error("ERROR: Incompatible pairing position in input file.");
        }

        clear(buffer);
        skipUntil(iter, IsNewline());
    }
    append(record.graph, graph);
    SEQAN_ASSERT_EQ(record.seqLen, length(record.sequence));
}


// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, Connect
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, Connect const & /*tag*/)     
{
    if (empty(record.sequence) && length(rows(record.align)) != 1)
        throw std::runtime_error("ERROR: Connect formatted file cannot contain an alignment.");
    if (length(record.graph) != 1)
        throw std::runtime_error("ERROR: Connect formatted file cannot contain multiple structure graphs.");

    Rna5String const sequence = empty(record.sequence) ? source(row(record.align, 0)) : record.sequence;
    
    //write old "header"
    appendNumber(target, record.seqLen);
    writeValue(target, '\t');
    write(target, "ENERGY = ");
    appendNumber(target, record.energy);
    writeValue(target, '\t');
    write(target, record.name);
    writeValue(target, '\n');

    //write "body"
    unsigned offset = record.offset > 0 ? record.offset : 1;
    for (unsigned i = 0; i < length(sequence); ++i)
    {
        writeValue(target, ' ');    //All records start with a space
        appendNumber(target, i + offset);
        writeValue(target, '\t');
        write(target, sequence[i]);
        writeValue(target, '\t');
        appendNumber(target, i + offset - 1);
        writeValue(target, '\t');
        appendNumber(target, i + offset + 1);
        writeValue(target, '\t');
        if (degree(record.graph[0], i) != 0)
        {
            TRnaAdjacencyIterator adjIter(record.graph[0], i);
            write(target, value(adjIter) + offset);
        }
        else
        {
            writeValue(target, '0');
        }
        writeValue(target, '\t');
        appendNumber(target, i + offset);
        writeValue(target, '\n');
    }
}

} //namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_CONNECT_READ_WRITE_H_
