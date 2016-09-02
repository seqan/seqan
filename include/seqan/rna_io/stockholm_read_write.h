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
// Author: Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// This file contains routines to read and write to connect format files (.ct)
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_STOCKHOLM_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_STOCKHOLM_READ_WRITE_H_

#include <seqan/stream.h>

/* IMPLEMENTATION NOTES

Stockholm FORMAT:

=> header with format version
=> tags specify semantic of line
        #=GR SS         means secondary structure
        #=GC SS_cons    means secondary structure of consensus
=> allowed symbols for RNA  .,;         unpaired
                            <>(){}[]    paired
                            AaBb...     pseudoknot
=> sequence lines consist of seqname and the aligned sequence
=> EOF is marked by //
=> example

    # STOCKHOLM 1.0
    #=GF ID    trna
    #=GF DE    Taken from Sprinzl alignment of 1415 tRNAs [Steinberg93]

    DF6280             GCGGAUUUAGCUCAGUUGGG.AGAGCGCCAGACUGAAGAUCUGGAGGUCC
    DE6280             UCCGAUAUAGUGUAAC.GGCUAUCACAUCACGCUUUCACCGUGGAGA.CC
    DD6280             UCCGUGAUAGUUUAAU.GGUCAGAAUGGGCGCUUGUCGCGUGCCAGA.UC
    DC6280             GCUCGUAUGGCGCAGU.GGU.AGCGCAGCAGAUUGCAAAUCUGUUGGUCC
    DA6280             GGGCACAUGGCGCAGUUGGU.AGCGCGCUUCCCUUGCAAGGAAGAGGUCA
    #=GC SS_cons       <<<<<<<..<<<<.........>>>>.<<<<<.......>>>>>.....<
    #=GC RF            xxxxxxxxxxxxxxxxxxxx.xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    DF6280             UGUGUUCGAUCCACAGAAUUCGCA
    DE6280             GGGGUUCGACUCCCCGUAUCGGAG
    DD6280             GGGGUUCAAUUCCCCGUCGCGGAG
    DC6280             UUAGUUCGAUCCUGAGUGCGAGCU
    DA6280             UCGGUUCGAUUCCGGUUGCGUCCA
    #=GC SS_cons       <<<<.......>>>>>>>>>>>>.
    #=GC RF            xxxxxxxxxxxxxxxxxxxxxxxx
    //

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

struct Stockholm_;
typedef Tag<Stockholm_> Stockholm;

template <typename T>
struct MagicHeader<Stockholm, T> :
    public MagicHeader<Nothing, T> {};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Stockholm, T>
{
    static char const * VALUE[1];
};
template <typename T>
char const * FileExtensions<Stockholm, T>::VALUE[1] =
{
    ".sth"      // default output extension
};


// ==========================================================================
// Functions
// ==========================================================================


// ----------------------------------------------------------------------------
// Function readRecord(); RnaRecord, Stockholm
// ----------------------------------------------------------------------------
template <typename TForwardIter>
inline void 
readRecord(RnaRecord & record, RnaIOContext & context, TForwardIter & iter, Stockholm const & /*tag*/)
{
    clear(record);
    clear(context.buffer);
    skipUntil(iter, NotFunctor<IsWhitespace>());
    readUntil(context.buffer, iter, IsNewline());               // Intro # STOCKHOLM 1.0
    if (std::string(toCString(context.buffer)).find("STOCKHOLM") == std::string::npos)
        throw ParseError("Expected STOCKHOLM identifier in the first line.");
    skipOne(iter);
    CharString bracket_str = "";
    StringSet<Rna5String, Owner<JournaledSet> > sequence;
    StringSet<String<unsigned> > gapPos;

    while (!atEnd(iter))
    {
        clear(context.buffer);
        readUntil(context.buffer, iter, IsWhitespace());        // read line until first whitespace

        if (context.buffer == "//")
        {                                                       // found terminal symbols
            break;
        }
        else if (context.buffer == "#=GC" || context.buffer == "#=GR")
        {                                                       // found a tag
            skipOne(iter);                                      // skip whitespace
            clear(context.buffer);
            readUntil(context.buffer, iter, IsWhitespace());    // read tag
            if (context.buffer == "SS" || context.buffer == "SS_cons")
            {                                                   // found secondary structure
                clear(context.buffer);
                skipUntil(iter, NotFunctor<IsWhitespace>());
                readUntil(context.buffer, iter, IsWhitespace());
                append(bracket_str, context.buffer);
            }
        }
        else if (context.buffer == "#=GF")
        {                                                       // found a tag
            skipOne(iter);                                      // skip whitespace
            clear(context.buffer);
            readUntil(context.buffer, iter, IsWhitespace());    // read tag
            if (context.buffer == "ID")
            {                                                   // found identification tag
                clear(context.buffer);
                skipUntil(iter, NotFunctor<IsWhitespace>());
                readUntil(context.buffer, iter, IsWhitespace());
                record.name = context.buffer;
            }
        }
        else if (length(context.buffer) > 0 && context.buffer[0] != '#')
        {                                                       // found sequence id
            unsigned idx = 0;
            while (idx < length(record.seq_id) && record.seq_id[idx] != context.buffer)
                ++idx;                                          // search sequence id

            unsigned const OFFSET = idx == length(record.seq_id) ? 0 : length(sequence[idx]);
            if (OFFSET == 0)
            {
                appendValue(record.seq_id, context.buffer);
                resize(rows(record.align), idx + 1);
                resize(gapPos, idx + 1);
            }

            clear(context.buffer);
            skipUntil(iter, NotFunctor<IsWhitespace>());
            readUntil(context.buffer, iter, IsWhitespace());    // read sequences

            unsigned gap = 0;
            while (gap < length(context.buffer))                // remove gap symbols
            {
                if (context.buffer[gap] != '.' && context.buffer[gap] != '-')
                {
                    ++gap;
                }
                else
                {
                    appendValue(gapPos[idx], gap + OFFSET);
                    erase(context.buffer, gap);
                }
            }

            if (OFFSET == 0)                                    // new sequence
                appendValue(sequence, context.buffer);
            else                                                // append to existing sequence
                append(sequence[idx], context.buffer);
        }
        skipUntil(iter, IsNewline());
        skipOne(iter);
    }
    clear(context.buffer);

    for (unsigned idx = 0; idx < length(gapPos); ++idx)
    {
        assignSource(row(record.align, idx), sequence[idx]);
        for (unsigned gap = length(gapPos[idx]); gap > 0u; --gap)
            insertGap(row(record.align, idx), gapPos[idx][gap - 1]);
    }

    record.amount = length(bracket_str);
    if (record.amount == 0)
        throw ParseError("Expected a secondary structure line.");
    if (empty(sequence))
        throw ParseError("Expected a sequence line.");
    bracket2graph(record.graph, bracket_str);
    record.begPos = 1;
    record.endPos = record.amount;
}


// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, Stockholm
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, Stockholm const & /*tag*/)     
{
    write(target, "# STOCKHOLM 1.0\n");                         // header
    if (record.name != " ")
    {                                                           // record name
        write(target, "#=GF ID      ");
        write(target, record.name);
        writeValue(target, '\n');
    }
    writeValue(target, '\n');
    if (!empty(record.sequence))
    {                                                           // single sequence
        write(target, "single-seq  \t");
        write(target, record.sequence);
        write(target, "\n#=GR SS     \t");
    }
    else
    {                                                           // alignment
        SEQAN_ASSERT_EQ(length(record.seq_id), length(rows(record.align)));
        unsigned maxlen = 12;                                   // determine indentation
        for (unsigned idx = 0; idx < length(record.seq_id); ++idx)
            maxlen = length(record.seq_id[idx]) > maxlen ? length(record.seq_id[idx]) : maxlen;
        for (unsigned idx = 0; idx < length(record.seq_id); ++idx)
        {
            CharString str = record.seq_id[idx];                // write id
            resize(str, maxlen, ' ');
            write(target, str);
            write(target, "\t");
            write(target, row(record.align, idx));              // write sequence
            writeValue(target, '\n');
        }
        write(target, "#=GC SS_cons\t");
    }

    CharString bracket_str("");
    graph2bracket(bracket_str, record.graph);
    write(target, bracket_str);
    write(target, "\n//\n");                                    // closing
}

} //namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_STOCKHOLM_READ_WRITE_H_
