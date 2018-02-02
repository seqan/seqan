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
// Author: Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// This file contains routines to read and write to connect format files (.ct)
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_STOCKHOLM_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_STOCKHOLM_READ_WRITE_H_

#include <seqan/stream.h>
// for bracket-graph transformation
#include <seqan/rna_io/dot_bracket_read_write.h>
#include "rna_record.h"

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

// --------------------------------------------------------------------------
// Tag Connect
// --------------------------------------------------------------------------

/*!
 * @tag FileFormats#Stockholm
 * @headerfile <seqan/rna_io.h>
 * @brief Stockholm format for RNA structures (*.sth).
 * @signature typedef Tag<Stockholm_> Stockholm;
 * @see FileFormats#RnaStruct
 */
struct Stockholm_;
typedef Tag<Stockholm_> Stockholm;

// --------------------------------------------------------------------------
// Class MagicHeader
// --------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Stockholm, T> :
    public MagicHeader<Nothing, T> {};

// ==========================================================================
// Metafunctions
// ==========================================================================

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
readRecord(RnaRecord & record, TForwardIter & iter, Stockholm const & /*tag*/)
{
    std::string buffer;
    clear(record);

    // read intro: # STOCKHOLM 1.0
    skipUntil(iter, NotFunctor<IsWhitespace>());
    readUntil(buffer, iter, IsNewline());
    if (buffer.find("STOCKHOLM") == std::string::npos)
        SEQAN_THROW(ParseError("Expected STOCKHOLM identifier in the first line."));
    skipOne(iter);
    clear(buffer);

    CharString bracketStr{};
    StringSet<Rna5String, Owner<JournaledSet> > sequence;
    StringSet<String<TSizeRna5String> > gapPos;

    while (!atEnd(iter))
    {
        // read tags and their values
        readUntil(buffer, iter, IsWhitespace());                // read line until first whitespace
        if (buffer == "//")
        {                                                       // found terminal symbols
            skipLine(iter);
            clear(buffer);
            break;
        }
        else if (buffer == "#=GC" || buffer == "#=GR")
        {                                                       // found a GC or GR tag
            clear(buffer);
            skipOne(iter);
            readUntil(buffer, iter, IsWhitespace());            // read tag
            if (buffer == "SS" || buffer == "SS_cons")
            {                                                   // found secondary structure
                clear(buffer);
                skipUntil(iter, NotFunctor<IsWhitespace>());
                readUntil(buffer, iter, IsWhitespace());
                append(bracketStr, buffer);
            }
        }
        else if (buffer == "#=GF")
        {                                                       // found a GF tag
            clear(buffer);
            skipOne(iter);
            readUntil(buffer, iter, IsWhitespace());            // read tag
            if (buffer == "ID")
            {                                                   // found identification
                clear(buffer);
                skipUntil(iter, NotFunctor<IsWhitespace>());
                readUntil(buffer, iter, IsWhitespace());
                record.name = buffer;
            }
            else if (buffer == "DE")
            {                                                   // found description
                clear(buffer);
                skipUntil(iter, NotFunctor<IsWhitespace>());
                readUntil(buffer, iter, IsWhitespace());
                record.comment = buffer;
            }
        }
        else if (length(buffer) > 0u && buffer[0] != '#')
        {                                                       // found a sequence id
            typename Size<StringSet<CharString> >::Type idx = 0u;  // search for index of this seq id
            while (idx < length(record.seqID) && record.seqID[idx] != buffer)
                ++idx;

            // determine length of stored sequence belonging to idx
            TSizeRna5String const offset = idx == length(record.seqID) ? 0u : length(sequence[idx]);
            if (offset == 0u)
            {                                                   // sequence id is new
                appendValue(record.seqID, buffer);
                resize(gapPos, idx + 1);
            }
            clear(buffer);

            // read sequence
            skipUntil(iter, NotFunctor<IsWhitespace>());
            readUntil(buffer, iter, IsWhitespace());

            // remove gap symbols from sequence and store their positions in gapPos list
            typename Size<std::string>::Type pos = 0u;
            while (pos < length(buffer))
            {
                if (buffer[pos] != '.' && buffer[pos] != '-')
                {
                    ++pos;
                }
                else
                {
                    appendValue(gapPos[idx], static_cast<TSizeRna5String>(pos) + offset);
                    erase(buffer, pos);
                }
            }

            if (offset == 0u)
                appendValue(sequence, buffer);                  // store new sequence
            else
                append(sequence[idx], buffer);                  // append to existing sequence
        }
        skipUntil(iter, IsNewline());
        skipOne(iter);
        clear(buffer);
    }
    if (empty(bracketStr))
        SEQAN_THROW(ParseError("Expected a secondary structure line."));
    if (empty(sequence))
        SEQAN_THROW(ParseError("Expected a sequence line."));

    SEQAN_ASSERT_EQ(length(sequence), length(gapPos));

    // store the alignment with gaps in record
    resize(rows(record.align), length(sequence));
    for (typename Size<StringSet<Rna5String, Owner<JournaledSet> > >::Type seq = 0; seq < length(sequence); ++seq)
    {
        assignSource(row(record.align, seq), sequence[seq]);
        for (typename Size<String<TSizeRna5String> >::Type pos = length(gapPos[seq]); pos > 0u; --pos)
        {
            insertGap(row(record.align, seq), gapPos[seq][pos - 1]);
        }
    }

    // store sequence length and secondary structure graph
    record.seqLen = length(bracketStr);
    bracket2graph(record.fixedGraphs, bracketStr);
}

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, RnaIOContext & /*context*/, TForwardIter & iter, Stockholm const & /*tag*/)
{
    readRecord(record, iter, Stockholm());
}

// ----------------------------------------------------------------------------
// Function writeRecord(); RnaRecord, Stockholm
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, Stockholm const & /*tag*/)
{
    if (length(record.fixedGraphs) != 1u)
        SEQAN_THROW(ParseError("ERROR: Cannot deal with multiple structure graphs."));

    write(target, "# STOCKHOLM 1.0\n");                         // header
    if (!empty(record.name))
    {                                                           // name
        write(target, "#=GF ID      ");
        write(target, record.name);
        writeValue(target, '\n');
    }
    if (!empty(record.comment))
    {                                                           // description
        write(target, "#=GF DE      ");
        write(target, record.comment);
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
        SEQAN_ASSERT_EQ(length(record.seqID), length(rows(record.align)));
        typename Size<CharString>::Type maxlen = 12;            // determine indentation
        forEach(record.seqID, [&](CharString const & sID)
        {
            maxlen = length(sID) > maxlen ? length(sID) : maxlen;
        });

        for (typename Size<StringSet<CharString> >::Type idx = 0; idx < length(record.seqID); ++idx)
        {
            CharString str = record.seqID[idx];                 // write id
            resize(str, maxlen, ' ');
            write(target, str);
            writeValue(target, '\t');
            write(target, row(record.align, idx));              // write sequence
            writeValue(target, '\n');
        }
        write(target, "#=GC SS_cons\t");
    }

    write(target, graph2bracket(record.fixedGraphs[0]));
    write(target, "\n//\n");                                    // closing
}

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, RnaIOContext & /*context*/, Stockholm const & /*tag*/)
{
    writeRecord(target, record, Stockholm());
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_STOCKHOLM_READ_WRITE_H_
