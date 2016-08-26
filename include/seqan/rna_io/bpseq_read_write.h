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
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_READ_BPSEQ_H_
#define SEQAN_INCLUDE_SEQAN_BPSEQ_READ_BPSEQ_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Bpseq
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Bpseq
 * @headerfile <seqan/bpseq_io.h>
 * @brief Variant callinf format file.
 *
 * @signature typedef Tag<Bpseq_> Bpseq;
 */
struct Bpseq_;
typedef Tag<Bpseq_> Bpseq;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BpseqHeader]
// ----------------------------------------------------------------------------


template <typename TForwardIter>
inline void
readHeader(RnaHeader & header,
           RnaIOContext & context,
           TForwardIter & iter,
           Bpseq const & /*tag*/)
{
    clear(header);
    CharString buffer;
    RnaHeaderRecord record;

    while (!atEnd(iter) && value(iter) == '#') // All the information stored in the # lines are saved in a single line
    {
        skipOne(iter);
        clear(buffer);
        // Write header key.
        record.key = "Info";
        // Read header value.
        readLine(record.value, iter);
    }
    appendValue(header, record);
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BpseqRecord]
// ----------------------------------------------------------------------------
// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record,
           RnaIOContext & context,
           TForwardIter & iter,
           Bpseq const & /*tag*/)
{
    typedef OrFunctor<IsSpace, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bpseq> > NextEntry;
    clear(record);
    CharString &buffer = context.buffer;
    CharString tmpStr="";
    unsigned counter = 0;
    while (!atEnd(iter) && value(iter) != '#')
    {
        // 
        clear(buffer);
        readUntil(buffer, iter, NextEntry());
        if (empty(buffer))
            SEQAN_THROW(EmptyFieldError("BEGPOS"));
        if(counter == 0)
            record.begPos = lexicalCast<__int32>(buffer);
        skipUntil(iter, NextEntry());

        // SEQUENCE
        clear(buffer);
        readUntil(buffer, iter, IsWhitespace());
        appendValue(record.sequence, buffer[0]);        //TO DO: APPEND VALUES TO SEQUENCE WITHOUT GETTING A SEG FAULT
        if (empty(record.sequence))
            SEQAN_THROW(EmptyFieldError("SEQUENCE"));     

        skipUntil(iter, NextEntry());

        // PAIR
        clear(buffer);
        readUntil(buffer, iter, OrFunctor<IsSpace, IsNewline>());
        if (empty(buffer))
            SEQAN_THROW(EmptyFieldError("PAIR"));
        appendValue(record.pair, lexicalCast<unsigned>(buffer));
        skipLine(iter);

        counter++;
    }
    if(record.begPos != 1)      //set beginning record position
        record.endPos = counter - record.begPos + 1;
    else
        record.endPos = counter;    //set end record position
    record.amount = record.endPos - record.begPos + 1;  //set amount of records

    return;
}

// ----------------------------------------------------------------------------
// Function writeHeader()                                           [BpseqHeader]
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeHeader(TTarget & target,
            RnaHeader const & header,
            RnaIOContext & context,
            Bpseq const & /*tag*/)
{
    for (unsigned i = 0; i < length(header); ++i)
    {
        write(target, "#");
        write(target, header[i].key);
        writeValue(target, '=');
        write(target, header[i].value);
        writeValue(target, '\n');
    }
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                           [BpseqRecord]
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target,
            RnaRecord const & record,
            RnaIOContext & context,
            Bpseq const & /*tag*/)
{
    for (unsigned i = 0; i < record.amount; ++i)
    {
        write(target, record.begPos + i);
        writeValue(target, ' ');
        write(target, record.sequence[0][i]);  //TO DO: PRINT SEQUENCE AT POSITION I WITHOUT SEG FAULTING
        writeValue(target, ' ');
        write(target, record.pair[i]);
        writeValue(target, '\n');
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_READ_BPSEQ_H_
