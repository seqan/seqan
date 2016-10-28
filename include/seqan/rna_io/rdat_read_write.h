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
// This file contains routines to read and write to Rdat format files (.rdat)
// ==========================================================================

#ifndef SEQAN_RDAT_FORMAT_READ_H
#define SEQAN_RDAT_FORMAT_READ_H

#include <seqan/stream.h>


/* IMPLEMENTATION NOTES



*/
namespace seqan{

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// ============================================================================
// Forwards
// ============================================================================
// --------------------------------------------------------------------------
// Tag Rdat
// --------------------------------------------------------------------------

struct Rdat_;
typedef Tag<Rdat_> Rdat;

template <typename T>
struct MagicHeader<Rdat, T> :
    public MagicHeader<Nothing, T> {};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Rdat, T>
{
    static char const * VALUE[1];
};
template <typename T>
char const * FileExtensions<Rdat, T>::VALUE[1] =
{
    ".rdat"      // default output extension
};


// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function readHeader(); RdatHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readHeader(RNAHeader & header, RNAIOContext & context, TForwardIter & iter, Rdat const & /*tag*/)
{
    typedef OrFunctor<IsTab, IsNewline> TNextEntry;

    //RDAT_VERSION    0.32
    //NAME        MedLoop

    clear(header);
    clear(context.buffer); 
    skipUntil(iter, NotFunctor<IsWhitespace>()); 

    readUntil(header.name,  iter, IsNewline());   
    
}

//SO ONE ONLY PASSES RECORD WHICH I PREFER, THIS ONE EXPLICITLY PASSES EVERY VALUE AND I DON'T UNDERSAND WHY NECESSARY

template <typename TForwardIter>
inline void 
readRecord(RNARecord & record, RNAIOContext & context, TForwardIter & iter, Rdat const & /*tag*/)
{
    typedef OrFunctor<IsTab, IsNewline> > NextEntry;

    clear(context);
    clear(record);

    skipUntil(iter, IsWhitespace());
    skipUntil(iter, NextEntry());
    Rna5String rec_base;
    readUntil(rec_base, iter, IsWhitespace());
    appendValue(record.sequence, rec_base);
    skipUntil(iter, IsNewline());

    skipUntil(iter, IsWhitespace());
    skipUntil(iter, NextEntry());
    readUntil(record.pair, iter, IsWhitespace()); // TODO: read connections properly, create graph
    skipUntil(iter, IsNewline());
    // OFFSET
    readUntil(context.buffer, iter, IsWhitespace());
    if (context.buffer == "OFFSET"){
        skipUntil(iter, NextEntry());
        readUntil(context.buffer, iter, IsWhitespace());
        if (!lexicalCast(record.offset, context.buffer))
            throw BadLexicalCast(record.offset, context.buffer);
        clear(context.buffer);
        skipUntil(iter, IsNewline());
    }
    else{
        SEQAN_THROW(EmptyFieldError("OFFSET"));
    }

    // SEQPOS
    readUntil(context.buffer, iter, NextEntry());
    if (context.buffer == "SEQPOS"){
        for(unsigned i = 0; i < length(base); ++i){
            skipUntil(iter, NotFunctor<IsWhitespace>());
            readUntil(context.buffer, iter, IsWhitespace());
            append(record.seqpos, context.buffer);
            clear(context.buffer);
            skipUntil(iter, IsNewline());
        }
    }
    else{
        SEQAN_THROW(EmptyFieldError("SEQPOS"));
    }
    //skipUntil(iter, IsNewline());


}


// ----------------------------------------------------------------------------
// Function writeRecord(Rdat);
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RNARecord const & record, Rdat const & /*tag*/)     
{

}

// ----------------------------------------------------------------------------
// Function writeHeader(); RdatHeader
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeHeader(TTarget & target, RNAHeader const & header, Rdat const & /*tag*/)
{

}


} //namespace seqan

#endif // SEQAN_Rdat_FORMAT_READ_H