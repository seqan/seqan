// ==========================================================================
//                         Mason - A Read Simulator
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// TSV format for variation sizes.  First column gives type, second column
// gives the variation size.
// ==========================================================================

#ifndef APPS_MASON2_VARIATION_SIZE_TSV_H_
#define APPS_MASON2_VARIATION_SIZE_TSV_H_

#include <seqan/sequence.h>
#include <seqan/stream.h>

// TODO(holtgrew): Consistently use exceptions instead of return values in mason.

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct VariationSizeTsv_;
typedef seqan::Tag<VariationSizeTsv_> VariationSizeTsv;

class VariationSizeRecord
{
public:
    enum Kind
    {
        INVALID = 0,
        INDEL = 1,
        INVERSION = 2,
        TRANSLOCATION = 3,
        DUPLICATION = 4
    };

    // Kind of the operation.
    Kind kind;

    // Size of the operation.  Always positive, also in case of deletions.
    int size;

    // The characters to insert into the genome.
    seqan::CharString seq;

    VariationSizeRecord() : kind(INVALID), size(-1)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TForwardIterator>
int readRecord(VariationSizeRecord & record, TForwardIterator & iter,
               VariationSizeTsv const & /*tag*/)
{
    seqan::CharString buffer, buffer2;

    // Read until tab, must not read file end.
    readUntil(buffer, iter, seqan::OrFunctor<seqan::IsTab, seqan::IsNewline>());

    // Read until tab, must not reach file end.
    skipOne(iter, seqan::IsTab());

    // Read until tab or line break, may reach end of file.
    readUntil(buffer2, iter, seqan::OrFunctor<seqan::IsTab, seqan::IsNewline>());

    lexicalCastWithException(record.size, buffer2);

    if (buffer == "INS" || buffer == "DEL")
        record.kind = VariationSizeRecord::INDEL;
    else if (buffer == "INV")
        record.kind = VariationSizeRecord::INVERSION;
    else if (buffer == "CTR")
        record.kind = VariationSizeRecord::TRANSLOCATION;
    else if (buffer == "DUP")
        record.kind = VariationSizeRecord::DUPLICATION;
    else
        return 1;  // Invalid record type.

    // Make sure that the record size is valid.
    if (record.size < 0)
        return 1;  // invalid size
    // Adjust record.size to type in case of indel.
    if (buffer == "DEL")
        record.size = -record.size;

    // In case of an insertion, check if we can read the sequence to be inserted.
    if (buffer == "INS" && !atEnd(iter) && *iter == '\t')
    {
        skipOne(iter, seqan::IsTab());

        clear(buffer2);
        readUntil(buffer2, iter, seqan::IsWhitespace());
        record.seq = buffer2;
        record.size = length(record.seq);
    }

    // Skip rest of the line, may reach end of file.
    skipLine(iter);

    return 0;
}

#endif  // #ifndef APPS_MASON2_VARIATION_SIZE_TSV_H_
